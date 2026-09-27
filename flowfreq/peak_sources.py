"""Pluggable sources of annual peak-flow data.

USGS is retiring the legacy NWIS web services (``nwis.waterdata.usgs.gov``,
``waterservices.usgs.gov``) in favour of the Water Data OGC APIs
(``api.waterdata.usgs.gov``). This module puts one interface in front of both,
so callers pick a backend by name and every backend returns the same frame.

Both backends work. The Water Data API backend was written against the live
``peaks`` collection (verified 2026-09-25), per this repository's rule
(``TODO.md``, ``docs/STREAMSTATS_NSS_ADDENDUM.md``) that no client is written
against an endpoint until it has been exercised live. It is the default: the
parity test against ``nwis-legacy`` (``tests/test_peak_backend_parity.py``)
passed live on 2026-09-26 -- same water years, flows and discharge codes, dates
within one day (the API's date is UTC; legacy's is local).

Roadmap: ``docs/MASTER_ROADMAP.md`` §1.1, issue #29.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Iterable, List, Optional, Protocol, Type, runtime_checkable

import pandas as pd
import requests

logger = logging.getLogger(__name__)

#: Columns every backend returns, in this order.
PEAK_COLUMNS = ("water_year", "peak_date", "peak_flow_cfs", "qualification_code")


@runtime_checkable
class PeakDataBackend(Protocol):
    """A source of annual peak-flow records for one site."""

    name: str

    def fetch_peaks(self, site_no: str) -> pd.DataFrame:
        """Return the site's annual peaks with columns :data:`PEAK_COLUMNS`."""
        ...  # pragma: no cover


def validate_peak_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Check a backend's output against the shared contract.

    Parameters
    ----------
    df : pandas.DataFrame
        Backend output.

    Returns
    -------
    pandas.DataFrame
        ``df`` restricted to :data:`PEAK_COLUMNS`, in order.

    Raises
    ------
    ValueError
        On missing columns, duplicate water years, or negative flows.
    """
    missing = [c for c in PEAK_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Peak frame missing columns: {missing}")
    if df["water_year"].duplicated().any():
        dups = sorted(df.loc[df["water_year"].duplicated(), "water_year"].unique())
        raise ValueError(f"Duplicate water years in peak frame: {dups}")
    if (df["peak_flow_cfs"] < 0).any():
        raise ValueError("Negative peak discharge in peak frame")
    return df.loc[:, list(PEAK_COLUMNS)]


class LegacyNwisBackend:
    """The existing NWIS RDB peak service, via :class:`flowfreq.usgs.USGSgage`."""

    name = "nwis-legacy"

    def fetch_peaks(self, site_no: str) -> pd.DataFrame:
        """Download peaks from the legacy NWIS peak service.

        Parameters
        ----------
        site_no : str
            USGS site number.

        Returns
        -------
        pandas.DataFrame
        """
        from flowfreq.usgs import USGSgage  # deferred: usgs imports plotting deps

        return validate_peak_frame(USGSgage(site_no).download_peak_flow())


#: Water Data OGC API ``peaks`` collection items endpoint.
WATERDATA_PEAKS_URL = "https://api.waterdata.usgs.gov/ogcapi/v1/collections/peaks/items"

#: Discharge. The ``peaks`` collection also holds gage-height (00065) rows for
#: the same water years, which would duplicate every year if not filtered out.
DISCHARGE_PARAMETER_CODE = "00060"

#: Water Data API ``qualifier`` tokens and the NWIS peak codes they stand for
#: (checked against the vendored WATSTORE files).
QUALIFIER_TOKEN_CODES: Dict[str, str] = {
    "MAXDAILYMEAN": "1",
    "ESTIMATED": "2",
    "DAMFAILURE": "3",
    "LESSTHAN": "4",
    "UNKNOWNREGULATION": "5",
    "REGULATED": "6",
    "HISTORIC": "7",
    "GREATERTHAN": "8",
    "EVENT": "9",
    "URBAN": "C",
    "OPPORTUNISTIC": "O",
    "REVISED": "R",
}

#: Tokens describing the precision of ``time``, not the discharge. The date
#: itself already carries the placeholder (day 01 / January), as peakfq does.
DATE_PRECISION_TOKENS = frozenset({"DAYUNKNOWN", "MONTHUNKNOWN"})

#: Gage-height flags that can appear on discharge rows; not discharge codes.
GAGE_HEIGHT_TOKENS = frozenset({"GHNOTASSCPKQ"})

_REQUEST_TIMEOUT_S = 60.0
_PAGE_LIMIT = 50_000
_MAX_PAGES = 1_000


def qualifiers_to_codes(
    tokens: Optional[Iterable[str]],
    site_no: str = "",
    water_year: Optional[int] = None,
) -> str:
    """Translate Water Data API ``qualifier`` tokens into an NWIS code string.

    The tokens must never reach :mod:`flowfreq.peak_codes` untranslated. As a
    list, no token matches a single-character code and censoring is silently
    lost. As a string such as ``"DAYUNKNOWN,HISTORIC"``, it is split into
    characters that include ``O`` and ``C``, and the historic peak is removed.

    Parameters
    ----------
    tokens : iterable of str or None
        The feature's ``qualifier`` property.
    site_no : str, optional
        Used only in log messages.
    water_year : int, optional
        Used only in log messages.

    Returns
    -------
    str
        Comma-separated NWIS codes in token order (``"1"``, ``"2,7"``), the
        format of the legacy RDB ``peak_cd`` column and of
        :func:`flowfreq.peak_codes.parse_codes`; ``""`` when there are none.
        Date-precision tokens are dropped, and gage-height tokens are dropped
        with a debug log. An unrecognised token is logged at warning level and
        left out: it has no code, and passing it through would have it read as
        single-character codes.
    """
    if tokens is None:
        return ""
    if isinstance(tokens, str):
        tokens = [tokens]
    codes: List[str] = []
    for raw in tokens:
        token = str(raw).strip().upper()
        if not token or token in DATE_PRECISION_TOKENS:
            continue
        if token in GAGE_HEIGHT_TOKENS:
            logger.debug(
                "Site %s WY %s: gage-height qualifier %s ignored on a discharge peak",
                site_no,
                water_year,
                token,
            )
            continue
        code = QUALIFIER_TOKEN_CODES.get(token)
        if code is None:
            logger.warning(
                "Site %s WY %s: unrecognised Water Data API peak qualifier %r has no "
                "NWIS code and is not carried into qualification_code",
                site_no,
                water_year,
                token,
            )
            continue
        if code not in codes:
            codes.append(code)
    return ",".join(codes)


class WaterDataApiBackend:
    """USGS Water Data OGC API ``peaks`` collection (issue #29).

    Live-verified 2026-09-25. Not the default backend until the parity test
    against ``nwis-legacy`` passes.

    Parameters
    ----------
    timeout : float, default 60
        Per-request timeout, seconds.
    """

    name = "waterdata-ogc"

    def __init__(self, timeout: float = _REQUEST_TIMEOUT_S) -> None:
        self.timeout = timeout

    def _get_page(self, url: str, params: Optional[Dict[str, Any]], site_no: str) -> Dict[str, Any]:
        """GET one page of items and return its decoded JSON."""
        try:
            response = requests.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            payload = response.json()
        except requests.RequestException as exc:
            raise requests.RequestException(
                f"Water Data API peak request failed for site {site_no} ({url}): {exc}"
            ) from exc
        except ValueError as exc:  # body is not JSON
            raise ValueError(
                f"Water Data API returned non-JSON for site {site_no} ({url}): {exc}"
            ) from exc
        if not isinstance(payload, dict):
            raise ValueError(f"Unexpected Water Data API response for site {site_no} ({url})")
        return payload

    def fetch_peaks(self, site_no: str) -> pd.DataFrame:
        """Download annual discharge peaks from the Water Data OGC API.

        Parameters
        ----------
        site_no : str
            USGS site number, with or without the ``USGS-`` prefix. The API
            requires the prefix (without it, it returns no features), so it
            is added here.

        Returns
        -------
        pandas.DataFrame
            Columns :data:`PEAK_COLUMNS`, sorted by water year (the API's
            order is arbitrary). ``water_year`` is the API's own field, never
            derived from the date. ``peak_date`` is the API's ``time``: the
            UTC date when the time of day is known, and a day-01 / January
            placeholder when the day or month is unknown, matching the legacy
            parser. ``qualification_code`` is an NWIS code string; see
            :func:`qualifiers_to_codes`.

        Raises
        ------
        requests.RequestException
            On an HTTP or network failure, naming the site and URL.
        ValueError
            If the site has no discharge peaks, or a response is malformed.
        """
        site = str(site_no).strip()
        if site.upper().startswith("USGS-"):
            site = site[5:]
        params: Optional[Dict[str, Any]] = {
            "f": "json",
            "monitoring_location_id": f"USGS-{site}",
            "parameter_code": DISCHARGE_PARAMETER_CODE,
            "skipGeometry": "true",
            "limit": _PAGE_LIMIT,
        }
        url: Optional[str] = WATERDATA_PEAKS_URL
        features: List[Dict[str, Any]] = []
        pages = 0
        while url is not None:
            pages += 1
            if pages > _MAX_PAGES:
                raise ValueError(
                    f"Water Data API paging for site {site} exceeded {_MAX_PAGES} pages ({url})"
                )
            payload = self._get_page(url, params, site)
            features.extend(payload.get("features") or [])
            # A next link already carries every query parameter plus the cursor.
            url = next(
                (lk.get("href") for lk in payload.get("links") or [] if lk.get("rel") == "next"),
                None,
            )
            params = None

        rows = []
        for feat in features:
            props = feat.get("properties") or {}
            pcode = props.get("parameter_code")
            if pcode is not None and pcode != DISCHARGE_PARAMETER_CODE:
                logger.warning(
                    "Site %s: skipping non-discharge peak row (parameter_code %s)", site, pcode
                )
                continue
            unit = props.get("unit_of_measure")
            if unit is not None and unit != "ft^3/s":
                raise ValueError(
                    f"Site {site}: unexpected peak discharge unit {unit!r} "
                    f"from {WATERDATA_PEAKS_URL}"
                )
            wy = props.get("water_year")
            rows.append(
                {
                    "water_year": wy,
                    "peak_date": props.get("time"),
                    "peak_flow_cfs": props.get("value"),
                    "qualification_code": qualifiers_to_codes(props.get("qualifier"), site, wy),
                }
            )

        if not rows:
            raise ValueError(
                f"No discharge peaks found for site {site} at {WATERDATA_PEAKS_URL} "
                f"(monitoring_location_id=USGS-{site}, "
                f"parameter_code={DISCHARGE_PARAMETER_CODE})"
            )

        df = pd.DataFrame(rows, columns=list(PEAK_COLUMNS))
        df["water_year"] = pd.to_numeric(df["water_year"], errors="coerce")
        df["peak_flow_cfs"] = pd.to_numeric(df["peak_flow_cfs"], errors="coerce")
        df["peak_date"] = pd.to_datetime(df["peak_date"], format="%Y-%m-%d", errors="coerce")
        incomplete = df["water_year"].isna() | df["peak_flow_cfs"].isna()
        if incomplete.any():
            logger.info(
                "Site %s: dropped %d peak row(s) with no water year or no value",
                site,
                int(incomplete.sum()),
            )
            df = df.loc[~incomplete].copy()
        df["water_year"] = df["water_year"].astype(int)
        df = df.sort_values("water_year", kind="stable").reset_index(drop=True)
        return validate_peak_frame(df)


_BACKENDS: Dict[str, Type[PeakDataBackend]] = {
    LegacyNwisBackend.name: LegacyNwisBackend,
    WaterDataApiBackend.name: WaterDataApiBackend,
}

#: Backend used when none is named. Switched from ``nwis-legacy`` once #29's
#: parity test passed live; USGS is retiring the legacy service.
DEFAULT_BACKEND = WaterDataApiBackend.name


def get_backend(name: str = DEFAULT_BACKEND) -> PeakDataBackend:
    """Return a backend instance by name.

    Parameters
    ----------
    name : str
        One of ``"nwis-legacy"`` or ``"waterdata-ogc"``.

    Returns
    -------
    PeakDataBackend

    Raises
    ------
    KeyError
        For an unknown backend name.
    """
    try:
        return _BACKENDS[name]()
    except KeyError:
        raise KeyError(f"Unknown peak backend {name!r}; choose from {sorted(_BACKENDS)}") from None
