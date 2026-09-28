"""USGS Water Data OGC API backend for instantaneous values.

USGS is retiring the legacy NWIS instantaneous-values service
(``waterservices.usgs.gov/nwis/iv``) in favour of the Water Data OGC APIs at
``api.waterdata.usgs.gov``. This module retrieves discharge (00060) and gage
height (00065) from the ``continuous`` collection and returns exactly the frame
the legacy path returns -- a tz-aware UTC index plus ``<value column>``,
``datetime_local``, ``tz_cd`` and ``qualification_code`` -- so
:mod:`flowfreq.subdaily` works on either unchanged.

It is the default backend of ``USGSgage.download_instantaneous_flow`` (and the
``_stage`` twin) since a live parity test against the legacy service
(``tests/test_iv_backend_parity.py``) passed; ``backend="nwis-legacy"`` still
selects the legacy service. Issue #29.

Every endpoint and behaviour this module relies on was exercised live on
2026-09-25 (captures under ``tests/fixtures/waterdata_ogc/``):

- ``monitoring_location_id`` needs the ``USGS-`` prefix; without it the API
  answers 200 with zero features, which would read as "no data".
- ``time`` is an RFC 3339 interval, **inclusive at both ends**, and may span at
  most :data:`MAX_TIME_ENVELOPE_DAYS`; a longer one is a hard 400. Omitting it
  returns only the most recent year, so a window is always sent.
- Pages hold at most :data:`PAGE_LIMIT` rows; ``links[rel=next]`` carries a
  cursor to the next page.
- An empty window is a 200 with zero features, not the legacy service's 400.
- Timestamps are UTC only. Local time is derived here from the monitoring
  location's ``time_zone_abbreviation`` and ``uses_daylight_savings``.
- A site measuring one parameter with several sensors returns them interleaved
  in long format at identical timestamps, told apart only by
  ``time_series_id``. Both of 03612600's 00065 series are flagged primary, so
  "primary" cannot be used to choose between them.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import timedelta, timezone
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import pandas as pd
import requests

from flowfreq.usgs import (
    _IV_COLUMNS,
    IV_PARAMETERS,
    NoInstantaneousDataError,
    _chunk_date_range,
    _empty_iv_frame,
    _iv_value_column,
    check_ts_id_form,
)

logger = logging.getLogger(__name__)

#: Root of the Water Data OGC API collections.
WATERDATA_BASE_URL = "https://api.waterdata.usgs.gov/ogcapi/v1/collections"

#: Longest ``time`` interval the ``continuous`` collection accepts, in days.
MAX_TIME_ENVELOPE_DAYS = 1100

#: Rows requested per page; the service's own cap.
PAGE_LIMIT = 50_000

#: Statistic code of an instantaneous (unit-value) series.
INSTANTANEOUS_STATISTIC = "00011"

#: IANA zones for the standard-time abbreviations the monitoring-locations
#: collection reports, used when the site observes daylight saving time.
DST_ZONES: Dict[str, str] = {
    "EST": "America/New_York",
    "CST": "America/Chicago",
    "MST": "America/Denver",
    "PST": "America/Los_Angeles",
    "AKST": "America/Anchorage",
    "HST": "Pacific/Honolulu",
    "AST": "America/Puerto_Rico",
}

#: UTC offsets in hours for standard-time abbreviations, used when the site
#: does *not* observe daylight saving time (e.g. MST without DST is Arizona).
STANDARD_OFFSETS: Dict[str, int] = {
    "UTC": 0,
    "GMT": 0,
    "AST": -4,
    "EST": -5,
    "CST": -6,
    "MST": -7,
    "PST": -8,
    "AKST": -9,
    "HST": -10,
    "SST": -11,
    "ChST": 10,
}

#: Approval status to the legacy NWIS approval code.
APPROVAL_CODES: Dict[str, str] = {"Approved": "A", "Provisional": "P"}

#: Qualifier tokens with a legacy NWIS code. Anything else is kept verbatim.
QUALIFIER_CODES: Dict[str, str] = {"ESTIMATED": "e"}

#: A local-time rule: an IANA zone name, or a fixed UTC offset plus its label.
LocalZone = Union[str, Tuple[timezone, str]]


class AmbiguousTimeSeriesError(ValueError):
    """Raised when a site reports one parameter from several sensors.

    The Water Data API returns every series for the parameter interleaved at
    identical timestamps, and more than one can be flagged primary. Picking one
    would hand back a plausible series from possibly the wrong instrument, and
    merging them would be worse, so the caller must choose with ``ts_id``. A
    ``ValueError`` subclass, as the legacy path's multi-sensor error is a plain
    ``ValueError``.
    """


@dataclass(frozen=True)
class TimeSeriesInfo:
    """One instantaneous series from the ``time-series-metadata`` collection."""

    time_series_id: str
    parameter_code: str
    sublocation_identifier: Optional[str]
    web_description: Optional[str]
    primary: Optional[str]
    begin: Optional[pd.Timestamp]
    end: Optional[pd.Timestamp]

    def describe(self) -> str:
        """One-line description used in error messages."""
        return (
            f"time_series_id={self.time_series_id} "
            f"sublocation={self.sublocation_identifier!r} "
            f"description={self.web_description!r} primary={self.primary!r} "
            f"begin={_fmt_ts(self.begin)} end={_fmt_ts(self.end)}"
        )


def _fmt_ts(ts: Optional[pd.Timestamp]) -> str:
    return "unknown" if ts is None else ts.strftime("%Y-%m-%dT%H:%M:%SZ")


def _parse_optional_ts(value: Any) -> Optional[pd.Timestamp]:
    if value is None or value == "":
        return None
    ts = pd.Timestamp(value)
    if pd.isna(ts):
        return None
    return ts.tz_convert("UTC") if ts.tzinfo is not None else ts.tz_localize("UTC")


def _location_id(site_no: str) -> str:
    """The ``USGS-`` prefixed id the API requires; a bare number matches nothing."""
    return f"USGS-{site_no}"


def _rfc3339(ts: pd.Timestamp) -> str:
    return ts.tz_convert("UTC").strftime("%Y-%m-%dT%H:%M:%SZ")


# ----------------------------------------------------------------------------
# HTTP
# ----------------------------------------------------------------------------


def _get_json(
    url: str, params: Optional[Mapping[str, Any]], timeout: int, what: str
) -> Dict[str, Any]:
    """GET a JSON document, naming the request in any failure."""
    try:
        response = requests.get(url, params=params, timeout=timeout)
        response.raise_for_status()
        payload = response.json()
    except (requests.RequestException, ValueError) as exc:
        raise requests.RequestException(f"Water Data API request failed ({what}): {exc}") from exc
    if not isinstance(payload, dict):
        raise requests.RequestException(
            f"Water Data API returned a non-object JSON document ({what})"
        )
    return payload


def _next_link(payload: Mapping[str, Any]) -> Optional[str]:
    for link in payload.get("links") or []:
        if link.get("rel") == "next" and link.get("href"):
            return str(link["href"])
    return None


def _get_all_features(
    url: str, params: Mapping[str, Any], timeout: int, what: str
) -> List[Dict[str, Any]]:
    """Fetch every page of an items request by following ``links[rel=next]``.

    The next link carries the full query plus a cursor, so it is requested
    as-is. A next link that repeats one already fetched raises rather than
    looping forever.
    """
    features: List[Dict[str, Any]] = []
    payload = _get_json(url, params, timeout, what)
    seen = set()
    page = 1
    while True:
        features.extend(payload.get("features") or [])
        nxt = _next_link(payload)
        if nxt is None:
            return features
        if nxt in seen:
            raise requests.RequestException(
                f"Water Data API paging loop ({what}): next link repeated after page {page}"
            )
        seen.add(nxt)
        page += 1
        payload = _get_json(nxt, None, timeout, f"{what}, page {page}")


# ----------------------------------------------------------------------------
# Monitoring location -> local time
# ----------------------------------------------------------------------------


def fetch_monitoring_location(site_no: str, timeout: int = 60) -> Dict[str, Any]:
    """Return the ``properties`` of a site's monitoring-locations record.

    Parameters
    ----------
    site_no : str
        USGS site number, without the ``USGS-`` prefix.
    timeout : int
        Request timeout in seconds.

    Returns
    -------
    dict
        The record's properties (``time_zone_abbreviation``,
        ``uses_daylight_savings``, ``drainage_area``, ...).

    Raises
    ------
    requests.RequestException
        The request failed.
    """
    payload = _get_json(
        f"{WATERDATA_BASE_URL}/monitoring-locations/items/{_location_id(site_no)}",
        {"f": "json"},
        timeout,
        f"monitoring location {site_no}",
    )
    return dict(payload.get("properties") or {})


def resolve_local_zone(tz_abbreviation: Optional[str], uses_dst: Optional[str]) -> LocalZone:
    """Turn a monitoring location's time-zone fields into a local-time rule.

    Parameters
    ----------
    tz_abbreviation : str
        ``time_zone_abbreviation`` -- always the *standard*-time abbreviation
        (``"CST"``, never ``"CDT"``).
    uses_dst : str
        ``uses_daylight_savings``, ``"Y"`` or ``"N"``.

    Returns
    -------
    str or (datetime.timezone, str)
        An IANA zone name when the site observes DST, so the abbreviation in
        effect is derived per instant; otherwise a fixed offset and the
        abbreviation to report for every row.

    Raises
    ------
    ValueError
        The abbreviation or DST flag is missing or not one this module knows.
        Local time is never guessed: a wrong zone shifts every local-day
        statistic in :mod:`flowfreq.subdaily`.
    """
    abbr = (tz_abbreviation or "").strip()
    flag = (uses_dst or "").strip().upper()
    if flag == "Y":
        if abbr not in DST_ZONES:
            raise ValueError(
                f"Unrecognized time_zone_abbreviation {tz_abbreviation!r} for a site that "
                f"observes daylight saving time; known: {sorted(DST_ZONES)}. Add it to "
                f"flowfreq.waterdata.DST_ZONES rather than guessing a zone."
            )
        return DST_ZONES[abbr]
    if flag == "N":
        if abbr not in STANDARD_OFFSETS:
            raise ValueError(
                f"Unrecognized time_zone_abbreviation {tz_abbreviation!r}; known: "
                f"{sorted(STANDARD_OFFSETS)}. Add it to "
                f"flowfreq.waterdata.STANDARD_OFFSETS rather than guessing an offset."
            )
        return (timezone(timedelta(hours=STANDARD_OFFSETS[abbr])), abbr)
    raise ValueError(
        f"uses_daylight_savings is {uses_dst!r}, expected 'Y' or 'N'; cannot derive "
        f"local time for this site"
    )


def local_time_columns(
    index: pd.DatetimeIndex, zone: LocalZone
) -> Tuple[pd.DatetimeIndex, List[str]]:
    """Naive local wall-clock times and zone codes for a UTC index.

    Parameters
    ----------
    index : pandas.DatetimeIndex
        tz-aware UTC timestamps.
    zone : str or (datetime.timezone, str)
        From :func:`resolve_local_zone`.

    Returns
    -------
    (pandas.DatetimeIndex, list of str)
        Naive local times (repeating across an autumn fall-back, exactly as the
        legacy service reports them) and the abbreviation in effect at each
        instant (``CST``/``CDT``).
    """
    if isinstance(zone, str):
        local = index.tz_convert(zone)
        codes = [str(c) for c in local.strftime("%Z")]
    else:
        fixed, label = zone
        local = index.tz_convert(fixed)
        codes = [label] * len(index)
    return local.tz_localize(None), codes


def _local_zone_arg(zone: LocalZone) -> Union[str, timezone]:
    return zone if isinstance(zone, str) else zone[0]


# ----------------------------------------------------------------------------
# Series discovery
# ----------------------------------------------------------------------------


def list_instantaneous_series(
    site_no: str, param_cd: str, timeout: int = 60
) -> List[TimeSeriesInfo]:
    """Instantaneous series a site reports for one parameter.

    Queries ``time-series-metadata`` filtered to the parameter, then keeps only
    statistic 00011 (instantaneous) -- the same parameter also carries daily
    means and annual-maximum series.

    Parameters
    ----------
    site_no : str
        USGS site number.
    param_cd : str
        ``"00060"`` or ``"00065"``.
    timeout : int
        Request timeout in seconds.

    Returns
    -------
    list of TimeSeriesInfo
        Possibly empty; sorted by ``time_series_id`` for a stable error message.
    """
    features = _get_all_features(
        f"{WATERDATA_BASE_URL}/time-series-metadata/items",
        {
            "f": "json",
            "monitoring_location_id": _location_id(site_no),
            "parameter_code": param_cd,
            "skipGeometry": "true",
            "limit": 1000,
        },
        timeout,
        f"time-series metadata for {site_no} {param_cd}",
    )
    series: List[TimeSeriesInfo] = []
    for feature in features:
        props = feature.get("properties") or {}
        if props.get("parameter_code") != param_cd:
            continue
        if props.get("statistic_id") != INSTANTANEOUS_STATISTIC:
            continue
        series.append(
            TimeSeriesInfo(
                time_series_id=str(props.get("id") or feature.get("id")),
                parameter_code=param_cd,
                sublocation_identifier=props.get("sublocation_identifier"),
                web_description=props.get("web_description"),
                primary=props.get("primary"),
                begin=_parse_optional_ts(props.get("begin")),
                end=_parse_optional_ts(props.get("end")),
            )
        )
    return sorted(series, key=lambda s: s.time_series_id)


def select_series(
    series: Sequence[TimeSeriesInfo],
    ts_id: Optional[str],
    site_no: str,
    param_cd: str,
    window: Optional[Tuple[pd.Timestamp, pd.Timestamp]] = None,
) -> TimeSeriesInfo:
    """Choose the one series to download, refusing when that is ambiguous.

    Parameters
    ----------
    series : sequence of TimeSeriesInfo
        From :func:`list_instantaneous_series`.
    ts_id : str, optional
        A ``time_series_id`` (32-hex UUID) chosen by the caller.
    site_no, param_cd : str
        For messages.
    window : (Timestamp, Timestamp), optional
        UTC request window. A series whose known ``begin``/``end`` lies wholly
        outside it cannot contribute rows and is not counted as a competitor,
        so a sensor retired years ago does not force every later request to
        name a ``ts_id``. Series with unknown bounds always count.

    Returns
    -------
    TimeSeriesInfo

    Raises
    ------
    NoInstantaneousDataError
        No instantaneous series for this parameter.
    ValueError
        ``ts_id`` names no instantaneous series for this parameter.
    AmbiguousTimeSeriesError
        Several series and no ``ts_id``.
    """
    description = IV_PARAMETERS[param_cd][1]
    if ts_id is not None:
        wanted = ts_id.strip().lower().replace("-", "")
        for info in series:
            if info.time_series_id.lower().replace("-", "") == wanted:
                return info
        available = "; ".join(s.describe() for s in series) or "none"
        raise ValueError(
            f"ts_id={ts_id!r} does not match any instantaneous {description} "
            f"({param_cd}) series at site {site_no}; available: {available}. On the "
            f"waterdata-ogc backend ts_id is the 32-hex time_series_id, not the "
            f"legacy NWIS DD number."
        )
    if not series:
        raise NoInstantaneousDataError(
            f"The Water Data API lists no instantaneous {description} ({param_cd}) "
            f"series for site {site_no}"
        )
    candidates = list(series)
    if window is not None and len(candidates) > 1:
        lo, hi = window
        candidates = [
            s
            for s in candidates
            if not ((s.end is not None and s.end < lo) or (s.begin is not None and s.begin > hi))
        ] or candidates
    if len(candidates) > 1:
        listing = "\n  ".join(s.describe() for s in candidates)
        raise AmbiguousTimeSeriesError(
            f"Site {site_no} reports {len(candidates)} separate instantaneous "
            f"{description} ({param_cd}) series:\n  {listing}\nPass ts_id=<time_series_id> "
            f"to choose one. Series are never merged, and 'primary' does not "
            f"disambiguate (several can be flagged primary)."
        )
    return candidates[0]


# ----------------------------------------------------------------------------
# Windows and chunks
# ----------------------------------------------------------------------------


def local_day_bounds(
    start_date: str, end_date: str, zone: LocalZone
) -> Tuple[pd.Timestamp, pd.Timestamp]:
    """UTC ``[start, end)`` covering local calendar days ``start_date``..``end_date``.

    The legacy service reads ``startDT``/``endDT`` as dates local to the gage,
    so the same window is reproduced here: from local midnight starting
    `start_date` to local midnight ending `end_date`, both converted to UTC.
    """
    tz = _local_zone_arg(zone)
    start = pd.Timestamp(start_date)
    end = pd.Timestamp(end_date)
    if pd.isna(start) or pd.isna(end):
        raise ValueError(f"Could not parse date range {start_date!r} to {end_date!r}")
    if start > end:
        raise ValueError(f"start_date {start_date} is after end_date {end_date}")
    lo = start.normalize().tz_localize(tz, nonexistent="shift_forward", ambiguous=True)
    hi = (end.normalize() + pd.Timedelta(days=1)).tz_localize(
        tz, nonexistent="shift_forward", ambiguous=True
    )
    return lo.tz_convert("UTC"), hi.tz_convert("UTC")


def plan_chunks(
    start_date: str, end_date: str, chunk_years: int, zone: LocalZone
) -> List[Tuple[pd.Timestamp, pd.Timestamp]]:
    """Half-open UTC intervals ``[lo, hi)`` for each request chunk.

    Chunks come from the legacy :func:`flowfreq.usgs._chunk_date_range`, so the
    two backends split a window identically; each chunk's ``hi`` is the next
    chunk's ``lo``, so together they cover the window with no gap and no
    overlap.

    Raises
    ------
    ValueError
        ``chunk_years`` < 1, or any chunk would exceed
        :data:`MAX_TIME_ENVELOPE_DAYS` (``chunk_years`` > 3). Checked before any
        request is sent, since the service rejects such a window outright.
    """
    if chunk_years < 1:
        raise ValueError(f"chunk_years must be >= 1, got {chunk_years}")
    chunks = [
        local_day_bounds(lo, hi, zone)
        for lo, hi in _chunk_date_range(start_date, end_date, chunk_years)
    ]
    limit = pd.Timedelta(days=MAX_TIME_ENVELOPE_DAYS)
    for lo, hi in chunks:
        if hi - lo > limit:
            raise ValueError(
                f"chunk_years={chunk_years} gives a request window of "
                f"{(hi - lo) / pd.Timedelta(days=1):.2f} days ({_rfc3339(lo)} to "
                f"{_rfc3339(hi)}), over the Water Data API's {MAX_TIME_ENVELOPE_DAYS}-day "
                f"limit. Use chunk_years <= 3."
            )
    return chunks


# ----------------------------------------------------------------------------
# Parsing
# ----------------------------------------------------------------------------


def map_qualification_code(approval_status: Optional[str], qualifiers: Any) -> str:
    """Legacy-style NWIS qualification code for one row.

    Approval first (``A``/``P``), then each qualifier token, joined with
    ``:`` as the legacy RDB ``_cd`` column does (``A:e``, ``P:Ice``). Mapped
    tokens take their legacy code (``ESTIMATED`` -> ``e``); unmapped tokens
    such as ``ICE`` are kept verbatim rather than dropped.

    Parameters
    ----------
    approval_status : str, optional
        ``"Approved"`` or ``"Provisional"``; an unknown status is kept verbatim.
    qualifiers : list of str, str or None
        The row's ``qualifier`` field.

    Returns
    -------
    str
    """
    parts: List[str] = []
    if approval_status:
        code = APPROVAL_CODES.get(approval_status)
        if code is None:
            logger.debug("Unmapped approval_status %r kept verbatim", approval_status)
            code = str(approval_status)
        parts.append(code)
    tokens: Iterable[Any]
    if qualifiers is None:
        tokens = ()
    elif isinstance(qualifiers, str):
        tokens = (qualifiers,)
    else:
        tokens = qualifiers
    for token in tokens:
        if token is None or str(token) == "":
            continue
        mapped = QUALIFIER_CODES.get(str(token))
        if mapped is None:
            logger.debug("Unmapped qualifier token %r kept verbatim", token)
            mapped = str(token)
        parts.append(mapped)
    return ":".join(parts)


def parse_continuous_features(
    features: Sequence[Mapping[str, Any]],
    param_cd: str,
    time_series_id: str,
) -> pd.DataFrame:
    """Parse ``continuous`` features for one series into a UTC-indexed frame.

    Parameters
    ----------
    features : sequence of dict
        GeoJSON features from the ``continuous`` collection.
    param_cd : str
        Parameter code; names the value column.
    time_series_id : str
        The series requested. Every row must carry it -- a row from another
        series means the response is not what was asked for, and it is never
        merged in.

    Returns
    -------
    pandas.DataFrame
        UTC index named ``datetime``; value column plus ``qualification_code``
        (local-time columns are added by the caller, which knows the zone).
        Rows with a null value (e.g. ice-affected, no estimate) are dropped, as
        the legacy parser drops non-numeric values.

    Raises
    ------
    ValueError
        A row from a different series or parameter, or conflicting duplicate
        timestamps within the series.
    """
    value_col = _iv_value_column(param_cd)
    if not features:
        frame = _empty_iv_frame(value_col)
        return frame[[value_col, "qualification_code"]]

    times: List[Any] = []
    values: List[Any] = []
    codes: List[str] = []
    for feature in features:
        props = feature.get("properties") or {}
        row_ts = props.get("time_series_id")
        if row_ts != time_series_id:
            raise ValueError(
                f"Water Data API returned a row from time series {row_ts!r} while "
                f"{time_series_id!r} was requested; refusing to merge series"
            )
        row_pc = props.get("parameter_code")
        if row_pc is not None and row_pc != param_cd:
            raise ValueError(f"Water Data API returned parameter {row_pc!r}, expected {param_cd}")
        times.append(props.get("time"))
        values.append(props.get("value"))
        codes.append(map_qualification_code(props.get("approval_status"), props.get("qualifier")))

    index = pd.DatetimeIndex(pd.to_datetime(times, utc=True), name="datetime")
    frame = pd.DataFrame(
        {
            value_col: pd.to_numeric(pd.Series(values, dtype=object), errors="coerce")
            .astype(float)
            .to_numpy(),
            "qualification_code": codes,
        },
        index=index,
    )
    frame = frame[frame.index.notna()]

    null_values = frame[value_col].isna()
    if null_values.any():
        logger.debug(
            "Dropping %d row(s) of series %s with no value (qualifiers: %s)",
            int(null_values.sum()),
            time_series_id,
            sorted(set(frame.loc[null_values, "qualification_code"])),
        )
        frame = frame[~null_values]

    return _drop_identical_duplicates(frame.sort_index(kind="stable"), time_series_id)


def _drop_identical_duplicates(frame: pd.DataFrame, time_series_id: str) -> pd.DataFrame:
    """Collapse exact duplicate rows of one series; raise on conflicting ones.

    A repeated timestamp *with the same value and code* carries no information
    (a re-sent row, e.g. across a page boundary) and is dropped with a warning.
    A repeated timestamp with *different* values means the source disagrees
    with itself, and keeping either would be an arbitrary choice presented as
    data -- the same reason several sensors are refused rather than merged.
    """
    dup = frame.index.duplicated(keep=False)
    if not dup.any():
        return frame
    block = frame[dup].reset_index()
    conflicting = block.groupby("datetime").nunique(dropna=False).gt(1).any(axis=1)
    if conflicting.any():
        stamps = [ts.isoformat() for ts in conflicting[conflicting].index[:5]]
        raise ValueError(
            f"Time series {time_series_id} has {int(conflicting.sum())} timestamp(s) "
            f"with conflicting values (e.g. {stamps}); refusing to pick one"
        )
    n_dropped = int(frame.index.duplicated(keep="first").sum())
    logger.warning(
        "Dropped %d exact duplicate row(s) from time series %s", n_dropped, time_series_id
    )
    return frame[~frame.index.duplicated(keep="first")]


# ----------------------------------------------------------------------------
# Retrieval
# ----------------------------------------------------------------------------


def fetch_continuous_chunk(
    site_no: str,
    param_cd: str,
    time_series_id: str,
    lo: pd.Timestamp,
    hi: pd.Timestamp,
    timeout: int = 60,
) -> pd.DataFrame:
    """All rows of one series in the half-open UTC window ``[lo, hi)``.

    The service's ``time`` interval is inclusive at both ends, so a row exactly
    at ``hi`` would be returned by this chunk *and* the next. It is requested
    and then discarded here, so each instant belongs to exactly one chunk. An
    empty window is an ordinary empty frame.

    Raises
    ------
    requests.RequestException
        A page failed; the message names the window.
    ValueError
        As :func:`parse_continuous_features`.
    """
    what = f"site {site_no}, parameter {param_cd}, window {_rfc3339(lo)}/{_rfc3339(hi)}"
    features = _get_all_features(
        f"{WATERDATA_BASE_URL}/continuous/items",
        {
            "f": "json",
            "monitoring_location_id": _location_id(site_no),
            "parameter_code": param_cd,
            "time_series_id": time_series_id,
            "time": f"{_rfc3339(lo)}/{_rfc3339(hi)}",
            "skipGeometry": "true",
            "limit": PAGE_LIMIT,
        },
        timeout,
        what,
    )
    frame = parse_continuous_features(features, param_cd, time_series_id)
    return frame[(frame.index >= lo) & (frame.index < hi)]


def download_instantaneous(
    site_no: str,
    param_cd: str,
    *,
    start_date: Optional[str] = None,
    end_date: Optional[str] = None,
    chunk_years: int = 1,
    ts_id: Optional[str] = None,
    timeout: int = 60,
) -> pd.DataFrame:
    """Instantaneous values for one parameter from the Water Data OGC API.

    The ``backend="waterdata-ogc"`` body of
    :meth:`flowfreq.usgs.USGSgage.download_instantaneous_flow` and
    :meth:`~flowfreq.usgs.USGSgage.download_instantaneous_stage`.

    Parameters
    ----------
    site_no : str
        USGS site number, without the ``USGS-`` prefix.
    param_cd : str
        ``"00060"`` (discharge) or ``"00065"`` (gage height).
    start_date, end_date : str, optional
        ``YYYY-MM-DD`` local calendar days, inclusive. Default to the chosen
        series' ``begin`` and ``end`` from ``time-series-metadata`` (its real
        period of record); ``end`` falls back to today.
    chunk_years : int
        Years per request window; must be 1-3 (the service caps a window at
        :data:`MAX_TIME_ENVELOPE_DAYS`). Each window is also paged.
    ts_id : str, optional
        ``time_series_id`` (32-hex UUID) selecting one series where the site
        has several. Required in that case. A legacy NWIS DD number is refused
        before any request (:func:`flowfreq.usgs.check_ts_id_form`).
    timeout : int
        Per-request timeout in seconds.

    Returns
    -------
    pandas.DataFrame
        The legacy instantaneous frame: UTC index named ``datetime``, sorted
        and unique; columns ``<value>``, ``datetime_local``, ``tz_cd``,
        ``qualification_code``.

    Raises
    ------
    NoInstantaneousDataError
        No instantaneous series, or no values in the window -- matching the
        legacy path, so callers need not know which backend ran.
    AmbiguousTimeSeriesError
        Several series and no ``ts_id``.
    ValueError
        Bad ``ts_id``, ``chunk_years`` out of range, an unmappable time zone,
        conflicting duplicate timestamps, or rows from an unrequested series.
    requests.RequestException
        Any request failed; no partial frame is returned.
    """
    value_col = _iv_value_column(param_cd)
    description = IV_PARAMETERS[param_cd][1]
    check_ts_id_form(ts_id, "waterdata-ogc", param_cd)
    if start_date is not None and end_date is not None:
        # Fail on an oversized chunk before any request; the check is repeated
        # below on the real local-day bounds.
        plan_chunks(start_date, end_date, chunk_years, (timezone.utc, "UTC"))
    elif chunk_years < 1:
        raise ValueError(f"chunk_years must be >= 1, got {chunk_years}")

    location = fetch_monitoring_location(site_no, timeout)
    zone = resolve_local_zone(
        location.get("time_zone_abbreviation"), location.get("uses_daylight_savings")
    )

    series = list_instantaneous_series(site_no, param_cd, timeout)
    window: Optional[Tuple[pd.Timestamp, pd.Timestamp]] = None
    if start_date is not None and end_date is not None:
        window = local_day_bounds(start_date, end_date, zone)
    chosen = select_series(series, ts_id, site_no, param_cd, window)

    local_tz = _local_zone_arg(zone)
    if start_date is None:
        if chosen.begin is None:
            raise NoInstantaneousDataError(
                f"Series {chosen.time_series_id} at site {site_no} has no recorded begin "
                f"date; pass start_date explicitly"
            )
        start_date = chosen.begin.tz_convert(local_tz).strftime("%Y-%m-%d")
    if end_date is None:
        end_ts = chosen.end if chosen.end is not None else pd.Timestamp.now(tz="UTC")
        end_date = end_ts.tz_convert(local_tz).strftime("%Y-%m-%d")

    chunks = plan_chunks(start_date, end_date, chunk_years, zone)
    frames = []
    for lo, hi in chunks:
        part = fetch_continuous_chunk(site_no, param_cd, chosen.time_series_id, lo, hi, timeout)
        if not part.empty:
            frames.append(part)

    if not frames:
        raise NoInstantaneousDataError(
            f"No instantaneous {description} (parameter {param_cd}) found for site "
            f"{site_no} between {start_date} and {end_date} in series "
            f"{chosen.time_series_id} (record {_fmt_ts(chosen.begin)} to "
            f"{_fmt_ts(chosen.end)})"
        )

    combined = pd.concat(frames).sort_index(kind="stable")
    combined = _drop_identical_duplicates(combined, chosen.time_series_id)

    index = pd.DatetimeIndex(combined.index, name="datetime")
    local, codes = local_time_columns(index, zone)
    out = pd.DataFrame(
        {
            value_col: combined[value_col].astype(float).to_numpy(),
            "datetime_local": local.to_numpy(),
            "tz_cd": pd.Series(codes, dtype=object).to_numpy(),
            "qualification_code": combined["qualification_code"].astype(object).to_numpy(),
        },
        index=index,
    )
    return out[[value_col, *_IV_COLUMNS[1:]]]
