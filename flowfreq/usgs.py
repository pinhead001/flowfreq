"""
flowfreq.usgs - USGS data retrieval
"""

from __future__ import annotations

import logging
import math
import re
from datetime import datetime, timezone
from functools import cached_property
from io import StringIO
from pathlib import Path
from typing import ClassVar, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import requests

from flowfreq.peak_sources import DEFAULT_BACKEND

logger = logging.getLogger(__name__)

NWIS_TZ_OFFSETS: Dict[str, int] = {
    "UTC": 0,
    "GMT": 0,
    "AST": -4,
    "ADT": -3,
    "EST": -5,
    "EDT": -4,
    "CST": -6,
    "CDT": -5,
    "MST": -7,
    "MDT": -6,
    "PST": -8,
    "PDT": -7,
    "AKST": -9,
    "AKDT": -8,
    "HST": -10,
    "HDT": -9,
    "SST": -11,
    "ChST": 10,
}
"""UTC offsets, in hours, for the time-zone abbreviations NWIS reports in the
``tz_cd`` column of instantaneous-value output.

NWIS returns wall-clock time local to the gage, so a single site's record mixes
standard and daylight codes across daylight-saving transitions. These offsets
are what let the series be expressed on a single monotonic UTC axis without
discarding the local time actually reported.
"""


class NoInstantaneousDataError(ValueError):
    """Raised when a site has no instantaneous-value record for the request.

    Distinct from a transport failure: the request succeeded, and NWIS has no
    unit-value discharge data for this site (or none within the requested
    window). Instantaneous records generally begin around 2007 and are often
    far shorter than the same gage's daily record.
    """


class GageAttributes:
    """Load and manage gage attributes from a local CSV file.

    The CSV file should have columns:
    - site_no: USGS site number (8-digit string)
    - site_name: Station name
    - drainage_area_sqmi: Drainage area in square miles
    - state: State abbreviation (optional)
    - huc8: HUC-8 watershed code (optional)
    """

    _instance: ClassVar[Optional["GageAttributes"]] = None
    _data: ClassVar[Optional[pd.DataFrame]] = None

    @classmethod
    def _find_data_file(cls) -> Optional[Path]:
        """Find the gage_attributes.csv file in various locations."""
        # Packaged location first: it is the only one that resolves for an
        # installed user. The cwd fallback lets a caller shadow the shipped
        # table with a local one without reinstalling.
        candidates = [
            Path(__file__).parent / "data" / "gage_attributes.csv",
            Path.cwd() / "data" / "gage_attributes.csv",
        ]

        for path in candidates:
            if path.exists():
                return path
        return None

    def __new__(cls, path: Optional[Path] = None):
        """Singleton pattern - only load the file once."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            if path:
                cls._load_data(path)
            else:
                data_file = cls._find_data_file()
                if data_file:
                    cls._load_data(data_file)
                else:
                    cls._data = pd.DataFrame()
        return cls._instance

    @classmethod
    def _load_data(cls, path: Path) -> None:
        """Load gage attributes from CSV file."""
        if path.exists():
            try:
                df = pd.read_csv(path, dtype={"site_no": str})
                # Ensure site_no is 8 digits with leading zeros
                df["site_no"] = df["site_no"].str.zfill(8)
                cls._data = df.set_index("site_no")
            except Exception:
                logger.warning("Failed to parse gage attributes file %s", path, exc_info=True)
                cls._data = pd.DataFrame()
        else:
            cls._data = pd.DataFrame()

    @classmethod
    def reload(cls, path: Optional[Path] = None) -> None:
        """Reload attributes from file (useful after file changes)."""
        if path:
            cls._load_data(path)
        else:
            data_file = cls._find_data_file()
            if data_file:
                cls._load_data(data_file)
            else:
                cls._data = pd.DataFrame()

    @classmethod
    def get_attributes(cls, site_no: str) -> Optional[Dict]:
        """Get attributes for a gage by site number.

        Returns dict with site_name, drainage_area_sqmi, etc. or None if not found.
        """
        if cls._data is None or cls._data.empty:
            cls()  # Initialize if needed

        site_no = str(site_no).zfill(8)
        if cls._data is not None and site_no in cls._data.index:
            row = cls._data.loc[site_no]
            return row.to_dict()
        return None

    @classmethod
    def get_drainage_area(cls, site_no: str) -> Optional[float]:
        """Get drainage area for a gage by site number."""
        attrs = cls.get_attributes(site_no)
        if attrs and "drainage_area_sqmi" in attrs:
            try:
                return float(attrs["drainage_area_sqmi"])
            except (ValueError, TypeError):
                return None
        return None

    @classmethod
    def get_site_name(cls, site_no: str) -> Optional[str]:
        """Get site name for a gage by site number."""
        attrs = cls.get_attributes(site_no)
        if attrs and "site_name" in attrs:
            return str(attrs["site_name"])
        return None

    @classmethod
    def status(cls) -> Dict:
        """Return status info about the loaded data file (for debugging)."""
        # Ensure singleton is initialized
        if cls._instance is None:
            cls()

        data_file = cls._find_data_file()
        return {
            "data_file": str(data_file) if data_file else None,
            "file_exists": data_file.exists() if data_file else False,
            "num_gages": len(cls._data) if cls._data is not None and not cls._data.empty else 0,
            "gages": list(cls._data.index) if cls._data is not None and not cls._data.empty else [],
        }


def _first_float(df: pd.DataFrame, column: str) -> Optional[float]:
    """First value of *column* as a float, or None if it is absent or blank.

    NWIS reports a missing numeric field as an empty column, which pandas reads
    as NaN. ``float(nan)`` succeeds, so a naive conversion silently stores NaN
    rather than None -- and NaN then propagates into summary tables and map
    coordinates without ever raising. Checking for it is the whole point of this
    helper.
    """
    if column not in df.columns or len(df) == 0:
        return None
    try:
        value = float(df[column].iloc[0])
    except (ValueError, TypeError):
        return None
    return None if math.isnan(value) else value


class USGSgage:
    """Class to handle USGS gage data retrieval and storage."""

    BASE_URL_DAILY: ClassVar[str] = "https://waterservices.usgs.gov/nwis/dv/"
    BASE_URL_IV: ClassVar[str] = "https://waterservices.usgs.gov/nwis/iv/"
    BASE_URL_PEAKS: ClassVar[str] = "https://nwis.waterdata.usgs.gov/nwis/peak"
    BASE_URL_SITE: ClassVar[str] = "https://waterservices.usgs.gov/nwis/site/"

    #: Start of the default daily-values window, preceding every USGS daily
    #: record. It exists so that a caller who supplies no range still gets one
    #: sent: NWIS reads a missing range as "the most recent day", not as "the
    #: whole record". See :meth:`download_daily_flow`.
    DEFAULT_START_DATE: ClassVar[str] = "1850-01-01"

    def __init__(self, site_no: str):
        self._site_no = str(site_no).zfill(8)
        self._site_name: Optional[str] = None
        self._drainage_area: Optional[float] = None
        self._daily_data: Optional[pd.DataFrame] = None
        self._peak_data: Optional[pd.DataFrame] = None
        self._daily_por_start: Optional[str] = None
        self._daily_por_end: Optional[str] = None
        self._latitude: Optional[float] = None
        self._longitude: Optional[float] = None
        self._instantaneous_data: Optional[pd.DataFrame] = None
        self._stage_data: Optional[pd.DataFrame] = None
        self._iv_por_start: Optional[str] = None
        self._iv_por_end: Optional[str] = None

    @property
    def site_no(self) -> str:
        return self._site_no

    @property
    def site_name(self) -> Optional[str]:
        return self._site_name

    @site_name.setter
    def site_name(self, value: str):
        self._site_name = value

    @property
    def drainage_area(self) -> Optional[float]:
        return self._drainage_area

    @drainage_area.setter
    def drainage_area(self, value: float):
        self._drainage_area = value

    @property
    def daily_data(self) -> Optional[pd.DataFrame]:
        return self._daily_data

    @daily_data.setter
    def daily_data(self, value: pd.DataFrame):
        self._daily_data = value

    @property
    def peak_data(self) -> Optional[pd.DataFrame]:
        return self._peak_data

    @peak_data.setter
    def peak_data(self, value: pd.DataFrame):
        self._peak_data = value

    @property
    def daily_por_start(self) -> Optional[str]:
        return self._daily_por_start

    @property
    def daily_por_end(self) -> Optional[str]:
        return self._daily_por_end

    @property
    def latitude(self) -> Optional[float]:
        """Decimal-degree latitude from the NWIS site service, or None."""
        return self._latitude

    @property
    def longitude(self) -> Optional[float]:
        """Decimal-degree longitude from the NWIS site service, or None.

        Negative in the western hemisphere, as NWIS reports it.
        """
        return self._longitude

    @property
    def instantaneous_data(self) -> Optional[pd.DataFrame]:
        """Most recently downloaded instantaneous (unit-value) *discharge* series.

        Gage height has its own :attr:`instantaneous_stage` rather than sharing
        this attribute: the two frames differ in their value column, so a stage
        download landing here would make ``gage.instantaneous_data["flow_cfs"]``
        raise KeyError at a call site that had no reason to expect it.
        """
        return self._instantaneous_data

    @instantaneous_data.setter
    def instantaneous_data(self, value: pd.DataFrame):
        self._instantaneous_data = value

    @property
    def instantaneous_stage(self) -> Optional[pd.DataFrame]:
        """Most recently downloaded instantaneous (unit-value) gage-height series."""
        return self._stage_data

    @instantaneous_stage.setter
    def instantaneous_stage(self, value: pd.DataFrame):
        self._stage_data = value

    @property
    def iv_por_start(self) -> Optional[str]:
        """First date of the instantaneous (unit-value) record, if known."""
        return self._iv_por_start

    @property
    def iv_por_end(self) -> Optional[str]:
        """Last date of the instantaneous (unit-value) record, if known."""
        return self._iv_por_end

    @cached_property
    def period_of_record(self) -> Optional[Tuple[int, int]]:
        if self._peak_data is not None:
            return (
                int(self._peak_data["water_year"].min()),
                int(self._peak_data["water_year"].max()),
            )
        return None

    def fetch_site_info(self, use_local_first: bool = True) -> None:
        """Fetch site information (name, drainage area, POR).

        Parameters
        ----------
        use_local_first : bool
            If True, check local gage_attributes.csv first for site name and
            drainage area before falling back to USGS API. Default True.
        """
        # First try to get attributes from local file
        if use_local_first:
            local_attrs = GageAttributes.get_attributes(self._site_no)
            if local_attrs:
                if "site_name" in local_attrs and pd.notna(local_attrs["site_name"]):
                    self._site_name = str(local_attrs["site_name"])
                if "drainage_area_sqmi" in local_attrs and pd.notna(
                    local_attrs["drainage_area_sqmi"]
                ):
                    try:
                        self._drainage_area = float(local_attrs["drainage_area_sqmi"])
                    except (ValueError, TypeError):
                        pass

        # Fetch from USGS Site Service API for POR dates and any missing info
        self._fetch_from_usgs_site_service()

    def _fetch_from_usgs_site_service(self) -> None:
        """Fetch site info from USGS Site Service API.

        Makes two separate API calls because siteOutput=expanded and
        seriesCatalogOutput=true cannot be combined in a single request.
        """
        self._last_api_error: Optional[str] = None

        # Call 1: site metadata (name, drainage area, lat/lon) via siteOutput=expanded
        if self._site_name is None or self._drainage_area is None or self._latitude is None:
            params_site = {
                "format": "rdb",
                "sites": self._site_no,
                "siteOutput": "expanded",
            }

            try:
                response = requests.get(self.BASE_URL_SITE, params=params_site, timeout=30)
                response.raise_for_status()

                lines = response.text.split("\n")
                data_lines = [l for l in lines if not l.startswith("#") and l.strip()]

                if len(data_lines) >= 2:
                    df = pd.read_csv(StringIO("\n".join(data_lines)), sep="\t", skiprows=[1])

                    if self._site_name is None and "station_nm" in df.columns and len(df) > 0:
                        self._site_name = df["station_nm"].iloc[0]

                    for attr, column in (
                        ("_drainage_area", "drain_area_va"),
                        ("_latitude", "dec_lat_va"),
                        ("_longitude", "dec_long_va"),
                    ):
                        if getattr(self, attr) is None:
                            setattr(self, attr, _first_float(df, column))
            except Exception as e:
                logger.warning("Site metadata request failed for %s: %s", self._site_no, e)
                self._last_api_error = str(e)

        # Call 2: Get period of record with seriesCatalogOutput=true
        params_por = {
            "format": "rdb",
            "sites": self._site_no,
            "seriesCatalogOutput": "true",
            "parameterCd": "00060",  # Discharge
        }

        try:
            response = requests.get(self.BASE_URL_SITE, params=params_por, timeout=30)
            response.raise_for_status()

            lines = response.text.split("\n")
            data_lines = [l for l in lines if not l.startswith("#") and l.strip()]

            if len(data_lines) >= 2:
                df = pd.read_csv(StringIO("\n".join(data_lines)), sep="\t", skiprows=[1])

                # Get daily value POR (data_type_cd == 'dv' for daily values)
                if "data_type_cd" in df.columns:
                    dv_rows = df[df["data_type_cd"] == "dv"]
                    if len(dv_rows) > 0:
                        if "begin_date" in df.columns:
                            self._daily_por_start = str(dv_rows["begin_date"].iloc[0])
                        if "end_date" in df.columns:
                            self._daily_por_end = str(dv_rows["end_date"].iloc[0])

                    # Instantaneous (unit-value) POR, data_type_cd == 'uv'.
                    # Usually starts around 2007 and is much shorter than 'dv'.
                    uv_rows = df[df["data_type_cd"] == "uv"]
                    if len(uv_rows) > 0:
                        if "begin_date" in df.columns:
                            self._iv_por_start = str(uv_rows["begin_date"].iloc[0])
                        if "end_date" in df.columns:
                            self._iv_por_end = str(uv_rows["end_date"].iloc[0])
        except Exception as e:
            logger.warning("Period-of-record request failed for %s: %s", self._site_no, e)
            if self._last_api_error:
                self._last_api_error += f"; {str(e)}"
            else:
                self._last_api_error = str(e)

    def download_daily_flow(
        self,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
        *,
        timeout: int = 60,
        backend: str = DEFAULT_BACKEND,
    ) -> pd.DataFrame:
        """Download mean daily streamflow data from USGS.

        Parameters
        ----------
        start_date, end_date : str, optional
            ISO dates bounding the request. Both default to the full period of
            record. If both are given, start_date must not be after end_date.
        timeout : int
            Per-request timeout in seconds. Default 60 -- the full
            period-of-record request this method sends by default can return
            tens of thousands of RDB rows for a long-running, high-frequency
            site, well past what a tight timeout is sized for.
        backend : str, default :data:`flowfreq.peak_sources.DEFAULT_BACKEND`
            ``"waterdata-ogc"`` (the default since issue #29) reads the USGS
            Water Data OGC API ``daily`` collection
            (:func:`flowfreq.waterdata.download_daily`); ``"nwis-legacy"`` the
            NWIS daily-values service this method always used before, which
            USGS is retiring. The frame is the same shape either way, and a
            live parity test (``tests/test_daily_backend_parity.py``) checks
            the two agree day for day. See Notes on ice.

        Returns
        -------
        pd.DataFrame
            Daily mean flows indexed by ``date`` (naive, the gage's local
            calendar day), one column ``flow_cfs``. Days with no published
            value are absent, not NaN.

        Raises
        ------
        ValueError
            start_date or end_date is not a parseable date, or start_date is
            after end_date; an unknown ``backend``; or no daily data for the
            site. On ``"waterdata-ogc"``, also
            :class:`flowfreq.waterdata.AmbiguousTimeSeriesError` (a
            ``ValueError``) if the site has more than one daily-mean series --
            call :func:`flowfreq.waterdata.download_daily` with ``ts_id``.
        requests.RequestException
            The request failed.

        Notes
        -----
        **Ice.** An ice-affected day is either an estimate or nothing. Where
        USGS published an estimate, both backends return the number (OGC
        qualifiers ``["ESTIMATED", "ICE"]``, legacy code ``A:e``). Where it
        published none -- typically provisional winter data -- legacy writes
        the text ``Ice`` in the value column and the OGC API a null value
        with qualifier ``["ICE"]``; both backends drop the day. So neither
        backend turns an ice day into NaN in the returned frame, and on every
        record compared live the two return the same days. The difference is
        only in what is kept alongside: the qualifier, which
        :func:`flowfreq.waterdata.download_daily` returns as
        ``qualification_code`` (``A:e:ICE``) and which this method, like the
        legacy path, leaves out.

        ``flow_cfs`` is always float on ``"waterdata-ogc"``. The legacy parser
        infers the dtype from the RDB text, so a window whose values are all
        whole numbers comes back as int64 there.

        A date range is **always** sent, even when the caller supplies none.
        The NWIS daily-values service answers a range-less request by
        returning only the most recent day, and a one-row frame is not an
        error -- it is a degenerate record that builds a plausible, meaningless
        curve wherever it is used. Measured on 12449500: one row with no range,
        26,958 with one. See ``DEFAULT_START_DATE``.

        The default end date is computed in UTC, not local wall-clock time: a
        host whose clock is behind UTC would otherwise silently request (and
        receive) a narrower range than intended, with no error to signal it.
        """
        _validate_daily_flow_range(start_date, end_date)
        start = start_date or self.DEFAULT_START_DATE
        end = end_date or datetime.now(timezone.utc).date().isoformat()

        if backend == "nwis-legacy":
            df = self._download_daily_flow_legacy(start, end, timeout)
        elif backend == "waterdata-ogc":
            # Deferred: flowfreq.waterdata imports from this module.
            from flowfreq.waterdata import download_daily

            daily = download_daily(
                self._site_no, "00060", start_date=start, end_date=end, timeout=timeout
            )
            df = daily[["flow_cfs"]]
        else:
            raise ValueError(
                f"Unknown daily-value backend {backend!r}; expected 'waterdata-ogc' or "
                f"'nwis-legacy'"
            )

        self._daily_data = df
        return df

    def _download_daily_flow_legacy(self, start: str, end: str, timeout: int) -> pd.DataFrame:
        """Daily means from the legacy NWIS daily-values RDB service, unchanged.

        The ``backend="nwis-legacy"`` body of :meth:`download_daily_flow`, which
        has already resolved and validated the range.
        """
        params = {
            "format": "rdb",
            "sites": self._site_no,
            "parameterCd": "00060",
            "statCd": "00003",
            "startDT": start,
            "endDT": end,
        }

        response = requests.get(self.BASE_URL_DAILY, params=params, timeout=timeout)
        response.raise_for_status()

        lines = response.text.split("\n")
        data_lines = [l for l in lines if not l.startswith("#") and l.strip()]

        if len(data_lines) < 2:
            raise ValueError(f"No daily data found for site {self._site_no}")

        header_idx = 0
        for i, line in enumerate(data_lines):
            if "datetime" in line.lower():
                header_idx = i
                break

        df = pd.read_csv(StringIO("\n".join(data_lines[header_idx:])), sep="\t", skiprows=[1])

        for line in lines:
            if "#" in line and "TS id" in line:
                name_start = line.find(self._site_no) + len(self._site_no)
                self._site_name = line[name_start:].strip()
                break

        date_col = [c for c in df.columns if "datetime" in c.lower()][0]
        flow_col = [c for c in df.columns if "00060" in c and "cd" not in c.lower()]

        if not flow_col:
            raise ValueError("Flow data column not found")

        flow_col = flow_col[0]

        df["date"] = pd.to_datetime(df[date_col])
        df["flow_cfs"] = pd.to_numeric(df[flow_col], errors="coerce")
        df = df[["date", "flow_cfs"]].dropna()
        return df.set_index("date")

    def download_instantaneous_flow(  # pylint: disable=too-many-arguments
        self,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
        *,
        tz: Optional[str] = None,
        chunk_years: int = 1,
        ts_id: Optional[str] = None,
        timeout: int = 60,
        backend: str = DEFAULT_BACKEND,
    ) -> pd.DataFrame:
        """Download instantaneous (unit-value) streamflow data from USGS.

        Retrieves parameter 00060 (discharge, cfs), by default from the USGS
        Water Data OGC API ``continuous`` collection (:mod:`flowfreq.waterdata`)
        and optionally from the legacy NWIS instantaneous-values service. There
        is no statistic code: every recorded value is returned, typically at a
        15-minute interval.

        Parameters
        ----------
        start_date : str, optional
            First date to retrieve, ``YYYY-MM-DD``, a calendar day local to the
            gage. Defaults to the start of the instantaneous period of record:
            the chosen series' ``begin`` from ``time-series-metadata`` on
            ``"waterdata-ogc"``, or the NWIS site service's ``uv`` period on
            ``"nwis-legacy"``.
        end_date : str, optional
            Last date to retrieve, ``YYYY-MM-DD``, inclusive. Defaults to the
            end of the instantaneous period of record.
        tz : str, optional
            IANA time-zone name (e.g. ``"America/Los_Angeles"``) for the
            returned index. Default ``None`` leaves the index in UTC. See
            Notes on what is and is not converted.
        chunk_years : int
            Number of years per HTTP request. Default 1. See Notes on the
            practical request limit.
        ts_id : str, optional
            Time-series identifier, used to disambiguate sites that report
            discharge from more than one sensor. If a site has multiple 00060
            series and this is not given, the call raises rather than silently
            picking one. **Its form depends on** ``backend``: on the default
            ``"waterdata-ogc"`` it is the 32-hex ``time_series_id`` (list them
            with :func:`flowfreq.waterdata.list_instantaneous_series`); on
            ``"nwis-legacy"`` it is the NWIS DD number (``"60629"``). Passing
            one form to the other backend raises a ``ValueError`` saying so,
            before any request is sent, rather than failing to match.
        timeout : int
            Per-request timeout in seconds. Default 60.
        backend : str, default :data:`flowfreq.peak_sources.DEFAULT_BACKEND`
            ``"waterdata-ogc"`` (the default since issue #29) for the USGS
            Water Data OGC API (:mod:`flowfreq.waterdata`), or
            ``"nwis-legacy"`` for the NWIS instantaneous-values service, which
            USGS is retiring. The returned frame has the same columns and UTC
            index either way, and a live parity test
            (``tests/test_iv_backend_parity.py``) checks the two agree row for
            row. On ``"waterdata-ogc"``, ``ts_id`` is the 32-hex
            ``time_series_id``, ``chunk_years`` may be at most 3 (the API caps a
            window at 1100 days), the default start is the series' own period
            of record from ``time-series-metadata``, and
            ``datetime_local``/``tz_cd`` are derived from the monitoring
            location's time zone because the API reports UTC only. Its
            ``qualification_code`` can also carry qualifier tokens the legacy
            RDB leaves out on rows that have a value (``A:e:ICE``,
            ``A:EQUIP:e``); the approval code and legacy's own tokens are the
            same.

        Returns
        -------
        pd.DataFrame
            Indexed by timezone-aware datetime (UTC unless ``tz`` is given),
            sorted ascending with duplicate timestamps removed. Columns:

            - ``flow_cfs`` : float, discharge in cubic feet per second
            - ``datetime_local`` : the naive local wall-clock time as reported
              by NWIS, unmodified
            - ``tz_cd`` : the NWIS time-zone abbreviation for that record
              (e.g. ``PST``, ``PDT``)
            - ``qualification_code`` : NWIS data qualifier (``P`` provisional,
              ``A`` approved, ``e`` estimated), ``:``-joined; see ``backend``

        Raises
        ------
        NoInstantaneousDataError
            The site has no instantaneous discharge record, or none within the
            requested window. Instantaneous records generally begin around
            2007 and are often far shorter than the same gage's daily record,
            so a site with a century of daily values may have fifteen years of
            unit values or none at all.
        ValueError
            The response carried data rows but no discharge column, the site
            reports multiple 00060 series and ``ts_id`` was not supplied, or
            NWIS reported a time-zone abbreviation this module cannot map to a
            UTC offset.
        requests.RequestException
            A chunk request failed. The error names the failed window. The
            partial result is discarded rather than returned, so a network
            failure can never be mistaken for a short record.

        Notes
        -----
        **Practical request limit.** A 15-minute series is roughly 35,000
        values per year. NWIS becomes unreliable well before a decade of unit
        values in a single request, so this method splits the window into
        chunks of ``chunk_years`` (default one calendar year) and issues one
        request each. Chunks that legitimately contain no data are skipped;
        a chunk whose request *fails* raises, and no truncated frame is
        returned.

        **Time zone.** NWIS reports wall-clock time local to the gage together
        with a ``tz_cd`` abbreviation, so a single site's record mixes standard
        and daylight offsets across daylight-saving transitions. Leaving that
        as a naive index makes it non-monotonic — the autumn transition repeats
        an hour — which quietly corrupts any statistic computed per day. The
        index is therefore placed on a single UTC axis, and the local time and
        its zone code are preserved verbatim as columns; nothing is discarded.
        Pass ``tz`` to get the index in a named zone instead.

        Note that diel statistics must be grouped on *local* calendar days: at
        a Pacific gage, UTC days are offset seven to eight hours and cut across
        the daily cycle. The functions in :mod:`flowfreq.subdaily` take an
        explicit time zone for this reason.

        **Gage height** for the same site and window comes from
        :meth:`download_instantaneous_stage`, which is this method with
        parameter 00065; the two are separate calls returning separate frames,
        joinable on the index.

        **Storage.** For a series of this size Parquet is the format worth
        reaching for — it round-trips the tz-aware index and float dtypes
        exactly, where CSV loses both, and it is several times smaller. See
        :func:`flowfreq.flowio.save_flow_frame`.

        Examples
        --------
        >>> gage = USGSgage("12449950")
        >>> iv = gage.download_instantaneous_flow("2022-06-01", "2022-09-30")
        >>> iv.index.tz is not None
        True
        """
        combined = self._download_instantaneous_backend(
            "00060",
            backend=backend,
            start_date=start_date,
            end_date=end_date,
            tz=tz,
            chunk_years=chunk_years,
            ts_id=ts_id,
            timeout=timeout,
        )
        self._instantaneous_data = combined
        return combined

    def download_instantaneous_stage(  # pylint: disable=too-many-arguments
        self,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
        *,
        tz: Optional[str] = None,
        chunk_years: int = 1,
        ts_id: Optional[str] = None,
        timeout: int = 60,
        backend: str = DEFAULT_BACKEND,
    ) -> pd.DataFrame:
        """Download instantaneous (unit-value) gage height from USGS.

        Retrieves parameter 00065 (gage height, feet), by default from the
        Water Data OGC API. Every argument, the backend, the chunking, the UTC
        axis and the multi-sensor ``ts_id`` rule are exactly
        :meth:`download_instantaneous_flow`'s -- see that method's Notes, all
        of which apply here -- and the only difference in the returned frame is
        that the value column is ``gage_height_ft`` rather than ``flow_cfs``.

        Stage is the series most FERC and HCP ramping-rate conditions are
        actually written against, and the one that drives fish stranding, so it
        feeds :func:`flowfreq.subdaily.ramping_rates` with
        ``value_col="gage_height_ft"``.

        Parameters
        ----------
        start_date, end_date : str, optional
            ``YYYY-MM-DD`` window bounds, inclusive. Default to the site's
            instantaneous period of record.
        tz : str, optional
            IANA time zone for the returned index; default ``None`` leaves it
            in UTC.
        chunk_years : int
            Years per HTTP request. Default 1.
        ts_id : str, optional
            Time-series identifier, to disambiguate a site that reports gage
            height from more than one sensor: the 32-hex ``time_series_id`` on
            ``"waterdata-ogc"``, the NWIS DD number on ``"nwis-legacy"`` (see
            :meth:`download_instantaneous_flow`). Multiple 00065 series are
            *more* common than multiple 00060 ones -- a site can carry separate
            headwater and tailwater, or primary and backup, stage sensors -- so
            this is worth expecting rather than treating as exotic.
        timeout : int
            Per-request timeout in seconds. Default 60.
        backend : str, default :data:`flowfreq.peak_sources.DEFAULT_BACKEND`
            ``"waterdata-ogc"`` (default) or ``"nwis-legacy"``; see
            :meth:`download_instantaneous_flow`.

        Returns
        -------
        pd.DataFrame
            Indexed by tz-aware datetime (UTC unless `tz` is given), with
            columns ``gage_height_ft``, ``datetime_local``, ``tz_cd``,
            ``qualification_code``.

        Raises
        ------
        NoInstantaneousDataError
            The site has no instantaneous gage-height record, or none in the
            requested window. **This is a good deal more common than for
            discharge**: a site with fifteen years of unit-value discharge may
            carry a partial stage record or none at all, since stage is an
            intermediate measurement a gage is not obliged to publish. The
            message names the parameter, so this is not mistaken for the site
            having no unit values whatsoever.
        ValueError
            As :meth:`download_instantaneous_flow`.
        requests.RequestException
            A chunk request failed; the partial result is discarded.

        Notes
        -----
        **Gage height is not comparable between gages, and is not depth.** It
        is measured from an arbitrary local datum, so 4.0 ft at one gage and
        4.0 ft at another say nothing about which is deeper or which carries
        more water, and neither is the depth of the water. Only *changes*
        within one gage's record are meaningful, which is why
        :func:`flowfreq.subdaily.ramping_rates` refuses to compute a
        percent-per-hour rate on a stage series.

        **A datum can be reset.** When a gage is rebuilt or resurveyed its
        datum may change, which puts a step in this series that is not a
        hydrologic event and will read as an enormous ramping rate. Nothing
        here detects that; check the site's NWIS history before analyzing a
        long stage record across a station rebuild.

        Examples
        --------
        >>> gage = USGSgage("12449950")
        >>> stage = gage.download_instantaneous_stage("2022-06-01", "2022-09-30")
        >>> "gage_height_ft" in stage.columns
        True
        """
        combined = self._download_instantaneous_backend(
            "00065",
            backend=backend,
            start_date=start_date,
            end_date=end_date,
            tz=tz,
            chunk_years=chunk_years,
            ts_id=ts_id,
            timeout=timeout,
        )
        self._stage_data = combined
        return combined

    def _download_instantaneous_backend(  # pylint: disable=too-many-arguments
        self,
        param_cd: str,
        *,
        backend: str,
        start_date: Optional[str],
        end_date: Optional[str],
        tz: Optional[str],
        chunk_years: int,
        ts_id: Optional[str],
        timeout: int,
    ) -> pd.DataFrame:
        """Route an instantaneous download to the named backend.

        ``"nwis-legacy"`` goes to :meth:`_download_instantaneous`, unchanged;
        ``"waterdata-ogc"`` to :func:`flowfreq.waterdata.download_instantaneous`.
        A ``ts_id`` in the other backend's form is refused first.
        """
        if backend in IV_BACKENDS:
            check_ts_id_form(ts_id, backend, param_cd)
        if backend == "nwis-legacy":
            return self._download_instantaneous(
                param_cd,
                start_date=start_date,
                end_date=end_date,
                tz=tz,
                chunk_years=chunk_years,
                ts_id=ts_id,
                timeout=timeout,
            )
        if backend == "waterdata-ogc":
            # Deferred: flowfreq.waterdata imports from this module.
            from flowfreq.waterdata import download_instantaneous

            combined = download_instantaneous(
                self._site_no,
                param_cd,
                start_date=start_date,
                end_date=end_date,
                chunk_years=chunk_years,
                ts_id=ts_id,
                timeout=timeout,
            )
            if tz is not None:
                combined.index = combined.index.tz_convert(tz)
            return combined
        raise ValueError(
            f"Unknown instantaneous-value backend {backend!r}; expected "
            f"'nwis-legacy' or 'waterdata-ogc'"
        )

    def _download_instantaneous(  # pylint: disable=too-many-arguments
        self,
        param_cd: str,
        *,
        start_date: Optional[str],
        end_date: Optional[str],
        tz: Optional[str],
        chunk_years: int,
        ts_id: Optional[str],
        timeout: int,
    ) -> pd.DataFrame:
        """Retrieve one instantaneous parameter, chunked, on a UTC axis.

        The shared body of :meth:`download_instantaneous_flow` and
        :meth:`download_instantaneous_stage`. Parameterized rather than copied
        because the chunking, the no-data-is-an-HTTP-400 handling, the
        tz-abbreviation mapping and the duplicate-timestamp rule would
        otherwise exist twice and drift apart on the next fix to any of them.
        """
        value_col = _iv_value_column(param_cd)
        description = IV_PARAMETERS[param_cd][1]

        if chunk_years < 1:
            raise ValueError(f"chunk_years must be >= 1, got {chunk_years}")

        start_date, end_date = self._resolve_iv_window(start_date, end_date)
        chunks = _chunk_date_range(start_date, end_date, chunk_years)

        frames: List[pd.DataFrame] = []
        for chunk_start, chunk_end in chunks:
            chunk = self._request_iv_chunk(chunk_start, chunk_end, ts_id, timeout, param_cd)
            if not chunk.empty:
                frames.append(chunk)

        combined = pd.concat(frames) if frames else _empty_iv_frame(value_col)

        if combined.empty:
            raise NoInstantaneousDataError(
                f"No instantaneous {description} (parameter {param_cd}) found for site "
                f"{self._site_no} between {start_date} and {end_date}"
                + (
                    f"; the site's instantaneous record runs "
                    f"{self._iv_por_start} to {self._iv_por_end}"
                    if self._iv_por_start
                    else f"; NWIS lists no instantaneous {description} series for this site"
                )
            )

        combined = combined[~combined.index.duplicated(keep="first")].sort_index()

        if tz is not None:
            combined.index = combined.index.tz_convert(tz)

        return combined

    def _resolve_iv_window(
        self, start_date: Optional[str], end_date: Optional[str]
    ) -> Tuple[str, str]:
        """Fill in a missing start or end date from the instantaneous POR.

        Looks the period of record up from the NWIS site service only when it
        is actually needed, so that a fully specified window costs no extra
        request. Raises rather than guessing when the site has no unit-value
        series at all, which avoids issuing a request that cannot succeed.
        """
        if start_date is not None and end_date is not None:
            return start_date, end_date

        if self._iv_por_start is None and self._iv_por_end is None:
            self.fetch_site_info()

        if self._iv_por_start is None and start_date is None:
            raise NoInstantaneousDataError(
                f"NWIS lists no instantaneous (unit-value) discharge series for site "
                f"{self._site_no}. Instantaneous records generally begin around 2007; "
                f"pass explicit start_date and end_date to request anyway."
            )

        resolved_start = start_date if start_date is not None else self._iv_por_start
        resolved_end = end_date if end_date is not None else self._iv_por_end
        if resolved_end is None:
            resolved_end = pd.Timestamp.today().strftime("%Y-%m-%d")

        return str(resolved_start), str(resolved_end)

    def _request_iv_chunk(  # pylint: disable=too-many-arguments
        self,
        chunk_start: str,
        chunk_end: str,
        ts_id: Optional[str],
        timeout: int,
        param_cd: str = "00060",
    ) -> pd.DataFrame:
        """Request and parse one date-range chunk of instantaneous values.

        An empty window is a normal outcome — instantaneous records have gaps —
        and NWIS signals it with HTTP 400 rather than an empty body, so that
        case is translated to an empty frame. Any other failure propagates with
        the failed window named, so a partial record is never silently
        returned as if it were complete.
        """
        params = {
            "format": "rdb",
            "sites": self._site_no,
            "parameterCd": param_cd,
            "startDT": chunk_start,
            "endDT": chunk_end,
        }

        try:
            response = requests.get(self.BASE_URL_IV, params=params, timeout=timeout)
            if response.status_code == 400 and _is_no_data_response(response.text):
                return _empty_iv_frame(_iv_value_column(param_cd))
            response.raise_for_status()
        except requests.RequestException as exc:
            raise requests.RequestException(
                f"Instantaneous-value request failed for site {self._site_no}, "
                f"parameter {param_cd}, window {chunk_start} to {chunk_end}: {exc}"
            ) from exc

        return _parse_iv_rdb(response.text, ts_id=ts_id, param_cd=param_cd)

    def download_peak_flow(self, backend: str = DEFAULT_BACKEND) -> pd.DataFrame:
        """Download annual peak streamflow data from USGS.

        Parameters
        ----------
        backend : str, default :data:`flowfreq.peak_sources.DEFAULT_BACKEND`
            Which USGS service to read, by :mod:`flowfreq.peak_sources` name.
            The default, ``"waterdata-ogc"``, is the Water Data OGC API
            ``peaks`` collection. ``"nwis-legacy"`` is the NWIS peak RDB
            service this method always used before, which USGS is retiring.

        Returns
        -------
        pandas.DataFrame
            ``water_year``, ``peak_date``, ``peak_flow_cfs``,
            ``qualification_code``, whichever backend, also stored as
            :attr:`peak_data`. On ``"waterdata-ogc"``, ``peak_date`` is the
            UTC date when the time of day is known, so an evening peak can
            read one day later than the legacy local date; ``water_year`` is
            the API's own and matches legacy. Legacy's date-precision codes
            ``Bd``/``Bm`` are not in the API's ``qualification_code``; the
            placeholder date carries that information on both.

        Notes
        -----
        Site metadata is filled in as a side effect, as before. The legacy
        service parses :attr:`site_name` and :attr:`drainage_area` from the
        RDB header. The OGC backend reads them from the ``monitoring-locations``
        collection (:func:`flowfreq.waterdata.fetch_monitoring_location`),
        whose ``monitoring_location_name`` is upper-case (``"BIG SANDY RIVER
        AT BRUCETON, TN"``). If that one extra request fails, the peaks are
        still returned and the two attributes are left as they were, with a
        warning logged.

        Raises
        ------
        ValueError
            For an unknown ``backend``, or a site with no peak data.
        requests.RequestException
            If the peak request fails.
        """
        if backend == "nwis-legacy":
            df, name, area = _download_peak_flow_rdb(self._site_no)
            if name is not None:
                self._site_name = name
            if area is not None:
                self._drainage_area = area
        else:
            from flowfreq.peak_sources import get_backend

            try:
                source = get_backend(backend)
            except KeyError as exc:
                raise ValueError(str(exc.args[0])) from None
            df = source.fetch_peaks(self._site_no)
            if backend == "waterdata-ogc":
                self._site_metadata_from_waterdata()

        self._peak_data = df.reset_index(drop=True)

        if "period_of_record" in self.__dict__:
            del self.__dict__["period_of_record"]

        return self._peak_data

    def _site_metadata_from_waterdata(self) -> None:
        """Set site name and drainage area from the OGC ``monitoring-locations`` record.

        The counterpart of the legacy RDB header's ``Station name`` and
        ``Drainage area`` lines. A failure is logged, not raised: the peaks
        are what was asked for.
        """
        # Deferred: flowfreq.waterdata imports from this module.
        from flowfreq.waterdata import fetch_monitoring_location

        try:
            location = fetch_monitoring_location(self._site_no, timeout=30)
        except requests.RequestException as exc:
            logger.warning(
                "Site %s: monitoring-location request failed (%s); site name and drainage "
                "area left unset",
                self._site_no,
                exc,
            )
            return
        name = location.get("monitoring_location_name")
        if isinstance(name, str) and name.strip():
            self._site_name = name.strip()
        area = location.get("drainage_area")
        try:
            area_f = float(area) if area is not None else float("nan")
        except (TypeError, ValueError):
            area_f = float("nan")
        if math.isfinite(area_f):
            self._drainage_area = area_f

    def __repr__(self) -> str:
        return f"USGSgage(site_no='{self._site_no}', name='{self._site_name}')"


#: Column layout of an instantaneous-value frame, in order.
_IV_COLUMNS: Tuple[str, ...] = (
    "flow_cfs",
    "datetime_local",
    "tz_cd",
    "qualification_code",
)

#: NWIS instantaneous parameter codes this module retrieves, mapped to the
#: value-column name and human description used for each. Discharge (00060)
#: and gage height (00065) share the whole retrieval path -- chunking, the UTC
#: axis, the tz mapping, multi-sensor disambiguation -- so the path is
#: parameterized over this table rather than duplicated per parameter, which
#: would let the two copies drift.
IV_PARAMETERS: Dict[str, Tuple[str, str]] = {
    "00060": ("flow_cfs", "discharge"),
    "00065": ("gage_height_ft", "gage height"),
}


#: Instantaneous-value backends :meth:`USGSgage.download_instantaneous_flow`
#: accepts. The default is :data:`flowfreq.peak_sources.DEFAULT_BACKEND`.
IV_BACKENDS: Tuple[str, ...] = ("waterdata-ogc", "nwis-legacy")

_LEGACY_TS_ID = re.compile(r"^\d+$")
_WATERDATA_TS_ID = re.compile(r"^[0-9a-f]{32}$")


def check_ts_id_form(ts_id: Optional[str], backend: str, param_cd: str = "00060") -> None:
    """Refuse a ``ts_id`` written for the other instantaneous backend.

    The two backends name a series differently: legacy NWIS by its DD number
    (``"60629"``), the Water Data OGC API by a 32-hex ``time_series_id``
    (``"15beb94252164ce0b9a77a90edce7528"``, dashes optional). Since the
    default backend changed from legacy to OGC (issue #29), a caller still
    passing a DD number would otherwise meet a "does not match" error only
    after several requests, or -- worse -- be tempted to drop ``ts_id`` and
    take whatever a single-sensor window returns. This check runs before any
    request and says which backend the value belongs to.

    Parameters
    ----------
    ts_id : str, optional
        The caller's ``ts_id``; ``None`` always passes.
    backend : str
        ``"waterdata-ogc"`` or ``"nwis-legacy"``.
    param_cd : str
        For the message.

    Raises
    ------
    ValueError
        ``ts_id`` is a DD number on ``"waterdata-ogc"``, or a 32-hex
        ``time_series_id`` on ``"nwis-legacy"``.
    """
    if ts_id is None:
        return
    text = str(ts_id).strip()
    if backend == "waterdata-ogc" and _LEGACY_TS_ID.match(text):
        raise ValueError(
            f"ts_id={ts_id!r} is a legacy NWIS DD number, but the backend is "
            f"'waterdata-ogc' (the default since issue #29), where ts_id is the 32-hex "
            f"time_series_id. Either pass backend='nwis-legacy' to keep using DD numbers, "
            f"or look the series up with "
            f"flowfreq.waterdata.list_instantaneous_series(site_no, {param_cd!r}) -- its "
            f"sublocation_identifier/web_description match the description legacy NWIS "
            f"prints beside the DD number."
        )
    if backend == "nwis-legacy" and _WATERDATA_TS_ID.match(text.lower().replace("-", "")):
        raise ValueError(
            f"ts_id={ts_id!r} is a Water Data OGC time_series_id, but the backend is "
            f"'nwis-legacy', where ts_id is the NWIS DD number (e.g. '60629'). Pass "
            f"backend='waterdata-ogc' to use it."
        )


def _iv_value_column(param_cd: str) -> str:
    """Value-column name for an NWIS instantaneous parameter code."""
    if param_cd not in IV_PARAMETERS:
        raise ValueError(
            f"Unsupported instantaneous parameter code {param_cd!r}; known: "
            f"{sorted(IV_PARAMETERS)}"
        )
    return IV_PARAMETERS[param_cd][0]


def _empty_iv_frame(value_col: str = "flow_cfs") -> pd.DataFrame:
    """Return a correctly typed, empty instantaneous-value frame.

    Used for chunks that legitimately contain no records, so that concatenating
    chunks never has to special-case dtype or index-tz mismatches.
    """
    return pd.DataFrame(
        {
            value_col: pd.Series(dtype=float),
            "datetime_local": pd.Series(dtype="datetime64[ns]"),
            "tz_cd": pd.Series(dtype=object),
            "qualification_code": pd.Series(dtype=object),
        },
        index=pd.DatetimeIndex([], tz="UTC", name="datetime"),
    )


def _is_no_data_response(text: str) -> bool:
    """Detect the NWIS 400 response that means 'nothing here', not 'broken'.

    The instantaneous-value service answers a window it has no records for with
    HTTP 400 and an explanatory body rather than an empty 200, so the status
    code alone cannot distinguish a gap in the record from a real failure.
    """
    lowered = text.lower()
    return "no sites" in lowered or "no data" in lowered


def _validate_daily_flow_range(start_date: Optional[str], end_date: Optional[str]) -> None:
    """Raise a clear ValueError for a malformed or reversed caller-supplied range.

    Otherwise a bad date string or a reversed range reaches NWIS exactly as given,
    and the resulting failure is indistinguishable from a network problem.
    """
    start_ts = pd.Timestamp(start_date) if start_date is not None else None
    end_ts = pd.Timestamp(end_date) if end_date is not None else None

    if (start_ts is not None and pd.isna(start_ts)) or (end_ts is not None and pd.isna(end_ts)):
        raise ValueError(f"Could not parse date range {start_date!r} to {end_date!r}")
    if start_ts is not None and end_ts is not None and start_ts > end_ts:
        raise ValueError(f"start_date {start_date} is after end_date {end_date}")


def _chunk_date_range(start_date: str, end_date: str, chunk_years: int) -> List[Tuple[str, str]]:
    """Split an inclusive date range into consecutive chunks of whole years.

    Parameters
    ----------
    start_date, end_date : str
        Inclusive bounds, ``YYYY-MM-DD``.
    chunk_years : int
        Length of each chunk in years.

    Returns
    -------
    list of (str, str)
        Inclusive, non-overlapping, gap-free ``(start, end)`` pairs covering
        exactly the requested range.
    """
    start_ts = pd.Timestamp(start_date)
    end_ts = pd.Timestamp(end_date)

    if pd.isna(start_ts) or pd.isna(end_ts):
        raise ValueError(f"Could not parse date range {start_date!r} to {end_date!r}")
    if start_ts > end_ts:
        raise ValueError(f"start_date {start_date} is after end_date {end_date}")

    chunks: List[Tuple[str, str]] = []
    current = start_ts
    while current <= end_ts:
        chunk_end = min(
            current + pd.DateOffset(years=chunk_years) - pd.Timedelta(days=1),
            end_ts,
        )
        chunks.append((current.strftime("%Y-%m-%d"), chunk_end.strftime("%Y-%m-%d")))
        current = chunk_end + pd.Timedelta(days=1)

    return chunks


def _resolve_value_column(df: pd.DataFrame, ts_id: Optional[str], param_cd: str) -> str:
    """Pick the single value column for one parameter out of a parsed RDB table.

    Value columns are named ``<DD>_<param_cd>``; their qualifier twins end
    ``_<param_cd>_cd``. A site that reports the parameter from more than one
    sensor yields more than one such column, and choosing between them silently
    would hand back a plausible-looking series from the wrong instrument.

    Parameters
    ----------
    df : pd.DataFrame
        Parsed RDB table, all columns as strings.
    ts_id : str, optional
        NWIS time-series (DD) identifier selecting one series.
    param_cd : str
        NWIS parameter code, e.g. ``"00060"`` for discharge or ``"00065"`` for
        gage height.

    Returns
    -------
    str
        Name of the value column to read.

    Raises
    ------
    ValueError
        No matching column, several with no ``ts_id`` to choose between them,
        or a ``ts_id`` matching none of them.
    """
    description = IV_PARAMETERS.get(param_cd, ("", param_cd))[1]
    value_cols = [c for c in df.columns if c.endswith(f"_{param_cd}")]

    if ts_id is not None:
        wanted = f"{ts_id}_{param_cd}"
        if wanted not in value_cols:
            raise ValueError(
                f"ts_id={ts_id!r} does not match any {description} series in this "
                f"response; available: {value_cols or 'none'}"
            )
        return wanted

    if not value_cols:
        raise ValueError(
            f"NWIS returned {len(df)} instantaneous records but no {description} "
            f"({param_cd}) column; columns were: {list(df.columns)}"
        )

    if len(value_cols) > 1:
        example = value_cols[0].split("_")[0]
        raise ValueError(
            f"Site reports {len(value_cols)} separate {param_cd} time series "
            f"({', '.join(value_cols)}). Pass ts_id to choose one "
            f"(e.g. ts_id={example!r}) rather than accepting an arbitrary sensor."
        )

    return value_cols[0]


def _qualifier_column(df: pd.DataFrame, value_col: str) -> pd.Series:
    """Return the NWIS data-qualifier column paired with a value column.

    Qualifier columns are the value column's name plus ``_cd``. They are not
    guaranteed to be present, so an all-empty column stands in when absent.
    """
    qual_col = f"{value_col}_cd"
    if qual_col in df.columns:
        return df[qual_col].fillna("").astype(str)
    return pd.Series("", index=df.index)


def _parse_iv_rdb(text: str, ts_id: Optional[str] = None, param_cd: str = "00060") -> pd.DataFrame:
    """Parse an NWIS instantaneous-value RDB payload into a UTC-indexed frame.

    Separated from the HTTP call so the parsing rules — column discovery,
    time-zone mapping, multi-sensor detection — can be tested against captured
    payloads without a network round trip.

    Parameters
    ----------
    text : str
        Raw RDB text as returned by the NWIS instantaneous-values service.
    ts_id : str, optional
        NWIS time-series (DD) identifier, to select one series at a site that
        reports the parameter from more than one sensor.
    param_cd : str
        NWIS parameter code to extract, one of :data:`IV_PARAMETERS`. Default
        ``"00060"`` (discharge); ``"00065"`` is gage height. The output column
        is named accordingly.

    Returns
    -------
    pd.DataFrame
        Indexed by tz-aware UTC datetime, with the columns described by
        :meth:`USGSgage.download_instantaneous_flow` -- except that the value
        column is named for `param_cd`. Empty if the payload carried no
        records.

    Raises
    ------
    ValueError
        Records are present but no column for `param_cd` was found, more than
        one such series is present and ``ts_id`` did not resolve it, a
        time-zone abbreviation could not be mapped to a UTC offset, or
        `param_cd` is not a supported parameter.
    """
    value_col = _iv_value_column(param_cd)
    lines = text.split("\n")
    data_lines = [line for line in lines if not line.startswith("#") and line.strip()]

    # An RDB payload needs a header, a format-spec row, and at least one record.
    if len(data_lines) < 3:
        return _empty_iv_frame(value_col)

    header_idx = next((i for i, line in enumerate(data_lines) if "datetime" in line.lower()), None)
    if header_idx is None:
        return _empty_iv_frame(value_col)

    df = pd.read_csv(
        StringIO("\n".join(data_lines[header_idx:])), sep="\t", skiprows=[1], dtype=str
    )
    if df.empty:
        return _empty_iv_frame(value_col)

    date_col = next((c for c in df.columns if c.lower() == "datetime"), None)
    if date_col is None:
        return _empty_iv_frame(value_col)

    raw_col = _resolve_value_column(df, ts_id, param_cd)

    if "tz_cd" not in df.columns:
        raise ValueError(
            "NWIS instantaneous response has no tz_cd column, so its local "
            "timestamps cannot be placed on a UTC axis"
        )

    tz_codes = df["tz_cd"].astype(str).str.strip()
    offsets = tz_codes.map(NWIS_TZ_OFFSETS)
    unknown = sorted(set(tz_codes[offsets.isna()]))
    if unknown:
        raise ValueError(
            f"Unrecognized NWIS time-zone code(s) {unknown} for this site. Add them to "
            f"flowfreq.usgs.NWIS_TZ_OFFSETS; records are not dropped silently."
        )

    frame = pd.DataFrame(
        {
            value_col: pd.to_numeric(df[raw_col], errors="coerce").astype(float),
            "datetime_local": pd.to_datetime(df[date_col], errors="coerce"),
            "tz_cd": tz_codes,
            "qualification_code": _qualifier_column(df, raw_col),
            "_offset_hours": offsets.astype(float),
        }
    )
    frame = frame.dropna(subset=[value_col, "datetime_local"])
    if frame.empty:
        return _empty_iv_frame(value_col)

    utc = frame["datetime_local"] - pd.to_timedelta(frame["_offset_hours"], unit="h")
    index = pd.DatetimeIndex(utc).tz_localize("UTC")
    index.name = "datetime"

    frame = frame.drop(columns=["_offset_hours"])
    frame.index = index
    return frame[[value_col, *_IV_COLUMNS[1:]]]


def _parse_peak_dt(peak_dt: pd.Series, site_no: str = "") -> pd.Series:
    """Parse NWIS peak dates, including partial dates, the way peakfq does.

    NWIS encodes an unknown day or month as ``00`` (``1897-03-00``,
    ``1968-00-00``), which :func:`pandas.to_datetime` rejects. Coercing those to
    NaT, as this parser once did, dropped the row -- and with it, typically, the
    historic peaks a Bulletin 17C analysis depends on most. Follows peakfq
    8.1.0's reader (``vendor/peakfqr/R/DataReaderFunctions_shinyapp.R``): an
    unknown month becomes January, which keeps the peak in the water year equal
    to its calendar year, and an unknown day becomes the 1st, which cannot
    change the water year. The returned date is therefore a placeholder where
    the day or month was unknown; the water year derived from it is exact.

    Parameters
    ----------
    peak_dt : pandas.Series
        ``peak_dt`` strings from the NWIS peak RDB (``YYYY-MM-DD``).
    site_no : str, optional
        Used only in the log message.

    Returns
    -------
    pandas.Series
        Datetimes; NaT only where the string is not a date at all.
    """
    parts = peak_dt.astype(str).str.strip().str.extract(r"^(\d{4})-(\d{2})-(\d{2})$")
    year = pd.to_numeric(parts[0], errors="coerce")
    month = pd.to_numeric(parts[1], errors="coerce")
    day = pd.to_numeric(parts[2], errors="coerce")

    partial = (month == 0) | (day == 0)
    if partial.any():
        logger.info(
            "Site %s: %d peak date(s) with unknown day or month; placeholder "
            "month 1 / day 1 used, as peakfq does",
            site_no,
            int(partial.sum()),
        )

    return pd.to_datetime(
        pd.DataFrame({"year": year, "month": month.replace(0, 1), "day": day.replace(0, 1)}),
        errors="coerce",
    )


def _download_peak_flow_rdb(
    site_no: str, timeout: float = 30
) -> Tuple[pd.DataFrame, Optional[str], Optional[float]]:
    """Fetch and parse a site's annual peaks from the legacy NWIS peak RDB service.

    The one implementation behind both ``USGSgage.download_peak_flow(
    backend="nwis-legacy")`` and :class:`flowfreq.peak_sources.LegacyNwisBackend`,
    so neither calls the other.

    Parameters
    ----------
    site_no : str
        USGS site number.
    timeout : float
        Request timeout, seconds.

    Returns
    -------
    tuple
        ``(frame, station_name, drainage_area_sqmi)``. The frame has columns
        ``water_year``, ``peak_date``, ``peak_flow_cfs``, ``qualification_code``;
        the name and area come from the RDB header and are ``None`` when it
        does not give them.

    Raises
    ------
    ValueError
        If the response holds no peak rows.
    requests.RequestException
        If the request fails.
    """
    params = {
        "site_no": site_no,
        "agency_cd": "USGS",
        "format": "rdb",
    }

    response = requests.get(USGSgage.BASE_URL_PEAKS, params=params, timeout=timeout)
    response.raise_for_status()

    lines = response.text.split("\n")
    data_lines = [l for l in lines if not l.startswith("#") and l.strip()]

    if len(data_lines) < 2:
        raise ValueError(f"No peak flow data found for site {site_no}")

    station_name: Optional[str] = None
    drainage_area: Optional[float] = None
    for line in lines:
        if "#" in line:
            if "DRAINAGE AREA" in line.upper():
                try:
                    parts = line.split(":")[-1].strip()
                    drainage_area = float(parts.split()[0])
                except (ValueError, IndexError):
                    pass
            if "STATION NAME" in line.upper():
                station_name = line.split(":")[-1].strip()

    # Codes as strings: a column holding only numeric codes (e.g. "7") and
    # blanks is otherwise inferred as float, and "7" comes back as "7.0".
    df = pd.read_csv(
        StringIO("\n".join(data_lines)),
        sep="\t",
        skiprows=[1],
        dtype={"peak_dt": str, "peak_cd": str, "gage_ht_cd": str},
    )

    df = df[df["agency_cd"] == "USGS"].copy()
    df["peak_date"] = _parse_peak_dt(df["peak_dt"], site_no)
    df["peak_flow_cfs"] = pd.to_numeric(df["peak_va"], errors="coerce")

    df["water_year"] = df["peak_date"].apply(
        lambda x: x.year + 1 if x.month >= 10 else x.year if pd.notna(x) else np.nan
    )

    if "peak_cd" in df.columns:
        df["qualification_code"] = df["peak_cd"].fillna("")
    else:
        df["qualification_code"] = ""

    df = df[["water_year", "peak_date", "peak_flow_cfs", "qualification_code"]].dropna(
        subset=["water_year", "peak_flow_cfs"]
    )
    df["water_year"] = df["water_year"].astype(int)
    return df.reset_index(drop=True), station_name, drainage_area


def fetch_nwis_peaks(site_no: str, backend: str = DEFAULT_BACKEND) -> List[Dict]:
    """
    Fetch peak flow records for a single USGS site.

    Parameters
    ----------
    site_no : str
        USGS site number
    backend : str, default :data:`flowfreq.peak_sources.DEFAULT_BACKEND`
        Passed to :meth:`USGSgage.download_peak_flow`. The default is now the
        Water Data OGC API; ``"nwis-legacy"`` is the service this always used
        before. The name ``fetch_nwis_peaks`` is kept for compatibility.

    Returns
    -------
    list of dict
        Peak flow records for the site
    """
    gage = USGSgage(site_no)
    gage.download_peak_flow(backend=backend)
    records = []
    for _, row in gage.peak_data.iterrows():
        records.append(
            {
                "year": int(row["water_year"]),
                "flow": float(row["peak_flow_cfs"]),
                "source": "USGS",
            }
        )
    return records


def fetch_nwis_batch(
    sites: List[str], workers: int = 6, backend: str = DEFAULT_BACKEND
) -> Tuple[Dict[str, List[Dict]], Dict[str, str]]:
    """
    Fetch peak flow records for multiple USGS sites in parallel.

    Parameters
    ----------
    sites : list of str
        USGS site numbers
    workers : int
        Number of parallel workers (default: 6)
    backend : str, default :data:`flowfreq.peak_sources.DEFAULT_BACKEND`
        Passed to :func:`fetch_nwis_peaks` for every site.

    Returns
    -------
    tuple
        (successful_results, errors) where:
        - successful_results: dict mapping site_no to list of records
        - errors: dict mapping site_no to error message
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    results: Dict[str, List[Dict]] = {}
    errors: Dict[str, str] = {}

    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_site = {executor.submit(fetch_nwis_peaks, site, backend): site for site in sites}

        for future in as_completed(future_to_site):
            site = future_to_site[future]
            try:
                results[site] = future.result()
            except Exception as e:
                errors[site] = str(e)

    return results, errors
