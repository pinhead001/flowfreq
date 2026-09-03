"""
hydrolib.usgs - USGS data retrieval
"""

from __future__ import annotations

import math
from functools import cached_property
from io import StringIO
from pathlib import Path
from typing import ClassVar, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import requests

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
        """Most recently downloaded instantaneous (unit-value) series."""
        return self._instantaneous_data

    @instantaneous_data.setter
    def instantaneous_data(self, value: pd.DataFrame):
        self._instantaneous_data = value

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
            if self._last_api_error:
                self._last_api_error += f"; {str(e)}"
            else:
                self._last_api_error = str(e)

    def download_daily_flow(self, start_date: str = None, end_date: str = None) -> pd.DataFrame:
        """Download mean daily streamflow data from USGS."""
        params = {
            "format": "rdb",
            "sites": self._site_no,
            "parameterCd": "00060",
            "statCd": "00003",
        }
        if start_date:
            params["startDT"] = start_date
        if end_date:
            params["endDT"] = end_date

        response = requests.get(self.BASE_URL_DAILY, params=params)
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
        df = df.set_index("date")

        self._daily_data = df
        return df

    def download_instantaneous_flow(  # pylint: disable=too-many-arguments
        self,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
        *,
        tz: Optional[str] = None,
        chunk_years: int = 1,
        ts_id: Optional[str] = None,
        timeout: int = 60,
    ) -> pd.DataFrame:
        """Download instantaneous (unit-value) streamflow data from USGS NWIS.

        Retrieves parameter 00060 (discharge, cfs) from the NWIS instantaneous
        values service. Unlike the daily service there is no statistic code:
        every recorded value is returned, typically at a 15-minute interval.

        Parameters
        ----------
        start_date : str, optional
            First date to retrieve, ``YYYY-MM-DD``. Defaults to the start of
            the site's instantaneous period of record, which is looked up from
            the NWIS site service if not already known.
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
            NWIS time-series (DD) identifier, used to disambiguate sites that
            report discharge from more than one sensor. If a site returns
            multiple 00060 series and this is not given, the call raises
            rather than silently picking one.
        timeout : int
            Per-request timeout in seconds. Default 60.

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
              ``A`` approved, ``e`` estimated)

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
        the daily cycle. The functions in :mod:`hydrolib.regime` take an
        explicit time zone for this reason.

        **Storage.** For a series of this size Parquet is the format worth
        reaching for — it round-trips the tz-aware index and float dtypes
        exactly, where CSV loses both, and it is several times smaller. See
        :func:`hydrolib.flowio.save_flow_frame`.

        Examples
        --------
        >>> gage = USGSgage("12449950")
        >>> iv = gage.download_instantaneous_flow("2022-06-01", "2022-09-30")
        >>> iv.index.tz is not None
        True
        """
        if chunk_years < 1:
            raise ValueError(f"chunk_years must be >= 1, got {chunk_years}")

        start_date, end_date = self._resolve_iv_window(start_date, end_date)
        chunks = _chunk_date_range(start_date, end_date, chunk_years)

        frames: List[pd.DataFrame] = []
        for chunk_start, chunk_end in chunks:
            chunk = self._request_iv_chunk(chunk_start, chunk_end, ts_id, timeout)
            if not chunk.empty:
                frames.append(chunk)

        combined = pd.concat(frames) if frames else _empty_iv_frame()

        if combined.empty:
            raise NoInstantaneousDataError(
                f"No instantaneous discharge (parameter 00060) found for site "
                f"{self._site_no} between {start_date} and {end_date}"
                + (
                    f"; the site's instantaneous record runs "
                    f"{self._iv_por_start} to {self._iv_por_end}"
                    if self._iv_por_start
                    else "; NWIS lists no instantaneous discharge series for this site"
                )
            )

        combined = combined[~combined.index.duplicated(keep="first")].sort_index()

        if tz is not None:
            combined.index = combined.index.tz_convert(tz)

        self._instantaneous_data = combined
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

    def _request_iv_chunk(
        self, chunk_start: str, chunk_end: str, ts_id: Optional[str], timeout: int
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
            "parameterCd": "00060",
            "startDT": chunk_start,
            "endDT": chunk_end,
        }

        try:
            response = requests.get(self.BASE_URL_IV, params=params, timeout=timeout)
            if response.status_code == 400 and _is_no_data_response(response.text):
                return _empty_iv_frame()
            response.raise_for_status()
        except requests.RequestException as exc:
            raise requests.RequestException(
                f"Instantaneous-value request failed for site {self._site_no}, "
                f"window {chunk_start} to {chunk_end}: {exc}"
            ) from exc

        return _parse_iv_rdb(response.text, ts_id=ts_id)

    def download_peak_flow(self) -> pd.DataFrame:
        """Download annual peak streamflow data from USGS."""
        params = {
            "site_no": self._site_no,
            "agency_cd": "USGS",
            "format": "rdb",
        }

        response = requests.get(self.BASE_URL_PEAKS, params=params)
        response.raise_for_status()

        lines = response.text.split("\n")
        data_lines = [l for l in lines if not l.startswith("#") and l.strip()]

        if len(data_lines) < 2:
            raise ValueError(f"No peak flow data found for site {self._site_no}")

        for line in lines:
            if "#" in line:
                if "DRAINAGE AREA" in line.upper():
                    try:
                        parts = line.split(":")[-1].strip()
                        self._drainage_area = float(parts.split()[0])
                    except (ValueError, IndexError):
                        pass
                if "STATION NAME" in line.upper():
                    self._site_name = line.split(":")[-1].strip()

        df = pd.read_csv(StringIO("\n".join(data_lines)), sep="\t", skiprows=[1])

        df = df[df["agency_cd"] == "USGS"].copy()
        df["peak_date"] = pd.to_datetime(df["peak_dt"], errors="coerce")
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

        self._peak_data = df.reset_index(drop=True)

        if "period_of_record" in self.__dict__:
            del self.__dict__["period_of_record"]

        return self._peak_data

    def __repr__(self) -> str:
        return f"USGSgage(site_no='{self._site_no}', name='{self._site_name}')"


#: Column layout of an instantaneous-value frame, in order.
_IV_COLUMNS: Tuple[str, ...] = (
    "flow_cfs",
    "datetime_local",
    "tz_cd",
    "qualification_code",
)


def _empty_iv_frame() -> pd.DataFrame:
    """Return a correctly typed, empty instantaneous-value frame.

    Used for chunks that legitimately contain no records, so that concatenating
    chunks never has to special-case dtype or index-tz mismatches.
    """
    return pd.DataFrame(
        {
            "flow_cfs": pd.Series(dtype=float),
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


def _resolve_flow_column(df: pd.DataFrame, ts_id: Optional[str]) -> str:
    """Pick the single discharge column out of a parsed RDB table.

    Value columns are named ``<DD>_00060``; their qualifier twins end
    ``_00060_cd``. A site that reports discharge from more than one sensor
    yields more than one such column, and choosing between them silently would
    hand back a plausible-looking series from the wrong instrument.

    Parameters
    ----------
    df : pd.DataFrame
        Parsed RDB table, all columns as strings.
    ts_id : str, optional
        NWIS time-series (DD) identifier selecting one series.

    Returns
    -------
    str
        Name of the discharge column to read.

    Raises
    ------
    ValueError
        No discharge column, several with no ``ts_id`` to choose between them,
        or a ``ts_id`` matching none of them.
    """
    flow_cols = [c for c in df.columns if c.endswith("_00060")]

    if ts_id is not None:
        wanted = f"{ts_id}_00060"
        if wanted not in flow_cols:
            raise ValueError(
                f"ts_id={ts_id!r} does not match any discharge series in this response; "
                f"available: {flow_cols or 'none'}"
            )
        return wanted

    if not flow_cols:
        raise ValueError(
            f"NWIS returned {len(df)} instantaneous records but no discharge (00060) "
            f"column; columns were: {list(df.columns)}"
        )

    if len(flow_cols) > 1:
        example = flow_cols[0].split("_")[0]
        raise ValueError(
            f"Site reports {len(flow_cols)} separate 00060 time series "
            f"({', '.join(flow_cols)}). Pass ts_id to choose one "
            f"(e.g. ts_id={example!r}) rather than accepting an arbitrary sensor."
        )

    return flow_cols[0]


def _qualifier_column(df: pd.DataFrame, flow_col: str) -> pd.Series:
    """Return the NWIS data-qualifier column paired with a discharge column.

    Qualifier columns are the value column's name plus ``_cd``. They are not
    guaranteed to be present, so an all-empty column stands in when absent.
    """
    qual_col = f"{flow_col}_cd"
    if qual_col in df.columns:
        return df[qual_col].fillna("").astype(str)
    return pd.Series("", index=df.index)


def _parse_iv_rdb(text: str, ts_id: Optional[str] = None) -> pd.DataFrame:
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
        reports discharge from more than one sensor.

    Returns
    -------
    pd.DataFrame
        Indexed by tz-aware UTC datetime, with the columns described by
        :meth:`USGSgage.download_instantaneous_flow`. Empty if the payload
        carried no records.

    Raises
    ------
    ValueError
        Records are present but no discharge column was found, more than one
        discharge series is present and ``ts_id`` did not resolve it, or a
        time-zone abbreviation could not be mapped to a UTC offset.
    """
    lines = text.split("\n")
    data_lines = [line for line in lines if not line.startswith("#") and line.strip()]

    # An RDB payload needs a header, a format-spec row, and at least one record.
    if len(data_lines) < 3:
        return _empty_iv_frame()

    header_idx = next((i for i, line in enumerate(data_lines) if "datetime" in line.lower()), None)
    if header_idx is None:
        return _empty_iv_frame()

    df = pd.read_csv(
        StringIO("\n".join(data_lines[header_idx:])), sep="\t", skiprows=[1], dtype=str
    )
    if df.empty:
        return _empty_iv_frame()

    date_col = next((c for c in df.columns if c.lower() == "datetime"), None)
    if date_col is None:
        return _empty_iv_frame()

    flow_col = _resolve_flow_column(df, ts_id)

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
            f"hydrolib.usgs.NWIS_TZ_OFFSETS; records are not dropped silently."
        )

    frame = pd.DataFrame(
        {
            "flow_cfs": pd.to_numeric(df[flow_col], errors="coerce").astype(float),
            "datetime_local": pd.to_datetime(df[date_col], errors="coerce"),
            "tz_cd": tz_codes,
            "qualification_code": _qualifier_column(df, flow_col),
            "_offset_hours": offsets.astype(float),
        }
    )
    frame = frame.dropna(subset=["flow_cfs", "datetime_local"])
    if frame.empty:
        return _empty_iv_frame()

    utc = frame["datetime_local"] - pd.to_timedelta(frame["_offset_hours"], unit="h")
    index = pd.DatetimeIndex(utc).tz_localize("UTC")
    index.name = "datetime"

    frame = frame.drop(columns=["_offset_hours"])
    frame.index = index
    return frame[list(_IV_COLUMNS)]


def fetch_nwis_peaks(site_no: str) -> List[Dict]:
    """
    Fetch peak flow records for a single USGS site.

    Parameters
    ----------
    site_no : str
        USGS site number

    Returns
    -------
    list of dict
        Peak flow records for the site
    """
    gage = USGSgage(site_no)
    gage.download_peak_flow()
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
    sites: List[str], workers: int = 6
) -> Tuple[Dict[str, List[Dict]], Dict[str, str]]:
    """
    Fetch peak flow records for multiple USGS sites in parallel.

    Parameters
    ----------
    sites : list of str
        USGS site numbers
    workers : int
        Number of parallel workers (default: 6)

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
        future_to_site = {executor.submit(fetch_nwis_peaks, site): site for site in sites}

        for future in as_completed(future_to_site):
            site = future_to_site[future]
            try:
                results[site] = future.result()
            except Exception as e:
                errors[site] = str(e)

    return results, errors
