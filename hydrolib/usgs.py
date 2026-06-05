"""
hydrolib.usgs - USGS data retrieval
"""

from __future__ import annotations

from functools import cached_property
from io import StringIO
from pathlib import Path
from typing import ClassVar, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import requests


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
        # Try multiple possible locations
        candidates = [
            # Relative to this module (for editable installs)
            Path(__file__).parent.parent / "data" / "gage_attributes.csv",
            # Relative to current working directory
            Path.cwd() / "data" / "gage_attributes.csv",
            # In hydrolib package data
            Path(__file__).parent / "data" / "gage_attributes.csv",
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


class USGSgage:
    """Class to handle USGS gage data retrieval and storage."""

    BASE_URL_DAILY: ClassVar[str] = "https://waterservices.usgs.gov/nwis/dv/"
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

        # Call 1: Get site metadata (name, drainage area) with siteOutput=expanded
        if self._site_name is None or self._drainage_area is None:
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

                    if (
                        self._drainage_area is None
                        and "drain_area_va" in df.columns
                        and len(df) > 0
                    ):
                        try:
                            self._drainage_area = float(df["drain_area_va"].iloc[0])
                        except (ValueError, TypeError):
                            pass
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
                    except:
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
