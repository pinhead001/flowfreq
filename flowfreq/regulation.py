"""Regulation / urbanization screen per gage, and the B17C refusal it drives.

Roadmap ``docs/MASTER_ROADMAP.md`` §1.1, issue #32. A gage is classed
``reference``, ``regulated``, ``urban`` or ``unknown`` from three kinds of
evidence:

1. **Basin attributes** (offline, packaged): ``data/regulation_screen.csv.gz``,
   built from GAGES-II (Falcone 2011, https://doi.org/10.5066/P96CPHOT) by
   ``tools/build_regulation_screen.py``. The rules and their sources are in
   that tool's docstring: GAGES-II ``Ref`` class -> reference; NID 2009 dam
   storage normalized by mean annual runoff above 127.8 days (Dudley and
   others 2018, https://doi.org/10.5066/P9AEGXY0, 75th percentile) ->
   regulated; NLCD 2006 impervious cover above 5 percent (Mastin and others
   2016, SIR 2016-5118 p. 23) -> urban. 9,322 gages; a gage GAGES-II does not
   cover has no table entry.
2. **Peak codes** (per record, at run time): any peak with NWIS code 6,
   "Discharge affected by regulation or diversion", makes the record
   ``regulated``. Mastin and others (2016, p. 20) treat code-6 gages as
   regulated and code-5 gages ("affected to unknown degree") as unregulated;
   this screen does the same, so code 5 is recorded but never decisive.
3. **Current basin data** (opt-in, online, per site): the same two rules
   re-evaluated on today's data instead of GAGES-II's 2009/2006 snapshot, for
   any gage StreamStats can delineate, in GAGES-II or not. See
   :func:`classify_site`'s ``use_current_nid`` and ``use_current_impervious``.

   - *Dam storage* is the current National Inventory of Dams
     (https://nid.sec.usace.army.mil/api/nation/csv, :func:`load_nid`) summed
     over the dams inside the StreamStats watershed polygon
     (:func:`nid_storage_in_basin`), normalized exactly as the packaged rule
     is: days = storage per unit area / mean annual runoff * 365.
   - *Impervious cover* is the newest NLCD impervious epoch StreamStats
     computes for the basin (:data:`NLCD_IMPERVIOUS_CODES`), with its year.
     Which epochs exist depends on the StreamStats region; see that constant.
   - StreamStats refuses to delineate at some gages with the warning "River is
     regulated, streamflow characteristics not valid". That is USGS's own
     regulated-reach layer, and it is taken as decisive.

Bulletin 17C itself sets no numeric regulation threshold (England and others,
2019, p. 36, "Regulated Flow Frequency"); it says regulated-flow methods need
national guidance it does not give. That is why a regulated record is refused
rather than fitted: :func:`require_unregulated` raises
:class:`RegulatedRecordError` unless the caller passes an explicit override,
and the override is recorded in the returned provenance. A gage with no
classification at all proceeds, with a warning.
"""

from __future__ import annotations

import gzip
import json
import logging
import math
import os
import time
from dataclasses import dataclass, field
from datetime import date, datetime, timezone
from enum import Enum
from functools import lru_cache
from importlib import resources
from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Iterable, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import requests

from .peak_codes import parse_codes

if TYPE_CHECKING:  # pragma: no cover
    from .streamstats import StreamStatsCache, WatershedCharacteristics

logger = logging.getLogger(__name__)

#: Peak code that, on any peak, makes a record regulated (SIR 2016-5118 p. 20).
REGULATED_PEAK_CODE = "6"
#: Recorded as evidence only; SIR 2016-5118 p. 20 treats it as unregulated.
UNKNOWN_DEGREE_PEAK_CODE = "5"

#: Dudley and others (2018), 75th percentile of normalized dam storage, days.
#: The same threshold ``tools/build_regulation_screen.py`` applies.
REGULATED_STORAGE_DAYS = 127.8
#: Mastin and others (2016), SIR 2016-5118 p. 23, percent impervious.
URBAN_IMPERVIOUS_PCT = 5.0

SCREEN_FILE = "regulation_screen.csv.gz"
SCREEN_COLUMNS = (
    "site_no",
    "site_class",
    "gagesii_class",
    "ndams_2009",
    "stor_nid_2009_ml_km2",
    "norm_storage_days",
    "imperv_pct_2006",
    "runave7100_mm",
    "basis",
    "source",
    "source_date",
)


class SiteClass(str, Enum):
    """Regulation / urbanization class of a gage."""

    REFERENCE = "reference"
    REGULATED = "regulated"
    URBAN = "urban"
    UNKNOWN = "unknown"


class RegulatedRecordError(ValueError):
    """Bulletin 17C was asked to fit a record classified as regulated."""


class NIDUnavailableError(RuntimeError):
    """The National Inventory of Dams could not be downloaded and no cached copy exists."""


@dataclass(frozen=True)
class SiteClassification:
    """A gage's screen result and the evidence behind it.

    Attributes
    ----------
    site_no : str
    site_class : SiteClass
    basis : tuple of str
        Human-readable reasons, attribute evidence first, then current basin
        data, then peak codes.
    attributes : dict
        The screen-table row (empty when the gage is not in GAGES-II).
    n_peaks_code6, n_peaks_code5, n_peaks : int
        Peak-code counts, when peak codes were supplied.
    current : dict
        Current-data evidence from the opt-in refinement (empty when it was
        not requested): NID storage, impervious time series, their dates and
        sources, and any reason a piece could not be computed.
    """

    site_no: str
    site_class: SiteClass
    basis: Tuple[str, ...] = ()
    attributes: Dict[str, Any] = field(default_factory=dict)
    n_peaks_code6: int = 0
    n_peaks_code5: int = 0
    n_peaks: int = 0
    current: Dict[str, Any] = field(default_factory=dict)

    @property
    def in_table(self) -> bool:
        """Whether GAGES-II basin-attribute evidence was available."""
        return bool(self.attributes)

    @property
    def has_current_evidence(self) -> bool:
        """Whether the opt-in refinement produced a storage or impervious value."""
        return any(
            k in self.current
            for k in ("nid_norm_storage_days", "imperv_pct", "streamstats_regulated")
        )

    def as_dict(self) -> Dict[str, Any]:
        """Plain-dict form for provenance."""
        return {
            "site_no": self.site_no,
            "site_class": self.site_class.value,
            "basis": list(self.basis),
            "attributes": dict(self.attributes),
            "n_peaks_code6": self.n_peaks_code6,
            "n_peaks_code5": self.n_peaks_code5,
            "n_peaks": self.n_peaks,
            "current": dict(self.current),
        }


def _default_path() -> Path:
    return Path(str(resources.files("flowfreq") / "data" / SCREEN_FILE))


@lru_cache(maxsize=4)
def _load(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, dtype={"site_no": str, "source_date": str})
    missing = [c for c in SCREEN_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"regulation screen missing columns: {missing}")
    bad = set(df["site_class"]) - {c.value for c in SiteClass}
    if bad:
        raise ValueError(f"unknown site_class value(s): {sorted(bad)}")
    if df["site_no"].duplicated().any():
        raise ValueError("duplicate site_no in regulation screen")
    return df.set_index("site_no", drop=False)


def load_screen(path: Union[str, Path, None] = None) -> pd.DataFrame:
    """Load the packaged regulation screen table.

    Parameters
    ----------
    path : str or Path, optional
        Read this file instead of the packaged one.

    Returns
    -------
    pandas.DataFrame
        One row per GAGES-II gage, indexed by ``site_no``, with
        :data:`SCREEN_COLUMNS`.

    Raises
    ------
    ValueError
        On missing columns, unknown classes or duplicate site numbers.
    """
    return _load(str(path if path is not None else _default_path())).copy()


def _clean(value: Any) -> Any:
    if isinstance(value, float) and pd.isna(value):
        return None
    if hasattr(value, "item"):
        return value.item()
    return value


# ------------------------------------------------------- current NID storage --

#: National Inventory of Dams, national CSV download (USACE). Verified live
#: 2026-10-03: 67 MB, first line ``Data Last Updated:,2026-9-23``, then a header
#: and 92,766 rows.
NID_CSV_URL = "https://nid.sec.usace.army.mil/api/nation/csv"
#: Columns kept from the NID CSV, and the names they are stored under.
NID_SOURCE_COLUMNS: Dict[str, str] = {
    "NID ID": "nid_id",
    "Latitude": "latitude",
    "Longitude": "longitude",
    "NID Storage (Acre-Ft)": "nid_storage_acft",
    "Year Completed": "year_completed",
    "Is Associated Structure?": "associated",
}
NID_COLUMNS = ("nid_id", "latitude", "longitude", "nid_storage_acft", "year_completed")
#: Re-download the cached NID after this many days (USACE updates it continually).
NID_MAX_AGE_DAYS = 90
NID_CACHE_FILE = "nid_national.csv.gz"
NID_META_FILE = "nid_national.meta.json"

#: Exact conversions: 1 acre-foot = 43,560 ft^3 = 1233.48183754752 m^3; 1 ML = 1000 m^3.
ACRE_FT_TO_ML = 1.23348183754752
CFS_TO_M3S = 0.028316846592
SQ_MI_TO_KM2 = 2.589988110336


@dataclass(frozen=True)
class NIDInventory:
    """The national dam inventory, reduced to what the storage sum needs.

    Attributes
    ----------
    dams : pandas.DataFrame
        :data:`NID_COLUMNS`; one row per primary structure (associated
        structures, which repeat their parent dam's NID ID and storage, are
        dropped), with coordinates present.
    data_last_updated : str
        The NID's own "Data Last Updated" date (``YYYY-MM-DD``).
    retrieved : str
        UTC date the CSV was downloaded.
    """

    dams: pd.DataFrame
    data_last_updated: str
    retrieved: str


def nid_cache_dir() -> Path:
    """Where :func:`load_nid` caches the NID: ``$FLOWFREQ_CACHE/nid`` or ``~/.flowfreq/nid``."""
    root = os.environ.get("FLOWFREQ_CACHE")
    return (Path(root) if root else Path.home() / ".flowfreq") / "nid"


def parse_nid_csv(text: str) -> Tuple[pd.DataFrame, str]:
    """Parse the NID national CSV.

    Parameters
    ----------
    text : str
        The whole download: a ``Data Last Updated:,YYYY-M-D`` line, then a
        header row, then one row per structure.

    Returns
    -------
    (pandas.DataFrame, str)
        The dams (:data:`NID_COLUMNS`) and the last-updated date as ``YYYY-MM-DD``.

    Raises
    ------
    ValueError
        If the first line is not the last-updated line, or a needed column is missing.
    """
    first, _, rest = text.partition("\n")
    key, _, value = first.strip().partition(",")
    if not key.lower().startswith("data last updated"):
        raise ValueError(f"NID CSV does not start with its last-updated line: {first[:80]!r}")
    updated = datetime.strptime(value.strip().strip('"'), "%Y-%m-%d").date().isoformat()
    raw = pd.read_csv(StringIO(rest), dtype=str, keep_default_na=False)
    missing = [c for c in NID_SOURCE_COLUMNS if c not in raw.columns]
    if missing:
        raise ValueError(f"NID CSV missing columns: {missing}")
    df = raw[list(NID_SOURCE_COLUMNS)].rename(columns=NID_SOURCE_COLUMNS)
    df = df[df["associated"].str.strip().str.lower() != "yes"].drop(columns="associated")
    for col in ("latitude", "longitude", "nid_storage_acft", "year_completed"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["latitude", "longitude"])
    df["nid_storage_acft"] = df["nid_storage_acft"].fillna(0.0)
    df = df.drop_duplicates("nid_id").sort_values("nid_id").reset_index(drop=True)
    return df.loc[:, list(NID_COLUMNS)], updated


def _download_nid(timeout: float, retries: int) -> str:
    last: Optional[BaseException] = None
    for attempt in range(retries):
        try:
            r = requests.get(NID_CSV_URL, timeout=timeout)
            r.raise_for_status()
            return r.content.decode("utf-8-sig", errors="replace")
        except requests.RequestException as exc:
            # The NID host has been seen to fail DNS resolution intermittently.
            last = exc
            logger.warning("NID download failed (attempt %d/%d): %s", attempt + 1, retries, exc)
            if attempt < retries - 1:
                time.sleep(5.0 * (attempt + 1))
    raise NIDUnavailableError(f"could not download {NID_CSV_URL}: {last}") from last


def load_nid(
    cache_dir: Union[str, Path, None] = None,
    *,
    refresh: bool = False,
    max_age_days: int = NID_MAX_AGE_DAYS,
    timeout: float = 600.0,
    retries: int = 3,
) -> NIDInventory:
    """The current National Inventory of Dams, cached in the user cache.

    The first call downloads the 67 MB national CSV (:data:`NID_CSV_URL`) and
    keeps a 1.3 MB extract, ``nid_national.csv.gz``, under :func:`nid_cache_dir`;
    later calls read that. Nothing is written into the package.

    Parameters
    ----------
    cache_dir : str or Path, optional
        Defaults to :func:`nid_cache_dir`.
    refresh : bool
        Download even if the cache is fresh.
    max_age_days : int
        Re-download a cache older than this.
    timeout : float
        Download timeout, seconds.
    retries : int
        Download attempts before giving up.

    Returns
    -------
    NIDInventory

    Raises
    ------
    NIDUnavailableError
        The download failed and there is no cached copy. A stale cache is used,
        with a warning, rather than raising.
    """
    folder = Path(cache_dir) if cache_dir is not None else nid_cache_dir()
    data_path, meta_path = folder / NID_CACHE_FILE, folder / NID_META_FILE
    meta: Dict[str, Any] = {}
    if data_path.exists() and meta_path.exists():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        age = (date.today() - date.fromisoformat(meta["retrieved"])).days
        if not refresh and age <= max_age_days:
            return _read_nid_cache(data_path, meta)
    try:
        text = _download_nid(timeout, retries)
    except NIDUnavailableError:
        if meta:
            logger.warning("Using the stale NID cache retrieved %s", meta["retrieved"])
            return _read_nid_cache(data_path, meta)
        raise
    dams, updated = parse_nid_csv(text)
    folder.mkdir(parents=True, exist_ok=True)
    with gzip.GzipFile(data_path, "wb", mtime=0) as fh:
        fh.write(dams.to_csv(index=False).encode("utf-8"))
    meta = {
        "source": NID_CSV_URL,
        "data_last_updated": updated,
        "retrieved": datetime.now(timezone.utc).date().isoformat(),
        "rows": len(dams),
    }
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return NIDInventory(dams, updated, meta["retrieved"])


def _read_nid_cache(path: Path, meta: Dict[str, Any]) -> NIDInventory:
    dams = pd.read_csv(path, dtype={"nid_id": str})
    return NIDInventory(dams, str(meta["data_last_updated"]), str(meta["retrieved"]))


def points_in_geometry(
    lon: np.ndarray, lat: np.ndarray, geometry: Dict[str, Any], chunk: int = 4096
) -> np.ndarray:
    """Even-odd point-in-polygon for many points against a GeoJSON (Multi)Polygon.

    Ray casting in plain longitude/latitude, as
    :func:`flowfreq.streamstats._point_in_ring` does for one point: exact for
    the polygon as drawn in that plane, whose edges are tens of metres long at
    StreamStats' 10 m DEM resolution, so treating them as straight lines in
    lon/lat rather than geodesics moves the boundary by well under a metre.
    Holes are honoured.

    Parameters
    ----------
    lon, lat : array_like
        Point coordinates, decimal degrees.
    geometry : dict
        GeoJSON ``Polygon`` or ``MultiPolygon`` geometry (``[lon, lat]`` order).
    chunk : int
        Ring edges processed per vectorized step (bounds memory).

    Returns
    -------
    numpy.ndarray of bool
    """
    from .streamstats import _polygons_of

    x = np.asarray(lon, dtype=float)
    y = np.asarray(lat, dtype=float)
    out = np.zeros(x.shape, dtype=bool)

    def ring_mask(ring: List[List[float]]) -> np.ndarray:
        r = np.asarray([p[:2] for p in ring], dtype=float)
        inside = np.zeros(x.shape, dtype=bool)
        x1, y1, x2, y2 = r[:-1, 0], r[:-1, 1], r[1:, 0], r[1:, 1]
        for s in range(0, len(x1), chunk):
            a, b, c, d = x1[s : s + chunk], y1[s : s + chunk], x2[s : s + chunk], y2[s : s + chunk]
            straddle = (b[None, :] > y[:, None]) != (d[None, :] > y[:, None])
            with np.errstate(divide="ignore", invalid="ignore"):
                xcross = (c - a)[None, :] * (y[:, None] - b[None, :]) / (d - b)[None, :] + a[
                    None, :
                ]
            crossings = np.count_nonzero(straddle & (x[:, None] < xcross), axis=1)
            inside ^= (crossings % 2).astype(bool)
        return inside

    for polygon in _polygons_of(geometry):
        if not polygon:
            continue
        hit = ring_mask(polygon[0])
        for hole in polygon[1:]:
            hit &= ~ring_mask(hole)
        out |= hit
    return out


@dataclass(frozen=True)
class BasinDamStorage:
    """NID dams inside a watershed polygon and their summed storage.

    Attributes
    ----------
    n_dams : int
    storage_acft : float
        Sum of ``NID Storage`` (the NID's maximum storage, the field GAGES-II's
        ``STOR_NID_2009`` sums), acre-feet.
    area_km2 : float
        Geodesic area of the polygon.
    storage_ml_km2 : float
        ``storage_acft`` per unit area, megaliters per km^2 (= mm of depth), the
        unit of GAGES-II ``STOR_NID_2009``.
    nid_ids : tuple of str
        The dams counted, largest storage first.
    data_last_updated : str
        The NID release used.
    """

    n_dams: int
    storage_acft: float
    area_km2: float
    storage_ml_km2: float
    nid_ids: Tuple[str, ...]
    data_last_updated: str


def nid_storage_in_basin(geometry: Dict[str, Any], nid: NIDInventory) -> BasinDamStorage:
    """Sum NID storage over the dams whose coordinates fall inside a basin.

    A dam is counted when its NID point is inside the polygon
    (:func:`points_in_geometry`). NID points mark the dam, which sits on the
    stream, so a dam within a few tens of metres of the outlet can fall on
    either side of a boundary drawn from a 10 m DEM; GAGES-II's 2009 sums
    share that limitation.

    Parameters
    ----------
    geometry : dict
        GeoJSON ``Polygon``/``MultiPolygon`` of the watershed.
    nid : NIDInventory

    Returns
    -------
    BasinDamStorage
    """
    from .streamstats import _polygons_of, geojson_area_sq_mi

    ext = np.asarray(
        [p[:2] for poly in _polygons_of(geometry) if poly for p in poly[0]], dtype=float
    )
    dams = nid.dams
    box = (
        dams["longitude"].between(ext[:, 0].min(), ext[:, 0].max())
        & dams["latitude"].between(ext[:, 1].min(), ext[:, 1].max())
    ).to_numpy()
    cand = dams[box]
    inside = points_in_geometry(cand["longitude"].to_numpy(), cand["latitude"].to_numpy(), geometry)
    hit = cand[inside].sort_values("nid_storage_acft", ascending=False)
    area_km2 = geojson_area_sq_mi(geometry) * SQ_MI_TO_KM2
    storage = float(hit["nid_storage_acft"].sum())
    return BasinDamStorage(
        n_dams=len(hit),
        storage_acft=storage,
        area_km2=area_km2,
        storage_ml_km2=storage * ACRE_FT_TO_ML / area_km2 if area_km2 > 0 else math.nan,
        nid_ids=tuple(hit["nid_id"].astype(str)),
        data_last_updated=nid.data_last_updated,
    )


def runoff_mm_from_flow(mean_flow_cfs: float, area_km2: float) -> float:
    """Mean annual runoff depth (mm/yr, 365-day year) from a mean flow and an area.

    With this, ``storage_ml_km2 / runoff * 365`` is exactly storage divided by
    one day's mean-flow volume, Dudley and others' (2018) normalization.
    """
    return mean_flow_cfs * CFS_TO_M3S * 86400.0 * 365.0 / (area_km2 * 1.0e6) * 1000.0


def mean_flow_from_daily(site_no: str, *, min_days: int = 365) -> Tuple[float, int]:
    """Mean daily discharge (cfs) over a gage's complete water years.

    Dudley and others (2018) normalized storage by observed mean flow; this is
    that, over the whole daily record in water years with at least
    ``min_days`` daily values (they used 1966-2015). One Water Data API query
    (:func:`flowfreq.waterdata.download_daily`).

    Returns
    -------
    (float, int)
        Mean flow and the number of complete water years used.

    Raises
    ------
    ValueError
        No complete water year of daily discharge.
    """
    from .waterdata import download_daily

    daily = download_daily(site_no)
    q = pd.to_numeric(daily["flow_cfs"], errors="coerce").dropna()
    if q.empty:
        raise ValueError(f"USGS {site_no} has no daily discharge")
    idx = pd.DatetimeIndex(q.index)
    wy = idx.year + (idx.month >= 10)
    counts = q.groupby(wy).count()
    full = counts[counts >= min_days].index
    q = q[np.isin(wy, full)]
    if q.empty:
        raise ValueError(f"USGS {site_no} has no complete water year of daily discharge")
    return float(q.mean()), len(full)


# ------------------------------------------------- NLCD impervious (current) --

#: StreamStats basin characteristics that are NLCD percent impervious, and the
#: NLCD epoch of each. Verified live 2026-10-03 against every region's
#: ``ss-hydro/v1/basin-characteristics/{region}`` list and by delineations that
#: computed them (GA: LC06IMP 9.33 / LC11IMP 10.56 / LC19IMP 12.94 at 02337000,
#: where GAGES-II's IMPNLCD06 is 9.41; CO, UT LC11IMP; WY LC16IMP; NV LC21IMP;
#: ID IMPNLCD01). Which codes a region offers varies: e.g. GA, NC, SC and VA
#: have three or more epochs, most states only LC11IMP, and WA, MT, NM, IL,
#: WI, NJ, VT, ME, LA, WV and AK none. ``IMPERV`` (DE, DC, MD, OR) is left out
#: because it names no NLCD epoch.
NLCD_IMPERVIOUS_CODES: Dict[str, int] = {
    "IMPNLCD01": 2001,
    "LC01IMP": 2001,
    "LC06IMP": 2006,
    "LC11IMP": 2011,
    "LC16IMP": 2016,
    "LC19IMP": 2019,
    "LC21IMP": 2021,
    "LC23IMP": 2023,
}

#: Regions whose impervious characteristic failed the live cross-check
#: against GAGES-II ``IMPNLCD06`` (2026-10-03), so it is recorded but not used.
#: Elsewhere it agreed: GA LC06IMP 9.33 vs 9.41 (02337000); TN LC11IMP 35.97 vs
#: 32.35 (07031692) and 25.57 vs 25.07 (03431300).
IMPERVIOUS_UNVERIFIED_REGIONS: Dict[str, str] = {
    "CO": (
        "its LC11IMP is about 2.6 times GAGES-II's NLCD 2006 impervious at both "
        "GAGES-II gages checked, 54.9 vs 21.2 percent at 06710150 and 39.3 vs 14.7 at "
        "07105600, more than five years of land-cover change explains"
    ),
}


def impervious_series(watershed: "WatershedCharacteristics") -> Dict[int, float]:
    """NLCD percent impervious by epoch year from a StreamStats result.

    Parameters
    ----------
    watershed : flowfreq.streamstats.WatershedCharacteristics
        Requested with the :data:`NLCD_IMPERVIOUS_CODES` it should report.

    Returns
    -------
    dict of int to float
        Year to percent impervious, ascending; empty when the region computes none.
    """
    out: Dict[int, float] = {}
    for code, year in NLCD_IMPERVIOUS_CODES.items():
        if code in watershed.characteristics:
            out[year] = float(watershed.characteristics[code].value)
    return dict(sorted(out.items()))


# ------------------------------------------------------ current-data refine --


def _locate(site_no: str) -> Tuple[str, float, float]:
    """(StreamStats region, lat, lon): the packaged catalog, else the Water Data API."""
    from .catalog import load_catalog

    try:
        cat = load_catalog()
        row = cat[cat["site_no"] == site_no]
    except (OSError, ValueError):
        row = pd.DataFrame()
    if len(row) and pd.notna(row.iloc[0]["latitude"]):
        r = row.iloc[0]
        return str(r["state"]), float(r["latitude"]), float(r["longitude"])
    from .regional_skew import STATE_FIPS
    from .waterdata import fetch_monitoring_location_feature, site_attributes

    attrs = site_attributes(fetch_monitoring_location_feature(site_no))
    state = STATE_FIPS.get(str(attrs.get("state_code") or "").zfill(2))
    if state is None or attrs["latitude"] is None or attrs["longitude"] is None:
        raise ValueError(f"USGS {site_no}: no state or coordinates to delineate from")
    return state, float(attrs["latitude"]), float(attrs["longitude"])


def current_basin_evidence(
    site_no: str,
    *,
    use_current_nid: bool = True,
    use_current_impervious: bool = True,
    watershed: Optional["WatershedCharacteristics"] = None,
    region: Optional[str] = None,
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    runoff_mm: Optional[float] = None,
    mean_flow_cfs: Optional[float] = None,
    nid: Optional[NIDInventory] = None,
    streamstats_cache: Optional["StreamStatsCache"] = None,
) -> Dict[str, Any]:
    """Current dam storage and impervious cover for one gage's basin.

    Parameters
    ----------
    site_no : str
    use_current_nid, use_current_impervious : bool
        Which pieces to compute.
    watershed : WatershedCharacteristics, optional
        A StreamStats result with its polygon and, for impervious cover, the
        :data:`NLCD_IMPERVIOUS_CODES`. Delineated here when absent.
    region, lat, lon : optional
        Pour point for that delineation; defaults to the gage's location
        (packaged catalog, else the Water Data API) with its state as region.
    runoff_mm : float, optional
        Mean annual runoff, mm/yr, to normalize storage by (GAGES-II
        ``RUNAVE7100`` for a GAGES-II gage).
    mean_flow_cfs : float, optional
        Used when ``runoff_mm`` is absent; otherwise the mean of the gage's
        complete water years of daily flow (:func:`mean_flow_from_daily`).
    nid : NIDInventory, optional
        Defaults to :func:`load_nid`.
    streamstats_cache : StreamStatsCache, optional
        Passed to the delineation.

    Returns
    -------
    dict
        Keys present only when computed: ``streamstats_regulated`` (the
        StreamStats warning text), ``nid_*`` (``n_dams``, ``storage_acft``,
        ``storage_ml_km2``, ``area_km2``, ``runoff_mm``, ``runoff_source``,
        ``norm_storage_days``, ``data_last_updated``, ``retrieved``),
        ``imperv_series`` (year to percent), ``imperv_pct``, ``imperv_year``,
        ``imperv_source``; and ``unavailable``, a list of reasons a requested
        piece could not be computed. Never raises for a service failure.
    """
    from . import streamstats as ss

    ev: Dict[str, Any] = {"unavailable": []}
    unavailable: List[str] = ev["unavailable"]
    if watershed is None:
        try:
            if region is None or lat is None or lon is None:
                r0, la0, lo0 = _locate(site_no)
                region = region or r0
                lat = la0 if lat is None else lat
                lon = lo0 if lon is None else lon
            watershed = ss.delineate_and_get_characteristics(
                region,
                lat,
                lon,
                characteristic_codes=["DRNAREA", *NLCD_IMPERVIOUS_CODES],
                cache=streamstats_cache,
            )
        except ss.DegenerateDelineationError as exc:
            if "regulated" in str(exc).lower():
                ev["streamstats_regulated"] = str(exc)
            else:
                unavailable.append(f"StreamStats delineation failed: {exc}")
            return ev
        except (ValueError, requests.RequestException) as exc:
            unavailable.append(f"StreamStats delineation failed: {exc}")
            return ev
    ev["streamstats_region"] = watershed.region

    geometry = (watershed.polygon_geojson or {}).get("geometry")
    if use_current_nid:
        if geometry is None:
            unavailable.append("NID: the StreamStats result carries no watershed polygon")
        else:
            try:
                inv = nid if nid is not None else load_nid()
            except NIDUnavailableError as exc:
                unavailable.append(f"NID: {exc}")
                inv = None
            if inv is not None:
                st = nid_storage_in_basin(geometry, inv)
                ev.update(
                    nid_n_dams=st.n_dams,
                    nid_storage_acft=round(st.storage_acft, 1),
                    nid_storage_ml_km2=round(st.storage_ml_km2, 2),
                    nid_area_km2=round(st.area_km2, 2),
                    nid_data_last_updated=inv.data_last_updated,
                    nid_retrieved=inv.retrieved,
                    nid_ids=list(st.nid_ids),
                )
                source = "GAGES-II RUNAVE7100"
                if runoff_mm is None:
                    try:
                        if mean_flow_cfs is None:
                            mean_flow_cfs, nyears = mean_flow_from_daily(site_no)
                            source = f"mean daily flow, {nyears} complete water years"
                        else:
                            source = "mean flow supplied by caller"
                        runoff_mm = runoff_mm_from_flow(mean_flow_cfs, st.area_km2)
                    except (ValueError, requests.RequestException) as exc:
                        unavailable.append(f"NID: no mean runoff to normalize by ({exc})")
                if runoff_mm is not None and runoff_mm > 0:
                    ev["nid_runoff_mm"] = round(float(runoff_mm), 1)
                    ev["nid_runoff_source"] = source
                    ev["nid_norm_storage_days"] = round(
                        st.storage_ml_km2 / float(runoff_mm) * 365.0, 1
                    )
    if use_current_impervious:
        series = impervious_series(watershed)
        suspect = IMPERVIOUS_UNVERIFIED_REGIONS.get(watershed.region)
        if series and suspect:
            ev["imperv_series_unverified"] = {str(k): v for k, v in series.items()}
            unavailable.append(
                f"impervious: StreamStats region {watershed.region} not used ({suspect})"
            )
        elif series:
            year = max(series)
            ev.update(
                imperv_series={str(k): v for k, v in series.items()},
                imperv_pct=round(series[year], 2),
                imperv_year=year,
                imperv_source=f"NLCD {year} via StreamStats region {watershed.region}",
            )
        else:
            unavailable.append(
                f"impervious: StreamStats region {watershed.region} computes no NLCD "
                "impervious characteristic"
            )
    return ev


def _attribute_rule(
    storage_days: Optional[float], imperv: Optional[float]
) -> Tuple[SiteClass, Optional[str]]:
    if storage_days is not None and storage_days > REGULATED_STORAGE_DAYS:
        return SiteClass.REGULATED, "storage"
    if imperv is not None and imperv > URBAN_IMPERVIOUS_PCT:
        return SiteClass.URBAN, "impervious"
    return SiteClass.UNKNOWN, None


def classify_site(
    site_no: str,
    peak_codes: Optional[Iterable[object]] = None,
    *,
    table: Optional[pd.DataFrame] = None,
    use_current_nid: bool = False,
    use_current_impervious: bool = False,
    watershed: Optional["WatershedCharacteristics"] = None,
    region: Optional[str] = None,
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    mean_flow_cfs: Optional[float] = None,
    nid: Optional[NIDInventory] = None,
    streamstats_cache: Optional["StreamStatsCache"] = None,
) -> SiteClassification:
    """Classify a gage from the screen table, current basin data and peak codes.

    By default this is offline: the GAGES-II table plus peak codes.

    With ``use_current_nid`` and/or ``use_current_impervious`` the attribute
    rules are re-evaluated on current data (:func:`current_basin_evidence`):
    today's NID storage inside the StreamStats watershed instead of NID 2009,
    and the newest NLCD impervious epoch StreamStats has for the region
    instead of NLCD 2006. A current value replaces its GAGES-II counterpart; a
    piece that cannot be computed falls back to it. The rules and thresholds
    are unchanged (first match wins):

    1. GAGES-II ``Ref`` stays ``reference`` (Dudley and others 2018 let their
       reference class stand whatever the storage); a current value over a
       threshold is noted in ``basis``.
    2. StreamStats' "River is regulated" delineation warning -> ``regulated``.
    3. Normalized storage > :data:`REGULATED_STORAGE_DAYS` -> ``regulated``.
    4. Impervious > :data:`URBAN_IMPERVIOUS_PCT` -> ``urban``.
    5. Otherwise ``unknown`` (a gage outside GAGES-II can never be shown to
       be ``reference`` this way).

    Peak code 6 makes the record ``regulated`` in every case.

    Parameters
    ----------
    site_no : str
        USGS site number.
    peak_codes : iterable, optional
        One raw NWIS ``peak_cd`` / ``qualification_code`` field per peak.
    table : pandas.DataFrame, optional
        A screen table (defaults to :func:`load_screen`).
    use_current_nid : bool, default False
        Opt in to current NID storage. Downloads the NID once (cached, see
        :func:`load_nid`) and delineates the basin with StreamStats.
    use_current_impervious : bool, default False
        Opt in to current NLCD impervious cover from StreamStats.
    watershed, region, lat, lon, mean_flow_cfs, nid, streamstats_cache
        Passed to :func:`current_basin_evidence`.

    Returns
    -------
    SiteClassification
        ``unknown`` with an empty ``attributes`` when the gage is not in the
        table and nothing else decides it.
    """
    site_no = str(site_no).strip()
    tab = table if table is not None else _load(str(_default_path()))
    attributes: Dict[str, Any] = {}
    basis = []
    site_class = SiteClass.UNKNOWN
    if site_no in tab.index:
        row = tab.loc[site_no]
        attributes = {c: _clean(row[c]) for c in SCREEN_COLUMNS if c != "site_no"}
        site_class = SiteClass(row["site_class"])
        basis.append(f"{row['basis']} ({row['source']} {row['source_date']})")
    else:
        basis.append("not in the GAGES-II screen table")

    current: Dict[str, Any] = {}
    if use_current_nid or use_current_impervious:
        current = current_basin_evidence(
            site_no,
            use_current_nid=use_current_nid,
            use_current_impervious=use_current_impervious,
            watershed=watershed,
            region=region,
            lat=lat,
            lon=lon,
            runoff_mm=attributes.get("runave7100_mm"),
            mean_flow_cfs=mean_flow_cfs,
            nid=nid,
            streamstats_cache=streamstats_cache,
        )
        site_class = _refine(site_class, attributes, current, basis)

    n6 = n5 = n = 0
    if peak_codes is not None:
        for raw in peak_codes:
            codes = parse_codes(raw)
            n += 1
            n6 += REGULATED_PEAK_CODE in codes
            n5 += UNKNOWN_DEGREE_PEAK_CODE in codes
        if n6:
            site_class = SiteClass.REGULATED
            basis.append(f"peak code 6 (regulation or diversion) on {n6} of {n} peaks")
        if n5:
            basis.append(f"peak code 5 (unknown degree) on {n5} of {n} peaks; not decisive")
    return SiteClassification(
        site_no=site_no,
        site_class=site_class,
        basis=tuple(basis),
        attributes=attributes,
        n_peaks_code6=n6,
        n_peaks_code5=n5,
        n_peaks=n,
        current=current,
    )


def _refine(
    site_class: SiteClass,
    attributes: Dict[str, Any],
    current: Dict[str, Any],
    basis: List[str],
) -> SiteClass:
    """Apply :func:`classify_site`'s rules 1-5 with current values; extend ``basis``."""
    for reason in current.get("unavailable", []):
        basis.append(f"current data: {reason}")
    days = current.get("nid_norm_storage_days")
    imperv = current.get("imperv_pct")
    if days is not None:
        basis.append(
            f"current NID ({current['nid_data_last_updated']}): {current['nid_n_dams']} dams, "
            f"storage {days:.1f} d of mean runoff ({current['nid_runoff_source']})"
        )
    if imperv is not None:
        basis.append(f"impervious {imperv:.2f}% ({current['imperv_source']})")
    if "streamstats_regulated" in current:
        basis.append(f"StreamStats: {current['streamstats_regulated']}")

    eff_days = days if days is not None else attributes.get("norm_storage_days")
    eff_imperv = imperv if imperv is not None else attributes.get("imperv_pct_2006")
    if site_class is SiteClass.REFERENCE:
        cls, why = _attribute_rule(eff_days, eff_imperv)
        if why or "streamstats_regulated" in current:
            basis.append(
                f"current {why or 'StreamStats warning'} would class it {cls.value}; "
                "GAGES-II Ref class retained"
            )
        return site_class
    if "streamstats_regulated" in current:
        return SiteClass.REGULATED
    if days is None and imperv is None:
        return site_class  # nothing current: the table's class (or unknown) stands
    cls, why = _attribute_rule(eff_days, eff_imperv)
    if why == "storage":
        basis.append(f"storage {eff_days:.1f} d > {REGULATED_STORAGE_DAYS} d")
    elif why == "impervious":
        basis.append(f"impervious {eff_imperv:.2f}% > {URBAN_IMPERVIOUS_PCT:g}%")
    return cls


def require_unregulated(
    classification: SiteClassification,
    *,
    allow_regulated: bool = False,
) -> Dict[str, Any]:
    """Refuse Bulletin 17C on a regulated record unless explicitly overridden.

    Parameters
    ----------
    classification : SiteClassification
    allow_regulated : bool, default False
        Fit anyway. Logged as a warning and recorded in the provenance.

    Returns
    -------
    dict
        Provenance: :meth:`SiteClassification.as_dict` plus
        ``regulation_override`` (True only when a regulated record was fitted
        because of ``allow_regulated``).

    Raises
    ------
    RegulatedRecordError
        If the record is ``regulated`` and ``allow_regulated`` is False.
    """
    prov = classification.as_dict()
    prov["regulation_override"] = False
    cls = classification.site_class
    if cls is SiteClass.REGULATED:
        reasons = "; ".join(classification.basis)
        if not allow_regulated:
            raise RegulatedRecordError(
                f"USGS {classification.site_no} is classified regulated ({reasons}). "
                "Bulletin 17C assumes unregulated flows and gives no regulated-flow method "
                "(England and others 2019, p. 36). Pass allow_regulated=True to fit it anyway; "
                "the override is recorded in the result."
            )
        logger.warning(
            "USGS %s is classified regulated (%s); fitting B17C because allow_regulated=True",
            classification.site_no,
            reasons,
        )
        prov["regulation_override"] = True
    elif (
        not classification.in_table
        and classification.n_peaks_code6 == 0
        and not classification.has_current_evidence
    ):
        logger.warning(
            "USGS %s has no regulation classification (not in GAGES-II); proceeding "
            "without the regulation screen",
            classification.site_no,
        )
    elif cls is SiteClass.URBAN:
        logger.info(
            "USGS %s is classified urban (%s)", classification.site_no, classification.basis
        )
    return prov
