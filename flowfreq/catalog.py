"""National gage catalog: schema, loader, and the pieces the build tool uses.

The catalog lists every USGS peak-flow site with at least 10 water years of
discharge peaks, with the attributes the rest of the roadmap selects on:
name, location, drainage area, state, HUC8, years of record, regulation
class, and (Wave 1 states) regression region. It is rebuilt by
``tools/build_gage_catalog.py`` from the Water Data OGC API ``peaks`` and
``monitoring-locations`` collections and is never edited by hand. The packaged
file is ``data/gage_catalog.csv.gz``, with its build metadata in
``data/gage_catalog.meta.json``.

Columns:

- ``regulation_class``: :func:`flowfreq.regulation.classify_site` on the
  site's peak codes plus the packaged GAGES-II screen: ``reference``,
  ``regulated``, ``urban`` or ``unknown``.
- ``regression_region``: the NSS region code(s) (``GC1750``) of the
  StreamStats peak-flow region polygon containing the **gage location**,
  ``;``-joined when the polygon names several (western Oregon's 2A/2B split
  is by basin mean elevation, not location). Filled for WA, OR, ID and MT
  only and blank elsewhere. A basin that crosses a region boundary needs the
  area-weighted lookup, :func:`flowfreq.streamstats.locate_regression_regions`.

The 3-row ``gage_attributes.csv`` seed that :class:`flowfreq.usgs.GageAttributes`
reads is unchanged.

Roadmap: ``docs/MASTER_ROADMAP.md`` §1.1, issue #33.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Union

import pandas as pd

from flowfreq.regression.jurisdictions import BY_CODE

SEED_PATH = Path(__file__).parent / "data" / "gage_attributes.csv"
CATALOG_PATH = Path(__file__).parent / "data" / "gage_catalog.csv.gz"
META_PATH = Path(__file__).parent / "data" / "gage_catalog.meta.json"

#: The catalog's inclusion rule: water years with a discharge peak.
MIN_YEARS = 10

REQUIRED_COLUMNS = ("site_no", "site_name", "state")
OPTIONAL_COLUMNS = (
    "latitude",
    "longitude",
    "drainage_area_sqmi",
    "huc8",
    "n_peaks",
    "first_water_year",
    "last_water_year",
    "regulation_class",
    "regression_region",
)
COLUMNS = REQUIRED_COLUMNS + OPTIONAL_COLUMNS
REGULATION_CLASSES = frozenset({"reference", "regulated", "urban", "unknown"})


def load_catalog(path: Union[str, Path, None] = None) -> pd.DataFrame:
    """Load and validate the gage catalog.

    Parameters
    ----------
    path : str or Path, optional
        Defaults to the built catalog if present, otherwise the seed file.
        ``.csv`` and ``.csv.gz`` are both read.

    Returns
    -------
    pandas.DataFrame
        All :data:`COLUMNS`, with ``site_no`` and ``huc8`` as strings.

    Raises
    ------
    ValueError
        On missing required columns, duplicate or non-numeric site numbers,
        unknown state codes, non-positive drainage areas, a malformed HUC8,
        or an unknown regulation class.
    """
    if path is None:
        path = CATALOG_PATH if CATALOG_PATH.exists() else SEED_PATH
    df = pd.read_csv(
        path,
        dtype={"site_no": str, "huc8": str, "state": str, "regression_region": str},
        keep_default_na=True,
    )
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"catalog missing required columns: {missing}")
    for col in OPTIONAL_COLUMNS:
        if col not in df.columns:
            df[col] = pd.NA
    df["site_no"] = df["site_no"].str.strip()
    if not df["site_no"].str.fullmatch(r"\d{8,15}").all():
        raise ValueError("site_no must be 8-15 digits")
    if df["site_no"].duplicated().any():
        raise ValueError(
            f"duplicate site_no: {sorted(df.loc[df['site_no'].duplicated(), 'site_no'])}"
        )
    unknown = sorted(set(df["state"].str.upper()) - BY_CODE.keys())
    if unknown:
        raise ValueError(f"unknown state code(s): {unknown}")
    da = pd.to_numeric(df["drainage_area_sqmi"], errors="coerce")
    if (da <= 0).any():
        raise ValueError("drainage_area_sqmi must be positive")
    huc = df["huc8"].dropna().astype(str)
    if not huc.str.fullmatch(r"\d{8}").all():
        raise ValueError("huc8 must be 8 digits")
    cls = set(df["regulation_class"].dropna()) - REGULATION_CLASSES
    if cls:
        raise ValueError(f"unknown regulation_class value(s): {sorted(cls)}")
    return df.loc[:, list(COLUMNS)]


def load_catalog_meta(path: Union[str, Path, None] = None) -> Dict[str, Any]:
    """Build metadata of the packaged catalog (sources, date, counts)."""
    return dict(json.loads(Path(path or META_PATH).read_text(encoding="utf-8")))


def select(
    catalog: pd.DataFrame,
    *,
    state: Optional[str] = None,
    min_years: Optional[int] = None,
    exclude_regulated: bool = False,
) -> pd.DataFrame:
    """Filter the catalog.

    Parameters
    ----------
    catalog : pandas.DataFrame
        Output of :func:`load_catalog`.
    state : str, optional
        Two-letter postal code.
    min_years : int, optional
        Minimum ``n_peaks``. Rows with unknown ``n_peaks`` are dropped.
    exclude_regulated : bool, default False
        Drop rows whose ``regulation_class`` is ``"regulated"``.

    Returns
    -------
    pandas.DataFrame
    """
    out = catalog
    if state is not None:
        out = out[out["state"].str.upper() == state.upper()]
    if min_years is not None:
        out = out[pd.to_numeric(out["n_peaks"], errors="coerce") >= min_years]
    if exclude_regulated:
        out = out[out["regulation_class"] != "regulated"]
    return out.reset_index(drop=True)


def catalog_row(
    site_no: str,
    state: str,
    peaks: pd.DataFrame,
    *,
    site_name: str = "",
    latitude: Optional[float] = None,
    longitude: Optional[float] = None,
    drainage_area_sqmi: Optional[float] = None,
    huc8: Optional[str] = None,
    regression_region: Optional[str] = None,
) -> dict:
    """Build one catalog row from a site's peaks and metadata.

    Parameters
    ----------
    site_no, state : str
    peaks : pandas.DataFrame
        A peak frame in the :mod:`flowfreq.peak_sources` contract
        (``water_year``, ``peak_flow_cfs``, ``qualification_code``). Only
        rows with a discharge count toward the years of record.
    site_name, latitude, longitude, drainage_area_sqmi, huc8, regression_region : optional
        Site metadata.

    Returns
    -------
    dict
        Keyed by :data:`COLUMNS`. ``regulation_class`` is
        :func:`flowfreq.regulation.classify_site` of the site and its codes.
    """
    from flowfreq.regulation import classify_site

    flows = pd.to_numeric(peaks["peak_flow_cfs"], errors="coerce")
    with_q = peaks.loc[flows.notna()]
    years = pd.to_numeric(with_q["water_year"], errors="coerce").dropna().astype(int)
    codes = with_q["qualification_code"].tolist() if "qualification_code" in with_q else None
    reg = classify_site(site_no, codes)
    return {
        "site_no": site_no,
        "site_name": site_name,
        "state": state.upper(),
        "latitude": latitude,
        "longitude": longitude,
        "drainage_area_sqmi": drainage_area_sqmi,
        "huc8": huc8,
        "n_peaks": int(years.nunique()),
        "first_water_year": int(years.min()) if len(years) else None,
        "last_water_year": int(years.max()) if len(years) else None,
        "regulation_class": reg.site_class.value,
        "regression_region": regression_region,
    }


# ----------------------------------------------------------------------------
# Bulk Water Data API features -> catalog rows (used by the build tool)
# ----------------------------------------------------------------------------


def _strip_prefix(location_id: str) -> str:
    s = str(location_id).strip()
    return s[5:] if s.upper().startswith("USGS-") else s


def _num(value: Any) -> float:
    """A float, or NaN for None / non-numeric text (the API sends strings)."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def peaks_by_site(features: Iterable[Mapping[str, Any]]) -> Dict[str, pd.DataFrame]:
    """Group Water Data API ``peaks`` features into per-site peak frames.

    Parameters
    ----------
    features : iterable of dict
        GeoJSON features with ``monitoring_location_id``, ``water_year``,
        ``value`` and ``qualifier`` properties (other properties are ignored).
        Non-USGS agencies' sites are skipped.

    Returns
    -------
    dict of str to pandas.DataFrame
        Site number to a frame with ``water_year``, ``peak_flow_cfs`` and
        ``qualification_code`` (NWIS code string, via
        :func:`flowfreq.peak_sources.qualifiers_to_codes`).
    """
    from flowfreq.peak_sources import qualifiers_to_codes

    rows: Dict[str, List[Dict[str, Any]]] = {}
    for feat in features:
        p = feat.get("properties") or {}
        loc = str(p.get("monitoring_location_id") or "")
        if not loc.upper().startswith("USGS-"):
            continue
        site = _strip_prefix(loc)
        wy = p.get("water_year")
        rows.setdefault(site, []).append(
            {
                "water_year": int(wy) if wy is not None else None,
                "peak_flow_cfs": _num(p.get("value")),
                "qualification_code": qualifiers_to_codes(p.get("qualifier"), site, wy),
            }
        )
    return {s: pd.DataFrame(r) for s, r in rows.items()}


def location_fields(feature: Mapping[str, Any]) -> Dict[str, Any]:
    """Catalog metadata from one ``monitoring-locations`` feature.

    Returns
    -------
    dict
        ``site_no``, ``site_name``, ``state`` (postal code, or None),
        ``latitude``, ``longitude``, ``drainage_area_sqmi`` (None unless
        positive) and ``huc8`` (first 8 digits of ``hydrologic_unit_code``).
    """
    from flowfreq.regional_skew import STATE_FIPS

    p = feature.get("properties") or {}
    geom = feature.get("geometry") or {}
    coords = geom.get("coordinates") if geom.get("type") == "Point" else None
    huc = str(p.get("hydrologic_unit_code") or "")
    da = _num(p.get("drainage_area"))
    return {
        "site_no": _strip_prefix(str(feature.get("id") or p.get("monitoring_location_number"))),
        "site_name": str(p.get("monitoring_location_name") or "").strip(),
        "state": STATE_FIPS.get(str(p.get("state_code") or "").zfill(2)),
        "latitude": float(coords[1]) if coords else None,
        "longitude": float(coords[0]) if coords else None,
        "drainage_area_sqmi": float(da) if pd.notna(da) and da > 0 else None,
        "huc8": huc[:8] if len(huc) >= 8 and huc[:8].isdigit() else None,
    }


def build_rows(
    peaks: Mapping[str, pd.DataFrame],
    locations: Mapping[str, Mapping[str, Any]],
    *,
    min_years: int = MIN_YEARS,
    regions: Optional[Mapping[str, str]] = None,
    default_state: Optional[str] = None,
) -> pd.DataFrame:
    """Assemble catalog rows for sites with at least ``min_years`` of peaks.

    Parameters
    ----------
    peaks : mapping
        Site number to peak frame (:func:`peaks_by_site`).
    locations : mapping
        Site number to :func:`location_fields` output. A site with peaks but
        no location record is kept with ``default_state`` and null metadata.
    min_years : int, default :data:`MIN_YEARS`
    regions : mapping, optional
        Site number to regression region code(s).
    default_state : str, optional
        Postal code for sites whose location lacks one.

    Returns
    -------
    pandas.DataFrame
        :data:`COLUMNS`, sorted by state and site.
    """
    out: List[dict] = []
    for site, frame in peaks.items():
        loc = dict(locations.get(site) or {})
        state = loc.get("state") or default_state
        if not state:
            continue
        row = catalog_row(
            site,
            state,
            frame,
            site_name=loc.get("site_name", ""),
            latitude=loc.get("latitude"),
            longitude=loc.get("longitude"),
            drainage_area_sqmi=loc.get("drainage_area_sqmi"),
            huc8=loc.get("huc8"),
            regression_region=(regions or {}).get(site),
        )
        if row["n_peaks"] >= min_years:
            out.append(row)
    df = pd.DataFrame(out, columns=list(COLUMNS))
    if len(df):
        df = df.sort_values(["state", "site_no"]).reset_index(drop=True)
    return df


def write_catalog(df: pd.DataFrame, path: Union[str, Path]) -> None:
    """Write the catalog deterministically (gzip without a timestamp for ``.gz``)."""
    p = Path(path)
    df = df.loc[:, list(COLUMNS)].copy()
    for col in ("latitude", "longitude"):
        df[col] = pd.to_numeric(df[col], errors="coerce").round(6)
    compression: Optional[Dict[str, Any]] = (
        {"method": "gzip", "mtime": 0, "compresslevel": 9} if p.suffix == ".gz" else None
    )
    df.to_csv(p, index=False, compression=compression)


def sites_needing_location(
    peaks: Mapping[str, pd.DataFrame], locations: Mapping[str, Any], min_years: int = MIN_YEARS
) -> Sequence[str]:
    """Sites with enough peaks that the bulk location query did not return."""
    need = []
    for site, frame in peaks.items():
        n = pd.to_numeric(frame["water_year"], errors="coerce")[
            pd.to_numeric(frame["peak_flow_cfs"], errors="coerce").notna()
        ].nunique()
        if n >= min_years and site not in locations:
            need.append(site)
    return sorted(need)
