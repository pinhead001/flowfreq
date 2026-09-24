"""National gage catalog: schema and loader.

The catalog lists the peak-flow gages available for analysis, with the
attributes the rest of the roadmap selects on: state, drainage area, years
of record, regulation class, and regression region. It is rebuilt by
``tools/build_gage_catalog.py`` and is never edited by hand.

Until that tool has been run with network access, the packaged seed is the
existing ``gage_attributes.csv`` (3 rows). The loader accepts it, and fills
the columns it lacks with nulls, so :class:`flowfreq.usgs.GageAttributes`
keeps working unchanged.

Roadmap: ``docs/MASTER_ROADMAP.md`` §1.1, issue #33.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import pandas as pd

from flowfreq.regression.jurisdictions import BY_CODE

SEED_PATH = Path(__file__).parent / "data" / "gage_attributes.csv"
CATALOG_PATH = Path(__file__).parent / "data" / "gage_catalog.csv"

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


def load_catalog(path: Union[str, Path, None] = None) -> pd.DataFrame:
    """Load and validate the gage catalog.

    Parameters
    ----------
    path : str or Path, optional
        Defaults to the built catalog if present, otherwise the seed file.

    Returns
    -------
    pandas.DataFrame
        All :data:`COLUMNS`, with ``site_no`` as a zero-padded string.

    Raises
    ------
    ValueError
        On missing required columns, duplicate or non-numeric site numbers,
        unknown state codes, or non-positive drainage areas.
    """
    if path is None:
        path = CATALOG_PATH if CATALOG_PATH.exists() else SEED_PATH
    df = pd.read_csv(path, dtype={"site_no": str, "huc8": str, "state": str})
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
    return df.loc[:, list(COLUMNS)]


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
) -> dict:
    """Build one catalog row from a site's peaks and metadata.

    The regulation class comes from peak-code evidence only
    (:func:`flowfreq.peak_codes.classify_from_codes`). ``regression_region``
    is left null until the Phase C region polygons exist.

    Parameters
    ----------
    site_no, state : str
    peaks : pandas.DataFrame
        A peak frame in the :mod:`flowfreq.peak_sources` contract.
    site_name, latitude, longitude, drainage_area_sqmi, huc8 : optional
        Site metadata.

    Returns
    -------
    dict
        Keyed by :data:`COLUMNS`.
    """
    from flowfreq.peak_codes import classify_from_codes

    years = pd.to_numeric(peaks["water_year"], errors="coerce").dropna()
    reg_class, _ = classify_from_codes(peaks["qualification_code"].tolist())
    return {
        "site_no": site_no,
        "site_name": site_name,
        "state": state.upper(),
        "latitude": latitude,
        "longitude": longitude,
        "drainage_area_sqmi": drainage_area_sqmi,
        "huc8": huc8,
        "n_peaks": int(len(years)),
        "first_water_year": int(years.min()) if len(years) else None,
        "last_water_year": int(years.max()) if len(years) else None,
        "regulation_class": reg_class.value,
        "regression_region": None,
    }
