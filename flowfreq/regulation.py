"""Regulation / urbanization screen per gage, and the B17C refusal it drives.

Roadmap ``docs/MASTER_ROADMAP.md`` §1.1, issue #32. A gage is classed
``reference``, ``regulated``, ``urban`` or ``unknown`` from two kinds of
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

Bulletin 17C itself sets no numeric regulation threshold (England and others,
2019, p. 36, "Regulated Flow Frequency"); it says regulated-flow methods need
national guidance it does not give. That is why a regulated record is refused
rather than fitted: :func:`require_unregulated` raises
:class:`RegulatedRecordError` unless the caller passes an explicit override,
and the override is recorded in the returned provenance. A gage with no
classification at all proceeds, with a warning.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from enum import Enum
from functools import lru_cache
from importlib import resources
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple, Union

import pandas as pd

from .peak_codes import parse_codes

logger = logging.getLogger(__name__)

#: Peak code that, on any peak, makes a record regulated (SIR 2016-5118 p. 20).
REGULATED_PEAK_CODE = "6"
#: Recorded as evidence only; SIR 2016-5118 p. 20 treats it as unregulated.
UNKNOWN_DEGREE_PEAK_CODE = "5"

SCREEN_FILE = "regulation_screen.csv.gz"
SCREEN_COLUMNS = (
    "site_no",
    "site_class",
    "gagesii_class",
    "ndams_2009",
    "stor_nid_2009_ml_km2",
    "norm_storage_days",
    "imperv_pct_2006",
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


@dataclass(frozen=True)
class SiteClassification:
    """A gage's screen result and the evidence behind it.

    Attributes
    ----------
    site_no : str
    site_class : SiteClass
    basis : tuple of str
        Human-readable reasons, attribute evidence first, then peak codes.
    attributes : dict
        The screen-table row (empty when the gage is not in GAGES-II).
    n_peaks_code6, n_peaks_code5, n_peaks : int
        Peak-code counts, when peak codes were supplied.
    """

    site_no: str
    site_class: SiteClass
    basis: Tuple[str, ...] = ()
    attributes: Dict[str, Any] = field(default_factory=dict)
    n_peaks_code6: int = 0
    n_peaks_code5: int = 0
    n_peaks: int = 0

    @property
    def in_table(self) -> bool:
        """Whether basin-attribute evidence was available."""
        return bool(self.attributes)

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


def classify_site(
    site_no: str,
    peak_codes: Optional[Iterable[object]] = None,
    *,
    table: Optional[pd.DataFrame] = None,
) -> SiteClassification:
    """Classify a gage from the screen table and, optionally, its peak codes.

    Parameters
    ----------
    site_no : str
        USGS site number.
    peak_codes : iterable, optional
        One raw NWIS ``peak_cd`` / ``qualification_code`` field per peak. Any
        code 6 makes the record ``regulated`` whatever the table says.
    table : pandas.DataFrame, optional
        A screen table (defaults to :func:`load_screen`).

    Returns
    -------
    SiteClassification
        ``unknown`` with an empty ``attributes`` when the gage is not in the
        table and no peak code decides it.
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
    )


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
    elif not classification.in_table and classification.n_peaks_code6 == 0:
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
