"""NWIS peak-flow qualification codes and their Bulletin 17C treatment.

The treatment rules are peakfq 8.1.0's, not a reinterpretation of them:
``vendor/peakfqr/R/readInputs.R::siteQT`` is the authority, and this module
reproduces its interval construction for the codes it acts on (4, 8, 3, O, 6,
C). Every other code is carried through as information only, which is also
what ``siteQT`` does with it.

``siteQT`` matches codes by substring (``grepl("4", peak_cd)``), so a
multi-code string such as ``"4,C"`` gets both treatments, applied in the fixed
order below. :func:`peak_interval` applies them in that same order so a
combination resolves identically.

Roadmap: ``docs/MASTER_ROADMAP.md`` §1.1, issues #30 (codes) and #32
(regulation screen, step 1).
"""

from __future__ import annotations

import logging
import math
import re
from dataclasses import dataclass
from enum import Enum
from typing import FrozenSet, Iterable, List, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

#: peakfq's own stand-ins for zero and infinity (``vendor/peakfqr/R/main.R``).
Q_MIN: float = 1e-20
Q_MAX: float = 1e20

#: USGS peak-flow file qualification codes. Descriptions follow the USGS
#: peak-flow data documentation; the codes ``siteQT`` acts on are the ones
#: that matter numerically.
PEAK_CODE_DESCRIPTIONS: dict[str, str] = {
    "1": "Discharge is a maximum daily average",
    "2": "Discharge is an estimate",
    "3": "Discharge affected by dam failure",
    "4": "Discharge less than indicated value, which is minimum recordable discharge at this site",
    "5": "Discharge affected to unknown degree by regulation or diversion",
    "6": "Discharge affected by regulation or diversion",
    "7": "Discharge is an historic peak",
    "8": "Discharge actually greater than indicated value",
    "9": "Discharge due to snowmelt, hurricane, ice-jam or debris dam breakup",
    "A": "Year of occurrence is unknown or not exact",
    "B": "Month or day of occurrence is unknown or not exact",
    # The legacy NWIS peak RDB header lists these two-character codes rather
    # than a bare B; WATSTORE files use them too.
    "Bd": "Day of occurrence is unknown or not exact",
    "Bm": "Month of occurrence is unknown or not exact",
    "C": "All or part of the record affected by urbanization, mining, agricultural changes, "
    "channelization, or other",
    "F": "Peak supplied by another agency",
    "O": "Opportunistic value not from systematic data collection",
    "R": "Revised",
}

LESS_THAN_CODES: FrozenSet[str] = frozenset({"4"})
GREATER_THAN_CODES: FrozenSet[str] = frozenset({"8"})
ALWAYS_REMOVED_CODES: FrozenSet[str] = frozenset({"3", "O"})
URBAN_REGULATED_CODES: FrozenSet[str] = frozenset({"6", "C"})
HISTORIC_CODES: FrozenSet[str] = frozenset({"7"})
#: Regulation evidence for the gage screen. Code 5 is not acted on by siteQT,
#: but it is still evidence of regulation for classification purposes.
REGULATION_EVIDENCE_CODES: FrozenSet[str] = frozenset({"5", "6"})
ALTERATION_EVIDENCE_CODES: FrozenSet[str] = frozenset({"C"})


class PeakTreatment(str, Enum):
    """How a single peak enters EMA."""

    SYSTEMATIC = "systematic"
    """Exact value: lower = upper = peak."""
    LESS_THAN = "less_than"
    """Code 4: interval [0, peak]."""
    GREATER_THAN = "greater_than"
    """Code 8: interval [peak, inf)."""
    REMOVED = "removed"
    """Codes 3/O, or 6/C unless urban/regulated peaks are kept: interval [0, inf)."""


#: A float-formatted numeric code, as pandas renders an all-numeric code column
#: it inferred as float ("7.0").
_FLOAT_CODE = re.compile(r"^\s*(\d+)\.0*\s*$")


def _float_code(value: float) -> str:
    """``7.0`` -> ``"7"``; a non-integral float is kept as its repr and logged."""
    if float(value).is_integer():
        return str(int(value))
    return repr(value)


#: One code within a field: the two-character date-precision codes ``Bd``/``Bm``
#: (any case), else a single non-space character.
_CODE_TOKEN = re.compile(r"[Bb][DdMm]|\S")


def _split_token(token: str) -> List[str]:
    """One comma-separated piece of a ``peak_cd`` field, as individual codes.

    ``"7.0"`` (a float-inferred code column) is code ``"7"``. ``Bd`` and
    ``Bm`` are whole codes, normalised to exactly that spelling; single
    characters are upper-cased. Splitting ``Bd`` into characters first would
    invent a code D (and ``Bm`` a code M), which is what this used to do.
    """
    match = _FLOAT_CODE.match(token)
    if match:
        token = match.group(1)
    out = []
    for code in _CODE_TOKEN.findall(token):
        out.append("B" + code[1].lower() if len(code) == 2 else code.upper())
    return out


def parse_codes(peak_cd: object) -> FrozenSet[str]:
    """Split an NWIS ``peak_cd`` field into individual codes.

    Accepts comma-separated strings (``"2,6"``), run-together strings
    (``"26"``, WATSTORE's ``"4Bm"``), ``None``/NaN, and iterables of codes. The
    date-precision codes ``Bd`` (day unknown) and ``Bm`` (month unknown) are
    read whole, as the NWIS peak RDB and WATSTORE both write them. A
    float-formatted numeric code (``"7.0"`` or ``7.0``, from a code column
    pandas inferred as float) is read as the integer code. Unknown codes are
    kept and logged once per call so they are never silently dropped.

    Parameters
    ----------
    peak_cd : object
        The raw field value.

    Returns
    -------
    frozenset of str
        Upper-cased single-character codes, plus ``"Bd"``/``"Bm"`` as spelled.
    """
    if peak_cd is None:
        return frozenset()
    if isinstance(peak_cd, float):
        if math.isnan(peak_cd):  # NaN from pandas
            return frozenset()
        # A code column pandas inferred as float: 7.0 is code "7".
        peak_cd = _float_code(peak_cd)
    if isinstance(peak_cd, str):
        chars = [c for t in peak_cd.split(",") for c in _split_token(t)]
    else:
        chars = [c for item in peak_cd for c in _split_token(str(item))]  # type: ignore[attr-defined]
    codes = frozenset(chars)
    unknown = codes - PEAK_CODE_DESCRIPTIONS.keys()
    if unknown:
        logger.warning(
            "Unrecognised NWIS peak code(s) %s; carried as information only", sorted(unknown)
        )
    return codes


@dataclass(frozen=True)
class PeakInterval:
    """A peak's EMA flow interval and the reason for it."""

    lower: float
    upper: float
    treatment: PeakTreatment
    codes: FrozenSet[str]

    @property
    def is_historic(self) -> bool:
        """Whether the peak carries the historic-peak code (7)."""
        return bool(self.codes & HISTORIC_CODES)


def peak_interval(
    value: float,
    peak_cd: object = None,
    *,
    include_urban_regulated: bool = False,
) -> PeakInterval:
    """Build a peak's flow interval the way ``siteQT`` does.

    Applied in ``siteQT``'s order, so later rules win: code 4 sets the lower
    bound to zero, code 8 sets the upper bound to infinity, codes 3/O remove
    the peak, and codes 6/C remove it unless ``include_urban_regulated``
    (peakfq's ``Urb/Reg = Yes``).

    Parameters
    ----------
    value : float
        Peak discharge, cfs. Must be non-negative.
    peak_cd : object, optional
        Raw NWIS qualification code field.
    include_urban_regulated : bool, default False
        Keep code 6/C peaks as systematic observations.

    Returns
    -------
    PeakInterval

    Raises
    ------
    ValueError
        If ``value`` is negative (``siteQT`` stops on this too).
    """
    if value < 0:
        raise ValueError(f"Negative discharge {value}")
    codes = parse_codes(peak_cd)
    lower, upper = float(value), float(value)
    treatment = PeakTreatment.SYSTEMATIC
    if codes & LESS_THAN_CODES:
        lower, treatment = Q_MIN, PeakTreatment.LESS_THAN
    if codes & GREATER_THAN_CODES:
        upper, treatment = Q_MAX, PeakTreatment.GREATER_THAN
    if codes & ALWAYS_REMOVED_CODES:
        lower, upper, treatment = Q_MIN, Q_MAX, PeakTreatment.REMOVED
    if not include_urban_regulated and codes & URBAN_REGULATED_CODES:
        lower, upper, treatment = Q_MIN, Q_MAX, PeakTreatment.REMOVED
    return PeakInterval(lower=lower, upper=upper, treatment=treatment, codes=codes)


class RegulationClass(str, Enum):
    """Gage classification from peak-code evidence alone.

    There is deliberately no ``REFERENCE`` member: the absence of a code is
    not evidence that a basin is unregulated. Reference status needs NID/
    GAGES-II data (issue #32, steps 2-3).
    """

    REGULATED = "regulated"
    ALTERED = "altered"
    NO_CODE_EVIDENCE = "no_code_evidence"


def classify_from_codes(peak_cds: Iterable[object]) -> Tuple[RegulationClass, float]:
    """Classify a gage from its peaks' qualification codes.

    Parameters
    ----------
    peak_cds : iterable
        One raw ``peak_cd`` field per peak.

    Returns
    -------
    (RegulationClass, float)
        The class and the fraction of peaks carrying the evidence for it.
        Regulation outranks alteration when both are present.
    """
    parsed = [parse_codes(c) for c in peak_cds]
    n = len(parsed)
    if n == 0:
        return RegulationClass.NO_CODE_EVIDENCE, 0.0
    n_reg = sum(1 for c in parsed if c & REGULATION_EVIDENCE_CODES)
    if n_reg:
        return RegulationClass.REGULATED, n_reg / n
    n_alt = sum(1 for c in parsed if c & ALTERATION_EVIDENCE_CODES)
    if n_alt:
        return RegulationClass.ALTERED, n_alt / n
    return RegulationClass.NO_CODE_EVIDENCE, 0.0


#: Every code ``siteQT`` changes a peak's EMA treatment for.
ACTED_ON_CODES: FrozenSet[str] = (
    LESS_THAN_CODES
    | GREATER_THAN_CODES
    | ALWAYS_REMOVED_CODES
    | URBAN_REGULATED_CODES
    | HISTORIC_CODES
)


def count_acted_on_codes(peak_cds: Iterable[object]) -> dict[str, int]:
    """How many peaks carry each code that changes their EMA treatment.

    For an analysis path that does not apply codes, to say what it ignored.

    Parameters
    ----------
    peak_cds : iterable
        One raw ``peak_cd`` field per peak.

    Returns
    -------
    dict of str to int
        Code to peak count, only for codes in :data:`ACTED_ON_CODES` that occur.
    """
    counts: dict[str, int] = {}
    for raw in peak_cds:
        for code in parse_codes(raw) & ACTED_ON_CODES:
            counts[code] = counts.get(code, 0) + 1
    return dict(sorted(counts.items()))


#: Columns :func:`peak_frame_intervals` adds to a peak frame.
INTERVAL_COLUMNS: Tuple[str, ...] = ("lower", "upper", "treatment", "is_historic", "codes")


def peak_frame_intervals(
    peaks: pd.DataFrame,
    *,
    include_urban_regulated: bool = False,
) -> pd.DataFrame:
    """Apply :func:`peak_interval` to every row of a peak frame.

    The bridge from a peak frame -- the ``water_year`` / ``peak_flow_cfs`` /
    ``qualification_code`` shape :meth:`flowfreq.usgs.USGSgage.download_peak_flow`,
    :mod:`flowfreq.peak_sources` and :func:`flowfreq.watstore.read_watstore`
    all return -- to the per-peak EMA treatment ``siteQT`` gives it.

    Rows with no discharge (gage-height-only peaks) are dropped, as ``siteQT``
    drops them, and the count is logged.

    Parameters
    ----------
    peaks : pandas.DataFrame
        Must have ``water_year`` and ``peak_flow_cfs``. ``qualification_code``
        is optional; without it every peak is systematic.
    include_urban_regulated : bool, default False
        peakfq's ``Urb/Reg = Yes``: keep code 6/C peaks.

    Returns
    -------
    pandas.DataFrame
        A copy of ``peaks`` (discharge-less rows removed) with the columns in
        :data:`INTERVAL_COLUMNS` appended: ``lower``/``upper`` in cfs (with
        :data:`Q_MIN`/:data:`Q_MAX` for zero/infinity), ``treatment`` as a
        :class:`PeakTreatment`, ``is_historic`` (code 7) and ``codes`` (a
        frozenset).

    Raises
    ------
    ValueError
        On a missing required column or a negative discharge.
    """
    missing = [c for c in ("water_year", "peak_flow_cfs") if c not in peaks.columns]
    if missing:
        raise ValueError(f"Peak frame missing columns: {missing}")
    flows = pd.to_numeric(peaks["peak_flow_cfs"], errors="coerce")
    keep = flows.notna()
    if not keep.all():
        logger.warning(
            "Dropping %d peak(s) with no discharge (gage-height only), as siteQT does",
            int((~keep).sum()),
        )
    out = peaks.loc[keep].copy()
    out["peak_flow_cfs"] = flows[keep].astype(float)
    raw_codes = (
        out["qualification_code"].tolist()
        if "qualification_code" in out.columns
        else [None] * len(out)
    )
    rows = [
        peak_interval(float(v), c, include_urban_regulated=include_urban_regulated)
        for v, c in zip(out["peak_flow_cfs"], raw_codes)
    ]
    out["lower"] = [r.lower for r in rows]
    out["upper"] = [r.upper for r in rows]
    out["treatment"] = [r.treatment for r in rows]
    out["is_historic"] = [r.is_historic for r in rows]
    out["codes"] = pd.Series([r.codes for r in rows], index=out.index, dtype=object)
    return out
