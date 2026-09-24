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
from dataclasses import dataclass
from enum import Enum
from typing import FrozenSet, Iterable, Tuple

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


def parse_codes(peak_cd: object) -> FrozenSet[str]:
    """Split an NWIS ``peak_cd`` field into individual codes.

    Accepts comma-separated strings (``"2,6"``), run-together strings
    (``"26"``), ``None``/NaN, and iterables of codes. Unknown codes are kept
    and logged once per call so they are never silently dropped.

    Parameters
    ----------
    peak_cd : object
        The raw field value.

    Returns
    -------
    frozenset of str
        Upper-cased single-character codes.
    """
    if peak_cd is None:
        return frozenset()
    if isinstance(peak_cd, float):  # NaN from pandas
        return frozenset()
    if isinstance(peak_cd, str):
        chars = [c for c in peak_cd.upper() if not c.isspace() and c != ","]
    else:
        chars = [str(c).strip().upper() for c in peak_cd if str(c).strip()]  # type: ignore[attr-defined]
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
