"""Matching annual exceedance probabilities that different sources write differently.

A Bulletin 17C fit tabulates the 1.5-year flood at AEP 0.67 and the 2.33-year flood at
0.43. A regression report prints 0.667 and 0.429 for the same floods, and a table of change
factors may print either. Compared exactly, the two never meet: the row is dropped, or
labelled extrapolated, or left unweighted, without a word. Everything that pairs AEPs from
two sources goes through this module instead.
"""

from __future__ import annotations

from typing import Iterable, Optional

#: How far apart, as a fraction of the larger, two AEPs may be and still be the same flood.
#: Relative, so that 0.005 and 0.002 are never taken for each other: 0.667 and 0.67 differ
#: by 0.45 percent, 0.429 and 0.43 by 0.23 percent, and the nearest pair of standard AEPs
#: (0.01 and 0.005) by 50 percent.
AEP_MATCH: float = 0.0075


def same_aep(a: float, b: float, tolerance: float = AEP_MATCH) -> bool:
    """True when two AEPs are the same flood written to different precision.

    >>> same_aep(0.67, 0.667), same_aep(0.005, 0.002), same_aep(0.01, 0.01)
    (True, False, True)
    """
    return abs(float(a) - float(b)) <= tolerance * max(abs(float(a)), abs(float(b)))


def match_aep(
    candidates: Iterable[float], aep: float, tolerance: float = AEP_MATCH
) -> Optional[float]:
    """The candidate that is the same flood as ``aep``, nearest first, or None.

    >>> match_aep([0.8, 0.667, 0.5], 0.67)
    0.667
    >>> match_aep([0.005, 0.01], 0.002) is None
    True
    """
    best: Optional[float] = None
    for candidate in candidates:
        value = float(candidate)
        if same_aep(value, aep, tolerance) and (best is None or abs(value - aep) < abs(best - aep)):
            best = value
    return best


__all__ = ["AEP_MATCH", "match_aep", "same_aep"]
