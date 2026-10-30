"""Flood-frequency estimates at ungaged sites from published regression equations.

This module holds :class:`SiteEstimate`, the result of an estimate at an ungaged
location; :func:`assemble`, the pure assembly of that estimate from inputs a caller
has already fetched (the located NSS regression regions, the basin characteristics,
and a state's equation library); and :func:`estimate_at_site`, the entry point that
fetches those inputs -- the library from the packaged data, the delineation and the
region location from StreamStats -- and assembles them, weighting in a caller's
at-site fit when the basin lies in a single region.

Assembly rules
--------------
- Each located region is mapped to a library region with
  :func:`region_for_nss_code`. A region with no library match is reported in
  ``unmatched`` and the estimate proceeds on the rest; if none matches,
  :class:`~flowfreq.regression.library.EquationsUnavailable` names the codes.
- A matched region with ``percent_weight`` 0 (NSS rounds to a whole percent, so a sliver
  under 0.5% arrives as 0) contributes nothing and is dropped before anything is
  evaluated, so a basin effectively in one region keeps that region's interval.
- A matched region whose equations need a characteristic the basin does not carry is
  skipped and the missing codes are reported by its NSS code. The equations checked are
  those at the requested AEPs when the region publishes every one of them, else (no
  request, or a requested AEP it does not publish) every AEP it publishes.
- The default AEPs are those every region that remains publishes, largest first; none in
  common raises :class:`~flowfreq.regression.library.EquationsUnavailable`.
- NSS serves superseded and current equation sets for one area side by side, each at
  its full ``percent_weight``. When the remaining regions' weights sum to more than
  100 % (beyond NSS's whole-percent rounding) they are alternatives, not a partition:
  ``region=`` chooses one and the others are reported as ``alternatives``; without a
  choice :class:`OverlappingRegionsError` lists the candidates. ``region=`` may instead be
  a sequence of codes, one per side of a region boundary, whose ``percent_weight`` sums
  to 100 % (within the same tolerance): they are used together as a partition and the
  rest are ``alternatives``.
- Area fractions come from NSS's ``percent_weight`` (what NSS itself weights by),
  renormalised over the regions actually used.
- One region gives its own prediction interval; several give an area-weighted flow with
  no interval and ``interval_method="none: area-weighted across N regions"``.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple, Union

import numpy as np
import pandas as pd

from flowfreq.core import FrequencyResults
from flowfreq.gage_weighting import require_variance_log, weight_with_regression
from flowfreq.regression.equations import RegressionEstimate, evaluate, evaluate_weighted
from flowfreq.regression.jurisdictions import jurisdiction
from flowfreq.regression.library import (
    OPT_IN_SUBDIRS,
    EquationsUnavailable,
    StateLibrary,
    available_states,
    load_state,
    nss_snapshot_states,
)
from flowfreq.streamstats import (
    PERCENT_WEIGHT_TOL,
    RegressionRegionLocation,
    StreamStatsCache,
    delineate_and_get_characteristics,
    locate_regression_regions,
)

logger = logging.getLogger(__name__)

__all__ = [
    "Assembly",
    "OverlappingRegionsError",
    "SiteEstimate",
    "assemble",
    "estimate_at_site",
    "region_for_nss_code",
]

#: Columns of :attr:`SiteEstimate.quantiles`, in order.
QUANTILE_COLUMNS: Tuple[str, ...] = (
    "aep",
    "flow_cfs",
    "lower",
    "upper",
    "interval_level",
    "interval_method",
    "extrapolated",
    "regions",
)

#: ``(fraction, estimate)`` pairs for one AEP, one per region used.
RegionEstimates = List[Tuple[float, RegressionEstimate]]


def region_for_nss_code(lib: StateLibrary, code: str) -> Optional[str]:
    """Map an NSS region code to a library ``region_code``.

    Equality with an equation's ``region_code`` wins; otherwise the first
    ``region_code`` (in the library's own order) whose ``region_name`` contains the NSS
    code as a whole token. Verified files write names such as
    ``"Flood region 1 (SIR 2014-5211; NSS GC1618)"``; the token rule keeps ``GC161``
    from matching that name, and ``GC1`` from matching ``GC1618``.

    Parameters
    ----------
    lib : StateLibrary
    code : str
        The NSS region code, as :class:`~flowfreq.streamstats.RegressionRegionLocation`
        carries it.

    Returns
    -------
    str or None
        The library region code, or ``None`` when nothing matches.
    """
    if any(e.region_code == code for e in lib.equations):
        return code
    token = re.compile(rf"(?<![A-Za-z0-9]){re.escape(code)}(?![A-Za-z0-9])")
    for e in lib.equations:
        if e.region_name and token.search(e.region_name):
            return e.region_code
    return None


def _needed_codes(lib: StateLibrary, region_code: str, aeps: Sequence[float]) -> List[str]:
    """Every characteristic code the region's equations at ``aeps`` require."""
    codes: List[str] = []
    for aep in aeps:
        eq = lib.equation(region_code, aep)
        for v in (*eq.variables, *eq.applicability):
            if v.code not in codes:
                codes.append(v.code)
    return codes


def _codes(regions: Sequence[RegressionRegionLocation]) -> str:
    return ", ".join(r.code for r in regions) or "none"


#: ``(NSS code, library region name, citation year)`` for one candidate region.
Candidate = Tuple[str, str, Optional[int]]


class OverlappingRegionsError(EquationsUnavailable):
    """The located NSS regions overlap: they are alternative equation sets, not a partition.

    NSS serves superseded and current equation sets for the same area side by side,
    each at its full ``percent_weight``; averaging them would blend two estimates of the
    same thing. A caller chooses one with ``region=`` (API) or ``--region`` (CLI).

    Parameters
    ----------
    state : str
        The jurisdiction.
    candidates : list of (str, str, int or None)
        ``(NSS code, library region name, citation year)`` for every region that could
        be used, in the located order.
    region : str, optional
        A choice that matched none of the candidates, named in the message.
    reason : str, optional
        Replaces the message's opening when the problem is more specific than "no
        choice" or "an unknown choice": a choice that NSS located but that was skipped,
        or a choice that covers only part of the basin (the regions overlap on a region
        boundary). The candidates are still listed after it.

    Attributes
    ----------
    candidates : list of (str, str, int or None)
    """

    def __init__(
        self,
        state: str,
        candidates: Sequence[Candidate],
        region: Optional[str] = None,
        *,
        reason: Optional[str] = None,
    ) -> None:
        self.candidates: List[Candidate] = list(candidates)
        listing = "\n".join(
            f"  {code}: {name} ({year if year is not None else 'year unknown'})"
            for code, name, year in self.candidates
        )
        if reason is not None:
            super().__init__(f"{state}: {reason}:\n{listing}")
            return
        head = (
            f"{state}: region {region!r} is not one of the overlapping NSS regions"
            if region is not None
            else f"{state}: the located NSS regions overlap (alternative equation sets for "
            "the same area, not a partition)"
        )
        super().__init__(
            f"{head}; choose one with region= (API) or --region (CLI), by NSS code or "
            f"library region code:\n{listing}"
        )


def _region_name(lib: StateLibrary, lib_code: str) -> str:
    """The library's name for a region (its first named equation's), else the code."""
    for e in lib.equations:
        if e.region_code == lib_code and e.region_name:
            return e.region_name
    return lib_code


def _region_citation(lib: StateLibrary, lib_code: str) -> Tuple[str, Optional[int]]:
    """``(citation text, year)`` of a region's first equation, as ``to_markdown`` cites it."""
    for e in lib.equations:
        if e.region_code == lib_code:
            c = e.citation
            text = c.publication + (f" ({c.year})" if c.year else "") + f", table {c.table}"
            return text, c.year
    return "", None


def _region_choices(region: Union[str, Sequence[str], None]) -> List[str]:
    """The ``region`` argument as a list of codes (empty when no choice was made).

    Raises
    ------
    ValueError
        If a sequence is empty or repeats a code (case-insensitively).
    """
    if region is None:
        return []
    if isinstance(region, str):
        return [region]
    codes = list(region)
    if not codes:
        raise ValueError("region must not be an empty sequence (use None for no choice)")
    seen: Set[str] = set()
    for code in codes:
        if code.upper() in seen:
            raise ValueError(f"region repeats the code {code!r}")
        seen.add(code.upper())
    return codes


def _why_not_in_use(
    choice: str,
    regions: Sequence[RegressionRegionLocation],
    region_codes: Mapping[str, str],
    unmatched: Sequence[RegressionRegionLocation],
    zero_weight: Sequence[RegressionRegionLocation],
    missing: Mapping[str, List[str]],
) -> Optional[str]:
    """Why a located region that ``choice`` names (NSS or library code) is not in use.

    ``None`` when ``choice`` names no located region at all.
    """
    wanted = choice.upper()
    for loc in regions:
        if loc.code.upper() != wanted and region_codes.get(loc.code, "").upper() != wanted:
            continue
        if any(u is loc for u in unmatched):
            return f"NSS located {loc.code}, but there is no library region for it"
        if any(z is loc for z in zero_weight):
            return f"NSS located {loc.code} at 0 % of the basin, so it was dropped"
        if loc.code in missing:
            return (
                f"NSS located {loc.code}, but it was skipped because the basin lacks "
                + ", ".join(missing[loc.code])
            )
    return None


@dataclass
class Assembly:
    """What :func:`assemble` produces: the estimate and what was left out of it.

    Attributes
    ----------
    quantiles : pandas.DataFrame
        One row per AEP with :data:`QUANTILE_COLUMNS`. ``regions`` is the NSS codes
        used, comma-separated; ``lower``, ``upper`` and ``interval_level`` are NaN where
        there is no interval.
    estimates : dict
        AEP to the ``(fraction, RegressionEstimate)`` pairs used for that row, in the
        order of ``used``. A single region carries fraction 1.0.
    used : list of RegressionRegionLocation
        The located regions the estimate rests on, in their original order.
    region_codes : dict
        NSS code to library ``region_code`` for every matched region (used or skipped).
    unmatched : list of RegressionRegionLocation
        Located regions with no library match.
    zero_weight : list of RegressionRegionLocation
        Matched regions with ``percent_weight`` 0, dropped unevaluated.
    missing : dict
        NSS region code to the characteristic codes the basin lacks, for every matched
        region skipped for that reason.
    alternatives : list of RegressionRegionLocation
        Overlapping regions set aside when ``region`` chose another (empty when the
        regions partition the basin).
    region_choice_source : str, optional
        Set when ``region`` named a library region that NSS did not locate (NSS located no
        peak-flow region mapping to the library); ``None`` on every other path.
    """

    quantiles: pd.DataFrame
    estimates: Dict[float, RegionEstimates]
    used: List[RegressionRegionLocation]
    region_codes: Dict[str, str]
    unmatched: List[RegressionRegionLocation]
    zero_weight: List[RegressionRegionLocation]
    missing: Dict[str, List[str]]
    alternatives: List[RegressionRegionLocation] = field(default_factory=list)
    region_choice_source: Optional[str] = None


def _library_region(lib: StateLibrary, code: str) -> Optional[str]:
    """Resolve a caller's code to a library ``region_code`` (case-insensitively), else None.

    Equal to a ``region_code`` in ``lib.regions``, or a whole-token match inside some
    equation's ``region_name`` (the token rule of :func:`region_for_nss_code`).
    """
    wanted = code.upper()
    for lib_code in lib.regions:
        if lib_code.upper() == wanted:
            return lib_code
    token = re.compile(rf"(?<![A-Za-z0-9]){re.escape(code)}(?![A-Za-z0-9])", re.IGNORECASE)
    for e in lib.equations:
        if e.region_name and token.search(e.region_name):
            return e.region_code
    return None


def assemble(
    lib: StateLibrary,
    regions: Sequence[RegressionRegionLocation],
    characteristics: Mapping[str, float],
    *,
    aeps: Optional[Sequence[float]],
    allow_extrapolation: bool,
    weighting_space: str,
    region: Union[str, Sequence[str], None] = None,
) -> Assembly:
    """Assemble the regression estimate for a basin from pre-fetched inputs.

    Parameters
    ----------
    lib : StateLibrary
        The jurisdiction's equations.
    regions : sequence of RegressionRegionLocation
        The NSS regression regions the basin overlaps, with their ``percent_weight``.
    characteristics : mapping of str to float
        Basin characteristics keyed by StreamStats code.
    aeps : sequence of float, optional
        AEPs to estimate, kept in the caller's order. ``None`` means every AEP all the
        regions used publish, largest first. A region is skipped for a characteristic
        the basin lacks only if an equation it would be evaluated with needs it: the
        equations at these AEPs when the region publishes every one of them, else
        (``None``, or a requested AEP the region does not publish) every equation it
        publishes, so a region that lacks a characteristic is skipped before an
        unpublished requested AEP could be raised for it.
    allow_extrapolation : bool
        Evaluate outside an equation's calibrated range and flag the row, instead of
        raising :class:`~flowfreq.regression.equations.OutOfRangeError`.
    weighting_space : {"log", "linear"}
        How a multi-region basin's estimates are area-averaged
        (:func:`~flowfreq.regression.equations.evaluate_weighted`).
    region : str, optional
        The region to use when the regions overlap, by NSS code or library
        ``region_code``. Overlap means the ``percent_weight`` of the regions remaining
        after the zero-weight and missing-characteristic steps sums to more than 100 %
        plus NSS's rounding tolerance (the 1.5 points ``streamstats`` allows, or half a
        point per region when that is more, since each weight is rounded to a whole
        percent); the chosen region is then used alone (fraction 1.0) and the others are
        returned in ``alternatives``; it must cover the whole basin (``percent_weight`` within the
        tolerance of 100). A partition is unaffected; ``region`` naming one of its
        regions is accepted as a no-op. A ``region`` naming no region in use is always
        refused. May also be a sequence of str: several codes that together cover the
        basin (one per side of a region boundary whose sides each have overlapping
        sets). Their ``percent_weight`` must sum to within the tolerance of 100; they
        are area-weighted (fractions renormalised) in their located order and every
        other region goes to ``alternatives``. An empty sequence, a repeated code
        (case-insensitively) or one region named twice under its NSS and library codes
        is a ``ValueError``; a one-element sequence is the same as that string. When NSS
        located no region (``regions`` empty) or none maps to the library, a single code
        that resolves against the library (a ``region_code``, or a whole token in an
        equation's ``region_name``, case-insensitively) is used alone as a region of
        ``percent_weight`` 100, recorded in ``region_codes``, with
        ``region_choice_source`` set; otherwise the refusal also says the code is not a
        library region (a sequence of two or more codes never resolves this way).

    Returns
    -------
    Assembly

    Raises
    ------
    OverlappingRegionsError
        If the regions overlap and ``region`` is not given, matches none of them, or
        names one covering only part of the basin (overlap on a region boundary), or
        names a library region that more than one located NSS region reaches.
    EquationsUnavailable
        If ``region`` names no region in use on a basin without overlap (the message
        says why when NSS located it: no library region, zero weight, or a missing
        characteristic).

        If NSS located no region, or no located region matches the library; if every
        matched region has zero weight or lacks a basin characteristic; if the regions
        used publish no AEP in common and none was requested; or if a requested AEP is
        not published for a region used.
    OutOfRangeError
        If an input is outside its range and ``allow_extrapolation`` is False.
    ValueError
        An empty ``aeps``; or, from :func:`evaluate_weighted`, an unknown
        ``weighting_space`` or a printed ``Q=0`` weighted in log space.
    """
    if aeps is not None and len(aeps) == 0:
        raise ValueError("aeps must be non-empty (use None for the regions' common AEPs)")
    choices = _region_choices(region)
    choice = choices[0] if choices else None
    matched: List[Tuple[RegressionRegionLocation, str]] = []
    unmatched: List[RegressionRegionLocation] = []
    for loc in regions:
        lib_code = region_for_nss_code(lib, loc.code)
        if lib_code is None:
            unmatched.append(loc)
        else:
            matched.append((loc, lib_code))
    region_choice_source: Optional[str] = None
    extra_refusal = ""
    pre_matched: Dict[str, str] = {}
    if region is not None and (not regions or not matched):
        resolved = _library_region(lib, choice) if len(choices) == 1 and choice else None
        if resolved is not None and choice is not None:
            located = RegressionRegionLocation(
                region_id=0,
                code=choice.upper(),
                name=_region_name(lib, resolved),
                percent_weight=100.0,
                area_sq_mi=0.0,
            )
            matched = [(located, resolved)]
            pre_matched = {located.code: resolved}
            region_choice_source = "caller; NSS located no peak-flow region for this basin"
        else:
            shown = choice if len(choices) == 1 else region
            extra_refusal = f"; region {shown!r} is not a library region either"
    if not regions and not matched:
        raise EquationsUnavailable(
            f"{lib.state}: NSS located no regression region for this basin in "
            f"{lib.state}; check the state and that the point is on a stream NSS covers"
            + extra_refusal
        )
    if not matched:
        raise EquationsUnavailable(
            f"{lib.state}: none of the located NSS regions ({_codes(regions)}) matches a "
            f"library region ({', '.join(lib.regions) or 'none'})" + extra_refusal
        )
    if unmatched:
        logger.warning(
            "%s: no library region for located NSS region(s) %s; estimating on %s",
            lib.state,
            _codes(unmatched),
            ", ".join(r.code for r, _ in matched),
        )
    region_codes = {loc.code: lib_code for loc, lib_code in matched}
    region_codes.update(pre_matched)

    zero_weight = [loc for loc, _ in matched if not loc.percent_weight > 0]
    weighted = [(loc, lib_code) for loc, lib_code in matched if loc.percent_weight > 0]
    if not weighted:
        raise EquationsUnavailable(
            f"{lib.state}: every matched NSS region ({_codes(zero_weight)}) has zero weight"
        )
    if zero_weight:
        logger.info("%s: dropping zero-weight NSS region(s) %s", lib.state, _codes(zero_weight))

    missing: Dict[str, List[str]] = {}
    used: List[Tuple[RegressionRegionLocation, str]] = []
    for loc, lib_code in weighted:
        # The equations at the requested AEPs when the region publishes every one of
        # them, so an equation the request never reaches cannot skip the region. Else
        # (no request, or a requested AEP it does not publish) every equation it
        # publishes: a region that lacks a characteristic is then skipped before an
        # unpublished requested AEP could be raised for it.
        published = list(lib.aeps(lib_code))
        check = published
        # The same 1e-12 match lib.equation uses.
        if aeps is not None and all(
            any(abs(p - float(a)) < 1e-12 for p in published) for a in aeps
        ):
            check = [float(a) for a in aeps]
        absent = [c for c in _needed_codes(lib, lib_code, check) if c not in characteristics]
        if absent:
            missing[loc.code] = absent
            logger.warning(
                "%s: skipping region %s, basin characteristics missing %s",
                lib.state,
                loc.code,
                absent,
            )
        else:
            used.append((loc, lib_code))
    if not used:
        raise EquationsUnavailable(
            f"{lib.state}: every matched region needs a basin characteristic the site "
            "lacks: "
            + "; ".join(f"{code} needs {', '.join(codes)}" for code, codes in missing.items())
        )

    # A choice must name a region in use (an NSS code is unambiguous; a library code
    # is the fallback), overlap or not: otherwise the estimate would rest on regions
    # the caller did not choose, with the choice recorded as if it had been applied.
    chosen: Optional[Tuple[RegressionRegionLocation, str]] = None
    ambiguous = False
    why: Optional[str] = None
    if choice is not None:
        wanted = choice.upper()
        chosen = next((pair for pair in used if pair[0].code.upper() == wanted), None)
        if chosen is None:
            by_library = [pair for pair in used if pair[1].upper() == wanted]
            ambiguous = len(by_library) > 1
            chosen = by_library[0] if by_library else None
        if chosen is None:
            why = _why_not_in_use(choice, regions, region_codes, unmatched, zero_weight, missing)

    # Further codes of a sequence resolve the same way; the first one that does not (or
    # that is ambiguous) becomes the choice reported.
    extra: List[Tuple[RegressionRegionLocation, str]] = []
    if chosen is not None and not ambiguous:
        for code in choices[1:]:
            wanted = code.upper()
            found = next((pair for pair in used if pair[0].code.upper() == wanted), None)
            by_library = []
            if found is None:
                by_library = [pair for pair in used if pair[1].upper() == wanted]
                found = by_library[0] if by_library else None
            if found is None:
                choice, chosen, extra = code, None, []
                why = _why_not_in_use(code, regions, region_codes, unmatched, zero_weight, missing)
                break
            if len(by_library) > 1:
                choice, chosen, extra, ambiguous = code, found, [], True
                break
            if any(found[0] is r for r, _ in [chosen, *extra]):
                # An NSS code and the same region's library code pass the string check
                # in _region_choices; counted twice, a 50 % side would "cover" the basin.
                raise ValueError(
                    f"region {code!r} names the same region as an earlier code "
                    f"({found[0].code}): list each side of the basin once"
                )
            extra.append(found)

    # NSS rounds each percent_weight to a whole percent, so an N-region partition sums
    # to 100 +/- 0.5 N; a real overlap sums to about 200, so scaling cannot mask one.
    overlap_tol = max(PERCENT_WEIGHT_TOL, 0.5 * len(used))
    overlap = sum(r.percent_weight for r, _ in used) > 100.0 + overlap_tol
    if choice is not None and chosen is None and not overlap:
        raise EquationsUnavailable(
            f"{lib.state}: region {choice!r} is not in use: "
            + (why if why is not None else f"NSS located only {_codes(regions)}")
            + f"; the estimate would rest on {', '.join(r.code for r, _ in used)}"
        )

    alternatives: List[RegressionRegionLocation] = []
    if overlap:
        candidates = [
            (r.code, _region_name(lib, lib_code), _region_citation(lib, lib_code)[1])
            for r, lib_code in used
        ]
        if choice is None:
            raise OverlappingRegionsError(lib.state, candidates)
        if chosen is None:
            if why is None:
                raise OverlappingRegionsError(lib.state, candidates, choice)
            raise OverlappingRegionsError(
                lib.state,
                candidates,
                choice,
                reason=f"region {choice!r} is not in use ({why}); the located NSS regions "
                "overlap, so choose one of these with region= (API) or --region (CLI)",
            )
        if ambiguous:
            raise OverlappingRegionsError(
                lib.state,
                candidates,
                choice,
                reason=f"region {choice!r} is ambiguous: library region {chosen[1]} is "
                f"reached by more than one located NSS region; choose by NSS code with "
                "region= (API) or --region (CLI)",
            )
        picked = [chosen] + extra
        if extra:
            picked_sum = sum(r.percent_weight for r, _ in picked)
            picked_tol = max(PERCENT_WEIGHT_TOL, 0.5 * len(picked))
            picked_codes = ", ".join(r.code for r, _ in picked)
            if picked_sum < 100.0 - picked_tol:
                raise OverlappingRegionsError(
                    lib.state,
                    candidates,
                    choice,
                    reason=f"region choice {picked_codes} covers only {picked_sum:g} % of the "
                    "basin: the located NSS regions overlap on a region boundary, so the "
                    "chosen equation sets do not cover the whole basin and choosing them "
                    "would drop the rest of it",
                )
            if picked_sum > 100.0 + picked_tol:
                raise OverlappingRegionsError(
                    lib.state,
                    candidates,
                    choice,
                    reason=f"the chosen regions ({picked_codes}) cover {picked_sum:g} % of "
                    "the basin, i.e. they overlap each other",
                )
        elif chosen[0].percent_weight < 100.0 - PERCENT_WEIGHT_TOL:
            # Overlapping sets on a region boundary: each side has its own
            # alternatives, and one set covers only its side of the basin.
            raise OverlappingRegionsError(
                lib.state,
                candidates,
                choice,
                reason=f"region {choice!r} ({chosen[0].code}) covers only "
                f"{chosen[0].percent_weight:g} % of the basin: the located NSS regions "
                "overlap on a region boundary, so no single equation set covers the "
                "whole basin and choosing one would drop the rest of it",
            )
        picked_ids = {id(r) for r, _ in picked}
        alternatives = [r for r, _ in used if id(r) not in picked_ids]
        logger.info(
            "%s: overlapping NSS regions; using %s, alternatives %s",
            lib.state,
            ", ".join(r.code for r, _ in picked),
            _codes(alternatives),
        )
        used = [pair for pair in used if id(pair[0]) in picked_ids]
    # A choice in use on a partition is a no-op: there is nothing to choose between.

    if aeps is None:
        common = set.intersection(*(set(lib.aeps(lib_code)) for _, lib_code in used))
        aep_list = sorted(common, reverse=True)
        if not aep_list:
            raise EquationsUnavailable(
                f"{lib.state}: the regions used publish no AEP in common: "
                + "; ".join(
                    f"{loc.code} at {', '.join(f'{a:g}' for a in lib.aeps(lib_code))}"
                    for loc, lib_code in used
                )
            )
    else:
        aep_list = [float(a) for a in aeps]

    total = sum(loc.percent_weight for loc, _ in used)
    fractions = [loc.percent_weight / total for loc, _ in used]
    region_label = ", ".join(loc.code for loc, _ in used)

    estimates: Dict[float, RegionEstimates] = {}
    rows: List[Dict[str, Any]] = []
    for aep in aep_list:
        per_region: RegionEstimates = [
            (
                fraction,
                evaluate(
                    lib.equation(lib_code, aep),
                    characteristics,
                    allow_extrapolation=allow_extrapolation,
                ),
            )
            for fraction, (_, lib_code) in zip(fractions, used)
        ]
        estimates[aep] = per_region
        if len(per_region) == 1:
            est = per_region[0][1]
            lower, upper = est.interval if est.interval is not None else (np.nan, np.nan)
            rows.append(
                {
                    "aep": aep,
                    "flow_cfs": est.flow_cfs,
                    "lower": lower,
                    "upper": upper,
                    "interval_level": est.interval_level if est.interval_level else np.nan,
                    "interval_method": est.interval_method,
                    "extrapolated": est.extrapolated,
                    "regions": region_label,
                }
            )
        else:
            rows.append(
                {
                    "aep": aep,
                    "flow_cfs": evaluate_weighted(per_region, space=weighting_space),
                    "lower": np.nan,
                    "upper": np.nan,
                    "interval_level": np.nan,
                    "interval_method": f"none: area-weighted across {len(per_region)} regions",
                    "extrapolated": any(e.extrapolated for _, e in per_region),
                    "regions": region_label,
                }
            )
    return Assembly(
        quantiles=pd.DataFrame(rows, columns=list(QUANTILE_COLUMNS)),
        estimates=estimates,
        used=[loc for loc, _ in used],
        region_codes=region_codes,
        unmatched=unmatched,
        zero_weight=zero_weight,
        missing=missing,
        alternatives=alternatives,
        region_choice_source=region_choice_source,
    )


@dataclass
class SiteEstimate:
    """A flood-frequency estimate at an ungaged site, with its provenance.

    Attributes
    ----------
    lat, lon : float
        The site, decimal degrees.
    state : str
        The jurisdiction whose equation library was used.
    method : str
        ``"regression"``, or ``"at-site weighted with regression"`` when a caller's
        at-site fit was weighted in.
    quantiles : pandas.DataFrame
        One row per AEP with :data:`QUANTILE_COLUMNS` (see :class:`Assembly`).
    regions : list of RegressionRegionLocation
        The located regions the estimate used (:attr:`Assembly.used`).
    unmatched_regions : list of RegressionRegionLocation
        Located regions with no library match, excluded and renormalised away.
    characteristics : dict
        Basin characteristics by StreamStats code.
    library_status : str
        The library's status (``verified``, ``provisional``, ...).
    advisory : str
        The library's advisory text, if any.
    estimates : dict
        AEP to ``(fraction, RegressionEstimate)`` pairs, one per region used.
    region_codes : dict
        NSS code to library ``region_code`` (:attr:`Assembly.region_codes`); it is how
        :meth:`to_markdown` finds each located region's citation.
    weighted : pandas.DataFrame, optional
        Bulletin 17C Appendix 9 weighting of an at-site fit with the regression, when
        the caller supplied one and the basin lies in a single region.
    provenance : dict
        How the inputs were obtained: snap, delineation, library file status,
        out-of-range codes, timestamps. :meth:`to_markdown` reads two keys the caller
        fills from the :class:`Assembly`: ``"missing_characteristics"``
        (:attr:`Assembly.missing`) for the skipped-region caveat and
        ``"zero_weight_regions"`` (the NSS codes of :attr:`Assembly.zero_weight`) for
        the dropped-region caveat, and ``"alternative_regions"`` (dicts with ``code``,
        ``name``, ``region_code``, ``citation``, ``year``) for the alternative
        equation sets not used.
    """

    lat: float
    lon: float
    state: str
    method: str
    quantiles: pd.DataFrame
    regions: List[RegressionRegionLocation]
    unmatched_regions: List[RegressionRegionLocation]
    characteristics: Dict[str, float]
    library_status: str
    advisory: str
    estimates: Dict[float, RegionEstimates]
    region_codes: Dict[str, str] = field(default_factory=dict)
    weighted: Optional[pd.DataFrame] = None
    provenance: Dict[str, Any] = field(default_factory=dict)

    def _citations(self) -> Dict[str, str]:
        """Library region code to the publication its equations cite."""
        cited: Dict[str, str] = {}
        for per_region in self.estimates.values():
            for _, e in per_region:
                c = e.equation.citation
                cited.setdefault(
                    e.equation.region_code,
                    c.publication + (f" ({c.year})" if c.year else "") + f", table {c.table}",
                )
        return cited

    def to_markdown(self) -> str:
        """A report-ready account of the estimate: regions, inputs, quantiles, caveats.

        Returns
        -------
        str
            Markdown.
        """
        lines = [
            f"# Regression estimate at {self.lat:.4f}, {self.lon:.4f} ({self.state})",
            "",
            f"- Method: {self.method}",
            f"- Equation library: {self.state}, status {self.library_status}",
        ]
        if self.advisory:
            lines.append(f"- Advisory: {self.advisory}")
        lines.extend(["", "## Regression regions", ""])
        if self.regions:
            lines.extend(
                [
                    "| NSS code | Name | Area fraction | Source |",
                    "|---|---|---:|---|",
                ]
            )
            citations = self._citations()
            # The fractions are what assemble() used: percent_weight renormalised over
            # the regions used. Each located region finds its library equation, and so
            # its citation, through region_codes rather than by position.
            total = sum(r.percent_weight for r in self.regions)
            for r in self.regions:
                fraction = r.percent_weight / total if total > 0 else float("nan")
                source = citations.get(self.region_codes.get(r.code, r.code), "")
                lines.append(f"| {r.code} | {r.name} | {fraction:.3f} | {source} |")
        else:
            lines.append("No regression region was used.")
        alternatives = self.provenance.get("alternative_regions")
        if alternatives:
            lines.extend(["", "Not used (alternative NSS equation sets):", ""])
            lines.extend(
                f"- {a['code']} ({a['name']}): {a['citation'] or 'citation unknown'}"
                for a in alternatives
            )
        if self.unmatched_regions:
            lines.append("")
            lines.append(
                "Located but not in the library, so excluded and the fractions renormalised: "
                + ", ".join(
                    f"{r.code} ({r.name}, {r.percent_weight:g}%)" for r in self.unmatched_regions
                )
                + "."
            )

        lines.extend(["", "## Basin characteristics", ""])
        if self.characteristics:
            lines.extend(["| Code | Value |", "|---|---:|"])
            lines.extend(f"| {code} | {value:g} |" for code, value in self.characteristics.items())
        else:
            lines.append("None.")

        lines.extend(["", "## Quantiles", ""])
        if len(self.quantiles):
            lines.extend(
                [
                    "| AEP (%) | Flow (cfs) | Lower | Upper | Level | Interval | Extrapolated "
                    "| Regions |",
                    "|---:|---:|---:|---:|---:|---|---|---|",
                ]
            )
            for _, row in self.quantiles.iterrows():
                lines.append(
                    f"| {row['aep'] * 100:g} | {row['flow_cfs']:,.0f} | {_cfs(row['lower'])} | "
                    f"{_cfs(row['upper'])} | {_level(row['interval_level'])} | "
                    f"{row['interval_method']} | {'yes' if row['extrapolated'] else 'no'} | "
                    f"{row['regions']} |"
                )
        else:
            lines.append("None.")

        if self.weighted is not None:
            lines.extend(["", "## At-site weighted with regression", ""])
            lines.extend(_frame_table(self.weighted))

        lines.extend(["", "## Caveats", ""])
        extrapolated = (
            self.quantiles["aep"][self.quantiles["extrapolated"].astype(bool)]
            if len(self.quantiles)
            else pd.Series(dtype=float)
        )
        if len(extrapolated):
            lines.append(
                f"- **{len(extrapolated)} of {len(self.quantiles)} quantiles are extrapolated**, "
                "at AEP " + ", ".join(f"{a:g}" for a in extrapolated) + ": a basin characteristic "
                "lies outside the calibrated range of the equation, so the published "
                "standard error does not apply there."
            )
        else:
            lines.append("- Every input lies within the calibrated range of the equations used.")
        if len(self.regions) > 1:
            lines.append(
                f"- The basin spans {len(self.regions)} regression regions; the flows are "
                "area-weighted and carry no prediction interval."
            )
        missing = self.provenance.get("missing_characteristics")
        if missing:
            lines.append(
                "- Regions skipped for missing basin characteristics: "
                + "; ".join(f"{code} needs {', '.join(codes)}" for code, codes in missing.items())
                + "."
            )
        zero_weight = self.provenance.get("zero_weight_regions")
        if zero_weight:
            lines.append(
                "- Regions NSS located at 0% of the basin (an overlap under half a percent), "
                "dropped: " + ", ".join(zero_weight) + "."
            )
        return "\n".join(lines)


def _nss_served(data_dir: Union[str, Path, None]) -> Set[str]:
    """Jurisdictions NSS serves peak-flow regions for: a published or snapshot file.

    Computed at call time from the library directory, so it follows the packaged data
    (or a ``data_dir`` override) rather than a list that would go stale.
    """
    return set(available_states(data_dir)) | set(nss_snapshot_states(data_dir))


def _load_library(
    code: str,
    data_dir: Union[str, Path, None],
    *,
    allow_nss_snapshot: bool,
    allow_provisional: bool,
) -> StateLibrary:
    """``load_state`` whose no-file error says whether NSS serves the jurisdiction."""
    try:
        return load_state(
            code,
            data_dir,
            allow_nss_snapshot=allow_nss_snapshot,
            allow_provisional=allow_provisional,
        )
    except EquationsUnavailable as err:
        served = code in _nss_served(data_dir)
        message = (
            f"{code}: no regression library file. NSS {'serves' if served else 'serves no'} "
            "peak-flow regions for this jurisdiction"
        )
        # A served jurisdiction without a published file has only a snapshot, which
        # load_state returns on opt-in alone: say so.
        if code in nss_snapshot_states(data_dir):
            message += (
                "; its harvested NSS set is available with allow_nss_snapshot=True "
                "(not report-verified)"
            )
        raise EquationsUnavailable(message) from err


def estimate_at_site(
    lat: float,
    lon: float,
    *,
    state: str,
    aeps: Optional[Sequence[float]] = None,
    at_site: Optional[FrequencyResults] = None,
    allow_nss_snapshot: bool = False,
    allow_provisional: bool = False,
    allow_extrapolation: bool = False,
    weighting_space: str = "log",
    cache: Optional[StreamStatsCache] = None,
    timeout: float = 60.0,
    data_dir: Union[str, Path, None] = None,
    region: Union[str, Sequence[str], None] = None,
) -> SiteEstimate:
    """Estimate flood quantiles at an ungaged site from the state's regression equations.

    Loads the jurisdiction's equation library, delineates the basin and fetches its
    characteristics from StreamStats
    (:func:`~flowfreq.streamstats.delineate_and_get_characteristics`), locates the NSS
    regression regions the basin overlaps
    (:func:`~flowfreq.streamstats.locate_regression_regions`), and hands the three to
    :func:`assemble`. The library is loaded first, so a jurisdiction without one is
    refused before any network call.

    Parameters
    ----------
    lat, lon : float
        The pour point, decimal degrees WGS84. StreamStats snaps it to the stream
        network; the snap is recorded in ``provenance["snap"]``.
    state : str
        Two-letter jurisdiction code: the StreamStats region and the equation library.
        Never inferred from the coordinate.
    aeps : sequence of float, optional
        AEPs to estimate, in the caller's order; ``None`` means every AEP all the regions
        used publish, largest first.
    at_site : FrequencyResults, optional
        A Bulletin 17C fit at the site. When the basin lies in a single regression
        region, it is weighted with the regression by inverse variance
        (:func:`~flowfreq.gage_weighting.weight_with_regression`, Bulletin 17C
        Appendix 9) and the result carries the weighted frame; a basin spanning several
        regions has an area-weighted flow with no variance of prediction, so the fit is
        returned unweighted and ``provenance["weighting"]`` says so.
    allow_nss_snapshot, allow_provisional : bool
        Opt in to the jurisdiction's ``nss_snapshot`` or ``provisional`` equation set
        (:func:`~flowfreq.regression.library.load_state`). The set's advisory is
        returned in ``advisory``.
    allow_extrapolation : bool
        Evaluate outside an equation's calibrated range and flag the row, instead of
        raising :class:`~flowfreq.regression.equations.OutOfRangeError`. Also lets an
        extrapolated regression take part in the at-site weighting.
    weighting_space : {"log", "linear"}
        How a multi-region basin's estimates are area-averaged.
    cache : StreamStatsCache, optional
        Passed to the delineation.
    timeout : float
        Per-request timeout, seconds, for both StreamStats calls.
    data_dir : str or Path, optional
        Override the packaged library directory (a test and override hook, passed to
        :func:`~flowfreq.regression.library.load_state`; the NSS-coverage check in the
        no-file error reads the same directory).
    region : str, optional
        Where NSS returns overlapping regions (alternative equation sets for the same
        area, such as a superseded report and its replacement), the one to use, by NSS
        code or library ``region_code`` (:func:`assemble`). The others are listed in
        ``provenance["alternative_regions"]``. May also be a sequence of codes that
        together cover the basin (see :func:`assemble`). Where NSS locates no peak-flow
        region for the state (Wyoming, Nevada), a library region code chosen from the
        report's map is used alone.

    Returns
    -------
    SiteEstimate
        ``method`` is ``"at-site weighted with regression"`` when the fit was weighted
        in, else ``"regression"``. ``provenance`` carries ``snap``, ``delineation``
        (polygon area, unavailable characteristics, the StreamStats provenance),
        ``library`` (status, file, retrieval date), ``regions_located`` (every region
        NSS returned, before matching), ``missing_characteristics``,
        ``zero_weight_regions``, ``alternative_regions`` (the overlapping regions
        set aside by ``region``, with citation and year), ``out_of_range`` (AEP to
        NSS code to the inputs outside their range), ``weighting``, ``options`` and
        ``estimated_at_utc``.

    Raises
    ------
    OverlappingRegionsError
        The located regions overlap and ``region`` does not choose one of them; the
        message and ``candidates`` list them.
    EquationsUnavailable
        If the jurisdiction has no library file; the message says whether NSS serves
        peak-flow regions for it at all (DC and four territories have none). Also the
        :func:`assemble` cases: no located region in the library, every region dropped,
        no AEP in common, or a requested AEP not published.
    OutOfRangeError
        An input outside its calibrated range without ``allow_extrapolation``.
    ValueError
        An unknown ``weighting_space``; an empty ``aeps``; a weighted ``Q=0`` in log space; or an
        ``at_site`` whose ``confidence_limits`` lack ``variance_log``.
    UnsnappablePointError, DegenerateDelineationError, StreamStatsResponseError, \
    StreamStatsTransportError
        From StreamStats, unchanged.
    """
    if weighting_space not in ("log", "linear"):
        raise ValueError(f"weighting_space must be 'log' or 'linear', got {weighting_space!r}")
    if aeps is not None and len(aeps) == 0:
        raise ValueError("aeps must be non-empty (use None for the regions' common AEPs)")
    if at_site is not None:
        # weight_with_regression's own refusal, raised before the StreamStats calls.
        require_variance_log(at_site)
    code = jurisdiction(state).code
    lib = _load_library(
        code,
        data_dir,
        allow_nss_snapshot=allow_nss_snapshot,
        allow_provisional=allow_provisional,
    )
    ws = delineate_and_get_characteristics(code, lat, lon, cache=cache, timeout=timeout)
    if ws.polygon_geojson is None:
        # The delineation is asked for its polygon; a cache entry saved before polygons
        # were returned is the one way to arrive here without one.
        raise ValueError(
            f"{code}: StreamStats returned no watershed polygon at {lat}, {lon}; the NSS "
            "regions cannot be located without one"
        )
    located = locate_regression_regions(code, ws.polygon_geojson, timeout=timeout)
    characteristics = {k: c.value for k, c in ws.characteristics.items()}

    assembly = assemble(
        lib,
        located,
        characteristics,
        aeps=aeps,
        allow_extrapolation=allow_extrapolation,
        weighting_space=weighting_space,
        region=region,
    )

    method = "regression"
    weighted: Optional[pd.DataFrame] = None
    if at_site is None:
        weighting = "no at-site fit supplied"
    elif len(assembly.used) == 1:
        regression = [per_region[0][1] for per_region in assembly.estimates.values()]
        weighted = weight_with_regression(
            at_site, regression, allow_extrapolated=allow_extrapolation
        )
        method = "at-site weighted with regression"
        weighting = f"at-site fit weighted with regression region {assembly.used[0].code}"
    else:
        weighting = (
            f"at-site fit not weighted: the basin spans {len(assembly.used)} regression "
            f"regions ({_codes(assembly.used)}), whose area-weighted flow carries no "
            "variance of prediction"
        )

    out_of_range: Dict[float, Dict[str, List[str]]] = {}
    for aep, per_region in assembly.estimates.items():
        for used_region, (_, est) in zip(assembly.used, per_region):
            if est.out_of_range:
                out_of_range.setdefault(aep, {})[used_region.code] = list(est.out_of_range)

    subdir = OPT_IN_SUBDIRS.get(lib.status)
    provenance: Dict[str, Any] = {
        "snap": ws.snap.to_dict(),
        "delineation": {
            "polygon_area_sq_mi": ws.polygon_area_sq_mi,
            "unavailable": dict(ws.unavailable),
            "streamstats": ws.provenance.to_dict() if ws.provenance else None,
        },
        "library": {
            "status": lib.status,
            "file": f"{subdir}/{code}.json" if subdir else f"{code}.json",
            "retrieved_utc": lib.retrieved_utc,
        },
        "regions_located": [r.to_dict() for r in located],
        "missing_characteristics": assembly.missing,
        "zero_weight_regions": [r.code for r in assembly.zero_weight],
        "alternative_regions": [
            _alternative(lib, r, assembly.region_codes[r.code]) for r in assembly.alternatives
        ],
        "out_of_range": out_of_range,
        "weighting": weighting,
        "options": {
            "allow_nss_snapshot": allow_nss_snapshot,
            "allow_provisional": allow_provisional,
            "allow_extrapolation": allow_extrapolation,
            "weighting_space": weighting_space,
            "region": region,
        },
        "estimated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    if assembly.region_choice_source is not None:
        provenance["region_choice_source"] = assembly.region_choice_source
    return SiteEstimate(
        lat=lat,
        lon=lon,
        state=code,
        method=method,
        quantiles=assembly.quantiles,
        regions=assembly.used,
        unmatched_regions=assembly.unmatched,
        characteristics=characteristics,
        library_status=lib.status,
        advisory=lib.advisory,
        estimates=assembly.estimates,
        region_codes=assembly.region_codes,
        weighted=weighted,
        provenance=provenance,
    )


def _alternative(lib: StateLibrary, r: RegressionRegionLocation, lib_code: str) -> Dict[str, Any]:
    """A provenance entry for an overlapping region set aside."""
    citation, year = _region_citation(lib, lib_code)
    return {
        "code": r.code,
        "name": r.name,
        "region_code": lib_code,
        "citation": citation,
        "year": year,
    }


def _cfs(value: float) -> str:
    return "" if pd.isna(value) else f"{value:,.0f}"


def _level(value: float) -> str:
    return "" if pd.isna(value) else f"{value:.0%}"


def _frame_table(frame: pd.DataFrame) -> List[str]:
    """A DataFrame as Markdown table rows, without the optional tabulate dependency."""
    cols = [str(c) for c in frame.columns]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(_cell(v) for v in row) + " |")
    return lines


def _cell(value: Any) -> str:
    if isinstance(value, (int, float)) and not isinstance(value, bool) and not pd.isna(value):
        return f"{value:,.6g}"
    return str(value)
