"""Flood-frequency estimates at ungaged sites from published regression equations.

This module holds :class:`SiteEstimate`, the result of an estimate at an ungaged
location, and :func:`assemble`, the pure assembly of that estimate from inputs a
caller has already fetched: the located NSS regression regions, the basin
characteristics, and a state's equation library. Nothing here touches the network;
the delineation, region location and library loading that feed it belong to the
``estimate_at_site`` entry point built on top.

Assembly rules
--------------
- Each located region is mapped to a library region with
  :func:`region_for_nss_code`. A region with no library match is reported in
  ``unmatched`` and the estimate proceeds on the rest; if none matches,
  :class:`~flowfreq.regression.library.EquationsUnavailable` names the codes.
- A matched region with ``percent_weight`` 0 (NSS rounds to a whole percent, so a sliver
  under 0.5% arrives as 0) contributes nothing and is dropped before anything is
  evaluated, so a basin effectively in one region keeps that region's interval.
- A matched region whose equations (at any AEP it publishes) need a characteristic the
  basin does not carry is skipped and the missing codes are reported by its NSS code.
- The default AEPs are those every region that remains publishes, largest first; none in
  common raises :class:`~flowfreq.regression.library.EquationsUnavailable`.
- Area fractions come from NSS's ``percent_weight`` (what NSS itself weights by),
  renormalised over the regions actually used.
- One region gives its own prediction interval; several give an area-weighted flow with
  no interval and ``interval_method="none: area-weighted across N regions"``.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from flowfreq.regression.equations import RegressionEstimate, evaluate, evaluate_weighted
from flowfreq.regression.library import EquationsUnavailable, StateLibrary
from flowfreq.streamstats import RegressionRegionLocation

logger = logging.getLogger(__name__)

__all__ = ["Assembly", "SiteEstimate", "assemble", "region_for_nss_code"]

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
    """

    quantiles: pd.DataFrame
    estimates: Dict[float, RegionEstimates]
    used: List[RegressionRegionLocation]
    region_codes: Dict[str, str]
    unmatched: List[RegressionRegionLocation]
    zero_weight: List[RegressionRegionLocation]
    missing: Dict[str, List[str]]


def assemble(
    lib: StateLibrary,
    regions: Sequence[RegressionRegionLocation],
    characteristics: Mapping[str, float],
    *,
    aeps: Optional[Sequence[float]],
    allow_extrapolation: bool,
    weighting_space: str,
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
        regions used publish, largest first.
    allow_extrapolation : bool
        Evaluate outside an equation's calibrated range and flag the row, instead of
        raising :class:`~flowfreq.regression.equations.OutOfRangeError`.
    weighting_space : {"log", "linear"}
        How a multi-region basin's estimates are area-averaged
        (:func:`~flowfreq.regression.equations.evaluate_weighted`).

    Returns
    -------
    Assembly

    Raises
    ------
    EquationsUnavailable
        If no located region matches the library; if every matched region has zero
        weight or lacks a basin characteristic; if the regions used publish no AEP in
        common and none was requested; or if a requested AEP is not published for a
        region used.
    OutOfRangeError
        If an input is outside its range and ``allow_extrapolation`` is False.
    ValueError
        From :func:`evaluate_weighted`: an unknown ``weighting_space``, or a printed
        ``Q=0`` weighted in log space.
    """
    matched: List[Tuple[RegressionRegionLocation, str]] = []
    unmatched: List[RegressionRegionLocation] = []
    for region in regions:
        lib_code = region_for_nss_code(lib, region.code)
        if lib_code is None:
            unmatched.append(region)
        else:
            matched.append((region, lib_code))
    if not matched:
        raise EquationsUnavailable(
            f"{lib.state}: none of the located NSS regions ({_codes(regions)}) matches a "
            f"library region ({', '.join(lib.regions) or 'none'})"
        )
    if unmatched:
        logger.warning(
            "%s: no library region for located NSS region(s) %s; estimating on %s",
            lib.state,
            _codes(unmatched),
            ", ".join(r.code for r, _ in matched),
        )
    region_codes = {region.code: lib_code for region, lib_code in matched}

    zero_weight = [region for region, _ in matched if not region.percent_weight > 0]
    weighted = [(region, lib_code) for region, lib_code in matched if region.percent_weight > 0]
    if not weighted:
        raise EquationsUnavailable(
            f"{lib.state}: every matched NSS region ({_codes(zero_weight)}) has zero weight"
        )
    if zero_weight:
        logger.info("%s: dropping zero-weight NSS region(s) %s", lib.state, _codes(zero_weight))

    missing: Dict[str, List[str]] = {}
    used: List[Tuple[RegressionRegionLocation, str]] = []
    for region, lib_code in weighted:
        # Every equation the region publishes, whatever AEPs were requested: a region
        # that lacks a characteristic is skipped before an unpublished requested AEP
        # could be raised for it.
        check = lib.aeps(lib_code)
        absent = [c for c in _needed_codes(lib, lib_code, check) if c not in characteristics]
        if absent:
            missing[region.code] = absent
            logger.warning(
                "%s: skipping region %s, basin characteristics missing %s",
                lib.state,
                region.code,
                absent,
            )
        else:
            used.append((region, lib_code))
    if not used:
        raise EquationsUnavailable(
            f"{lib.state}: every matched region needs a basin characteristic the site "
            "lacks: "
            + "; ".join(f"{code} needs {', '.join(codes)}" for code, codes in missing.items())
        )

    if aeps is None:
        common = set.intersection(*(set(lib.aeps(lib_code)) for _, lib_code in used))
        aep_list = sorted(common, reverse=True)
        if not aep_list:
            raise EquationsUnavailable(
                f"{lib.state}: the regions used publish no AEP in common: "
                + "; ".join(
                    f"{region.code} at {', '.join(f'{a:g}' for a in lib.aeps(lib_code))}"
                    for region, lib_code in used
                )
            )
    else:
        aep_list = [float(a) for a in aeps]

    total = sum(region.percent_weight for region, _ in used)
    fractions = [region.percent_weight / total for region, _ in used]
    region_label = ", ".join(region.code for region, _ in used)

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
        used=[region for region, _ in used],
        region_codes=region_codes,
        unmatched=unmatched,
        zero_weight=zero_weight,
        missing=missing,
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
        the dropped-region caveat.
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
                    "| AEP | Flow (cfs) | Lower | Upper | Level | Interval | Extrapolated | Regions |",
                    "|---:|---:|---:|---:|---:|---|---|---|",
                ]
            )
            for _, row in self.quantiles.iterrows():
                lines.append(
                    f"| {row['aep']:g} | {row['flow_cfs']:,.0f} | {_cfs(row['lower'])} | "
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
