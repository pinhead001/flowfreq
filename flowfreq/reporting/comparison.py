"""The ``comparison`` section: several estimates at one site, side by side."""

from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from ..aep import same_aep
from ._base import SectionContent, _demote, _version, cfs
from .blocks import Block, Heading, Paragraph, Table
from .gaged import _gaged_ffa
from .regional import _ungaged
from .spec import ReportError
from .transposed import _transposed
from .weighted import _nested

# -- several estimates of the flood quantiles at one site, side by side --------------------

#: The kinds of estimate a ``comparison`` section can set side by side.
_COMPARED_KINDS = ("gaged_ffa", "transposed", "ungaged")


_ESTIMATE_KEYS = ("label", "kind", "options")


#: The caption of a ``comparison`` section's main table.
COMPARISON_CAPTION = "Flood quantiles at the site by each method."


def _compared(
    index: int, given: Any, section: str = "comparison"
) -> Tuple[str, SectionContent, pd.DataFrame, str]:
    """One estimate: its label, its section, its flows with limits, and what the limits are."""
    owner = f"{section}: estimate {index + 1}"
    if not isinstance(given, dict):
        raise ReportError(f"{owner} must be an object with label, kind and options")
    extra = sorted(set(given) - set(_ESTIMATE_KEYS))
    if extra:
        raise ReportError(
            f"{owner}: unknown key(s) {', '.join(extra)}; allowed: {', '.join(_ESTIMATE_KEYS)}"
        )
    label = str(given.get("label") or "").strip()
    kind = given.get("kind")
    if not label:
        raise ReportError(f"{owner} needs a label, the name its column and subsection take")
    if kind not in _COMPARED_KINDS:
        raise ReportError(
            f"{owner}: kind must be one of {', '.join(_COMPARED_KINDS)}; got {kind!r}"
        )
    options = _nested(str(kind), "options", given.get("options"), owner)
    if kind == "gaged_ffa":
        content, b17c = _gaged_ffa(options)
        flows = b17c.results.quantiles[["aep", "flow_cfs"]].merge(
            b17c.results.confidence_limits[["aep", "lower_5pct", "upper_5pct"]],
            on="aep",
            how="left",
        )
        flows = flows.rename(columns={"lower_5pct": "lower", "upper_5pct": "upper"})
        return label, content, flows, "90-percent confidence limits"
    if kind == "transposed":
        content, results = _transposed(options)
        flows = results.quantiles[["aep", "flow_cfs"]].copy()
        limits = results.confidence_limits
        if limits.empty:
            flows["lower"], flows["upper"] = np.nan, np.nan
        else:
            flows = flows.merge(limits[["aep", "lower_5pct", "upper_5pct"]], on="aep", how="left")
            flows = flows.rename(columns={"lower_5pct": "lower", "upper_5pct": "upper"})
        return label, content, flows, "the donor's 90-percent confidence limits, scaled"
    content, estimate = _ungaged(options)
    flows = estimate.quantiles[["aep", "flow_cfs", "lower", "upper"]].copy()
    levels = sorted({float(v) for v in estimate.quantiles["interval_level"].dropna()})
    said = (
        f"{levels[0]:.0%} prediction limits"
        if len(levels) == 1
        else "prediction limits" if levels else "no limits"
    )
    return label, content, flows, said.replace("%", "-percent")


def _flow_table(flows: pd.DataFrame) -> Dict[float, Any]:
    """An estimate's rows that have a flow, by AEP."""
    return {float(r["aep"]): r for _, r in flows[flows["flow_cfs"].notna()].iterrows()}


def _row_at(table: Mapping[float, Any], aep: float) -> Any:
    """The estimate's row at ``aep``, or at the AEP nearest it within ``AEP_MATCH``; else None."""
    if not table:
        return None
    nearest = min(table, key=lambda a: abs(a - aep))
    return table[nearest] if same_aep(nearest, aep) else None


def _matched_aeps(
    owner: str,
    tables: Sequence[Mapping[float, Any]],
    asked: Any,
    default: Optional[Sequence[float]] = None,
) -> Tuple[List[float], List[str]]:
    """The AEPs to tabulate across several estimates, and which of them join unlike AEPs.

    ``asked`` is the section's ``aeps`` option. Without it the AEPs are those of ``default``
    that every estimate has, or, with no default, every AEP all the estimates have. An AEP
    is "had" within :data:`AEP_MATCH`. The second value lists, in words, the rows whose
    estimates sit at slightly different AEPs.
    """
    if asked is None:
        if default is None:
            candidates: List[float] = []
            for aep in sorted({a for table in tables for a in table}, reverse=True):
                if not any(same_aep(aep, c) for c in candidates):
                    candidates.append(aep)
        else:
            candidates = [float(a) for a in default]
        aeps = [a for a in candidates if all(_row_at(t, a) is not None for t in tables)]
        if not aeps:
            raise ReportError(f"{owner}: the estimates have no AEP in common")
    else:
        if not isinstance(asked, list) or not asked:
            raise ReportError(f"{owner}: aeps must be a non-empty list of AEPs")
        try:
            aeps = [float(a) for a in asked]
        except (TypeError, ValueError) as exc:
            raise ReportError(f"{owner}: aeps must be a list of numbers") from exc
        missing = [a for a in aeps if any(_row_at(t, a) is None for t in tables)]
        if missing:
            raise ReportError(
                f"{owner}: not every estimate has a quantile at AEP "
                + ", ".join(f"{a:g}" for a in missing)
            )
    shown: List[float] = []
    apart: List[str] = []
    for aep in aeps:
        used = sorted({float(_row_at(table, aep)["aep"]) for table in tables})
        # Where the estimates differ, the row takes the value written with the most digits:
        # 0.667, not 0.67, for the 1.5-year flood.
        label = (
            aep
            if asked is not None or default is not None
            else max(used, key=lambda a: len(format(a, "g")))
        )
        shown.append(label)
        if len(used) > 1:
            apart.append(
                format(label * 100, ".3g")
                + " percent ("
                + ", ".join(format(value, "g") for value in used)
                + ")"
            )
    return shown, apart


def _build_comparison(given: Dict[str, Any]) -> SectionContent:
    options = dict(given)
    listed = options["estimates"]
    if not isinstance(listed, list) or len(listed) < 2:
        raise ReportError("comparison: estimates must be a list of two or more estimates")
    built = [_compared(i, each) for i, each in enumerate(listed)]
    labels = [label for label, _, _, _ in built]
    if len(set(labels)) != len(labels):
        raise ReportError("comparison: two estimates have the same label")
    tables = [_flow_table(flows) for _, _, flows, _ in built]
    aeps, apart = _matched_aeps("comparison", tables, options.get("aeps"))

    site = str(options["site_name"]).strip() if options.get("site_name") else ""
    where = site or "the site"
    flow_rows = []
    limit_rows = []
    for aep in aeps:
        matched = [_row_at(table, aep) for table in tables]
        flows = [float(row["flow_cfs"]) for row in matched]
        low, high = min(flows), max(flows)
        flow_rows.append(
            (
                f"{aep * 100:.3g}",
                *(cfs(v) for v in flows),
                f"{high / low:.2f}" if low > 0 else "-",
            )
        )
        cells = [f"{aep * 100:.3g}"]
        for row in matched:
            lower, upper = row["lower"], row["upper"]
            cells.append(
                "-" if pd.isna(lower) or pd.isna(upper) else f"{cfs(lower)} to {cfs(upper)}"
            )
        limit_rows.append(tuple(cells))
    blocks: List[Block] = [
        Paragraph(
            f"The flood quantiles at {where} were estimated by {len(built)} methods, set side by "
            "side below at the annual exceedance probabilities all of them give. Each "
            "estimate is computed as its own subsection describes"
            + ("." if bool(options.get("detail", True)) else ", in this report's other sections.")
            + " The last column is the largest estimate over the smallest. The estimates "
            "are compared, not combined: no weighting is applied here, and they are not "
            "independent where they share a gage record or a regression. The computation "
            f"was made with flowfreq {_version()}."
        ),
        Table(
            (
                "Annual exceedance probability (percent)",
                *(f"{label} (cfs)" for label in labels),
                "Largest over smallest",
            ),
            tuple(flow_rows),
            ("right",) * (len(labels) + 2),
            COMPARISON_CAPTION,
        ),
    ]
    if apart:
        blocks.append(
            Paragraph(
                "The methods tabulate some rows at slightly different probabilities, which "
                "are shown together: " + "; ".join(apart) + "."
            )
        )
    if any(cell != "-" for row in limit_rows for cell in row[1:]):
        blocks.append(
            Paragraph(
                "Limits, where a method gives them: "
                + "; ".join(f"{label}, {said}" for label, _, _, said in built)
                + ". A transposed estimate's limits carry the donor gage's sampling "
                "uncertainty only, not the uncertainty of the transposition."
            )
        )
        blocks.append(
            Table(
                (
                    "Annual exceedance probability (percent)",
                    *(f"{label} (cfs)" for label in labels),
                ),
                tuple(limit_rows),
                ("right",) * (len(labels) + 1),
                "Lower and upper limits of each estimate.",
            )
        )
    references: List[str] = []
    resolved = []
    for (label, content, _, _), each in zip(built, listed):
        if bool(options.get("detail", True)):
            blocks.append(Heading(2, label))
            blocks.extend(_demote(content.blocks))
        references.extend(r for r in content.references if r not in references)
        resolved.append({"label": label, "kind": each["kind"], "options": content.options})
    options["estimates"] = resolved
    return SectionContent(
        title=(
            f"Flood quantiles at {site} by {len(built)} methods"
            if site
            else "Comparison of flood estimates"
        ),
        blocks=blocks,
        options=options,
        references=references,
    )
