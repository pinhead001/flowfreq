"""The ``comparison`` section: several estimates at one site, side by side."""

from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from ..aep import same_aep
from ._base import SUPERSEDED, SectionContent, _curve_figure, _demote, _version, cfs, with_label
from .b17b import _b17b_ffa
from .blocks import Block, Figure, Heading, Paragraph, Table
from .gaged import _gaged_ffa
from .regional import _caution_paragraphs, _ungaged
from .spec import ReportError
from .transposed import _transposed
from .weighted import _nested

# -- several estimates of the flood quantiles at one site, side by side --------------------

#: The kinds of estimate a ``comparison`` section can set side by side.
_COMPARED_KINDS = ("gaged_ffa", "transposed", "ungaged", "b17b_ffa")


#: Why a ``site_summary`` section takes no ``b17b_ffa`` estimate.
B17B_NOT_SUMMARIZED = (
    f"{SUPERSEDED} is shown in a b17b_ffa section or a comparison, not in a summary"
)


_ESTIMATE_KEYS = ("label", "kind", "options")


#: The caption of a ``comparison`` section's main table.
COMPARISON_CAPTION = "Flood quantiles at the site by each method."


def _compared(
    index: int,
    given: Any,
    section: str = "comparison",
    kinds: Sequence[str] = _COMPARED_KINDS,
) -> Tuple[str, SectionContent, pd.DataFrame, str]:
    """One estimate: its label, its section, its flows with limits, and what the limits are.

    ``kinds`` are those the section takes. The flows are empty for a ``b17b_ffa`` estimate
    when Bulletin 17B gives no answer.
    """
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
    if kind not in kinds:
        message = f"{owner}: kind must be one of {', '.join(kinds)}; got {kind!r}"
        if kind == "b17b_ffa":
            message += f". {B17B_NOT_SUMMARIZED}"
        raise ReportError(message)
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
    if kind == "b17b_ffa":
        content, fit = _b17b_ffa(options)
        columns = ["aep", "flow_cfs", "lower", "upper"]
        flows = pd.DataFrame(columns=columns) if fit is None else fit.quantiles[columns].copy()
        return label, content, flows, "90-percent confidence limits (Bulletin 17B)"
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


def _b17b_label(label: str) -> str:
    """A ``b17b_ffa`` estimate's label as shown: it says that Bulletin 17B is superseded.

    A label that names Bulletin 17B gets the one word it lacks, any other the whole phrase.
    """
    words = label.lower()
    if "bulletin 17b" not in words:
        return with_label(label, SUPERSEDED)
    return label if "superseded" in words else f"{label} (superseded)"


def _comparison_label(options: Mapping[str, Any]) -> str:
    """What every title of the section has to say: that a Bulletin 17B estimate is in it."""
    present = any(each["kind"] == "b17b_ffa" for each in options["estimates"])
    return f"with {SUPERSEDED}" if present else ""


def _within(flows: pd.DataFrame, aeps: Sequence[float]) -> pd.DataFrame:
    """An estimate's rows from the most common of ``aeps`` to the rarest."""
    top, bottom = max(aeps), min(aeps)
    keep = [
        bottom <= a <= top or same_aep(a, top) or same_aep(a, bottom)
        for a in flows["aep"].astype(float)
    ]
    return flows[keep]


def _curves_block(
    estimates: Sequence[Tuple[str, str, pd.DataFrame, str]], aeps: Sequence[float], where: str
) -> Figure:
    """The curves of ``(label, kind, flows, limits)`` estimates; those with no flows are named."""
    shown = [e for e in estimates if not e[2].empty]
    kinds = {kind for _, kind, _, _ in shown}
    named = []
    if kinds & {"gaged_ffa", "transposed"}:
        named.append("Bulletin 17C")
    if "b17b_ffa" in kinds:
        named.append(SUPERSEDED)
    caption = " and ".join(named) + " frequency curves." if named else "Frequency curves."
    superseded = [label for label, kind, _, _ in shown if kind == "b17b_ffa"]
    if superseded:
        caption += f" Superseded, dashed: {'; '.join(superseded)}."
    missing = [label for label, _, flows, _ in estimates if flows.empty]
    if missing:
        caption += f" Not drawn: {'; '.join(missing)}."
    limits = "; ".join(f"{label}, {said}" for label, _, _, said in shown)
    caption += f" Shaded bands are each estimate's limits ({limits})."
    png = _curve_figure(
        [(label, _within(flows, aeps)) for label, _, flows, _ in shown],
        aeps,
        where,
        dashed=superseded,
        legend_below=True,
    )
    return Figure(png, caption)


def _opening(where: str, n: int, n_refused: int, detail: bool, has_b17b: bool) -> str:
    """The paragraph above the table; with every estimate refused there is no table."""
    if n_refused == n:
        return (
            f"The flood quantiles at {where} were asked for by {n} methods, and none of them "
            "gives an answer, so there is no table. The reasons are below. The computation "
            f"was made with flowfreq {_version()}."
        )
    if n_refused:
        said = (
            f"The flood quantiles at {where} were asked for by {n} methods, and {n_refused} of "
            f"them {'gives' if n_refused == 1 else 'give'} no answer. The rest are set side by "
            "side below at the annual exceedance probabilities all of those give. Each "
        )
    else:
        said = (
            f"The flood quantiles at {where} were estimated by {n} methods, set side by "
            "side below at the annual exceedance probabilities all of them give. Each "
        )
    said += (
        "estimate is computed as its own subsection describes"
        + ("." if detail else ", in this report's other sections.")
        + " The last column is the largest estimate over the smallest. The estimates "
        "are compared, not combined: no weighting is applied here, and they are not "
        "independent where they share a gage record or a regression. The computation "
        f"was made with flowfreq {_version()}."
    )
    if has_b17b:
        said += (
            " Bulletin 17B was superseded by Bulletin 17C; its estimates are here for "
            "comparison only."
        )
    return said


def _comparison_tables(
    built: Sequence[Tuple[str, SectionContent, pd.DataFrame, str]],
    tables: Sequence[Mapping[float, Any]],
    refused: Sequence[bool],
    aeps: Sequence[float],
    apart: Sequence[str],
) -> List[Block]:
    """The table of flows, the table of limits and what they need said.

    An estimate that is ``refused`` has a column of dashes and no part in the ratio.
    """
    labels = [label for label, _, _, _ in built]
    flow_rows = []
    limit_rows = []
    for aep in aeps:
        matched = [_row_at(table, aep) for table in tables]
        flows = [float(row["flow_cfs"]) for row in matched if row is not None]
        low, high = min(flows), max(flows)
        # One flow left by the refusals is its own largest and smallest.
        compared = low > 0 and (len(flows) > 1 or not any(refused))
        flow_rows.append(
            (
                f"{aep * 100:.3g}",
                *("-" if row is None else cfs(float(row["flow_cfs"])) for row in matched),
                f"{high / low:.2f}" if compared else "-",
            )
        )
        cells = [f"{aep * 100:.3g}"]
        for row in matched:
            if row is None or pd.isna(row["lower"]) or pd.isna(row["upper"]):
                cells.append("-")
            else:
                cells.append(f"{cfs(row['lower'])} to {cfs(row['upper'])}")
        limit_rows.append(tuple(cells))
    blocks: List[Block] = [
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
                + "; ".join(
                    f"{label}, {said}"
                    for (label, _, _, said), no_answer in zip(built, refused)
                    if not no_answer
                )
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
    return blocks


def _build_comparison(given: Dict[str, Any]) -> SectionContent:
    options = dict(given)
    listed = options["estimates"]
    if not isinstance(listed, list) or len(listed) < 2:
        raise ReportError("comparison: estimates must be a list of two or more estimates")
    built = []
    for i, each in enumerate(listed):
        label, content, estimate, said = _compared(i, each)
        shown = _b17b_label(label) if each["kind"] == "b17b_ffa" else label
        built.append((shown, content, estimate, said))
    kinds = [each["kind"] for each in listed]
    labels = [label for label, _, _, _ in built]
    if len(set(labels)) != len(labels):
        raise ReportError("comparison: two estimates have the same label")
    detail = bool(options.get("detail", True))
    has_b17b = "b17b_ffa" in kinds
    # Bulletin 17B refuses some records: such an estimate has no flows to match or compare.
    refused = [k == "b17b_ffa" and e.empty for k, (_, _, e, _) in zip(kinds, built)]
    tables = [_flow_table(flows) for _, _, flows, _ in built]
    answered = [table for table, no_answer in zip(tables, refused) if not no_answer]
    asked = options.get("aeps")
    aeps: List[float] = []
    apart: List[str] = []
    if answered or asked is not None:
        aeps, apart = _matched_aeps("comparison", answered, asked)

    site = str(options["site_name"]).strip() if options.get("site_name") else ""
    where = site or "the site"
    blocks: List[Block] = [Paragraph(_opening(where, len(built), sum(refused), detail, has_b17b))]
    if answered:
        blocks.extend(_comparison_tables(built, tables, refused, aeps, apart))
    if not detail:
        blocks.extend(
            Paragraph(f"{label}: {text}")
            for (label, content, _, _), no_answer in zip(built, refused)
            if no_answer
            for text in content.cautions
        )
    if has_b17b and answered:
        curves = [(label, kind, e, said) for (label, _, e, said), kind in zip(built, kinds)]
        blocks.append(_curves_block(curves, aeps, where))
    references: List[str] = []
    resolved = []
    for (label, content, _, _), each, no_answer in zip(built, listed, refused):
        if detail:
            blocks.append(Heading(2, label))
            blocks.extend(_demote(content.blocks))
        elif not no_answer:
            blocks.extend(_caution_paragraphs(label, content))
        references.extend(r for r in content.references if r not in references)
        # The caller's label: the words Bulletin 17B adds to it are added again.
        own = str(each["label"]).strip()
        resolved.append({"label": own, "kind": each["kind"], "options": content.options})
    options["estimates"] = resolved
    title = (
        f"Flood quantiles at {site} by {len(built)} methods"
        if site
        else "Comparison of flood estimates"
    )
    return SectionContent(
        title=with_label(title, _comparison_label(options)),
        blocks=blocks,
        options=options,
        references=references,
        cautions=[f"{label}: {text}" for label, c, _, _ in built for text in c.cautions],
    )
