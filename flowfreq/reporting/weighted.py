"""The ``weighted`` section: an at-site fit weighted with the regression."""

from __future__ import annotations

import math
from typing import Any, Dict, List

from ._base import SECTION_KINDS, SectionContent, _demote, _version, cfs
from .blocks import Block, Heading, Paragraph, Table
from .gaged import _gaged_ffa
from .regional import _ungaged
from .spec import ReportError

# -- at-site fit weighted with the regression ---------------------------------------------


def _nested(kind: str, name: str, value: Any, owner: str = "weighted") -> Dict[str, Any]:
    """A nested section's options, checked against the kind they belong to."""
    if not isinstance(value, dict):
        raise ReportError(f"{owner}: {name} must be an object of {kind} options")
    SECTION_KINDS[kind].check(value)
    return dict(value)


def _build_weighted(given: Dict[str, Any]) -> SectionContent:
    from ..gage_weighting import weight_with_regression

    options = dict(given)
    gage, b17c = _gaged_ffa(_nested("gaged_ffa", "gage", options["gage"]))
    regression, estimate = _ungaged(_nested("ungaged", "regression", options["regression"]))
    if len(estimate.regions) != 1:
        raise ReportError(
            f"weighted: the basin spans {len(estimate.regions)} regression regions, whose "
            "area-weighted flow has no variance of prediction, so it cannot be weighted with "
            "the at-site fit"
        )
    extrapolate = bool(options["regression"].get("allow_extrapolation", False))
    frame = weight_with_regression(
        b17c.results,
        [per_region[0][1] for per_region in estimate.estimates.values()],
        allow_extrapolated=extrapolate,
    )
    # The at-site fit has more AEPs than any regression publishes; those rows say nothing.
    frame = frame[frame["regression_flow_cfs"].notna()]
    if frame.empty:
        raise ReportError("weighted: the regression and the at-site fit share no AEP")
    if frame["weighted_flow_cfs"].isna().all():
        raise ReportError(
            "weighted: no AEP could be weighted: "
            + "; ".join(sorted({str(note) for note in frame["note"] if note}))
        )
    noted = bool(frame["note"].astype(bool).any())
    header = [
        "Annual exceedance probability (percent)",
        "At-site flow (cfs)",
        "Regression flow (cfs)",
        "Weight on the at-site flow",
        "Weighted flow (cfs)",
        "Lower 90-percent limit (cfs)",
        "Upper 90-percent limit (cfs)",
    ]
    rows = []
    for _, row in frame.iterrows():
        weight = row["weight_at_site"]
        cells = [
            f"{row['aep'] * 100:.3g}",
            cfs(row["flow_cfs"]),
            cfs(row["regression_flow_cfs"]),
            f"{weight:.2f}" if math.isfinite(weight) else "-",
            cfs(row["weighted_flow_cfs"]),
            cfs(row["weighted_lower_cfs"]),
            cfs(row["weighted_upper_cfs"]),
        ]
        if noted:
            cells.append(str(row["note"]))
        rows.append(tuple(cells))
    if noted:
        header.append("Note")
    blocks: List[Block] = [
        Paragraph(
            "The at-site flood quantiles were weighted with the regional regression estimates "
            "by the inverse of their variances, in log10 of flow, following Bulletin 17C "
            "Appendix 9 (England and others, 2019). The at-site variance is the variance of "
            "the fitted quantile and the regression's is its variance of prediction, so the "
            "estimate with the smaller variance carries the larger weight. The weighted limits "
            "treat the two estimates as independent. The computation was made with flowfreq "
            f"{_version()}."
        ),
        Table(
            tuple(header),
            tuple(rows),
            (*("right",) * 7, *(("left",) if noted else ())),
            "At-site and regression flood quantiles and their weighted estimate.",
        ),
    ]
    if bool(options.get("detail", True)):
        blocks.append(Heading(2, "At-site analysis"))
        blocks.extend(_demote(gage.blocks))
        blocks.append(Heading(2, "Regression estimate"))
        blocks.extend(_demote(regression.blocks))
    references = list(gage.references)
    references.extend(r for r in regression.references if r not in references)
    options["gage"], options["regression"] = gage.options, regression.options
    title = gage.title.replace("Flood frequency", "Weighted flood frequency", 1)
    return SectionContent(title=title, blocks=blocks, options=options, references=references)
