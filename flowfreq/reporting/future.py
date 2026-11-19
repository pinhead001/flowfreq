"""The ``future_flow`` section: flood quantiles under a future condition."""

from __future__ import annotations

from typing import Any, Dict, List

from ..aep import match_aep, same_aep
from ._base import SectionContent, _demote, _version, cfs, short_citation
from .blocks import Block, Heading, Paragraph, Table
from .gaged import _gaged_ffa
from .regional import _ungaged
from .spec import ReportError
from .weighted import _nested

# -- flood quantiles under a future condition ---------------------------------------------

_FACTOR_SET_KEYS = (
    "name",
    "aeps",
    "factors",
    "scenario",
    "horizon",
    "citation",
    "level",
    "geography",
    "legal_status",
)


#: The caption of a ``future_flow`` section's main table.
FUTURE_FLOW_CAPTION = "Flood quantiles under the future condition."


def _factor_set(given: Any) -> Any:
    """The analyst's change factors as a ``ChangeFactorSet``."""
    from ..future_flow import factor_set_from_dict

    if not isinstance(given, dict):
        raise ReportError("future_flow: factor_set must be an object")
    extra = sorted(set(given) - set(_FACTOR_SET_KEYS))
    if extra:
        raise ReportError(
            f"future_flow: unknown factor_set key(s) {', '.join(extra)}; "
            f"allowed: {', '.join(_FACTOR_SET_KEYS)}"
        )
    try:
        return factor_set_from_dict(given)
    except KeyError as exc:
        raise ReportError(f"future_flow: {exc.args[0]}") from exc
    except (TypeError, ValueError) as exc:
        raise ReportError(f"future_flow: factor_set: {exc}") from exc


def _build_future_flow(given: Dict[str, Any]) -> SectionContent:
    from ..future_flow import apply_change_factors

    options = dict(given)
    has_gage = options.get("gage") is not None
    if has_gage == (options.get("regression") is not None):
        raise ReportError(
            "future_flow: give the present-condition analysis as gage (gaged_ffa options) or "
            "regression (ungaged options), one of the two"
        )
    factors = _factor_set(options["factor_set"])
    if has_gage:
        base, b17c = _gaged_ffa(_nested("gaged_ffa", "gage", options["gage"], "future_flow"))
        present = b17c.results.quantiles[["aep", "flow_cfs"]]
        source = "the Bulletin 17C fit at the gage"
        options["gage"] = base.options
    else:
        base, estimate = _ungaged(
            _nested("ungaged", "regression", options["regression"], "future_flow")
        )
        present = estimate.quantiles[["aep", "flow_cfs"]]
        source = "the regional regression estimate"
        options["regression"] = base.options
    present = present[present["flow_cfs"].notna()].sort_values("aep", ascending=False)
    have = [float(a) for a in present["aep"]]
    published = [float(a) for a in factors.aeps]
    asked = options.get("aeps")
    if asked is None:
        # Only the AEPs the factors cover: nothing is held beyond them unless asked for.
        low, high = min(published), max(published)
        keep = [a for a in have if low <= a <= high or same_aep(a, low) or same_aep(a, high)]
        if not keep:
            raise ReportError(
                "future_flow: the present-condition quantiles have no AEP within the factor "
                f"set's range ({min(published):g} to {max(published):g}); name the AEPs to "
                "scale with aeps, and the end factor is held for them"
            )
    else:
        if not isinstance(asked, list) or not asked:
            raise ReportError("future_flow: aeps must be a non-empty list of AEPs")
        keep = []
        for value in asked:
            try:
                wanted = float(value)
            except (TypeError, ValueError) as exc:
                raise ReportError("future_flow: aeps must be a list of numbers") from exc
            match = [a for a in have if abs(a - wanted) < 1e-9]
            if not match:
                raise ReportError(
                    f"future_flow: the present-condition analysis has no quantile at AEP "
                    f"{wanted:g}; it has " + ", ".join(f"{a:g}" for a in have)
                )
            keep.append(match[0])
    present = present[present["aep"].isin(keep)].reset_index(drop=True)
    future = apply_change_factors(present, factors)

    def basis(aep: float, extrapolated: bool) -> str:
        if extrapolated:
            return "held beyond the published AEPs"
        if match_aep(published, aep) is not None:
            return "published"
        return "interpolated"

    cited = short_citation(factors.citation)
    blocks: List[Block] = [
        Paragraph(
            f"The flood quantiles from {source} were multiplied by change factors for "
            f"{factors.scenario}, {factors.horizon}, from {cited} (factor set "
            f"{factors.name!r}; {factors.level}-level, {factors.legal_status}). A factor "
            "between published AEPs is interpolated on the standard normal variate of the "
            "AEP. The factors and their citation were supplied by the analyst: flowfreq "
            "ships no change factors and did not check these against their source. A "
            "factor scales the flow only, so the present-condition confidence or "
            "prediction limits do not carry over and none are given. The computation was "
            f"made with flowfreq {_version()}."
        ),
        Table(
            ("Annual exceedance probability (percent)", "Change factor"),
            tuple(
                (f"{aep * 100:.3g}", f"{factor:.3f}")
                for aep, factor in sorted(zip(published, factors.factors), reverse=True)
            ),
            ("right", "right"),
            "Change factors as supplied.",
        ),
        Table(
            (
                "Annual exceedance probability (percent)",
                "Present flow (cfs)",
                "Change factor",
                "Future flow (cfs)",
                "Factor",
            ),
            tuple(
                (
                    f"{row['aep'] * 100:.3g}",
                    cfs(row["regulatory_flow_cfs"]),
                    f"{row['change_factor']:.3f}",
                    cfs(row["future_flow_cfs"]),
                    basis(float(row["aep"]), bool(row["extrapolated"])),
                )
                for _, row in future.table.iterrows()
            ),
            ("right", "right", "right", "right", "left"),
            FUTURE_FLOW_CAPTION,
        ),
    ]
    if bool(future.table["extrapolated"].any()):
        blocks.append(
            Paragraph(
                "**A factor held beyond the published AEPs is the nearest published factor, "
                "not a published value for that AEP.**"
            )
        )
    if bool(options.get("detail", True)):
        blocks.append(Heading(2, "Present-condition analysis"))
        blocks.extend(_demote(base.blocks))
    references = list(base.references)
    if factors.citation not in references:
        references.append(factors.citation)
    title = base.title[:1].lower() + base.title[1:] if base.title else "flood quantiles"
    return SectionContent(
        title=f"Future-condition {title}", blocks=blocks, options=options, references=references
    )
