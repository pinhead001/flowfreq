"""The ``transposed`` section: a gage's quantiles moved to a site by area ratio."""

from __future__ import annotations

import math
from typing import Any, Dict, List, Tuple

from ._base import SectionContent, _demote, _version, cfs, short_citation
from .blocks import Block, Bullets, Heading, Paragraph, Table, parse_markdown, without_title
from .gaged import _gaged_ffa
from .spec import ReportError
from .weighted import _nested

# -- a gage's quantiles transposed to an ungaged site -------------------------------------

_EXPONENT_KEYS = ("citation", "exponent", "aeps", "exponents", "region", "valid_area_range_sqmi")


def _positive(owner: str, name: str, value: Any) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        number = math.nan
    if isinstance(value, bool) or not math.isfinite(number) or number <= 0.0:
        raise ReportError(f"{owner}: {name} must be a positive number, got {value!r}")
    return number


def _regression_exponents(given: Any) -> Any:
    """The spec's ``exponents`` object as a :class:`~flowfreq.transpose.RegressionExponents`."""
    from ..transpose import RegressionExponents

    if not isinstance(given, dict):
        raise ReportError(
            "transposed: exponents must be an object with a citation and either exponent "
            "(one value for every AEP) or aeps and exponents (the published values)"
        )
    extra = sorted(set(given) - set(_EXPONENT_KEYS))
    if extra:
        raise ReportError(
            f"transposed: unknown exponents key(s) {', '.join(extra)}; "
            f"allowed: {', '.join(_EXPONENT_KEYS)}"
        )
    citation = given.get("citation")
    if not isinstance(citation, str) or not citation.strip():
        raise ReportError(
            "transposed: exponents needs a citation naming the report and table the "
            "exponents came from"
        )
    one = given.get("exponent")
    by_aep = given.get("aeps") is not None or given.get("exponents") is not None
    if (one is None) != by_aep:
        raise ReportError(
            "transposed: exponents takes either exponent (one value for every AEP) or aeps "
            "and exponents (the published values), not both and not neither"
        )
    area_range = given.get("valid_area_range_sqmi")
    if area_range is not None:
        if not isinstance(area_range, (list, tuple)) or len(area_range) != 2:
            raise ReportError(
                "transposed: exponents valid_area_range_sqmi must be the smallest and the "
                "largest drainage area, in square miles"
            )
        area_range = (
            _positive("transposed", "exponents valid_area_range_sqmi", area_range[0]),
            _positive("transposed", "exponents valid_area_range_sqmi", area_range[1]),
        )
    region = str(given.get("region") or "")
    if one is not None:
        return RegressionExponents.constant(
            _positive("transposed", "exponents exponent", one),
            citation,
            region=region,
            valid_area_range_sqmi=area_range,
        )
    if given.get("aeps") is None or given.get("exponents") is None:
        raise ReportError("transposed: exponents by AEP needs both aeps and exponents")
    return RegressionExponents(
        aeps=given["aeps"],
        exponents=given["exponents"],
        citation=citation,
        region=region,
        valid_area_range_sqmi=area_range,
    )


def _site_class(name: str, value: Any) -> Any:
    from ..regulation import SiteClass

    if value is None:
        return None
    try:
        return SiteClass(str(value))
    except ValueError:
        raise ReportError(
            f"transposed: {name} must be one of {', '.join(c.value for c in SiteClass)}; "
            f"got {value!r}"
        ) from None


def _build_transposed(given: Dict[str, Any]) -> SectionContent:
    return _transposed(given)[0]


def _transposed(given: Dict[str, Any]) -> Tuple[SectionContent, Any]:
    """The transposed section and the ``TransposedFrequency`` behind it."""
    from ..transpose import transpose_frequency

    options = dict(given)
    donor, b17c = _gaged_ffa(_nested("gaged_ffa", "donor", options["donor"], "transposed"))
    target_area = _positive("transposed", "target_area_sqmi", options["target_area_sqmi"])
    donor_area = options.get("donor_area_sqmi")
    if donor_area is None:
        donor_area = donor.options.get("drainage_area_sqmi")
    if donor_area is None:
        raise ReportError(
            "transposed: the donor has no drainage area; give donor_area_sqmi, or "
            "drainage_area_sqmi in the donor's options"
        )
    donor_area = _positive("transposed", "donor_area_sqmi", donor_area)
    exponents = _regression_exponents(options["exponents"])
    max_aep = options.get("max_aep")
    if max_aep is not None:
        max_aep = _positive("transposed", "max_aep", max_aep)
    site_no = donor.options.get("site_no")
    target_name = options.get("target_name")
    results = transpose_frequency(
        b17c.results,
        donor_area,
        target_area,
        exponents,
        allow_out_of_range=bool(options.get("allow_out_of_range", False)),
        max_aep=max_aep,
        donor_site_no=f"USGS {site_no}" if site_no else None,
        target_name=str(target_name) if target_name else None,
        donor_regulation=_site_class("donor_site_class", options.get("donor_site_class")),
        target_regulation=_site_class("target_site_class", options.get("target_site_class")),
        allow_regulation_mismatch=bool(options.get("allow_regulation_mismatch", False)),
    )

    limits = results.confidence_limits
    limit_at = {} if limits.empty else {float(r["aep"]): r for _, r in limits.iterrows()}
    header = [
        "Annual exceedance probability (percent)",
        "Donor flow (cfs)",
        "Drainage-area exponent",
        "Factor applied",
        "Transposed flow (cfs)",
    ]
    if limit_at:
        header += ["Lower 90-percent limit (cfs)", "Upper 90-percent limit (cfs)"]
    header.append("Exponent basis")
    rows = []
    for _, row in results.quantiles.iterrows():
        cells = [
            f"{row['aep'] * 100:.3g}",
            cfs(row["donor_flow_cfs"]),
            f"{row['exponent']:.3f}",
            f"{row['area_ratio_factor']:.4f}",
            cfs(row["flow_cfs"]),
        ]
        if limit_at:
            limit = limit_at.get(float(row["aep"]))
            cells += [
                cfs(limit["lower_5pct"]) if limit is not None else "-",
                cfs(limit["upper_5pct"]) if limit is not None else "-",
            ]
        cells.append(str(row["source"]))
        rows.append(tuple(cells))
    table = Table(
        tuple(header),
        tuple(rows),
        (*("right",) * (len(header) - 1), "left"),
        "Flood quantiles transposed from the gage by the drainage-area ratio.",
    )

    where = str(target_name) if target_name else "the ungaged site"
    gage_label = " ".join(
        x for x in (f"USGS {site_no}" if site_no else None, donor.options.get("site_name")) if x
    )
    blocks: List[Block] = [
        Paragraph(
            f"Flood quantiles at {where} were transposed from the donor gage"
            f"{', ' + gage_label + ',' if gage_label else ''} by the drainage-area ratio "
            "method: the gage's quantile at each annual exceedance probability, multiplied "
            "by the ratio of the site's drainage area to the gage's raised to the "
            "drainage-area exponent of a regional regression at that probability. The "
            f"exponents are from {short_citation(exponents.citation)}. No distribution was "
            "fitted at the site, so there are no moments to report there. The computation "
            f"was made with flowfreq {_version()}."
        )
    ]
    _, body = without_title(parse_markdown(results.to_markdown()))
    for block in body:
        if isinstance(block, Paragraph) and "Q_target" in block.text:
            continue  # the formula, which the paragraph above gives in words
        if isinstance(block, Heading) and block.text == "Quantiles":
            continue
        blocks.append(table if isinstance(block, Table) else block)
    # transpose_frequency only logs this one, and a report has no log.
    developed_over = exponents.valid_area_range_sqmi
    if developed_over is not None:
        low, high = developed_over
        outside = tuple(
            f"**The {label}'s drainage area, {area:,.4g} square miles, is outside the range "
            f"the cited regression was developed over** ({low:,.4g} to {high:,.4g} square "
            "miles), so the exponent is applied beyond its basis."
            for label, area in (("donor gage", donor_area), ("site", target_area))
            if not low <= area <= high
        )
        if outside and isinstance(blocks[-1], Bullets):
            blocks[-1] = Bullets((*outside, *blocks[-1].items))
        elif outside:
            blocks.append(Bullets(outside))
    if bool(options.get("detail", True)):
        blocks.append(Heading(2, "Analysis at the donor gage"))
        blocks.extend(_demote(donor.blocks))
    references = list(donor.references)
    if exponents.citation.strip() not in references:
        references.append(exponents.citation.strip())
    options["donor"] = donor.options
    options["donor_area_sqmi"] = donor_area
    title = f"Transposed flood frequency at {target_name}" if target_name else ""
    content = SectionContent(
        title=title or "Transposed flood frequency",
        blocks=blocks,
        options=options,
        references=references,
    )
    return content, results
