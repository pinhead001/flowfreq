"""The ``ungaged`` section: the regional regression at a site."""

from __future__ import annotations

import datetime as _dt
from dataclasses import replace
from typing import Any, Dict, List, Mapping, Optional, Tuple

import numpy as np
import pandas as pd

from ._base import SectionContent, _version
from .blocks import Block, Bullets, Paragraph, Table, parse_markdown, without_title
from .spec import ReportError

# -- ungaged: regression from characteristics ---------------------------------------------


_STREAMSTATS_KEYS = (
    "retrieved",
    "snapped_lat",
    "snapped_lon",
    "snap_distance_m",
    "polygon_area_sq_mi",
    "regions_located",
    "region_choice",
    "supplied",
    "notes",
)


def _streamstats_notes(estimate: Any) -> List[str]:
    """What the StreamStats run set aside or warned of, which the stored inputs alone lose."""
    provenance = estimate.provenance
    notes: List[str] = []
    if provenance.get("region_choice_source"):
        notes.append(f"Region choice: {provenance['region_choice_source']}.")
    alternatives = provenance.get("alternative_regions") or []
    if alternatives:
        notes.append(
            "Alternative NSS equation sets for the same area, not used: "
            + "; ".join(
                f"{a['code']} ({a['name']}): {a['citation'] or 'citation unknown'}"
                for a in alternatives
            )
            + "."
        )
    if estimate.unmatched_regions:
        notes.append(
            "Located by NSS but not in the library, so excluded and the area fractions "
            "renormalised: "
            + ", ".join(
                f"{r.code} ({r.name}, {r.percent_weight:g} percent)"
                for r in estimate.unmatched_regions
            )
            + "."
        )
    missing = provenance.get("missing_characteristics") or {}
    if missing:
        notes.append(
            "Regions skipped because StreamStats returned no value for a characteristic "
            "their equations need: "
            + "; ".join(f"{code} needs {', '.join(codes)}" for code, codes in missing.items())
            + "."
        )
    zero_weight = provenance.get("zero_weight_regions") or []
    if zero_weight:
        notes.append(
            "Regions NSS located at 0 percent of the basin (an overlap under half a percent), "
            "dropped: " + ", ".join(zero_weight) + "."
        )
    advisories = (provenance.get("delineation") or {}).get("advisories") or []
    if advisories:
        notes.append(
            "StreamStats warned about this basin: "
            + "; ".join(str(a) for a in advisories)
            + ". The regression equations may not represent it; check against local data."
        )
    withheld = sorted(
        code
        for code, why in ((provenance.get("delineation") or {}).get("unavailable") or {}).items()
        if str(why).startswith("withheld")
    )
    if withheld:
        notes.append(
            "Characteristics StreamStats computed and said may be inaccurate for this basin, "
            "withheld from the equations: " + ", ".join(withheld) + "."
        )
    covered = provenance.get("covered_percent")
    if covered is not None and float(covered) < 98.0:
        notes.append(
            f"The regression regions used cover {float(covered):.0f} percent of the basin; "
            "their shares were scaled up to the whole, so the rest (in another state or "
            "country, or in a region with no equations here) is treated as like the part "
            "covered."
        )
    return notes


def _same_quantiles(a: pd.DataFrame, b: pd.DataFrame) -> bool:
    if len(a) != len(b):
        return False
    return all(
        np.allclose(
            a[column].to_numpy(dtype=float),
            b[column].to_numpy(dtype=float),
            rtol=1e-9,
            atol=0.0,
            equal_nan=True,
        )
        for column in ("aep", "flow_cfs", "lower", "upper")
    )


def _freeze_streamstats(options: Dict[str, Any], lat: float, lon: float) -> None:
    """Delineate the site with StreamStats and write what it returned into the options.

    Afterwards the options hold the basin characteristics, the regions with their percent
    of the basin and a ``streamstats`` record of the run, so the section rebuilds from
    them without StreamStats. The estimate the stored inputs give is checked against the
    one StreamStats gave, since the report shows the former.
    """
    from ..ungaged import estimate_at_site, estimate_from_characteristics

    chosen = options.get("region")
    if chosen is not None and not isinstance(chosen, (str, list)):
        raise ReportError(
            "ungaged: without characteristics, region is the choice among overlapping NSS "
            "regions: a code or a list of codes"
        )
    snapshot = bool(options.get("allow_nss_snapshot", False))
    provisional = bool(options.get("allow_provisional", False))
    extrapolate = bool(options.get("allow_extrapolation", False))
    coverage: Dict[str, Any] = {}
    if options.get("min_coverage_percent") is not None:
        coverage["min_coverage_percent"] = options["min_coverage_percent"]
    supplied = options.pop("supplied_characteristics", None)
    extra: Dict[str, Any] = {}
    if supplied is not None:
        if not isinstance(supplied, dict) or not supplied:
            raise ReportError(
                "ungaged: supplied_characteristics must be an object of code to value"
            )
        extra["supplied_characteristics"] = supplied
    live = estimate_at_site(
        lat,
        lon,
        state=str(options["state"]),
        aeps=options.get("aeps"),
        weighting_space=options.get("weighting_space"),
        region=chosen,
        allow_nss_snapshot=snapshot,
        allow_provisional=provisional,
        allow_extrapolation=extrapolate,
        allow_inaccurate=bool(options.get("allow_inaccurate", False)),
        **coverage,
        **extra,
    )
    total = sum(r.percent_weight for r in live.regions)
    if not live.regions or total <= 0:
        raise ReportError("ungaged: StreamStats located no regression region the library has")
    regions: Dict[str, float] = {}
    for r in live.regions:
        code = live.region_codes.get(r.code, r.code)
        regions[code] = regions.get(code, 0.0) + 100.0 * r.percent_weight / total
    characteristics = {str(code): float(value) for code, value in live.characteristics.items()}
    stored = estimate_from_characteristics(
        str(options["state"]),
        characteristics,
        region=regions,
        lat=lat,
        lon=lon,
        aeps=options.get("aeps"),
        weighting_space=options.get("weighting_space"),
        allow_nss_snapshot=snapshot,
        allow_provisional=provisional,
        allow_extrapolation=extrapolate,
    )
    if not _same_quantiles(live.quantiles, stored.quantiles):
        raise ReportError(
            "ungaged: the estimate from the characteristics and regions StreamStats returned "
            "does not reproduce the StreamStats estimate itself, so the section cannot be "
            "stored faithfully. Give characteristics and region yourself instead."
        )
    snap = live.provenance.get("snap") or {}
    delineation = live.provenance.get("delineation") or {}
    record: Dict[str, Any] = {
        "retrieved": _dt.date.today().isoformat(),
        "snapped_lat": snap.get("snapped_lat"),
        "snapped_lon": snap.get("snapped_lon"),
        "snap_distance_m": snap.get("distance_m"),
        "polygon_area_sq_mi": delineation.get("polygon_area_sq_mi"),
        "regions_located": [
            {"code": r.code, "name": r.name, "percent": round(float(r.percent_weight), 6)}
            for r in live.regions
        ],
        "notes": _streamstats_notes(live),
    }
    if chosen is not None:
        record["region_choice"] = chosen
    if supplied:
        record["supplied"] = sorted(supplied)
    options["characteristics"] = characteristics
    options["region"] = regions
    options["streamstats"] = record


def _streamstats_record(options: Mapping[str, Any]) -> Optional[Dict[str, Any]]:
    record = options.get("streamstats")
    if record is None:
        return None
    if not isinstance(record, dict):
        raise ReportError("ungaged: streamstats must be an object")
    extra = sorted(set(record) - set(_STREAMSTATS_KEYS))
    if extra:
        raise ReportError(
            f"ungaged: unknown streamstats key(s) {', '.join(extra)}; "
            f"allowed: {', '.join(_STREAMSTATS_KEYS)}"
        )
    notes = record.get("notes") or []
    if not isinstance(notes, list) or any(not isinstance(n, str) for n in notes):
        raise ReportError("ungaged: streamstats notes must be a list of sentences")
    supplied = record.get("supplied") or []
    if not isinstance(supplied, list) or any(not isinstance(c, str) for c in supplied):
        raise ReportError("ungaged: streamstats supplied must be a list of codes")
    return record


def _used_codes(estimate: Any) -> List[str]:
    """The characteristic codes the equations used read, in the order first met."""
    codes: List[str] = []
    for per_region in estimate.estimates.values():
        for _, each in per_region:
            for variable in (*each.equation.variables, *each.equation.applicability):
                if variable.code not in codes:
                    codes.append(variable.code)
    return codes


#: The quantile table's columns as an ungaged section prints them.
_UNGAGED_HEADERS = {
    "AEP (%)": "Annual exceedance probability (percent)",
    "Lower": "Lower prediction limit (cfs)",
    "Upper": "Upper prediction limit (cfs)",
    "Level": "Level of the limits",
    "Interval": "How the limits were computed",
    "Extrapolated": "Outside the calibrated range",
    "Regions": "Regression regions",
}


def _ungaged(given: Dict[str, Any]) -> Tuple[SectionContent, Any]:
    """The ungaged section and the estimate behind it."""
    from ..ungaged import estimate_from_characteristics

    options = dict(given)
    try:
        lat, lon = float(options["lat"]), float(options["lon"])
    except (TypeError, ValueError) as exc:
        raise ReportError(f"ungaged: lat and lon must be numbers: {exc}") from exc
    if options.get("characteristics") is None:
        if options.get("streamstats") is not None:
            raise ReportError("ungaged: a streamstats record needs the characteristics it gave")
        _freeze_streamstats(options, lat, lon)
    elif options.get("supplied_characteristics") is not None:
        raise ReportError(
            "ungaged: supplied_characteristics is for a section without characteristics, whose "
            "basin StreamStats gives; put the values in characteristics"
        )
    record = _streamstats_record(options)
    characteristics = options["characteristics"]
    if not isinstance(characteristics, dict) or not characteristics:
        raise ReportError("ungaged: characteristics must be an object of code to value")
    region = options.get("region")
    if not isinstance(region, (str, dict)):
        raise ReportError(
            "ungaged: with characteristics, region must be a region code or an object of "
            "code to percent"
        )
    try:
        values = {str(code): float(value) for code, value in characteristics.items()}
    except (TypeError, ValueError) as exc:
        raise ReportError(f"ungaged: characteristics must be numbers: {exc}") from exc
    estimate = estimate_from_characteristics(
        str(options["state"]),
        values,
        region=region,
        lat=lat,
        lon=lon,
        aeps=options.get("aeps"),
        allow_nss_snapshot=bool(options.get("allow_nss_snapshot", False)),
        allow_provisional=bool(options.get("allow_provisional", False)),
        allow_extrapolation=bool(options.get("allow_extrapolation", False)),
        weighting_space=options.get("weighting_space"),
    )
    analysts = list((record or {}).get("supplied") or [])
    by_analyst = " and ".join(filter(None, [", ".join(analysts[:-1]), *analysts[-1:]]))
    heading, body = without_title(parse_markdown(estimate.to_markdown()))
    site_name = options.get("site_name")
    if record is None:
        opening = (
            "Flood quantiles were estimated from published regional regression equations, "
            "evaluated on basin characteristics supplied by the analyst. The computation was "
            f"made with flowfreq {_version()}; the characteristics were not measured by it."
        )
    else:
        opening = (
            "Flood quantiles were estimated from published regional regression equations, "
            "evaluated on basin characteristics from USGS StreamStats"
            + (f", retrieved {record['retrieved']}" if record.get("retrieved") else "")
            + f". StreamStats delineated the basin above {lat:.5f}, {lon:.5f}"
        )
        if record.get("snapped_lat") is not None and record.get("snapped_lon") is not None:
            opening += (
                f", which it moved to {float(record['snapped_lat']):.5f}, "
                f"{float(record['snapped_lon']):.5f} on its stream network"
            )
            if record.get("snap_distance_m") is not None:
                opening += f" ({float(record['snap_distance_m']):,.0f} m away)"
        if record.get("polygon_area_sq_mi") is not None:
            opening += (
                f"; the delineated basin is {float(record['polygon_area_sq_mi']):,.4g} "
                "square miles"
            )
        if by_analyst:
            verb = "was" if len(analysts) == 1 else "were"
            opening += f". {by_analyst} {verb} supplied by the analyst, not StreamStats"
        opening += (
            f". The computation was made with flowfreq {_version()} from the stored values, "
            "so it does not change if StreamStats later does."
        )
    blocks: List[Block] = [Paragraph(opening)]
    used = _used_codes(estimate)
    show_all = record is None or bool(options.get("all_characteristics", False))
    for block in body:
        if isinstance(block, Table) and block.header and block.header[0] == "AEP (%)":
            # The estimate's own report abbreviates; a section spells its columns out
            # with their units, as the other kinds do.
            spelled = tuple(_UNGAGED_HEADERS.get(name, name) for name in block.header)
            block = replace(
                block,
                header=spelled,
                caption="Regression flood quantiles with prediction limits.",
            )
        elif isinstance(block, Table) and block.header and block.header[0] == "NSS code":
            block = replace(
                block, caption="Regression regions the basin is in, and the equations' source."
            )
        elif isinstance(block, Table) and block.header == ("Code", "Value") and show_all:
            block = replace(block, caption="Basin characteristics.")
        elif isinstance(block, Table) and block.header == ("Code", "Value") and not show_all:
            # StreamStats returns every characteristic it computes for the state.
            block = replace(
                block,
                rows=tuple(row for row in block.rows if row[0] in used),
                caption="Basin characteristics the equations use, from StreamStats"
                + (f" except {by_analyst}, supplied by the analyst." if by_analyst else "."),
            )
        elif isinstance(block, Bullets) and block is body[-1] and record is not None:
            # The caveats list is the last block; the stored notes belong in it.
            block = replace(block, items=(*block.items, *(record.get("notes") or ())))
        blocks.append(block)

    references: List[str] = []
    for per_region in estimate.estimates.values():
        for _, each in per_region:
            c = each.equation.citation
            cited = c.publication.rstrip(". ")
            if c.year and str(c.year) not in cited:
                cited += f" ({c.year})"
            cited += f". {c.url}" if c.url else "."
            if cited not in references:
                references.append(cited)
    title = (
        f"Regression estimate at {site_name}" if site_name else (heading or "Regression estimate")
    )
    content = SectionContent(
        title=title,
        blocks=blocks,
        options=options,
        references=references,
        cautions=[estimate.advisory] if estimate.advisory else [],
    )
    return content, estimate


def _caution_paragraphs(label: str, content: SectionContent) -> List[Block]:
    """An estimate's cautions, for a section that shows its flows without its own text."""
    return [Paragraph(f"Cautions on {label}: {text}") for text in content.cautions]


def _build_ungaged(given: Dict[str, Any]) -> SectionContent:
    return _ungaged(given)[0]
