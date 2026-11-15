"""Report section kinds: each computes its tables and figures from flowfreq, and nothing else.

A section builder takes the options a :class:`~flowfreq.reporting.spec.SectionSpec` carries
and returns the section's blocks, its references and its options *as resolved*: whatever the
builder fetched (a gage's annual peaks, a regional skew looked up by state) is written back
as plain values, so the same options build the same section again without the network.
"""

from __future__ import annotations

import datetime as _dt
import io
import math
import re
from dataclasses import dataclass, replace
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .blocks import Block, Bullets, Figure, Heading, Paragraph, Table, parse_markdown, without_title
from .spec import ReportError

#: The reference every Bulletin 17C section cites.
BULLETIN_17C = (
    "England, J.F., Jr., Cohn, T.A., Faber, B.A., Stedinger, J.R., Thomas, W.O., Jr., "
    "Veilleux, A.G., Kiang, J.E., and Mason, R.R., Jr., 2019, Guidelines for determining "
    "flood flow frequency, Bulletin 17C (ver. 1.1, May 2019): U.S. Geological Survey "
    "Techniques and Methods, book 4, chap. B5, 148 p., https://doi.org/10.3133/tm4B5."
)

_PEAK_KEYS = (
    "water_year",
    "peak_flow_cfs",
    "qualification_code",
    "peak_date",
    "source",
    "retrieved",
)
_FIGURE_DPI = 200
#: The authors and year a USGS-style citation opens with: "Eash, D.A., and Veilleux, A.G., 2013,".
_AUTHORS_YEAR = re.compile(r"^(?P<authors>.+?),\s*(?P<year>(?:18|19|20)\d{2})[a-z]?,")
_INITIALS = re.compile(r"(?:[A-Z]\.[ -]?)+|Jr\.?|Sr\.?|II|III")


@dataclass(frozen=True)
class SectionContent:
    """What a section builder returns.

    Attributes
    ----------
    title : str
        The section's own heading, used when the spec gives none.
    blocks : list of Block
        Tables, figures and boilerplate text, all computed.
    options : dict
        The options as resolved: enough to build the same section offline.
    references : list of str
        Full citations for the report's reference list.
    """

    title: str
    blocks: List[Block]
    options: Dict[str, Any]
    references: List[str]


@dataclass(frozen=True)
class SectionKind:
    """A registered section kind: its name, what it takes and the function that builds it."""

    name: str
    summary: str
    required: Tuple[str, ...]
    optional: Tuple[str, ...]
    build: Callable[[Dict[str, Any]], SectionContent]

    def check(self, options: Mapping[str, Any]) -> None:
        """Raise :class:`ReportError` for an unknown or a missing option."""
        allowed = (*self.required, *self.optional)
        extra = sorted(set(options) - set(allowed))
        if extra:
            raise ReportError(
                f"{self.name}: unknown option(s) {', '.join(extra)}; "
                f"allowed: {', '.join(allowed) or 'none'}"
            )
        missing = [key for key in self.required if options.get(key) is None]
        if missing:
            raise ReportError(f"{self.name}: missing option(s) {', '.join(missing)}")


def _version() -> str:
    from .. import __version__

    return str(__version__)


def cfs(value: Any) -> str:
    """A flow for a table cell: whole cfs at 100 and above, three significant figures below."""
    if value is None:
        return "-"
    number = float(value)
    if not math.isfinite(number):
        return "-"
    return f"{number:,.0f}" if abs(number) >= 100 else f"{number:.3g}"


def short_citation(citation: str) -> str:
    """``Eash and others (2013)`` for a citation that opens with authors and a year.

    One author gives ``Eash (2013)`` and two ``Southard and Veilleux (2014)``. A citation in
    any other form is returned whole, so nothing is guessed.
    """
    match = _AUTHORS_YEAR.match(citation.strip())
    if match is None:
        return citation.strip()
    # "Eash, D.A., Barnes, K.K., and Veilleux, A.G.": the surnames are the parts that are
    # not initials.
    parts = [p.strip().removeprefix("and ").strip() for p in match["authors"].split(",")]
    surnames = [p for p in parts if p and not _INITIALS.fullmatch(p)]
    if not surnames:
        return citation.strip()
    if len(surnames) == 1:
        who = surnames[0]
    elif len(surnames) == 2:
        who = f"{surnames[0]} and {surnames[1]}"
    else:
        who = f"{surnames[0]} and others"
    return f"{who} ({match['year']})"


def _figure_png(fig: Any) -> bytes:
    """A pyplot figure as PNG bytes; the figure is closed."""
    import matplotlib.pyplot as plt

    buffer = io.BytesIO()
    # No Software entry: it carries the matplotlib version, which is not report content.
    fig.savefig(
        buffer, format="png", dpi=_FIGURE_DPI, bbox_inches="tight", metadata={"Software": None}
    )
    plt.close(fig)
    return buffer.getvalue()


def _demote(blocks: Sequence[Block], by: int = 1) -> List[Block]:
    """Push every heading ``by`` levels down, to nest a parsed page under a subsection."""
    return [replace(b, level=b.level + by) if isinstance(b, Heading) else b for b in blocks]


# -- text ---------------------------------------------------------------------------------


def _build_text(options: Dict[str, Any]) -> SectionContent:
    return SectionContent(title="", blocks=[], options=dict(options), references=[])


# -- gaged flood frequency ----------------------------------------------------------------


def download_peaks(site_no: str) -> Dict[str, Any]:
    """A gage's annual peaks and site attributes, in the form a spec stores.

    Parameters
    ----------
    site_no : str
        USGS site number.

    Returns
    -------
    dict
        ``peaks`` (the inline record: ``water_year``, ``peak_flow_cfs``,
        ``qualification_code``, ``peak_date``, ``source``, ``retrieved``), ``site_name`` and
        ``drainage_area_sqmi`` (either may be None).

    Raises
    ------
    ReportError
        The download failed or returned no peaks.
    """
    import requests

    from ..usgs import USGSgage

    gage = USGSgage(site_no)
    try:
        frame = gage.download_peak_flow()
    except (requests.RequestException, ValueError, KeyError) as exc:
        raise ReportError(
            f"annual peaks for USGS {gage.site_no} could not be downloaded: {exc}"
        ) from exc
    frame = frame.dropna(subset=["water_year", "peak_flow_cfs"]).sort_values("water_year")
    if frame.empty:
        raise ReportError(f"USGS {gage.site_no} has no annual peaks")

    def text(value: Any) -> Optional[str]:
        return None if pd.isna(value) or str(value).strip() == "" else str(value).strip()

    return {
        "peaks": {
            "water_year": [int(y) for y in frame["water_year"]],
            "peak_flow_cfs": [float(q) for q in frame["peak_flow_cfs"]],
            "qualification_code": [text(c) for c in frame["qualification_code"]],
            "peak_date": [None if pd.isna(d) else str(d)[:10] for d in frame["peak_date"]],
            "source": f"USGS Water Data API annual peaks, site {gage.site_no}",
            "retrieved": _dt.date.today().isoformat(),
        },
        "site_name": gage.site_name,
        "drainage_area_sqmi": gage.drainage_area,
    }


def _peak_frame(peaks: Any) -> pd.DataFrame:
    """The inline peaks as a frame, validated."""
    if not isinstance(peaks, dict):
        raise ReportError("gaged_ffa: peaks must be an object of columns")
    extra = sorted(set(peaks) - set(_PEAK_KEYS))
    if extra:
        raise ReportError(
            f"gaged_ffa: unknown peaks key(s) {', '.join(extra)}; allowed: {', '.join(_PEAK_KEYS)}"
        )
    years, flows = peaks.get("water_year"), peaks.get("peak_flow_cfs")
    if not isinstance(years, list) or not isinstance(flows, list) or not years:
        raise ReportError("gaged_ffa: peaks needs water_year and peak_flow_cfs lists")
    n = len(years)
    columns: Dict[str, List[Any]] = {"water_year": years, "peak_flow_cfs": flows}
    for key in ("qualification_code", "peak_date"):
        values = peaks.get(key)
        columns[key] = [None] * n if values is None else values
    if any(not isinstance(v, list) or len(v) != n for v in columns.values()):
        raise ReportError(f"gaged_ffa: every peaks column must be a list of {n} values")
    try:
        frame = pd.DataFrame(
            {
                "water_year": [int(y) for y in years],
                "peak_flow_cfs": [float(q) for q in flows],
                "qualification_code": columns["qualification_code"],
                "peak_date": columns["peak_date"],
            }
        )
    except (TypeError, ValueError) as exc:
        raise ReportError(f"gaged_ffa: peaks are not numeric: {exc}") from exc
    if not np.all(np.isfinite(frame["peak_flow_cfs"].to_numpy())):
        raise ReportError("gaged_ffa: peak_flow_cfs holds a non-finite value")
    if frame["water_year"].duplicated().any():
        raise ReportError("gaged_ffa: a water year appears more than once in peaks")
    return frame


def _resolve_skew(options: Dict[str, Any]) -> None:
    """Turn a regional-skew lookup into the values it found, in place."""
    by_state = options.get("regional_skew_state")
    by_site = bool(options.get("regional_skew_from_site"))
    if not by_state and not by_site:
        if options.get("regional_skew_region") is not None:
            raise ReportError("gaged_ffa: regional_skew_region needs regional_skew_state")
        return
    if by_state and by_site:
        raise ReportError(
            "gaged_ffa: give regional_skew_state or regional_skew_from_site, not both"
        )
    if options.get("regional_skew") is not None:
        raise ReportError(
            "gaged_ffa: give a regional_skew value or a lookup (regional_skew_state, "
            "regional_skew_from_site), not both"
        )
    from ..regional_skew import RegionalSkewUnavailable, regional_skew_for, regional_skew_for_site

    provisional = bool(options.get("allow_provisional_skew", False))
    try:
        if by_state:
            found = regional_skew_for(
                str(by_state), options.get("regional_skew_region"), allow_provisional=provisional
            )
        else:
            if not options.get("site_no"):
                raise ReportError("gaged_ffa: regional_skew_from_site needs site_no")
            found = regional_skew_for_site(
                str(options["site_no"]),
                allow_provisional=provisional,
                drainage_area_sqmi=options.get("drainage_area_sqmi"),
            )
    except RegionalSkewUnavailable as exc:
        raise ReportError(f"gaged_ffa: no regional skew: {exc}") from exc
    for key in ("regional_skew_state", "regional_skew_region", "regional_skew_from_site"):
        options.pop(key, None)
    options["regional_skew"] = float(found.skew)
    # The standard error is what is stored and what run_ffa is given, on the first build and
    # on every rebuild, so both square the same float.
    options["regional_skew_se"] = math.sqrt(float(found.skew_mse))
    options["regional_skew_citation"] = found.citation
    options["regional_skew_note"] = (found.skew_region or found.state) + (
        f"; {found.table}" if found.table else ""
    )
    if found.advisory:
        options["regional_skew_advisory"] = found.advisory


def _skew_sentence(options: Mapping[str, Any], source: str) -> str:
    from ..workflow import B17C_DEFAULT_SKEW, B17C_DEFAULT_SKEW_SE

    if source == "station":
        return "The station skew was used without regional weighting."
    if source == "default":
        return (
            f"The station skew was weighted with flowfreq's fallback regional skew of "
            f"{B17C_DEFAULT_SKEW:.3f} (standard error {B17C_DEFAULT_SKEW_SE:.2f}), which has no "
            "located source; a published regional skew should replace it."
        )
    cited = options.get("regional_skew_citation")
    note = options.get("regional_skew_note")
    sentence = (
        f"The station skew was weighted with a regional skew of "
        f"{float(options['regional_skew']):.3f} (standard error "
        f"{float(options['regional_skew_se']):.3f})"
    )
    if cited:
        sentence += f" from {short_citation(str(cited))}" + (f": {note}." if note else ".")
    else:
        sentence += ", supplied by the analyst without a citation."
    if options.get("regional_skew_advisory"):
        sentence += f" {options['regional_skew_advisory']}"
    return sentence


def _threshold_rows(thresholds: Any) -> List[Dict[str, Any]]:
    if thresholds is None:
        return []
    if not isinstance(thresholds, list):
        raise ReportError("gaged_ffa: perception_thresholds must be a list")
    for t in thresholds:
        if not isinstance(t, dict) or not {"start_year", "end_year", "threshold_cfs"} <= set(t):
            raise ReportError(
                "gaged_ffa: each perception threshold needs start_year, end_year, threshold_cfs"
            )
    return list(thresholds)


def _build_gaged_ffa(given: Dict[str, Any]) -> SectionContent:
    return _gaged_ffa(given)[0]


def _gaged_ffa(given: Dict[str, Any]) -> Tuple[SectionContent, Any]:
    """The gaged section and the fitted ``Bulletin17C`` behind it."""
    from ..freq_plot import plot_frequency_curve, plot_peak_flows_with_thresholds
    from ..workflow import run_ffa

    options = dict(given)
    site_no = str(options["site_no"]).strip() if options.get("site_no") else None
    if options.get("peaks") is None:
        if site_no is None:
            raise ReportError("gaged_ffa: give site_no (to download the peaks) or peaks")
        fetched = download_peaks(site_no)
        options["peaks"] = fetched["peaks"]
        for key in ("site_name", "drainage_area_sqmi"):
            if options.get(key) is None and fetched.get(key) is not None:
                options[key] = fetched[key]
    frame = _peak_frame(options["peaks"])
    _resolve_skew(options)
    thresholds = _threshold_rows(options.get("perception_thresholds"))
    apply_codes = bool(options.get("apply_peak_codes", True))
    codes = frame["qualification_code"].tolist()
    has_codes = any(c is not None for c in codes)

    result = run_ffa(
        frame["peak_flow_cfs"].to_numpy(dtype=float),
        frame["water_year"].to_numpy(dtype=int),
        regional_skew=options.get("regional_skew"),
        regional_skew_se=options.get("regional_skew_se"),
        perception_thresholds=thresholds or None,
        low_outlier_threshold_override=options.get("low_outlier_threshold"),
        use_default_skew=bool(options.get("use_default_skew", False)),
        station_skew_only=bool(options.get("station_skew_only", False)),
        peak_codes=codes if has_codes else None,
        apply_peak_codes=apply_codes,
        site_no=site_no,
        allow_regulated=bool(options.get("allow_regulated", False)),
        stationarity=bool(options.get("stationarity", True)),
        stationarity_figure=bool(options.get("stationarity", True)),
    )
    if result["error"]:
        raise ReportError(f"gaged_ffa: the frequency analysis failed: {result['error']}")
    p = result["parameters"]
    b17c = result["b17c"]
    site_name = options.get("site_name")
    label = " ".join(x for x in (f"USGS {site_no}" if site_no else None, site_name) if x)
    title = f"Flood frequency at {label}" if label else "Flood frequency"
    blocks: List[Block] = []

    # Station and record
    year_values = frame["water_year"].to_numpy(dtype=int)
    flow_values = frame["peak_flow_cfs"].to_numpy(dtype=float)
    largest, smallest = int(np.argmax(flow_values)), int(np.argmin(flow_values))
    record: List[Tuple[str, str]] = []
    if site_no:
        record.append(("USGS station number", site_no))
    if site_name:
        record.append(("Station name", str(site_name)))
    if options.get("drainage_area_sqmi") is not None:
        record.append(
            ("Drainage area (square miles)", f"{float(options['drainage_area_sqmi']):,.4g}")
        )
    record.extend(
        [
            ("Annual peaks", f"{len(frame)}"),
            ("Water years", f"{year_values.min()} to {year_values.max()}"),
            (
                "Largest peak (cfs)",
                f"{cfs(flow_values[largest])} (water year {year_values[largest]})",
            ),
            (
                "Smallest peak (cfs)",
                f"{cfs(flow_values[smallest])} (water year {year_values[smallest]})",
            ),
        ]
    )
    if options["peaks"].get("source"):
        retrieved = options["peaks"].get("retrieved")
        record.append(
            (
                "Data source",
                str(options["peaks"]["source"]) + (f", retrieved {retrieved}" if retrieved else ""),
            )
        )
    classification = p.get("site_classification")
    if classification:
        basis = "; ".join(str(b) for b in classification.get("basis") or [])
        record.append(
            (
                "Regulation screen",
                str(classification["site_class"]) + (f": {basis}" if basis else ""),
            )
        )
    blocks.append(
        Table(("Item", "Value"), tuple(record), ("left", "left"), "Station and peak-flow record.")
    )
    if classification and classification.get("regulation_override"):
        blocks.append(
            Paragraph(
                "**This record is classified as regulated.** Bulletin 17C describes "
                "unregulated flood peaks; the fit below was made at the analyst's direction "
                "and describes the regulated record as observed."
            )
        )

    # Method
    method = [
        "The annual peaks were fitted with a log-Pearson Type III distribution following "
        "Bulletin 17C (England and others, 2019): the Expected Moments Algorithm (EMA) for "
        "the moments and the Multiple Grubbs-Beck Test (MGBT) for potentially influential "
        f"low floods. The computation was made with flowfreq {_version()}.",
        _skew_sentence(options, str(p["regional_skew_source"])),
    ]
    threshold = float(p["low_outlier_threshold"] or 0.0)
    censored = int(p["n_low_outliers"] or 0)
    if p["low_outlier_source"] == "override":
        method.append(
            f"A low-outlier threshold of {cfs(threshold)} cfs was set by the analyst; "
            f"{censored} value(s) below it were treated as censored."
        )
    elif threshold > 0 and censored:
        method.append(
            f"MGBT set the low-outlier threshold at {cfs(threshold)} cfs; {censored} value(s) "
            "below it were treated as censored."
        )
    else:
        method.append("MGBT identified no potentially influential low floods.")
    if thresholds:
        method.append(
            "Perception thresholds: "
            + "; ".join(
                f"{int(t['start_year'])} to {int(t['end_year'])}, {cfs(t['threshold_cfs'])} cfs"
                for t in thresholds
            )
            + ". A year in those periods without a recorded peak was treated as below its "
            "threshold."
        )
    applied = p.get("peak_codes_applied") or {}
    if applied:
        method.append(
            "NWIS qualification codes were applied as peakfq applies them ("
            + ", ".join(f"code {code} on {count} peak(s)" for code, count in applied.items())
            + ")."
        )
    elif has_codes and not apply_codes:
        method.append(
            "NWIS qualification codes were not applied: every peak was fitted as an exact "
            "systematic value."
        )
    if result["method"] == "mom":
        method.append(
            "**EMA did not converge on this record, so the moments are from the method of "
            "moments.**"
        )
    elif not result["converged"]:
        method.append(
            "**EMA did not converge on this record; the results are its last iterate and "
            "should be checked before use.**"
        )
    blocks.append(Paragraph(" ".join(method)))

    # Fitted parameters
    fitted: List[Tuple[str, str]] = [
        ("Mean of log10 flow", f"{p['mean_log']:.4f}"),
        ("Standard deviation of log10 flow", f"{p['std_log']:.4f}"),
        ("Station skew", f"{p['skew_station']:.3f}"),
    ]
    if p["regional_skew"] is not None:
        fitted.append(("Regional skew", f"{p['regional_skew']:.3f}"))
        fitted.append(("Weighted skew", f"{p['skew_weighted']:.3f}"))
    fitted.extend(
        [
            ("Skew used", f"{p['skew_used']:.3f}"),
            ("Low-outlier threshold (cfs)", cfs(threshold) if threshold > 0 else "none"),
            ("Values censored below the threshold", f"{censored}"),
            ("Peaks flagged by MGBT", f"{int(p['n_mgbt_outliers'] or 0)}"),
        ]
    )
    blocks.append(
        Table(
            ("Parameter", "Value"),
            tuple(fitted),
            ("left", "right"),
            "Fitted log-Pearson Type III parameters and low-outlier screening.",
        )
    )

    # Quantiles
    quantiles = result["quantile_df"]
    blocks.append(
        Table(
            (
                "Annual exceedance probability (percent)",
                "Return interval (years)",
                "Flow (cfs)",
                "Lower 90-percent limit (cfs)",
                "Upper 90-percent limit (cfs)",
            ),
            tuple(
                (
                    f"{row['AEP (%)'] * 100:.3g}",
                    f"{row['Return Interval (yr)']:g}",
                    cfs(row["Flow (cfs)"]),
                    cfs(row["Lower 90% CI"]),
                    cfs(row["Upper 90% CI"]),
                )
                for _, row in quantiles.iterrows()
            ),
            ("right",) * 5,
            "Flood quantiles with 90-percent confidence limits.",
        )
    )

    # Figures
    blocks.append(
        Figure(
            _figure_png(plot_frequency_curve(b17c, site_name=site_name, site_no=site_no)),
            "Fitted frequency curve with 90-percent confidence limits and the annual peaks.",
        )
    )
    blocks.append(
        Figure(
            _figure_png(
                plot_peak_flows_with_thresholds(
                    frame,
                    site_name=site_name,
                    site_no=site_no,
                    thresholds=[
                        {
                            "start_year": int(t["start_year"]),
                            "end_year": int(t["end_year"]),
                            "lower_cfs": float(t["threshold_cfs"]),
                        }
                        for t in thresholds
                    ]
                    or None,
                    mgbt_threshold=threshold if threshold > 0 else None,
                    mgbt_threshold_source=str(p["low_outlier_source"]),
                )
            ),
            "Annual peak flows.",
        )
    )

    # Stationarity
    page = result.get("stationarity_report")
    screen = result.get("stationarity")
    if page is not None:
        _, body = without_title(parse_markdown(page.markdown))
        blocks.append(Heading(2, "Stationarity screen"))
        blocks.extend(_demote(body))
        blocks.append(
            Figure(
                _figure_png(page.figure),
                "Stationarity screen of the annual peaks: trend, change point and segment "
                "medians.",
            )
        )
    elif isinstance(screen, dict) and screen.get("error"):
        blocks.append(Heading(2, "Stationarity screen"))
        blocks.append(Paragraph(f"The stationarity screen was not run: {screen['error']}"))

    # The record itself
    if bool(options.get("peaks_table", True)):
        # Date and codes appear only when the record carries them.
        columns: List[Tuple[str, str, List[str]]] = [
            ("Water year", "right", [f"{int(y)}" for y in frame["water_year"]]),
        ]
        dates = ["" if d is None else str(d) for d in frame["peak_date"]]
        if any(dates):
            columns.append(("Date", "left", dates))
        columns.append(("Peak flow (cfs)", "right", [cfs(q) for q in frame["peak_flow_cfs"]]))
        if has_codes:
            columns.append(
                ("Qualification codes", "left", ["" if c is None else str(c) for c in codes])
            )
        blocks.append(Heading(2, "Annual peaks"))
        blocks.append(
            Table(
                tuple(name for name, _, _ in columns),
                tuple(zip(*(cells for _, _, cells in columns))),
                tuple(side for _, side, _ in columns),
                "Annual peak flows used in the analysis.",
            )
        )

    references = [BULLETIN_17C]
    if p["regional_skew_source"] == "user" and options.get("regional_skew_citation"):
        references.append(str(options["regional_skew_citation"]))
    content = SectionContent(title=title, blocks=blocks, options=options, references=references)
    return content, b17c


# -- ungaged: regression from characteristics ---------------------------------------------


_STREAMSTATS_KEYS = (
    "retrieved",
    "snapped_lat",
    "snapped_lon",
    "snap_distance_m",
    "polygon_area_sq_mi",
    "regions_located",
    "region_choice",
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
        opening += (
            f". The computation was made with flowfreq {_version()} from the stored values, "
            "so it does not change if StreamStats later does."
        )
    blocks: List[Block] = [Paragraph(opening)]
    used = _used_codes(estimate)
    show_all = record is None or bool(options.get("all_characteristics", False))
    for block in body:
        if isinstance(block, Table) and block.header and block.header[0] == "AEP (%)":
            block = replace(block, caption="Regression flood quantiles with prediction limits.")
        elif isinstance(block, Table) and block.header == ("Code", "Value") and not show_all:
            # StreamStats returns every characteristic it computes for the state.
            block = replace(
                block,
                rows=tuple(row for row in block.rows if row[0] in used),
                caption="Basin characteristics the equations use, from StreamStats.",
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
    content = SectionContent(title=title, blocks=blocks, options=options, references=references)
    return content, estimate


def _build_ungaged(given: Dict[str, Any]) -> SectionContent:
    return _ungaged(given)[0]


# -- at-site fit weighted with the regression ---------------------------------------------


def _nested(kind: str, name: str, value: Any) -> Dict[str, Any]:
    """A nested section's options, checked against the kind they belong to."""
    if not isinstance(value, dict):
        raise ReportError(f"weighted: {name} must be an object of {kind} options")
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


#: The section kinds a spec may name.
SECTION_KINDS: Dict[str, SectionKind] = {
    kind.name: kind
    for kind in (
        SectionKind(
            name="text",
            summary="A titled section of narrative only.",
            required=(),
            optional=(),
            build=_build_text,
        ),
        SectionKind(
            name="gaged_ffa",
            summary=(
                "Bulletin 17C flood frequency at a gage: record, method, fitted parameters, "
                "quantiles with confidence limits, frequency curve, stationarity screen."
            ),
            required=(),
            optional=(
                "site_no",
                "site_name",
                "drainage_area_sqmi",
                "peaks",
                "regional_skew",
                "regional_skew_se",
                "regional_skew_citation",
                "regional_skew_note",
                "regional_skew_advisory",
                "regional_skew_state",
                "regional_skew_region",
                "regional_skew_from_site",
                "allow_provisional_skew",
                "station_skew_only",
                "use_default_skew",
                "perception_thresholds",
                "low_outlier_threshold",
                "apply_peak_codes",
                "allow_regulated",
                "stationarity",
                "peaks_table",
            ),
            build=_build_gaged_ffa,
        ),
        SectionKind(
            name="ungaged",
            summary=(
                "Regional regression estimate at an ungaged site: regions, inputs, quantiles "
                "with prediction limits, caveats. The basin characteristics and region are "
                "the analyst's, or come from StreamStats when none are given and are then "
                "stored in the specification."
            ),
            required=("state", "lat", "lon"),
            optional=(
                "characteristics",
                "region",
                "streamstats",
                "all_characteristics",
                "site_name",
                "aeps",
                "allow_nss_snapshot",
                "allow_provisional",
                "allow_extrapolation",
                "weighting_space",
            ),
            build=_build_ungaged,
        ),
        SectionKind(
            name="weighted",
            summary=(
                "A gage's Bulletin 17C quantiles weighted with the regional regression at the "
                "gage (Bulletin 17C Appendix 9): the weighted table, then both analyses."
            ),
            required=("gage", "regression"),
            optional=("detail",),
            build=_build_weighted,
        ),
    )
}


def build_section(kind: str, options: Mapping[str, Any]) -> SectionContent:
    """Build one section.

    Parameters
    ----------
    kind : str
        A key of :data:`SECTION_KINDS`.
    options : mapping
        The kind's options (see ``SECTION_KINDS[kind].required`` and ``.optional``).

    Returns
    -------
    SectionContent

    Raises
    ------
    ReportError
        An unknown kind or option, a missing option, or an analysis that failed. Refusals
        the analyses raise themselves (a regulated record, an input outside an equation's
        range) pass through unchanged.
    """
    if kind not in SECTION_KINDS:
        raise ReportError(
            f"unknown section kind {kind!r}; choose from {', '.join(sorted(SECTION_KINDS))}"
        )
    section = SECTION_KINDS[kind]
    section.check(options)
    return section.build(dict(options))


__all__ = [
    "BULLETIN_17C",
    "SECTION_KINDS",
    "SectionContent",
    "SectionKind",
    "build_section",
    "cfs",
    "download_peaks",
    "short_citation",
]
