"""Report section kinds: each computes its tables and figures from flowfreq, and nothing else.

A section builder takes the options a :class:`~flowfreq.reporting.spec.SectionSpec` carries
and returns the section's blocks, its references and its options *as resolved*: whatever the
builder fetched (a gage's annual peaks, a regional skew looked up by state) is written back
as plain values, so the same options build the same section again without the network.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import io
import json
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


_EXTENSION_KEYS = ("index", "n_extend", "years")
_INDEX_KEYS = ("site_no", "site_name", "peaks")


def _extended_frame(options: Dict[str, Any], frame: pd.DataFrame) -> Tuple[Any, pd.DataFrame]:
    """The record extended from an index gage (Bulletin 17C Appendix 8), and the extension.

    Downloads the index gage's peaks when the options carry none, and writes them back.
    The frame returned has the observed peaks and the generated ones, with a ``source``
    column saying which.
    """
    from ..peak_extension import extend_peaks

    given = options["extension"]
    if not isinstance(given, dict) or not isinstance(given.get("index"), dict):
        raise ReportError("gaged_ffa: extension must be an object with an index gage")
    extra = sorted(set(given) - set(_EXTENSION_KEYS)) + sorted(
        set(given["index"]) - set(_INDEX_KEYS)
    )
    if extra:
        raise ReportError(
            f"gaged_ffa: unknown extension key(s) {', '.join(extra)}; an extension takes "
            f"{', '.join(_EXTENSION_KEYS)} and its index {', '.join(_INDEX_KEYS)}"
        )
    extension = dict(given)
    index = dict(extension["index"])
    if index.get("peaks") is None:
        if not index.get("site_no"):
            raise ReportError("gaged_ffa: the extension's index gage needs site_no or peaks")
        fetched = download_peaks(str(index["site_no"]).strip())
        index["peaks"] = fetched["peaks"]
        if index.get("site_name") is None and fetched.get("site_name") is not None:
            index["site_name"] = fetched["site_name"]
    index_frame = _peak_frame(index["peaks"])
    try:
        result = extend_peaks(
            pd.Series(frame["peak_flow_cfs"].to_numpy(), index=frame["water_year"].to_numpy()),
            pd.Series(
                index_frame["peak_flow_cfs"].to_numpy(), index=index_frame["water_year"].to_numpy()
            ),
            n_extend=extension.get("n_extend"),
            years=extension.get("years"),
        )
    except (TypeError, ValueError) as exc:
        raise ReportError(f"gaged_ffa: the record could not be extended: {exc}") from exc
    extension["index"] = index
    options["extension"] = extension
    added = pd.DataFrame(
        {
            "water_year": result.extended["water_year"].astype(int),
            "peak_flow_cfs": result.extended["flow_cfs"].astype(float),
            "qualification_code": None,
            "peak_date": None,
            "source": "extended",
        }
    )
    both = pd.concat([frame.assign(source="observed"), added], ignore_index=True)
    return result, both.sort_values("water_year").reset_index(drop=True)


def _extension_blocks(result: Any, options: Dict[str, Any]) -> List[Block]:
    """What the record extension did, in a paragraph, its statistics and any cautions."""
    index = options["extension"]["index"]
    label = " ".join(
        str(x)
        for x in (
            f"USGS {index['site_no']}" if index.get("site_no") else None,
            index.get("site_name"),
        )
        if x
    )
    years = [int(y) for y in result.extended["water_year"]]
    blocks: List[Block] = [
        Heading(2, "Record extension"),
        Paragraph(
            f"The record was extended from the index gage{' ' + label if label else ''} "
            "following Bulletin 17C Appendix 8 (England and others, 2019). The mean and "
            "variance of the logarithms of the peaks were first improved with the index "
            "gage's longer record (the Matalas-Jacobs estimators), and the effective record "
            f"that improvement is worth, {result.effective_length_variance:.1f} years in the "
            f"variance, set the number of peaks to add: {result.n_extended}. They were "
            f"generated from the index gage's peaks of {min(years)} to {max(years)} by the "
            "line below, whose intercept and slope make the observed and generated peaks "
            "together reproduce the improved mean and variance. The generated peaks stand "
            "for the information the index gage adds and are not estimates of what "
            "happened in those years; they enter the frequency analysis as if observed. "
            "The stationarity screen is not run on an extended record."
        ),
        Table(
            ("Item", "Value"),
            (
                ("Concurrent years", f"{result.n_concurrent}"),
                ("Index-gage years outside them", f"{result.n_additional}"),
                (
                    "Years at this gage the index gage lacks (not in the fit)",
                    f"{result.n_short_only}",
                ),
                ("Correlation of the concurrent log10 peaks", f"{result.correlation:.4f}"),
                ("Mean of log10 flow, concurrent years", f"{result.mean_short:.4f}"),
                ("Mean of log10 flow, improved", f"{result.mean_augmented:.4f}"),
                ("Standard deviation of log10 flow, concurrent years", f"{result.sd_short:.4f}"),
                (
                    "Standard deviation of log10 flow, improved",
                    f"{math.sqrt(result.variance_augmented):.4f}",
                ),
                ("Effective record length, mean (years)", f"{result.effective_length_mean:.1f}"),
                (
                    "Effective record length, variance (years)",
                    f"{result.effective_length_variance:.1f}",
                ),
                ("Peaks added", f"{result.n_extended}"),
                (
                    "Extension line",
                    f"log10 Q = {result.intercept:.4f} + {result.slope:.4f} "
                    f"(log10 Q_index - {result.mean_index_extended:.4f})",
                ),
            ),
            ("left", "right"),
            "Record extension from the index gage.",
        ),
    ]
    if result.notes:
        blocks.append(Bullets(tuple(f"**{note}**" for note in result.notes)))
    return blocks


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
    observed = frame
    extension = None
    if options.get("extension") is not None:
        extension, frame = _extended_frame(options, frame)
    screened = bool(options.get("stationarity", True)) and extension is None
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
        stationarity=screened,
        stationarity_figure=screened,
    )
    if result["error"]:
        raise ReportError(f"gaged_ffa: the frequency analysis failed: {result['error']}")
    p = result["parameters"]
    b17c = result["b17c"]
    site_name = options.get("site_name")
    label = " ".join(x for x in (f"USGS {site_no}" if site_no else None, site_name) if x)
    title = f"Flood frequency at {label}" if label else "Flood frequency"
    if extension is not None:
        title += ", record extended"
    blocks: List[Block] = []

    # Station and record
    year_values = observed["water_year"].to_numpy(dtype=int)
    flow_values = observed["peak_flow_cfs"].to_numpy(dtype=float)
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
            ("Annual peaks", f"{len(observed)}"),
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
    if extension is not None:
        record.append(("Peaks added by record extension", f"{extension.n_extended}"))
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
    if extension is not None:
        blocks.extend(_extension_blocks(extension, options))

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
        if extension is not None:
            columns.append(("Source", "left", [str(s) for s in frame["source"]]))
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


# -- engine_check: the fit against a reference computation ---------------------------------


_REFERENCE_KEYS = (
    "source",
    "computed",
    "record",
    "parameters",
    "low_outlier_count",
    "low_outlier_threshold_cfs",
    "quantiles",
)
_REFERENCE_QUANTILE_KEYS = ("aep", "flow_cfs", "lower_cfs", "upper_cfs")
#: The fitted parameters a reference may give, with the name each prints under.
_REFERENCE_PARAMETERS = {
    "mean_log": "Mean of log10 flow",
    "std_log": "Standard deviation of log10 flow",
    "skew_at_site": "Station skew",
    "skew_weighted": "Skew used",
}
_TOLERANCES = {"tolerance_pct": 1.0, "parameter_tolerance_pct": 0.5, "ci_tolerance_pct": 2.0}
_NO_EXTENSION = (
    "engine_check: the reference fit is made by flowfreq's Fortran extension (USGS peakfq "
    "8.1.0), which is not built in this environment. Build it with `python "
    "build_fortran/build.py` (gfortran and meson), or give the numbers of a reference "
    "computation in `reference`."
)


def _finite(owner: str, name: str, value: Any) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        number = math.nan
    if isinstance(value, bool) or not math.isfinite(number):
        raise ReportError(f"{owner}: {name} must be a number, got {value!r}")
    return number


def _record_digest(gage: Mapping[str, Any]) -> str:
    """A digest of the gage inputs that decide the fit, so a stored reference names its record."""
    frame = _peak_frame(gage["peaks"])
    low = gage.get("low_outlier_threshold")
    held = {
        "years": [int(y) for y in frame["water_year"]],
        "flows": [float(q) for q in frame["peak_flow_cfs"]],
        "codes": [None if c is None else str(c) for c in frame["qualification_code"]],
        "skew": [
            None if gage.get(key) is None else float(gage[key])
            for key in ("regional_skew", "regional_skew_se")
        ],
        "flags": [
            bool(gage.get(key, default))
            for key, default in (
                ("station_skew_only", False),
                ("use_default_skew", False),
                ("apply_peak_codes", True),
            )
        ],
        "thresholds": [
            [int(t["start_year"]), int(t["end_year"]), float(t["threshold_cfs"])]
            for t in _threshold_rows(gage.get("perception_thresholds"))
        ],
        "low_outlier_threshold": None if low is None else float(low),
    }
    text = json.dumps(held, sort_keys=True)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def _live_reference(b17c: Any, digest: str) -> Dict[str, Any]:
    """The record through peakfq's ``emafitpr``, as the plain values the options store."""
    try:
        live = b17c.fortran_reference()
    except ImportError as error:
        raise ReportError(_NO_EXTENSION) from error
    rows = []
    for aep in sorted(live.quantiles, reverse=True):
        lower, upper = live.confidence_intervals[aep]
        rows.append(
            {
                "aep": float(aep),
                "flow_cfs": float(live.quantiles[aep]),
                "lower_cfs": float(lower),
                "upper_cfs": float(upper),
            }
        )
    return {
        "source": (
            "USGS peakfq 8.1.0, its Fortran routine emafitpr called on the same peaks, "
            f"thresholds and regional skew by flowfreq {_version()}"
        ),
        "computed": _dt.date.today().isoformat(),
        "record": digest,
        "parameters": {name: float(live.parameters[name]) for name in _REFERENCE_PARAMETERS},
        "low_outlier_count": int(live.low_outlier_count),
        "low_outlier_threshold_cfs": float(live.low_outlier_threshold),
        "quantiles": rows,
    }


def _reference_result(given: Any) -> Any:
    """A stored or analyst-given reference as a ``ReferenceResult``, validated."""
    from ..validation.reference import ReferenceResult

    owner = "engine_check"
    if not isinstance(given, dict):
        raise ReportError(f"{owner}: reference must be an object")
    extra = sorted(set(given) - set(_REFERENCE_KEYS))
    if extra:
        raise ReportError(
            f"{owner}: unknown reference key(s) {', '.join(extra)}; "
            f"allowed: {', '.join(_REFERENCE_KEYS)}"
        )
    source = given.get("source")
    if not isinstance(source, str) or not source.strip():
        raise ReportError(
            f"{owner}: reference needs a source: the program, its version and where the "
            "numbers were taken from"
        )
    parameters = given.get("parameters") or {}
    if not isinstance(parameters, dict):
        raise ReportError(f"{owner}: reference parameters must be an object")
    unknown = sorted(set(parameters) - set(_REFERENCE_PARAMETERS))
    if unknown:
        raise ReportError(
            f"{owner}: unknown reference parameter(s) {', '.join(unknown)}; "
            f"allowed: {', '.join(_REFERENCE_PARAMETERS)}"
        )
    rows = given.get("quantiles")
    if not isinstance(rows, list) or not rows:
        raise ReportError(
            f"{owner}: reference needs quantiles: a list of objects with aep and flow_cfs, "
            "and lower_cfs and upper_cfs when it has 90-percent confidence limits"
        )
    quantiles: Dict[float, float] = {}
    limits: Dict[float, Tuple[float, float]] = {}
    for row in rows:
        if not isinstance(row, dict) or not {"aep", "flow_cfs"} <= set(row):
            raise ReportError(f"{owner}: each reference quantile needs aep and flow_cfs")
        stray = sorted(set(row) - set(_REFERENCE_QUANTILE_KEYS))
        if stray:
            raise ReportError(
                f"{owner}: unknown reference quantile key(s) {', '.join(stray)}; "
                f"allowed: {', '.join(_REFERENCE_QUANTILE_KEYS)}"
            )
        aep = _positive(owner, "a reference aep", row["aep"])
        if aep >= 1.0:
            raise ReportError(
                f"{owner}: a reference aep is a fraction between 0 and 1, got {row['aep']!r}"
            )
        if aep in quantiles:
            raise ReportError(f"{owner}: the reference gives AEP {aep:g} more than once")
        quantiles[aep] = _positive(owner, f"the reference flow at AEP {aep:g}", row["flow_cfs"])
        bounds = [row.get("lower_cfs"), row.get("upper_cfs")]
        if any(b is not None for b in bounds):
            if any(b is None for b in bounds):
                raise ReportError(
                    f"{owner}: the reference at AEP {aep:g} needs lower_cfs and upper_cfs "
                    "together"
                )
            limits[aep] = (
                _positive(owner, f"the reference lower limit at AEP {aep:g}", bounds[0]),
                _positive(owner, f"the reference upper limit at AEP {aep:g}", bounds[1]),
            )
    count = given.get("low_outlier_count")
    if count is not None and (isinstance(count, bool) or not isinstance(count, int) or count < 0):
        raise ReportError(
            f"{owner}: reference low_outlier_count must be a whole number, got {count!r}"
        )
    threshold = given.get("low_outlier_threshold_cfs")
    if (
        threshold is not None
        and _finite(owner, "reference low_outlier_threshold_cfs", threshold) < 0
    ):
        raise ReportError(
            f"{owner}: reference low_outlier_threshold_cfs must not be negative, got {threshold!r}"
        )
    return ReferenceResult(
        source=source.strip(),
        low_outlier_count=-1 if count is None else int(count),
        low_outlier_threshold=-1.0 if threshold is None else float(threshold),
        parameters={
            name: _finite(owner, f"the reference {name}", value)
            for name, value in parameters.items()
        },
        quantiles=quantiles,
        confidence_intervals=limits,
    )


def _difference(percent: float) -> str:
    """A percent difference to three decimals, without printing a small one as zero."""
    return "<0.001" if 0.0 < percent < 0.0005 else f"{percent:.3f}"


def _build_engine_check(given: Dict[str, Any]) -> SectionContent:
    from ..validation.comparisons import FrequencyComparator

    owner = "engine_check"
    options = dict(given)
    gage, b17c = _gaged_ffa(_nested("gaged_ffa", "gage", options["gage"], owner))
    tolerances = {
        name: default if options.get(name) is None else _positive(owner, name, options[name])
        for name, default in _TOLERANCES.items()
    }
    digest = _record_digest(gage.options)
    stored = options.get("reference")
    if stored is None:
        stored = _live_reference(b17c, digest)
    reference = _reference_result(stored)
    if stored.get("record") is not None and stored["record"] != digest:
        raise ReportError(
            f"{owner}: the stored reference was computed for other peaks, thresholds or "
            "regional skew than the gage has now. Remove reference to compute it again."
        )

    # flowfreq at the reference's AEPs, so the two are compared at the same probabilities.
    aeps = np.array(sorted(reference.quantiles, reverse=True), dtype=float)
    flows = b17c.compute_quantiles(aep=aeps)["flow_cfs"].to_numpy(dtype=float)
    bounds = b17c.compute_confidence_limits(aep=aeps)
    lower = bounds["lower_5pct"].to_numpy(dtype=float)
    upper = bounds["upper_5pct"].to_numpy(dtype=float)
    native_flow = {float(a): float(q) for a, q in zip(aeps, flows)}
    native_limits = {float(a): (float(lo), float(hi)) for a, lo, hi in zip(aeps, lower, upper)}
    native_parameters = b17c.to_comparison_dict()["parameters"]
    comparator = FrequencyComparator(
        tolerance_pct=tolerances["tolerance_pct"],
        parameter_tolerance_pct=tolerances["parameter_tolerance_pct"],
        ci_tolerance_pct=tolerances["ci_tolerance_pct"],
    )
    comparison = comparator.compare(
        {
            "parameters": native_parameters,
            "quantiles": native_flow,
            "confidence_intervals": native_limits,
        },
        reference,
    )
    results = b17c.results
    native_count = int(results.n_low_outliers or 0)
    native_threshold = float(results.low_outlier_threshold or 0.0)
    counts_differ = reference.low_outlier_count >= 0 and reference.low_outlier_count != native_count

    def within(ok: bool) -> str:
        return "yes" if ok else "no"

    # The verdict
    measured: List[str] = []
    outside: List[str] = []
    if comparison.quantile_diffs:
        worst = max(comparison.quantile_diffs.values())
        measured.append(
            f"{_difference(worst)} percent among the {len(comparison.quantile_diffs)} flood "
            f"quantiles (tolerance {comparator.tolerance_pct:g} percent)"
        )
        if worst > comparator.tolerance_pct:
            outside.append("the flood quantiles")
    if comparison.ci_diffs:
        worst = max(comparison.ci_diffs.values())
        measured.append(
            f"{_difference(worst)} percent among their 90-percent confidence limits "
            f"(tolerance {comparator.ci_tolerance_pct:g} percent)"
        )
        if worst > comparator.ci_tolerance_pct:
            outside.append("the confidence limits")
    if comparison.parameter_diffs:
        worst = max(comparison.parameter_diffs.values())
        measured.append(
            f"{_difference(worst)} percent in the mean and standard deviation of log10 flow "
            f"(tolerance {comparator.parameter_tolerance_pct:g} percent)"
        )
        if worst > comparator.parameter_tolerance_pct:
            outside.append("the mean or standard deviation")
    if comparison.skew_diffs:
        worst = max(comparison.skew_diffs.values())
        measured.append(
            f"{worst:.4f} in the skew coefficients (tolerance "
            f"{comparator.skew_tolerance_abs:g}, in skew units)"
        )
        if worst > comparator.skew_tolerance_abs:
            outside.append("the skew")
    if counts_differ:
        outside.append(
            f"the number of values censored as low outliers (flowfreq {native_count}, the "
            f"reference {reference.low_outlier_count})"
        )
    agree = comparison.passed and not counts_differ
    prefix = "Flood frequency at "
    label = gage.title[len(prefix) :] if gage.title.startswith(prefix) else ""
    computed = stored.get("computed")
    sentences = [
        f"The Bulletin 17C fit{' at ' + label if label else ''} was checked against a "
        f"reference computation of the same record: {reference.source}"
        + (f" ({computed})." if computed else "."),
        (
            "The two agree within the tolerances set."
            if agree
            else "**The two do not agree within the tolerances set.**"
        ),
        "The largest differences are " + "; ".join(measured) + ".",
    ]
    if outside:
        sentences.append("Outside tolerance: " + "; ".join(outside) + ".")
    elif reference.low_outlier_count >= 0:
        sentences.append(
            f"Both censor {native_count} value(s) as low outliers."
            if native_count
            else "Neither censors a value as a low outlier."
        )
    blocks: List[Block] = [Paragraph(" ".join(sentences))]

    # Fitted parameters
    fitted: List[Tuple[str, ...]] = []
    for name, printed in _REFERENCE_PARAMETERS.items():
        if name not in reference.parameters or name not in native_parameters:
            continue
        if name in comparison.skew_diffs:
            gap = comparison.skew_diffs[name]
            tolerance = comparator.skew_tolerance_abs
            shown, allowed = f"{gap:.4f}", f"{tolerance:g}"
        else:
            gap = comparison.parameter_diffs[name]
            tolerance = comparator.parameter_tolerance_pct
            shown, allowed = f"{_difference(gap)} percent", f"{tolerance:g} percent"
        fitted.append(
            (
                printed,
                f"{native_parameters[name]:.4f}",
                f"{reference.parameters[name]:.4f}",
                shown,
                allowed,
                within(gap <= tolerance),
            )
        )
    if reference.low_outlier_count >= 0:
        fitted.append(
            (
                "Values censored as low outliers",
                f"{native_count}",
                f"{reference.low_outlier_count}",
                f"{abs(native_count - reference.low_outlier_count)}",
                "0",
                within(not counts_differ),
            )
        )
    if reference.low_outlier_threshold >= 0:
        fitted.append(
            (
                "Low-outlier threshold (cfs)",
                cfs(native_threshold) if native_threshold > 0 else "none",
                (
                    cfs(reference.low_outlier_threshold)
                    if reference.low_outlier_threshold > 0
                    else "none"
                ),
                "",
                "",
                "not judged",
            )
        )
    if fitted:
        blocks.append(
            Table(
                (
                    "Parameter",
                    "flowfreq",
                    "Reference",
                    "Difference",
                    "Tolerance",
                    "Within tolerance",
                ),
                tuple(fitted),
                ("left", "right", "right", "right", "right", "left"),
                "Fitted parameters from flowfreq and from the reference.",
            )
        )

    # Quantiles
    blocks.append(
        Table(
            (
                "Annual exceedance probability (percent)",
                "flowfreq (cfs)",
                "Reference (cfs)",
                "Difference (percent)",
                "Within tolerance",
            ),
            tuple(
                (
                    f"{aep * 100:.4g}",
                    cfs(native_flow[aep]),
                    cfs(reference.quantiles[aep]),
                    _difference(comparison.quantile_diffs[aep]),
                    within(comparison.quantile_diffs[aep] <= comparator.tolerance_pct),
                )
                for aep in native_flow
            ),
            ("right", "right", "right", "right", "left"),
            "Flood quantiles from flowfreq and from the reference.",
        )
    )
    if comparison.ci_diffs:
        blocks.append(
            Table(
                (
                    "Annual exceedance probability (percent)",
                    "flowfreq lower limit (cfs)",
                    "flowfreq upper limit (cfs)",
                    "Reference lower limit (cfs)",
                    "Reference upper limit (cfs)",
                    "Larger difference (percent)",
                    "Within tolerance",
                ),
                tuple(
                    (
                        f"{aep * 100:.4g}",
                        cfs(native_limits[aep][0]),
                        cfs(native_limits[aep][1]),
                        cfs(reference.confidence_intervals[aep][0]),
                        cfs(reference.confidence_intervals[aep][1]),
                        _difference(comparison.ci_diffs[aep]),
                        within(comparison.ci_diffs[aep] <= comparator.ci_tolerance_pct),
                    )
                    for aep in native_flow
                    if aep in comparison.ci_diffs
                ),
                ("right", "right", "right", "right", "right", "right", "left"),
                "90-percent confidence limits from flowfreq and from the reference.",
            )
        )
    if bool(options.get("detail", True)):
        blocks.append(Heading(2, "Analysis at the gage"))
        blocks.extend(_demote(gage.blocks))
    options["gage"] = gage.options
    options["reference"] = stored
    return SectionContent(
        title=f"Computation check at {label}" if label else "Computation check",
        blocks=blocks,
        options=options,
        references=list(gage.references),
    )


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


#: How far apart, as a fraction of the larger, two AEPs may be and still be one row. A
#: Bulletin 17C fit tabulates the 1.5-year flood at 0.67 and the 2.33-year flood at 0.43,
#: where a regression report has 0.667 and 0.429. Relative, so that 0.005 and 0.002 are
#: never taken for each other.
AEP_MATCH = 0.0075


def _flow_table(flows: pd.DataFrame) -> Dict[float, Any]:
    """An estimate's rows that have a flow, by AEP."""
    return {float(r["aep"]): r for _, r in flows[flows["flow_cfs"].notna()].iterrows()}


def _row_at(table: Mapping[float, Any], aep: float) -> Any:
    """The estimate's row at ``aep``, or at the AEP nearest it within ``AEP_MATCH``; else None."""
    if not table:
        return None
    nearest = min(table, key=lambda a: abs(a - aep))
    return table[nearest] if abs(nearest - aep) <= AEP_MATCH * max(nearest, aep) else None


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
                if not any(abs(aep - c) <= AEP_MATCH * max(aep, c) for c in candidates):
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
        keep = [a for a in have if min(published) - 1e-12 <= a <= max(published) + 1e-12]
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
        if any(abs(aep - p) < 1e-9 for p in published):
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
                "extension",
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
        SectionKind(
            name="transposed",
            summary=(
                "A gage's Bulletin 17C quantiles transposed to an ungaged site on the same "
                "stream by the drainage-area ratio, with a cited regression exponent: the "
                "arithmetic, the caveats, then the analysis at the gage."
            ),
            required=("donor", "target_area_sqmi", "exponents"),
            optional=(
                "donor_area_sqmi",
                "target_name",
                "max_aep",
                "allow_out_of_range",
                "donor_site_class",
                "target_site_class",
                "allow_regulation_mismatch",
                "detail",
            ),
            build=_build_transposed,
        ),
        SectionKind(
            name="engine_check",
            summary=(
                "A gage's Bulletin 17C fit checked against a reference computation of the "
                "same record: USGS peakfq 8.1.0 through flowfreq's Fortran extension, whose "
                "numbers are then stored in the specification, or numbers the analyst gives. "
                "Parameters, quantiles and confidence limits side by side, with the "
                "differences and whether each is within tolerance."
            ),
            required=("gage",),
            optional=(
                "reference",
                "tolerance_pct",
                "parameter_tolerance_pct",
                "ci_tolerance_pct",
                "detail",
            ),
            build=_build_engine_check,
        ),
        SectionKind(
            name="future_flow",
            summary=(
                "Flood quantiles under a future condition: a gage's Bulletin 17C quantiles or "
                "a regression estimate multiplied by change factors the analyst supplies and "
                "cites, with which factors were published, interpolated or held."
            ),
            required=("factor_set",),
            optional=("gage", "regression", "aeps", "detail"),
            build=_build_future_flow,
        ),
        SectionKind(
            name="comparison",
            summary=(
                "Several estimates of the flood quantiles at one site, side by side: each a "
                "gage's Bulletin 17C fit, a transposition from a gage or a regional "
                "regression, with its limits, the ratio of the largest to the smallest, and "
                "each estimate's own analysis below."
            ),
            required=("estimates",),
            optional=("site_name", "aeps", "detail"),
            build=_build_comparison,
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
    "COMPARISON_CAPTION",
    "FUTURE_FLOW_CAPTION",
    "SECTION_KINDS",
    "SectionContent",
    "SectionKind",
    "build_section",
    "cfs",
    "download_peaks",
    "short_citation",
]
