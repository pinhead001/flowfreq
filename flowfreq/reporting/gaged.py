"""The ``gaged_ffa`` section: a Bulletin 17C fit of a gage's annual peaks."""

from __future__ import annotations

import copy
import datetime as _dt
import json
import math
import threading
from collections import OrderedDict
from dataclasses import replace
from typing import Any, Dict, List, Mapping, Optional, Tuple

import numpy as np
import pandas as pd

from ._base import (
    _PEAK_KEYS,
    BULLETIN_17C,
    SectionContent,
    _demote,
    _figure_png,
    _version,
    cfs,
    short_citation,
)
from .blocks import Block, Bullets, Figure, Heading, Paragraph, Table, parse_markdown, without_title
from .spec import ReportError

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
    values = frame["peak_flow_cfs"].to_numpy()
    if not np.all(np.isfinite(values)):
        raise ReportError("gaged_ffa: peak_flow_cfs holds a non-finite value")
    if np.any(values < 0):
        raise ReportError("gaged_ffa: peak_flow_cfs holds a negative value")
    if n < 3 or int(np.sum(values > 0)) < 3:
        raise ReportError(
            f"gaged_ffa: a frequency analysis needs at least 3 peaks above zero; the record "
            f"has {int(np.sum(values > 0))} of {n}"
        )
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


_EXTENSION_KEYS = ("index", "n_extend", "years", "split_year")


_INDEX_KEYS = ("site_no", "site_name", "peaks")


def _extended_frame(options: Dict[str, Any], frame: pd.DataFrame) -> Tuple[Any, pd.DataFrame]:
    """The record extended from an index gage (Bulletin 17C Appendix 8), and the extension.

    Downloads the index gage's peaks when the options carry none, and writes them back.
    The frame returned has the observed peaks and the generated ones, with a ``source``
    column saying which.
    """
    from ..peak_extension import concurrent_relation, extend_peaks

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
    count = extension.get("n_extend")
    if count is not None and (isinstance(count, bool) or not isinstance(count, int) or count < 2):
        raise ReportError("gaged_ffa: extension n_extend must be a whole number, 2 or more")
    chosen = extension.get("years")
    if chosen is not None and (
        not isinstance(chosen, list)
        or not chosen
        or any(isinstance(y, bool) or not isinstance(y, int) for y in chosen)
    ):
        raise ReportError("gaged_ffa: extension years must be a list of water years")
    split = extension.get("split_year")
    if split is not None and (isinstance(split, bool) or not isinstance(split, int)):
        raise ReportError("gaged_ffa: extension split_year must be a water year")
    index = dict(extension["index"])
    if index.get("peaks") is None:
        if not index.get("site_no"):
            raise ReportError("gaged_ffa: the extension's index gage needs site_no or peaks")
        fetched = download_peaks(str(index["site_no"]).strip())
        index["peaks"] = fetched["peaks"]
        if index.get("site_name") is None and fetched.get("site_name") is not None:
            index["site_name"] = fetched["site_name"]
    index_frame = _peak_frame(index["peaks"])
    short = pd.Series(frame["peak_flow_cfs"].to_numpy(), index=frame["water_year"].to_numpy())
    long = pd.Series(
        index_frame["peak_flow_cfs"].to_numpy(), index=index_frame["water_year"].to_numpy()
    )
    try:
        result = extend_peaks(
            short, long, n_extend=extension.get("n_extend"), years=extension.get("years")
        )
        relation = concurrent_relation(short, long, split_year=extension.get("split_year"))
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
    return (result, relation), both.sort_values("water_year").reset_index(drop=True)


def _extension_blocks(result: Any, relation: pd.DataFrame, options: Dict[str, Any]) -> List[Block]:
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
    chosen = options["extension"].get("split_year") is not None
    blocks.append(
        Paragraph(
            "Record extension assumes one fixed relation between the two gages. The table "
            "below sets the two parts of the concurrent period side by side, divided at "
            + ("the year the analyst gave" if chosen else "its middle")
            + ": the median of this gage's peak over the index gage's in the same year, "
            "the spread of that ratio, and the correlation of the logarithms. A ratio or a "
            "correlation that differs between the parts means the relation moved."
        )
    )
    blocks.append(
        Table(
            (
                "Water years",
                "Concurrent years",
                "Median ratio to the index gage",
                "Standard deviation of log10 ratio",
                "Correlation of log10 peaks",
            ),
            tuple(
                (
                    str(row["period"]).replace("-", " to "),
                    f"{int(row['n'])}",
                    f"{row['median_ratio']:.3f}",
                    "-" if pd.isna(row["ratio_log_sd"]) else f"{row['ratio_log_sd']:.4f}",
                    "-" if pd.isna(row["correlation"]) else f"{row['correlation']:.3f}",
                )
                for _, row in relation.iterrows()
            ),
            ("left", "right", "right", "right", "right"),
            "Relation of the two gages' concurrent peaks, by part of the concurrent period.",
        )
    )
    return blocks


def _build_gaged_ffa(given: Dict[str, Any]) -> SectionContent:
    return _gaged_ffa(given)[0]


#: How many gaged fits :func:`_gaged_ffa` keeps, least recently used dropped first.
FIT_CACHE_SIZE: int = 32

_FIT_CACHE: OrderedDict[str, Tuple[SectionContent, Any]] = OrderedDict()
_FIT_CACHE_LOCK = threading.Lock()


def clear_fit_cache() -> None:
    """Forget every memoised gaged fit.

    A report fits each gage once: :func:`_gaged_ffa` keeps the last
    :data:`FIT_CACHE_SIZE` fits, keyed on their options. Clearing makes the next build of
    any gaged section fit again.

    Returns
    -------
    None
    """
    with _FIT_CACHE_LOCK:
        _FIT_CACHE.clear()


def _fit_key(options: Mapping[str, Any]) -> Optional[str]:
    """The canonical form of a gage's options, or None when JSON cannot encode them."""
    try:
        return json.dumps(options, sort_keys=True)
    except (TypeError, ValueError):
        return None


def _cache_get(key: Optional[str]) -> Optional[Tuple[SectionContent, Any]]:
    if key is None:
        return None
    with _FIT_CACHE_LOCK:
        hit = _FIT_CACHE.get(key)
        if hit is None:
            return None
        _FIT_CACHE.move_to_end(key)
        content, b17c = hit
        return copy.deepcopy(content), b17c


def _cache_put(keys: List[Optional[str]], entry: Tuple[SectionContent, Any]) -> None:
    with _FIT_CACHE_LOCK:
        for key in keys:
            if key is None:
                continue
            _FIT_CACHE[key] = entry
            _FIT_CACHE.move_to_end(key)
        while len(_FIT_CACHE) > FIT_CACHE_SIZE:
            _FIT_CACHE.popitem(last=False)


def _gaged_ffa(given: Dict[str, Any]) -> Tuple[SectionContent, Any]:
    """The gaged section and the fitted ``Bulletin17C`` behind it.

    Fits are memoised (see :func:`clear_fit_cache`) under the options as resolved, which carry
    the peaks, and under those given when they carry ``peaks`` too, so the sections nested from
    a previous build's resolved options reuse the fit. Options without peaks are downloaded and
    fitted each time. The section returned is the caller's own copy; the
    ``Bulletin17C`` object may be shared between callers and must be treated as read-only.
    ``given`` is not changed.
    """
    cacheable = given.get("peaks") is not None
    key = _fit_key(given) if cacheable else None
    hit = _cache_get(key)
    if hit is not None:
        return hit
    content, b17c = _fit_gaged_ffa(given)
    resolved = _fit_key(content.options)
    if key is not None or resolved is not None:
        _cache_put([key, resolved], (copy.deepcopy(content), b17c))
    return content, b17c


def _fit_gaged_ffa(given: Dict[str, Any]) -> Tuple[SectionContent, Any]:
    """Build the gaged section and fit it, uncached."""
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
    relation = None
    if options.get("extension") is not None:
        (extension, relation), frame = _extended_frame(options, frame)
    screened = bool(options.get("stationarity", True)) and extension is None
    _resolve_skew(options)
    thresholds = _threshold_rows(options.get("perception_thresholds"))
    apply_codes = bool(options.get("apply_peak_codes", True))
    codes = [None if pd.isna(c) else c for c in frame["qualification_code"].tolist()]
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
    if extension is not None and relation is not None:
        blocks.extend(_extension_blocks(extension, relation, options))

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
                "Recurrence interval (years)",
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
        blocks.extend(
            (
                replace(b, caption="Stationarity tests of the annual peaks.")
                if isinstance(b, Table) and not b.caption and b.header[:1] == ("Test",)
                else b
            )
            for b in _demote(body)
        )
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
        dates = ["" if pd.isna(d) else str(d) for d in frame["peak_date"]]
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


__all__ = ["FIT_CACHE_SIZE", "clear_fit_cache", "download_peaks"]
