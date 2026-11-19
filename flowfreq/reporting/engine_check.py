"""The ``engine_check`` section: the fit against a reference computation."""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import math
from typing import Any, Dict, List, Mapping, Tuple

import numpy as np
import pandas as pd

from ._base import SectionContent, _demote, _version, cfs
from .blocks import Block, Heading, Paragraph, Table
from .gaged import _gaged_ffa, _peak_frame, _threshold_rows
from .spec import ReportError
from .transposed import _positive
from .weighted import _nested

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


def _extension_digest(extension: Any) -> Any:
    """What of a record extension decides the fit: the index peaks and which years are used."""
    if not isinstance(extension, dict) or not isinstance(extension.get("index"), dict):
        return None
    peaks = extension["index"].get("peaks") or {}
    return {
        "years": [int(y) for y in peaks.get("water_year") or []],
        "flows": [float(q) for q in peaks.get("peak_flow_cfs") or []],
        "n_extend": extension.get("n_extend"),
        "chosen": extension.get("years"),
    }


def _record_digest(gage: Mapping[str, Any]) -> str:
    """A digest of the gage inputs that decide the fit, so a stored reference names its record."""
    frame = _peak_frame(gage["peaks"])
    low = gage.get("low_outlier_threshold")
    held = {
        "years": [int(y) for y in frame["water_year"]],
        "flows": [float(q) for q in frame["peak_flow_cfs"]],
        "codes": [None if pd.isna(c) else str(c) for c in frame["qualification_code"]],
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
    extension = _extension_digest(gage.get("extension"))
    if extension is not None:
        # Only when there is one, so a reference stored before extensions existed still matches.
        held["extension"] = extension
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
