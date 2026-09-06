"""
flowfreq.workflow - High-level Bulletin 17C analysis entry points.

One call from annual peaks to a fitted frequency curve, plus the skew-variant
helpers that go with it. This is the layer a consumer wants when it does not
want to assemble :class:`~flowfreq.bulletin17c.Bulletin17C`, quantiles and
confidence limits by hand.

Everything here returns plain numbers and DataFrames of numbers. Turning those
into labelled, rounded, comma-separated strings is presentation and belongs to
the caller.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .bulletin17c import Bulletin17C
from .core import EMAParameters, FrequencyResults, kfactor_array

if TYPE_CHECKING:  # pragma: no cover - annotations only
    # Same layering reason as flowfreq.bulletin17c.Bulletin17C.validate: a
    # runtime import here would make every `import flowfreq.workflow` pull in
    # the validation subsystem, which flowfreq/__init__ leaves opt-in.
    # compare_engines imports what it needs lazily at call time.
    from .validation.comparisons import ComparisonResult
    from .validation.reference import ReferenceResult

logger = logging.getLogger(__name__)

#: Nationwide B17C default generalized skew (England et al. 2019).
B17C_DEFAULT_SKEW: float = -0.302

#: Return intervals reported by :func:`run_ffa` and :func:`compute_skew_tables`.
DEFAULT_RETURN_INTERVALS: List[float] = [1.5, 2, 5, 10, 25, 50, 100, 200, 500]

#: The same series as annual exceedance probabilities.
#: = [0.667, 0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002]
DEFAULT_AEP: List[float] = [1 / ri for ri in DEFAULT_RETURN_INTERVALS]

#: Canonical skew option labels, in the order a report presents them.
SKEW_OPTIONS: List[str] = ["Station Skew", "Weighted Skew", "Regional Skew"]


def _low_outlier_source(override: Optional[float]) -> str:
    """Describe where the reported PILF threshold came from.

    Parameters
    ----------
    override : float or None
        The user's requested threshold, if any.

    Returns
    -------
    str
        ``"MGBT"`` or ``"override"``. Both EMA and MOM censor on a supplied
        threshold now, so the label needs no method-specific caveat.
    """
    if override is None:
        return "MGBT"
    return "override"


def run_ffa(
    peak_flows: np.ndarray,
    water_years: np.ndarray,
    regional_skew: float = B17C_DEFAULT_SKEW,
    regional_skew_se: float = 0.55,
    perception_thresholds: Optional[List[dict]] = None,
    low_outlier_threshold_override: Optional[float] = None,
) -> dict:
    """Run Bulletin 17C flood frequency analysis.

    Fits by EMA, falling back to MOM only when EMA fails to converge *and* no
    perception thresholds are in play -- MOM cannot represent censored
    intervals, so a threshold-bearing record keeps its non-converged EMA fit
    rather than silently losing the thresholds.

    Errors are returned, not raised: any failure comes back as a string under
    ``error`` with the other keys left at their empty defaults. That suits an
    interactive caller that wants to show the message rather than crash. A
    caller that would rather have an exception should check ``error`` and
    raise its own.

    Parameters
    ----------
    peak_flows : np.ndarray
        Annual peak flows in cfs.
    water_years : np.ndarray
        Corresponding water years.
    regional_skew : float
        Regional skew coefficient. Defaults to the nationwide B17C value.
    regional_skew_se : float
        Regional skew standard error.
    perception_thresholds : list of dict, optional
        Each dict has keys ``start_year``, ``end_year``, ``threshold_cfs``
        (legacy) or ``lower_cfs`` / ``upper_cfs``.  Converts to the
        ``Dict[Tuple[int,int], float]`` format expected by
        :class:`~flowfreq.bulletin17c.Bulletin17C` and passed to EMA so that
        years in each period without a recorded peak are treated as
        left-censored observations (peak < threshold).
    low_outlier_threshold_override : float, optional
        User-supplied PILF threshold (cfs).  When > 0, overrides the MGBT
        result and censors all peaks below this value.  The threshold actually
        applied, its source and the resulting PILF count come back under
        ``parameters`` so a caller can show which cut produced the fit.

    Returns
    -------
    dict
        Keys: b17c, converged, method, parameters, quantile_df, error.

    Examples
    --------
    >>> result = run_ffa(peak_flows, water_years)
    >>> result["quantile_df"]["Flow (cfs)"]

    >>> result = run_ffa(peak_flows, water_years, regional_skew=-0.05,
    ...                  low_outlier_threshold_override=500.0)
    """
    result = {
        "b17c": None,
        "converged": False,
        "method": None,
        "parameters": {},
        "quantile_df": pd.DataFrame(),
        "error": None,
    }

    try:
        # Convert list of threshold dicts → Dict[Tuple[int,int], float]
        pt_dict: Optional[Dict[Tuple[int, int], float]] = None
        if perception_thresholds:
            pt_dict = {
                (int(t["start_year"]), int(t["end_year"])): float(t["threshold_cfs"])
                for t in perception_thresholds
                if float(t.get("threshold_cfs", 0)) > 0
            } or None

        lo_override = (
            float(low_outlier_threshold_override)
            if low_outlier_threshold_override and low_outlier_threshold_override > 0
            else None
        )

        b17c = Bulletin17C(
            peak_flows=peak_flows,
            water_years=water_years,
            regional_skew=regional_skew,
            regional_skew_mse=regional_skew_se**2,
            perception_thresholds=pt_dict,
            user_low_outlier_threshold=lo_override,
        )

        b17c.run_analysis(method="ema")
        method = "ema"
        converged = bool(b17c.results.ema_converged)

        # Only fall back to MOM when no perception thresholds are in play — MOM has no
        # mechanism to incorporate censored intervals, so we keep the (non-converged)
        # EMA result when thresholds extend the record.
        if not converged and not pt_dict:
            logger.warning("EMA did not converge, falling back to MOM")
            b17c.run_analysis(method="mom")
            method = "mom"
            converged = True

        aep = np.array(DEFAULT_AEP)
        quantiles_df = b17c.compute_quantiles(aep=aep)
        ci_df = b17c.compute_confidence_limits(aep=aep)

        quantile_df = pd.DataFrame(
            {
                "Return Interval (yr)": DEFAULT_RETURN_INTERVALS,
                "AEP (%)": aep,
                "Flow (cfs)": quantiles_df["flow_cfs"].values,
                "Lower 90% CI": ci_df["lower_5pct"].values,
                "Upper 90% CI": ci_df["upper_5pct"].values,
            }
        )

        r = b17c.results
        result.update(
            {
                "b17c": b17c,
                "converged": converged,
                "method": method,
                "parameters": {
                    "mean_log": r.mean_log,
                    "std_log": r.std_log,
                    "skew_station": r.skew_station,
                    "skew_weighted": r.skew_weighted,
                    "skew_used": r.skew_used,
                    "regional_skew": regional_skew,
                    # The low-outlier cut and where it came from. Without these
                    # a caller could offer the override but never show its effect.
                    # Both EMA and MOM censor on it now, so the source is just
                    # whether it came from MGBT or the user.
                    "low_outlier_threshold": r.low_outlier_threshold,
                    "n_low_outliers": r.n_low_outliers,
                    "low_outlier_source": _low_outlier_source(lo_override),
                },
                "quantile_df": quantile_df,
            }
        )

    except Exception as e:
        logger.exception("FFA analysis failed")
        result["error"] = str(e)

    return result


def _skew_values_from_result(ffa_result: dict) -> Dict[str, Optional[float]]:
    """Return the three skew values stored in an ffa_result dict."""
    p = ffa_result.get("parameters", {})
    return {
        "Station Skew": p.get("skew_station"),
        "Weighted Skew": p.get("skew_weighted"),
        "Regional Skew": p.get("regional_skew"),
    }


def compute_skew_tables(
    ffa_result: dict,
    selected_labels: List[str],
) -> Dict[str, pd.DataFrame]:
    """Compute a raw quantile+CI table for each selected skew option.

    Uses the LP3 moments (mean_log, std_log) already fitted by EMA/MOM and
    substitutes the requested skew value to produce separate frequency tables
    without re-running the full analysis.

    Parameters
    ----------
    ffa_result : dict
        Output from :func:`run_ffa`.
    selected_labels : list[str]
        Subset of :data:`SKEW_OPTIONS`.

    Returns
    -------
    dict[str, pd.DataFrame]
        Maps label → DataFrame with columns:
        ``Return Interval (yr)``, ``AEP (%)``,
        ``Flow (cfs)``, ``Lower 90% CI``, ``Upper 90% CI``.
        Returns an empty dict if ffa_result has an error.
    """
    if ffa_result.get("error") or ffa_result.get("b17c") is None:
        return {}

    r = ffa_result["b17c"].results
    mean_log = r.mean_log
    std_log = r.std_log
    n = r.n_systematic or r.n_peaks

    skew_map = _skew_values_from_result(ffa_result)
    aep = np.array(DEFAULT_AEP)
    z_alpha = 1.6449  # norm.ppf(0.95): two-sided 90% CI

    tables: Dict[str, pd.DataFrame] = {}
    for label in selected_labels:
        skew_val = skew_map.get(label)
        if skew_val is None:
            continue

        K = kfactor_array(skew_val, aep)
        log_Q = mean_log + K * std_log
        Q = 10.0**log_Q

        var_factor = 1 / n + K**2 * (1 + 0.75 * skew_val**2) / (2 * (n - 1))
        se_log = std_log * np.sqrt(var_factor)
        lower = 10.0 ** (log_Q - z_alpha * se_log)
        upper = 10.0 ** (log_Q + z_alpha * se_log)

        tables[label] = pd.DataFrame(
            {
                "Return Interval (yr)": DEFAULT_RETURN_INTERVALS,
                "AEP (%)": aep,
                "Flow (cfs)": Q,
                "Lower 90% CI": lower,
                "Upper 90% CI": upper,
            }
        )

    return tables


def build_skew_curves_dict(
    ffa_result: dict,
    selected_labels: List[str],
) -> Dict[str, float]:
    """Return ``{label: skew_value}`` for the selected skew options.

    Intended for passing directly to
    :func:`flowfreq.freq_plot.plot_frequency_curve` as the
    ``skew_curves`` argument.

    Parameters
    ----------
    ffa_result : dict
        Output from :func:`run_ffa`.
    selected_labels : list[str]
        Skew labels the caller has selected.

    Returns
    -------
    dict[str, float]
        Empty dict (fall back to default) when no valid labels are found.
    """
    skew_map = _skew_values_from_result(ffa_result)
    return {lbl: skew_map[lbl] for lbl in selected_labels if skew_map.get(lbl) is not None}


@dataclass
class EngineComparisonReport:
    """The result of running one record through both engines and comparing them.

    ``docs/FORTRAN_ENGINE_DESIGN.md`` section 5: "the comparison ... is the
    feature; ``engine=`` alone leaves the user diffing two analyses by hand."
    Built by :func:`compare_engines`, which reuses the already-existing
    :meth:`~flowfreq.bulletin17c.Bulletin17C.validate` /
    :class:`~flowfreq.validation.comparisons.FrequencyComparator` machinery
    rather than duplicating comparison logic.

    Attributes
    ----------
    native : FrequencyResults
        The native engine's fitted result, evaluated at the same AEPs the
        Fortran side was.
    reference : ReferenceResult
        The live ``emafitpr`` output.
    comparison : ComparisonResult
        Per-field differences and pass/fail, from
        :class:`~flowfreq.validation.comparisons.FrequencyComparator`.
    site_name : str
        Carried through for the markdown title only.
    """

    native: FrequencyResults
    reference: "ReferenceResult"
    comparison: "ComparisonResult"
    site_name: str = ""

    @property
    def max_quantile_deviation_pct(self) -> float:
        """Largest quantile percent difference across every AEP compared.

        Deliberately not :attr:`ComparisonResult.max_diff_pct` -- that also
        folds in parameter and confidence-interval differences, which answer
        a different question than "how far apart are the flood quantiles
        themselves."
        """
        return max(self.comparison.quantile_diffs.values(), default=0.0)

    def to_markdown(self) -> str:
        """A markdown report suitable for a submittal appendix.

        Parameters, skews (in skew units, not percent -- see
        :class:`~flowfreq.validation.comparisons.ComparisonResult`), quantiles
        and confidence limits, each as a table of native vs. Fortran with the
        percent (or, for skew, absolute) difference already computed by
        :class:`~flowfreq.validation.comparisons.FrequencyComparator`.
        """
        n, r, c = self.native, self.reference, self.comparison
        lines: List[str] = []

        title = "Engine comparison"
        if self.site_name:
            title += f": {self.site_name}"
        lines.append(f"# {title}")
        lines.append("")
        status = "PASS" if c.passed else "FAIL"
        lines.append(
            f"**{status}** -- native EMA vs. USGS peakfq (`emafitpr`, the vendored Fortran "
            "reference)"
        )
        lines.append("")
        lines.append(
            f"- Max quantile deviation: **{self.max_quantile_deviation_pct:.3f}%** "
            f"(tolerance {c.tolerance_pct:.2f}%)"
        )
        lines.append(f"- {c.summary}")
        lines.append("")

        native_params = {"mean_log": n.mean_log, "std_log": n.std_log}
        if c.parameter_diffs:
            lines.append("## Parameters")
            lines.append("")
            lines.append("| Parameter | Native | Fortran | Diff % |")
            lines.append("|---|---:|---:|---:|")
            for key in sorted(c.parameter_diffs):
                nat_val = native_params.get(key)
                ref_val = r.parameters.get(key)
                nat_str = f"{nat_val:.6g}" if nat_val is not None else "n/a"
                ref_str = f"{ref_val:.6g}" if ref_val is not None else "n/a"
                lines.append(f"| {key} | {nat_str} | {ref_str} | {c.parameter_diffs[key]:.3f} |")
            lines.append("")

        if c.skew_diffs:
            native_skew = {"skew_at_site": n.skew_station, "skew_weighted": n.skew_weighted}
            lines.append("## Skew (skew units, not percent)")
            lines.append("")
            lines.append("| Parameter | Native | Fortran | Abs diff |")
            lines.append("|---|---:|---:|---:|")
            for key in sorted(c.skew_diffs):
                nat_val = native_skew.get(key)
                ref_val = r.parameters.get(key)
                nat_str = f"{nat_val:.4f}" if nat_val is not None else "n/a"
                ref_str = f"{ref_val:.4f}" if ref_val is not None else "n/a"
                lines.append(f"| {key} | {nat_str} | {ref_str} | {c.skew_diffs[key]:.4f} |")
            lines.append("")

        if c.quantile_diffs:
            native_q: Dict[float, float] = {}
            if n.quantiles is not None and not n.quantiles.empty:
                native_q = dict(zip(n.quantiles["aep"], n.quantiles["flow_cfs"]))
            lines.append("## Quantiles (cfs)")
            lines.append("")
            lines.append("| AEP | Return period (yr) | Native | Fortran | Diff % |")
            lines.append("|---:|---:|---:|---:|---:|")
            for aep in sorted(c.quantile_diffs, reverse=True):
                nat_val = native_q.get(aep)
                ref_val = r.quantiles.get(aep)
                nat_str = f"{nat_val:,.0f}" if nat_val is not None else "n/a"
                ref_str = f"{ref_val:,.0f}" if ref_val is not None else "n/a"
                lines.append(
                    f"| {aep:.4g} | {1 / aep:,.1f} | {nat_str} | {ref_str} | "
                    f"{c.quantile_diffs[aep]:.3f} |"
                )
            lines.append("")

        if c.ci_diffs:
            native_ci: Dict[float, Tuple[float, float]] = {}
            if n.confidence_limits is not None and not n.confidence_limits.empty:
                for _, row in n.confidence_limits.iterrows():
                    native_ci[float(row["aep"])] = (row["lower_5pct"], row["upper_5pct"])
            lines.append("## Confidence intervals (cfs)")
            lines.append("")
            lines.append(
                "| AEP | Native lower | Native upper | Fortran lower | Fortran upper | Diff % |"
            )
            lines.append("|---:|---:|---:|---:|---:|---:|")
            for aep in sorted(c.ci_diffs, reverse=True):
                nat_lo, nat_hi = native_ci.get(aep, (None, None))
                ref_lo, ref_hi = r.confidence_intervals.get(aep, (None, None))
                nat_lo_s = f"{nat_lo:,.0f}" if nat_lo is not None else "n/a"
                nat_hi_s = f"{nat_hi:,.0f}" if nat_hi is not None else "n/a"
                ref_lo_s = f"{ref_lo:,.0f}" if ref_lo is not None else "n/a"
                ref_hi_s = f"{ref_hi:,.0f}" if ref_hi is not None else "n/a"
                lines.append(
                    f"| {aep:.4g} | {nat_lo_s} | {nat_hi_s} | {ref_lo_s} | {ref_hi_s} | "
                    f"{c.ci_diffs[aep]:.3f} |"
                )
            lines.append("")

        return "\n".join(lines)


def compare_engines(
    peak_flows: np.ndarray,
    water_years: Optional[np.ndarray] = None,
    regional_skew: float = B17C_DEFAULT_SKEW,
    regional_skew_se: float = 0.55,
    historical_peaks: Optional[List[Tuple[int, float]]] = None,
    perception_thresholds: Optional[Dict[Tuple[int, int], float]] = None,
    user_low_outlier_threshold: Optional[float] = None,
    ema_params: Optional[EMAParameters] = None,
    aeps: Optional[np.ndarray] = None,
    site_name: str = "",
    tolerance_pct: float = 1.0,
    parameter_tolerance_pct: float = 0.5,
    ci_tolerance_pct: float = 2.0,
) -> EngineComparisonReport:
    """Run one record through both engines and compare them.

    ``docs/FORTRAN_ENGINE_DESIGN.md`` section 5-- this is the feature the
    Fortran bridge exists for: a per-run comparison against USGS peakfq 8.1.0
    on the caller's own data, not just the four sites already committed as
    parity goldens.

    **Requires the built f2py extension, and raises the same actionable
    ``ImportError`` :mod:`flowfreq.peakfqr` raises when it is absent -- there
    is no golden-file fallback.** Decided, not left open (design doc section
    9, question 4): a fallback would only ever cover the four sites already
    proven to agree, which is exactly where this comparison proves the least;
    on a caller's own record there is no golden to fall back to, so the
    feature would silently work on some inputs and not others.

    Parameters
    ----------
    peak_flows, water_years, historical_peaks, perception_thresholds,
    user_low_outlier_threshold, ema_params : see :class:`~flowfreq.bulletin17c.Bulletin17C`.
    regional_skew, regional_skew_se : float
        As in :func:`run_ffa`. Defaults to the nationwide B17C generalized
        skew, so a caller who only wants a quick parity check does not have
        to look one up.
    aeps : array-like, optional
        Annual exceedance probabilities to compare at. Defaults to
        :attr:`~flowfreq.bulletin17c.FloodFrequencyAnalysis.STANDARD_AEP`.
        Both engines are evaluated at exactly this list -- the native side's
        own quantiles/confidence limits are recomputed at *aeps* even when it
        matches the default, so the two are always compared like for like
        rather than at whatever AEPs the native fit happened to be run at.
    site_name : str, optional
        Carried through to :meth:`EngineComparisonReport.to_markdown`'s title.
    tolerance_pct, parameter_tolerance_pct, ci_tolerance_pct : float
        As in :meth:`~flowfreq.bulletin17c.Bulletin17C.validate`.

    Returns
    -------
    EngineComparisonReport

    Raises
    ------
    ImportError
        The f2py extension is not built; run
        ``python build_fortran/build.py`` (needs gfortran and meson).
    """
    import flowfreq.peakfqr  # noqa: F401 -- raise before doing any native work if absent

    from .bulletin17c import FloodFrequencyAnalysis
    from .fortran_engine import run_fortran_reference

    if aeps is None:
        aeps = FloodFrequencyAnalysis.STANDARD_AEP
    aeps = np.asarray(aeps, dtype=float)

    native = Bulletin17C(
        peak_flows=peak_flows,
        water_years=water_years,
        regional_skew=regional_skew,
        regional_skew_mse=regional_skew_se**2,
        historical_peaks=historical_peaks,
        perception_thresholds=perception_thresholds,
        ema_params=ema_params,
        user_low_outlier_threshold=user_low_outlier_threshold,
    )
    native.run_analysis(method="ema", engine="native")
    # Recompute at *aeps* explicitly: run_analysis() always fits at
    # STANDARD_AEP internally, so a caller-supplied aeps list would otherwise
    # leave native.results.quantiles at a different AEP set than the
    # reference below, and FrequencyComparator would silently compare nothing
    # (no aep keys in common) rather than raise.
    native.results.quantiles = native.compute_quantiles(aep=aeps)
    native.results.confidence_limits = native.compute_confidence_limits(aep=aeps)

    reference, _arrays = run_fortran_reference(
        peak_flows,
        water_years=water_years,
        historical_peaks=historical_peaks,
        perception_thresholds=perception_thresholds,
        user_low_outlier_threshold=user_low_outlier_threshold,
        ema_params=ema_params,
        regional_skew=regional_skew,
        regional_skew_mse=regional_skew_se**2,
        aeps=aeps,
        station_name=site_name,
    )

    comparison = native.validate(
        reference,
        tolerance_pct=tolerance_pct,
        parameter_tolerance_pct=parameter_tolerance_pct,
        ci_tolerance_pct=ci_tolerance_pct,
    )
    return EngineComparisonReport(
        native=native.results, reference=reference, comparison=comparison, site_name=site_name
    )
