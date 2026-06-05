"""
Comparison engine for native vs reference frequency analysis results.

Compares HydroLib native EMA output against PeakfqSA/reference results
with configurable tolerance thresholds per output category.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

from hydrolib.peakfqsa.parsers import PeakfqSAResult

logger = logging.getLogger(__name__)


@dataclass
class ComparisonResult:
    """Result of comparing native and reference frequency analyses.

    Parameters
    ----------
    passed : bool
        True if all differences are within tolerance.
    tolerance_pct : float
        Tolerance threshold used for the comparison.
    parameter_diffs : dict[str, float]
        Parameter name to percent difference mapping.
    quantile_diffs : dict[float, float]
        AEP to percent difference mapping for quantiles.
    ci_diffs : dict[float, float]
        AEP to percent difference mapping for CI bounds.
    max_diff_pct : float
        Maximum percent difference across all comparisons.
    summary : str
        Human-readable one-line summary of comparison.
    """

    passed: bool = False
    tolerance_pct: float = 1.0
    parameter_diffs: dict[str, float] = field(default_factory=dict)
    quantile_diffs: dict[float, float] = field(default_factory=dict)
    ci_diffs: dict[float, float] = field(default_factory=dict)
    max_diff_pct: float = 0.0
    summary: str = ""


def _pct_diff(native_val: float, ref_val: float) -> float:
    """Compute percent difference between two values.

    Parameters
    ----------
    native_val : float
        Value from the native analysis.
    ref_val : float
        Reference value.

    Returns
    -------
    float
        Absolute percent difference. Returns 0.0 if both values are zero.
    """
    if ref_val == 0.0:
        if native_val == 0.0:
            return 0.0
        return 100.0
    return abs((native_val - ref_val) / ref_val) * 100.0


class FrequencyComparator:
    """Compare native HydroLib analysis against a reference result.

    Parameters
    ----------
    tolerance_pct : float
        Default tolerance for quantile comparisons.
    parameter_tolerance_pct : float
        Tolerance for LP3 parameter comparisons.
    ci_tolerance_pct : float
        Tolerance for confidence interval comparisons.
    """

    def __init__(
        self,
        tolerance_pct: float = 1.0,
        parameter_tolerance_pct: float = 0.5,
        ci_tolerance_pct: float = 2.0,
    ) -> None:
        self.tolerance_pct = tolerance_pct
        self.parameter_tolerance_pct = parameter_tolerance_pct
        self.ci_tolerance_pct = ci_tolerance_pct

    def compare(
        self,
        native: dict[str, Any],
        reference: PeakfqSAResult,
    ) -> ComparisonResult:
        """Compare native analysis output against a reference result.

        Parameters
        ----------
        native : dict
            Output from HydroLib native analysis. Expected keys:
            'parameters', 'quantiles', 'confidence_intervals'.
        reference : PeakfqSAResult
            Reference result to compare against.

        Returns
        -------
        ComparisonResult
            Detailed comparison with per-field differences.
        """
        param_diffs = self.compare_parameters(native, reference)
        quant_diffs = self.compare_quantiles(native, reference)
        ci_diffs = self.compare_ci(native, reference)

        # Determine max difference across all categories
        all_diffs: list[float] = []
        all_diffs.extend(param_diffs.values())
        all_diffs.extend(quant_diffs.values())
        all_diffs.extend(ci_diffs.values())
        max_diff = max(all_diffs) if all_diffs else 0.0

        # Check pass/fail per category
        params_ok = all(d <= self.parameter_tolerance_pct for d in param_diffs.values())
        quants_ok = all(d <= self.tolerance_pct for d in quant_diffs.values())
        cis_ok = all(d <= self.ci_tolerance_pct for d in ci_diffs.values())
        passed = params_ok and quants_ok and cis_ok

        # Build summary
        n_params = len(param_diffs)
        n_quants = len(quant_diffs)
        n_cis = len(ci_diffs)
        status = "PASS" if passed else "FAIL"
        summary = (
            f"{status}: max diff {max_diff:.3f}% "
            f"(params={n_params}, quantiles={n_quants}, CIs={n_cis})"
        )

        return ComparisonResult(
            passed=passed,
            tolerance_pct=self.tolerance_pct,
            parameter_diffs=param_diffs,
            quantile_diffs=quant_diffs,
            ci_diffs=ci_diffs,
            max_diff_pct=max_diff,
            summary=summary,
        )

    def compare_parameters(self, native: dict[str, Any], ref: PeakfqSAResult) -> dict[str, float]:
        """Compare LP3 parameters between native and reference.

        Parameters
        ----------
        native : dict
            Native analysis output with 'parameters' key.
        ref : PeakfqSAResult
            Reference result.

        Returns
        -------
        dict[str, float]
            Parameter name to percent difference.
        """
        diffs: dict[str, float] = {}
        native_params = native.get("parameters", {})

        for key, ref_val in ref.parameters.items():
            if key in native_params:
                diffs[key] = _pct_diff(native_params[key], ref_val)
            else:
                logger.debug("Parameter '%s' not in native output, skipping", key)

        return diffs

    def compare_quantiles(self, native: dict[str, Any], ref: PeakfqSAResult) -> dict[float, float]:
        """Compare quantile estimates between native and reference.

        Parameters
        ----------
        native : dict
            Native analysis output with 'quantiles' key.
        ref : PeakfqSAResult
            Reference result.

        Returns
        -------
        dict[float, float]
            AEP to percent difference.
        """
        diffs: dict[float, float] = {}
        native_quantiles = native.get("quantiles", {})

        for aep, ref_val in ref.quantiles.items():
            if aep in native_quantiles:
                diffs[aep] = _pct_diff(native_quantiles[aep], ref_val)
            else:
                logger.debug("AEP %s not in native quantiles, skipping", aep)

        return diffs

    def compare_ci(self, native: dict[str, Any], ref: PeakfqSAResult) -> dict[float, float]:
        """Compare confidence intervals between native and reference.

        The percent difference is the maximum of lower and upper bound
        differences for each AEP.

        Parameters
        ----------
        native : dict
            Native analysis output with 'confidence_intervals' key.
        ref : PeakfqSAResult
            Reference result.

        Returns
        -------
        dict[float, float]
            AEP to percent difference on bounds.
        """
        diffs: dict[float, float] = {}
        native_cis = native.get("confidence_intervals", {})

        for aep, (ref_lo, ref_hi) in ref.confidence_intervals.items():
            if aep in native_cis:
                nat_lo, nat_hi = native_cis[aep]
                diff_lo = _pct_diff(nat_lo, ref_lo)
                diff_hi = _pct_diff(nat_hi, ref_hi)
                diffs[aep] = max(diff_lo, diff_hi)
            else:
                logger.debug("AEP %s not in native CIs, skipping", aep)

        return diffs
