"""
Comparison engine for native vs reference frequency analysis results.

Compares HydroLib native EMA output against a
:class:`~hydrolib.validation.reference.ReferenceResult` with configurable
tolerance thresholds per output category.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

from hydrolib.validation.reference import ReferenceResult

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
    skew_diffs : dict[str, float]
        Skew parameter name to **absolute** difference, in skew units.
        Skew is the one LP3 parameter that legitimately passes through zero --
        Big Sandy's at-site skew is 0.0066 -- so a percent difference on it
        measures the denominator, not the disagreement: a gap of 0.016 reads as
        249%. These are judged against ``skew_tolerance_abs`` instead, and are
        excluded from ``parameter_diffs`` and from ``max_diff_pct``.
    max_diff_pct : float
        Maximum percent difference across all comparisons. Skews are not in it,
        because they are not measured in percent.
    summary : str
        Human-readable one-line summary of comparison.
    """

    passed: bool = False
    tolerance_pct: float = 1.0
    parameter_diffs: dict[str, float] = field(default_factory=dict)
    quantile_diffs: dict[float, float] = field(default_factory=dict)
    ci_diffs: dict[float, float] = field(default_factory=dict)
    skew_diffs: dict[str, float] = field(default_factory=dict)
    max_diff_pct: float = 0.0
    summary: str = ""


def _is_skew(parameter: str) -> bool:
    """Whether a parameter name denotes a skew coefficient.

    Skews are compared in their own units rather than by percent; see
    :meth:`FrequencyComparator.compare_skews`.
    """
    return "skew" in parameter.lower()


def _pct_diff(native_val: float, ref_val: float, scale_floor: float = 0.0) -> float:
    """Compute percent difference between two values.

    Parameters
    ----------
    native_val : float
        Value from the native analysis.
    ref_val : float
        Reference value.
    scale_floor : float, optional
        Lower bound on the denominator. Percent difference is meaningless for a
        quantity that legitimately passes through zero: Big Sandy's at-site
        skew is 0.0066 under peakfq 8.1.0, so an absolute gap of 0.016 -- small,
        in skew units -- divides out to 249%, which then dominates
        ``max_diff_pct`` and hides everything else in the report. With a floor
        of 0.1 the same gap reads 16%. Left at 0.0 the behaviour is unchanged,
        which is what discharges want: they are in the thousands and never
        approach any sane floor.

    Returns
    -------
    float
        Absolute percent difference. Returns 0.0 if both values are zero.
    """
    denominator = max(abs(ref_val), scale_floor)
    if denominator == 0.0:
        return 0.0 if native_val == 0.0 else 100.0
    return abs(native_val - ref_val) / denominator * 100.0


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
    parameter_scale_floor : float
        Denominator floor for parameter percent differences; see
        :func:`_pct_diff`. Applies to non-skew parameters only, not to
        quantiles or confidence intervals.
    skew_tolerance_abs : float
        Tolerance for skew parameters, as an **absolute** difference in skew
        units. Skew is compared this way rather than by percent because it
        legitimately passes through zero; see :meth:`compare_skews`.
    """

    def __init__(
        self,
        tolerance_pct: float = 1.0,
        parameter_tolerance_pct: float = 0.5,
        ci_tolerance_pct: float = 2.0,
        parameter_scale_floor: float = 0.1,
        skew_tolerance_abs: float = 0.05,
    ) -> None:
        self.tolerance_pct = tolerance_pct
        self.parameter_tolerance_pct = parameter_tolerance_pct
        self.ci_tolerance_pct = ci_tolerance_pct
        self.parameter_scale_floor = parameter_scale_floor
        self.skew_tolerance_abs = skew_tolerance_abs

    def compare(
        self,
        native: dict[str, Any],
        reference: ReferenceResult,
    ) -> ComparisonResult:
        """Compare native analysis output against a reference result.

        Parameters
        ----------
        native : dict
            Output from HydroLib native analysis. Expected keys:
            'parameters', 'quantiles', 'confidence_intervals'.
        reference : ReferenceResult
            Reference result to compare against.

        Returns
        -------
        ComparisonResult
            Detailed comparison with per-field differences.
        """
        param_diffs = self.compare_parameters(native, reference)
        skew_diffs = self.compare_skews(native, reference)
        quant_diffs = self.compare_quantiles(native, reference)
        ci_diffs = self.compare_ci(native, reference)

        # Skews are deliberately absent here: they are absolute differences in
        # skew units, so mixing them into a percent maximum would compare
        # unlike quantities and let a near-zero denominator dominate.
        all_diffs: list[float] = []
        all_diffs.extend(param_diffs.values())
        all_diffs.extend(quant_diffs.values())
        all_diffs.extend(ci_diffs.values())
        max_diff = max(all_diffs) if all_diffs else 0.0

        # Check pass/fail per category
        params_ok = all(d <= self.parameter_tolerance_pct for d in param_diffs.values())
        skews_ok = all(d <= self.skew_tolerance_abs for d in skew_diffs.values())
        quants_ok = all(d <= self.tolerance_pct for d in quant_diffs.values())
        cis_ok = all(d <= self.ci_tolerance_pct for d in ci_diffs.values())
        passed = params_ok and skews_ok and quants_ok and cis_ok

        # Build summary
        n_params = len(param_diffs)
        n_quants = len(quant_diffs)
        n_cis = len(ci_diffs)
        status = "PASS" if passed else "FAIL"
        worst_skew = max(skew_diffs.values(), default=0.0)
        skew_note = f", skews={len(skew_diffs)} (worst {worst_skew:.4f})" if skew_diffs else ""
        summary = (
            f"{status}: max diff {max_diff:.3f}% "
            f"(params={n_params}, quantiles={n_quants}, CIs={n_cis}{skew_note})"
        )

        return ComparisonResult(
            passed=passed,
            tolerance_pct=self.tolerance_pct,
            parameter_diffs=param_diffs,
            quantile_diffs=quant_diffs,
            ci_diffs=ci_diffs,
            skew_diffs=skew_diffs,
            max_diff_pct=max_diff,
            summary=summary,
        )

    def compare_parameters(self, native: dict[str, Any], ref: ReferenceResult) -> dict[str, float]:
        """Compare LP3 parameters between native and reference.

        Parameters
        ----------
        native : dict
            Native analysis output with 'parameters' key.
        ref : ReferenceResult
            Reference result.

        Returns
        -------
        dict[str, float]
            Parameter name to percent difference.
        """
        diffs: dict[str, float] = {}
        native_params = native.get("parameters", {})

        for key, ref_val in ref.parameters.items():
            if _is_skew(key):
                continue  # compared in skew units by compare_skews
            if key in native_params:
                diffs[key] = _pct_diff(
                    native_params[key], ref_val, scale_floor=self.parameter_scale_floor
                )
            else:
                logger.debug("Parameter '%s' not in native output, skipping", key)

        return diffs

    def compare_skews(self, native: dict[str, Any], ref: ReferenceResult) -> dict[str, float]:
        """Compare skew parameters by absolute difference, in skew units.

        Percent difference is the wrong metric for a quantity that passes
        through zero. Big Sandy's reference at-site skew is 0.0066, so an
        absolute gap of 0.016 -- small, in skew units -- divides out to 249%
        and swamps every other number in the report.

        Parameters
        ----------
        native : dict
            Native analysis output with a 'parameters' key.
        ref : ReferenceResult
            Reference result.

        Returns
        -------
        dict[str, float]
            Skew parameter name to absolute difference.
        """
        diffs: dict[str, float] = {}
        native_params = native.get("parameters", {})

        for key, ref_val in ref.parameters.items():
            if not _is_skew(key):
                continue
            if key in native_params:
                diffs[key] = abs(float(native_params[key]) - float(ref_val))
            else:
                logger.debug("Skew parameter '%s' not in native output, skipping", key)

        return diffs

    def compare_quantiles(self, native: dict[str, Any], ref: ReferenceResult) -> dict[float, float]:
        """Compare quantile estimates between native and reference.

        Parameters
        ----------
        native : dict
            Native analysis output with 'quantiles' key.
        ref : ReferenceResult
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

    def compare_ci(self, native: dict[str, Any], ref: ReferenceResult) -> dict[float, float]:
        """Compare confidence intervals between native and reference.

        The percent difference is the maximum of lower and upper bound
        differences for each AEP.

        Parameters
        ----------
        native : dict
            Native analysis output with 'confidence_intervals' key.
        ref : ReferenceResult
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
