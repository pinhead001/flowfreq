"""
Tolerance-based result validation for PeakfqSA outputs.

Provides functions to compare PeakfqSA results against expected values
with configurable percentage tolerances.
"""

from __future__ import annotations

import logging
from typing import Any

from hydrolib.peakfqsa.parsers import PeakfqSAResult
from hydrolib.validation.comparisons import FrequencyComparator

logger = logging.getLogger(__name__)


def validate_against_expected(
    result: PeakfqSAResult,
    expected_quantiles: dict[float, float],
    expected_parameters: dict[str, float] | None = None,
    expected_ci: dict[float, tuple[float, float]] | None = None,
    tolerance_pct: float = 1.0,
) -> dict[str, Any]:
    """Validate a PeakfqSA result against expected values.

    Parameters
    ----------
    result : PeakfqSAResult
        Result to validate.
    expected_quantiles : dict[float, float]
        Expected AEP-to-discharge mapping.
    expected_parameters : dict[str, float] or None
        Expected LP3 parameters.
    expected_ci : dict[float, tuple[float, float]] or None
        Expected confidence intervals.
    tolerance_pct : float
        Maximum allowed percent difference.

    Returns
    -------
    dict[str, Any]
        Validation report with keys: passed, max_diff_pct, details.
    """
    # Build a reference PeakfqSAResult from expected values
    reference = PeakfqSAResult(
        parameters=expected_parameters or {},
        quantiles=expected_quantiles,
        confidence_intervals=expected_ci or {},
    )

    # Build native-style dict from the result being validated
    native = {
        "parameters": dict(result.parameters),
        "quantiles": dict(result.quantiles),
        "confidence_intervals": dict(result.confidence_intervals),
    }

    comparator = FrequencyComparator(
        tolerance_pct=tolerance_pct,
        parameter_tolerance_pct=tolerance_pct,
        ci_tolerance_pct=tolerance_pct * 2,
    )
    comparison = comparator.compare(native, reference)

    return {
        "passed": comparison.passed,
        "max_diff_pct": comparison.max_diff_pct,
        "summary": comparison.summary,
        "parameter_diffs": comparison.parameter_diffs,
        "quantile_diffs": comparison.quantile_diffs,
        "ci_diffs": comparison.ci_diffs,
    }
