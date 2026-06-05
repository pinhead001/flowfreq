"""
Bulletin 17C benchmark test cases.

Each benchmark is a named, self-contained test scenario with known inputs
and expected outputs from official USGS examples. Used for validating
HydroLib native EMA against published reference values.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np

from hydrolib.peakfqsa.parsers import PeakfqSAResult
from hydrolib.validation.comparisons import ComparisonResult, FrequencyComparator

logger = logging.getLogger(__name__)


@dataclass
class Benchmark:
    """A single benchmark test case.

    Parameters
    ----------
    name : str
        Short identifier for the benchmark.
    description : str
        Human-readable description.
    peaks : dict[int, float]
        Systematic peak flows {water_year: discharge}.
    historical : dict[int, float]
        Historical peaks {water_year: discharge}.
    thresholds : list[dict]
        Perception thresholds.
    expected_parameters : dict[str, float]
        Expected LP3 parameters from reference.
    expected_quantiles : dict[float, float]
        Expected AEP to discharge mapping.
    expected_confidence_intervals : dict[float, tuple[float, float]]
        Expected AEP to (lower, upper) CI mapping.
    tolerance_pct : float
        Tolerance for passing.
    regional_skew : float
        Regional skew coefficient.
    regional_skew_sd : float
        Standard deviation of regional skew.
    begyear : int
        Analysis start year.
    endyear : int
        Analysis end year.
    """

    name: str = ""
    description: str = ""
    peaks: dict[int, float] = field(default_factory=dict)
    historical: dict[int, float] = field(default_factory=dict)
    thresholds: list[dict[str, Any]] = field(default_factory=list)
    expected_parameters: dict[str, float] = field(default_factory=dict)
    expected_quantiles: dict[float, float] = field(default_factory=dict)
    expected_confidence_intervals: dict[float, tuple[float, float]] = field(default_factory=dict)
    tolerance_pct: float = 1.0
    regional_skew: float = -0.302
    regional_skew_sd: float = 0.55
    begyear: int = 0
    endyear: int = 0

    def run_native(self) -> dict[str, Any]:
        """Run benchmark using HydroLib native EMA.

        Returns
        -------
        dict
            Native analysis result in comparison dict format.
        """
        from hydrolib.bulletin17c import Bulletin17C

        peak_flows = np.array(list(self.peaks.values()))
        water_years = np.array(list(self.peaks.keys()))

        historical = [(year, q) for year, q in self.historical.items()] or None

        b17c = Bulletin17C(
            peak_flows=peak_flows,
            water_years=water_years,
            regional_skew=self.regional_skew,
            regional_skew_mse=self.regional_skew_sd**2,
            historical_peaks=historical,
        )
        b17c.run_analysis(method="ema")
        return b17c.to_comparison_dict()

    def validate_against_expected(self) -> ComparisonResult:
        """Validate native results against expected values.

        Returns
        -------
        ComparisonResult
            Comparison between native and expected results.
        """
        native = self.run_native()

        reference = PeakfqSAResult(
            parameters=dict(self.expected_parameters),
            quantiles=dict(self.expected_quantiles),
            confidence_intervals={
                aep: (lo, hi) for aep, (lo, hi) in self.expected_confidence_intervals.items()
            },
        )

        comparator = FrequencyComparator(
            tolerance_pct=self.tolerance_pct,
            parameter_tolerance_pct=self.tolerance_pct,
            ci_tolerance_pct=self.tolerance_pct * 2,
        )
        return comparator.compare(native, reference)


def _create_big_sandy_benchmark() -> Benchmark:
    """Create the Big Sandy River benchmark from fixture data."""
    from tests.peakfqsa.fixtures.big_sandy import (
        BEGYEAR,
        ENDYEAR,
        EXPECTED_CONFIDENCE_INTERVALS,
        EXPECTED_PARAMETERS,
        EXPECTED_QUANTILES,
        HISTORICAL_PEAKS,
        REGIONAL_SKEW,
        REGIONAL_SKEW_SD,
        SYSTEMATIC_PEAKS,
        THRESHOLDS,
        TOLERANCE_PERCENT,
    )

    return Benchmark(
        name="big_sandy",
        description="Big Sandy River at Bruceton, TN (USGS 03606500) - PeakfqSA manual example",
        peaks=dict(SYSTEMATIC_PEAKS),
        historical=dict(HISTORICAL_PEAKS),
        thresholds=[dict(t) for t in THRESHOLDS],
        expected_parameters=dict(EXPECTED_PARAMETERS),
        expected_quantiles=dict(EXPECTED_QUANTILES),
        expected_confidence_intervals=dict(EXPECTED_CONFIDENCE_INTERVALS),
        tolerance_pct=TOLERANCE_PERCENT,
        regional_skew=REGIONAL_SKEW,
        regional_skew_sd=REGIONAL_SKEW_SD,
        begyear=BEGYEAR,
        endyear=ENDYEAR,
    )


# Registry of available benchmarks
BENCHMARKS: dict[str, Benchmark] = {}


def _create_fortran_respec_benchmark() -> Benchmark:
    """Create benchmark from Respec Fortran test fixture (moms_p3/p3est_ema)."""
    from tests.peakfqsa.fixtures.fortran_respec import (
        MOMS_P3_EXPECTED,
        P3EST_EMA_EXPECTED,
        TRUTH_MOMS_P3_EXPECTED,
    )

    return Benchmark(
        name="fortran_respec",
        description=(
            "Respec EMA Fortran tests — moms_p3, p3est_ema, var_mom, qP3 "
            "from peakfqr test-fortran.R"
        ),
        expected_parameters={
            "moms_p3_orig_mean": MOMS_P3_EXPECTED["orig"][0],
            "moms_p3_orig_var": MOMS_P3_EXPECTED["orig"][1],
            "moms_p3_orig_skew": MOMS_P3_EXPECTED["orig"][2],
            "p3est_ema_orig_mean": P3EST_EMA_EXPECTED["orig"][0],
            "p3est_ema_orig_var": P3EST_EMA_EXPECTED["orig"][1],
            "p3est_ema_orig_skew": P3EST_EMA_EXPECTED["orig"][2],
            "truth_mean": TRUTH_MOMS_P3_EXPECTED[0],
            "truth_var": TRUTH_MOMS_P3_EXPECTED[1],
            "truth_skew": TRUTH_MOMS_P3_EXPECTED[2],
        },
        tolerance_pct=0.01,
    )


def _create_skew_weighting_benchmark() -> Benchmark:
    """Create benchmark from Halloween skew weighting tests."""
    return Benchmark(
        name="skew_weighting",
        description=(
            "Halloween determinant ratio skew weighting — "
            "detrat() vs Greg Schwarz SAS output (312 cases)"
        ),
        tolerance_pct=0.1,
    )


def register_benchmarks() -> None:
    """Populate the BENCHMARKS registry with all available benchmarks."""
    if "big_sandy" not in BENCHMARKS:
        BENCHMARKS["big_sandy"] = _create_big_sandy_benchmark()
    if "fortran_respec" not in BENCHMARKS:
        BENCHMARKS["fortran_respec"] = _create_fortran_respec_benchmark()
    if "skew_weighting" not in BENCHMARKS:
        BENCHMARKS["skew_weighting"] = _create_skew_weighting_benchmark()


def run_all_benchmarks() -> dict[str, ComparisonResult]:
    """Run all registered benchmarks against expected values.

    Returns
    -------
    dict[str, ComparisonResult]
        Benchmark name to comparison result mapping.
    """
    register_benchmarks()
    results: dict[str, ComparisonResult] = {}

    for name, benchmark in BENCHMARKS.items():
        logger.info("Running benchmark: %s", name)
        try:
            results[name] = benchmark.validate_against_expected()
        except Exception as e:
            logger.error("Benchmark '%s' failed: %s", name, e)
            results[name] = ComparisonResult(
                passed=False,
                summary=f"ERROR: {e}",
            )

    return results


def print_benchmark_report(results: dict[str, ComparisonResult]) -> None:
    """Print a formatted report of benchmark results.

    Parameters
    ----------
    results : dict[str, ComparisonResult]
        Results from run_all_benchmarks.
    """
    print("\n" + "=" * 60)
    print("  HydroLib Benchmark Report")
    print("=" * 60)

    n_pass = sum(1 for r in results.values() if r.passed)
    n_total = len(results)

    for name, result in results.items():
        status = "PASS" if result.passed else "FAIL"
        print(f"\n  [{status}] {name}")
        print(f"         Max diff: {result.max_diff_pct:.3f}%")
        print(f"         {result.summary}")

    print("\n" + "-" * 60)
    print(f"  Total: {n_pass}/{n_total} passed")
    print("=" * 60 + "\n")
