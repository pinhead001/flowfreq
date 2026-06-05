"""Tests for benchmark framework."""

from __future__ import annotations

import pytest

from hydrolib.validation.benchmarks import (
    BENCHMARKS,
    Benchmark,
    register_benchmarks,
    run_all_benchmarks,
)
from tests.peakfqsa.fixtures.big_sandy import EXPECTED_QUANTILES, SYSTEMATIC_PEAKS


class TestBenchmark:
    """Tests for Benchmark class."""

    def test_register_benchmarks(self) -> None:
        """register_benchmarks populates the registry."""
        register_benchmarks()
        assert "big_sandy" in BENCHMARKS

    def test_big_sandy_benchmark_fields(self) -> None:
        """Big Sandy benchmark has expected fields."""
        register_benchmarks()
        bm = BENCHMARKS["big_sandy"]
        assert bm.name == "big_sandy"
        assert len(bm.peaks) > 40
        assert len(bm.historical) > 0
        assert len(bm.expected_quantiles) > 0

    def test_simple_benchmark_run_native(self) -> None:
        """run_native succeeds on a simple (systematic-only) benchmark."""
        bm = Benchmark(
            name="simple_test",
            peaks=dict(SYSTEMATIC_PEAKS),
            regional_skew=-0.5,
            regional_skew_sd=0.55,
        )
        result = bm.run_native()
        assert "parameters" in result
        assert "quantiles" in result
        assert "mean_log" in result["parameters"]

    def test_simple_benchmark_validate(self) -> None:
        """validate_against_expected returns a ComparisonResult for simple case."""
        # Use a small subset of expected quantiles with generous tolerance
        bm = Benchmark(
            name="simple_test",
            peaks=dict(SYSTEMATIC_PEAKS),
            expected_parameters={"mean_log": 3.717272, "std_log": 0.289200},
            expected_quantiles={0.01: 23158.65},
            regional_skew=-0.5,
            regional_skew_sd=0.55,
            tolerance_pct=50.0,  # Generous since systematic-only vs full analysis
        )
        result = bm.validate_against_expected()
        assert result.summary is not None
        assert result.max_diff_pct >= 0

    def test_run_all_benchmarks(self) -> None:
        """run_all_benchmarks processes all registered benchmarks."""
        results = run_all_benchmarks()
        assert "big_sandy" in results
        # Big Sandy with historical data may not converge, but should not crash
        assert results["big_sandy"].summary is not None
