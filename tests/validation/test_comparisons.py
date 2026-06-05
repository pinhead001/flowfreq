"""Tests for the comparison engine."""

from __future__ import annotations

import pytest

from hydrolib.peakfqsa.parsers import PeakfqSAResult
from hydrolib.validation.comparisons import ComparisonResult, FrequencyComparator, _pct_diff


class TestPctDiff:
    """Tests for percent difference helper."""

    def test_identical_values(self) -> None:
        """Identical values give 0% difference."""
        assert _pct_diff(100.0, 100.0) == 0.0

    def test_both_zero(self) -> None:
        """Both zero gives 0% difference."""
        assert _pct_diff(0.0, 0.0) == 0.0

    def test_ref_zero(self) -> None:
        """Reference zero with nonzero native gives 100%."""
        assert _pct_diff(5.0, 0.0) == 100.0

    def test_ten_percent_diff(self) -> None:
        """10% difference."""
        assert abs(_pct_diff(110.0, 100.0) - 10.0) < 0.01

    def test_negative_values(self) -> None:
        """Handles negative values correctly."""
        assert abs(_pct_diff(-0.12, -0.10) - 20.0) < 0.01


class TestComparisonResult:
    """Tests for ComparisonResult dataclass."""

    def test_default_is_fail(self) -> None:
        """Default result is a failure."""
        cr = ComparisonResult()
        assert cr.passed is False

    def test_fields(self) -> None:
        """All fields are accessible."""
        cr = ComparisonResult(
            passed=True,
            tolerance_pct=1.0,
            max_diff_pct=0.5,
            summary="PASS",
        )
        assert cr.passed is True
        assert cr.max_diff_pct == 0.5


class TestFrequencyComparator:
    """Tests for FrequencyComparator."""

    def _make_reference(self) -> PeakfqSAResult:
        """Create a reference result for testing."""
        return PeakfqSAResult(
            parameters={
                "mean_log": 3.717272,
                "std_log": 0.289200,
                "skew_weighted": -0.118702,
            },
            quantiles={
                0.01: 23158.65,
                0.02: 19617.73,
                0.10: 12134.65,
                0.50: 5284.36,
            },
            confidence_intervals={
                0.01: (17388.03, 37986.08),
                0.02: (15154.99, 29124.18),
                0.10: (9766.00, 15218.32),
            },
        )

    def test_identical_results_pass(self) -> None:
        """Identical native and reference results pass comparison."""
        ref = self._make_reference()
        native = {
            "parameters": dict(ref.parameters),
            "quantiles": dict(ref.quantiles),
            "confidence_intervals": dict(ref.confidence_intervals),
        }

        comp = FrequencyComparator()
        result = comp.compare(native, ref)

        assert result.passed is True
        assert result.max_diff_pct == 0.0
        assert "PASS" in result.summary

    def test_within_tolerance_passes(self) -> None:
        """Results within tolerance pass."""
        ref = self._make_reference()
        native = {
            "parameters": {
                "mean_log": 3.717272 * 1.003,  # 0.3% off
                "std_log": 0.289200 * 1.004,  # 0.4% off
                "skew_weighted": -0.118702 * 1.004,
            },
            "quantiles": {
                0.01: 23158.65 * 1.005,  # 0.5% off
                0.02: 19617.73 * 1.008,  # 0.8% off
                0.10: 12134.65 * 1.003,
                0.50: 5284.36 * 1.002,
            },
            "confidence_intervals": {
                0.01: (17388.03 * 1.01, 37986.08 * 1.015),  # 1.5% off
                0.02: (15154.99 * 1.01, 29124.18 * 1.01),
                0.10: (9766.00 * 1.005, 15218.32 * 1.005),
            },
        }

        comp = FrequencyComparator(
            tolerance_pct=1.0, parameter_tolerance_pct=0.5, ci_tolerance_pct=2.0
        )
        result = comp.compare(native, ref)

        assert result.passed is True

    def test_parameter_exceeds_tolerance_fails(self) -> None:
        """Parameter exceeding tolerance causes failure."""
        ref = self._make_reference()
        native = {
            "parameters": {
                "mean_log": 3.717272 * 1.01,  # 1% off, exceeds 0.5% tolerance
                "std_log": 0.289200,
                "skew_weighted": -0.118702,
            },
            "quantiles": dict(ref.quantiles),
            "confidence_intervals": dict(ref.confidence_intervals),
        }

        comp = FrequencyComparator(parameter_tolerance_pct=0.5)
        result = comp.compare(native, ref)

        assert result.passed is False
        assert "FAIL" in result.summary

    def test_quantile_exceeds_tolerance_fails(self) -> None:
        """Quantile exceeding tolerance causes failure."""
        ref = self._make_reference()
        native = {
            "parameters": dict(ref.parameters),
            "quantiles": {
                0.01: 23158.65 * 1.02,  # 2% off, exceeds 1% tolerance
                0.02: 19617.73,
                0.10: 12134.65,
                0.50: 5284.36,
            },
            "confidence_intervals": dict(ref.confidence_intervals),
        }

        comp = FrequencyComparator(tolerance_pct=1.0)
        result = comp.compare(native, ref)

        assert result.passed is False
        assert result.quantile_diffs[0.01] > 1.0

    def test_ci_exceeds_tolerance_fails(self) -> None:
        """CI exceeding tolerance causes failure."""
        ref = self._make_reference()
        native = {
            "parameters": dict(ref.parameters),
            "quantiles": dict(ref.quantiles),
            "confidence_intervals": {
                0.01: (17388.03, 37986.08 * 1.05),  # 5% off, exceeds 2%
                0.02: (15154.99, 29124.18),
                0.10: (9766.00, 15218.32),
            },
        }

        comp = FrequencyComparator(ci_tolerance_pct=2.0)
        result = comp.compare(native, ref)

        assert result.passed is False
        assert result.ci_diffs[0.01] > 2.0

    def test_empty_native_passes_vacuously(self) -> None:
        """Empty native result passes (no comparisons to fail)."""
        ref = self._make_reference()
        native: dict = {"parameters": {}, "quantiles": {}, "confidence_intervals": {}}

        comp = FrequencyComparator()
        result = comp.compare(native, ref)

        assert result.passed is True
        assert result.max_diff_pct == 0.0

    def test_compare_parameters_only(self) -> None:
        """compare_parameters works independently."""
        ref = self._make_reference()
        native = {"parameters": {"mean_log": 3.72, "std_log": 0.29}}

        comp = FrequencyComparator()
        diffs = comp.compare_parameters(native, ref)

        assert "mean_log" in diffs
        assert "std_log" in diffs
        assert diffs["mean_log"] > 0

    def test_compare_quantiles_only(self) -> None:
        """compare_quantiles works independently."""
        ref = self._make_reference()
        native = {"quantiles": {0.01: 23200.0, 0.10: 12000.0}}

        comp = FrequencyComparator()
        diffs = comp.compare_quantiles(native, ref)

        assert 0.01 in diffs
        assert 0.10 in diffs

    def test_max_diff_tracked(self) -> None:
        """max_diff_pct reflects the largest difference."""
        ref = self._make_reference()
        native = {
            "parameters": {"mean_log": 3.717272},
            "quantiles": {0.01: 23158.65 * 1.03},  # 3% off
            "confidence_intervals": {},
        }

        comp = FrequencyComparator()
        result = comp.compare(native, ref)

        assert result.max_diff_pct > 2.5
