"""Tests for benchmark framework."""

from __future__ import annotations

import pathlib

import pytest

import hydrolib.validation
from hydrolib.validation.benchmarks import (
    BENCHMARKS,
    DATA_DIR,
    Benchmark,
    load_case,
    register_benchmarks,
    run_all_benchmarks,
)
from tests.fixtures.big_sandy import SYSTEMATIC_PEAKS


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
        assert results["big_sandy"].summary is not None

    def test_big_sandy_passes_against_peakfq_810(self) -> None:
        """The benchmark is only worth shipping if it can pass.

        It could not before: it validated against the 2012 PeakfqSA manual,
        which peakfq 8.1.0 does not reproduce, so a correct implementation was
        guaranteed a FAIL. It now compares against the 8.1.0 reference.
        """
        assert run_all_benchmarks()["big_sandy"].passed

    def test_every_registered_benchmark_can_be_fitted(self) -> None:
        """Two entries used to be registered with no peaks at all.

        fortran_respec and skew_weighting are routine-level checks of moms_p3,
        p3est_ema and detrat, not frequency analyses, so every run of
        `hydrolib benchmark` ended in "zero-size array to reduction operation
        minimum". They are covered by tests/test_r_fixtures.py instead.
        """
        register_benchmarks()
        assert BENCHMARKS
        for name, bm in BENCHMARKS.items():
            assert bm.peaks, f"benchmark {name} has no peaks to fit"

    def test_run_native_rejects_a_benchmark_with_no_peaks(self) -> None:
        with pytest.raises(ValueError, match="no systematic peaks"):
            Benchmark(name="empty").run_native()

    def test_known_deviations_are_excluded_from_the_comparison(self) -> None:
        """The open HWN skew defect is reported, not compared and not hidden."""
        register_benchmarks()
        bm = BENCHMARKS["big_sandy"]
        assert "skew_weighted" in bm.known_deviations
        assert "TODO.md" in bm.known_deviations["skew_weighted"]
        compared = bm.validate_against_expected().parameter_diffs
        assert "skew_weighted" not in compared
        assert "mean_log" in compared


class TestPackagedCases:
    """Benchmark data ships with the package; tests/ does not."""

    def test_load_case_reads_packaged_data(self) -> None:
        case = load_case("big_sandy_03606500")
        assert case["name"] == "big_sandy"
        assert len(case["systematic_peaks"]) == 44

    def test_unknown_case_names_what_is_available(self) -> None:
        with pytest.raises(FileNotFoundError, match="big_sandy_03606500"):
            load_case("no_such_site")

    def test_case_data_is_not_reachable_only_from_the_repo_root(self) -> None:
        """The defect this replaced: benchmarks.py imported from tests/.

        tests/ is not in [tool.setuptools.packages.find], so `hydrolib
        benchmark` raised ModuleNotFoundError for every installed user. The
        data directory has to resolve relative to the package itself.
        """
        assert DATA_DIR.is_dir()
        assert DATA_DIR.is_relative_to(pathlib.Path(hydrolib.validation.__file__).parent)
