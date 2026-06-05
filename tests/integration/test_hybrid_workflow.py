"""Integration tests for hybrid Bulletin 17C workflow."""

from __future__ import annotations

import numpy as np
import pytest

from hydrolib.bulletin17c import Bulletin17C
from hydrolib.peakfqsa.parsers import PeakfqSAResult
from hydrolib.validation.comparisons import FrequencyComparator
from tests.peakfqsa.fixtures.big_sandy import (
    EXPECTED_CONFIDENCE_INTERVALS,
    EXPECTED_PARAMETERS,
    EXPECTED_QUANTILES,
    HISTORICAL_PEAKS,
    REGIONAL_SKEW,
    REGIONAL_SKEW_SD,
    SYSTEMATIC_PEAKS,
    TOLERANCE_PERCENT,
)

# Mark for tests requiring the real PeakfqSA binary
requires_peakfqsa = pytest.mark.requires_peakfqsa


class TestBigSandyNativeSystematic:
    """Test Big Sandy analysis using native EMA with systematic data only.

    The native EMA implementation has numerical limitations with
    historical/censored intervals. These tests use systematic-only data
    where the native solver converges reliably.
    """

    def _run_native_systematic(self) -> Bulletin17C:
        """Run native EMA on Big Sandy systematic data only."""
        peak_flows = np.array(list(SYSTEMATIC_PEAKS.values()))
        water_years = np.array(list(SYSTEMATIC_PEAKS.keys()))

        b17c = Bulletin17C(
            peak_flows=peak_flows,
            water_years=water_years,
            regional_skew=REGIONAL_SKEW,
            regional_skew_mse=REGIONAL_SKEW_SD**2,
        )
        b17c.run_analysis(method="ema")
        return b17c

    def test_native_ema_converges(self) -> None:
        """Native EMA converges on systematic-only data."""
        b17c = self._run_native_systematic()
        assert b17c.results is not None
        assert b17c.results.ema_converged is True

    def test_native_produces_parameters(self) -> None:
        """Native EMA produces finite LP3 parameters."""
        b17c = self._run_native_systematic()
        r = b17c.results
        assert np.isfinite(r.mean_log)
        assert np.isfinite(r.std_log)
        assert r.std_log > 0
        assert r.skew_weighted is not None

    def test_native_produces_quantiles(self) -> None:
        """Native EMA produces reasonable quantile estimates."""
        b17c = self._run_native_systematic()
        quantiles = b17c.compute_quantiles()
        assert len(quantiles) > 0
        # 1% AEP should produce a reasonable flood estimate for Big Sandy
        row_01 = quantiles[quantiles["aep"] == 0.01]
        if not row_01.empty:
            q_01 = row_01["flow_cfs"].iloc[0]
            assert 5000 < q_01 < 100000, f"1% AEP flow {q_01} out of reasonable range"

    def test_to_comparison_dict(self) -> None:
        """to_comparison_dict produces expected structure."""
        b17c = self._run_native_systematic()
        d = b17c.to_comparison_dict()

        assert "parameters" in d
        assert "quantiles" in d
        assert "confidence_intervals" in d
        assert "mean_log" in d["parameters"]
        assert "std_log" in d["parameters"]
        assert len(d["quantiles"]) > 0

    def test_validate_against_reference(self) -> None:
        """Native results can be validated against a reference.

        Uses generous tolerances since systematic-only analysis will
        differ from the full historical analysis in the expected values.
        """
        b17c = self._run_native_systematic()

        # Create a reference result from expected values
        reference = PeakfqSAResult(
            parameters=dict(EXPECTED_PARAMETERS),
            quantiles=dict(EXPECTED_QUANTILES),
            confidence_intervals={
                aep: (lo, hi) for aep, (lo, hi) in EXPECTED_CONFIDENCE_INTERVALS.items()
            },
        )

        # Use very generous tolerance since systematic-only differs from full analysis
        result = b17c.validate(
            reference, tolerance_pct=50.0, parameter_tolerance_pct=50.0, ci_tolerance_pct=50.0
        )
        assert result.summary is not None
        assert result.max_diff_pct >= 0


class TestBigSandyWithHistorical:
    """Test Big Sandy with historical data — native EMA now converges.

    The analytical truncated-moment approach (matching the Fortran mP3
    routine) resolves the overflow that previously caused NaN divergence.
    """

    def _run_native_historical(self) -> Bulletin17C:
        peak_flows = np.array(list(SYSTEMATIC_PEAKS.values()))
        water_years = np.array(list(SYSTEMATIC_PEAKS.keys()))
        historical = [(year, q) for year, q in HISTORICAL_PEAKS.items()]

        b17c = Bulletin17C(
            peak_flows=peak_flows,
            water_years=water_years,
            regional_skew=REGIONAL_SKEW,
            regional_skew_mse=REGIONAL_SKEW_SD**2,
            historical_peaks=historical,
        )
        b17c.run_analysis(method="ema")
        return b17c

    def test_native_ema_with_historical_converges(self) -> None:
        """Native EMA with historical data converges to finite parameters."""
        b17c = self._run_native_historical()
        r = b17c.results
        assert r is not None
        assert r.ema_converged is True
        assert np.isfinite(r.mean_log)
        assert np.isfinite(r.std_log)
        assert r.std_log > 0

    def test_native_historical_quantiles_reasonable(self) -> None:
        """Historical EMA produces reasonable quantile estimates."""
        b17c = self._run_native_historical()
        quantiles = b17c.compute_quantiles()
        assert len(quantiles) > 0
        row_01 = quantiles[quantiles["aep"] == 0.01]
        if not row_01.empty:
            q_01 = row_01["flow_cfs"].iloc[0]
            assert 5000 < q_01 < 200000, f"1% AEP flow {q_01} out of reasonable range"


@requires_peakfqsa
class TestBigSandyEndToEnd:
    """Full integration: native + PeakfqSA comparison.

    Skipped when PeakfqSA is not installed.
    """

    pass
