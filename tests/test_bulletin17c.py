"""Tests for Bulletin 17C flood frequency analysis."""

import numpy as np
import pytest

from hydrolib import (
    AnalysisMethod,
    Bulletin17C,
    ExpectedMomentsAlgorithm,
    MethodOfMoments,
    grubbs_beck_critical_value,
    kfactor,
)


# Fixtures
@pytest.fixture
def synthetic_peaks():
    """Generate synthetic peak flow data."""
    np.random.seed(42)
    n = 50
    mean_log, std_log, skew = 4.5, 0.25, 0.3

    alpha = 4 / skew**2
    z = (np.random.gamma(alpha, 1, n) - alpha) / np.sqrt(alpha)
    return 10 ** (mean_log + std_log * z)


@pytest.fixture
def water_years():
    """Water years for synthetic data."""
    return np.arange(1971, 2021)


# Core utility tests
class TestUtilities:
    def test_kfactor_zero_skew(self):
        """K-factor with zero skew should equal standard normal quantile."""
        from scipy.special import ndtri

        aep = 0.01
        assert abs(kfactor(0.0, aep) - ndtri(1 - aep)) < 0.001

    def test_kfactor_positive_skew(self):
        """K-factor with positive skew for 1% AEP."""
        K = kfactor(0.5, 0.01)
        assert 2.0 < K < 3.0  # Reasonable range

    def test_kfactor_cached(self):
        """Verify caching works (same result on repeated calls)."""
        k1 = kfactor(0.3, 0.01)
        k2 = kfactor(0.3, 0.01)
        assert k1 == k2

    def test_grubbs_beck_critical_value(self):
        """Test Grubbs-Beck critical value calculation."""
        k_50 = grubbs_beck_critical_value(50)
        assert 2.7 < k_50 < 2.8

        k_100 = grubbs_beck_critical_value(100)
        assert k_100 > k_50  # Should increase with n

    def test_kfactor_clips_extreme_skew(self):
        """K-factor should clamp implausible skew rather than blow up.

        Skew coefficients well outside the Bulletin 17B/17C Appendix 3
        table domain (|skew| <= 3) should be treated the same as the
        boundary value, not produce runaway K-factors.
        """
        k_bound = kfactor(-3.0, 0.01)
        k_extreme = kfactor(-7.0, 0.01)
        assert k_extreme == pytest.approx(k_bound)


# Method of Moments tests
class TestMethodOfMoments:
    def test_basic_analysis(self, synthetic_peaks):
        """Test basic MOM analysis runs without error."""
        mom = MethodOfMoments(synthetic_peaks)
        results = mom.run_analysis()

        assert results.n_peaks == len(synthetic_peaks)
        assert results.method == AnalysisMethod.MOM
        assert results.mean_log > 0
        assert results.std_log > 0

    def test_with_regional_skew(self, synthetic_peaks):
        """Test MOM with regional skew weighting."""
        mom = MethodOfMoments(synthetic_peaks, regional_skew=0.0, regional_skew_mse=0.15)
        results = mom.run_analysis()

        assert results.skew_weighted is not None
        assert results.skew_regional == 0.0
        # Weighted skew should be between station and regional
        assert (
            min(results.skew_station, 0.0)
            <= results.skew_weighted
            <= max(results.skew_station, 0.0)
        )

    def test_quantiles_computed(self, synthetic_peaks):
        """Test that quantiles are computed."""
        mom = MethodOfMoments(synthetic_peaks)
        results = mom.run_analysis()

        assert not results.quantiles.empty
        assert "aep" in results.quantiles.columns
        assert "flow_cfs" in results.quantiles.columns

        # 100-year flow should be greater than 10-year
        q10 = results.quantiles[results.quantiles["aep"] == 0.10]["flow_cfs"].values[0]
        q100 = results.quantiles[results.quantiles["aep"] == 0.01]["flow_cfs"].values[0]
        assert q100 > q10

    def test_confidence_limits(self, synthetic_peaks):
        """Test confidence limit computation."""
        mom = MethodOfMoments(synthetic_peaks)
        results = mom.run_analysis()

        assert not results.confidence_limits.empty
        assert "lower_5pct" in results.confidence_limits.columns
        assert "upper_5pct" in results.confidence_limits.columns

        # Upper limit should be greater than estimate
        for _, row in results.confidence_limits.iterrows():
            assert row["lower_5pct"] < row["flow_cfs"] < row["upper_5pct"]


# Expected Moments Algorithm tests
class TestEMA:
    def test_basic_ema(self, synthetic_peaks, water_years):
        """Test basic EMA analysis."""
        ema = ExpectedMomentsAlgorithm(synthetic_peaks, water_years=water_years)
        results = ema.run_analysis()

        assert results.method == AnalysisMethod.EMA
        assert results.ema_iterations is not None
        assert results.ema_converged is not None

    def test_ema_convergence(self, synthetic_peaks, water_years):
        """Test EMA converges."""
        ema = ExpectedMomentsAlgorithm(synthetic_peaks, water_years=water_years)
        results = ema.run_analysis()

        assert results.ema_converged == True
        assert results.ema_iterations < 100

    def test_ema_with_historical(self, synthetic_peaks, water_years):
        """Test EMA with historical peaks."""
        historical = [(1936, 150000), (1955, 120000)]

        ema = ExpectedMomentsAlgorithm(
            synthetic_peaks, water_years=water_years, historical_peaks=historical
        )
        results = ema.run_analysis()

        assert results.n_historical == 2
        assert results.n_peaks > len(synthetic_peaks)

    def test_ema_vs_mom_similar_without_historical(self, synthetic_peaks, water_years):
        """EMA and MOM should give similar results without historical data."""
        mom = MethodOfMoments(synthetic_peaks)
        mom_results = mom.run_analysis()

        ema = ExpectedMomentsAlgorithm(synthetic_peaks, water_years=water_years)
        ema_results = ema.run_analysis()

        # Mean and std should be close (within 5%)
        assert abs(ema_results.mean_log - mom_results.mean_log) / mom_results.mean_log < 0.05
        assert abs(ema_results.std_log - mom_results.std_log) / mom_results.std_log < 0.10

    def test_ema_low_outliers_do_not_blow_up_quantiles(self):
        """Regression test: a handful of low outliers should not cause the
        EMA moment-matching iteration to diverge to an implausible skew and
        wildly overestimate the flood quantiles.

        Reproduces a reported bug where a record with peaks in the ~20,000
        cfs range (a few low outliers censored by the Multiple Grubbs-Beck
        test) produced a Q100 estimate above 550,000 cfs -- ~27x the largest
        observed peak -- because the EMA skew fixed point ran away to an
        unphysical value (~-7) instead of stabilizing near the sample skew.
        """
        np.random.seed(1)
        years = np.arange(1990, 2021)
        peaks = np.random.lognormal(mean=np.log(20000), sigma=0.25, size=len(years))
        peaks[2] = 3000
        peaks[5] = 4500
        peaks[10] = 2000

        b17c = Bulletin17C(
            peak_flows=peaks, water_years=years, regional_skew=-0.302, regional_skew_mse=0.55**2
        )
        results = b17c.run_analysis(method="ema")

        assert results.n_low_outliers > 0
        assert abs(results.skew_station) <= 3.0

        q100 = b17c.compute_quantiles(aep=np.array([0.01]))["flow_cfs"].iloc[0]
        # Q100 should stay within a physically reasonable multiple of the
        # largest observed peak, not blow up to 10-25x it.
        assert q100 < 5 * peaks.max()


# Unified interface tests
class TestBulletin17C:
    def test_default_method_is_ema(self, synthetic_peaks):
        """Default method should be EMA."""
        b17c = Bulletin17C(synthetic_peaks)
        results = b17c.run_analysis()
        assert results.method == AnalysisMethod.EMA

    def test_mom_method_selection(self, synthetic_peaks):
        """Test MOM method selection."""
        b17c = Bulletin17C(synthetic_peaks)
        results = b17c.run_analysis(method="mom")
        assert results.method == AnalysisMethod.MOM

    def test_ema_method_selection(self, synthetic_peaks):
        """Test EMA method selection."""
        b17c = Bulletin17C(synthetic_peaks)
        results = b17c.run_analysis(method="ema")
        assert results.method == AnalysisMethod.EMA

    def test_property_access(self, synthetic_peaks):
        """Test convenience property access."""
        b17c = Bulletin17C(synthetic_peaks)
        b17c.run_analysis()

        assert b17c.mean_log is not None
        assert b17c.std_log is not None
        assert b17c.skew_station is not None
        assert b17c.quantiles is not None


# Edge cases
class TestEdgeCases:
    def test_small_sample(self):
        """Test with minimum viable sample size."""
        peaks = np.array([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
        mom = MethodOfMoments(peaks)
        results = mom.run_analysis()
        assert results.n_peaks == 10

    def test_handles_zeros(self):
        """Test that zeros are filtered out."""
        peaks = np.array([0, 100, 200, 0, 300, 400])
        mom = MethodOfMoments(peaks)
        results = mom.run_analysis()
        assert results.n_peaks == 4

    def test_handles_nans(self):
        """Test that NaNs are filtered out."""
        peaks = np.array([100, np.nan, 200, 300, np.nan, 400])
        mom = MethodOfMoments(peaks)
        results = mom.run_analysis()
        assert results.n_peaks == 4


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
