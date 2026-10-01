"""Tests for Bulletin 17C flood frequency analysis."""

import logging

import numpy as np
import pytest

from flowfreq import (
    AnalysisMethod,
    Bulletin17C,
    ExpectedMomentsAlgorithm,
    MethodOfMoments,
    grubbs_beck_critical_value,
    kfactor,
)
from flowfreq.bulletin17c import _b17b_skew_mse


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

    def test_historical_perception_threshold_overrides_peak_max(self):
        """Regression test: an explicit perception threshold for the
        historical period must be used as the EMA historical censoring
        threshold, not the maximum of the listed historical peak values.

        These are different quantities: "the biggest flood we happen to
        have a record of" (max of historical peak values) is not "the
        smallest flood that would have been noticed at all" (the actual
        perception threshold). Silently substituting the former censors
        every unlisted historical year against too weak a bound and biases
        every downstream moment. Uses the Big Sandy River fixture values
        (PeakfqSA User Manual, Cohn 2012): historical peaks 1897/1919/1927
        max out at 25,000 cfs, but the documented perception threshold for
        that period is 18,000 cfs.
        """
        water_years = np.arange(1930, 1974)
        peak_flows = np.full(len(water_years), 5000.0)
        historical_peaks = [(1897, 25000.0), (1919, 21000.0), (1927, 18500.0)]
        perception_thresholds = {(1890, 1929): 18000.0}

        ema = ExpectedMomentsAlgorithm(
            peak_flows,
            water_years=water_years,
            historical_peaks=historical_peaks,
            perception_thresholds=perception_thresholds,
        )

        assert ema._ema_params.historical_threshold == 18000.0

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


# MGBT reference validation
class TestMGBTOrestimba:
    """Validate MGBT against the USGS B17C Appendix 10 PILF example.

    Reference: Bulletin 17C Appendix 10 — Orestimba Creek near Newman, CA
    (USGS 11274500, WY 1932-2013).  The MGBT should identify 782 cfs as the
    PILF threshold, censoring 30 peaks (12 zero-flow years + 18 non-zero
    peaks < 782 cfs) with a significance level ≈ 0.0007.
    """

    @pytest.fixture
    def orestimba_data(self):
        """Actual annual peak flows for USGS 11274500 (WY 1932-2013)."""
        peaks = {
            1932: 4260,
            1933: 345,
            1934: 516,
            1935: 1320,
            1936: 1200,
            1937: 2180,
            1938: 3230,
            1939: 115,
            1940: 3440,
            1941: 3070,
            1942: 1880,
            1943: 6450,
            1944: 1290,
            1945: 5970,
            1946: 782,
            1947: 0,
            1948: 0,
            1949: 335,
            1950: 175,
            1951: 2920,
            1952: 3660,
            1953: 147,
            1954: 0,
            1955: 16,
            1956: 5620,
            1957: 1440,
            1958: 10200,
            1959: 5380,
            1960: 448,
            1961: 0,
            1962: 1740,
            1963: 8300,
            1964: 156,
            1965: 560,
            1966: 128,
            1967: 4200,
            1968: 0,
            1969: 5080,
            1970: 1010,
            1971: 584,
            1972: 0,
            1973: 1510,
            1974: 922,
            1975: 1010,
            1976: 0,
            1977: 0,
            1978: 4360,
            1979: 1270,
            1980: 5210,
            1981: 1130,
            1982: 5550,
            1983: 6360,
            1984: 991,
            1985: 50,
            1986: 6990,
            1987: 112,
            1988: 0,
            1989: 0,
            1990: 4,
            1991: 1260,
            1992: 888,
            1993: 4190,
            1994: 12,
            1995: 12000,
            1996: 3130,
            1997: 3320,
            1998: 9470,
            1999: 833,
            2000: 2550,
            2001: 958,
            2002: 425,
            2003: 2790,
            2004: 2990,
            2005: 1820,
            2006: 1630,
            2007: 0,
            2008: 2110,
            2009: 310,
            2010: 4400,
            2011: 4440,
            2012: 0,
            2013: 6250,
        }
        wys = np.array(sorted(peaks.keys()))
        flows = np.array([peaks[y] for y in wys], dtype=float)
        return flows, wys

    def test_mgbt_threshold_782(self, orestimba_data):
        """MGBT threshold must equal 782 cfs (B17C Appendix 10 reference)."""
        flows, wys = orestimba_data
        b = Bulletin17C(
            peak_flows=flows,
            water_years=wys,
            regional_skew=-0.302,
            regional_skew_mse=0.302,
        )
        b.run_analysis(method="ema")
        assert b.results.low_outlier_threshold == pytest.approx(782.0, abs=1.0)

    def test_mgbt_n_low_outliers_30(self, orestimba_data):
        """MGBT must censor exactly 30 peaks (12 zeros + 18 non-zero < 782)."""
        flows, wys = orestimba_data
        b = Bulletin17C(
            peak_flows=flows,
            water_years=wys,
            regional_skew=-0.302,
            regional_skew_mse=0.302,
        )
        b.run_analysis(method="ema")
        assert b.results.n_low_outliers == 30

    def test_mgbt_pilf_includes_12_zeros(self, orestimba_data):
        """PILF list must include all 12 zero-flow years."""
        flows, wys = orestimba_data
        b = Bulletin17C(
            peak_flows=flows,
            water_years=wys,
            regional_skew=-0.302,
            regional_skew_mse=0.302,
        )
        b.run_analysis(method="ema")
        zeros_in_pilf = sum(1 for f in b.results.pilf_flows if f == 0.0)
        assert zeros_in_pilf == 12


class TestMGBTMemoization:
    """Guards on the lru_cache around ExpectedMomentsAlgorithm._mgbt_pvalue.

    The cache is what takes the suite from ~75 s to ~14 s, so it is worth a
    little insurance: one test that it is still in place and reachable, and one
    that it has not changed an answer.
    """

    @staticmethod
    def _descriptor():
        return ExpectedMomentsAlgorithm.__dict__["_mgbt_pvalue"]

    def test_decorator_order_survives(self):
        """staticmethod must stay outermost.

        Reversed, this breaks only on Python 3.9, where a staticmethod object
        is not callable and lru_cache cannot wrap it -- a red matrix job that a
        green local 3.11 run would not predict.
        """
        assert isinstance(self._descriptor(), staticmethod)
        assert hasattr(ExpectedMomentsAlgorithm._mgbt_pvalue, "cache_info")

    def test_cached_value_matches_uncached(self):
        """The cache must be transparent: same key, bit-identical result."""
        cached = ExpectedMomentsAlgorithm._mgbt_pvalue
        for args in ((50, 3, -2.1), (88, 12, -1.4), (30, 1, 0.5)):
            assert cached(*args) == cached.__wrapped__(*args)

    def test_repeated_key_is_served_from_the_cache(self):
        cached = ExpectedMomentsAlgorithm._mgbt_pvalue
        cached(97, 7, -1.9)
        before = cached.cache_info().hits
        cached(97, 7, -1.9)
        assert cached.cache_info().hits == before + 1

    def test_unhashable_argument_is_a_type_error(self):
        """A numpy array would silently defeat the cache; it must raise instead."""
        with pytest.raises(TypeError):
            ExpectedMomentsAlgorithm._mgbt_pvalue(50, 3, np.array([-2.1, -1.0]))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


class TestRegionalSkewInTheFixedPoint:
    """The regional skew enters the EMA loop, it is not averaged in afterwards.

    peakfq 8.1.0 implements Bulletin 17C eq. 7-10 inside moms_p3
    (vendor/peakfqr/src/emafit.f:1344): the regional skew arrives as nG
    pseudo-observations at value rG, evaluated every iteration. Because the
    skew feeds the P3 distribution that produces the next round of expected
    moments, that also moves the mean and the variance -- which a post-hoc
    weighted average of two skews cannot do.
    """

    @staticmethod
    def _big_sandy(**kwargs):
        from tests.fixtures.big_sandy import (
            HISTORICAL_PEAKS,
            REGIONAL_SKEW,
            REGIONAL_SKEW_SD,
            SYSTEMATIC_PEAKS,
            THRESHOLDS,
        )

        params = dict(
            peak_flows=list(SYSTEMATIC_PEAKS.values()),
            water_years=list(SYSTEMATIC_PEAKS.keys()),
            regional_skew=REGIONAL_SKEW,
            regional_skew_mse=REGIONAL_SKEW_SD**2,
            historical_peaks=[(y, f) for y, f in HISTORICAL_PEAKS.items()],
            perception_thresholds={(t["start"], t["end"]): t["lower"] for t in THRESHOLDS},
        )
        params.update(kwargs)
        b = Bulletin17C(**params)
        b.run_analysis()
        return b

    def test_weighting_moves_the_mean_and_variance_too(self):
        """Not just the skew. This is what distinguishes the two formulations."""
        weighted = self._big_sandy().results
        at_site = self._big_sandy(regional_skew=None, regional_skew_mse=None).results
        assert weighted.mean_log != at_site.mean_log
        assert weighted.std_log != at_site.std_log

    def test_equivalent_years_follows_griffis_2004(self):
        """nG = n * Wd * as_G_mse / r_G_mse, emafit.f:1194.

        as_G_mse is ADJE's censoring-adjusted skew MSE
        (``_adje_skew_mse``), not ``_b17b_skew_mse`` alone -- see
        TODO.md P3's "skew weighting" item.
        """
        ema = self._big_sandy()._analyzer
        n = len(ema.intervals)
        mean_log, std_log = ema.results.mean_log, ema.results.std_log
        n_g = ema._regional_skew_equivalent_years(mean_log, std_log, 0.0066, n)
        expected = (
            n * 1.0 * ema._adje_skew_mse(mean_log, std_log, 0.0066, n) / ema._regional_skew_mse
        )
        assert n_g == pytest.approx(expected)

    def test_no_regional_skew_means_no_pseudo_observations(self):
        ema = self._big_sandy(regional_skew=None, regional_skew_mse=None)._analyzer
        assert ema._regional_skew_equivalent_years(3.7, 0.09, 0.0066, 84) == 0.0

    @pytest.mark.parametrize("mse", [0.0, 1e11])
    def test_generalized_and_station_only_get_no_pseudo_observations(self, mse):
        """MSE 0 is 'generalized, no error'; >= 1e10 is station-only."""
        ema = self._big_sandy(regional_skew_mse=mse)._analyzer
        assert ema._regional_skew_equivalent_years(3.7, 0.09, 0.0066, 84) == 0.0

    def test_bias_correction_applies_to_exact_peaks_only(self):
        """c2 and c3 correct the exact peaks; censored intervals go in raw.

        With an uncensored record the two formulations coincide, which is what
        makes this checkable: the Fortran split must reproduce the plain
        sample moments exactly.
        """
        flows = np.array([9100.0, 2060, 7820, 3220, 5580, 17000, 6740, 13800, 4270, 5940])
        years = np.arange(1960, 1970)
        b = Bulletin17C(peak_flows=flows, water_years=years)
        b.run_analysis(method="ema")
        logs = np.log10(flows)
        n = len(logs)
        assert b.results.mean_log == pytest.approx(logs.mean())
        assert b.results.std_log == pytest.approx(logs.std(ddof=1))
        m3 = np.sum((logs - logs.mean()) ** 3)
        expected_skew = n * m3 / ((n - 1) * (n - 2) * logs.std(ddof=1) ** 3)
        assert b.results.skew_station == pytest.approx(expected_skew)


class TestB17BSkewMse:
    """mseg(), emafit.f:1707."""

    def test_matches_the_fortran_formula(self):
        # a = -0.33 + 0.08|g|, b = 0.94 - 0.26|g|, mseg = 10**(a - b*log10(n/10))
        expected = 10 ** ((-0.33 + 0.08 * 0.5) - (0.94 - 0.26 * 0.5) * np.log10(8.4))
        assert _b17b_skew_mse(84, 0.5) == pytest.approx(expected)

    def test_large_skew_switches_coefficients(self):
        """|g| > 0.9 changes a; |g| > 1.5 pins b at 0.55."""
        expected = 10 ** ((-0.52 + 0.30 * 1.6) - 0.55 * np.log10(8.4))
        assert _b17b_skew_mse(84, 1.6) == pytest.approx(expected)

    def test_sign_of_skew_is_irrelevant(self):
        assert _b17b_skew_mse(84, -0.5) == pytest.approx(_b17b_skew_mse(84, 0.5))

    def test_record_length_is_capped_at_150(self):
        """mseg_all passes min(n, 150); a longer record buys no more certainty."""
        assert _b17b_skew_mse(400, 0.2) == pytest.approx(_b17b_skew_mse(150, 0.2))

    def test_mse_falls_with_record_length(self):
        assert _b17b_skew_mse(100, 0.2) < _b17b_skew_mse(20, 0.2)


class TestMomUserLowOutlierThreshold:
    """MOM used to drop a user PILF threshold silently, and never censored on it.

    Bulletin17C.run_analysis did not pass user_low_outlier_threshold down the
    MOM path at all, so a user who set one and got a MOM fit -- which is what
    happens when EMA does not converge and flowfreq.workflow falls back -- saw a
    Grubbs-Beck value they had not asked for, with no indication their setting
    had been ignored. MOM also computed its moments from every peak regardless
    of the threshold.

    Both are fixed: the override is passed down, and MOM now applies the
    Bulletin 17B conditional-probability adjustment -- peaks below the
    threshold are dropped from the moments, and quantiles are evaluated at an
    adjusted exceedance probability that accounts for the omitted peaks.
    """

    FLOWS = np.array(
        [
            9100.0,
            2060,
            7820,
            3220,
            5580,
            17000,
            6740,
            13800,
            4270,
            5940,
            1680,
            1200,
            10100,
            3780,
            5340,
            5630,
            12000,
            3980,
            6130,
            4740,
        ]
    )
    YEARS = np.arange(1930, 1950)

    def _run(self, override=None):
        b17c = Bulletin17C(
            peak_flows=self.FLOWS,
            water_years=self.YEARS,
            regional_skew=-0.5,
            regional_skew_mse=0.3025,
            user_low_outlier_threshold=override,
        )
        b17c.run_analysis(method="mom")
        return b17c.results

    def test_without_an_override_it_reports_grubbs_beck(self):
        results = self._run()
        assert results.low_outlier_threshold > 0
        assert results.low_outlier_threshold == pytest.approx(
            10
            ** (
                np.mean(np.log10(self.FLOWS))
                - grubbs_beck_critical_value(len(self.FLOWS)) * np.std(np.log10(self.FLOWS), ddof=1)
            )
        )

    def test_override_is_reported_instead_of_grubbs_beck(self):
        assert self._run(override=4000.0).low_outlier_threshold == pytest.approx(4000.0)

    def test_override_changes_the_reported_pilf_count(self):
        assert self._run(override=4000.0).n_low_outliers > self._run().n_low_outliers

    def test_override_logs_that_mom_censors(self, caplog):
        """The fit now acts on the override; say so."""
        with caplog.at_level(logging.INFO, logger="flowfreq.bulletin17c"):
            self._run(override=4000.0)
        assert any("censors" in r.message for r in caplog.records)

    def test_moments_shift_with_the_override(self):
        """The conditional-probability adjustment means the override changes the fit."""
        base, forced = self._run(), self._run(override=4000.0)
        assert forced.mean_log != pytest.approx(base.mean_log)
        assert forced.std_log != pytest.approx(base.std_log)
        assert forced.n_systematic == len(self.FLOWS) - forced.n_low_outliers

    def test_zero_and_none_both_mean_grubbs_beck(self):
        assert self._run(override=0.0).low_outlier_threshold == pytest.approx(
            self._run().low_outlier_threshold
        )

    def test_conditional_moments_come_from_the_censored_sample_only(self):
        """mean_log/std_log/skew_station must match a fit on the surviving peaks alone."""
        results = self._run(override=4000.0)
        conditional = np.log10(self.FLOWS[self.FLOWS >= 4000.0])
        n_c = len(conditional)
        assert results.n_systematic == n_c
        assert results.mean_log == pytest.approx(np.mean(conditional))
        assert results.std_log == pytest.approx(np.std(conditional, ddof=1))

    def test_quantiles_use_the_conditional_probability(self):
        """Pc = P * n / n_conditional should reproduce the K-factor used at each AEP."""
        b17c = Bulletin17C(
            peak_flows=self.FLOWS,
            water_years=self.YEARS,
            regional_skew=-0.5,
            regional_skew_mse=0.3025,
            user_low_outlier_threshold=4000.0,
        )
        b17c.run_analysis(method="mom")
        results = b17c.results
        n, n_c = len(self.FLOWS), results.n_systematic

        quantiles = b17c.compute_quantiles(aep=np.array([0.10, 0.02]))
        for aep, row in zip([0.10, 0.02], quantiles.itertuples()):
            pc = aep * n / n_c
            expected_K = kfactor(results.skew_used, pc)
            assert row.K_factor == pytest.approx(expected_K)

    def test_high_pilf_fraction_yields_nan_for_undefined_aep(self):
        """An AEP whose conditional probability would exceed 1 has no answer."""
        results = self._run(override=8000.0)
        n, n_c = len(self.FLOWS), results.n_systematic
        assert n_c < n
        # aep=0.995 will fail the Pc < 1 test whenever n / n_c > ~1.005.
        quantiles = MethodOfMoments(
            self.FLOWS,
            regional_skew=-0.5,
            regional_skew_mse=0.3025,
            user_low_outlier_threshold=8000.0,
        ).compute_quantiles(aep=np.array([0.995]))
        if 0.995 * n / n_c >= 1.0:
            assert np.isnan(quantiles["flow_cfs"].iloc[0])

    def test_too_few_conditional_peaks_raises(self):
        """A threshold that censors nearly everything cannot fit a skew."""
        with pytest.raises(ValueError):
            self._run(override=1_000_000.0)


class TestEngineParameter:
    """``Bulletin17C.run_analysis(engine=...)``'s guard clauses.

    The fortran engine itself needs the built f2py extension
    (``tests/fortran_parity/test_interval_builder_live.py`` and
    ``test_compare_engines.py`` exercise that live); everything here checks
    validation that happens before any Fortran call is attempted, so it runs
    with or without the extension built.
    """

    FLOWS = np.array([1000, 1200, 1500, 900, 1100, 1300, 1050, 1400, 950, 1250], dtype=float)
    YEARS = np.arange(2010, 2020)

    def test_default_engine_is_native(self):
        b17c = Bulletin17C(peak_flows=self.FLOWS, water_years=self.YEARS)
        b17c.run_analysis(method="ema")
        assert isinstance(b17c._analyzer, ExpectedMomentsAlgorithm)

    def test_unknown_engine_raises(self):
        b17c = Bulletin17C(peak_flows=self.FLOWS, water_years=self.YEARS)
        with pytest.raises(ValueError, match="engine"):
            b17c.run_analysis(engine="bogus")

    def test_fortran_engine_rejects_method_of_moments(self):
        """peakfq 8.1.0's emafitpr has no MOM path -- this must raise before
        ever trying to import the extension, so it is checkable without one
        built."""
        b17c = Bulletin17C(peak_flows=self.FLOWS, water_years=self.YEARS)
        with pytest.raises(ValueError, match="Method-of-Moments"):
            b17c.run_analysis(method="mom", engine="fortran")


class TestEMAZeroFlowYears:
    """Zero (and NaN) flows are dropped from the fit; their years must go with them.

    Regression test: the base class filtered ``_peak_flows`` but EMA kept every
    water year, so ``zip(flows, years)`` paired each peak after the first zero
    with the wrong year. Found by the .psf converter's row check (PR #52), which
    refused 15 of 24 WY/MT stations for it.
    """

    FLOWS = [0.0, 100.0, 200.0, 0.0, 300.0, np.nan, 400.0]
    YEARS = [1990, 1991, 1992, 1993, 1994, 1995, 1996]

    def _ema(self, **kwargs):
        from flowfreq.bulletin17c import ExpectedMomentsAlgorithm

        return ExpectedMomentsAlgorithm(np.array(self.FLOWS), np.array(self.YEARS), **kwargs)

    def test_each_peak_keeps_its_own_year(self):
        intervals = self._ema()._build_flow_intervals()
        assert {(i.year, i.lower) for i in intervals if i.lower > 1.0} == {
            (1991, 100.0),
            (1992, 200.0),
            (1994, 300.0),
            (1996, 400.0),
        }

    def test_zero_years_get_a_row_like_sitqt(self):
        """siteQT records a zero exactly at Qmin (1e-20); the NaN year gets none."""
        from flowfreq.bulletin17c import _PERCEPTION_QMIN

        rows = {i.year: i for i in self._ema()._build_flow_intervals()}
        assert set(rows) == {1990, 1991, 1992, 1993, 1994, 1996}
        for year in (1990, 1993):
            assert rows[year].lower == rows[year].upper == _PERCEPTION_QMIN

    def test_mgbt_censors_zero_years_below_the_low_threshold(self):
        """gbtest: a zero's Qmin row lies below any cutoff, so it becomes (0, cutoff)."""
        rows = {i.year: i for i in self._ema()._build_flow_intervals(low_threshold=150.0)}
        for year in (1990, 1993):
            assert rows[year].is_censored
            assert (rows[year].lower, rows[year].upper) == (1e-6, 150.0)  # (gbtmin, cutoff)
        assert rows[1991].is_censored  # 100 < 150: an ordinary PILF
        assert not rows[1992].is_censored

    def test_a_fit_with_zeros_counts_them_as_rows(self):
        """End to end: MGBT treats zeros as low outliers and each gets a row."""
        from flowfreq.bulletin17c import ExpectedMomentsAlgorithm

        rng = np.random.default_rng(3)
        flows = 10 ** rng.normal(3.0, 0.3, 30)
        flows[[4, 11, 20]] = 0.0
        years = np.arange(1990, 2020)
        ema = ExpectedMomentsAlgorithm(flows, years)
        results = ema.run_analysis()
        assert results.n_peaks == 30
        zero_rows = [i for i in ema._intervals if i.year in (1994, 2001, 2010)]
        assert len(zero_rows) == 3 and all(i.is_censored for i in zero_rows)
        assert results.n_low_outliers >= 3

    def test_zero_flow_years_are_recorded_not_gaps(self):
        ema = self._ema()
        # Only the NaN year is missing; the zero years were observed.
        assert ema._ema_params.systematic_start == 1990
        assert ema._ema_params.systematic_end == 1996
        assert ema._ema_params.historical_threshold is not None  # gap at 1995 only
        assert ema._ema_params.historical_end == 1994

    def test_zero_years_are_not_censored_as_historical(self):
        ema = self._ema(perception_thresholds={(1985, 1996): 50.0})
        censored_years = {i.year for i in ema._build_flow_intervals() if i.is_censored}
        assert 1990 not in censored_years
        assert 1993 not in censored_years

    def test_length_mismatch_raises(self):
        from flowfreq.bulletin17c import ExpectedMomentsAlgorithm

        with pytest.raises(ValueError, match="same length"):
            ExpectedMomentsAlgorithm(np.array([1.0, 2.0]), np.array([2000]))


class TestAtSiteSkewMseOption:
    """emafit.f:707-711: B17B skew MSE when MGBT computes the cut and finds PILFs."""

    @staticmethod
    def _flows(with_pilfs: bool) -> np.ndarray:
        rng = np.random.default_rng(7)
        flows = 10 ** rng.normal(3.0, 0.25, 40)
        if with_pilfs:
            flows[:4] = [5.0, 6.0, 7.0, 8.0]
        return flows

    def _ema(self, flows, **kwargs):
        from flowfreq.bulletin17c import ExpectedMomentsAlgorithm

        ema = ExpectedMomentsAlgorithm(
            flows,
            np.arange(1980, 1980 + len(flows)),
            regional_skew=-0.1,
            regional_skew_mse=0.12,
            **kwargs,
        )
        ema.run_analysis()
        return ema

    def test_mgbt_with_low_outliers_uses_b17b(self):
        ema = self._ema(self._flows(with_pilfs=True))
        assert ema._results.n_low_outliers > 0
        assert ema._at_site_option == "B17B"

    def test_mgbt_without_low_outliers_uses_adje(self):
        ema = self._ema(self._flows(with_pilfs=False))
        assert ema._results.n_low_outliers == 0
        assert ema._at_site_option == "ADJE"

    def test_user_threshold_keeps_adje_even_with_low_outliers(self):
        """A supplied threshold is peakfq's FIXED (gbthrsh0 > -6): no switch."""
        ema = self._ema(self._flows(with_pilfs=True), user_low_outlier_threshold=10.0)
        assert ema._results.n_low_outliers > 0
        assert ema._at_site_option == "ADJE"

    def test_b17b_mode_uses_uncapped_mseg(self):
        from flowfreq.bulletin17c import _mseg

        ema = self._ema(self._flows(with_pilfs=True))
        mse = ema._at_site_skew_mse(3.0, 0.25, -0.3, 200)
        assert mse == pytest.approx(_mseg(200, -0.3))

    def test_mseg_caps_only_in_the_adje_helper(self):
        from flowfreq.bulletin17c import _b17b_skew_mse, _mseg

        assert _b17b_skew_mse(200, 0.2) == pytest.approx(_mseg(150, 0.2))
        assert _mseg(200, 0.2) < _mseg(150, 0.2)


class TestNearZeroSkewConfidenceBounds:
    """emafit.f section 4.2: interpolate the bounds between var_emab at +/-skewmin."""

    def _ema(self):
        from flowfreq.bulletin17c import ExpectedMomentsAlgorithm

        rng = np.random.default_rng(1)
        flows = 10 ** rng.normal(3.0, 0.25, 50)
        ema = ExpectedMomentsAlgorithm(
            flows, np.arange(1970, 2020), regional_skew=0.0, regional_skew_mse=0.02
        )
        ema.run_analysis()
        return ema

    def test_bounds_are_the_linear_blend_of_the_two_edges(self):
        from flowfreq.bulletin17c import _CI_SKEWMIN

        ema = self._ema()
        g = ema._results.skew_used
        assert abs(g) < _CI_SKEWMIN
        aep = np.array([0.1, 0.01])
        ci = ema.compute_confidence_limits(aep)

        nobs, tl, tu = ema._perception_threshold_groups()
        args = (
            tuple(nobs),
            tuple(tl),
            tuple(tu),
            float(ema._results.mean_log),
            float(ema._results.std_log) ** 2,
        )
        rest = (tuple(1.0 - aep), 0.90, 0.02, ema._at_site_option)
        lo_n, hi_n = ema._cohn_confidence_bounds(*args, -_CI_SKEWMIN, *rest)
        lo_p, hi_p = ema._cohn_confidence_bounds(*args, _CI_SKEWMIN, *rest)
        wt = (g + _CI_SKEWMIN) / (2 * _CI_SKEWMIN)
        np.testing.assert_allclose(
            np.log10(ci["lower_5pct"].to_numpy()), (1 - wt) * lo_n + wt * lo_p, rtol=1e-12
        )
        np.testing.assert_allclose(
            np.log10(ci["upper_5pct"].to_numpy()), (1 - wt) * hi_n + wt * hi_p, rtol=1e-12
        )

    def test_quantiles_stay_exact_at_the_skew_itself(self):
        """Only the bounds are interpolated; the quantile uses the real skew."""
        from flowfreq.core import kfactor

        ema = self._ema()
        r = ema._results
        q = ema.compute_quantiles(np.array([0.01]))["flow_cfs"].iloc[0]
        expected = 10 ** (r.mean_log + kfactor(r.skew_used, 0.01) * r.std_log)
        assert q == pytest.approx(expected, rel=1e-10)


class TestLowOutlierCounts:
    """``n_low_outliers`` is peakfq's ``gbnlow``; ``n_mgbt_outliers`` the flagged peaks."""

    @staticmethod
    def _record():
        rng = np.random.default_rng(7)
        years = np.arange(1960, 2010)
        flows = np.round(10 ** rng.normal(3.0, 0.35, len(years)), 0)
        gaps = [1970, 1971, 1972]  # no peak, censored below a 50 cfs perception threshold
        keep = ~np.isin(years, gaps)
        return flows[keep], years[keep], gaps

    def test_gap_years_below_the_cutoff_are_low_outlier_rows(self):
        flows, years, gaps = self._record()
        # Above the smallest peak, so gbtest does not put the gap years into
        # the test's own sample; below the 300 cfs cutoff, so it censors them.
        perception = (float(flows.min()) + 300.0) / 2.0
        b = Bulletin17C(
            flows,
            years,
            user_low_outlier_threshold=300.0,
            perception_thresholds={(1960, 2009): perception},
        )
        r = b.run_analysis(method="ema")
        below = int(np.sum(flows < 300.0))
        assert below > 0
        assert r.n_mgbt_outliers == below
        assert r.n_low_outliers == below + len(gaps)
        assert b.n_mgbt_outliers == below

    def test_less_than_gap_years_are_in_the_tests_sample(self):
        """A gap year censored no higher than the smallest peak joins MGBT's sample."""
        flows, years, gaps = self._record()
        r = Bulletin17C(
            flows,
            years,
            user_low_outlier_threshold=300.0,
            perception_thresholds={(1960, 2009): float(flows.min()) / 2.0},
        ).run_analysis(method="ema")
        below = int(np.sum(flows < 300.0))
        assert r.n_mgbt_outliers == r.n_low_outliers == below + len(gaps)

    def test_gap_years_above_the_cutoff_are_not(self):
        flows, years, _gaps = self._record()
        r = Bulletin17C(
            flows,
            years,
            user_low_outlier_threshold=300.0,
            perception_thresholds={(1960, 2009): 400.0},
        ).run_analysis(method="ema")
        assert r.n_low_outliers == r.n_mgbt_outliers == int(np.sum(flows < 300.0))

    def test_method_of_moments_reports_both_alike(self):
        flows, years, _gaps = self._record()
        r = Bulletin17C(flows, years, user_low_outlier_threshold=300.0).run_analysis(method="mom")
        assert r.n_low_outliers == r.n_mgbt_outliers == int(np.sum(flows < 300.0))


class TestP3estEmaFixedPoint:
    """``_ema_fixed_point`` follows ``p3est_ema``'s start and stopping rule."""

    @staticmethod
    def _ema(**extra):
        rng = np.random.default_rng(7)
        years = np.arange(1960, 2010)
        flows = np.round(10 ** rng.normal(3.0, 0.35, len(years)), 0)
        ema = ExpectedMomentsAlgorithm(flows, years, user_low_outlier_threshold=300.0, **extra)
        ema.run_analysis()
        return ema

    @staticmethod
    def _dist_p3(prev, cur):
        # emafit.f dist_p3, on (mean, variance, skew)
        return (
            (prev[0] - cur[0]) ** 2 / cur[1]
            + (prev[1] - cur[1]) ** 2 / cur[1] ** 2 / 10.0
            + (prev[2] - cur[2]) ** 2 / 100.0
        )

    def test_converges_to_the_dist_p3_floor(self):
        ema = self._ema()
        mean, std, skew, converged, iterations = ema._ema_fixed_point(9.0, 9.0, 9.0)
        assert converged
        assert 1 < iterations < ema._ema_params.max_iterations
        nxt = ema._ema_iteration(mean, std, skew)
        assert self._dist_p3((mean, std**2, skew), (nxt[0], nxt[1] ** 2, nxt[2])) <= 1e-10

    def test_start_is_ignored(self):
        """p3est_ema always starts from (0, 1, 0), so the given start cannot matter."""
        ema = self._ema()
        a = ema._ema_fixed_point(0.0, 1.0, 0.0)
        b = ema._ema_fixed_point(2.9, 0.3, 0.05)
        assert a == b

    def test_iteration_cap_reports_non_convergence(self):
        ema = self._ema()
        ema._ema_params.max_iterations = 2
        *_moments, converged, iterations = ema._ema_fixed_point(0.0, 1.0, 0.0)
        assert not converged
        assert iterations == 2

    def test_defaults_are_p3est_emas(self):
        params = self._ema()._ema_params
        assert params.max_iterations == 20000
        assert params.tolerance == 1e-10
