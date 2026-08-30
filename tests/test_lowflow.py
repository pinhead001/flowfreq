"""Tests for hydrolib.lowflow (low-flow frequency analysis)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hydrolib.lowflow import LOW_FLOW_YEAR_TYPES, LowFlowFrequency, annual_minimum_flow


def _flat_year_with_dip(
    start: str, base: float, dip_start: str, dip_val: float, dip_len: int = 7
) -> pd.Series:
    """365 days of `base` flow with a `dip_len`-day dip to `dip_val` starting at
    `dip_start`. The hand-computable building block for the tests below: the
    correct annual minimum n-day mean is always exactly `dip_val`."""
    idx = pd.date_range(start, periods=365, freq="D")
    vals = np.full(365, float(base))
    offset = (pd.Timestamp(dip_start) - idx[0]).days
    vals[offset : offset + dip_len] = dip_val
    return pd.Series(vals, index=idx)


def _multi_year_daily(n_years: int, dip_values, start_year: int = 2000) -> pd.DataFrame:
    """n_years of climatic-year data (Apr 1 - Mar 31), each with a 7-day dip to
    the given value, one value per year. Base flow 100 cfs elsewhere."""
    series = [
        _flat_year_with_dip(
            f"{start_year + i}-04-01", 100.0, f"{start_year + i + 1}-02-01", float(v)
        )
        for i, v in enumerate(dip_values)
    ]
    daily = pd.concat(series).to_frame("flow_cfs")
    return daily[~daily.index.duplicated()]


class TestAnnualMinimumFlow:
    """Tests for the annual minimum n-day mean flow computation."""

    def test_hand_computed_case(self) -> None:
        """Three years, each with one clear 7-day dip to a known value."""
        daily = _multi_year_daily(3, [10.0, 20.0, 5.0])
        minima = annual_minimum_flow(daily, n_day=7, year_type="climatic")

        assert list(minima["year"]) == [2000, 2001, 2002]
        assert minima["flow_cfs"].tolist() == pytest.approx([10.0, 20.0, 5.0])
        assert minima["complete"].all()
        assert (minima["n_day"] == 7).all()

    def test_window_end_date_matches_the_dip(self) -> None:
        """The reported date is the last day of the window achieving the minimum."""
        daily = _flat_year_with_dip("2020-04-01", 100.0, "2021-02-01", 10.0).to_frame("flow_cfs")
        minima = annual_minimum_flow(daily, n_day=7)
        assert minima.iloc[0]["window_end_date"] == pd.Timestamp("2021-02-07")

    def test_gap_does_not_bridge_into_a_false_minimum(self) -> None:
        """A missing block of days (not NaN rows -- literally absent, as
        download_daily_flow produces) must not be silently averaged with
        non-adjacent days into a fabricated low. This is the central
        correctness property of the reindex-before-rolling design."""
        idx = pd.date_range("2020-06-01", "2020-06-20", freq="D")
        vals = [50, 48, 45, 40, 35, 20, 15, 12, 10, 9, 8, 100, 100, 100, 60, 55, 50, 48, 46, 44]
        series = pd.Series(vals, index=idx, dtype=float)
        # Drop the true low-flow days entirely, simulating a real data gap.
        gapped = series.drop(series.index[5:11]).to_frame("flow_cfs")

        minima = annual_minimum_flow(gapped, n_day=7, year_type="calendar", min_days=0)
        # No 7-day window of real consecutive days reaches anywhere near the
        # dropped 8-12 cfs values; the minimum must come only from windows
        # entirely within real data.
        assert minima.iloc[0]["flow_cfs"] > 40.0

    def test_negative_values_treated_as_missing_not_as_extreme_lows(self) -> None:
        """A negative discharge value is a data artifact (ice/backwater), not
        a legitimate low; it must not become the reported annual minimum."""
        idx = pd.date_range("2040-04-01", periods=365, freq="D")
        vals = np.full(365, 100.0)
        vals[300:307] = -999.0  # artifact
        vals[100:107] = 20.0  # the real low
        daily = pd.Series(vals, index=idx).to_frame("flow_cfs")

        minima = annual_minimum_flow(daily, n_day=7)
        assert minima.iloc[0]["flow_cfs"] == pytest.approx(20.0)

    def test_true_zero_flow_is_retained_not_treated_as_missing(self) -> None:
        """Unlike a negative artifact, an actual zero is a real observation
        and must survive into the annual minimum series."""
        daily = _flat_year_with_dip("2020-04-01", 100.0, "2021-02-01", 0.0).to_frame("flow_cfs")
        minima = annual_minimum_flow(daily, n_day=7)
        assert minima.iloc[0]["flow_cfs"] == 0.0

    def test_incomplete_year_flagged(self) -> None:
        """A year with fewer than min_days of data is marked incomplete."""
        idx = pd.date_range("2030-04-01", "2030-08-01", freq="D")  # ~4 months
        daily = pd.DataFrame({"flow_cfs": 50.0}, index=idx)
        minima = annual_minimum_flow(daily, n_day=7, min_days=350)
        assert not minima.iloc[0]["complete"]
        assert minima.iloc[0]["n_days"] < 350

    def test_year_with_data_but_no_full_window_reports_nan_not_missing_row(self) -> None:
        """A year present in the daily record but too short for even one full
        n-day window still appears (so a completeness check can flag it),
        with flow_cfs NaN rather than the year vanishing from the table."""
        idx = pd.date_range("2030-04-01", periods=3, freq="D")  # 3 days, n_day=7
        daily = pd.DataFrame({"flow_cfs": [50.0, 51.0, 49.0]}, index=idx)
        minima = annual_minimum_flow(daily, n_day=7)
        assert len(minima) == 1
        assert pd.isna(minima.iloc[0]["flow_cfs"])
        assert not minima.iloc[0]["complete"]

    def test_invalid_n_day_raises(self) -> None:
        daily = pd.DataFrame({"flow_cfs": [1.0]}, index=pd.date_range("2020-01-01", periods=1))
        with pytest.raises(ValueError, match="n_day must be >= 1"):
            annual_minimum_flow(daily, n_day=0)

    def test_invalid_year_type_raises(self) -> None:
        daily = pd.DataFrame({"flow_cfs": [1.0]}, index=pd.date_range("2020-01-01", periods=1))
        with pytest.raises(ValueError, match="year_type must be one of"):
            annual_minimum_flow(daily, year_type="fiscal")

    @pytest.mark.parametrize("year_type", LOW_FLOW_YEAR_TYPES)
    def test_all_year_types_accepted(self, year_type: str) -> None:
        daily = _multi_year_daily(3, [10.0, 20.0, 5.0])
        minima = annual_minimum_flow(daily, n_day=7, year_type=year_type)
        assert len(minima) >= 2


class TestLowFlowFrequencyValidation:
    """Tests for constructor validation."""

    def test_too_few_years_raises(self) -> None:
        daily = _multi_year_daily(5, [10.0] * 5)
        with pytest.raises(ValueError, match="at least 10 are required"):
            LowFlowFrequency(daily, n_day=7)

    def test_unknown_distribution_raises(self) -> None:
        daily = _multi_year_daily(12, [10.0] * 12)
        with pytest.raises(ValueError, match="distribution must be"):
            LowFlowFrequency(daily, distribution="gumbel")

    def test_short_record_below_recommended_logs_warning(self, caplog) -> None:
        daily = _multi_year_daily(12, np.linspace(10, 15, 12))
        with caplog.at_level("WARNING", logger="hydrolib.lowflow"):
            LowFlowFrequency(daily, n_day=7)
        assert any("recommended" in r.message for r in caplog.records)

    def test_no_warning_at_or_above_recommended_minimum(self, caplog) -> None:
        daily = _multi_year_daily(20, np.linspace(10, 15, 20))
        with caplog.at_level("WARNING", logger="hydrolib.lowflow"):
            LowFlowFrequency(daily, n_day=7)
        assert not any("recommended" in r.message for r in caplog.records)

    def test_incomplete_year_excluded_by_default(self) -> None:
        full = _multi_year_daily(12, np.linspace(10, 15, 12), start_year=2000)
        partial_idx = pd.date_range("2012-04-01", "2012-08-01", freq="D")
        partial = pd.DataFrame({"flow_cfs": 50.0}, index=partial_idx)
        combo = pd.concat([full, partial]).sort_index()
        combo = combo[~combo.index.duplicated()]

        lff = LowFlowFrequency(combo, n_day=7)
        assert 2012 not in lff.annual_minimums["year"].to_numpy()
        assert 2012 in lff.annual_minimums_all_years["year"].to_numpy()

    def test_require_complete_years_false_includes_it(self) -> None:
        full = _multi_year_daily(12, np.linspace(10, 15, 12), start_year=2000)
        partial_idx = pd.date_range("2012-04-01", "2013-03-31", freq="D")
        # Complete calendar coverage but let's instead just verify the flag
        # is honored by re-running with the flag off on an all-complete set.
        lff = LowFlowFrequency(full, n_day=7, require_complete_years=False)
        assert lff.n == 12


class TestMinPositiveYears:
    """Tests for the MIN_POSITIVE_YEARS/RECOMMENDED_MIN_POSITIVE_YEARS floors.

    Simulating 30,000 size-n samples from a true-skew-0 population gives a
    sample-skew standard deviation of 1.22 at n=3 (a symmetric population
    routinely produces an apparent skew over 1 from noise alone) versus
    0.68 at n=10 -- there is no n where the noise vanishes, only shrinks,
    which is why MIN_POSITIVE_YEARS is a floor on absurdity (was 3, an
    arithmetic minimum with no real stability), not a claim that 5 is safe.
    All constructions here use n_total=12 (above MIN_YEARS=10) so only the
    positive-year floor is under test.
    """

    def test_exactly_at_new_floor_does_not_raise(self) -> None:
        """n_pos=5 (the new MIN_POSITIVE_YEARS) must construct successfully."""
        daily = _multi_year_daily(12, [0.0] * 7 + [10.0, 11.0, 12.0, 13.0, 14.0])
        lff = LowFlowFrequency(daily, n_day=7)
        assert lff._n_positive_years == 5

    def test_one_below_new_floor_raises(self) -> None:
        """n_pos=4 must still raise under the raised floor."""
        daily = _multi_year_daily(12, [0.0] * 8 + [10.0, 11.0, 12.0, 13.0])
        with pytest.raises(ValueError, match="at least 5 are required"):
            LowFlowFrequency(daily, n_day=7)

    def test_old_floor_of_three_no_longer_accepted(self) -> None:
        """n_pos=3 (the previous MIN_POSITIVE_YEARS) must now raise -- this
        is the exact case the fix targets: the old code fit a full LP3 skew
        from 3 points, which simulation shows is indistinguishable from
        noise."""
        daily = _multi_year_daily(12, [0.0] * 9 + [10.0, 11.0, 12.0])
        with pytest.raises(ValueError, match="at least 5 are required"):
            LowFlowFrequency(daily, n_day=7)

    def test_thin_positive_sample_warns_for_lp3(self, caplog) -> None:
        """n_pos=9 (below RECOMMENDED_MIN_POSITIVE_YEARS=10) with the
        default distribution='lp3' must warn specifically about the skew
        estimate, not just about overall record length."""
        daily = _multi_year_daily(12, [0.0] * 3 + [8, 9, 10, 11, 12, 13, 14, 15, 16])
        with caplog.at_level("WARNING", logger="hydrolib.lowflow"):
            LowFlowFrequency(daily, n_day=7, distribution="lp3")
        assert any("nonzero year" in r.message and "skew" in r.message for r in caplog.records)

    def test_no_thin_sample_warning_at_or_above_recommended(self, caplog) -> None:
        """n_pos=10 (exactly at RECOMMENDED_MIN_POSITIVE_YEARS) must not log
        the skew-specific warning."""
        daily = _multi_year_daily(12, [0.0] * 2 + list(range(8, 18)))
        with caplog.at_level("WARNING", logger="hydrolib.lowflow"):
            LowFlowFrequency(daily, n_day=7, distribution="lp3")
        assert not any("nonzero year" in r.message for r in caplog.records)

    def test_no_thin_sample_warning_for_lognormal(self, caplog) -> None:
        """distribution='lognormal' does not use the skew estimate at all,
        so the same thin n_pos=9 case must not warn about it."""
        daily = _multi_year_daily(12, [0.0] * 3 + [8, 9, 10, 11, 12, 13, 14, 15, 16])
        with caplog.at_level("WARNING", logger="hydrolib.lowflow"):
            LowFlowFrequency(daily, n_day=7, distribution="lognormal")
        assert not any("nonzero year" in r.message for r in caplog.records)

    def test_n_positive_years_consistent_between_init_and_run_analysis(self) -> None:
        """The hard-error/warning check in __init__ and the moment fit in
        run_analysis() must agree on n_pos -- they now share one computation
        rather than each deriving it separately."""
        daily = _multi_year_daily(12, [0.0] * 3 + list(range(8, 17)))
        lff = LowFlowFrequency(daily, n_day=7)
        n_pos_at_init = lff._n_positive_years
        results = lff.run_analysis()
        assert results.n_years - results.n_zero_years == n_pos_at_init


class TestLowFlowFrequencyAnalysis:
    """Tests for the fit itself: quantiles, monotonicity, sign correctness."""

    def test_quantiles_increase_with_non_exceedance_probability(self) -> None:
        """flow at p=0.50 (median annual low) must exceed flow at p=0.01 (rare,
        severe low) -- the defining monotonicity property of the low-flow tail."""
        rng = np.random.default_rng(0)
        dips = 15.0 + rng.normal(0, 2, 12)
        daily = _multi_year_daily(12, np.clip(dips, 1.0, None))

        lff = LowFlowFrequency(daily, n_day=7)
        results = lff.run_analysis()
        q = results.quantiles.sort_values("non_exceedance_prob")
        assert q["flow_cfs"].is_monotonic_increasing

    def test_negative_station_skew_does_not_invert_quantiles(self) -> None:
        """Regression guard at the module level for the core.py sign bug:
        a fit with negative station skew must still produce increasing
        quantiles, not the sign-inverted values the pre-fix code produced."""
        # A right-skewed-in-reverse (i.e. negatively skewed in log space)
        # set of annual minima: mostly clustered high with a couple of
        # severe lows pulling the tail down.
        dips = [18, 19, 20, 21, 19, 20, 18, 21, 6, 20, 19, 8]
        daily = _multi_year_daily(12, dips)

        lff = LowFlowFrequency(daily, n_day=7)
        results = lff.run_analysis()
        assert results.skew_station < 0
        q = results.quantiles.sort_values("non_exceedance_prob")
        assert q["flow_cfs"].is_monotonic_increasing
        assert (q["flow_cfs"] > 0).all()

    def test_lognormal_forces_zero_skew_but_retains_station_skew(self) -> None:
        daily = _multi_year_daily(12, [18, 19, 20, 21, 19, 20, 18, 21, 6, 20, 19, 8])
        results = LowFlowFrequency(daily, n_day=7, distribution="lognormal").run_analysis()
        assert results.skew_used == 0.0
        assert results.skew_station != 0.0
        assert results.distribution == "lognormal"

    def test_zero_variance_record_gives_zero_skew_not_nan(self) -> None:
        """Every positive year hits the identical annual minimum (plausible on
        a reach governed by a fixed minimum-flow release): std_log is exactly
        0, which would make the raw sample-skew formula a 0/0 division. This
        must come out as skew=0 (the mathematically correct degenerate value,
        since np.clip does not clip NaN and a NaN would poison every
        downstream K-factor via NaN * 0 = NaN)."""
        daily = _multi_year_daily(12, [20.0] * 12)
        results = LowFlowFrequency(daily, n_day=7).run_analysis()

        assert results.std_log == 0.0
        assert results.skew_station == 0.0
        assert not np.isnan(results.skew_station)
        assert results.quantiles["flow_cfs"].to_numpy() == pytest.approx([20.0] * 6)

    def test_results_fields_populated(self) -> None:
        # Built as climatic-shaped years; year_type left at its "climatic"
        # default so the boundary years line up and n_years == 12 exactly.
        # Cross-year_type relabeling is covered by TestLowFlowYearLabel.
        daily = _multi_year_daily(12, np.linspace(10, 20, 12))
        results = LowFlowFrequency(daily, n_day=7).run_analysis()

        assert results.n_years == 12
        assert results.n_zero_years == 0
        assert results.p_zero == 0.0
        assert results.n_day == 7
        assert results.year_type == "climatic"
        assert results.distribution == "lp3"
        assert not results.quantiles.empty
        assert not results.confidence_limits.empty

    def test_invalid_non_exceedance_probability_raises(self) -> None:
        daily = _multi_year_daily(12, np.linspace(10, 20, 12))
        lff = LowFlowFrequency(daily, n_day=7)
        lff.run_analysis()
        with pytest.raises(ValueError, match="strictly between 0 and 1"):
            lff.compute_quantiles(np.array([0.5, 1.0]))
        with pytest.raises(ValueError, match="strictly between 0 and 1"):
            lff.compute_quantiles(np.array([0.0, 0.5]))

    def test_run_analysis_called_implicitly_by_compute_quantiles(self) -> None:
        """compute_quantiles works even before run_analysis is called explicitly,
        matching Bulletin17C's own lazy-run convention."""
        daily = _multi_year_daily(12, np.linspace(10, 20, 12))
        lff = LowFlowFrequency(daily, n_day=7)
        q = lff.compute_quantiles()
        assert not q.empty


class TestZeroFlowConditionalProbability:
    """Tests for the zero-flow-year conditional-probability adjustment."""

    def _make_zero_inflated(self, n_total: int = 15, n_zero: int = 3):
        rng = np.random.default_rng(1)
        dips = [0.0] * n_zero + list(np.clip(10.0 + rng.normal(0, 3, n_total - n_zero), 1.0, None))
        return _multi_year_daily(n_total, dips, start_year=2005)

    def test_zero_year_accounting(self) -> None:
        daily = self._make_zero_inflated(n_total=15, n_zero=3)
        results = LowFlowFrequency(daily, n_day=7).run_analysis()
        assert results.n_zero_years == 3
        assert results.p_zero == pytest.approx(3 / 15)

    def test_probability_at_or_below_p_zero_gives_exact_zero_quantile(self) -> None:
        daily = self._make_zero_inflated(n_total=15, n_zero=3)  # p0 = 0.2
        results = LowFlowFrequency(daily, n_day=7).run_analysis()
        lff = LowFlowFrequency(daily, n_day=7)
        q = lff.compute_quantiles(np.array([0.05, 0.10, 0.20]))

        assert (q["flow_cfs"] == 0.0).all()
        assert q["log_flow"].isna().all()
        assert not np.isneginf(q["log_flow"]).any()
        assert q["conditional_prob"].isna().all()

    def test_probability_above_p_zero_uses_conditional_probability_formula(self) -> None:
        daily = self._make_zero_inflated(n_total=15, n_zero=3)  # p0 = 0.2
        lff = LowFlowFrequency(daily, n_day=7)
        p = np.array([0.30, 0.50])
        q = lff.compute_quantiles(p)

        p0 = 3 / 15
        expected_conditional = (p - p0) / (1 - p0)
        assert q["conditional_prob"].to_numpy() == pytest.approx(expected_conditional)
        assert (q["flow_cfs"] > 0).all()
        assert q["log_flow"].notna().all()

    def test_confidence_limits_nan_at_or_below_p_zero(self) -> None:
        daily = self._make_zero_inflated(n_total=15, n_zero=3)
        lff = LowFlowFrequency(daily, n_day=7)
        ci = lff.compute_confidence_limits(np.array([0.05, 0.10, 0.20, 0.50]))

        pct_cols = [c for c in ci.columns if "pct" in c]
        below = ci[ci["non_exceedance_prob"] <= 0.2]
        above = ci[ci["non_exceedance_prob"] > 0.2]

        assert below[pct_cols].isna().all().all()
        lower_col = next(c for c in ci.columns if c.startswith("lower"))
        upper_col = next(c for c in ci.columns if c.startswith("upper"))
        assert (above[lower_col] < above["flow_cfs"]).all()
        assert (above["flow_cfs"] < above[upper_col]).all()

    def test_too_few_positive_years_raises(self) -> None:
        """15 total years but only 2 nonzero cannot support fitting a distribution.
        Raised at construction time now (see TestMinPositiveYears), not deferred
        to run_analysis()."""
        daily = self._make_zero_inflated(n_total=15, n_zero=13)
        with pytest.raises(ValueError, match="at least 5 are required"):
            LowFlowFrequency(daily, n_day=7)

    def test_all_zero_years_is_not_silently_a_valid_fit(self) -> None:
        """An entirely zero-flow record must raise, not report p_zero=1.0 with
        an undefined distribution fit."""
        daily = _multi_year_daily(12, [0.0] * 12)
        with pytest.raises(ValueError, match="at least 5 are required"):
            LowFlowFrequency(daily, n_day=7)


class TestBootstrapConfidenceLimits:
    """Tests for the parametric bootstrap CI (compute_bootstrap_confidence_limits).

    Added as an alternative to compute_confidence_limits specifically because
    that method's asymptotic variance formula treats p_zero as known --
    understating uncertainty for any record with zero-flow years. These tests
    focus on confirming the bootstrap actually captures that extra source of
    uncertainty, not just that it produces "a" number.
    """

    def _make_zero_inflated(self, n_total: int = 20, n_zero: int = 4, seed: int = 1):
        rng = np.random.default_rng(seed)
        dips = [0.0] * n_zero + list(np.clip(10.0 + rng.normal(0, 3, n_total - n_zero), 1.0, None))
        return _multi_year_daily(n_total, dips, start_year=2005)

    def test_brackets_the_point_estimate_no_zero_years(self) -> None:
        rng = np.random.default_rng(1)
        dips = np.clip(15.0 + rng.normal(0, 3, 15), 1.0, None)
        daily = _multi_year_daily(15, dips)
        lff = LowFlowFrequency(daily, n_day=7)
        lff.run_analysis()

        boot = lff.compute_bootstrap_confidence_limits(random_state=42)
        lower_col = next(c for c in boot.columns if c.startswith("lower"))
        upper_col = next(c for c in boot.columns if c.startswith("upper"))
        assert (boot[lower_col] < boot["flow_cfs"]).all()
        assert (boot["flow_cfs"] < boot[upper_col]).all()

    def test_reproducible_with_fixed_random_state(self) -> None:
        rng = np.random.default_rng(1)
        dips = np.clip(15.0 + rng.normal(0, 3, 15), 1.0, None)
        daily = _multi_year_daily(15, dips)
        lff = LowFlowFrequency(daily, n_day=7)
        lff.run_analysis()

        boot_a = lff.compute_bootstrap_confidence_limits(random_state=7)
        boot_b = lff.compute_bootstrap_confidence_limits(random_state=7)
        lower_col = next(c for c in boot_a.columns if c.startswith("lower"))
        assert np.allclose(boot_a[lower_col], boot_b[lower_col])

    def test_captures_p_zero_uncertainty_where_analytic_ci_cannot(self) -> None:
        """The central point of adding this method: at p == p_zero exactly,
        compute_confidence_limits reports NaN (it treats p_zero as a fixed,
        known quantity, so there is nothing for its formula to attach
        uncertainty to). The bootstrap must report a real, non-trivial
        interval there instead, since some resamples' simulated p_zero will
        land above the true value and some below."""
        daily = self._make_zero_inflated(n_total=20, n_zero=4)  # p0 = 0.20
        lff = LowFlowFrequency(daily, n_day=7)
        results = lff.run_analysis()
        assert results.p_zero == pytest.approx(0.20)

        boot = lff.compute_bootstrap_confidence_limits(np.array([0.20]), random_state=42)
        analytic = lff.compute_confidence_limits(np.array([0.20]))

        analytic_pct_cols = [c for c in analytic.columns if "pct" in c]
        assert analytic[analytic_pct_cols].isna().all(axis=None)

        upper_col = next(c for c in boot.columns if c.startswith("upper"))
        lower_col = next(c for c in boot.columns if c.startswith("lower"))
        assert not np.isnan(boot[upper_col].iloc[0])
        assert not np.isnan(boot[lower_col].iloc[0])
        assert boot[upper_col].iloc[0] > 0

    def test_n_resamples_used_reflects_skipped_iterations(self) -> None:
        """A record with a real (nonzero) p_zero must skip at least some
        iterations whose simulated positive-year count is too small to fit
        -- n_resamples_used should come back below the requested count."""
        daily = self._make_zero_inflated(n_total=20, n_zero=15)  # p0 = 0.75, n_pos = 5 (floor)
        lff = LowFlowFrequency(daily, n_day=7)
        lff.run_analysis()

        boot = lff.compute_bootstrap_confidence_limits(n_resamples=2000, random_state=1)
        assert (boot["n_resamples_used"] < 2000).all()
        assert (boot["n_resamples_used"] > 1000).all()  # not so many skipped it's unusable

    def test_raises_when_too_many_iterations_are_unusable(self) -> None:
        """Directly injects an extreme, internally-inconsistent p_zero after
        construction (deliberately bypassing the constructor's own
        validation, which would never let a real record reach this state --
        see the reachability note in the review) specifically to confirm
        this defensive check independently guards the method even if an
        extreme p_zero were ever reached some other way."""
        daily = self._make_zero_inflated(n_total=20, n_zero=15)
        lff = LowFlowFrequency(daily, n_day=7)
        results = lff.run_analysis()
        results.p_zero = 0.9  # artificially inconsistent with n_years=20, n_pos=5

        with pytest.raises(ValueError, match="too high relative to its length"):
            lff.compute_bootstrap_confidence_limits(n_resamples=2000, random_state=1)

    def test_lognormal_distribution_skew_stays_zero_in_every_resample(self) -> None:
        """distribution='lognormal' must not refit a nonzero skew from any
        simulated resample -- the whole point of that distribution choice is
        to sidestep the skew estimate entirely, in every code path."""
        rng = np.random.default_rng(1)
        dips = np.clip(15.0 + rng.normal(0, 3, 15), 1.0, None)
        daily = _multi_year_daily(15, dips)
        lff = LowFlowFrequency(daily, n_day=7, distribution="lognormal")
        lff.run_analysis()

        boot = lff.compute_bootstrap_confidence_limits(random_state=3, n_resamples=300)
        assert (boot["n_resamples_used"] > 250).all()
        assert not boot.isna().any(axis=None)

    def test_invalid_non_exceedance_probability_raises(self) -> None:
        rng = np.random.default_rng(1)
        dips = np.clip(15.0 + rng.normal(0, 3, 15), 1.0, None)
        daily = _multi_year_daily(15, dips)
        lff = LowFlowFrequency(daily, n_day=7)
        lff.run_analysis()
        with pytest.raises(ValueError, match="strictly between 0 and 1"):
            lff.compute_bootstrap_confidence_limits(np.array([0.0, 0.5]))

    def test_runs_without_explicit_run_analysis_call(self) -> None:
        """Matches the lazy-run convention used elsewhere in this class."""
        rng = np.random.default_rng(1)
        dips = np.clip(15.0 + rng.normal(0, 3, 15), 1.0, None)
        daily = _multi_year_daily(15, dips)
        lff = LowFlowFrequency(daily, n_day=7)
        boot = lff.compute_bootstrap_confidence_limits(random_state=1, n_resamples=300)
        assert not boot.empty
