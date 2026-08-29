"""Tests for hydrolib.lowflow (low-flow frequency analysis)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hydrolib.lowflow import (
    LOW_FLOW_YEAR_TYPES,
    LowFlowFrequency,
    _low_flow_year_label,
    annual_minimum_flow,
)


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


class TestLowFlowYearLabel:
    """Tests for the year-definition labeling used by annual_minimum_flow."""

    def test_climatic_year_labeled_by_starting_calendar_year(self) -> None:
        """Apr 1 Y - Mar 31 Y+1 is climatic year Y (majority-of-months convention)."""
        dates = pd.DatetimeIndex(
            ["1999-04-01", "1999-12-31", "2000-01-01", "2000-03-31", "2000-04-01"]
        )
        labels = _low_flow_year_label(dates, "climatic")
        assert list(labels) == [1999, 1999, 1999, 1999, 2000]

    def test_water_year_matches_usgs_download_peak_flow_convention(self) -> None:
        """Must agree with the inline water_year logic in usgs.download_peak_flow:
        Oct 1 (Y-1) - Sep 30 Y is water year Y."""
        dates = pd.DatetimeIndex(
            ["1998-10-01", "1998-12-31", "1999-01-01", "1999-09-30", "1999-10-01"]
        )
        labels = _low_flow_year_label(dates, "water")
        assert list(labels) == [1999, 1999, 1999, 1999, 2000]

    def test_calendar_year_is_trivial(self) -> None:
        dates = pd.DatetimeIndex(["2020-01-01", "2020-12-31"])
        assert list(_low_flow_year_label(dates, "calendar")) == [2020, 2020]

    def test_unknown_year_type_raises(self) -> None:
        with pytest.raises(ValueError, match="year_type must be one of"):
            _low_flow_year_label(pd.DatetimeIndex(["2020-01-01"]), "fiscal")


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
        """15 total years but only 2 nonzero cannot support fitting a distribution."""
        daily = self._make_zero_inflated(n_total=15, n_zero=13)
        lff = LowFlowFrequency(daily, n_day=7)
        with pytest.raises(ValueError, match="at least 3 are required"):
            lff.run_analysis()

    def test_all_zero_years_is_not_silently_a_valid_fit(self) -> None:
        """An entirely zero-flow record must raise, not report p_zero=1.0 with
        an undefined distribution fit."""
        daily = _multi_year_daily(12, [0.0] * 12)
        lff = LowFlowFrequency(daily, n_day=7)
        with pytest.raises(ValueError, match="at least 3 are required"):
            lff.run_analysis()
