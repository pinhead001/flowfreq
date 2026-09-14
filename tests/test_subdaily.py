"""Tests for flowfreq.subdaily (sub-daily timing and ramping-rate metrics).

Two of these carry most of the weight, because they are the errors that look
correct in the output:

- ``test_circular_mean_near_midnight_is_not_the_arithmetic_mean``: hours 23.5
  and 0.5 must summarize to ~0, not to noon. A future "simplification" back to
  ``.mean()`` passes every other test in this file.
- ``test_dst_fall_back_day_has_no_negative_intervals``: a local wall-clock
  ``dt`` goes backwards across the autumn transition, which turns a modest ramp
  into a large one of the wrong sign. Only this test notices.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from flowfreq.subdaily import (
    DEFAULT_MAX_GAP_MULTIPLE,
    KNOWN_VALUE_UNITS,
    TRUE_ZERO_UNITS,
    _count_reversals,
    _local_day_hours,
    _resolve_units,
    circular_hour_statistics,
    daily_extreme_timing,
    diel_variation,
    extreme_timing_summary,
    ramping_rate_summary,
    ramping_rates,
)

PACIFIC = "America/Los_Angeles"


def _pacific_diel_series(
    n_days: int = 11, amplitude: float = 30.0, mean: float = 100.0
) -> pd.DataFrame:
    """A UTC-indexed instantaneous series with a clean diel sinusoid in LOCAL
    (Pacific) time: trough at 04:00 local, peak at 16:00 local.

    Deliberately the same shape as ``tests/test_regime.py``'s fixture of the
    same name, so the timing functions are exercised on a series whose phase is
    already established by the diel-variation tests.
    """
    local_times = pd.date_range("2020-01-10 00:00", periods=96 * n_days, freq="15min", tz=PACIFIC)
    hour_local = np.array([t.hour + t.minute / 60 for t in local_times])
    flow = mean + amplitude * np.sin(2 * np.pi * (hour_local - 4) / 24 - np.pi / 2)
    return pd.DataFrame({"flow_cfs": flow}, index=local_times.tz_convert("UTC"))


def _flat_day(value: float = 42.0, date: str = "2021-07-01") -> pd.DataFrame:
    return pd.DataFrame(
        {"flow_cfs": np.full(96, value)},
        index=pd.date_range(date, periods=96, freq="15min", tz="UTC"),
    )


class TestCircularHourStatistics:
    """Tests for circular_hour_statistics."""

    def test_circular_mean_near_midnight_is_not_the_arithmetic_mean(self) -> None:
        """The whole reason this function exists.

        A site peaking at 23:30 one day and 00:30 the next peaks at midnight.
        The arithmetic mean of 23.5 and 0.5 is 12.0 -- noon, the one time of
        day it never peaks -- and nothing about that number looks wrong.
        """
        stats = circular_hour_statistics([23.5, 0.5])
        assert stats["mean_hour"] == pytest.approx(0.0, abs=1e-9)
        assert np.mean([23.5, 0.5]) == 12.0  # the answer this must not give
        assert stats["concentration"] > 0.99

    def test_mean_hour_stays_in_range(self) -> None:
        """A mean direction at the wrap point must come back as 0.0, not 24.0.

        atan2 returns a hair below zero there, which `% 2*pi` lifts to almost
        exactly 2*pi -- hour 24.0, outside the documented [0, 24) range.
        """
        for hours in ([23.5, 0.5], [0.0], [23.999], [0.0, 0.0, 23.9999]):
            mean_hour = float(circular_hour_statistics(hours)["mean_hour"])
            assert 0.0 <= mean_hour < 24.0, f"{hours} gave {mean_hour}"

    def test_identical_hours_give_full_concentration_and_zero_spread(self) -> None:
        stats = circular_hour_statistics([16.0] * 11)
        assert stats["mean_hour"] == pytest.approx(16.0)
        assert stats["concentration"] == pytest.approx(1.0)
        # Not -0.0: sqrt(-2*log(1)) is sqrt(-0.0), which prints as a negative
        # standard deviation for the most consistent timing possible.
        assert stats["circular_std_hours"] == 0.0
        assert not np.signbit(stats["circular_std_hours"])
        assert stats["rayleigh_p"] < 0.01

    def test_uniform_hours_give_low_concentration_and_no_significance(self) -> None:
        """Timing spread round the clock must announce itself through `r`.

        The mean direction of a uniform sample is arbitrary; concentration is
        what tells the caller not to read it as a peak time.
        """
        stats = circular_hour_statistics(np.arange(0, 24, 0.25))
        assert stats["concentration"] < 1e-10
        assert stats["rayleigh_p"] > 0.9
        assert stats["circular_std_hours"] > 24.0

    def test_nans_dropped_and_counted(self) -> None:
        stats = circular_hour_statistics([16.0, np.nan, 16.0])
        assert stats["n"] == 2
        assert stats["mean_hour"] == pytest.approx(16.0)

    def test_hours_outside_range_wrapped_not_rejected(self) -> None:
        assert circular_hour_statistics([25.0])["mean_hour"] == pytest.approx(1.0)

    def test_empty_input_gives_nan_not_error(self) -> None:
        stats = circular_hour_statistics([])
        assert stats["n"] == 0
        assert np.isnan(stats["mean_hour"])
        assert np.isnan(stats["rayleigh_p"])

    def test_rayleigh_p_never_exceeds_one(self) -> None:
        """Zar's approximation can overshoot 1.0; a p-value cannot."""
        for hours in ([0.0, 12.0], [0.0, 8.0, 16.0], np.arange(0, 24, 0.5)):
            assert circular_hour_statistics(hours)["rayleigh_p"] <= 1.0


class TestDailyExtremeTiming:
    """Tests for daily_extreme_timing."""

    def test_recovers_known_local_phase(self) -> None:
        """The central correctness property: a sinusoid with its trough at
        04:00 and peak at 16:00 local must be read back at those hours, which
        fails outright if local-day grouping is dropped for the UTC index."""
        timing = daily_extreme_timing(_pacific_diel_series(n_days=11), tz=PACIFIC)
        middle = timing.iloc[2:-2]
        assert middle["hour_of_max"].between(15.5, 16.5).all()
        assert middle["hour_of_min"].between(3.5, 4.5).all()
        assert middle["complete"].all()

    def test_grouping_by_utc_would_give_a_different_answer(self) -> None:
        """Confirms the tz conversion does real work for this input rather
        than being an accidental no-op."""
        iv = _pacific_diel_series(n_days=5)
        utc_timing = daily_extreme_timing(iv, tz="UTC")
        assert not utc_timing["hour_of_max"].between(15.5, 16.5).all()

    def test_flat_day_gives_nan_timing_not_midnight(self) -> None:
        """A constant-release day has no peak hour, and 0.0 would be a
        fabricated midnight one that drags the circular mean."""
        timing = daily_extreme_timing(_flat_day(), tz="UTC")
        assert np.isnan(timing["hour_of_max"].iloc[0])
        assert np.isnan(timing["hour_of_min"].iloc[0])
        assert timing["n_tied_max"].iloc[0] == 96

    def test_tie_reports_first_occurrence_and_counts_it(self) -> None:
        """Two equal crests: the first is returned, and n_tied_max says so."""
        idx = pd.date_range("2021-08-01 00:00", periods=8, freq="3h", tz="UTC")
        values = np.array([10.0, 20.0, 50.0, 30.0, 25.0, 50.0, 20.0, 15.0])
        timing = daily_extreme_timing(pd.DataFrame({"flow_cfs": values}, index=idx), tz="UTC")
        assert timing["hour_of_max"].iloc[0] == pytest.approx(6.0)  # the first 50
        assert timing["n_tied_max"].iloc[0] == 2
        assert timing["n_tied_min"].iloc[0] == 1

    def test_fractional_hour_keeps_sub_hour_resolution(self) -> None:
        """Rounding a 15-minute record to the hour throws away three quarters
        of the timing resolution it carries."""
        idx = pd.date_range("2021-08-01 00:00", periods=96, freq="15min", tz="UTC")
        values = np.full(96, 10.0)
        values[57] = 99.0  # 14:15
        timing = daily_extreme_timing(pd.DataFrame({"flow_cfs": values}, index=idx), tz="UTC")
        assert timing["hour_of_max"].iloc[0] == pytest.approx(14.25)

    def test_negative_values_treated_as_missing(self) -> None:
        idx = pd.date_range("2021-06-01", periods=8, freq="3h", tz="UTC")
        values = np.array([50.0, 55.0, -999.0, 60.0, 65.0, 45.0, 40.0, 52.0])
        timing = daily_extreme_timing(pd.DataFrame({"flow_cfs": values}, index=idx), tz="UTC")
        assert timing["n_obs"].iloc[0] == 7
        assert timing["min_value"].iloc[0] == 40.0

    def test_expected_obs_uses_actual_local_day_length(self) -> None:
        """Unlike diel_variation's fixed 1440 minutes: a 25-hour fall-back day
        expects 100 quarter-hours, not 96, so it is not mismarked."""
        idx = pd.date_range("2020-10-31", "2020-11-02 23:45", freq="15min", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": 100 + np.arange(len(idx)) * 0.1}, index=idx)
        timing = daily_extreme_timing(iv, tz=PACIFIC).set_index("date")
        assert timing.loc[pd.Timestamp("2020-11-01").date(), "expected_obs"] == pytest.approx(100.0)
        assert timing.loc[pd.Timestamp("2020-10-31").date(), "expected_obs"] == pytest.approx(96.0)

    def test_spring_forward_day_expects_fewer_observations(self) -> None:
        idx = pd.date_range("2021-03-13", "2021-03-15 23:45", freq="15min", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": 100 + np.arange(len(idx)) * 0.1}, index=idx)
        timing = daily_extreme_timing(iv, tz=PACIFIC).set_index("date")
        assert timing.loc[pd.Timestamp("2021-03-14").date(), "expected_obs"] == pytest.approx(92.0)

    def test_incomplete_day_still_reports_timing(self) -> None:
        """Matching diel_variation: a gappy day usually still catches its peak,
        so the value is reported and `complete` grades it."""
        idx = pd.date_range("2021-06-01 06:00", periods=4, freq="1h", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": [10.0, 30.0, 20.0, 15.0]}, index=idx)
        timing = daily_extreme_timing(iv, tz="UTC")
        assert timing["hour_of_max"].iloc[0] == pytest.approx(7.0)
        assert not timing["complete"].iloc[0]

    def test_works_on_a_stage_series(self) -> None:
        idx = pd.date_range("2021-06-01", periods=96, freq="15min", tz="UTC")
        stage = pd.DataFrame({"gage_height_ft": 4.0 + np.sin(np.arange(96) / 15)}, index=idx)
        timing = daily_extreme_timing(stage, tz="UTC", value_col="gage_height_ft")
        assert timing["max_value"].iloc[0] > timing["min_value"].iloc[0]

    def test_too_few_timestamps_raises(self) -> None:
        iv = pd.DataFrame(
            {"flow_cfs": [50.0]}, index=pd.date_range("2020-01-01", periods=1, tz="UTC")
        )
        with pytest.raises(ValueError, match="at least 2 timestamps"):
            daily_extreme_timing(iv, tz="UTC")

    def test_missing_value_column_raises_naming_it(self) -> None:
        iv = pd.DataFrame(
            {"flow_cfs": [1.0, 2.0]},
            index=pd.date_range("2020-01-01", periods=2, freq="1h", tz="UTC"),
        )
        with pytest.raises(KeyError, match="gage_height_ft"):
            daily_extreme_timing(iv, tz="UTC", value_col="gage_height_ft")

    def test_naive_index_raises(self) -> None:
        iv = pd.DataFrame(
            {"flow_cfs": 50.0}, index=pd.date_range("2020-01-01", periods=10, freq="1h")
        )
        with pytest.raises(TypeError, match="timezone-aware"):
            daily_extreme_timing(iv, tz="UTC")


class TestExtremeTimingSummary:
    """Tests for extreme_timing_summary."""

    def test_summary_matches_known_phase(self) -> None:
        timing = daily_extreme_timing(_pacific_diel_series(n_days=11), tz=PACIFIC)
        summary = extreme_timing_summary(timing, which="max")
        assert summary["mean_hour_of_max"] == pytest.approx(16.0, abs=0.3)
        assert summary["concentration_of_max"] > 0.99
        assert summary["rayleigh_p_of_max"] < 0.01

    def test_min_summary_recovers_the_other_extreme(self) -> None:
        """The pair is the point: a morning minimum with an afternoon maximum
        is the snowmelt signature, and either half alone is ambiguous."""
        timing = daily_extreme_timing(_pacific_diel_series(n_days=11), tz=PACIFIC)
        summary = extreme_timing_summary(timing, which="min")
        assert summary["mean_hour_of_min"] == pytest.approx(4.0, abs=0.3)

    def test_only_complete_days_pooled(self) -> None:
        timing = daily_extreme_timing(_pacific_diel_series(n_days=5), tz=PACIFIC)
        assert extreme_timing_summary(timing)["n_days"] == int(timing["complete"].sum())

    def test_flat_days_excluded_from_n_days(self) -> None:
        """A flat day is complete but has no defined timing, so n_days counts
        fewer rows than `complete` does."""
        timing = daily_extreme_timing(_pacific_diel_series(n_days=5), tz=PACIFIC)
        timing.loc[2, ["hour_of_max", "hour_of_min"]] = np.nan
        assert extreme_timing_summary(timing)["n_days"] == int(timing["complete"].sum()) - 1

    def test_all_incomplete_gives_nan_not_error(self) -> None:
        timing = daily_extreme_timing(_pacific_diel_series(n_days=5), tz=PACIFIC)
        timing["complete"] = False
        summary = extreme_timing_summary(timing)
        assert summary["n_days"] == 0
        assert np.isnan(summary["mean_hour_of_max"])

    def test_bad_which_raises(self) -> None:
        timing = daily_extreme_timing(_pacific_diel_series(n_days=3), tz=PACIFIC)
        with pytest.raises(ValueError, match="must be 'max' or 'min'"):
            extreme_timing_summary(timing, which="peak")


class TestRampingRates:
    """Tests for ramping_rates."""

    def test_known_linear_ramp_recovered_exactly(self) -> None:
        """A series rising 10 cfs every 15 minutes ramps at 40 cfs/hr."""
        idx = pd.date_range("2021-06-01", periods=96, freq="15min", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": 100.0 + 10.0 * np.arange(96)}, index=idx)
        ramps = ramping_rates(iv, tz="UTC")
        assert ramps["max_up_ramp_cfs_per_hr"].iloc[0] == pytest.approx(40.0)
        assert ramps["mean_abs_ramp_cfs_per_hr"].iloc[0] == pytest.approx(40.0)

    def test_down_ramp_is_signed_and_never_positive(self) -> None:
        """A magnitude compared against a limit written as a negative number is
        a bug waiting to happen; a value that cannot be positive is not."""
        idx = pd.date_range("2021-06-01", periods=96, freq="15min", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": 1000.0 - 5.0 * np.arange(96)}, index=idx)
        ramps = ramping_rates(iv, tz="UTC")
        assert ramps["max_down_ramp_cfs_per_hr"].iloc[0] == pytest.approx(-20.0)
        assert (ramps["max_down_ramp_cfs_per_hr"] <= 0).all()

    def test_monotone_fall_gives_zero_up_ramp_not_nan(self) -> None:
        """0.0 means "measured, and nothing rose"; NaN is reserved for "no
        valid interval at all". A recession is the first, not the second."""
        idx = pd.date_range("2021-06-01", periods=96, freq="15min", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": 1000.0 - 5.0 * np.arange(96)}, index=idx)
        ramps = ramping_rates(iv, tz="UTC")
        assert ramps["max_up_ramp_cfs_per_hr"].iloc[0] == 0.0
        assert ramps["n_intervals"].iloc[0] > 0

    def test_dst_fall_back_day_has_no_negative_intervals(self) -> None:
        """The assertion that catches a dt measured on local wall-clock time.

        Across the autumn transition the local clock runs 01:59 -> 01:00, so a
        local-time difference is negative for a real forward hour: a 15-minute
        interval spanning it computes as -45 minutes, turning a modest rise into
        a large fall. The input here rises monotonically, so any negative
        down-ramp is that bug.
        """
        idx = pd.date_range("2020-10-31", "2020-11-02 23:45", freq="15min", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": 100.0 + np.arange(len(idx)) * 0.1}, index=idx)
        ramps = ramping_rates(iv, tz=PACIFIC)
        assert (ramps["max_down_ramp_cfs_per_hr"] == 0.0).all()
        assert (ramps["max_up_ramp_cfs_per_hr"] > 0).all()

    def test_gap_excluded_and_counted(self) -> None:
        """A change across a hole is not a ramping rate: it reads as a gentle
        ramp however violently flow moved inside the gap."""
        first = pd.date_range("2021-06-01 00:00", periods=20, freq="15min", tz="UTC")
        second = pd.date_range("2021-06-01 09:00", periods=20, freq="15min", tz="UTC")
        values = np.concatenate([np.full(20, 100.0), np.full(20, 500.0)])
        iv = pd.DataFrame({"flow_cfs": values}, index=first.append(second))
        ramps = ramping_rates(iv, tz="UTC")
        assert ramps["n_intervals_gapped"].iloc[0] == 1
        assert ramps["n_intervals"].iloc[0] == 38
        # The 400-cfs jump is excluded, so nothing ramped at all.
        assert ramps["max_up_ramp_cfs_per_hr"].iloc[0] == 0.0
        assert ramps["mean_abs_ramp_cfs_per_hr"].iloc[0] == pytest.approx(0.0)

    def test_explicit_max_gap_hours_admits_the_gap(self) -> None:
        """The same data with a wide enough threshold does include the jump,
        confirming the exclusion above is the threshold's doing."""
        first = pd.date_range("2021-06-01 00:00", periods=20, freq="15min", tz="UTC")
        second = pd.date_range("2021-06-01 09:00", periods=20, freq="15min", tz="UTC")
        values = np.concatenate([np.full(20, 100.0), np.full(20, 500.0)])
        iv = pd.DataFrame({"flow_cfs": values}, index=first.append(second))
        ramps = ramping_rates(iv, tz="UTC", max_gap_hours=6.0)
        assert ramps["n_intervals_gapped"].iloc[0] == 0
        assert ramps["max_up_ramp_cfs_per_hr"].iloc[0] == pytest.approx(400.0 / 4.25)

    def test_default_gap_threshold_scales_with_the_record(self) -> None:
        """An hourly record must not have every interval judged a gap by a
        threshold calibrated for 15-minute data."""
        idx = pd.date_range("2021-06-01", periods=48, freq="1h", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": 100.0 + np.arange(48)}, index=idx)
        ramps = ramping_rates(iv, tz="UTC")
        assert (ramps["n_intervals_gapped"] == 0).all()
        assert DEFAULT_MAX_GAP_MULTIPLE == 3.0

    def test_interval_label_shifts_midnight_straddling_intervals(self) -> None:
        """The convention is a real choice, so the two settings must differ."""
        idx = pd.date_range("2021-06-01 23:00", periods=8, freq="15min", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": 100.0 + 10.0 * np.arange(8)}, index=idx)
        by_end = ramping_rates(iv, tz="UTC", interval_label="end").set_index("date")
        by_start = ramping_rates(iv, tz="UTC", interval_label="start").set_index("date")
        june1 = pd.Timestamp("2021-06-01").date()
        assert by_end.loc[june1, "n_intervals"] != by_start.loc[june1, "n_intervals"]

    def test_percent_rate_relative_to_preceding_value(self) -> None:
        idx = pd.date_range("2021-06-01", periods=5, freq="1h", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": [100.0, 110.0, 121.0, 133.1, 146.41]}, index=idx)
        ramps = ramping_rates(iv, tz="UTC")
        assert ramps["max_up_ramp_pct_per_hr"].iloc[0] == pytest.approx(10.0)

    def test_percent_rate_nan_from_a_zero_preceding_value(self) -> None:
        """A percentage change from zero is not a number, and an intermittent
        reach produces these. The absolute rate is still reported.

        The NaN is per *interval*, so a day's percent rate is only NaN when
        every one of its intervals started from zero -- which is the case here.
        """
        idx = pd.date_range("2021-06-01", periods=4, freq="1h", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": [0.0, 0.0, 0.0, 5.0]}, index=idx)
        ramps = ramping_rates(iv, tz="UTC")
        assert np.isnan(ramps["max_up_ramp_pct_per_hr"].iloc[0])
        assert ramps["max_up_ramp_cfs_per_hr"].iloc[0] == pytest.approx(5.0)

    def test_one_defined_percent_interval_is_enough_for_a_value(self) -> None:
        """The other side of the rule above: a day with a single interval whose
        preceding value was positive gets a percent rate from that interval,
        and 0.0 for the direction nothing moved in."""
        idx = pd.date_range("2021-06-01", periods=4, freq="1h", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": [0.0, 5.0, 0.0, 0.0]}, index=idx)
        ramps = ramping_rates(iv, tz="UTC")
        assert ramps["max_down_ramp_pct_per_hr"].iloc[0] == pytest.approx(-100.0)
        assert ramps["max_up_ramp_pct_per_hr"].iloc[0] == 0.0

    def test_stage_series_gets_ft_columns_and_no_percent(self) -> None:
        idx = pd.date_range("2021-06-01", periods=96, freq="15min", tz="UTC")
        stage = pd.DataFrame({"gage_height_ft": 4.0 + 0.01 * np.arange(96)}, index=idx)
        ramps = ramping_rates(stage, tz="UTC", value_col="gage_height_ft")
        assert "max_up_ramp_ft_per_hr" in ramps.columns
        assert not [c for c in ramps.columns if "pct" in c]
        assert ramps["max_up_ramp_ft_per_hr"].iloc[0] == pytest.approx(0.04)

    def test_relative_on_stage_raises_naming_the_datum(self) -> None:
        """Gage height is measured from an arbitrary local datum, so a percent
        change in it has no referent while looking exactly like a regulatory
        ramping rate. Refused rather than returned."""
        idx = pd.date_range("2021-06-01", periods=96, freq="15min", tz="UTC")
        stage = pd.DataFrame({"gage_height_ft": 4.0 + 0.01 * np.arange(96)}, index=idx)
        with pytest.raises(ValueError, match="datum"):
            ramping_rates(stage, tz="UTC", value_col="gage_height_ft", relative=True)

    def test_reversals_counted(self) -> None:
        """Sign changes in the rate sequence -- the standard hydropeaking count."""
        idx = pd.date_range("2021-06-01", periods=7, freq="1h", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": [10.0, 20.0, 10.0, 20.0, 10.0, 20.0, 10.0]}, index=idx)
        assert ramping_rates(iv, tz="UTC")["n_reversals"].iloc[0] == 5

    def test_negative_values_do_not_create_intervals(self) -> None:
        """A reading masked as an artifact leaves a hole no different from an
        outage; the pair spans it and is then judged against the gap rule."""
        idx = pd.date_range("2021-06-01", periods=6, freq="1h", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": [100.0, 110.0, -999.0, 130.0, 140.0, 150.0]}, index=idx)
        ramps = ramping_rates(iv, tz="UTC")
        assert ramps["n_obs"].iloc[0] == 5
        assert ramps["n_intervals"].iloc[0] + ramps["n_intervals_gapped"].iloc[0] == 4

    def test_unknown_value_column_units_must_be_given(self) -> None:
        idx = pd.date_range("2021-06-01", periods=10, freq="1h", tz="UTC")
        iv = pd.DataFrame({"temp_c": 10.0 + np.arange(10)}, index=idx)
        with pytest.raises(ValueError, match="Cannot infer units"):
            ramping_rates(iv, tz="UTC", value_col="temp_c")
        ramps = ramping_rates(iv, tz="UTC", value_col="temp_c", units="degC")
        assert "max_up_ramp_degC_per_hr" in ramps.columns

    def test_bad_interval_label_raises(self) -> None:
        idx = pd.date_range("2021-06-01", periods=10, freq="1h", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": np.arange(10.0)}, index=idx)
        with pytest.raises(ValueError, match="interval_label must be one of"):
            ramping_rates(iv, tz="UTC", interval_label="middle")

    def test_non_positive_max_gap_hours_raises(self) -> None:
        idx = pd.date_range("2021-06-01", periods=10, freq="1h", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": np.arange(10.0)}, index=idx)
        with pytest.raises(ValueError, match="max_gap_hours must be positive"):
            ramping_rates(iv, tz="UTC", max_gap_hours=0.0)

    def test_too_few_timestamps_raises(self) -> None:
        iv = pd.DataFrame(
            {"flow_cfs": [50.0]}, index=pd.date_range("2020-01-01", periods=1, tz="UTC")
        )
        with pytest.raises(ValueError, match="at least 2 timestamps"):
            ramping_rates(iv, tz="UTC")


class TestRampingRateSummary:
    """Tests for ramping_rate_summary."""

    def test_summary_over_a_known_ramp(self) -> None:
        idx = pd.date_range("2021-06-01", periods=96 * 3, freq="15min", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": 500.0 + 20.0 * np.sin(np.arange(len(idx)) / 24)}, index=idx)
        summary = ramping_rate_summary(ramping_rates(iv, tz="UTC"))
        assert summary["n_days"] > 0
        assert summary["extreme_max_up_ramp_cfs_per_hr"] > 0
        assert summary["extreme_max_down_ramp_cfs_per_hr"] < 0

    def test_extreme_is_the_most_severe_not_the_largest(self) -> None:
        """For a down-ramp column the extreme is the minimum, so it is the
        worst value in the period rather than the numerically largest."""
        idx = pd.date_range("2021-06-01", periods=96 * 3, freq="15min", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": 500.0 + 20.0 * np.sin(np.arange(len(idx)) / 24)}, index=idx)
        ramps = ramping_rates(iv, tz="UTC")
        summary = ramping_rate_summary(ramps)
        complete = ramps[ramps["complete"]]
        assert summary["extreme_max_down_ramp_cfs_per_hr"] == pytest.approx(
            complete["max_down_ramp_cfs_per_hr"].min()
        )

    def test_negative_limit_counts_days_at_or_below_it(self) -> None:
        idx = pd.date_range("2021-06-01", periods=96 * 3, freq="15min", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": 500.0 - 1.0 * np.arange(len(idx))}, index=idx)
        ramps = ramping_rates(iv, tz="UTC")
        summary = ramping_rate_summary(ramps, limits={"max_down_ramp_cfs_per_hr": -1.0})
        assert summary["n_days_exceeding_max_down_ramp_cfs_per_hr"] == summary["n_days"]
        assert summary["frac_days_exceeding_max_down_ramp_cfs_per_hr"] == pytest.approx(1.0)

    def test_positive_limit_counts_days_at_or_above_it(self) -> None:
        idx = pd.date_range("2021-06-01", periods=96 * 3, freq="15min", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": 500.0 + 1.0 * np.arange(len(idx))}, index=idx)
        ramps = ramping_rates(iv, tz="UTC")
        summary = ramping_rate_summary(ramps, limits={"max_up_ramp_cfs_per_hr": 1e9})
        assert summary["n_days_exceeding_max_up_ramp_cfs_per_hr"] == 0

    def test_only_complete_days_pooled(self) -> None:
        idx = pd.date_range("2021-06-01 06:00", periods=96 * 2, freq="15min", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": 500.0 + np.arange(len(idx))}, index=idx)
        ramps = ramping_rates(iv, tz="UTC")
        assert ramping_rate_summary(ramps)["n_days"] == int(ramps["complete"].sum())

    def test_unknown_limit_column_raises(self) -> None:
        """Reporting zero exceedances for a limit that was never checked would
        be worse than raising."""
        idx = pd.date_range("2021-06-01", periods=96, freq="15min", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": 500.0 + np.arange(96)}, index=idx)
        ramps = ramping_rates(iv, tz="UTC")
        with pytest.raises(KeyError, match="max_down_ramp_ft_per_hr"):
            ramping_rate_summary(ramps, limits={"max_down_ramp_ft_per_hr": -0.5})


class TestModuleHelpers:
    """Tests for the private helpers that carry their own rules."""

    def test_reversals_skip_zero_rate_intervals(self) -> None:
        """A flat spell between a rise and a fall is one reversal, not two."""
        assert _count_reversals(pd.Series([5.0, 0.0, 0.0, -5.0])) == 1
        assert _count_reversals(pd.Series([5.0, -5.0, 5.0])) == 2
        assert _count_reversals(pd.Series([5.0])) == 0
        assert _count_reversals(pd.Series([np.nan, np.nan], dtype=float)) == 0

    def test_local_day_hours_handles_dst(self) -> None:
        dates = [pd.Timestamp(d).date() for d in ("2020-10-31", "2020-11-01", "2021-03-14")]
        hours = _local_day_hours(dates, PACIFIC)
        assert list(hours) == [24.0, 25.0, 23.0]

    def test_resolve_units_known_and_override(self) -> None:
        assert _resolve_units("flow_cfs", None) == "cfs"
        assert _resolve_units("gage_height_ft", None) == "ft"
        assert _resolve_units("anything", "m3s") == "m3s"
        assert KNOWN_VALUE_UNITS["flow_cfs"] == "cfs"
        assert "ft" not in TRUE_ZERO_UNITS


class TestDielVariationStillWorksFromHere:
    """diel_variation moved into this module; regime.py re-exports it.

    tests/test_regime.py covers its behaviour and was deliberately left
    importing it from flowfreq.regime, which is the back-compatibility test.
    This is only the direct-import check.
    """

    def test_importable_from_subdaily(self) -> None:
        daily = diel_variation(_pacific_diel_series(n_days=5), tz=PACIFIC)
        assert daily["range_cfs"].iloc[2] == pytest.approx(60.0, abs=1.0)

    def test_re_exported_from_regime_and_top_level(self) -> None:
        import flowfreq
        from flowfreq.regime import diel_variation as from_regime

        assert from_regime is diel_variation
        assert flowfreq.diel_variation is diel_variation
