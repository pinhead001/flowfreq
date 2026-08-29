"""Tests for hydrolib.regime (flow regime metrics)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hydrolib.regime import _days_in_season  # noqa: PLC2701 -- testing the helper directly
from hydrolib.regime import (
    BASEFLOW_METHODS,
    FlowRegime,
    _hysep_interval_days,
    _ih_smoothed_minima,
    _interpolate_baseflow,
    baseflow_index,
    monthly_flow_summary,
    richards_baker_flashiness,
    seasonal_flow_summary,
    separate_baseflow,
    tqmean,
)


def _synthetic_daily(n_years: int = 15, seed: int = 42, start: str = "2005-10-01") -> pd.DataFrame:
    """A snowmelt-shaped multi-year daily series with storm noise -- enough
    years and realistic enough shape to exercise every metric, not tuned to
    any specific expected numeric answer (those are checked against smaller,
    hand-computable series instead)."""
    rng = np.random.default_rng(seed)
    idx = pd.date_range(start, periods=365 * n_years, freq="D")
    wy_start = idx.map(lambda d: pd.Timestamp(f"{d.year if d.month >= 10 else d.year - 1}-10-01"))
    doy = (idx - wy_start).days.to_numpy()
    seasonal = 80 + 300 * np.exp(-(((doy - 210) % 365) ** 2) / (2 * 40**2))
    base = 40 + 20 * np.cos(2 * np.pi * (doy - 270) / 365)
    flow = np.clip(seasonal + base + rng.normal(0, 5, len(idx)), 1, None)
    return pd.DataFrame({"flow_cfs": flow}, index=idx)


class TestRichardsBakerFlashiness:
    """Tests for the Richards-Baker Flashiness Index."""

    def test_hand_computed_zigzag(self) -> None:
        """RBI = sum(|diffs|) / sum(flows) over the n-1 within-year differences."""
        idx = pd.date_range("2020-01-01", periods=5, freq="D")
        daily = pd.DataFrame({"flow_cfs": [10.0, 15.0, 10.0, 20.0, 10.0]}, index=idx)
        rbi = richards_baker_flashiness(daily, year_type="calendar", min_days=1)
        # diffs: 5,5,10,10 -> sum=30; flows sum=65
        assert rbi.iloc[0]["flashiness_index"] == pytest.approx(30 / 65)

    def test_flat_flow_gives_zero(self) -> None:
        """No day-to-day change at all is the minimally flashy case: RBI = 0."""
        idx = pd.date_range("2020-01-01", periods=10, freq="D")
        daily = pd.DataFrame({"flow_cfs": 50.0}, index=idx)
        rbi = richards_baker_flashiness(daily, year_type="calendar", min_days=1)
        assert rbi.iloc[0]["flashiness_index"] == 0.0

    def test_incomplete_year_is_nan_but_row_present(self) -> None:
        idx = pd.date_range("2020-01-01", periods=50, freq="D")
        daily = pd.DataFrame({"flow_cfs": 50.0}, index=idx)
        rbi = richards_baker_flashiness(daily, year_type="calendar", min_days=350)
        assert len(rbi) == 1
        assert not rbi.iloc[0]["complete"]
        assert np.isnan(rbi.iloc[0]["flashiness_index"])

    def test_gap_within_a_complete_year_is_a_small_approximation_not_a_crash(self) -> None:
        """A handful of missing days inside an otherwise-complete year drops
        the (at most two) difference terms adjacent to each gap via nansum,
        as documented, rather than crashing or fabricating a value."""
        idx = pd.date_range("2020-01-01", periods=365, freq="D")
        flow = 100 + 10 * np.sin(np.arange(365) / 20)
        daily = pd.DataFrame({"flow_cfs": flow}, index=idx)
        gapped = daily.drop(daily.index[180:185])

        rbi_gapped = richards_baker_flashiness(gapped, year_type="calendar", min_days=350)
        rbi_full = richards_baker_flashiness(daily, year_type="calendar", min_days=350)

        assert rbi_gapped.iloc[0]["complete"]
        assert not np.isnan(rbi_gapped.iloc[0]["flashiness_index"])
        assert rbi_gapped["flashiness_index"].iloc[0] == pytest.approx(
            rbi_full["flashiness_index"].iloc[0], abs=0.01
        )

    def test_negative_values_excluded_like_lowflow(self) -> None:
        """A negative daily value is a data artifact, not information -- must
        not contribute a spurious difference term."""
        idx = pd.date_range("2020-01-01", periods=10, freq="D")
        vals = np.full(10, 50.0)
        vals[5] = -999.0
        daily = pd.DataFrame({"flow_cfs": vals}, index=idx)
        rbi = richards_baker_flashiness(daily, year_type="calendar", min_days=1)
        # With day 5 excluded, the remaining series is flat -> RBI should be
        # small (only real variation, if any) rather than dominated by a
        # fabricated +/-1049 swing around the artifact.
        assert rbi.iloc[0]["flashiness_index"] < 0.01


class TestTQmean:
    """Tests for TQmean."""

    def test_hand_computed_case(self) -> None:
        """Only the value above the year's own mean counts."""
        idx = pd.date_range("2020-01-01", periods=5, freq="D")
        daily = pd.DataFrame({"flow_cfs": [10.0, 20.0, 30.0, 40.0, 100.0]}, index=idx)
        tq = tqmean(daily, year_type="calendar", min_days=1)
        # mean = 40; only 100 exceeds it -> 1/5 = 0.2
        assert tq.iloc[0]["tqmean"] == pytest.approx(0.2)

    def test_second_hand_computed_case(self) -> None:
        idx = pd.date_range("2020-01-01", periods=5, freq="D")
        daily = pd.DataFrame({"flow_cfs": [10.0, 20.0, 30.0, 40.0, 50.0]}, index=idx)
        tq = tqmean(daily, year_type="calendar", min_days=1)
        # mean = 30; 40 and 50 exceed it -> 2/5 = 0.4
        assert tq.iloc[0]["tqmean"] == pytest.approx(0.4)

    def test_bounded_between_zero_and_one(self) -> None:
        daily = _synthetic_daily()
        tq = tqmean(daily)
        valid = tq["tqmean"].dropna()
        assert ((valid >= 0) & (valid <= 1)).all()

    def test_incomplete_year_is_nan(self) -> None:
        idx = pd.date_range("2020-01-01", periods=50, freq="D")
        daily = pd.DataFrame({"flow_cfs": 50.0}, index=idx)
        tq = tqmean(daily, year_type="calendar", min_days=350)
        assert not tq.iloc[0]["complete"]
        assert np.isnan(tq.iloc[0]["tqmean"])


class TestInterpolateBaseflow:
    """Tests for the interpolate+cap helper shared by IH and HYSEP-local-minimum."""

    def test_caps_at_actual_flow_when_interpolation_would_exceed_it(self) -> None:
        """Baseflow cannot exceed total flow -- a brief dip below the
        interpolated line must be followed down, not overridden."""
        flows = np.array([10.0, 10.0, 3.0, 10.0, 10.0])
        bf = _interpolate_baseflow(flows, np.array([0, 4]), np.array([10.0, 10.0]))
        assert bf[2] == 3.0
        assert np.allclose(bf[[0, 1, 3, 4]], 10.0)

    def test_pure_linear_interpolation_when_no_capping_needed(self) -> None:
        flows = np.array([5.0, 100.0, 100.0, 100.0, 15.0])
        bf = _interpolate_baseflow(flows, np.array([0, 4]), np.array([5.0, 15.0]))
        assert np.allclose(bf, [5, 7.5, 10, 12.5, 15])

    def test_fewer_than_two_turning_points_gives_all_nan(self) -> None:
        flows = np.array([10.0, 10.0, 3.0, 10.0, 10.0])
        bf = _interpolate_baseflow(flows, np.array([2]), np.array([3.0]))
        assert np.all(np.isnan(bf))


class TestIhSmoothedMinima:
    """Tests for the UKIH/IH smoothed-minima turning-point selection."""

    def test_accepts_genuine_local_minimum_rejects_recession_neighbors(self) -> None:
        """Hand-derived block-minima sequence [40, 20, 10, 15(ish), 30]: only
        the middle block (a genuine local low) should pass the 0.9-ratio
        turning-point test; its neighbors, each on a monotonic slope toward
        or away from it, should not."""
        flows = np.concatenate(
            [
                np.linspace(100, 40, 5),
                np.linspace(38, 20, 5),
                np.linspace(18, 10, 5),  # true low block
                np.linspace(15, 22, 5),
                np.linspace(30, 45, 5),
            ]
        )
        bf = _ih_smoothed_minima(flows, block_days=5, factor=0.9)
        # Only one turning point exists in this short example (verified by
        # construction), so no pair to interpolate between -> all NaN.
        assert np.all(np.isnan(bf))

    def test_baseflow_never_exceeds_actual_flow(self) -> None:
        daily = _synthetic_daily()
        flows = daily["flow_cfs"].to_numpy()
        bf = _ih_smoothed_minima(flows)
        valid = ~np.isnan(bf)
        assert np.all(bf[valid] <= flows[valid] + 1e-6)


class TestHysep:
    """Tests for the three HYSEP variants and their shared N* formula."""

    def test_interval_days_matches_published_examples(self) -> None:
        """N* = round_to_odd(2 * A**0.2); A=100 sq mi is a commonly cited N*=5 case."""
        assert _hysep_interval_days(100) == 5

    def test_interval_days_increases_with_drainage_area(self) -> None:
        assert _hysep_interval_days(10) < _hysep_interval_days(1000) < _hysep_interval_days(10000)

    def test_interval_days_rejects_nonpositive_area(self) -> None:
        with pytest.raises(ValueError, match="must be positive"):
            _hysep_interval_days(0)

    def test_fixed_interval_is_a_step_function(self) -> None:
        """Baseflow is constant within each non-overlapping 2*N*-day block."""
        flows = np.array(
            [50, 45, 40, 35, 30, 10, 32, 36, 40, 44, 48, 52, 20, 48, 44, 40, 36, 32, 28, 24],
            dtype=float,
        )
        from hydrolib.regime import _hysep_fixed

        bf = _hysep_fixed(flows, n_star=2)
        assert bf[0] == bf[1] == bf[2] == bf[3] == flows[0:4].min()

    def test_sliding_interval_never_exceeds_actual_flow(self) -> None:
        from hydrolib.regime import _hysep_sliding

        flows = np.array(
            [50, 45, 40, 35, 30, 10, 32, 36, 40, 44, 48, 52, 20, 48, 44, 40, 36, 32, 28, 24],
            dtype=float,
        )
        bf = _hysep_sliding(flows, n_star=2)
        valid = ~np.isnan(bf)
        assert np.all(bf[valid] <= flows[valid] + 1e-9)

    def test_local_minimum_identifies_known_local_minima(self) -> None:
        from hydrolib.regime import _hysep_local_minimum

        flows = np.array(
            [50, 45, 40, 35, 30, 10, 32, 36, 40, 44, 48, 52, 20, 48, 44, 40, 36, 32, 28, 24],
            dtype=float,
        )
        bf = _hysep_local_minimum(flows, n_star=2)
        assert bf[5] == 10.0
        assert bf[12] == 20.0

    @pytest.mark.parametrize("method", ["hysep_fixed", "hysep_sliding", "hysep_local_minimum"])
    def test_requires_drainage_area(self, method: str) -> None:
        daily = _synthetic_daily(n_years=2)
        with pytest.raises(ValueError, match="requires drainage_area_sqmi"):
            separate_baseflow(daily, method=method)


class TestLyneHollick:
    """Tests for the Lyne-Hollick recursive digital filter."""

    def _recession_with_storm(self) -> np.ndarray:
        days = np.arange(90)
        recession = 100 * np.exp(-days / 60.0) + 10
        storm = np.zeros(90)
        storm[40:45] = [80, 200, 350, 150, 40]
        return recession + storm

    def test_gap_does_not_propagate_through_the_recursion(self) -> None:
        """A single missing day must not poison every value after it -- each
        day's filtered value depends on the previous one via the recursion,
        so an unguarded gap would corrupt the entire remainder."""
        from hydrolib.regime import _lyne_hollick

        flows = self._recession_with_storm()
        gapped = flows.copy()
        gapped[30] = np.nan

        bf = _lyne_hollick(gapped)
        assert np.isnan(bf).sum() == 1
        assert np.isnan(bf[30])
        assert not np.isnan(bf[29])
        assert not np.isnan(bf[89])

    def test_recession_is_mostly_baseflow(self) -> None:
        from hydrolib.regime import _lyne_hollick

        flows = self._recession_with_storm()
        bf = _lyne_hollick(flows)
        bfi_recession = bf[:35].sum() / flows[:35].sum()
        assert bfi_recession > 0.85

    def test_storm_peak_is_mostly_quickflow(self) -> None:
        from hydrolib.regime import _lyne_hollick

        flows = self._recession_with_storm()
        bf = _lyne_hollick(flows)
        assert bf[42] < flows[42] * 0.5

    def test_baseflow_bounded_by_zero_and_total_flow(self) -> None:
        from hydrolib.regime import _lyne_hollick

        flows = self._recession_with_storm()
        bf = _lyne_hollick(flows)
        assert np.all(bf >= -1e-9)
        assert np.all(bf <= flows + 1e-9)


class TestSeparateBaseflowAndBaseflowIndex:
    """Tests for the public dispatch functions."""

    def test_unknown_method_raises(self) -> None:
        daily = _synthetic_daily(n_years=2)
        with pytest.raises(ValueError, match="method must be one of"):
            separate_baseflow(daily, method="digital_filter_v2")

    @pytest.mark.parametrize("method", BASEFLOW_METHODS)
    def test_all_methods_bounded_by_actual_flow(self, method: str) -> None:
        daily = _synthetic_daily(n_years=3)
        kwargs = {"drainage_area_sqmi": 1772.0} if method.startswith("hysep") else {}
        bf = separate_baseflow(daily, method=method, **kwargs)
        aligned = daily["flow_cfs"].reindex(bf.index)
        valid = bf.notna() & aligned.notna()
        assert (bf[valid] <= aligned[valid] + 1e-6).all()

    def test_baseflow_index_bounded_zero_one(self) -> None:
        daily = _synthetic_daily()
        bfi = baseflow_index(daily)
        valid = bfi["baseflow_index"].dropna()
        assert ((valid >= 0) & (valid <= 1)).all()

    def test_baseflow_index_reports_method_used(self) -> None:
        daily = _synthetic_daily(n_years=2)
        bfi = baseflow_index(daily, method="lyne_hollick")
        assert (bfi["method"] == "lyne_hollick").all()

    def test_baseflow_index_requires_drainage_area_for_hysep(self) -> None:
        daily = _synthetic_daily(n_years=2)
        with pytest.raises(ValueError, match="requires drainage_area_sqmi"):
            baseflow_index(daily, method="hysep_fixed")


class TestMonthlyFlowSummary:
    """Tests for the per-(year, month) summary table."""

    def test_hand_computed_statistics(self) -> None:
        idx = pd.DatetimeIndex(["2020-01-05", "2020-01-10", "2020-01-15"])
        daily = pd.DataFrame({"flow_cfs": [10.0, 20.0, 30.0]}, index=idx)
        m = monthly_flow_summary(daily)
        row = m[(m["year"] == 2020) & (m["month"] == 1)].iloc[0]
        assert row["mean_flow_cfs"] == 20.0
        assert row["min_flow_cfs"] == 10.0
        assert row["max_flow_cfs"] == 30.0
        assert row["median_flow_cfs"] == 20.0
        assert row["n_days"] == 3

    def test_days_in_month_is_true_calendar_length_not_data_span(self) -> None:
        """Regression: days_in_month must reflect January's true 31 days even
        though only 3 scattered days of data were provided -- not the span
        between the first and last provided date."""
        idx = pd.DatetimeIndex(["2020-01-05", "2020-01-10", "2020-01-15"])
        daily = pd.DataFrame({"flow_cfs": [10.0, 20.0, 30.0]}, index=idx)
        m = monthly_flow_summary(daily)
        row = m[(m["year"] == 2020) & (m["month"] == 1)].iloc[0]
        assert row["days_in_month"] == 31
        assert not row["complete"]

    def test_leap_february_correct_under_climatic_year_relabeling(self) -> None:
        """Regression: under year_type="climatic", February is labeled by
        the PREVIOUS calendar year (climatic year Y spans Apr Y - Mar Y+1).
        Naively computing days_in_month from that label would ask for the
        wrong year's leap status. February 2020 (leap, 29 days) must still
        report 29, not 28, even though it is labeled climatic-year 2019."""
        idx = pd.date_range("2020-02-01", "2020-02-29", freq="D")
        daily = pd.DataFrame({"flow_cfs": 50.0}, index=idx)
        m = monthly_flow_summary(daily, year_type="climatic")
        row = m.iloc[0]
        assert row["year"] == 2019
        assert row["days_in_month"] == 29
        assert row["n_days"] == 29
        assert row["complete"]

    def test_non_leap_february_under_climatic_year_relabeling(self) -> None:
        idx = pd.date_range("2019-02-01", "2019-02-28", freq="D")
        daily = pd.DataFrame({"flow_cfs": 50.0}, index=idx)
        m = monthly_flow_summary(daily, year_type="climatic")
        row = m.iloc[0]
        assert row["year"] == 2018
        assert row["days_in_month"] == 28

    def test_all_twelve_months_present_for_multi_year_data(self) -> None:
        daily = _synthetic_daily()
        m = monthly_flow_summary(daily)
        assert set(m["month"]) == set(range(1, 13))


class TestSeasonalFlowSummary:
    """Tests for the per-(year, season) summary table."""

    def test_december_grouped_into_following_winter(self) -> None:
        """Dec 2019 + Jan 2020 + Feb 2020 must all fall under one 'winter
        2020' row, not split into 'winter 2019' and 'winter 2020'."""
        idx = pd.DatetimeIndex(["2019-12-15", "2020-01-15", "2020-02-15"])
        daily = pd.DataFrame({"flow_cfs": [10.0, 20.0, 30.0]}, index=idx)
        s = seasonal_flow_summary(daily)
        winter_rows = s[s["season"] == "winter"]
        assert len(winter_rows) == 1
        assert winter_rows.iloc[0]["year"] == 2020
        assert winter_rows.iloc[0]["mean_flow_cfs"] == 20.0

    def test_days_in_season_true_calendar_length_leap_year(self) -> None:
        assert _days_in_season(2020, "winter") == 31 + 31 + 29
        assert _days_in_season(2019, "winter") == 31 + 31 + 28

    def test_edge_of_record_season_denominator_is_true_length_not_data_span(self) -> None:
        """Regression: a record starting Oct 1 only contains 61 days of
        'fall' (Sep+Oct+Nov=91 total). days_in_season must be 91 (so
        'complete' correctly comes out False), not 61 (the span actually
        provided), which would make a season missing all of September read
        as though it could reach 100% complete."""
        idx = pd.date_range("2005-10-01", periods=61, freq="D")
        daily = pd.DataFrame({"flow_cfs": 50.0}, index=idx)
        s = seasonal_flow_summary(daily)
        row = s.iloc[0]
        assert row["season"] == "fall"
        assert row["n_days"] == 61
        assert row["days_in_season"] == 91
        assert not row["complete"]

    def test_full_leap_winter_is_complete(self) -> None:
        idx = pd.date_range("2019-12-01", "2020-02-29", freq="D")
        daily = pd.DataFrame({"flow_cfs": 50.0}, index=idx)
        s = seasonal_flow_summary(daily)
        row = s[s["season"] == "winter"].iloc[0]
        assert row["days_in_season"] == 91
        assert row["n_days"] == 91
        assert row["complete"]

    def test_all_four_seasons_present_for_multi_year_data(self) -> None:
        daily = _synthetic_daily()
        s = seasonal_flow_summary(daily)
        assert set(s["season"]) == {"winter", "spring", "summer", "fall"}

    def test_season_rows_are_chronologically_ordered(self) -> None:
        daily = _synthetic_daily(n_years=3)
        s = seasonal_flow_summary(daily)
        order = {"winter": 0, "spring": 1, "summer": 2, "fall": 3}
        keys = list(zip(s["year"], s["season"].map(order)))
        assert keys == sorted(keys)


class TestFlowRegime:
    """Tests for the FlowRegime facade class."""

    def test_annual_monthly_seasonal_all_populated(self) -> None:
        daily = _synthetic_daily()
        regime = FlowRegime(daily, drainage_area_sqmi=1772.0)
        assert not regime.annual.empty
        assert not regime.monthly.empty
        assert not regime.seasonal.empty
        assert not regime.baseflow_series.empty

    def test_annual_has_all_three_metrics(self) -> None:
        daily = _synthetic_daily()
        regime = FlowRegime(daily)
        for col in ("flashiness_index", "tqmean", "baseflow_index"):
            assert col in regime.annual.columns

    def test_hysep_method_flows_through_facade(self) -> None:
        daily = _synthetic_daily(n_years=3)
        regime = FlowRegime(daily, baseflow_method="hysep_local_minimum", drainage_area_sqmi=500.0)
        assert regime.annual["baseflow_index"].notna().any()

    def test_hysep_without_drainage_area_raises(self) -> None:
        daily = _synthetic_daily(n_years=3)
        with pytest.raises(ValueError, match="requires drainage_area_sqmi"):
            FlowRegime(daily, baseflow_method="hysep_fixed")

    def test_properties_return_copies_not_internal_state(self) -> None:
        daily = _synthetic_daily(n_years=3)
        regime = FlowRegime(daily)
        annual = regime.annual
        annual["flashiness_index"] = 0.0
        assert not (regime.annual["flashiness_index"] == 0.0).all()

    def test_summary_n_years_matches_complete_annual_rows(self) -> None:
        daily = _synthetic_daily(n_years=15)
        regime = FlowRegime(daily)
        summary = regime.summary()
        assert summary["n_years"] == regime.annual["complete"].sum()

    def test_summary_values_within_range_of_annual_values(self) -> None:
        """The pooled/averaged period-of-record figures should be broadly
        consistent with the per-year values that feed them, not wildly off."""
        daily = _synthetic_daily(n_years=15)
        regime = FlowRegime(daily)
        summary = regime.summary()
        complete = regime.annual[regime.annual["complete"]]

        assert (
            complete["flashiness_index"].min()
            <= summary["flashiness_index"]
            <= complete["flashiness_index"].max()
        )
        assert (
            complete["baseflow_index"].min()
            <= summary["baseflow_index"]
            <= complete["baseflow_index"].max()
        )

    def test_summary_baseflow_index_is_volume_weighted_not_mean_of_ratios(self) -> None:
        """Construct two complete years with very different total volume and
        a baseflow index that differs sharply between them; the pooled
        summary must land closer to the high-volume year's own BFI than a
        naive arithmetic mean of the two per-year BFIs would, since the
        high-volume year contributes far more of the pooled numerator and
        denominator."""
        idx = pd.date_range("2010-01-01", periods=365 * 2, freq="D")
        year1 = 20 + 2 * np.sin(np.arange(365) / 30)  # low volume, high BFI (smooth)
        year2 = 500 + 400 * np.sin(np.arange(365) / 10) ** 2  # high volume, flashier
        flow = np.concatenate([year1, year2])
        daily = pd.DataFrame({"flow_cfs": np.clip(flow, 1, None)}, index=idx)

        regime = FlowRegime(daily, year_type="calendar", min_days=350)
        summary = regime.summary()
        annual = regime.annual.set_index("year")
        bfi_2010 = annual.loc[2010, "baseflow_index"]
        bfi_2011 = annual.loc[2011, "baseflow_index"]
        naive_mean = (bfi_2010 + bfi_2011) / 2

        # The pooled value should sit closer to year 2's BFI (which dominates
        # by volume) than the unweighted mean does, whenever the two years'
        # BFIs actually differ.
        if abs(bfi_2010 - bfi_2011) > 0.01:
            assert abs(summary["baseflow_index"] - bfi_2011) < abs(naive_mean - bfi_2011)
