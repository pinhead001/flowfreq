"""Tests for flowfreq.regime (flow regime metrics)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from flowfreq.regime import _days_in_season  # noqa: PLC2701 -- testing the helper directly
from flowfreq.regime import (
    BASEFLOW_METHODS,
    FlowRegime,
    _hysep_interval_days,
    _ih_smoothed_minima,
    _interpolate_baseflow,
    _lyne_hollick,
    _water_year_blocks,
    baseflow_index,
    diel_variation,
    diel_variation_summary,
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


class TestGapsInTheRecord:
    """Interpolated baseflow does not cross a gap too long for a complete year."""

    METHODS = (
        "ih_smoothed_minima",
        "ih_smoothed_minima_log",
        "ih_smoothed_minima_log_wy",
        "hysep_local_minimum",
    )

    @pytest.mark.parametrize("method", METHODS)
    def test_a_long_gap_separates_as_two_records(self, method: str) -> None:
        daily = _synthetic_daily(n_years=6, start="1990-10-01").iloc[:2000]
        early, late = daily.iloc[:800], daily.iloc[1500:]
        whole = separate_baseflow(pd.concat([early, late]), method, 42.0)
        for part in (early, late):
            alone = separate_baseflow(part, method, 42.0)
            pd.testing.assert_series_equal(whole.loc[alone.index], alone, check_freq=False)

    @pytest.mark.parametrize("method", METHODS)
    def test_a_short_gap_is_bridged(self, method: str) -> None:
        daily = _synthetic_daily(n_years=3, start="1990-10-01")
        gapped = daily.drop(daily.index[500:503])
        bf = separate_baseflow(gapped, method, 42.0)
        assert bf.iloc[500:503].isna().all()
        assert bf.iloc[480:500].notna().all() and bf.iloc[503:523].notna().all()

    def test_the_longest_gap_bridged(self) -> None:
        flows = np.full(60, 10.0)
        positions, values = np.array([5, 55]), np.array([10.0, 10.0])
        for missing, bridged in ((15, True), (16, False)):
            gapped = flows.copy()
            gapped[20 : 20 + missing] = np.nan
            bf = _interpolate_baseflow(gapped, positions, values)
            assert np.isnan(bf[20 : 20 + missing]).all()
            assert (not np.isnan(bf[10])) is bridged


class TestLogInterpolation:
    """``ih_smoothed_minima_log``: the same turning points, joined in logarithms."""

    def test_the_line_between_turning_points_is_geometric(self) -> None:
        flows = np.full(5, 100.0)
        bf = _interpolate_baseflow(flows, np.array([0, 4]), np.array([1.0, 16.0]), log=True)
        np.testing.assert_allclose(bf, [1.0, 2.0, 4.0, 8.0, 16.0])

    def test_a_turning_point_of_zero_is_joined_straight(self) -> None:
        flows = np.full(9, 100.0)
        positions, values = np.array([0, 2, 4, 6, 8]), np.array([1.0, 9.0, 0.0, 9.0, 16.0])
        bf = _interpolate_baseflow(flows, positions, values, log=True)
        np.testing.assert_allclose(bf, [1.0, 3.0, 9.0, 4.5, 0.0, 4.5, 9.0, 12.0, 16.0])

    def test_never_above_the_streamflow(self) -> None:
        flows = np.array([1.0, 1.5, 100.0, 100.0, 16.0])
        bf = _interpolate_baseflow(flows, np.array([0, 4]), np.array([1.0, 16.0]), log=True)
        np.testing.assert_allclose(bf, [1.0, 1.5, 4.0, 8.0, 16.0])

    def test_below_the_straight_lines_and_defined_on_the_same_days(self) -> None:
        daily = _synthetic_daily(n_years=3)
        straight = separate_baseflow(daily, "ih_smoothed_minima")
        log = separate_baseflow(daily, "ih_smoothed_minima_log")
        assert (straight.isna() == log.isna()).all()
        assert (log.dropna() <= straight.dropna() + 1e-9).all()
        assert (log.dropna() < straight.dropna() - 1e-6).any()

    def test_block_days_and_factor_are_used(self) -> None:
        daily = _synthetic_daily(n_years=3)
        usual = separate_baseflow(daily, "ih_smoothed_minima_log")
        assert not usual.equals(separate_baseflow(daily, "ih_smoothed_minima_log", block_days=7))
        assert not usual.equals(separate_baseflow(daily, "ih_smoothed_minima_log", factor=0.7))

    def test_against_the_usgs_programs_published_indexes(self) -> None:
        """OFR 03-146's BFI_AVE, and what each method differs from it by
        (``tests/fixtures/bfi/README.md``, where one of the gages is the method's worst of 23)."""
        published = {"04131500": 0.544, "12065500": 0.525, "13077700": 0.660}
        off = {}
        for method in BASEFLOW_METHODS[:3]:
            for site, index in published.items():
                path = Path(__file__).parent / "fixtures" / "bfi" / f"{site}_daily.csv.gz"
                daily = pd.read_csv(path, parse_dates=["date"], index_col="date")
                table = baseflow_index(daily, method=method)
                mean = float(table.loc[table["complete"], "baseflow_index"].mean())
                off[method, site] = round(mean - index, 3)
        assert [off["ih_smoothed_minima", s] for s in published] == [0.015, 0.012, 0.049]
        assert [off["ih_smoothed_minima_log", s] for s in published] == [-0.033, -0.002, 0.015]
        assert [off["ih_smoothed_minima_log_wy", s] for s in published] == [-0.003, -0.003, 0.023]


class TestWaterYearBlocks:
    """``ih_smoothed_minima_log_wy``: the 5-day blocks start again each October 1."""

    METHOD = "ih_smoothed_minima_log_wy"

    @staticmethod
    def _blocks(start: str, days: int, block_days: int = 5) -> list:
        index = pd.date_range(start, periods=days, freq="D")
        edges = _water_year_blocks(index, block_days)
        return [(index[a].strftime("%m-%d"), b - a) for a, b in edges]

    def test_every_year_has_73_blocks_and_the_last_takes_february_29(self) -> None:
        blocks = self._blocks("1990-10-01", 365 + 366 + 365)
        assert len(blocks) == 3 * 73
        assert [blocks[i][0] for i in (0, 73, 146)] == ["10-01"] * 3
        assert [n for _, n in blocks[:73]] == [5] * 73
        assert [n for _, n in blocks[73:146]] == [5] * 72 + [6]

    def test_the_last_block_takes_what_is_left_over(self) -> None:
        blocks = self._blocks("1990-10-01", 365, block_days=7)
        assert len(blocks) == 52 and blocks[-1] == ("09-23", 8)

    def test_partial_blocks_at_the_ends_of_the_record_are_left_out(self) -> None:
        blocks = self._blocks("1990-10-04", 20)
        assert blocks == [("10-06", 5), ("10-11", 5), ("10-16", 5)]
        assert self._blocks("1991-09-27", 3) == []
        # The leap year's last block needs its sixth day.
        assert self._blocks("1992-09-25", 5) == [] and len(self._blocks("1992-09-25", 6)) == 1

    def test_the_first_date_of_the_record_does_not_matter(self) -> None:
        daily = _synthetic_daily(n_years=3, start="1990-10-01")
        whole = separate_baseflow(daily, self.METHOD)
        usual = separate_baseflow(daily, "ih_smoothed_minima_log")
        moved = False
        for shift in (1, 2, 3, 4):
            part = separate_baseflow(daily.iloc[shift:], self.METHOD).iloc[90:]
            pd.testing.assert_series_equal(part, whole.loc[part.index], check_freq=False)
            other = separate_baseflow(daily.iloc[shift:], "ih_smoothed_minima_log").iloc[90:]
            both = other.notna() & usual.loc[other.index].notna()
            moved = moved or not np.allclose(other[both], usual.loc[other.index][both])
        assert moved

    def test_the_first_water_year_is_as_without_the_restart(self) -> None:
        daily = _synthetic_daily(n_years=1, start="1990-10-01").iloc[:365]
        pd.testing.assert_series_equal(
            separate_baseflow(daily, self.METHOD),
            separate_baseflow(daily, "ih_smoothed_minima_log"),
        )

    def test_turning_points_compare_blocks_across_october_1(self) -> None:
        """The last block of a year and the first of the next are neighbors."""
        days = pd.date_range("1990-10-01", "1992-09-30", freq="D")
        flow = pd.Series(50.0, index=days)
        flow["1991-10-01":"1991-10-05"] = 1.0
        bf = separate_baseflow(flow.to_frame("flow_cfs"), self.METHOD)
        # 0.9 x 50 is over the 1 next to it, so the block before is no turning point, and the
        # line runs from September 21, the first day of the block before that, to October 1.
        assert bf["1991-09-21"] == pytest.approx(50.0) and bf["1991-10-01"] == pytest.approx(1.0)
        assert bf["1991-09-26"] == pytest.approx(50.0**0.5)

    def test_block_days_and_factor_are_used(self) -> None:
        daily = _synthetic_daily(n_years=3)
        usual = separate_baseflow(daily, self.METHOD)
        assert not usual.equals(separate_baseflow(daily, self.METHOD, block_days=7))
        assert not usual.equals(separate_baseflow(daily, self.METHOD, factor=0.7))

    def test_never_above_the_streamflow(self) -> None:
        daily = _synthetic_daily(n_years=3)
        bf = separate_baseflow(daily, self.METHOD).dropna()
        assert (bf <= daily["flow_cfs"].loc[bf.index] + 1e-9).all()


class TestHysep:
    """Tests for the three HYSEP variants and their shared N* formula."""

    @pytest.mark.parametrize(
        "area,days",
        [
            (0.5, 3),
            (32.0, 3),
            (59.1, 5),
            (100.0, 5),
            (1000.0, 7),
            (3125.0, 9),
            (4000.0, 11),
            (1e6, 11),
        ],
    )
    def test_interval_days(self, area: float, days: int) -> None:
        """2N* is the odd integer from 3 to 11 nearest 2 * A**0.2 (WRIR 96-4040); the report
        gives 5 for French Creek's 59.1 square miles. A tie (32) goes to the smaller."""
        assert _hysep_interval_days(area) == days

    def test_interval_days_rejects_nonpositive_area(self) -> None:
        with pytest.raises(ValueError, match="must be positive"):
            _hysep_interval_days(0)

    def test_fixed_interval_is_a_step_function(self) -> None:
        """Baseflow is constant within each non-overlapping 2N*-day block."""
        flows = np.array(
            [50, 45, 40, 35, 30, 10, 32, 36, 40, 44, 48, 52, 20, 48, 44, 40, 36, 32, 28, 24],
            dtype=float,
        )
        from flowfreq.regime import _hysep_fixed

        bf = _hysep_fixed(flows, interval_days=3)
        assert bf[0] == bf[1] == bf[2] == flows[0:3].min()
        assert bf[3] == flows[3:6].min()

    def test_sliding_interval_never_exceeds_actual_flow(self) -> None:
        from flowfreq.regime import _hysep_sliding

        flows = np.array(
            [50, 45, 40, 35, 30, 10, 32, 36, 40, 44, 48, 52, 20, 48, 44, 40, 36, 32, 28, 24],
            dtype=float,
        )
        bf = _hysep_sliding(flows, interval_days=3)
        valid = ~np.isnan(bf)
        assert np.all(bf[valid] <= flows[valid] + 1e-9)

    def test_local_minimum_identifies_known_local_minima(self) -> None:
        from flowfreq.regime import _hysep_local_minimum

        flows = np.array(
            [50, 45, 40, 35, 30, 10, 32, 36, 40, 44, 48, 52, 20, 48, 44, 40, 36, 32, 28, 24],
            dtype=float,
        )
        bf = _hysep_local_minimum(flows, interval_days=3)
        assert bf[5] == 10.0
        assert bf[12] == 20.0

    @pytest.mark.parametrize("method", ["hysep_fixed", "hysep_sliding", "hysep_local_minimum"])
    def test_requires_drainage_area(self, method: str) -> None:
        daily = _synthetic_daily(n_years=2)
        with pytest.raises(ValueError, match="requires drainage_area_sqmi"):
            separate_baseflow(daily, method=method)


class TestHysepFrenchCreek:
    """The worked example of USGS WRIR 96-4040: French Creek, 59.1 square miles, 5 days.

    ``tests/fixtures/hysep/README.md`` says where the record and the 11-day offset come from.
    """

    AREA = 59.1
    STREAMFLOW = [42.87, 42.18, 57.45, 46.60, 28.55, 19.30]
    STREAMFLOW += [94.13, 30.39, 21.70, 22.48, 43.10, 50.48]
    BASEFLOW = [31.56, 36.37, 41.56, 36.88, 21.99, 14.53]
    BASEFLOW += [28.06, 24.18, 19.02, 19.96, 19.72, 30.48]

    @pytest.fixture
    def record(self) -> pd.DataFrame:
        path = Path(__file__).parent / "fixtures" / "hysep" / "01472157_daily.csv"
        return pd.read_csv(path, index_col="date", parse_dates=True)

    def test_fixed_interval(self, record: pd.DataFrame) -> None:
        """April 5 to 9, 1992: every day takes the lowest, 49 cfs."""
        bf = separate_baseflow(record.loc["1992-04-05":], "hysep_fixed", self.AREA)
        assert bf.loc["1992-04-05":"1992-04-09"].tolist() == [49.0] * 5

    def test_sliding_interval(self, record: pd.DataFrame) -> None:
        """April 16 to 20, 1992: the lowest, 42 cfs, goes to the median day, April 18."""
        bf = separate_baseflow(record.loc["1992-03-01":], "hysep_sliding", self.AREA)
        assert bf.loc["1992-04-18"] == 42.0
        assert bf.loc["1992-04-19"] > 42.0  # April 17 to 21 does not reach back to it

    def test_local_minimum(self, record: pd.DataFrame) -> None:
        """Figure 4: local minimums on April 9, 15, 21 and 24, 1992."""
        spring = record.loc["1992-03-01":]
        bf = separate_baseflow(spring, "hysep_local_minimum", self.AREA)
        for day in (9, 15, 21, 24):
            date = pd.Timestamp(1992, 4, day)
            assert bf.loc[date] == spring["flow_cfs"].loc[date]

    def test_local_minimum_monthly_summary(self, record: pd.DataFrame) -> None:
        """Appendix 2, 1969, whose dates run 11 days early."""
        listed = record.loc["1968-12-01":"1970-02-28"]
        listed = listed.set_axis(listed.index - pd.Timedelta(days=11))
        flow = listed["flow_cfs"].loc["1969"]
        np.testing.assert_allclose(
            flow.groupby(flow.index.month).mean(), self.STREAMFLOW, atol=0.0051
        )
        bf = separate_baseflow(listed, "hysep_local_minimum", self.AREA).loc["1969"]
        np.testing.assert_allclose(bf.groupby(bf.index.month).mean(), self.BASEFLOW, atol=0.3)
        assert bf.mean() == pytest.approx(27.00, rel=0.005)


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
        from flowfreq.regime import _lyne_hollick

        flows = self._recession_with_storm()
        gapped = flows.copy()
        gapped[30] = np.nan

        bf = _lyne_hollick(gapped)
        assert np.isnan(bf).sum() == 1
        assert np.isnan(bf[30])
        assert not np.isnan(bf[29])
        assert not np.isnan(bf[89])

    def test_recession_is_mostly_baseflow(self) -> None:
        from flowfreq.regime import _lyne_hollick

        flows = self._recession_with_storm()
        bf = _lyne_hollick(flows)
        bfi_recession = bf[:35].sum() / flows[:35].sum()
        assert bfi_recession > 0.85

    def test_storm_peak_is_mostly_quickflow(self) -> None:
        from flowfreq.regime import _lyne_hollick

        flows = self._recession_with_storm()
        bf = _lyne_hollick(flows)
        assert bf[42] < flows[42] * 0.5

    def test_baseflow_bounded_by_zero_and_total_flow(self) -> None:
        from flowfreq.regime import _lyne_hollick

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

    def test_n_days_counts_matched_days_not_raw_flow_days(self) -> None:
        """Regression: n_days must reflect days where BOTH flow_cfs and the
        separated baseflow are valid, not just flow_cfs -- ih_smoothed_minima
        always leaves a multi-day warm-up region NaN before its first
        turning point (and, for an isolated year, another before its last),
        even on perfectly gap-free daily data."""
        rng = np.random.default_rng(3)
        idx = pd.date_range("2010-01-01", periods=365, freq="D")
        flow = np.clip(50 + 5 * np.sin(np.arange(365) / 30) + rng.normal(0, 1, 365), 20, None)
        daily = pd.DataFrame({"flow_cfs": flow}, index=idx)

        bf = separate_baseflow(daily, method="ih_smoothed_minima")
        expected_n_days = int(bf.notna().sum())
        assert expected_n_days < 365, "the reproduction must actually have a warm-up gap"

        bfi = baseflow_index(daily, year_type="calendar", min_days=350)
        row = bfi.iloc[0]
        assert row["n_days"] < 365, "n_days must be less than the full flow record"
        assert row["n_days"] == expected_n_days

    def test_baseflow_index_matches_manually_computed_matched_day_ratio(self) -> None:
        """Regression for the core defect: baseflow_index must equal
        sum(baseflow)/sum(flow) over the SAME (matched) day-set, not
        sum(baseflow over available days)/sum(flow over ALL days). Before
        the fix this reported 0.9604 here; the matched-day-set ratio is
        0.9765 -- a 1.6-point difference from mixing denominators, with no
        warning since 'complete' was (wrongly) True either way."""
        rng = np.random.default_rng(3)
        idx = pd.date_range("2010-01-01", periods=365 * 4, freq="D")
        flow = np.clip(
            50 + 5 * np.sin(np.arange(len(idx)) / 30) + rng.normal(0, 1, len(idx)), 20, None
        )
        daily = pd.DataFrame({"flow_cfs": flow}, index=idx)

        bf = separate_baseflow(daily, method="ih_smoothed_minima")
        aligned_flow = daily["flow_cfs"].reindex(bf.index)
        matched = bf.notna() & aligned_flow.notna()
        year_labels = pd.DatetimeIndex(bf.index).year
        year0 = year_labels[0]
        mask = matched & (year_labels == year0)
        expected = bf[mask].sum() / aligned_flow[mask].sum()

        bfi = baseflow_index(daily, year_type="calendar", min_days=350)
        reported = bfi[bfi["year"] == year0].iloc[0]["baseflow_index"]
        assert reported == pytest.approx(expected, abs=1e-9)
        assert reported == pytest.approx(0.9765, abs=0.001)

    def test_a_storm_inside_the_warmup_window_no_longer_biases_the_ratio(self) -> None:
        """The dramatic case: a large event landing entirely inside the
        method's warm-up window used to inflate the denominator (all of it
        counted in total flow) while contributing nothing to the numerator
        (excluded from the baseflow sum) -- reported BFI came out 0.886
        against a correct 0.951. Must now match the matched-day value."""
        rng = np.random.default_rng(7)
        idx = pd.date_range("2010-01-01", periods=365 * 3, freq="D")
        flow = np.clip(
            50 + 5 * np.sin(np.arange(len(idx)) / 30) + rng.normal(0, 1, len(idx)), 20, None
        )
        flow[0:25] += 300  # spike squarely inside the warm-up window
        daily = pd.DataFrame({"flow_cfs": flow}, index=idx)

        bfi = baseflow_index(daily, year_type="calendar", min_days=350)
        row = bfi[bfi["year"] == 2010].iloc[0]
        assert row["baseflow_index"] == pytest.approx(0.951, abs=0.001)

    def test_complete_true_never_accompanies_a_mismatched_ratio(self) -> None:
        """No dataset should be able to reach complete=True while its
        reported ratio still mixes day-sets -- complete must gate on the
        same matched set the ratio itself uses."""
        for seed in range(5):
            rng = np.random.default_rng(seed)
            idx = pd.date_range("2010-01-01", periods=365 * 2, freq="D")
            flow = np.clip(
                50 + 5 * np.sin(np.arange(len(idx)) / 30) + rng.normal(0, 1, len(idx)), 20, None
            )
            daily = pd.DataFrame({"flow_cfs": flow}, index=idx)
            bfi = baseflow_index(daily, year_type="calendar", min_days=350)

            bf = separate_baseflow(daily, method="ih_smoothed_minima")
            aligned_flow = daily["flow_cfs"].reindex(bf.index)
            matched = bf.notna() & aligned_flow.notna()
            year_labels = pd.DatetimeIndex(bf.index).year

            for _, row in bfi.iterrows():
                if not row["complete"]:
                    continue
                mask = matched & (year_labels == row["year"])
                expected = bf[mask].sum() / aligned_flow[mask].sum()
                assert row["baseflow_index"] == pytest.approx(expected, abs=1e-9)


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

    def test_baseflow_complete_column_present_and_can_diverge_from_complete(self) -> None:
        """Regression: complete (flashiness_index/tqmean's own gate) and
        baseflow_complete (baseflow_index's own, independent gate) are
        different questions and must be reported separately. A short record
        can have a fully-valid daily flow series (complete=True) while
        hysep_local_minimum's long warm-up window still leaves too little
        matched data for baseflow_index (baseflow_complete=False) --
        without the separate column this would show up as an unexplained
        NaN baseflow_index next to complete=True."""
        idx = pd.date_range("2020-01-01", periods=355, freq="D")
        flow = 50 + 5 * np.sin(np.arange(355) / 20)
        daily = pd.DataFrame({"flow_cfs": flow}, index=idx)

        regime = FlowRegime(
            daily,
            year_type="calendar",
            min_days=350,
            baseflow_method="hysep_local_minimum",
            drainage_area_sqmi=500.0,
        )
        assert "baseflow_complete" in regime.annual.columns
        row = regime.annual.iloc[0]
        assert row["complete"]
        assert not row["baseflow_complete"]
        assert np.isnan(row["baseflow_index"])

    def test_pooled_baseflow_index_matches_manual_matched_day_computation(self) -> None:
        """Regression: FlowRegime.summary()'s pooled BFI must match summing
        baseflow and flow over the SAME matched day-set across complete
        years, not baseflow-available-days over all-days -- the same defect
        as baseflow_index() itself, but inside the pooling loop."""
        rng = np.random.default_rng(3)
        idx = pd.date_range("2010-01-01", periods=365 * 4, freq="D")
        flow = np.clip(
            50 + 5 * np.sin(np.arange(len(idx)) / 30) + rng.normal(0, 1, len(idx)), 20, None
        )
        daily = pd.DataFrame({"flow_cfs": flow}, index=idx)

        regime = FlowRegime(daily, year_type="calendar", min_days=350)
        summary = regime.summary()

        bf = separate_baseflow(daily, method="ih_smoothed_minima")
        aligned_flow = daily["flow_cfs"].reindex(bf.index)
        matched = bf.notna() & aligned_flow.notna()
        complete_years = set(regime.annual.loc[regime.annual["baseflow_complete"], "year"])
        year_labels = pd.DatetimeIndex(bf.index).year
        include = matched & pd.Series(year_labels, index=bf.index).isin(complete_years).to_numpy()
        expected = bf[include].sum() / aligned_flow[include].sum()

        assert summary["baseflow_index"] == pytest.approx(expected, abs=1e-9)

    def test_summary_reports_baseflow_n_years_separately(self) -> None:
        daily = _synthetic_daily(n_years=15)
        regime = FlowRegime(daily)
        summary = regime.summary()
        assert summary["baseflow_n_years"] == regime.annual["baseflow_complete"].sum()


def _pacific_diel_series(
    n_days: int = 11, amplitude: float = 30.0, mean: float = 100.0
) -> pd.DataFrame:
    """A UTC-indexed instantaneous series (as download_instantaneous_flow
    returns) with a clean, known diel sinusoid in LOCAL (Pacific) time:
    trough near 4am local, peak near 4pm local, true range = 2*amplitude."""
    local_times = pd.date_range(
        "2020-01-10 00:00", periods=96 * n_days, freq="15min", tz="America/Los_Angeles"
    )
    hour_local = np.array([t.hour + t.minute / 60 for t in local_times])
    flow = mean + amplitude * np.sin(2 * np.pi * (hour_local - 4) / 24 - np.pi / 2)
    return pd.DataFrame({"flow_cfs": flow}, index=local_times.tz_convert("UTC"))


class TestDielVariation:
    """Tests for diel_variation."""

    def test_local_day_grouping_recovers_true_diel_range(self) -> None:
        """The central correctness property: grouping on local calendar days
        (not the UTC index the data is stored in) must recover the true
        diel amplitude on every fully-bracketed day."""
        iv = _pacific_diel_series(n_days=11, amplitude=30.0)
        daily = diel_variation(iv, tz="America/Los_Angeles")
        middle = daily.iloc[2:-2]  # avoid the partial days at the UTC storage edges
        assert (middle["range_cfs"] > 55).all() and (middle["range_cfs"] < 65).all()
        assert middle["complete"].all()

    def test_grouping_by_utc_day_would_give_a_different_wrong_answer(self) -> None:
        """Guard against silently reverting to UTC-day grouping: manually
        regrouping the same data by the raw (UTC) index date must NOT
        recover the true range, confirming the tz conversion is doing
        real work rather than being a no-op for this input."""
        iv = _pacific_diel_series(n_days=5, amplitude=30.0)
        wrong = pd.DataFrame({"date": iv.index.date, "flow_cfs": iv["flow_cfs"].to_numpy()})
        wrong_ranges = wrong.groupby("date")["flow_cfs"].agg(lambda s: s.max() - s.min())
        # At least one UTC-labeled bucket must show a range that is NOT
        # close to the true 60 -- i.e. grouping naively really would corrupt it.
        assert not (wrong_ranges.between(55, 65)).all()

    def test_zero_mean_day_gives_nan_cv_not_error(self) -> None:
        idx = pd.date_range("2020-06-01", periods=4 * 24, freq="15min", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": np.zeros(4 * 24)}, index=idx)
        daily = diel_variation(iv, tz="UTC")
        assert daily.iloc[0]["mean_flow_cfs"] == 0.0
        assert np.isnan(daily.iloc[0]["cv"])
        assert daily.iloc[0]["range_cfs"] == 0.0

    def test_no_runtime_warning_on_zero_mean_division(self) -> None:
        """The cv computation must not emit a divide-by-zero warning even
        though it is vectorized across all days including zero-mean ones."""
        import warnings

        idx = pd.date_range("2020-06-01", periods=4 * 24, freq="15min", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": np.zeros(4 * 24)}, index=idx)
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            diel_variation(iv, tz="UTC")

    def test_negative_values_excluded_from_range(self) -> None:
        """A negative artifact must not become the day's minimum."""
        idx = pd.date_range("2020-06-01", periods=8, freq="3h", tz="UTC")
        vals = np.array([50.0, 55.0, -999.0, 60.0, 65.0, 45.0, 40.0, 52.0])
        iv = pd.DataFrame({"flow_cfs": vals}, index=idx)
        daily = diel_variation(iv, tz="UTC")
        assert daily.iloc[0]["n_obs"] == 7
        assert daily.iloc[0]["min_flow_cfs"] == 40.0
        assert daily.iloc[0]["max_flow_cfs"] == 65.0

    def test_single_observation_day_zero_range_nan_cv(self) -> None:
        """range is defined for n=1 (a point has no spread); cv (needing a
        sample std) is not."""
        idx = pd.date_range("2020-07-01 00:00", periods=8, freq="3h", tz="UTC").append(
            pd.date_range("2020-07-02 12:00", periods=1, freq="3h", tz="UTC")
        )
        vals = np.concatenate([np.linspace(40, 80, 8), [55.0]])
        iv = pd.DataFrame({"flow_cfs": vals}, index=idx)
        daily = diel_variation(iv, tz="UTC")
        single_row = daily[daily["n_obs"] == 1].iloc[0]
        assert single_row["range_cfs"] == 0.0
        assert np.isnan(single_row["cv"])
        assert not single_row["complete"]

    def test_sparse_day_still_reports_value_but_marked_incomplete(self) -> None:
        """Departure from the annual/monthly/seasonal metrics: an incomplete
        day's range/cv are reported as computed (still informative at
        partial coverage), not suppressed to NaN -- only the completeness
        flag signals reduced precision."""
        idx3 = pd.date_range("2020-07-01 00:00", periods=8, freq="3h", tz="UTC")
        sparse_day = pd.date_range("2020-07-02 06:00", periods=2, freq="6h", tz="UTC")
        full_idx = idx3.append(sparse_day)
        vals3 = np.concatenate([np.linspace(40, 80, 8), [55.0, 58.0]])
        iv = pd.DataFrame({"flow_cfs": vals3}, index=full_idx)
        daily = diel_variation(iv, tz="UTC")

        sparse_row = daily[daily["n_obs"] == 2].iloc[0]
        assert sparse_row["range_cfs"] == 3.0
        assert not np.isnan(sparse_row["cv"])
        assert not sparse_row["complete"]

    def test_expected_obs_inferred_from_median_step_not_hardcoded(self) -> None:
        """A record sampled hourly (not NWIS's usual 15-minute) must infer
        ~24 expected observations per day, not a hardcoded 96."""
        idx = pd.date_range("2020-01-01", periods=24 * 5, freq="1h", tz="UTC")
        iv = pd.DataFrame({"flow_cfs": 50.0 + 10 * np.sin(np.arange(len(idx)) / 4)}, index=idx)
        daily = diel_variation(iv, tz="UTC")
        assert daily["expected_obs"].iloc[0] == pytest.approx(24.0)

    def test_too_few_timestamps_raises(self) -> None:
        idx = pd.date_range("2020-01-01", periods=1, tz="UTC")
        iv = pd.DataFrame({"flow_cfs": [50.0]}, index=idx)
        with pytest.raises(ValueError, match="at least 2 timestamps"):
            diel_variation(iv, tz="UTC")

    def test_requires_tz_aware_index(self) -> None:
        """A naive (non-tz-aware) index cannot be tz_convert'd; this should
        surface as an error rather than silently misbehaving."""
        idx = pd.date_range("2020-01-01", periods=10, freq="1h")  # no tz
        iv = pd.DataFrame({"flow_cfs": 50.0}, index=idx)
        with pytest.raises(TypeError):
            diel_variation(iv, tz="UTC")


class TestDielVariationSummary:
    """Tests for diel_variation_summary."""

    def test_summary_matches_known_amplitude(self) -> None:
        iv = _pacific_diel_series(n_days=11, amplitude=30.0)
        daily = diel_variation(iv, tz="America/Los_Angeles")
        summary = diel_variation_summary(daily)
        assert summary["mean_diel_amplitude_cfs"] == pytest.approx(60.0, abs=0.5)
        assert summary["n_days"] == daily["complete"].sum()

    def test_only_complete_days_included(self) -> None:
        iv = _pacific_diel_series(n_days=5, amplitude=30.0)
        daily = diel_variation(iv, tz="America/Los_Angeles")
        daily_with_incomplete = daily.copy()
        daily_with_incomplete.loc[0, "complete"] = False
        daily_with_incomplete.loc[0, "range_cfs"] = 99999.0  # sentinel that must be excluded

        summary = diel_variation_summary(daily_with_incomplete)
        assert summary["n_days"] == len(daily_with_incomplete) - 1
        assert summary["mean_diel_amplitude_cfs"] < 1000  # sentinel correctly excluded

    def test_composes_with_an_arbitrary_caller_defined_slice(self) -> None:
        """No grouping logic of its own: filtering to a sub-period before
        calling is how 'a period' is chosen."""
        iv = _pacific_diel_series(n_days=11, amplitude=30.0)
        daily = diel_variation(iv, tz="America/Los_Angeles")
        first_half = daily.iloc[: len(daily) // 2]
        summary = diel_variation_summary(first_half)
        assert summary["n_days"] <= len(first_half)

    def test_empty_input_gives_nan_not_error(self) -> None:
        empty = pd.DataFrame(columns=["date", "range_cfs", "cv", "complete"])
        summary = diel_variation_summary(empty)
        assert summary["n_days"] == 0
        assert np.isnan(summary["mean_diel_amplitude_cfs"])
        assert np.isnan(summary["mean_diel_cv"])


class TestLyneHollickAgainstLadson:
    """Ladson and others (2013), the program BFI.R of 1 October 2013 (TonyLadson on GitHub).

    It reflects 30 days at each end and starts the filter at the first flow; ours does neither,
    so the two differ near the ends of a record and nowhere else.
    """

    #: Bass River at Loch, 67 daily flows (the paper's own sample).
    BASS = [
        5,
        7,
        108,
        117,
        57,
        36,
        26,
        95,
        1169,
        308,
        144,
        89,
        62,
        48,
        40,
        35,
        73,
        82,
        342,
        393,
        310,
        275,
        260,
        245,
        256,
        141,
        119,
        934,
        382,
        158,
        96,
        122,
        103,
        83,
        67,
        148,
        366,
        161,
        119,
        82,
        330,
        294,
        261,
        266,
        153,
        247,
        703,
        498,
        286,
        163,
        124,
        85,
        94,
        81,
        62,
        47,
        37,
        30,
        26,
        24,
        24,
        22,
        21,
        20,
        19,
        18,
        18,
    ]

    @staticmethod
    def _ladson(flows: np.ndarray, alpha: float, reflect: int = 30) -> np.ndarray:
        def one_pass(y: np.ndarray) -> np.ndarray:
            quick = np.empty(len(y))
            quick[0] = y[0]
            for i in range(1, len(y)):
                quick[i] = alpha * quick[i - 1] + 0.5 * (1 + alpha) * (y[i] - y[i - 1])
            return np.where(quick > 0, y - quick, y)

        padded = np.concatenate([flows[reflect:0:-1], flows, flows[-2 : -reflect - 2 : -1]])
        base = one_pass(one_pass(one_pass(padded)[::-1])[::-1])
        return np.maximum(base[reflect:-reflect], 0.0)

    def test_the_transcription_gives_the_programs_printed_index(self):
        flows = np.array(self.BASS, dtype=float)
        index = self._ladson(flows, 0.98).sum() / flows.sum()
        assert index == pytest.approx(0.1965759, abs=5e-8)  # printed in BFI_example.R

    def test_away_from_the_ends_of_a_record_ours_is_the_same(self):
        frame = pd.read_csv(Path(__file__).parent / "fixtures" / "bfi" / "13077700_daily.csv.gz")
        flows = frame["flow_cfs"].dropna().to_numpy(dtype=float)
        ours, theirs = _lyne_hollick(flows, 0.925), self._ladson(flows, 0.925)
        assert len(flows) == 11049
        apart = np.abs(ours - theirs)
        assert apart[200:-200].max() < 1e-6 * flows.mean()
        assert ours.sum() / theirs.sum() == pytest.approx(1.0, abs=3e-4)

    def test_a_record_of_two_months_is_not(self):
        flows = np.array(self.BASS, dtype=float)
        ours = _lyne_hollick(flows, 0.98).sum() / flows.sum()
        assert ours == pytest.approx(0.1022, abs=5e-5)  # against 0.1966
