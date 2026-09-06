"""Tests for flowfreq.hydrograph.Hydrograph.

Daily-series plotting and summary statistics. No network: daily data is built
by hand as a DatetimeIndex/`flow_cfs` DataFrame, the same shape
`tests/test_report.py` uses for its stub gage.
"""

from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from flowfreq.hydrograph import Hydrograph

matplotlib.use("Agg")


def _daily_frame(days: int = 400, start: str = "2000-10-01") -> pd.DataFrame:
    """A DatetimeIndex/flow_cfs DataFrame spanning more than a water year."""
    idx = pd.date_range(start, periods=days, freq="D")
    flow = 500 + 400 * np.sin(np.arange(days) / 58.0) ** 2
    return pd.DataFrame({"flow_cfs": flow}, index=idx)


@pytest.fixture
def daily():
    return _daily_frame()


class TestDayOfWaterYear:
    """Water year runs Oct 1 - Sep 30; day 1 is Oct 1."""

    def test_oct_1_is_day_1(self):
        assert Hydrograph.day_of_water_year(10, 1) == 1

    def test_sep_30_is_the_last_day(self):
        assert Hydrograph.day_of_water_year(9, 30) == 365

    def test_monotonic_across_the_water_year(self):
        """Every day of a non-leap water year, in Oct-Sep order, strictly increases."""
        months_days = [(10, d) for d in range(1, 32)]
        months_days += [(11, d) for d in range(1, 31)]
        months_days += [(12, d) for d in range(1, 32)]
        for m in (1, 2, 3, 4, 5, 6, 7, 8, 9):
            days_in_m = {1: 31, 2: 28, 3: 31, 4: 30, 5: 31, 6: 30, 7: 31, 8: 31, 9: 30}[m]
            months_days += [(m, d) for d in range(1, days_in_m + 1)]
        values = [Hydrograph.day_of_water_year(m, d) for m, d in months_days]
        assert values == sorted(values)
        assert values == list(range(1, 366))

    def test_non_numeric_month_raises(self):
        """No bounds checking at all -- month 13 or -1 silently wraps via `%`
        rather than raising, so the only reachable error is a bad type."""
        with pytest.raises(TypeError):
            Hydrograph.day_of_water_year("October", 1)


class TestPlotDailyTimeseries:
    def test_returns_figure(self, daily):
        fig = Hydrograph.plot_daily_timeseries(daily)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_site_name_and_no_appear_in_title(self, daily):
        fig = Hydrograph.plot_daily_timeseries(
            daily, site_name="Big Sandy River at Bruceton, TN", site_no="03606500"
        )
        title = fig.axes[0].get_title()
        assert "03606500" in title
        assert "Big Sandy" in title
        plt.close(fig)

    def test_log_yscale_sets_power_of_ten_ticks(self, daily):
        fig = Hydrograph.plot_daily_timeseries(daily, yscale="log")
        ax = fig.axes[0]
        assert ax.get_yscale() == "log"
        ticks = ax.get_yticks()
        assert len(ticks) > 0
        assert all(t > 0 for t in ticks)
        plt.close(fig)

    def test_save_path_writes_a_file(self, daily, tmp_path):
        out = tmp_path / "daily.png"
        fig = Hydrograph.plot_daily_timeseries(daily, save_path=str(out))
        assert out.is_file()
        plt.close(fig)

    def test_missing_flow_column_raises(self, daily):
        with pytest.raises(KeyError):
            Hydrograph.plot_daily_timeseries(daily.rename(columns={"flow_cfs": "q"}))


class TestPlotSummaryHydrograph:
    def test_returns_figure(self, daily):
        fig = Hydrograph.plot_summary_hydrograph(daily)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_x_axis_spans_the_water_year(self, daily):
        fig = Hydrograph.plot_summary_hydrograph(daily)
        ax = fig.axes[0]
        assert ax.get_xlim() == (1.0, 366.0)
        assert list(ax.get_xticklabels()[0].get_text() for _ in [0])[0] in Hydrograph.MONTH_LABELS
        plt.close(fig)

    def test_fewer_than_four_percentiles_omits_the_bands(self, daily):
        """The two fill_between bands are gated on len(percentiles) >= 4."""
        fig = Hydrograph.plot_summary_hydrograph(daily, percentiles=[10, 50, 90])
        assert len(fig.axes[0].collections) == 0
        plt.close(fig)

    def test_default_percentiles_produce_two_bands(self, daily):
        fig = Hydrograph.plot_summary_hydrograph(daily)
        assert len(fig.axes[0].collections) == 2
        plt.close(fig)

    def test_missing_flow_column_raises(self, daily):
        with pytest.raises(KeyError):
            Hydrograph.plot_summary_hydrograph(daily.rename(columns={"flow_cfs": "q"}))


class TestGetSummaryStats:
    def test_columns_and_index_range(self, daily):
        stats = Hydrograph.get_summary_stats(daily)
        assert list(stats.columns) == [
            "day_of_water_year",
            "mean",
            "min",
            "max",
            "p10",
            "p25",
            "p50",
            "p75",
            "p90",
        ]
        assert stats["day_of_water_year"].min() >= 1
        assert stats["day_of_water_year"].max() <= 366

    def test_min_le_median_le_max_every_day(self, daily):
        stats = Hydrograph.get_summary_stats(daily)
        assert (stats["min"] <= stats["p50"]).all()
        assert (stats["p50"] <= stats["max"]).all()

    def test_custom_percentiles_are_honoured(self, daily):
        stats = Hydrograph.get_summary_stats(daily, percentiles=[5, 95])
        assert list(stats.columns) == ["day_of_water_year", "mean", "min", "max", "p5", "p95"]

    def test_missing_flow_column_raises(self, daily):
        with pytest.raises(KeyError):
            Hydrograph.get_summary_stats(daily.rename(columns={"flow_cfs": "q"}))


class TestPlotFlowDurationCurve:
    def test_returns_figure_and_table(self, daily):
        fig, table = Hydrograph.plot_flow_duration_curve(daily)
        assert isinstance(fig, plt.Figure)
        assert isinstance(table, pd.DataFrame)
        plt.close(fig)

    def test_table_columns(self, daily):
        _, table = Hydrograph.plot_flow_duration_curve(daily)
        assert list(table.columns) == ["Percent Exceeded", "Flow (cfs)"]
        assert len(table) == 9
        plt.close(plt.gcf())

    def test_flow_decreases_as_percent_exceeded_increases(self, daily):
        """A flow exceeded 99% of the time must be smaller than one exceeded 1%."""
        _, table = Hydrograph.plot_flow_duration_curve(daily)
        flows = table["Flow (cfs)"].to_numpy()
        assert np.all(np.diff(flows) <= 0)
        plt.close(plt.gcf())

    def test_table_path_writes_csv(self, daily, tmp_path):
        out = tmp_path / "fdc.csv"
        fig, table = Hydrograph.plot_flow_duration_curve(daily, table_path=str(out))
        assert out.is_file()
        written = pd.read_csv(out)
        assert list(written["Flow (cfs)"]) == pytest.approx(list(table["Flow (cfs)"]))
        plt.close(fig)

    def test_empty_dataframe_raises(self):
        empty = pd.DataFrame({"flow_cfs": []}, index=pd.DatetimeIndex([]))
        with pytest.raises(ValueError):
            Hydrograph.plot_flow_duration_curve(empty)
