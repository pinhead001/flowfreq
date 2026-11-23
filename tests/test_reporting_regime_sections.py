"""The `monthly_statistics`, `summary_hydrograph` and `baseflow` report sections.

Nothing here makes a request: the daily record is synthetic and carried inline.
"""

from __future__ import annotations

from typing import Any, Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from flowfreq.hydrograph import Hydrograph
from flowfreq.regime import (
    BASEFLOW_METHODS,
    FlowRegime,
    baseflow_index,
    monthly_flow_summary,
    separate_baseflow,
)
from flowfreq.reporting import ReportError, build_section, regime_sections
from flowfreq.reporting.blocks import Figure, Table
from flowfreq.reporting.regime_sections import (
    _METHOD_WORDS,
    BASEFLOW_CAPTION,
    MONTHLY_CAPTION,
    SUMMARY_HYDROGRAPH_CAPTION,
    _typical_year,
)
from flowfreq.reporting.sections import cfs
from flowfreq.reporting.session import MAIN_TABLES
from tests.test_reporting_daily import FRAME, START, STATION, VALUES, _table, _text

MONTHS = [
    "October",
    "November",
    "December",
    "January",
    "February",
    "March",
    "April",
    "May",
    "June",
    "July",
    "August",
    "September",
]


def _station(values: Any, start: str = START) -> Dict[str, Any]:
    """The test station with another daily record."""
    flows = [None if v is None or np.isnan(v) else float(v) for v in values]
    return {**STATION, "daily": {"start": start, "flow_cfs": flows}}


def _figures(blocks: List[Any]) -> List[Figure]:
    return [b for b in blocks if isinstance(b, Figure)]


def _tables(blocks: List[Any]) -> List[Table]:
    return [b for b in blocks if isinstance(b, Table)]


@pytest.fixture
def drawn(monkeypatch: pytest.MonkeyPatch) -> List[Any]:
    """The axes of every figure a section draws, kept open to look at."""
    axes: List[Any] = []

    def keep(fig: Any) -> bytes:
        fig.canvas.draw()
        axes.append(fig.axes[0])
        plt.close(fig)
        return b""

    monkeypatch.setattr(regime_sections, "_figure_png", keep)
    return axes


class TestWhatIsDrawn:
    """Each figure shows its data: the limits hold it and the axis runs upward."""

    @staticmethod
    def _holds(ax: Any, low: float, high: float) -> bool:
        bottom, top = ax.get_ylim()
        return bottom <= low and high <= top

    @staticmethod
    def _tick_labels(ax: Any) -> List[str]:
        return [t.get_text() for t in ax.get_yticklabels() if t.get_text()]

    @pytest.mark.parametrize("floor", [0.0, 0.02, 25.0])
    def test_baseflow(self, drawn: List[Any], floor: float) -> None:
        values = np.maximum(VALUES - VALUES.min(), 0.0) + floor
        build_section("baseflow", _station(values))
        ax, bars = drawn
        line = ax.get_lines()[0].get_ydata()
        assert self._holds(ax, line.min(), line.max())
        assert ax.get_yscale() == ("log" if line.min() > 0 else "linear")
        low = ax.collections[0].get_paths()[0].vertices[:, 1].min()
        assert ax.get_ylim()[0] <= low * (1 + 1e-9)
        assert not [t for t in self._tick_labels(ax) if "10^" in t or t == "0" and floor]
        assert len(bars.patches) == 29

    def test_baseflow_of_an_intermittent_stream_is_linear(self, drawn: List[Any]) -> None:
        values = VALUES.copy()
        values[::50] = 0.0
        build_section("baseflow", _station(values))
        ax = drawn[0]
        assert ax.get_yscale() == "linear" and ax.get_ylim()[0] == 0.0
        assert ax.get_ylim()[1] >= ax.get_lines()[0].get_ydata().max()

    def test_monthly(self, drawn: List[Any]) -> None:
        build_section("monthly_statistics", dict(STATION))
        (ax,) = drawn
        means = monthly_flow_summary(FRAME).query("complete")["mean_flow_cfs"]
        assert self._holds(ax, means.min(), means.max())
        labels = self._tick_labels(ax)
        assert labels and not [t for t in labels if "10^" in t]

    def test_monthly_with_a_dry_month_is_linear(self, drawn: List[Any]) -> None:
        values = VALUES.copy()
        values[:30] = 0.0  # April 1990
        build_section("monthly_statistics", _station(values))
        assert drawn[0].get_yscale() == "linear"

    def test_summary(self, drawn: List[Any]) -> None:
        build_section("summary_hydrograph", dict(STATION))
        (ax,) = drawn
        assert self._holds(ax, VALUES.min(), VALUES.max())
        assert ax.get_yscale() == "log"
        stats = Hydrograph.get_summary_stats(FRAME)
        for line, name in zip(ax.get_lines(), ("p50", "min", "max")):
            np.testing.assert_array_equal(line.get_xdata(), stats["day_of_water_year"])
            np.testing.assert_allclose(line.get_ydata(), stats[name])

    def test_summary_of_an_intermittent_stream_is_linear(self, drawn: List[Any]) -> None:
        values = VALUES.copy()
        values[::50] = 0.0
        build_section("summary_hydrograph", _station(values))
        assert drawn[0].get_yscale() == "linear"

    def test_the_separated_year_is_the_one_named(self, drawn: List[Any]) -> None:
        built = build_section("baseflow", dict(STATION))
        ax = drawn[0]
        year = int(ax.get_xlabel().removeprefix("Water year "))
        flows = FRAME["flow_cfs"].loc[f"{year - 1}-10-01":f"{year}-09-30"]
        np.testing.assert_array_equal(ax.get_lines()[0].get_ydata(), flows)
        assert f"water year {year}" in _figures(built.blocks)[0].caption

    def test_missing_days_are_a_break_in_the_separated_hydrograph(self, drawn: List[Any]) -> None:
        values = VALUES.copy()
        missing = (FRAME.index.month == 5) & (FRAME.index.day <= 10)
        values[missing] = np.nan
        build_section("baseflow", {**_station(values), "min_days": 300})
        line = drawn[0].get_lines()[0]
        year = int(drawn[0].get_xlabel().removeprefix("Water year "))
        assert len(line.get_ydata()) == (366 if pd.Timestamp(year, 1, 1).is_leap_year else 365)
        assert np.isnan(line.get_ydata()).sum() == 10


class TestMonthlyStatistics:
    @pytest.fixture
    def built(self) -> Any:
        return build_section("monthly_statistics", dict(STATION))

    def test_statistics_are_those_of_the_complete_monthly_means(self, built: Any) -> None:
        monthly = monthly_flow_summary(FRAME)
        monthly = monthly[monthly["complete"]]
        table = _table(built.blocks, MONTHLY_CAPTION)
        assert MAIN_TABLES["monthly_statistics"] == MONTHLY_CAPTION
        assert [row[0] for row in table.rows] == MONTHS
        for row, number in zip(table.rows, [10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]):
            rows = monthly[monthly["month"] == number]
            means = rows["mean_flow_cfs"]
            most, least = rows.loc[means.idxmax()], rows.loc[means.idxmin()]
            assert row[1:7] == (
                f"{len(rows)}",
                cfs(means.mean()),
                cfs(means.median()),
                cfs(means.std(ddof=1)),
                f"{cfs(most['mean_flow_cfs'])} ({int(most['year'])})",
                f"{cfs(least['mean_flow_cfs'])} ({int(least['year'])})",
            )

    def test_percents_of_annual_volume_sum_to_100(self, built: Any) -> None:
        table = _table(built.blocks, MONTHLY_CAPTION)
        assert sum(float(row[7]) for row in table.rows) == pytest.approx(100.0, abs=0.6)

    def test_percents_of_annual_volume_weigh_each_month_by_its_length(self) -> None:
        # A constant flow: the share of a month is its share of the days.
        built = build_section("monthly_statistics", _station(np.full(1461, 50.0), "1991-10-01"))
        shares = {row[0]: row[7] for row in _table(built.blocks, MONTHLY_CAPTION).rows}
        assert shares["October"] == f"{100 * 31 / 365.25:.1f}"
        assert shares["February"] == f"{100 * 28.25 / 365.25:.1f}"

    def test_title_method_and_figure(self, built: Any) -> None:
        assert built.title == "Monthly flow statistics at USGS 01234567 Test Creek near Nowhere"
        assert "at least 90 percent of its days" in _text(built.blocks)
        assert len(_figures(built.blocks)) == 1
        assert len(_tables(built.blocks)) == 2

    def test_an_incomplete_month_is_left_out_and_counted(self) -> None:
        values = VALUES.copy()
        values[:20] = np.nan  # April 1990
        table = _table(
            build_section("monthly_statistics", _station(values)).blocks, MONTHLY_CAPTION
        )
        counts = {row[0]: row[1] for row in table.rows}
        assert counts["April"] == "29" and counts["May"] == "30"

    def test_one_year_has_no_standard_deviation(self) -> None:
        built = build_section("monthly_statistics", _station(VALUES[:365]))
        rows = _table(built.blocks, MONTHLY_CAPTION).rows
        assert {row[1] for row in rows} == {"1"} and {row[4] for row in rows} == {"-"}

    def test_a_month_with_no_complete_year_is_a_row_of_dashes(self) -> None:
        built = build_section("monthly_statistics", _station(VALUES[:60]))
        rows = {row[0]: row for row in _table(built.blocks, MONTHLY_CAPTION).rows}
        assert rows["October"][1:] == ("0", "-", "-", "-", "-", "-", "-")
        assert rows["April"][1] == "1"
        assert rows["April"][7] == "-"  # a share of the year needs all twelve months

    def test_annual_table_on_request(self) -> None:
        built = build_section("monthly_statistics", {**STATION, "annual_table": True})
        table = _tables(built.blocks)[-1]
        assert table.header == ("Water year", *[m[:3] for m in MONTHS])
        assert table.rows[0][:7] == ("1990", "-", "-", "-", "-", "-", "-")  # starts in April
        monthly = monthly_flow_summary(FRAME).set_index(["year", "month"])["mean_flow_cfs"]
        assert table.rows[1][1] == cfs(monthly[(1991, 10)])
        assert len(table.rows) == 31

    def test_annual_table_leaves_out_a_year_with_no_complete_month(self) -> None:
        values = VALUES.copy()
        values[(FRAME.index >= "1995-09-25") & (FRAME.index <= "1998-10-10")] = np.nan
        built = build_section("monthly_statistics", {**_station(values), "annual_table": True})
        years = [row[0] for row in _tables(built.blocks)[-1].rows]
        assert years == [str(y) for y in range(1990, 2021) if y not in (1996, 1997, 1998)]

    def test_negative_values_are_missing(self) -> None:
        values = VALUES.copy()
        values[:30] = -5.0
        built = build_section("monthly_statistics", _station(values))
        rows = {row[0]: row for row in _table(built.blocks, MONTHLY_CAPTION).rows}
        assert rows["April"][1] == "29"
        assert "30 negative value(s)" in _text(built.blocks)

    def test_zero_flows_draw(self) -> None:
        built = build_section("monthly_statistics", _station(np.zeros(800)))
        rows = _table(built.blocks, MONTHLY_CAPTION).rows
        assert {row[2] for row in rows} == {"0"} and {row[7] for row in rows} == {"-"}

    def test_rebuilds_from_its_own_options(self, built: Any) -> None:
        again = build_section("monthly_statistics", built.options)
        assert _tables(again.blocks) == _tables(built.blocks)

    @pytest.mark.parametrize(
        "change,said",
        [
            ({"annual_table": "yes"}, "annual_table must be true or false"),
            ({"n_day": 7}, r"unknown option\(s\) n_day"),
            ({"daily": {"start": START, "flow_cfs": [5.0] * 10}}, "no complete month"),
            ({"daily": {"start": START, "flow_cfs": [-5.0] * 40}}, "no value of zero or more"),
        ],
    )
    def test_refusals(self, change: Dict[str, Any], said: str) -> None:
        with pytest.raises(ReportError, match=said):
            build_section("monthly_statistics", {**STATION, **change})


class TestSummaryHydrograph:
    @pytest.fixture
    def built(self) -> Any:
        return build_section("summary_hydrograph", dict(STATION))

    def test_rows_are_the_hydrograph_own_statistics(self, built: Any) -> None:
        stats = Hydrograph.get_summary_stats(FRAME).set_index("day_of_water_year")
        table = _table(built.blocks, SUMMARY_HYDROGRAPH_CAPTION)
        assert MAIN_TABLES["summary_hydrograph"] == SUMMARY_HYDROGRAPH_CAPTION
        assert len(table.rows) == 24
        assert [row[0] for row in table.rows[:3]] == ["Oct 1", "Oct 15", "Nov 1"]
        for row, (month, day) in zip(table.rows, [(10, 1), (10, 15), (11, 1), (11, 15)]):
            here = stats.loc[Hydrograph.day_of_water_year(month, day)]
            names = ("min", "p10", "p25", "p50", "p75", "p90", "max")
            assert row[1:] == ("30", *(cfs(here[name]) for name in names))

    def test_percentiles_are_named_as_non_exceedance(self, built: Any) -> None:
        text = _text(built.blocks)
        assert "percentiles of non-exceedance" in text
        assert "q90" in text  # the 10th percentile by the name the duration section uses
        assert built.title == "Summary hydrograph at USGS 01234567 Test Creek near Nowhere"
        assert len(_figures(built.blocks)) == 1

    def test_a_date_without_a_value_is_a_row_of_dashes(self) -> None:
        built = build_section("summary_hydrograph", _station(VALUES[:30]))  # April 1990
        rows = {row[0]: row for row in _table(built.blocks, SUMMARY_HYDROGRAPH_CAPTION).rows}
        assert rows["Oct 1"][1:] == ("0", *["-"] * 7)
        assert rows["Apr 15"][1] == "1"

    def test_negative_values_are_missing(self) -> None:
        values = VALUES.copy()
        values[0] = -9.0  # April 1, 1990
        built = build_section("summary_hydrograph", _station(values))
        rows = {row[0]: row for row in _table(built.blocks, SUMMARY_HYDROGRAPH_CAPTION).rows}
        assert rows["Apr 1"][1] == "29" and "-9" not in rows["Apr 1"]
        assert "1 negative value(s)" in _text(built.blocks)

    def test_zero_flows_draw(self) -> None:
        built = build_section("summary_hydrograph", _station(np.zeros(800)))
        assert len(_figures(built.blocks)) == 1

    def test_rebuilds_from_its_own_options(self, built: Any) -> None:
        again = build_section("summary_hydrograph", built.options)
        assert _tables(again.blocks) == _tables(built.blocks)

    @pytest.mark.parametrize(
        "change,said",
        [
            ({"percentiles": [5, 95]}, r"unknown option\(s\) percentiles"),
            (
                {"daily": {"start": "1990-04-02", "flow_cfs": [5.0] * 10}},
                "no value on the 1st or 15th of any month",
            ),
            ({"daily": {"start": START, "flow_cfs": [-5.0] * 40}}, "no value of zero or more"),
        ],
    )
    def test_refusals(self, change: Dict[str, Any], said: str) -> None:
        with pytest.raises(ReportError, match=said):
            build_section("summary_hydrograph", {**STATION, **change})


class TestBaseflow:
    @pytest.fixture
    def built(self) -> Any:
        return build_section("baseflow", dict(STATION))

    def test_the_index_is_the_pooled_one(self, built: Any) -> None:
        regime = FlowRegime(FRAME)
        rows = dict(_table(built.blocks, BASEFLOW_CAPTION).rows)
        assert MAIN_TABLES["baseflow"] == BASEFLOW_CAPTION
        assert rows["Baseflow index"] == f"{regime.summary()['baseflow_index']:.3f}"
        assert rows["Complete water years"] == "29"
        assert rows["Separation method"] == "ih_smoothed_minima"

    def test_means_are_over_the_days_with_both_flows(self, built: Any) -> None:
        series = separate_baseflow(FRAME)
        both = pd.DataFrame({"flow": FRAME["flow_cfs"], "base": series}).dropna()
        both = both.loc["1990-10-01":"2019-09-30"]
        rows = dict(_table(built.blocks, BASEFLOW_CAPTION).rows)
        assert rows["Mean streamflow (cfs)"] == cfs(both["flow"].mean())
        assert rows["Mean baseflow (cfs)"] == cfs(both["base"].mean())
        inches = 365.25 * 86400 * 12 / (42.0 * 5280**2)
        assert rows["Mean streamflow (inches per year)"] == f"{both['flow'].mean() * inches:.2f}"
        assert rows["Mean baseflow (inches per year)"] == f"{both['base'].mean() * inches:.2f}"

    def test_no_inches_without_a_drainage_area(self) -> None:
        options = {k: v for k, v in STATION.items() if k != "drainage_area_sqmi"}
        table = _table(build_section("baseflow", options).blocks, BASEFLOW_CAPTION)
        assert not [row for row in table.rows if "inches" in row[0]]

    def test_annual_table(self, built: Any) -> None:
        annual = baseflow_index(FRAME)
        annual = annual[annual["complete"]]
        table = _tables(built.blocks)[-1]
        assert table.header == ("Water year", "Days used", "Baseflow index")
        assert table.rows == tuple(
            (f"{year}", f"{days}", f"{index:.3f}")
            for year, days, index in zip(annual["year"], annual["n_days"], annual["baseflow_index"])
        )
        assert "2 water year(s) with too few days" in _text(built.blocks)
        brief = build_section("baseflow", {**STATION, "annual_table": False})
        assert len(_tables(brief.blocks)) == len(_tables(built.blocks)) - 1

    def test_title_cautions_and_figures(self, built: Any) -> None:
        assert built.title == "Baseflow at USGS 01234567 Test Creek near Nowhere"
        text = _text(built.blocks)
        assert "not comparable" in text and "regulation" in text and "snowmelt" in text
        assert len(_figures(built.blocks)) == 2

    @pytest.mark.parametrize("method", BASEFLOW_METHODS)
    def test_every_method_is_named_and_cited(self, method: str) -> None:
        built = build_section("baseflow", {**STATION, "method": method})
        words, citation = _METHOD_WORDS[method]
        assert words in _text(built.blocks)
        assert built.references == [citation]
        regime = FlowRegime(FRAME, baseflow_method=method, drainage_area_sqmi=42.0)
        rows = dict(_table(built.blocks, BASEFLOW_CAPTION).rows)
        assert rows["Baseflow index"] == f"{regime.summary()['baseflow_index']:.3f}"
        assert ("interval of 5 days" in _text(built.blocks)) == method.startswith("hysep")

    def test_methods_have_words_for_each(self) -> None:
        assert tuple(_METHOD_WORDS) == BASEFLOW_METHODS

    def test_typical_year_is_nearest_the_median_and_earliest_on_a_tie(self) -> None:
        assert _typical_year(pd.Series({2001: 5.0, 2002: 9.0, 2003: 7.5})) == 2003
        assert _typical_year(pd.Series({2001: 5.0, 2002: 9.0})) == 2001

    def test_a_year_of_zero_flow_has_no_index(self) -> None:
        values = VALUES.copy()
        dry = (FRAME.index >= "1999-10-01") & (FRAME.index <= "2000-09-30")
        values[dry] = 0.0
        built = build_section("baseflow", _station(values))
        rows = {row[0]: row for row in _tables(built.blocks)[-1].rows}
        assert rows["2000"][2] == "-"
        assert len(_figures(built.blocks)) == 2

    def test_negative_values_are_missing(self) -> None:
        values = VALUES.copy()
        values[400:430] = -5.0
        built = build_section("baseflow", _station(values))
        assert "30 negative value(s)" in _text(built.blocks)
        years = [row[0] for row in _tables(built.blocks)[-1].rows]
        assert "1991" not in years and "1992" in years

    def test_rebuilds_from_its_own_options(self, built: Any) -> None:
        again = build_section("baseflow", built.options)
        assert _tables(again.blocks) == _tables(built.blocks)

    @pytest.mark.parametrize(
        "change,said",
        [
            ({"method": "eckhardt"}, "method must be one of ih_smoothed_minima"),
            ({"method": "hysep_fixed", "drainage_area_sqmi": None}, "give drainage_area_sqmi"),
            ({"drainage_area_sqmi": 0.0}, "drainage_area_sqmi must be a positive number"),
            ({"drainage_area_sqmi": "big"}, "drainage_area_sqmi must be a positive number"),
            ({"min_days": 0}, "min_days must be a whole number of days from 1 to 366"),
            ({"annual_table": 1}, "annual_table must be true or false"),
            ({"daily": {"start": START, "flow_cfs": [5.0] * 200}}, "no complete water year"),
            ({"daily": {"start": START, "flow_cfs": [0.0] * 900}}, "every flow.*is zero"),
            ({"n_day": 7}, r"unknown option\(s\) n_day"),
        ],
    )
    def test_refusals(self, change: Dict[str, Any], said: str) -> None:
        with pytest.raises(ReportError, match=said):
            build_section("baseflow", {**STATION, **change})
