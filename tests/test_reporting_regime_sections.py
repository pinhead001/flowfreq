"""The `monthly_statistics`, `summary_hydrograph` and `baseflow` report sections.

Nothing here makes a request: the daily record is synthetic and carried inline.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np
import pandas as pd
import pytest

from flowfreq.regime import monthly_flow_summary
from flowfreq.reporting import ReportError, build_section
from flowfreq.reporting.blocks import Figure, Table
from flowfreq.reporting.regime_sections import MONTHLY_CAPTION
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
