"""Naturalization add-back (methow item 11).

Natural flow = gaged flow + sum over upstream diversions of rate x (1 - return fraction)
on the days each is in effect. Inputs are the caller's: flowfreq does not estimate
diversion rates or return flows.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from flowfreq.naturalize import Diversion, naturalize

DAYS = pd.date_range("2020-10-01", "2021-09-30", freq="D")
GAGED = pd.Series(100.0, index=DAYS)


def test_constant_diversion_with_return_flow() -> None:
    out = naturalize(GAGED, [Diversion("canal", rate_cfs=20.0, return_fraction=0.25)])
    np.testing.assert_allclose(out.flow, 115.0)
    assert list(out.added.columns) == ["canal"]
    assert out.summary.loc["canal", "mean_added_cfs"] == pytest.approx(15.0)
    assert out.summary.loc["canal", "days_in_effect"] == len(DAYS)


def test_season_and_period_limit_the_add_back() -> None:
    d = Diversion("ditch", rate_cfs=10.0, months=(5, 6, 7, 8, 9), start="2021-06-01")
    out = naturalize(GAGED, [d])
    june_on = (DAYS >= "2021-06-01") & DAYS.month.isin([5, 6, 7, 8, 9])
    np.testing.assert_allclose(out.flow[june_on], 110.0)
    np.testing.assert_allclose(out.flow[~june_on], 100.0)
    assert out.summary.loc["ditch", "days_in_effect"] == int(june_on.sum())


def test_daily_rate_series_and_several_diversions() -> None:
    rate = pd.Series(np.linspace(0.0, 30.0, len(DAYS)), index=DAYS)
    out = naturalize(
        GAGED,
        [Diversion("a", rate_cfs=rate, return_fraction=0.5), Diversion("b", rate_cfs=4.0)],
    )
    np.testing.assert_allclose(out.flow, 100.0 + 0.5 * rate + 4.0)
    np.testing.assert_allclose(out.added["a"] + out.added["b"], out.flow - GAGED)


def test_a_rate_series_missing_days_adds_nothing_there_and_says_so() -> None:
    rate = pd.Series(8.0, index=DAYS[:100])
    out = naturalize(GAGED, [Diversion("short", rate_cfs=rate)])
    np.testing.assert_allclose(out.flow.iloc[:100], 108.0)
    np.testing.assert_allclose(out.flow.iloc[100:], 100.0)
    assert out.summary.loc["short", "days_without_rate"] == len(DAYS) - 100


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"rate_cfs": -1.0}, "rate"),
        ({"rate_cfs": 5.0, "return_fraction": 1.5}, "return_fraction"),
        ({"rate_cfs": 5.0, "months": (0, 6)}, "months"),
        ({"rate_cfs": 5.0, "start": "2021-06-01", "end": "2021-01-01"}, "start"),
    ],
)
def test_diversion_checks(kwargs: dict, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        Diversion("bad", **kwargs)


def test_gaged_series_and_names_are_checked() -> None:
    with pytest.raises(TypeError):
        naturalize(pd.Series([1.0, 2.0]), [Diversion("x", rate_cfs=1.0)])
    with pytest.raises(ValueError, match="duplicate"):
        naturalize(GAGED, [Diversion("x", rate_cfs=1.0), Diversion("x", rate_cfs=2.0)])
