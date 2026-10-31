"""Regression tests for the second independent review of 2026-10-31.

One test (or more) per confirmed defect in monthly/duration transfer, naturalize,
gls and stationarity.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from flowfreq.duration_transfer import duration_statistics
from flowfreq.monthly_transfer import build_monthly_pairs, monthly_leave_one_out, monthly_statistics
from flowfreq.naturalize import Diversion, naturalize
from flowfreq.regression.gls import _nearest_psd
from flowfreq.stationarity import pettitt

YEARS = list(range(2000, 2012))
F = np.array([1.0, 1.3, 0.8, 1.1, 0.6, 1.5, 0.9, 1.2, 0.7, 1.4, 1.05, 0.95])


def _daily(scale: float, f=F) -> pd.Series:
    days = pd.date_range("1999-10-01", "2011-09-30", freq="D")
    wy = days.year + (days.month >= 10)
    fy = dict(zip(YEARS, f))
    return pd.Series([scale * (50.0 + 10 * m) * fy[y] for m, y in zip(days.month, wy)], index=days)


def _network(n: int = 8):
    rng = np.random.default_rng(0)
    sites = [f"g{k}" for k in range(n)]
    daily = {
        s: _daily(float(rng.uniform(0.5, 3.0)), f=F * rng.uniform(0.8, 1.2, 12)) for s in sites
    }
    gages = pd.DataFrame(
        {"site_no": sites, "drainage_area_sqmi": rng.uniform(20, 2000, n),
         "runoff_depth": rng.uniform(10, 30, n)}
    )  # fmt: skip
    stats = monthly_statistics(daily, start_water_year=2000, end_water_year=2011)
    return stats, gages, [(sites[k], sites[k + 1]) for k in range(n - 1)]


# 1. duplicate or self-loop connections ------------------------------------------


@pytest.mark.parametrize("extra, match", [(("g3", "g2"), "listed twice"), (("g3", "g3"), "itself")])
def test_duplicate_and_self_connections_are_refused(extra, match) -> None:
    stats, gages, conns = _network()
    with pytest.raises(ValueError, match=match):
        build_monthly_pairs(stats, 6, gages, conns + [extra])
    with pytest.raises(ValueError, match=match):
        monthly_leave_one_out(stats, 6, gages, conns + [extra])


# 2. naturalize index alignment --------------------------------------------------


def test_rate_aligns_by_date_whatever_the_time_of_day_or_timezone() -> None:
    days = pd.date_range("2020-01-01", "2020-03-31", freq="D", tz="America/Denver")
    gaged = pd.Series(100.0, index=days)
    rate = pd.Series(10.0, index=pd.date_range("2020-01-01 12:00", periods=91, freq="D"))
    out = naturalize(gaged, [Diversion("c", rate_cfs=rate, start="2020-02-01")])
    feb_on = np.asarray(days >= pd.Timestamp("2020-02-01", tz="America/Denver"))
    np.testing.assert_allclose(out.flow[feb_on], 110.0)
    np.testing.assert_allclose(out.flow[~feb_on], 100.0)


def test_a_rate_covering_few_days_in_effect_warns(caplog) -> None:
    days = pd.date_range("2020-01-01", "2020-03-31", freq="D")
    monthly = pd.Series(10.0, index=pd.date_range("2020-01-01", periods=3, freq="MS"))
    with caplog.at_level("WARNING", logger="flowfreq.naturalize"):
        out = naturalize(pd.Series(100.0, index=days), [Diversion("m", rate_cfs=monthly)])
    assert out.summary.loc["m", "days_without_rate"] == 88
    assert any("rate" in r.message and "m" in r.message for r in caplog.records)


def test_a_rate_with_duplicate_dates_is_refused() -> None:
    days = pd.date_range("2020-01-01", "2020-01-10", freq="D")
    rate = pd.Series(1.0, index=pd.DatetimeIndex(["2020-01-02", "2020-01-02", "2020-01-03"]))
    with pytest.raises(ValueError, match="duplicate"):
        naturalize(pd.Series(100.0, index=days), [Diversion("d", rate_cfs=rate)])


# 3. completeness counts distinct days; non-daily input refused --------------------


def test_a_duplicated_day_does_not_hide_a_missing_one() -> None:
    d = _daily(1.0)
    d = d.drop(pd.Timestamp("2005-03-01"))
    d = pd.concat([d, d.loc[["2005-03-02"]]]).sort_index()
    with pytest.raises(ValueError, match="one value per day"):
        duration_statistics({"a": d}, start_water_year=2000, end_water_year=2011)
    with pytest.raises(ValueError, match="one value per day"):
        monthly_statistics({"a": d}, start_water_year=2000, end_water_year=2011)


def test_sub_daily_input_is_refused() -> None:
    hourly = pd.Series(1.0, index=pd.date_range("1999-10-01", "2011-09-30 23:00", freq="h"))
    with pytest.raises(ValueError, match="one value per day"):
        duration_statistics({"a": hourly}, start_water_year=2000, end_water_year=2011)


# 4. monthly messages unchanged by the refactor -----------------------------------


def test_monthly_messages_keep_their_wording() -> None:
    def dry(s: pd.Series) -> pd.Series:
        s = s.copy()
        s.loc[(s.index.month == 9) & (s.index.year == 2003)] = 0.0
        return s

    st = monthly_statistics({"a": dry(_daily(1.0))}, start_water_year=2000,
                            end_water_year=2011, statistic="median")  # fmt: skip
    reasons = set(st.excluded["reason"])
    assert "zero flow in a yearly mean (median works in log space)" in reasons
    stats, gages, conns = _network()
    with pytest.raises(KeyError, match="monthly statistics"):
        build_monthly_pairs(stats, 6, gages.assign(site_no=gages["site_no"].replace("g0", "zz")),
                            [("zz", "g1")])  # fmt: skip


# 5. projection keeps the known diagonal ------------------------------------------


def test_nearest_psd_keeps_the_diagonal() -> None:
    S = np.array([[1.0, 0.9, -0.9], [0.9, 1.0, 0.9], [-0.9, 0.9, 1.0]]) * 0.04
    assert np.linalg.eigvalsh(S)[0] < 0
    P = _nearest_psd(S)
    np.testing.assert_allclose(np.diag(P), np.diag(S))
    assert np.linalg.eigvalsh(P)[0] >= -1e-12
    np.testing.assert_allclose(P, P.T)


# 6. Pettitt shift direction --------------------------------------------------------


def test_pettitt_shift_from_the_statistic_not_medians() -> None:
    assert pettitt([1, 1, 1, 1, 1, 1, 9, 9, 9, 9, 9]).shift == "up"
    assert pettitt([9, 9, 9, 9, 9, 1, 1, 1, 1, 1, 1]).shift == "down"
    assert pettitt([1, 1, 1, 1, 1, 1, 9, 9, 1, 1, 1]).shift == "up"
