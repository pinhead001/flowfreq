"""Duration-quantile transfer variance (methow item 7, part 3).

For each complete water year, a gage's daily flows give the flow equaled or exceeded on p
percent of days (q10, q50, q90; Weibull plotting position). Over a window of years the
statistic is the median of those yearly values (the median annual flow-duration curve) or
their mean, with log10 sampling variance and concurrent-year covariance as in
``flowfreq.monthly_transfer``. Pairs, the GLS blend and leave-one-out follow the monthly
machinery, keyed by exceedance percentage instead of month.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from flowfreq.duration_transfer import (
    DurationStatistics,
    blend_duration,
    build_duration_pairs,
    duration_leave_one_out,
    duration_statistics,
    duration_transfers,
)
from flowfreq.transfer_variance import combine, fit_dissimilarity

YEARS = list(range(2000, 2012))
F = np.array([1.0, 1.3, 0.8, 1.1, 0.6, 1.5, 0.9, 1.2, 0.7, 1.4, 1.05, 0.95])
BASE = {m: 20.0 + 15.0 * m for m in range(1, 13)}  # within-year variation by month


def _daily(scale: float, f=F, years=YEARS, edit=None) -> pd.Series:
    days = pd.date_range(f"{years[0] - 1}-10-01", f"{years[-1]}-09-30", freq="D")
    wy = days.year + (days.month >= 10)
    fy = dict(zip(years, f))
    vals = np.array([scale * BASE[m] * fy[y] for m, y in zip(days.month, wy)], dtype=float)
    s = pd.Series(vals, index=days)
    return edit(s) if edit is not None else s


def _yearly_q(series: pd.Series, wy: int, p: int) -> float:
    s = series[(series.index >= f"{wy - 1}-10-01") & (series.index <= f"{wy}-09-30")]
    return float(np.quantile(s.to_numpy(), 1.0 - p / 100.0, method="weibull"))


def _stats(daily, statistic="median", **kw) -> DurationStatistics:
    return duration_statistics(
        daily, start_water_year=2000, end_water_year=2011, statistic=statistic, **kw
    )


def _row(st: DurationStatistics, site: str, p: int) -> pd.Series:
    t = st.table
    return t[(t["site_no"] == site) & (t["exceedance"] == p)].iloc[0]


def test_yearly_duration_quantiles_and_median_statistic() -> None:
    d = _daily(2.0)
    st = _stats({"a": d})
    assert st.exceedances == (10, 50, 90)
    for p in (10, 50, 90):
        x = np.array([_yearly_q(d, y, p) for y in YEARS])
        r = _row(st, "a", p)
        lx = np.log10(x)
        assert r["value_cfs"] == pytest.approx(np.median(x))
        assert r["log_q"] == pytest.approx(np.median(lx))
        assert r["n_years"] == 12
        assert r["variance"] == pytest.approx(math.pi / 2 * lx.var(ddof=1) / 12)
    q10, q90 = _row(st, "a", 10)["value_cfs"], _row(st, "a", 90)["value_cfs"]
    assert q10 > q90  # q10 is the high flow, exceeded on only 10 percent of days


def test_mean_statistic_and_covariance() -> None:
    fb = np.array([1.1, 1.2, 0.9, 1.0, 0.7, 1.3, 1.0, 1.1, 0.8, 1.2, 1.0, 0.9])
    da, db = _daily(1.0), _daily(3.0, f=fb)
    st = _stats({"a": da, "b": db}, statistic="mean")
    xa = np.array([_yearly_q(da, y, 50) for y in YEARS])
    xb = np.array([_yearly_q(db, y, 50) for y in YEARS])
    ln10_2 = math.log(10.0) ** 2
    expected = 12 * np.cov(xa, xb, ddof=1)[0, 1] / (12 * 12 * xa.mean() * xb.mean() * ln10_2)
    assert st.covariance[50].loc["a", "b"] == pytest.approx(expected)
    assert st.covariance[50].loc["a", "a"] == pytest.approx(_row(st, "a", 50)["variance"])


def test_incomplete_water_year_is_dropped() -> None:
    def gap(s: pd.Series) -> pd.Series:
        return s.drop(s.index[(s.index >= "2004-03-01") & (s.index <= "2004-03-05")])

    st = _stats({"a": _daily(1.0, edit=gap)})
    assert _row(st, "a", 50)["n_years"] == 11
    assert 2004 not in set(st.yearly.loc[st.yearly["site_no"] == "a", "water_year"])


def test_zero_low_flows_are_excluded_and_reported() -> None:
    def dry(s: pd.Series) -> pd.Series:
        s = s.copy()
        s.loc[s.index.month.isin([8, 9, 10])] = 0.0  # a quarter of each year dry
        return s

    st = _stats({"a": _daily(1.0, edit=dry)})
    assert math.isnan(_row(st, "a", 90)["log_q"])
    assert not math.isnan(_row(st, "a", 10)["log_q"])
    ex = st.excluded
    assert list(ex.loc[ex["site_no"] == "a", "exceedance"]) == [90]


@pytest.mark.parametrize(
    "kw, match",
    [
        ({"exceedances": (0, 50)}, "exceedance"),
        ({"exceedances": (50, 100)}, "exceedance"),
        ({"statistic": "mode"}, "statistic"),
        ({"min_years": 1}, "min_years"),
    ],
)
def test_argument_checks(kw: dict, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        _stats({"a": _daily(1.0)}, **kw)


GAGES = pd.DataFrame(
    {
        "site_no": ["a", "b", "c"],
        "drainage_area_sqmi": [100.0, 300.0, 50.0],
        "runoff_depth": [20.0, 15.0, 25.0],
    }
)


def test_pairs_use_the_runoff_volume_ratio_and_exceedance_key() -> None:
    st = _stats({"a": _daily(1.0), "b": _daily(2.0), "c": _daily(0.6)})
    pairs = build_duration_pairs(st, 90, GAGES, [("a", "b")])
    ab = pairs.iloc[0]
    la, lb = _row(st, "a", 90)["log_q"], _row(st, "b", 90)["log_q"]
    log_ar, log_rr = math.log10(3.0), math.log10(15.0 / 20.0)
    assert ab["exceedance"] == 90 and "month" not in pairs.columns
    assert ab["transfer"] == pytest.approx(la + log_ar + log_rr)
    assert ab["error"] == pytest.approx(la + log_ar + log_rr - lb)


def test_unknown_exceedance_is_refused() -> None:
    st = _stats({"a": _daily(1.0), "b": _daily(2.0)})
    with pytest.raises(KeyError, match="25"):
        build_duration_pairs(st, 25, GAGES, [("a", "b")])


def _network(n: int = 14, seed: int = 1):
    rng = np.random.default_rng(seed)
    sites = [f"g{k}" for k in range(n)]
    daily = {
        s: _daily(float(rng.uniform(0.5, 3.0)), f=F * rng.uniform(0.8, 1.2, 12)) for s in sites
    }
    gages = pd.DataFrame(
        {"site_no": sites, "drainage_area_sqmi": rng.uniform(20, 2000, n),
         "runoff_depth": rng.uniform(10, 30, n)}
    )  # fmt: skip
    conns = [(sites[k], sites[k + 1]) for k in range(n - 1)]
    return _stats(daily), gages, conns


def test_blend_matches_combine_and_leave_one_out_runs() -> None:
    st, gages, conns = _network()
    model = fit_dissimilarity(build_duration_pairs(st, 50, gages, conns))
    res = blend_duration(model, st, 50, gages, "g4", ["g3", "g5"], donor_correlation=0.1)
    transfers, dcov = duration_transfers(st, 50, gages, "g4", ["g3", "g5"])
    ref = combine(model, transfers, None, donor_cov=dcov, donor_correlation=0.1)
    assert res.log_q == pytest.approx(ref.log_q) and res.variance == pytest.approx(ref.variance)
    out = duration_leave_one_out(st, 50, gages, conns)
    assert list(out["site_no"]) == list(gages["site_no"])
    assert out["exceedance"].eq(50).all() and (out["note"] == "").all()
    assert np.allclose(out["blend_error"], out["blend"] - out["observed"])
