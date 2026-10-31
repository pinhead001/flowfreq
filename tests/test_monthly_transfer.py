"""Monthly transfer-variance inputs: per-month statistics and gage pairs (methow item 7).

Extends ``flowfreq.transfer_variance`` from flood quantiles to monthly flows. For each
calendar month, every gage's yearly monthly means over a water-year window give two
statistics: their **mean** (the period monthly mean) and their **median** (the "typical"
month, e.g. the September median a hydraulic model uses for least width). Each statistic
comes with its log10 sampling variance and, from concurrent water years, the sampling
covariance between gages. Pairs are transferred by the runoff-volume ratio
``(At Rt) / (Ag Rg)`` with exponent 1, so ``fit_dissimilarity`` fits ``s2_d`` per month
unchanged. Decisions recorded in TODO.md (2026-10-26): exponent 1 on runoff volume;
empirical concurrent-year covariance; zero flows excluded and reported.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from flowfreq.monthly_transfer import MonthlyStatistics, build_monthly_pairs, monthly_statistics
from flowfreq.transfer_variance import fit_dissimilarity

LN10_2 = math.log(10.0) ** 2
YEARS = list(range(2000, 2012))  # 12 water years
F = np.array([1.0, 1.3, 0.8, 1.1, 0.6, 1.5, 0.9, 1.2, 0.7, 1.4, 1.05, 0.95])
BASE = {m: 100.0 + 10.0 * m for m in range(1, 13)}


def _daily(scale: float, years=YEARS, f=F, edit=None) -> pd.Series:
    """Daily flow = scale * BASE[month] * f[water year], constant within each month."""
    days = pd.date_range(f"{years[0] - 1}-10-01", f"{years[-1]}-09-30", freq="D")
    wy = days.year + (days.month >= 10)
    fy = dict(zip(years, f))
    vals = np.array([scale * BASE[m] * fy[y] for m, y in zip(days.month, wy)], dtype=float)
    s = pd.Series(vals, index=days)
    if edit is not None:
        s = edit(s)
    return s


def _stats(daily, statistic="mean", **kw) -> MonthlyStatistics:
    return monthly_statistics(daily, start_water_year=2000, end_water_year=2011,
                              statistic=statistic, **kw)  # fmt: skip


def _row(st: MonthlyStatistics, site: str, month: int) -> pd.Series:
    t = st.table
    return t[(t["site_no"] == site) & (t["month"] == month)].iloc[0]


# --------------------------------------------------------------------------- statistics


def test_mean_statistic_value_and_variance() -> None:
    st = _stats({"a": _daily(1.0)})
    r = _row(st, "a", 9)
    x = BASE[9] * F
    assert r["value_cfs"] == pytest.approx(x.mean())
    assert r["log_q"] == pytest.approx(math.log10(x.mean()))
    assert r["n_years"] == 12
    expected_var = x.var(ddof=1) / (12 * x.mean() ** 2 * LN10_2)
    assert r["variance"] == pytest.approx(expected_var)
    assert st.covariance[9].loc["a", "a"] == pytest.approx(expected_var)
    assert set(st.table["month"]) == set(range(1, 13))


def test_median_statistic_value_and_variance() -> None:
    st = _stats({"a": _daily(1.0)}, statistic="median")
    r = _row(st, "a", 9)
    lx = np.log10(BASE[9] * F)
    assert r["value_cfs"] == pytest.approx(np.median(BASE[9] * F))
    assert r["log_q"] == pytest.approx(np.median(lx))
    assert r["variance"] == pytest.approx(math.pi / 2 * lx.var(ddof=1) / 12)


def test_yearly_series_is_kept() -> None:
    st = _stats({"a": _daily(2.0)})
    y = st.yearly
    assert list(y.columns) == ["site_no", "water_year", "month", "flow_cfs"]
    sep = y[(y["month"] == 9)].sort_values("water_year")
    np.testing.assert_allclose(sep["flow_cfs"], 2.0 * BASE[9] * F)


def test_perfectly_correlated_gages_have_zero_known_variance() -> None:
    for statistic in ("mean", "median"):
        st = _stats({"a": _daily(1.0), "b": _daily(3.0)}, statistic=statistic)
        cov = st.covariance[5]
        known = cov.loc["a", "a"] + cov.loc["b", "b"] - 2 * cov.loc["a", "b"]
        assert known == pytest.approx(0.0, abs=1e-15), statistic


def test_mean_covariance_over_concurrent_years() -> None:
    # Gage d starts 4 years later and has its own year factors.
    fd = np.array([0.9, 1.6, 0.7, 1.2, 1.0, 1.1, 0.8, 1.3])
    st = _stats({"a": _daily(1.0), "d": _daily(1.0, years=YEARS[4:], f=fd)}, min_years=5)
    xa, xd = BASE[7] * F, BASE[7] * fd
    xa_c = xa[4:]
    s_xy = np.cov(xa_c, xd, ddof=1)[0, 1]
    expected = 8 * s_xy / (12 * 8 * xa.mean() * xd.mean() * LN10_2)
    assert st.covariance[7].loc["a", "d"] == pytest.approx(expected)
    assert _row(st, "d", 7)["n_years"] == 8


def test_median_covariance_uses_arcsin_of_the_correlation() -> None:
    fb = np.array([1.1, 1.2, 0.9, 1.0, 0.7, 1.3, 1.0, 1.1, 0.8, 1.2, 1.0, 0.9])
    st = _stats({"a": _daily(1.0), "b": _daily(1.0, f=fb)}, statistic="median")
    la, lb = np.log10(BASE[3] * F), np.log10(BASE[3] * fb)
    rho = np.corrcoef(la, lb)[0, 1]
    expected = la.std(ddof=1) * lb.std(ddof=1) * math.asin(rho) / 12
    assert st.covariance[3].loc["a", "b"] == pytest.approx(expected)


def test_incomplete_month_drops_that_year() -> None:
    def gap(s: pd.Series) -> pd.Series:
        s = s.copy()
        s.loc["2005-09-10"] = np.nan
        return s

    st = _stats({"a": _daily(1.0, edit=gap)})
    assert _row(st, "a", 9)["n_years"] == 11
    assert _row(st, "a", 8)["n_years"] == 12


def test_too_few_years_is_excluded_and_reported() -> None:
    st = _stats({"a": _daily(1.0), "d": _daily(1.0, years=YEARS[4:], f=F[4:])}, min_years=10)
    ex = st.excluded
    assert list(ex.columns) == ["site_no", "month", "reason"]
    assert set(ex[ex["site_no"] == "d"]["month"]) == set(range(1, 13))
    assert ex["reason"].str.contains("years").all()
    assert np.isnan(_row(st, "d", 1)["log_q"])


def test_zero_flows_are_excluded_and_reported() -> None:
    def dry_one_september(s: pd.Series) -> pd.Series:
        s = s.copy()
        s.loc["2003-09-01":"2003-09-30"] = 0.0
        return s

    def dry_septembers(s: pd.Series) -> pd.Series:
        s = s.copy()
        s.loc[s.index.month == 9] = 0.0
        return s

    one = {"c": _daily(1.0, edit=dry_one_september)}
    # The median works in log space, so one zero year excludes the gage-month ...
    med = _stats(one, statistic="median")
    assert ((med.excluded["site_no"] == "c") & (med.excluded["month"] == 9)).any()
    assert med.excluded["reason"].str.contains("zero").any()
    # ... the mean only needs a positive period mean.
    mean = _stats(one)
    assert mean.excluded.empty
    assert _row(mean, "c", 9)["n_years"] == 12
    # Every September dry: the period mean is zero, excluded for the mean too.
    allzero = _stats({"c": _daily(1.0, edit=dry_septembers)})
    assert (
        (allzero.excluded["month"] == 9) & allzero.excluded["reason"].str.contains("zero")
    ).any()
    assert np.isnan(_row(allzero, "c", 9)["log_q"])


@pytest.mark.parametrize(
    "kw,match",
    [
        (dict(statistic="mode"), "statistic"),
        (dict(start_water_year=2011, end_water_year=2000), "water year"),
        (dict(min_years=1), "min_years"),
    ],
)
def test_statistics_argument_checks(kw: dict, match: str) -> None:
    args = dict(start_water_year=2000, end_water_year=2011, statistic="mean")
    args.update(kw)
    with pytest.raises(ValueError, match=match):
        monthly_statistics({"a": _daily(1.0)}, **args)


def test_non_daily_index_is_refused() -> None:
    s = _daily(1.0)
    s.index = range(len(s))
    with pytest.raises(TypeError, match="DatetimeIndex"):
        _stats({"a": s})


# --------------------------------------------------------------------------- pairs

GAGES = pd.DataFrame(
    {
        "site_no": ["a", "b", "c"],
        "drainage_area_sqmi": [100.0, 300.0, 50.0],
        "runoff_depth": [20.0, 15.0, 25.0],
        "precip_in": [40.0, 32.0, 46.0],
    }
)


def test_pairs_use_the_runoff_volume_ratio() -> None:
    fb = np.array([1.1, 1.2, 0.9, 1.0, 0.7, 1.3, 1.0, 1.1, 0.8, 1.2, 1.0, 0.9])
    st = _stats({"a": _daily(1.0), "b": _daily(2.0, f=fb), "c": _daily(0.6)})
    pairs = build_monthly_pairs(st, 9, GAGES, [("a", "b"), ("c", "a")], attributes=["precip_in"])
    ab = pairs.iloc[0]
    la, lb = _row(st, "a", 9)["log_q"], _row(st, "b", 9)["log_q"]
    log_ar, log_rr = math.log10(3.0), math.log10(15.0 / 20.0)
    assert (ab["donor"], ab["receiver"], ab["month"], ab["statistic"]) == ("a", "b", 9, "mean")
    assert ab["log_ar"] == pytest.approx(log_ar)
    assert ab["log_rr"] == pytest.approx(log_rr)
    assert ab["exponent"] == 1.0
    assert ab["transfer"] == pytest.approx(la + log_ar + log_rr)
    assert ab["error"] == pytest.approx(la + log_ar + log_rr - lb)
    cov = st.covariance[9]
    assert ab["known_variance"] == pytest.approx(
        cov.loc["a", "a"] + cov.loc["b", "b"] - 2 * cov.loc["a", "b"]
    )
    assert ab["d_precip_in"] == pytest.approx(-8.0)
    assert len(pairs) == 2


def test_runoff_depth_of_one_is_the_area_ratio() -> None:
    st = _stats({"a": _daily(1.0), "b": _daily(3.0)})
    flat = GAGES.assign(runoff_depth=1.0)
    p = build_monthly_pairs(st, 4, flat, [("a", "b")]).iloc[0]
    assert p["log_rr"] == 0.0
    assert p["error"] == pytest.approx(0.0, abs=1e-12)  # b = 3a, areas 100 -> 300


def test_pairs_with_an_excluded_gage_are_dropped() -> None:
    def dry(s: pd.Series) -> pd.Series:
        s = s.copy()
        s.loc[s.index.month == 9] = 0.0
        return s

    st = _stats({"a": _daily(1.0), "b": _daily(3.0), "c": _daily(1.0, edit=dry)})
    pairs = build_monthly_pairs(st, 9, GAGES, [("a", "b"), ("c", "a")])
    assert list(zip(pairs["donor"], pairs["receiver"])) == [("a", "b")]
    assert len(build_monthly_pairs(st, 8, GAGES, [("a", "b"), ("c", "a")])) == 2


def test_pairs_feed_fit_dissimilarity() -> None:
    rng = np.random.default_rng(0)
    sites = [f"g{k}" for k in range(16)]
    daily = {
        s: _daily(float(rng.uniform(0.5, 3.0)), f=F * rng.uniform(0.8, 1.2, 12)) for s in sites
    }
    gages = pd.DataFrame(
        {"site_no": sites, "drainage_area_sqmi": rng.uniform(20, 2000, 16),
         "runoff_depth": rng.uniform(10, 30, 16)}
    )  # fmt: skip
    st = _stats(daily)
    conns = [(sites[k], sites[k + 1]) for k in range(15)]
    model = fit_dissimilarity(build_monthly_pairs(st, 6, gages, conns))
    assert model.n_pairs == 15 and model.coefficients["const"] >= 0.0


@pytest.mark.parametrize("month", [0, 13])
def test_month_is_checked(month: int) -> None:
    st = _stats({"a": _daily(1.0), "b": _daily(3.0)})
    with pytest.raises(ValueError, match="month"):
        build_monthly_pairs(st, month, GAGES, [("a", "b")])


def test_runoff_depth_column_is_required() -> None:
    st = _stats({"a": _daily(1.0), "b": _daily(3.0)})
    with pytest.raises(ValueError, match="runoff_depth"):
        build_monthly_pairs(st, 4, GAGES.drop(columns="runoff_depth"), [("a", "b")])


def test_unknown_site_in_connections_is_refused() -> None:
    st = _stats({"a": _daily(1.0), "b": _daily(3.0)})
    with pytest.raises(KeyError, match="zz"):
        build_monthly_pairs(st, 4, GAGES, [("a", "zz")])


# ------------------------------------------------------------ part 2: transfers, blend, LOO
#
# Methow item 7, part 2: monthly leave-one-out and the GLS blend, mirroring
# flowfreq.transfer_variance.leave_one_out / combine for flood quantiles.


def _network(n: int = 16, seed: int = 0):
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


def test_monthly_transfers_frame_and_donor_covariance() -> None:
    from flowfreq.monthly_transfer import monthly_transfers

    st = _stats({"a": _daily(1.0), "b": _daily(2.0), "c": _daily(0.6)})
    transfers, dcov = monthly_transfers(st, 9, GAGES, "b", ["a", "c"], attributes=["precip_in"])
    la = _row(st, "a", 9)["log_q"]
    a = transfers.set_index("donor").loc["a"]
    assert a["log_ar"] == pytest.approx(math.log10(3.0))
    assert a["log_rr"] == pytest.approx(math.log10(15.0 / 20.0))
    assert a["transfer"] == pytest.approx(la + math.log10(3.0) + math.log10(15.0 / 20.0))
    assert a["donor_variance"] == pytest.approx(st.covariance[9].loc["a", "a"])
    assert a["d_precip_in"] == pytest.approx(-8.0)
    assert list(dcov.index) == ["a", "c"] and dcov.loc["a", "c"] == pytest.approx(
        st.covariance[9].loc["a", "c"]
    )


def test_monthly_transfers_refuse_an_excluded_donor() -> None:
    from flowfreq.monthly_transfer import monthly_transfers

    def dry(s: pd.Series) -> pd.Series:
        s = s.copy()
        s.loc[s.index.month == 9] = 0.0
        return s

    st = _stats({"a": _daily(1.0), "b": _daily(2.0), "c": _daily(1.0, edit=dry)})
    with pytest.raises(ValueError, match="c.*excluded"):
        monthly_transfers(st, 9, GAGES, "b", ["a", "c"])


def test_blend_month_matches_combine() -> None:
    from flowfreq.monthly_transfer import blend_month, monthly_transfers
    from flowfreq.transfer_variance import combine

    st, gages, conns = _network()
    model = fit_dissimilarity(build_monthly_pairs(st, 6, gages, conns))
    res = blend_month(model, st, 6, gages, "g5", ["g4", "g6"], donor_correlation=0.2)
    transfers, dcov = monthly_transfers(st, 6, gages, "g5", ["g4", "g6"])
    ref = combine(model, transfers, None, donor_cov=dcov, donor_correlation=0.2)
    assert res.log_q == pytest.approx(ref.log_q) and res.variance == pytest.approx(ref.variance)
    with_reg = blend_month(model, st, 6, gages, "g5", ["g4", "g6"], regression=(2.5, 0.02))
    assert with_reg.variance < res.variance or with_reg.variance == pytest.approx(res.variance)


def test_monthly_leave_one_out_scores_every_connected_gage() -> None:
    from flowfreq.monthly_transfer import monthly_leave_one_out

    st, gages, conns = _network()
    out = monthly_leave_one_out(st, 6, gages, conns)
    assert list(out["site_no"]) == list(gages["site_no"])
    assert set(out["n_donors"]) == {1, 2}
    assert out["month"].eq(6).all() and out["statistic"].eq("mean").all()
    ok = out["note"] == ""
    assert ok.all()
    assert np.allclose(out["blend_error"], out["blend"] - out["observed"])
    # Holding out g5 refits s2_d without its pairs, then blends g4's and g6's transfers.
    from flowfreq.monthly_transfer import blend_month

    pairs = build_monthly_pairs(st, 6, gages, conns)
    keep = (pairs["donor"] != "g5") & (pairs["receiver"] != "g5")
    model = fit_dissimilarity(pairs[keep])
    ref = blend_month(model, st, 6, gages, "g5", ["g4", "g6"])
    g5 = out.set_index("site_no").loc["g5"]
    assert g5["blend"] == pytest.approx(ref.log_q)
    assert g5["blend_variance"] == pytest.approx(ref.variance)
    assert g5["blend_covered"] in (True, False)


def test_monthly_leave_one_out_skips_excluded_gages_and_reports_refit_failures() -> None:
    from flowfreq.monthly_transfer import monthly_leave_one_out

    def dry(s: pd.Series) -> pd.Series:
        s = s.copy()
        s.loc[s.index.month == 9] = 0.0
        return s

    st = _stats({"a": _daily(1.0), "b": _daily(3.0), "c": _daily(1.0, edit=dry)})
    out = monthly_leave_one_out(st, 9, GAGES, [("a", "b"), ("c", "a")])
    # c is excluded in September; a and b each keep one pair, which leaves no pairs to
    # refit s2_d from when either is held out.
    assert "c" not in set(out["site_no"])
    assert out["note"].str.contains("s2_d refit failed").all()
    assert out["blend"].isna().all()
