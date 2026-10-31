"""Stationarity screening (roadmap section 6.1; methow item 8) and the mixed-population
seasonality screen (methow item 9).

Reference values are computed by hand from the textbook formulas (Mann 1945, Kendall 1975,
Sen 1968, Pettitt 1979, Hamed and Rao 1998) or cross-checked against scipy.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from flowfreq.stationarity import (
    mann_kendall,
    mixed_population_screen,
    pettitt,
    sens_slope,
    spearman_trend,
    stationarity_screen,
)


def test_mann_kendall_increasing_series_by_hand() -> None:
    r = mann_kendall([1.0, 2.0, 3.0, 4.0, 5.0])
    assert r.s == 10
    assert r.var_s == pytest.approx(5 * 4 * 15 / 18)
    z = (10 - 1) / math.sqrt(5 * 4 * 15 / 18)
    assert r.z == pytest.approx(z)
    assert r.p_value == pytest.approx(2 * (1 - stats.norm.cdf(z)))
    assert r.tau == pytest.approx(1.0)
    assert r.trend == "increasing"


def test_mann_kendall_tie_correction() -> None:
    x = [1.0, 2.0, 2.0, 3.0, 3.0, 3.0]
    r = mann_kendall(x)
    n = 6
    ties = (2 * 1 * 9) + (3 * 2 * 11)  # sum t (t-1) (2t+5) over tie groups of 2 and 3
    assert r.var_s == pytest.approx((n * (n - 1) * (2 * n + 5) - ties) / 18)
    # Without ties the statistic agrees with scipy's Kendall tau.
    y = np.random.default_rng(1).normal(size=30)
    tau, _ = stats.kendalltau(np.arange(30), y)
    assert mann_kendall(y).tau == pytest.approx(tau)


def test_mann_kendall_no_trend_and_decreasing() -> None:
    assert mann_kendall([3.0, 1.0, 3.0, 1.0, 3.0, 1.0]).trend == "no trend"
    assert mann_kendall(np.arange(20.0)[::-1]).trend == "decreasing"


def test_hamed_rao_inflates_the_variance_for_autocorrelated_series() -> None:
    rng = np.random.default_rng(4)
    e = rng.normal(size=200)
    x = np.zeros(200)
    for t in range(1, 200):
        x[t] = 0.7 * x[t - 1] + e[t]  # AR(1), no trend
    plain = mann_kendall(x)
    corrected = mann_kendall(x, method="hamed_rao")
    assert corrected.var_s > plain.var_s
    assert corrected.p_value > plain.p_value
    white = rng.normal(size=200)
    ratio = mann_kendall(white, method="hamed_rao").var_s / mann_kendall(white).var_s
    assert 0.7 < ratio < 1.5


def test_sens_slope_by_hand() -> None:
    years = np.array([2000, 2001, 2002, 2003, 2004])
    x = np.array([1.0, 3.0, 2.0, 5.0, 4.0])
    slopes = [(x[j] - x[i]) / (years[j] - years[i]) for i in range(5) for j in range(i + 1, 5)]
    slope, intercept = sens_slope(x, years)
    assert slope == pytest.approx(np.median(slopes))
    assert intercept == pytest.approx(np.median(x - slope * years))


def test_spearman_matches_scipy() -> None:
    y = np.random.default_rng(2).normal(size=25) + 0.05 * np.arange(25)
    rho, p = spearman_trend(y)
    ref = stats.spearmanr(np.arange(25), y)
    assert rho == pytest.approx(ref.statistic) and p == pytest.approx(ref.pvalue)


def test_pettitt_step_change_by_hand() -> None:
    x = [1.0] * 5 + [5.0] * 5
    r = pettitt(x, years=range(1990, 2000))
    assert r.index == 5 and r.year == 1995  # first year of the second segment
    assert r.k == 25  # t (n - t) for a clean step
    assert r.p_value == pytest.approx(2 * math.exp(-6 * 25**2 / (10**3 + 10**2)))
    assert r.shift == "up"


def test_screen_summarizes_and_flags() -> None:
    years = np.arange(1960, 2020)
    rng = np.random.default_rng(5)
    flat = pd.Series(np.exp(rng.normal(8, 0.4, 60)), index=years)
    out = stationarity_screen(flat)
    assert set(out) >= {"n", "mann_kendall", "sens_slope_log10_per_year", "spearman", "pettitt",
                        "flags"}  # fmt: skip
    assert out["n"] == 60 and out["flags"] == []
    stepped = flat.copy()
    stepped.loc[1990:] *= 3.0
    flagged = stationarity_screen(stepped)
    assert any("change point" in f for f in flagged["flags"])
    assert flagged["pettitt"]["year"] in range(1985, 1996)


@pytest.mark.parametrize("bad", [[1.0, 2.0], [1.0, float("nan"), 2.0, 3.0]])
def test_screen_refuses_short_or_missing_series(bad) -> None:
    with pytest.raises(ValueError):
        mann_kendall(bad)


def test_mixed_population_screen() -> None:
    snow = pd.to_datetime([f"{y}-05-20" for y in range(1980, 2010)])
    out = mixed_population_screen(snow)
    assert out["cool_season_share"] == 0.0 and not out["mixed"]
    mixed_dates = pd.to_datetime(
        [f"{y}-05-20" for y in range(1980, 1995)] + [f"{y}-11-15" for y in range(1995, 2010)]
    )
    out = mixed_population_screen(mixed_dates)
    assert out["cool_season_share"] == pytest.approx(0.5)
    assert out["mixed"] and "mixed" in out["message"]
    # The circular mean date's concentration drops when peaks split between seasons.
    assert mixed_population_screen(snow)["concentration"] > out["concentration"]
