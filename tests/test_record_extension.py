"""Tests for flowfreq.record_extension (MOVE.1)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from flowfreq.record_extension import MIN_CONCURRENT, Move1, extend_record, fit_move1


def _pair(n: int = 400, r_noise: float = 0.3, seed: int = 0):
    rng = np.random.default_rng(seed)
    dates = pd.date_range("2000-01-01", periods=n, freq="D")
    lx = rng.normal(2.0, 0.5, n)
    ly = 1.0 + 0.8 * (lx - 2.0) + rng.normal(0.0, r_noise, n)
    return pd.Series(10**ly, index=dates), pd.Series(10**lx, index=dates)


class TestFitMove1:
    def test_exact_relation_is_recovered(self):
        dates = pd.date_range("2000-01-01", periods=50, freq="D")
        x = pd.Series(np.linspace(10, 1000, 50), index=dates)
        y = 0.3 * x**1.2
        fit = fit_move1(y, x)
        assert fit.slope == pytest.approx(1.2)
        assert fit.correlation == pytest.approx(1.0)
        np.testing.assert_allclose(fit.predict(x.to_numpy()), y.to_numpy(), rtol=1e-10)

    def test_preserves_mean_and_variance_unlike_ols(self):
        y, x = _pair()
        fit = fit_move1(y, x)
        ly = np.log10(fit.predict(x.to_numpy()))
        assert ly.mean() == pytest.approx(np.log10(y).mean())
        assert ly.std(ddof=1) == pytest.approx(np.log10(y).std(ddof=1))
        # OLS would shrink the variance by r; MOVE.1 does not.
        assert fit.correlation < 0.9
        assert fit.slope == pytest.approx(np.log10(y).std(ddof=1) / np.log10(x).std(ddof=1))

    def test_nonpositive_pairs_are_excluded_and_counted(self):
        y, x = _pair()
        y.iloc[:5] = 0.0
        x.iloc[5:8] = -1.0
        fit = fit_move1(y, x)
        assert fit.n_excluded_nonpositive == 8
        assert fit.n_concurrent == len(y) - 8

    def test_only_concurrent_dates_are_paired(self):
        y, x = _pair()
        fit = fit_move1(y.iloc[100:200], x)
        assert fit.n_concurrent == 100
        assert fit.concurrent_start == y.index[100]
        assert fit.concurrent_end == y.index[199]

    def test_too_few_pairs_raise(self):
        y, x = _pair(n=MIN_CONCURRENT - 1)
        with pytest.raises(ValueError, match="concurrent positive pairs"):
            fit_move1(y, x)

    def test_negative_correlation_raises(self):
        y, x = _pair()
        with pytest.raises(ValueError, match="positively correlated"):
            fit_move1(pd.Series(1.0 / y.to_numpy(), index=y.index), x)

    def test_constant_record_raises(self):
        y, x = _pair()
        with pytest.raises(ValueError, match="variation"):
            fit_move1(pd.Series(5.0, index=y.index), x)

    def test_non_series_input_raises(self):
        y, x = _pair()
        with pytest.raises(TypeError):
            fit_move1(y.to_numpy(), x)
        with pytest.raises(TypeError, match="DatetimeIndex"):
            fit_move1(y.reset_index(drop=True), x)


class TestPredict:
    def test_zero_and_missing_index_give_nan(self):
        y, x = _pair()
        out = fit_move1(y, x).predict(np.array([0.0, -2.0, np.nan, 100.0]))
        assert np.isnan(out[:3]).all()
        assert out[3] > 0


class TestExtendRecord:
    def test_keeps_observed_and_fills_the_rest(self):
        y, x = _pair(n=600)
        short = y.iloc[200:300]
        out = extend_record(short, x)
        assert len(out) == len(x)
        assert (out.loc[short.index, "source"] == "observed").all()
        np.testing.assert_allclose(out.loc[short.index, "flow"], short.to_numpy())
        assert (out["source"] == "move1").sum() == len(x) - len(short)

    def test_index_zero_gives_no_estimate(self):
        y, x = _pair(n=600)
        x.iloc[0] = 0.0
        out = extend_record(y.iloc[200:300], x)
        assert out["source"].iloc[0] == "none"
        assert np.isnan(out["flow"].iloc[0])

    def test_supplied_fit_is_used(self):
        y, x = _pair(n=600)
        fit = Move1(0.0, 1.0, 0.0, 1.0, 1.0, 10, 0, x.index[0], x.index[9])
        out = extend_record(y.iloc[:10], x, fit=fit)
        np.testing.assert_allclose(out["flow"].iloc[10:], x.iloc[10:].to_numpy())

    def test_duplicate_dates_raise(self):
        y, x = _pair()
        with pytest.raises(ValueError, match="duplicate"):
            extend_record(pd.concat([y, y.iloc[:1]]), x)
