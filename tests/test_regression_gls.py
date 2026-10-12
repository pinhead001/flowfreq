"""Tests for flowfreq.regression.gls (WREG-equivalent GLS for flood quantiles)."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import norm

from flowfreq.regression.gls import (
    GLSResult,
    WregCorrelation,
    fit_correlation,
    fit_gls,
    lp3_quantile_sampling_variance,
    sampling_covariance,
)


def _design(n: int, rng: np.random.Generator) -> np.ndarray:
    return np.column_stack([np.ones(n), rng.normal(2.0, 0.6, n), rng.normal(0.0, 1.0, n)])


def test_zero_sampling_error_equals_ols():
    rng = np.random.default_rng(1)
    X = _design(60, rng)
    y = X @ np.array([1.0, 0.8, -0.2]) + rng.normal(0, 0.15, 60)
    res = fit_gls(y, X, np.zeros((60, 60)), names=["const", "logA", "P"])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    np.testing.assert_allclose(res.beta, beta, atol=1e-10)
    sse = float(np.sum((y - X @ beta) ** 2))
    assert res.model_error_variance == pytest.approx(sse / (60 - 3), abs=1e-10)
    assert isinstance(res, GLSResult) and res.names == ["const", "logA", "P"]


def test_known_model_error_variance_recovered():
    # Mean of 10 independent fits (n = 800 each): the estimator's SE is about
    # 0.0034 per fit, 0.0011 for the mean, so 0.003 is ~3 SE.
    mevs = []
    for seed in range(10):
        rng = np.random.default_rng(seed)
        n = 800
        X = _design(n, rng)
        samp = rng.uniform(0.005, 0.05, n)
        y = X @ np.array([2.0, 0.7, 0.1]) + rng.normal(0, 0.2, n) + rng.normal(0, np.sqrt(samp))
        res = fit_gls(y, X, np.diag(samp))
        # AVP = MEV + average sampling variance of the prediction; always >= MEV.
        assert res.avp > res.model_error_variance
        assert 0.0 < res.pseudo_r2 <= 1.0
        mevs.append(res.model_error_variance)
    assert float(np.mean(mevs)) == pytest.approx(0.04, abs=0.003)


def test_model_error_variance_floors_at_zero():
    # Residuals much smaller than the stated sampling error: no positive root.
    rng = np.random.default_rng(3)
    X = _design(40, rng)
    y = X @ np.array([1.0, 0.5, 0.0]) + rng.normal(0, 0.01, 40)
    res = fit_gls(y, X, np.eye(40) * 0.5)
    assert res.model_error_variance == 0.0


def test_predict_returns_variance_including_mev():
    rng = np.random.default_rng(5)
    X = _design(80, rng)
    y = X @ np.array([1.0, 0.5, 0.3]) + rng.normal(0, 0.1, 80)
    res = fit_gls(y, X, np.eye(80) * 0.001)
    yhat, var = res.predict(np.array([1.0, 2.0, 0.0]))
    assert yhat == pytest.approx(1.0 + 1.0, abs=0.1)
    assert var > res.model_error_variance


def test_singular_design_raises():
    rng = np.random.default_rng(2)
    X = _design(30, rng)
    X = np.column_stack([X, X[:, 1]])
    with pytest.raises(ValueError, match="rank"):
        fit_gls(rng.normal(size=30), X, np.zeros((30, 30)))


def test_fewer_sites_than_parameters_raises():
    X = np.ones((3, 3))
    X[:, 1] = [1, 2, 3]
    X[:, 2] = [0, 1, 0]
    with pytest.raises(ValueError, match="sites"):
        fit_gls(np.array([1.0, 2.0, 3.0]), X, np.zeros((3, 3)))


def test_mismatched_shapes_raise():
    with pytest.raises(ValueError, match="shape"):
        fit_gls(np.zeros(10), np.ones((9, 1)), np.zeros((10, 10)))


def test_variance_at_zero_skew_closed_form():
    k = norm.ppf(0.99)
    expected = 0.3**2 / 50 * (1 + k**2 / 2)
    got = lp3_quantile_sampling_variance(0.3, 0.0, 50, 0.01, 0.0)
    assert float(got) == pytest.approx(expected, rel=1e-8)


def test_skew_mse_term_increases_variance():
    a = float(lp3_quantile_sampling_variance(0.3, 0.2, 50, 0.01, 0.0))
    b = float(lp3_quantile_sampling_variance(0.3, 0.2, 50, 0.01, 0.12))
    assert b > a


def test_variance_rejects_bad_inputs():
    with pytest.raises(ValueError, match="aep"):
        lp3_quantile_sampling_variance(0.3, 0.0, 50, 1.5, 0.0)
    with pytest.raises(ValueError, match="record length"):
        lp3_quantile_sampling_variance(0.3, 0.0, 0, 0.01, 0.0)


def test_wreg_correlation_model():
    m = WregCorrelation(alpha=0.002, theta=0.98)
    r = m.rho(np.array([0.0, 50.0, 500.0]))
    assert r[0] == pytest.approx(1.0)
    assert 1.0 > r[1] > r[2] > 0.0


def test_fit_correlation_recovers_parameters():
    true = WregCorrelation(alpha=0.004, theta=0.985)
    d = np.linspace(1, 600, 400)
    r = true.rho(d) + np.random.default_rng(0).normal(0, 0.01, d.size)
    fit = fit_correlation(d, r)
    assert fit.theta == pytest.approx(true.theta, abs=0.005)
    assert fit.rho(np.array([100.0]))[0] == pytest.approx(true.rho(np.array([100.0]))[0], abs=0.03)


def test_sampling_covariance_zero_without_concurrent_record():
    sd = np.array([0.3, 0.25, 0.28])
    skew = np.array([0.0, 0.1, -0.1])
    n = np.array([40.0, 30.0, 50.0])
    conc = np.array([[40, 0, 20], [0, 30, 10], [20, 10, 50]], dtype=float)
    lat = np.array([47.0, 47.1, 47.2])
    lon = np.array([-120.0, -120.1, -120.2])
    S = sampling_covariance(sd, skew, n, 0.01, 0.18, lat, lon, conc, WregCorrelation(0.002, 0.98))
    assert S[0, 1] == 0.0 and S[1, 0] == 0.0
    assert S[0, 2] > 0.0 and np.allclose(S, S.T)
    np.testing.assert_allclose(np.diag(S), lp3_quantile_sampling_variance(sd, skew, n, 0.01, 0.18))
