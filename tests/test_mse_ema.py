"""Tests for hydrolib._mse_ema that do not require the Fortran extension.

Parity against the vendored Fortran lives in
tests/fortran_parity/test_fortran_oracles.py (skipped when the extension
isn't built). These check properties and closed-form cross-checks that
should hold regardless.
"""

from __future__ import annotations

import numpy as np
import pytest

from hydrolib._mse_ema import _chol33, _covw, _gammaquad, _normquad, mc2mnvb, mn2m_var, mse_ema


class TestQuadrature:
    def test_normquad_weights_sum_to_one_and_are_symmetric(self):
        t, w = _normquad(5)
        assert w.sum() == pytest.approx(1.0)
        assert np.allclose(sorted(t), sorted(-t))  # symmetric about 0
        assert np.average(t, weights=w) == pytest.approx(0.0, abs=1e-9)

    def test_normquad_reproduces_the_standard_normal_variance(self):
        t, w = _normquad(5)
        assert np.average(t**2, weights=w) == pytest.approx(1.0, rel=1e-6)

    def test_gammaquad_weights_sum_to_one(self):
        t, w = _gammaquad(4, alpha=6.0, beta=0.5)
        assert w.sum() == pytest.approx(1.0)

    def test_gammaquad_mean_matches_alpha_beta(self):
        t, w = _gammaquad(4, alpha=6.0, beta=0.5)
        assert np.average(t, weights=w) == pytest.approx(6.0 * 0.5, rel=1e-6)

    def test_gammaquad_nodes_are_nonnegative(self):
        """A Gamma distribution has support [0, infinity)."""
        t, _ = _gammaquad(4, alpha=6.0, beta=0.5)
        assert np.all(t >= 0.0)


class TestChol33:
    def test_reconstructs_the_covariance(self):
        s = np.array([[2.0, 0.3, -0.1], [0.3, 1.5, 0.2], [-0.1, 0.2, 0.8]])
        v = _chol33(s)
        assert v is not None
        assert np.allclose(v.T @ v, s)

    def test_has_the_documented_zero_pattern(self):
        s = np.array([[2.0, 0.3, -0.1], [0.3, 1.5, 0.2], [-0.1, 0.2, 0.8]])
        v = _chol33(s)
        assert v[0, 1] == 0.0
        assert v[2, 0] == 0.0
        assert v[2, 1] == 0.0

    def test_returns_none_for_non_positive_diagonal(self):
        s = np.array([[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        assert _chol33(s) is None

    def test_returns_none_for_non_positive_semi_definite(self):
        s = np.array([[1.0, 5.0, 0.0], [5.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        assert _chol33(s) is None


class TestCovw:
    def test_matches_plain_weighted_covariance_with_equal_weights(self):
        rng = np.random.default_rng(0)
        x = rng.normal(size=20)
        y = rng.normal(size=20)
        w = np.ones(20)
        assert _covw(x, y, w) == pytest.approx(np.cov(x, y, ddof=0)[0, 1])

    def test_variance_is_nonnegative(self):
        rng = np.random.default_rng(1)
        x = rng.normal(size=10)
        w = rng.uniform(0.1, 1.0, size=10)
        assert _covw(x, x, w) >= 0.0


class TestMc2mnvb:
    def test_returns_symmetric_matrix(self):
        mc = np.array([3.7, 0.09, 0.3])
        s_mc = np.array(
            [[0.001, 0.0001, 0.0002], [0.0001, 0.0002, 0.00005], [0.0002, 0.00005, 0.001]]
        )
        s_mn = mc2mnvb(mc, s_mc)
        assert np.allclose(s_mn, s_mn.T)

    def test_rejects_non_positive_semi_definite_covariance(self):
        mc = np.array([3.7, 0.09, 0.3])
        s_mc = np.array([[1.0, 5.0, 0.0], [5.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        with pytest.raises(ValueError):
            mc2mnvb(mc, s_mc)


class TestMn2mVar:
    def test_matches_m2mn_point_estimate(self):
        """mn2m_var's mc must be the exact mn2m(mn) -- it only linearizes

        the covariance, not the point estimate.
        """
        from hydrolib._p3_moments import mn2m

        mn = np.array([3.71, 13.8, 51.6])
        s_mn = np.eye(3) * 0.001
        mc, _ = mn2m_var(mn, s_mn)
        assert np.allclose(mc, mn2m(mn))

    def test_covariance_is_symmetric(self):
        mn = np.array([3.71, 13.8, 51.6])
        s_mn = np.array(
            [[0.001, 0.0002, 0.0004], [0.0002, 0.0005, 0.0003], [0.0004, 0.0003, 0.002]]
        )
        _, s_mc = mn2m_var(mn, s_mn)
        assert np.allclose(s_mc, s_mc.T)


class TestMseEma:
    def test_rejects_kmom_out_of_range(self):
        mc = np.array([3.7, 0.09, 0.3])
        with pytest.raises(ValueError):
            mse_ema(np.array([50.0]), np.array([-20.0]), np.array([20.0]), mc, kmom=4)
        with pytest.raises(ValueError):
            mse_ema(np.array([50.0]), np.array([-20.0]), np.array([20.0]), mc, kmom=0)

    def test_is_positive(self):
        mc = np.array([3.7, 0.09, 0.3])
        for kmom in (1, 2, 3):
            assert mse_ema(np.array([50.0]), np.array([-20.0]), np.array([20.0]), mc, kmom) > 0.0

    def test_more_observations_reduces_the_mse(self):
        mc = np.array([3.7, 0.09, 0.3])
        small = mse_ema(np.array([20.0]), np.array([-20.0]), np.array([20.0]), mc, kmom=3)
        large = mse_ema(np.array([200.0]), np.array([-20.0]), np.array([20.0]), mc, kmom=3)
        assert large < small

    def test_blends_continuously_across_skew_zero(self):
        """The skewmin blend (|skew| <= 0.06324555) should not produce a

        visible jump right at the boundary.
        """
        mc_at_floor = np.array([3.7, 0.09, 0.06324555])
        mc_just_inside = np.array([3.7, 0.09, 0.06])
        just_outside = mse_ema(
            np.array([50.0]), np.array([-20.0]), np.array([20.0]), mc_at_floor, 3
        )
        just_inside = mse_ema(
            np.array([50.0]), np.array([-20.0]), np.array([20.0]), mc_just_inside, 3
        )
        assert just_outside == pytest.approx(just_inside, rel=0.05)
