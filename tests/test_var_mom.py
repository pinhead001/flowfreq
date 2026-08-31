"""Tests for hydrolib._var_mom that do not require the Fortran extension.

Parity against the vendored Fortran lives in
tests/fortran_parity/test_fortran_oracles.py (skipped when the extension
isn't built). These check properties and closed-form cross-checks that
should hold regardless.
"""

from __future__ import annotations

import numpy as np
import pytest

from hydrolib._p3_moments import m2p, m_p3
from hydrolib._var_mom import d_est, expmomderiv, var_mom, varb, varc


class TestVarb:
    def test_zero_when_no_split_uncertainty(self):
        """p1 = p3 = 0 (all mass in the middle region) leaves no multinomial

        variance to propagate -- vn is identically zero.
        """
        mu_x = np.array([[1.0, 3.7, 1.0], [1.0, 14.0, 1.0], [1.0, 55.0, 1.0]])
        vb = varb(mu_x, 50.0, 0.0, 1.0, 0.0)
        assert np.allclose(vb, 0.0)

    def test_symmetric(self):
        mu_x = np.array([[1.0, 3.7, 6.0], [2.0, 14.0, 40.0], [4.0, 55.0, 260.0]])
        vb = varb(mu_x, 50.0, 0.2, 0.5, 0.3)
        assert np.allclose(vb, vb.T)


class TestVarc:
    def test_symmetric(self):
        e_x = np.array([3.7, 14.0, 55.0, 230.0, 1000.0, 4600.0])
        vc = varc(e_x, 50.0, 0.6)
        assert np.allclose(vc, vc.T)

    def test_diagonal_matches_variance_formula(self):
        """vc(i,i) = nh*p2*(E[X^2i] - E[X^i]**2): the ordinary variance of

        the i-th raw moment, scaled by the number of observations actually
        in the interval.
        """
        e_x = np.array([3.7, 14.0, 55.0, 230.0, 1000.0, 4600.0])
        nh, p2 = 50.0, 0.6
        vc = varc(e_x, nh, p2)
        mu_n = nh * p2
        for i in range(1, 4):
            expected = mu_n * (e_x[2 * i - 1] - e_x[i - 1] ** 2)
            assert vc[i - 1, i - 1] == pytest.approx(expected)


class TestDEst:
    def test_zero_below_probability_cutoff(self):
        mc = np.array([0.0, 1.0, 0.3])
        d = d_est(50.0, mc, pa=0.0, pb=1.0)
        assert np.allclose(d, 0.0)

    def test_nonzero_when_a_tail_has_real_probability(self):
        mc = np.array([0.0, 1.0, 0.3])
        d = d_est(50.0, mc, pa=0.3, pb=1.0)
        assert not np.allclose(d, 0.0)


class TestExpMomDeriv:
    def test_shape(self):
        mc = np.array([0.0, 1.0, 0.3])
        parms = m2p(mc)
        jac = expmomderiv(parms, -1e19, 2.0)
        assert jac.shape == (3, 3)

    def test_zero_for_an_empty_interval(self):
        mc = np.array([0.0, 1.0, 0.3])
        parms = m2p(mc)
        assert np.allclose(expmomderiv(parms, 5.0, 5.0), 0.0)
        assert np.allclose(expmomderiv(parms, 5.0, 4.0), 0.0)


class TestVarMom:
    def test_returns_symmetric_positive_definite_matrix(self):
        n = np.array([85.0])
        tl, tu = np.array([-20.0]), np.array([20.0])
        mc = np.array([3.7, 0.09, -0.2])
        v = var_mom(n, tl, tu, mc)
        assert v.shape == (3, 3)
        assert np.allclose(v, v.T)
        assert np.all(np.linalg.eigvalsh(v) > 0)

    def test_uncensored_case_matches_the_classical_moment_covariance(self):
        """With nothing censored, d_est and varb both vanish (see TestVarb),

        so var_mom reduces to Cov(raw moments)/n -- the ordinary
        delta-method covariance of moment estimators from n iid draws,
        computable independently from m_p3's own raw moments without going
        through var_mom's threshold-group machinery at all.
        """
        mc = np.array([3.7, 0.09, -0.2])
        n = 85.0
        v = var_mom(np.array([n]), np.array([-20.0]), np.array([20.0]), mc)

        # var_mom itself always works in mean-centred coordinates (mean 0) --
        # see hydrolib._var_mom.var_mom -- so the independent check must too.
        mc_centred = np.array([0.0, mc[1], mc[2]])
        e_x = m_p3(-999.0, 999.0, mc_centred, 6)
        expected = np.zeros((3, 3))
        for i in range(1, 4):
            for j in range(1, 4):
                expected[i - 1, j - 1] = (e_x[i + j - 1] - e_x[i - 1] * e_x[j - 1]) / n
        assert np.allclose(v, expected, rtol=1e-9)

    def test_rejects_inverted_threshold_pair(self):
        mc = np.array([3.7, 0.09, -0.2])
        with pytest.raises(ValueError):
            var_mom(np.array([50.0]), np.array([5.0]), np.array([2.0]), mc)

    def test_more_observations_shrinks_the_covariance(self):
        mc = np.array([3.7, 0.09, -0.2])
        v_small = var_mom(np.array([20.0]), np.array([-20.0]), np.array([20.0]), mc)
        v_large = var_mom(np.array([200.0]), np.array([-20.0]), np.array([20.0]), mc)
        assert np.all(np.diag(v_large) < np.diag(v_small))
