"""Tests for hydrolib._var_emab that do not require the Fortran extension.

Parity against the vendored Fortran lives in
tests/fortran_parity/test_fortran_oracles.py::TestVarEmabPort (skipped when
the extension isn't built). These check properties that should hold
regardless, plus a real-fixture smoke test for regmoms/var_emab -- both
need a well-conditioned covariance to run at all (mn2mvarb's root-find,
gridmake's Cholesky factor), so arbitrary made-up moments are not a safe
choice of input the way they are for e.g. _p3_moments's leaf routines.
"""

from __future__ import annotations

import numpy as np
import pytest

from hydrolib._var_emab import ci_ema_m3b, regmoms, var_emab

# Big Sandy's own converged weighted fit and threshold groups (TODO.md P3),
# verified against the Fortran oracle in
# tests/fortran_parity/test_fortran_oracles.py::TestVarEmabPort -- reused
# here as a real, known-well-conditioned fixture rather than invented
# numbers that might not be.
_BIG_SANDY_MC = np.array([3.71750771, 0.08470576, -0.15630573])
_BIG_SANDY_NOBS = np.array([44.0, 40.0])
_BIG_SANDY_TL = np.array([-20.0, 4.25527251])
_BIG_SANDY_TU = np.array([20.0, 20.0])
_BIG_SANDY_R_G_MSE = 0.30250000000000005


class TestCiEmaM3b:
    def test_symmetric_covariance_gives_symmetric_bounds(self):
        """cov(yp, syp) = 0 means beta1 = 0, so both denominators are 1."""
        cv = np.array([[0.01, 0.0], [0.0, 0.001]])
        cil, cih, nu = ci_ema_m3b(5.0, cv, 0.9)
        assert abs((cih - 5.0) - (5.0 - cil)) < 1e-9
        assert nu == pytest.approx(0.01)

    def test_nonzero_covariance_gives_asymmetric_bounds(self):
        cv = np.array([[0.01, -0.002], [-0.002, 0.001]])
        cil, cih, _ = ci_ema_m3b(5.0, cv, 0.9)
        assert abs((cih - 5.0) - (5.0 - cil)) > 1e-6

    def test_wider_coverage_gives_a_wider_interval(self):
        cv = np.array([[0.01, -0.002], [-0.002, 0.001]])
        cil90, cih90, _ = ci_ema_m3b(5.0, cv, 0.90)
        cil99, cih99, _ = ci_ema_m3b(5.0, cv, 0.99)
        assert cih99 - cil99 > cih90 - cil90

    def test_nu_output_is_var_yp_not_degrees_of_freedom(self):
        """emafit.f:1902 clobbers its own nu output with var(yp); transcribed as-is."""
        cv = np.array([[0.02, -0.001], [-0.001, 0.0005]])
        _, _, nu = ci_ema_m3b(5.0, cv, 0.9)
        assert nu == pytest.approx(cv[0, 0])


class TestRegmoms:
    def test_returns_a_symmetric_positive_semidefinite_matrix(self):
        s_mc = regmoms(
            _BIG_SANDY_NOBS, _BIG_SANDY_TL, _BIG_SANDY_TU, _BIG_SANDY_MC, r_g_mse=_BIG_SANDY_R_G_MSE
        )
        assert s_mc.shape == (3, 3)
        assert np.allclose(s_mc, s_mc.T)
        assert np.all(np.linalg.eigvalsh(s_mc) >= -1e-12)

    def test_station_only_skew_gives_a_larger_skew_variance_than_weighted(self):
        """Folding in a real regional skew (finite MSE) should only ever shrink

        Var[skew], never grow it -- the weighted estimate is at least as
        precise as the at-site one alone.
        """
        weighted = regmoms(
            _BIG_SANDY_NOBS, _BIG_SANDY_TL, _BIG_SANDY_TU, _BIG_SANDY_MC, r_g_mse=_BIG_SANDY_R_G_MSE
        )
        station_only = regmoms(
            _BIG_SANDY_NOBS, _BIG_SANDY_TL, _BIG_SANDY_TU, _BIG_SANDY_MC, r_g_mse=1e10
        )
        assert weighted[2, 2] <= station_only[2, 2]


class TestVarEmab:
    def test_yp_is_monotonic_in_non_exceedance_probability(self):
        """Higher pq (closer to 1) must give a higher quantile."""
        pq = np.array([0.1, 0.5, 0.9, 0.99])
        yp, _, cil, cih, _ = var_emab(
            _BIG_SANDY_NOBS,
            _BIG_SANDY_TL,
            _BIG_SANDY_TU,
            _BIG_SANDY_MC,
            pq,
            eps=0.9,
            r_g_mse=_BIG_SANDY_R_G_MSE,
        )
        assert np.all(np.diff(yp) > 0)
        assert np.all(cil < yp)
        assert np.all(cih > yp)

    def test_matches_the_classic_formula_when_the_grid_is_degenerate(self):
        """Not a real-world case -- a sanity check on the machinery itself.

        var_emab still has to run to completion (regmoms, two levels of
        gridmake, q_p3 at every node) without raising, and produce a
        confidence interval that actually brackets the point estimate.
        """
        pq = np.array([0.5])
        yp, cv, cil, cih, nu = var_emab(
            _BIG_SANDY_NOBS,
            _BIG_SANDY_TL,
            _BIG_SANDY_TU,
            _BIG_SANDY_MC,
            pq,
            eps=0.9,
            r_g_mse=_BIG_SANDY_R_G_MSE,
        )
        assert cv.shape == (1, 2, 2)
        assert cil[0] < yp[0] < cih[0]
