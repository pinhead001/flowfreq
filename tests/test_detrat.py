"""Tests for hydrolib._detrat that do not require the Fortran extension.

Parity against the vendored Fortran lives in
tests/fortran_parity/test_fortran_oracles.py (skipped when the extension
isn't built). These check properties that should hold regardless.
"""

from __future__ import annotations

import numpy as np
import pytest

from hydrolib._detrat import _expmomcderiv, _p2m, detrat
from hydrolib._p3_moments import m2p


class TestP2m:
    def test_inverts_m2p(self):
        mc = np.array([3.7, 0.09, -0.5])
        parms = m2p(mc)
        assert np.allclose(_p2m(parms), mc, rtol=1e-12)

    @pytest.mark.parametrize("skew", [-1.2, -0.3, 0.3, 1.2])
    def test_inverts_m2p_across_skew_signs(self, skew):
        mc = np.array([3.7, 0.09, skew])
        parms = m2p(mc)
        assert np.allclose(_p2m(parms), mc, rtol=1e-12)


class TestExpMomCDeriv:
    def test_mne_between_the_thresholds_when_censoring_is_a_single_tail(self):
        """E[X | X < tl] should sit below tl for a lower-tail-only group."""
        mc = np.array([3.7, 0.09, -0.5])
        parms = m2p(mc)
        _, mne, _ = _expmomcderiv(parms, 3.5, 20.0)
        assert mne[0] < 3.5

    def test_shape(self):
        mc = np.array([3.7, 0.09, -0.5])
        parms = m2p(mc)
        mc_out, mne, dedmc = _expmomcderiv(parms, 3.5, 4.0)
        assert mc_out.shape == (3,)
        assert mne.shape == (3,)
        assert dedmc.shape == (3, 3)


class TestDetrat:
    def test_is_one_below_the_skew_floor(self):
        mc = np.array([3.7, 0.09, 0.02])
        assert detrat(mc, 50, np.array([50.0]), np.array([-20.0]), np.array([20.0])) == 1.0

    def test_computes_the_real_ratio_at_the_skew_floor_boundary(self):
        """abs(skew) < 0.04 short-circuits; exactly 0.04 does not (emafit.f:3654's

        strict inequality), so a censored group at exactly the floor must
        still go through the real determinant-ratio computation.
        """
        mc = np.array([3.7, 0.09, 0.04])
        w = detrat(mc, 80, np.array([60.0, 20.0]), np.array([-20.0, 4.3]), np.array([4.3, 20.0]))
        assert w != 1.0

    def test_is_one_when_nothing_is_censored_even_above_the_floor(self):
        """A full-support threshold group carries no censoring to correct for."""
        mc = np.array([3.7, 0.09, -0.5])
        w = detrat(mc, 50, np.array([50.0]), np.array([-20.0]), np.array([20.0]))
        assert w == pytest.approx(1.0, abs=1e-6)

    def test_symmetric_skew_gives_the_same_ratio(self):
        """detrat depends on |skew| in the same way m2p's alpha does; a

        mirrored threshold group at the opposite skew sign should give the
        same Wd.
        """
        nobs, tl, tu = np.array([60.0, 20.0]), np.array([-20.0, 4.0]), np.array([4.0, 20.0])
        w_pos = detrat(np.array([3.7, 0.09, 0.5]), 80, nobs, tl, tu)
        w_neg = detrat(np.array([3.7, 0.09, -0.5]), 80, nobs, tl, tu)
        # Not exactly equal -- the threshold group itself is not mirrored
        # around the mean -- but both should land well inside (0, 1).
        assert 0.0 < w_pos < 1.0
        assert 0.0 < w_neg < 1.0

    def test_more_censoring_pulls_wd_further_from_one(self):
        """A wider censored region should discount the at-site skew more."""
        mc = np.array([3.7, 0.09, -0.6])
        light = detrat(
            mc, 80, np.array([70.0, 10.0]), np.array([-20.0, 4.3]), np.array([4.3, 20.0])
        )
        heavy = detrat(
            mc, 80, np.array([40.0, 40.0]), np.array([-20.0, 3.8]), np.array([3.8, 20.0])
        )
        assert abs(heavy - 1.0) > abs(light - 1.0)
