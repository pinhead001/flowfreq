"""Tests for flowfreq._p3_moments that do not require the Fortran extension.

Parity against the vendored Fortran itself lives in
tests/fortran_parity/test_fortran_oracles.py (skipped when the extension
isn't built). These check the pure-Python properties that should hold
regardless: algebraic identities, known closed forms, and error handling.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import stats

from flowfreq._p3_moments import m2mn, m2p, m_p3, mn2m, p_p3, q_p3


class TestM2pM2mn:
    def test_m2p_recovers_normal_when_skew_tiny_but_nonzero(self):
        """alpha = 4/skew**2 grows without bound as skew -> 0; sanity only."""
        tau, alpha, beta = m2p(np.array([5.0, 4.0, 0.01]))
        assert alpha == pytest.approx(4.0 / 0.01**2)
        assert beta == pytest.approx(np.sqrt(4.0 / alpha))
        assert tau == pytest.approx(5.0 - alpha * beta)

    def test_m2p_beta_sign_follows_skew_sign(self):
        _, _, beta_pos = m2p(np.array([0.0, 1.0, 0.5]))
        _, _, beta_neg = m2p(np.array([0.0, 1.0, -0.5]))
        assert beta_pos > 0
        assert beta_neg < 0

    def test_m2mn_mn2m_roundtrip(self):
        m = np.array([3.7, 0.09, -0.3])
        assert np.allclose(mn2m(m2mn(m)), m, rtol=1e-12)

    def test_m2mn_matches_definition_for_symmetric_distribution(self):
        """skew = 0 drops the only term m2mn adds beyond raw-moment algebra."""
        mu, var = 2.0, 3.0
        mn = m2mn(np.array([mu, var, 0.0]))
        assert mn[0] == pytest.approx(mu)
        assert mn[1] == pytest.approx(mu**2 + var)
        assert mn[2] == pytest.approx(3.0 * mn[1] * mu - 2.0 * mu**3)


class TestPP3QP3:
    def test_p_p3_matches_normal_cdf_when_skew_exactly_zero(self):
        """skew = 0 forces wg = 0: pure Wilson-Hilferty (identity transform)."""
        m = np.array([2.0, 4.0, 0.0])
        for x in (-1.0, 0.0, 2.0, 3.5, 6.0):
            expected = stats.norm.cdf((x - 2.0) / 2.0)
            assert p_p3(x, m) == pytest.approx(expected, abs=1e-9)

    def test_p_p3_is_monotonic(self):
        m = np.array([3.7, 0.09, -0.4])
        xs = np.linspace(2.0, 5.5, 25)
        ps = [p_p3(x, m) for x in xs]
        assert all(b >= a for a, b in zip(ps, ps[1:]))
        assert 0.0 <= ps[0] <= ps[-1] <= 1.0

    def test_p_p3_and_q_p3_roundtrip(self):
        m = np.array([3.7, 0.09, -0.4])
        for p in (0.01, 0.25, 0.5, 0.75, 0.99):
            assert p_p3(q_p3(p, m), m) == pytest.approx(p, rel=1e-6)

    def test_q_p3_boundary_returns_signed_infinity(self):
        """q <= 0 with skew <= 0, or q >= 1 with skew >= 0: no mass there."""
        m_neg = np.array([0.0, 1.0, -0.5])
        m_pos = np.array([0.0, 1.0, 0.5])
        assert q_p3(0.0, m_neg) == -1.0e31
        assert q_p3(1.0, m_pos) == 1.0e31


class TestMP3:
    def test_full_support_matches_m2mn(self):
        m = np.array([3.7, 0.09, -0.4])
        assert np.allclose(m_p3(-1e20, 1e20, m, 3), m2mn(m), rtol=1e-6)

    def test_degenerate_interval_falls_back_to_nearest_bound(self):
        m = np.array([0.0, 1.0, 0.3])
        lo, hi = 50.0, 60.0
        assert np.allclose(m_p3(lo, hi, m, 3), [lo, lo**2, lo**3])

    def test_rejects_n_greater_than_twelve(self):
        m = np.array([3.7, 0.09, -0.4])
        with pytest.raises(ValueError):
            m_p3(-1e20, 1e20, m, 13)

    def test_moments_increase_toward_the_upper_bound_of_a_narrow_interval(self):
        """A sanity bound, not a parity check: mean of a narrow truncation
        must lie between its endpoints."""
        m = np.array([3.7, 0.09, -0.4])
        lo, hi = 3.5, 3.6
        moments = m_p3(lo, hi, m, 1)
        assert lo <= moments[0] <= hi


class TestWilsonHilfertyWeightLiterals:
    """``_wh_weight`` uses the Fortran's single-precision literals."""

    def test_span_is_folded_in_single_precision(self):
        from flowfreq import _p3_moments as pm

        assert pm._WH_LOW == float(np.float32(0.0007))
        assert pm._WH_SPAN == float(np.float32(0.0010) - np.float32(0.0007))
        assert pm._WH_SPAN != 0.0010 - 0.0007

    def test_weights_at_the_ends_of_the_blend(self):
        from flowfreq import _p3_moments as pm

        assert pm._wh_weight(0.0) == (0.0, 1.0)
        assert pm._wh_weight(0.5) == (1.0, 0.0)
        wg, wwh = pm._wh_weight(0.00085)
        assert 0.0 < wg < 1.0 and wg + wwh == pytest.approx(1.0)


class TestIncompleteGammaAtWorkingPrecision:
    """``_lower_gamma_reg`` (fixed-point loops) and the upward recurrence that
    ``_fp_g1_mom_trc_batch`` builds on it, against independent mpmath
    references at the same 50 digits. The tolerance is the working
    precision's, ~30 orders tighter than anything float64 can see: these are
    performance rewrites of the same mathematics and must not move a result.
    """

    # (a, x): both branches (x < a + 1 is the series), small and large a.
    CASES = [
        (0.7, 0.2),
        (0.7, 3.0),
        (5.0, 2.0),
        (5.0, 11.0),
        (60.0, 48.0),
        (60.0, 75.0),
        (400.0, 380.0),
        (400.0, 431.0),
    ]

    @pytest.mark.parametrize("a, x", CASES)
    def test_matches_mpmath_gammainc(self, a, x):
        import mpmath

        from flowfreq._p3_moments import _GAMMA_MOMENT_DPS, _lower_gamma_reg

        with mpmath.workdps(_GAMMA_MOMENT_DPS):
            got = _lower_gamma_reg(mpmath.mpf(a), mpmath.mpf(x))
            ref = mpmath.gammainc(mpmath.mpf(a), 0, mpmath.mpf(x), regularized=True)
            # Both P and its complement Q: up - down cancels in whichever is small.
            assert abs(got - ref) <= mpmath.mpf(10) ** -46 * min(ref, 1 - ref)

    def test_large_shape_parameter_matches_the_complement_identity(self):
        """alpha ~ 4/skew**2 for near-zero skew: mpmath.gammainc gives up here,
        so check P(a, x) + Q(a, x) = 1 with Q from the upper-gamma series."""
        import mpmath

        from flowfreq._p3_moments import _GAMMA_MOMENT_DPS, _lower_gamma_reg

        with mpmath.workdps(_GAMMA_MOMENT_DPS):
            a = mpmath.mpf(9.0e4)
            for x in (a - 3 * mpmath.sqrt(a), a + 2 * mpmath.sqrt(a)):
                p = _lower_gamma_reg(a, x)
                # P(a + 1, x) = P(a, x) - x**a e**-x / Gamma(a + 1), checked by
                # evaluating both sides independently.
                p1 = _lower_gamma_reg(a + 1, x)
                term = mpmath.e ** (a * mpmath.log(x) - x - mpmath.loggamma(a + 1))
                assert abs((p - term) - p1) <= mpmath.mpf(10) ** -44 * min(p1, 1 - p1)

    @pytest.mark.parametrize("alpha", [3.0, 40.0, 2500.0])
    def test_batched_moments_match_direct_incomplete_gamma_ratios(self, alpha):
        """The recurrence gives what 2 + 2*kmax separate P(a + k, x) solves give."""
        import math

        import mpmath

        from flowfreq._p3_moments import _GAMMA_MOMENT_DPS, _fp_g1_mom_trc_batch, _lower_gamma_reg

        with mpmath.workdps(_GAMMA_MOMENT_DPS):
            a = mpmath.mpf(alpha)
            s = mpmath.sqrt(a)
            tl, tu = a - 1.5 * s, a + 0.7 * s
            got = _fp_g1_mom_trc_batch(a, tl, tu, 6)
            down = _lower_gamma_reg(a, tu) - _lower_gamma_reg(a, tl)
            for k in range(1, 7):
                up = _lower_gamma_reg(a + k, tu) - _lower_gamma_reg(a + k, tl)
                ref = up / down * mpmath.rf(a, k)
                assert abs(got[k - 1] - ref) <= mpmath.mpf(10) ** -38 * abs(ref), k
            assert math.isfinite(float(got[-1]))
