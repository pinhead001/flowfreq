"""Tests for hydrolib._p3_moments that do not require the Fortran extension.

Parity against the vendored Fortran itself lives in
tests/fortran_parity/test_fortran_oracles.py (skipped when the extension
isn't built). These check the pure-Python properties that should hold
regardless: algebraic identities, known closed forms, and error handling.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import stats

from hydrolib._p3_moments import m2mn, m2p, m_p3, mn2m, p_p3, q_p3


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
