"""Native port of peakfq's Pearson III moment/parameter primitives.

Phase 1 of the ``var_mom`` port (TODO.md P3). ``var_mom`` and the routines it
needs are not implemented anywhere in flowfreq yet; this module is the leaf
layer of that dependency tree -- the piece with no further flowfreq
dependencies of its own:

    var_mom -> mP3(k=6), pP3, varb, varc, d_est -> m2p, qP3, expmomderiv
                m2mn

Everything here is a direct, function-for-function transcription of
``vendor/peakfqr/src/emafit.f`` and ``vendor/peakfqr/src/probfun.f``, checked
routine-by-routine against the Fortran through the oracles
``build_fortran/_emafort.pyf`` exposes (``m2p``, ``m2mn``, ``mn2m``, ``pp3``,
``qp3sub``, ``mp3``) -- see ``tests/fortran_parity/test_fortran_oracles.py``.

Not ported here: ``varb``, ``varc``, ``d_est``, ``var_mom`` itself, or
anything above it (``mse_ema``, ``VAR_EMAB``, ``mn2mvarb``, ``ci_ema_m3b``).
Those are later phases. Nothing in this module is wired into
``Bulletin17C``/``ExpectedMomentsAlgorithm`` yet -- it exists to be verified
in isolation first.

Every function takes and returns *central* moments ``(mean, variance, skew)``
in the same convention as ``EMAParameters``/``bulletin17c.py``, log10 space
when used on flood data, exactly as the Fortran does.
"""

from __future__ import annotations

import math
from typing import Tuple

import mpmath
import numpy as np
from mpmath.libmp import from_man_exp, to_fixed
from scipy import stats
from scipy.special import ndtri

__all__ = [
    "m2p",
    "m2mn",
    "mn2m",
    "p_p3",
    "q_p3",
    "m_p3",
]


def _f32(value: float) -> float:
    """A default-kind (single-precision) Fortran real literal, as a double.

    peakfq is compiled without ``-fdefault-real-8``, so an unsuffixed literal
    such as ``0.0007`` is a single-precision real, widened when it meets a
    double.
    """
    return float(np.float32(value))


_WH_LOW = _f32(0.0007)
# ``(0.0010-0.0007)`` is two single-precision literals, folded in single precision.
_WH_SPAN = float(np.float32(0.0010) - np.float32(0.0007))
_WH_PI = _f32(3.14159265359)


def _wh_weight(skew: float) -> Tuple[float, float]:
    """(wg, wwh): weight on the incomplete-gamma vs. Wilson-Hilferty solution.

    ``emafit.f`` blends the two smoothly over ``|skew|`` in [0.0007, 0.0010]
    (e.g. ``pP3`` at emafit.f:3231) so that pure-gamma code near skew = 0,
    where alpha = 4/skew**2 blows up, is never evaluated.
    """
    w = max(0.0, (abs(skew) - _WH_LOW) / _WH_SPAN)
    wwh = 0.0 if w >= 1.0 else (1.0 + math.cos(_WH_PI * w)) / 2.0
    return 1.0 - wwh, wwh


def _whlp2z(x: float, skew: float) -> float:
    """Wilson-Hilferty standardizing transform (``emafit.f`` ``whlp2z``)."""
    if skew == 0.0:
        return x
    return (6.0 / skew) * (skew**2 / 36.0 - 1.0 + max(0.0, 1.0 + skew * x / 2.0) ** (1.0 / 3.0))


def _whz2lp(z: float, skew: float) -> float:
    """Inverse of ``_whlp2z`` (``emafit.f`` ``whz2lp``)."""
    if skew == 0.0:
        return z
    return (2.0 / skew) * (1.0 + skew * z / 6.0 - skew**2 / 36.0) ** 3 - 2.0 / skew


def m2p(m: np.ndarray) -> np.ndarray:
    """Central moments (mu, sigma^2, skew) -> Pearson III params (tau, alpha, beta).

    Pure algebra, ``probfun.f`` ``M2P`` (line 768). Undefined at skew = 0
    (alpha = 4/skew**2), which is fine: callers only reach this branch where
    ``_wh_weight`` gives it nonzero weight, and that excludes skew = 0.
    """
    mu, var, skew = float(m[0]), float(m[1]), float(m[2])
    alpha = 4.0 / skew**2
    beta = math.sqrt(var / alpha)
    if skew < 0.0:
        beta = -beta
    tau = mu - alpha * beta
    return np.array([tau, alpha, beta])


def m2mn(m: np.ndarray) -> np.ndarray:
    """Central moments -> noncentral (raw) moments about zero.

    ``probfun.f`` ``M2MN`` (line 818).
    """
    m1, m2, m3 = float(m[0]), float(m[1]), float(m[2])
    mn1 = m1
    mn2 = m2 + m1**2
    mn3 = m3 * m2**1.5 + 3.0 * mn2 * m1 - 2.0 * m1**3
    return np.array([mn1, mn2, mn3])


def mn2m(mn: np.ndarray) -> np.ndarray:
    """Noncentral moments -> central moments; inverse of ``m2mn``.

    ``probfun.f`` ``MN2M`` (line 830).
    """
    mn1, mn2, mn3 = float(mn[0]), float(mn[1]), float(mn[2])
    m1 = mn1
    m2 = mn2 - mn1**2
    m3 = (mn3 - 3.0 * mn2 * mn1 + 2.0 * mn1**3) / m2**1.5
    return np.array([m1, m2, m3])


def _fp_g3_cdf(x: float, parms: np.ndarray) -> float:
    """CDF of the 3-parameter (shifted) gamma at x. ``probfun.f`` FP_G3_CDF."""
    tau, alpha, beta = parms
    y = x - tau
    if beta > 0.0:
        return float(stats.gamma.cdf(y / beta, alpha))
    if beta < 0.0:
        return 1.0 - float(stats.gamma.cdf(y / beta, alpha))
    return 1.0 if y > 0.0 else 0.0


def p_p3(x: float, m: np.ndarray) -> float:
    """CDF at x of a Pearson III variate with central moments m.

    ``emafit.f`` ``pP3`` (line 3210): blends the incomplete-gamma CDF with a
    Wilson-Hilferty normal approximation via ``_wh_weight``.
    """
    skew = float(m[2])
    wg, wwh = _wh_weight(skew)
    p_g = 0.0
    p_wh = 0.0
    if wg > 0.0:
        p_g = _fp_g3_cdf(x, m2p(m))
    if wwh > 0.0:
        mu, s = float(m[0]), math.sqrt(float(m[1]))
        z = _whlp2z((x - mu) / s, skew)
        p_wh = float(stats.norm.cdf(z))
    return wg * p_g + wwh * p_wh


def _fp_g1_icdf(p: float, alpha: float) -> float:
    """Inverse CDF of the 1-parameter gamma(alpha). ``probfun.f`` FP_G1_ICDF.

    ``Gamma(alpha, scale=1)`` and ``chi2(2*alpha)/2`` are the same
    distribution, which is what the Fortran's ``DCHIIN`` call exploits;
    ``scipy.stats.gamma.ppf`` gives the identical quantile directly.
    """
    if p > 0.999999999:
        return 9.0e99
    if p < 0.000000001:
        return 0.0
    return float(stats.gamma.ppf(p, alpha))


def _fp_g3_icdf(q: float, parms: np.ndarray) -> float:
    """Inverse CDF of the shifted gamma. ``probfun.f`` FP_G3_ICDF/FP_G2_ICDF."""
    tau, alpha, beta = parms
    if beta > 0.0:
        return tau + beta * _fp_g1_icdf(q, alpha)
    if beta < 0.0:
        return tau + beta * _fp_g1_icdf(1.0 - q, alpha)
    raise ValueError("invalid parameters in fp_g3_icdf: beta = 0")


def q_p3(q: float, m: np.ndarray) -> float:
    """Inverse CDF (quantile) of a Pearson III variate with central moments m.

    ``emafit.f`` ``qP3`` (line 3266). The two boundary branches only
    short-circuit (return +/-infinity) when the requested tail and the sign
    of skew agree that the distribution has no mass there; otherwise, like
    the Fortran, this falls through to the ordinary blended computation.
    """
    skew = float(m[2])
    if q <= 0.0 and skew <= 0.0:
        return -1.0e31
    if q >= 1.0 and skew >= 0.0:
        return 1.0e31

    wg, wwh = _wh_weight(skew)
    q_g = 0.0
    q_wh = 0.0
    if wg > 0.0:
        q_g = _fp_g3_icdf(q, m2p(m))
    if wwh > 0.0:
        mu, s = float(m[0]), math.sqrt(float(m[1]))
        q_wh = mu + s * _whz2lp(float(ndtri(q)), skew)
    return wg * q_g + wwh * q_wh


#: Decimal digits of precision for the mpmath incomplete-gamma moment ratio
#: in _fp_g1_mom_trc. Real flood records have at-site skew near zero, which
#: makes alpha = 4/skew**2 enormous (Big Sandy: ~9e4; a skew of 0.001 gives
#: ~4e6) -- emafit.f's own author hit this: FP_G1_CDF/FP_G1_MOM_TRC were
#: patched in September 2024 to evaluate DGAMDF in real128 (~34 decimal
#: digits) specifically because float64 catastrophically cancels here. The
#: up/down ratio subtracts two incomplete-gamma values that agree in most of
#: their leading digits, then multiplies by alpha ~10**5 raised to the k-th
#: power (k up to 6) -- float64's ~16 digits are consumed by the
#: subtraction long before that amplification, verified against the
#: mp3/qp3sub/pp3 Fortran oracles (tests/fortran_parity/test_fortran_oracles.py
#: ``TestMP3Port``) where float64 alone was off by double digits of percent
#: at k=5-6 on Big Sandy's own censoring thresholds.
_GAMMA_MOMENT_DPS = 50

#: Extra bits, beyond the working precision, of _lower_gamma_reg's fixed-point
#: scale. Absorbs the floor rounding of a few hundred to a few thousand
#: series/continued-fraction steps with room to spare.
_FIXED_GUARD_BITS = 64


def _lower_gamma_reg(a: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """Regularized lower incomplete gamma P(a, x), a, x >= 0, a arbitrary.

    ``mpmath.gammainc(..., regularized=True)`` evaluates this through a
    hypergeometric series that, for the shape parameters this module needs
    (alpha = 4/skew**2, in the hundreds of thousands to millions for
    near-zero skew), fails outright with "Hypergeometric series converges
    too slowly" -- the series needs on the order of x terms, and x is
    itself ~alpha here. This is the classical combined series/continued-
    fraction algorithm (Numerical Recipes ``gammp``/``gammq``), which
    converges in O(sqrt(a)) terms in both regimes regardless of how large a
    is, evaluated at ``mpmath``'s working precision.

    The two loops run in fixed-point integer arithmetic, as ``mpmath``'s own
    series code does, rather than on ``mpf`` objects: hundreds of terms per
    call, and this function is the innermost cost of every EMA iteration,
    ``mse_ema`` and ``var_emab`` -- where per-operation ``mpf`` allocation,
    dispatch and renormalization made it the single largest cost in the test
    suite. The fixed-point scale carries ``_FIXED_GUARD_BITS`` beyond the
    working precision plus twice the integer bits of ``max(a, x)``, so the
    result agrees with the ``mpf`` evaluation to ~1e-48 relative at 50
    digits -- the working precision's own rounding level, some thirty orders
    below anything that survives the conversion to float64 (see
    ``tests/test_p3_moments.py``).
    """
    if x <= 0:
        return mpmath.mpf(0)
    a = mpmath.mpf(a)
    x = mpmath.mpf(x)
    prec, rnd = mpmath.mp._prec_rounding
    eps = mpmath.mpf(10) ** (-(mpmath.mp.dps + 5))
    wp = prec + _FIXED_GUARD_BITS + 2 * max(int(x), int(a), 1).bit_length()
    one = 1 << wp
    eps_f = to_fixed(eps._mpf_, wp)
    a_f = to_fixed(a._mpf_, wp)
    x_f = to_fixed(x._mpf_, wp)
    maxiter = 2_000_000
    if x < a + 1:
        s_f = (one << wp) // a_f
        d_f = s_f
        ap_f = a_f
        for _ in range(maxiter):
            ap_f += one
            d_f = (d_f * x_f) // ap_f
            s_f += d_f
            if d_f * one < s_f * eps_f:  # d < s * eps
                break
        else:
            raise RuntimeError("_lower_gamma_reg: series did not converge")
        s = mpmath.mp.make_mpf(from_man_exp(s_f, -wp, prec, rnd))
        return s * mpmath.e ** (-x + a * mpmath.log(x) - mpmath.loggamma(a))
    # Modified Lentz continued fraction for Q(a, x).
    fpmin_f = to_fixed((mpmath.mpf(10) ** (-(mpmath.mp.dps + 20)))._mpf_, wp)
    b_f = x_f + one - a_f
    c_f = (one << wp) // fpmin_f
    d_f = (one << wp) // b_f
    h_f = d_f
    for i in range(1, maxiter + 1):
        an_f = -i * (i * one - a_f)
        b_f += 2 * one
        d_f = ((an_f * d_f) >> wp) + b_f
        if abs(d_f) < fpmin_f:
            d_f = fpmin_f
        c_f = b_f + (an_f << wp) // c_f
        if abs(c_f) < fpmin_f:
            c_f = fpmin_f
        d_f = (one << wp) // d_f
        delta_f = (d_f * c_f) >> wp
        h_f = (h_f * delta_f) >> wp
        if abs(delta_f - one) < eps_f:
            break
    else:
        raise RuntimeError("_lower_gamma_reg: continued fraction did not converge")
    h = mpmath.mp.make_mpf(from_man_exp(h_f, -wp, prec, rnd))
    gammcf = mpmath.e ** (-x + a * mpmath.log(x) - mpmath.loggamma(a)) * h
    return 1 - gammcf


def _gamma_term(a: mpmath.mpf, x: mpmath.mpf) -> mpmath.mpf:
    """``x**a * exp(-x) / Gamma(a + 1)``, the step in P(a+1, x) = P(a, x) - term."""
    if x <= 0:
        return mpmath.mpf(0)
    return mpmath.e ** (a * mpmath.log(x) - x - mpmath.loggamma(a + 1))


def _fp_g1_mom_trc_batch(alpha, tl, tu, kmax: int) -> list:
    """[E[X^k | tl < X < tu] for k = 1..kmax], X ~ Gamma(alpha, scale=1).

    Batched form of ``FP_G1_MOM_TRC``: the incomplete-gamma ``down`` term
    does not depend on k, so computing k = 1..kmax together shares it
    instead of recomputing it kmax times -- the dominant cost once this is
    called from a numerical Jacobian (``flowfreq._var_mom._dexpect``),
    which needs several k at once, several times over for a finite
    difference.
    """
    tl1 = max(min(mpmath.mpf(0), tu), tl)
    if tl1 == tu:
        return [tl1**k for k in range(1, kmax + 1)]
    p_tu = _lower_gamma_reg(alpha, tu)
    p_tl = _lower_gamma_reg(alpha, tl1)
    down = p_tu - p_tl
    # P(alpha + k, x) for k = 1..kmax from P(alpha, x) by the exact upward
    # recurrence P(a + 1, x) = P(a, x) - x**a * exp(-x) / Gamma(a + 1)
    # (DLMF 8.8.5), instead of 2*kmax further _lower_gamma_reg solves. Each
    # of those is an O(sqrt(alpha))-term series or continued fraction in
    # mpmath, and alpha = 4/skew**2 runs to 1e5-1e6 on real records: this
    # was the single largest cost in the test suite (the EMA fixed point,
    # mse_ema and var_emab all bottom out here). The recurrence terms are
    # evaluated at the same working precision as P itself, so the k-th
    # truncated moment keeps the ~50-digit accuracy the direct solve had --
    # verified unchanged in float64 against the mp3 Fortran oracle
    # (tests/fortran_parity/test_fortran_oracles.py::TestMP3Port).
    term_tu = _gamma_term(alpha, tu)
    term_tl = _gamma_term(alpha, tl1)
    result = []
    for k in range(1, kmax + 1):
        p_tu -= term_tu
        p_tl -= term_tl
        term_tu = term_tu * tu / (alpha + k)
        term_tl = term_tl * tl1 / (alpha + k)
        if down > 0:
            up = p_tu - p_tl
            ans = up / down
            for j in range(k):
                ans *= alpha + j
            result.append(ans)
        else:
            result.append(tl1**k if tl1 > alpha else tu**k)
    return result


def _fp_g2_mom_trc_batch(alpha, beta, tl, tu, kmax: int) -> list:
    """Batched ``FP_G2_MOM_TRC``: [E[X^k | tl < X < tu] for k = 1..kmax]."""
    if beta == 0.0:
        raise ValueError("invalid parameters in fp_g2_mom_trc: beta = 0")
    alpha_mp, beta_mp = mpmath.mpf(alpha), mpmath.mpf(beta)
    tl_mp, tu_mp = mpmath.mpf(tl), mpmath.mpf(tu)
    if beta > 0.0:
        fp = _fp_g1_mom_trc_batch(alpha_mp, tl_mp / beta_mp, tu_mp / beta_mp, kmax)
    else:
        fp = _fp_g1_mom_trc_batch(alpha_mp, tu_mp / beta_mp, tl_mp / beta_mp, kmax)
    return [beta_mp**k * fp[k - 1] for k in range(1, kmax + 1)]


def _gamma_trunc_moments(tau, alpha, beta, tl, tu, kmax: int) -> list:
    """[E[X^k | tl < X < tu] for k = 1..kmax] for the 3-parameter shifted gamma.

    Batched counterpart to ``m_p3``'s gamma branch (same binomial
    expansion, ``emafit.f:3049``), expressed directly in ``(tau, alpha,
    beta)`` rather than central moments -- what ``flowfreq._var_mom``'s
    ``_dexpect`` differentiates numerically with respect to each
    parameter. Accepts float or ``mpmath.mpf`` inputs; always returns a
    list of ``mpmath.mpf``.
    """
    tau_mp = mpmath.mpf(tau)
    fp = [mpmath.mpf(1)] + _fp_g2_mom_trc_batch(alpha, beta, tl - tau, tu - tau, kmax)
    result = []
    for k in range(1, kmax + 1):
        mg = fp[k]
        for j in range(k):
            mg += math.comb(k, j) * tau_mp ** (k - j) * fp[j]
        result.append(mg)
    return result


def m_p3(tl: float, tu: float, m: np.ndarray, n: int) -> np.ndarray:
    """E[X^k | tl < X < tu] for k = 1..n, X ~ Pearson III with central moments m.

    ``emafit.f`` ``mP3`` (line 2983). This is the truncated-moment machinery
    ``var_mom`` needs for ``mu_x``/``e_x`` on censored intervals, and the
    piece TODO.md P3 names as the source of ``_ema_iteration``'s residual
    disagreement with ``moms_p3`` on censored rows (Cains Coulee: 0.70%
    variance, 4.94% skew) -- flowfreq's own
    ``_truncated_gamma_moment``/``_truncated_normal_moments`` are an
    approximation of what this function computes exactly.

    Like the Fortran, blends an incomplete-gamma solution (large |skew|)
    with a Wilson-Hilferty power-series solution (small |skew|).
    """
    if n > 12:
        raise ValueError("m_p3: n too large (greater than 12)")

    mu, var, skew = float(m[0]), float(m[1]), float(m[2])
    s = math.sqrt(var)
    zl = max(-1e20, min(1e20, _whlp2z((tl - mu) / s, skew)))
    zu = max(-1e20, min(1e20, _whlp2z((tu - mu) / s, skew)))
    cdfu = float(stats.norm.cdf(zu))
    cdfl = float(stats.norm.cdf(zl))

    if cdfu == cdfl:
        # Off the support of f -- use whichever bound is closer to the mean.
        t = tl if abs(tl - mu) < abs(tu - mu) else tu
        return np.array([t**k for k in range(1, n + 1)])

    wg, wwh = _wh_weight(skew)

    mg = np.zeros(n + 1)
    if wg > 0.0:
        tau, alpha, beta = m2p(m)
        with mpmath.workdps(_GAMMA_MOMENT_DPS):
            moments = _gamma_trunc_moments(tau, alpha, beta, tl, tu, n)
        mg[1:] = [float(v) for v in moments]

    mwh = np.zeros(n + 1)
    if wwh > 0.0:
        pdfu = float(stats.norm.pdf(zu))
        pdfl = float(stats.norm.pdf(zl))
        expz = np.zeros(3 * n + 1)
        expz[0] = 1.0
        expz[1] = (-pdfu + pdfl) / (cdfu - cdfl)
        zu1 = zu if cdfu < 1.0 else 0.0
        zl1 = zl if cdfl > 0.0 else 0.0
        for i in range(2, 3 * n + 1):
            expz[i] = (-(zu1 ** (i - 1)) * pdfu + (zl1 ** (i - 1)) * pdfl) / (cdfu - cdfl) + (
                i - 1
            ) * expz[i - 2]

        a = [
            -skew * (3888.0 - 108.0 * skew**2 + skew**4) / 23328.0,
            (1.0 - skew**2 / 36.0) ** 2,
            (skew / 6.0) * (1.0 - skew**2 / 36.0),
            skew**2 / 108.0,
        ]

        b = np.zeros((n + 1, 3 * n + 1))
        b[0, 0] = 1.0
        fwh = np.zeros(n + 1)
        fwh[0] = 1.0
        for i in range(1, n + 1):
            for j in range(0, 3 * i + 1):
                for k in range(max(0, j - 3 * (i - 1)), min(j, 3) + 1):
                    b[i, j] += b[i - 1, j - k] * a[k]
                fwh[i] += b[i, j] * expz[j]

        for i in range(1, n + 1):
            mwh[i] = s**i * fwh[i]
            for j in range(i):
                mwh[i] += math.comb(i, j) * mu ** (i - j) * s**j * fwh[j]

    return np.array([wg * mg[i] + wwh * mwh[i] for i in range(1, n + 1)])
