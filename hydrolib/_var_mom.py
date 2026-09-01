"""Native port of peakfq's var_mom: variance-covariance of the EMA moment estimators.

Phase 2 of the ``var_mom`` port (TODO.md P3). Phase 1 (``hydrolib._p3_moments``)
ported the leaf layer -- ``m2p``, ``m2mn``, ``mn2m``, ``p_p3``, ``q_p3``,
``m_p3``. This module composes the rest of the dependency tree on top of it:

    var_mom -> mP3(k=6), pP3, varb, varc, d_est -> m2p, qP3, expmomderiv
                m2mn

Everything here is a direct, function-for-function transcription of
``vendor/peakfqr/src/emafit.f`` and ``vendor/peakfqr/src/probfun.f``, checked
routine-by-routine against six Fortran oracles ``build_fortran/_emafort.pyf``
exposes (``varb``, ``varc``, ``d_est``, ``expmomderiv``, and Phase 1's
``m2p``/``mp3``/``pp3``/``qp3sub``) -- see
``tests/fortran_parity/test_fortran_oracles.py``.

Not ported here: ``mn2mvarb`` and anything above it (``mse_ema``, ``VAR_EMAB``,
``ci_ema_m3b``, the pseudo effective record length). Those are later phases.
Nothing in this module is wired into ``Bulletin17C``/``ExpectedMomentsAlgorithm``
yet -- Phase 3 decides how (or whether) it changes what those compute.

The 3x3 matrix bookkeeping the Fortran does through named helper subroutines
(``DSET``, ``DMSUM``, ``DMRRRR``, ``DMXYTF``, ``DMMULT``, ``DLGINV``) is just
generic linear algebra -- filling, summing, and multiplying small matrices,
and a matrix inverse -- with no numerical subtlety of its own, so it is done
here with plain ``numpy``/``numpy.linalg`` rather than transcribed routine by
routine.
"""

from __future__ import annotations

import math

import mpmath
import numpy as np

from hydrolib._p3_moments import _GAMMA_MOMENT_DPS, _gamma_trunc_moments, m2p, m_p3, p_p3, q_p3

__all__ = ["expmomderiv", "d_est", "varb", "varc", "var_mom"]

#: emafit.f:2485 -- below this tail probability, d_est skips that tail
#: entirely (the matching qP3 call is skipped too, not just zeroed).
_D_EST_PMIN = 1.0e-6

#: emafit.f:2485/3283 -- var_mom's and qP3's stand-ins for -infinity/+infinity.
_D_EST_INF = 1.0e19

#: emafit.f:2322/2324 -- var_mom's own sentinel for "unbounded" and the
#: [skewmin, skewmax] clamp it applies to the skew it works with internally
#: (see var_mom below). Distinct from _D_EST_INF: this is what var_mom's
#: *own* mP3/pP3 calls use for the open tails, not what d_est/qP3 use.
_VAR_MOM_XINF = 999.0
_SKEWMIN = 0.06324555
_SKEWMAX = 1.41


def _dpdm(parms: np.ndarray) -> np.ndarray:
    """Jacobian of (tau, alpha, beta) w.r.t. noncentral moments (m1, m2, m3).

    ``JACT[i, j] = d param_i / d moment_j``. ``probfun.f`` ``DPDM`` (line
    648): closed-form algebra, no truncation or quadrature involved, so
    plain float64 throughout.
    """
    tau, alpha, beta = float(parms[0]), float(parms[1]), float(parms[2])

    # P2MN (probfun.f:787): params -> noncentral moments of the *untruncated*
    # distribution those params describe.
    m1 = alpha * beta + tau
    m2 = alpha * (1.0 + alpha) * beta**2 + 2.0 * alpha * beta * tau + tau**2
    m3 = (
        alpha * (1.0 + alpha) * (2.0 + alpha) * beta**3
        + 3.0 * alpha * (1.0 + alpha) * beta**2 * tau
        + 3.0 * alpha * beta * tau**2
        + tau**3
    )

    v = m2 - m1**2
    dvdm1, dvdm2, dvdm3 = -2.0 * m1, 1.0, 0.0

    s = m3 - 3.0 * m2 * m1 + 2.0 * m1**3
    dsdm1, dsdm2, dsdm3 = -3.0 * m2 + 6.0 * m1**2, -3.0 * m1, 1.0

    a = 4.0 * v**3 / s**2
    dadm1 = 4.0 * (3.0 * v**2 * dvdm1 * s**2 - 2.0 * s * dsdm1 * v**3) / s**4
    dadm2 = 4.0 * (3.0 * v**2 * dvdm2 * s**2 - 2.0 * s * dsdm2 * v**3) / s**4
    dadm3 = 4.0 * (3.0 * v**2 * dvdm3 * s**2 - 2.0 * s * dsdm3 * v**3) / s**4

    d = v / a
    dddm1 = (dvdm1 * a - dadm1 * v) / a**2
    dddm2 = (dvdm2 * a - dadm2 * v) / a**2
    dddm3 = (dvdm3 * a - dadm3 * v) / a**2

    b = math.copysign(math.sqrt(d), s)
    dbdm1 = 0.5 * dddm1 / b
    dbdm2 = 0.5 * dddm2 / b
    dbdm3 = 0.5 * dddm3 / b

    dtdm1 = 1.0 - (dadm1 * b + a * dbdm1)
    dtdm2 = -(dadm2 * b + a * dbdm2)
    dtdm3 = -(dadm3 * b + a * dbdm3)

    # JAC[moment, param]; DPDM returns its transpose, JACT[param, moment].
    jac = np.array(
        [
            [dtdm1, dadm1, dbdm1],
            [dtdm2, dadm2, dbdm2],
            [dtdm3, dadm3, dbdm3],
        ]
    )
    return jac.T


#: Relative step for _dexpect's manual central difference. mpmath's working
#: precision (_GAMMA_MOMENT_DPS, 50 decimal digits) leaves so much headroom
#: over this that truncation error (~h**2) dominates cancellation error
#: (~eps/h) by many orders of magnitude at any reasonable h -- unlike in
#: float64, there is no delicate step-size tradeoff to make here.
_DEXPECT_STEP = mpmath.mpf("1e-12")


def _dexpect(tau: float, alpha: float, beta: float, tl: float, tu: float):
    """E[X^j] over (tl, tu), j=1..3, and its Jacobian w.r.t. (tau, alpha, beta).

    ``probfun.f`` ``DEXPECT`` (line 1136).

    The Fortran computes both the moments and the Jacobian from a shared
    closed-form chain (``ADJ``/``G1``/``B`` and their derivatives, plus
    ``DDGAM`` for d(incomplete gamma)/d(alpha)). This differentiates
    ``_gamma_trunc_moments`` numerically instead, at the same working
    precision Phase 1's ``m_p3`` uses for the same reason
    (``_GAMMA_MOMENT_DPS``): that function is already the one this port
    trusts (verified against the Fortran and independently against
    quadrature in ``tests/fortran_parity/test_fortran_oracles.py``), so
    reusing it here means one moment implementation to get right rather
    than two. A manual central difference, batched over j=1..3 per
    evaluation rather than ``mpmath.diff`` called once per (j, param) pair,
    matters in practice: this is on the hot path of ``mse_ema``, called
    from ``ExpectedMomentsAlgorithm`` once per analysis, and the naive
    9-``mpmath.diff`` version measured close to a second per call on Big
    Sandy's censored group -- batching the three moments together and
    sharing the incomplete-gamma work within each evaluation (see
    ``_gamma_trunc_moments``) cuts that by roughly an order of magnitude.

    Returns zeros, matching the Fortran's early exit, when the interval is
    empty (``tu <= tl``) or ``beta == 0``.
    """
    if tu <= tl or beta == 0.0:
        return np.zeros(3), np.zeros((3, 3))

    with mpmath.workdps(_GAMMA_MOMENT_DPS):
        tau_mp, alpha_mp, beta_mp = mpmath.mpf(tau), mpmath.mpf(alpha), mpmath.mpf(beta)
        tl_mp, tu_mp = mpmath.mpf(tl), mpmath.mpf(tu)

        m3 = np.array(
            [float(v) for v in _gamma_trunc_moments(tau_mp, alpha_mp, beta_mp, tl_mp, tu_mp, 3)]
        )

        dm3 = np.zeros((3, 3))
        for col, center in enumerate((tau_mp, alpha_mp, beta_mp)):
            h = max(abs(center), mpmath.mpf(1)) * _DEXPECT_STEP
            args_plus = [tau_mp, alpha_mp, beta_mp]
            args_minus = [tau_mp, alpha_mp, beta_mp]
            args_plus[col] += h
            args_minus[col] -= h
            plus = _gamma_trunc_moments(*args_plus, tl_mp, tu_mp, 3)
            minus = _gamma_trunc_moments(*args_minus, tl_mp, tu_mp, 3)
            for row in range(3):
                dm3[row, col] = float((plus[row] - minus[row]) / (2 * h))
    return m3, dm3


def expmomderiv(parms: np.ndarray, xmin: float, xmax: float) -> np.ndarray:
    """Jacobian of E[X^j | xmin < X < xmax], j=1..3, w.r.t. noncentral moments.

    ``emafit.f`` ``expmomderiv`` (line 841): ``_dexpect``'s Jacobian w.r.t.
    (tau, alpha, beta) chained through ``_dpdm``'s Jacobian of
    (tau, alpha, beta) w.r.t. those moments.
    """
    tau, alpha, beta = float(parms[0]), float(parms[1]), float(parms[2])
    _, jac1 = _dexpect(tau, alpha, beta, xmin, xmax)
    jac2 = _dpdm(parms)
    return jac1 @ jac2


def d_est(nh: float, mc: np.ndarray, pa: float, pb: float) -> np.ndarray:
    """The open-tail bias-correction matrix var_mom sums over threshold groups.

    ``emafit.f`` ``d_est`` (line 2473). Skips (zeros) a tail whose
    probability is below ``_D_EST_PMIN`` -- ``q_p3`` is not evaluated
    there either, matching the Fortran, which never calls ``qP3`` at a
    probability so close to 0 or 1 that its own boundary clamps would bite.
    """
    parms = m2p(mc)
    p1 = pa
    p3 = 1.0 - pb

    jacl = np.zeros((3, 3))
    if p1 > _D_EST_PMIN:
        t1 = q_p3(pa, mc)
        jacl = (nh * p1) * expmomderiv(parms, -_D_EST_INF, t1)

    jacg = np.zeros((3, 3))
    if p3 > _D_EST_PMIN:
        t3 = q_p3(pb, mc)
        jacg = (nh * p3) * expmomderiv(parms, t3, _D_EST_INF)

    return jacl + jacg


def varb(mu_x: np.ndarray, nh: float, p1: float, p2: float, p3: float) -> np.ndarray:
    """Multinomial-sampling contribution to the moment estimators' covariance.

    ``emafit.f`` ``varb`` (line 2416): the delta-method "sandwich"
    ``mu_x @ vn @ mu_x.T``, where ``vn`` is the covariance of a
    multinomial(nh; p1, p2, p3) split across the three regions
    (below tl, between tl and tu, above tu) and ``mu_x``'s columns are
    each region's raw moments 1..3.
    """
    vn = np.array(
        [
            [nh * p1 * (1.0 - p1), -nh * p1 * p2, -nh * p1 * p3],
            [-nh * p1 * p2, nh * p2 * (1.0 - p2), -nh * p2 * p3],
            [-nh * p1 * p3, -nh * p2 * p3, nh * p3 * (1.0 - p3)],
        ]
    )
    return mu_x @ vn @ mu_x.T


def varc(e_x: np.ndarray, nh: float, p2: float) -> np.ndarray:
    """Within-interval variance from not knowing exactly where inside (tl, tu)

    a censored observation falls. ``emafit.f`` ``varc`` (line 2444).
    ``e_x`` is ``m_p3(tl, tu, mc, 6)`` -- moments 1..6 of the exact
    interval, needed because this reaches ``e_x[i + j]`` up to index 6.
    """
    mu_n = p2 * nh
    vc = np.zeros((3, 3))
    for i in range(1, 4):
        for j in range(1, 4):
            if i > j:
                vc[i - 1, j - 1] = vc[j - 1, i - 1]
            else:
                vc[i - 1, j - 1] = mu_n * (e_x[i + j - 1] - e_x[i - 1] * e_x[j - 1])
    return vc


def var_mom(
    n_in: np.ndarray, tl_in: np.ndarray, tu_in: np.ndarray, mc_in: np.ndarray
) -> np.ndarray:
    """Variance-covariance matrix of the fitted (mean, variance, skew) estimators.

    ``emafit.f`` ``var_mom`` (line 2299). ``n_in``/``tl_in``/``tu_in``
    describe distinct perception-threshold *groups* -- ``n_in[i]``
    observations share the pair ``(tl_in[i], tu_in[i])`` -- the same
    convention ``mseg_all_sub``/``detratsub`` use (see
    ``tests/fortran_parity/test_fortran_oracles.py::_threshold_groups``).
    ``mc_in`` is the fitted at-site central moments (mean, variance, skew).

    Two things var_mom does to its own working moments before calling into
    the leaf layer, easy to miss reading the Fortran linearly:

    * Recentres to mean 0 (``tl``/``tu`` are shifted by ``mc_in[0]`` before
      any ``p_p3``/``m_p3`` call).
    * Clamps ``|skew|`` to ``[_SKEWMIN, _SKEWMAX]`` (0.0632 to 1.41) before
      using it internally -- so even a near-zero at-site skew (Big Sandy:
      0.0066) never drives alpha = 4/skew**2 past roughly 1000 here, well
      short of the regime where Phase 1 found the Fortran's own ``mP3``
      losing precision (see TODO.md P3): that finding used ``mP3``'s raw,
      unclamped skew directly, which is not what ``var_mom`` itself ever
      hands it.
    """
    n_in = np.asarray(n_in, dtype=float)
    tl_in = np.asarray(tl_in, dtype=float)
    tu_in = np.asarray(tu_in, dtype=float)
    mc_in = np.asarray(mc_in, dtype=float)

    skew = math.copysign(min(max(abs(mc_in[2]), _SKEWMIN), _SKEWMAX), mc_in[2])
    mc = np.array([0.0, mc_in[1], skew])

    vb_t = np.zeros((3, 3))
    vc_t = np.zeros((3, 3))
    d_t = np.zeros((3, 3))
    n_t = 0.0

    for nh, tl_group, tu_group in zip(n_in, tl_in, tu_in):
        n_t += nh
        tl = tl_group - mc_in[0]
        tu = tu_group - mc_in[0]
        if tl > tu:
            raise ValueError("var_mom: invalid perception threshold tl > tu")

        pa = p_p3(tl, mc)
        pb = p_p3(tu, mc)
        p1, p2, p3 = pa, pb - pa, 1.0 - pb

        mnouta = m_p3(-_VAR_MOM_XINF, tl, mc, 3)
        mnoutb = m_p3(tu, _VAR_MOM_XINF, mc, 3)
        e_x = m_p3(tl, tu, mc, 6)

        mu_x = np.array([mnouta, e_x[:3], mnoutb]).T

        vb_t += varb(mu_x, nh, p1, p2, p3)
        vc_t += varc(e_x, nh, p2)
        d_t += d_est(nh, mc, pa, pb)

    ainv = np.eye(3) - d_t / n_t
    a = np.linalg.inv(ainv)
    varm = a @ (vb_t + vc_t) @ a.T
    return varm / n_t**2
