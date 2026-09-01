"""Native port of peakfq's mse_ema: MSE of an at-site EMA moment under ADJE.

Phase 3 of the ``var_mom`` port (TODO.md P3). Phases 1 (``hydrolib._p3_moments``)
and 2 (``hydrolib._var_mom``) built ``var_mom``, the variance-covariance of the
*noncentral* moment estimators. ``mse_ema`` is the last piece the ADJE skew
weighting needs -- ``as_G_mse = bias_adj * mseg(...)`` where
``bias_adj = mse_ema(censored)/mse_ema(uncensored)`` (TODO.md P3) -- and it
needs the covariance of the *central* moments, which is a nonlinear function
of the noncentral-moment covariance:

    mseg_all(ADJE) -> mse_ema -> var_mom -> mP3(k=6), pP3, varb, varc,
                                             d_est -> m2p, qP3, expmomderiv
                                  m2mn
                                  mn2mvarb -> mn2m_var, mc2mnvb,
                                              chol33   (iterative solver)

Everything here is a direct, function-for-function transcription of
``vendor/peakfqr/src/emafit.f`` and ``vendor/peakfqr/src/probfun.f``, checked
routine-by-routine against seven Fortran oracles ``build_fortran/_emafort.pyf``
exposes (``normquad``, ``gammaquad``, ``chol33``, ``mc2mnvb``, ``mn2m_var``,
``mn2mvarb``, ``mse_ema``) -- see ``tests/fortran_parity/test_fortran_oracles.py``.

One deliberate departure from the Fortran, not a shortcut around anything
uncertain: ``mn2mvarb`` solves an inverse problem (find the central-moment
covariance whose quadrature-based forward map, ``mc2mnvb``, reproduces the
given noncentral-moment covariance) that the Fortran solves with a bespoke
step-halved Newton iteration, checking positive-semi-definiteness via
``chol33`` at every trial step. Reparametrizing the unknown as the six free
entries of ``chol33``'s own factor ``V`` (so the candidate covariance is
always ``V.T @ V``, positive-semi-definite for *any* real ``V`` with no
guard needed) turns that into an unconstrained root-finding problem, solved
with ``scipy.optimize.root`` from the same starting point (``mn2m_var``'s
linearized estimate) the Fortran uses. Both sides are finding a root of the
same smooth map from the same starting point, so they should -- and are
verified to -- land on the same answer; see ``TestMn2mvarbPort``.

Not ported here: ``VAR_EMAB``, ``regmoms``, ``gridmake``'s use in
``ci_ema_m3b`` (the confidence-interval shape). Those are a later phase.
Nothing in this module is wired into ``ExpectedMomentsAlgorithm``/
``Bulletin17C`` yet.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from scipy.optimize import root
from scipy.special import gamma as _gamma_fn
from scipy.special import roots_genlaguerre

from hydrolib._p3_moments import m2mn, mn2m
from hydrolib._var_mom import var_mom

__all__ = ["mc2mnvb", "mn2m_var", "mn2mvarb", "mse_ema"]

#: emafit.f:1668/2324 -- var_mom's own clamp floor on |skew|, reused here as
#: mse_ema's *second*, smoother blend across skew = 0 (emafit.f:1687-1705).
_SKEWMIN = 0.06324555

#: gammaquad's own stabilization threshold (probfun.f:2694): above this
#: shape parameter, Gauss-Laguerre quadrature nodes become numerically
#: unreliable, so a mean-preserving Gamma(160, .) proxy is used instead.
_GAMMAQUAD_ALPHAMAX = 160.0

#: mn2mvarb always calls mc2mnvb with this grid (emafit.f:2558's NND/2,2,2/).
_NND = (2, 2, 2)


def _normquad(n: int) -> Tuple[np.ndarray, np.ndarray]:
    """Gauss-Hermite quadrature for the standard normal. ``probfun.f`` NORMQUAD.

    ``integral phi(z) g(z) dz ~= sum(w(i) * g(t(i)))``.
    """
    x, w = np.polynomial.hermite.hermgauss(n)
    return x * np.sqrt(2.0), w / np.sqrt(np.pi)


def _gammaquad(n: int, alpha: float, beta: float) -> Tuple[np.ndarray, np.ndarray]:
    """Gauss quadrature for a Gamma(alpha, beta) distribution. ``probfun.f`` GAMMAQUAD."""
    if alpha <= _GAMMAQUAD_ALPHAMAX:
        a, b, c = alpha, beta, 0.0
    else:
        a = _GAMMAQUAD_ALPHAMAX
        b = beta * np.sqrt(alpha / a)
        c = alpha * beta - a * b  # mean-preserving: a*b + c == alpha*beta
    t, w = roots_genlaguerre(n, a - 1.0)
    return c + t * b, w / _gamma_fn(a)


def _chol33(s: np.ndarray) -> Optional[np.ndarray]:
    """Permuted (2,1,3) Cholesky-like factor, V.T @ V == S. ``probfun.f`` CHOL33.

    Returns ``None`` (the Fortran's ``iflag = 1``) when S is not positive
    semi-definite in the order this pivoting checks -- s(1,1), s(2,2),
    s(3,3), then the two derived terms.
    """
    if s[0, 0] <= 0.0 or s[1, 1] <= 0.0 or s[2, 2] <= 0.0:
        return None
    v = np.zeros((3, 3))
    v[1, 1] = np.sqrt(s[1, 1])
    v[1, 0] = s[1, 0] / v[1, 1]
    v[1, 2] = s[1, 2] / v[1, 1]
    t1 = s[0, 0] - v[1, 0] ** 2
    if t1 < 0.0:
        return None
    v[0, 0] = np.sqrt(t1)
    v[0, 2] = (s[2, 0] - v[1, 2] * v[1, 0]) / v[0, 0]
    t1 = s[2, 2] - v[1, 2] ** 2 - v[0, 2] ** 2
    if t1 < 0.0:
        return None
    v[2, 2] = np.sqrt(t1)
    return v


def _tri_in(s: np.ndarray) -> np.ndarray:
    """The 6 distinct entries of a symmetric 3x3 matrix. ``emafit.f`` TRI_IN."""
    return np.array([s[0, 0], s[1, 0], s[2, 0], s[1, 1], s[1, 2], s[2, 2]])


def _tri_out(t: np.ndarray) -> np.ndarray:
    """Inverse of ``_tri_in``. ``emafit.f`` TRI_OUT."""
    s = np.empty((3, 3))
    s[0, 0], s[1, 0], s[2, 0], s[1, 1], s[1, 2], s[2, 2] = t
    s[0, 1], s[0, 2], s[2, 1] = s[1, 0], s[2, 0], s[1, 2]
    return s


def _v_from_tri(x: np.ndarray) -> np.ndarray:
    """The chol33-shaped V (see module docstring) from its 6 free entries."""
    v11, v21, v13, v22, v23, v33 = x
    return np.array([[v11, 0.0, v13], [v21, v22, v23], [0.0, 0.0, v33]])


def _gridmake(
    mc: np.ndarray, s_mc: np.ndarray, nnd: Tuple[int, int, int] = _NND
) -> Tuple[np.ndarray, np.ndarray]:
    """Quadrature (weights, points) for E[f(mc-space X)] given Cov(X) = s_mc.

    ``emafit.f`` ``GRIDMAKE``. Each axis is quadrature for an independently
    *shaped* variable (mean: normal, variance: gamma matched by method of
    moments to its own marginal mean/variance, skew: normal), then the
    ``chol33`` factor injects the desired covariance ``s_mc`` across axes,
    and ``mc`` is added back as the location.
    """
    v = _chol33(s_mc)
    if v is None:
        raise ValueError("gridmake: covariance matrix is not positive semi-definite")

    t1, w1 = _normquad(nnd[0])
    alpha = mc[1] ** 2 / s_mc[1, 1]
    beta = 1.0 / np.sqrt(alpha)
    t2, w2 = _gammaquad(nnd[1], alpha, beta)
    t2 = t2 - alpha * beta
    t3, w3 = _normquad(nnd[2])

    z = np.array(np.meshgrid(t1, t2, t3, indexing="ij")).reshape(3, -1)
    w = np.multiply.outer(np.multiply.outer(w1, w2), w3).reshape(-1)
    points = (mc[:, None] + v.T @ z).T
    return w, points


def _covw(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    """Weighted covariance. ``probfun.f`` COVW."""
    mx = np.average(x, weights=w)
    my = np.average(y, weights=w)
    return float(np.sum(w * (x - mx) * (y - my)) / np.sum(w))


def mc2mnvb(mc: np.ndarray, s_mc: np.ndarray, nnd: Tuple[int, int, int] = _NND) -> np.ndarray:
    """Covariance of noncentral moments given central moments (mc, s_mc).

    ``emafit.f`` ``MC2MNVB`` (line 2594): the forward map ``mn2mvarb``
    inverts. Quadrature points from ``_gridmake``, mapped through ``m2mn``,
    combined with ``_covw``.
    """
    w, points = _gridmake(mc, s_mc, nnd)
    mn_points = np.array([m2mn(p) for p in points])
    s_mn = np.empty((3, 3))
    for i in range(3):
        for j in range(3):
            s_mn[i, j] = _covw(mn_points[:, i], mn_points[:, j], w)
    return s_mn


def mn2m_var(mn: np.ndarray, s_mn: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Linearized (delta-method) (mc, s_mc) given (mn, s_mn).

    ``emafit.f`` ``mn2m_var`` (line 2864): closed-form Jacobian of central
    w.r.t. noncentral moments, then the usual sandwich ``df @ s_mn @ df.T``.
    Used as ``mn2mvarb``'s starting point, not its answer -- linear in a
    problem that mostly is not.
    """
    m1, m2, m3 = mn
    mc = mn2m(mn)
    denom = -(m1**2) + m2
    df = np.array(
        [
            [1.0, 0.0, 0.0],
            [-2.0 * m1, 1.0, 0.0],
            [
                (6.0 * m1**2 - 3.0 * m2) / denom**1.5
                + (3.0 * m1 * (2.0 * m1**3 - 3.0 * m1 * m2 + m3)) / denom**2.5,
                -3.0 * m1 / denom**1.5
                - (3.0 * (2.0 * m1**3 - 3.0 * m1 * m2 + m3)) / (2.0 * denom**2.5),
                denom ** (-1.5),
            ],
        ]
    )
    s_mc = df @ s_mn @ df.T
    return mc, s_mc


def mn2mvarb(mn: np.ndarray, s_mn: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """(mc, s_mc) such that ``mc2mnvb(mc, s_mc)`` reproduces ``s_mn``.

    ``emafit.f`` ``MN2MVARB`` (line 2514). See the module docstring for how
    the solve itself differs from the Fortran's step-halved Newton
    iteration while targeting the same root from the same starting point.
    """
    mc, s_mc0 = mn2m_var(mn, s_mn)
    v0 = _chol33(s_mc0)
    if v0 is None:
        raise ValueError("mn2mvarb: mn2m_var's initial covariance is not positive semi-definite")
    x0 = np.array([v0[0, 0], v0[1, 0], v0[0, 2], v0[1, 1], v0[1, 2], v0[2, 2]])
    target = _tri_in(s_mn)

    def residual(x: np.ndarray) -> np.ndarray:
        s_mc = _v_from_tri(x).T @ _v_from_tri(x)
        return _tri_in(mc2mnvb(mc, s_mc)) - target

    solution = root(residual, x0, method="hybr", tol=1e-10)
    if not solution.success:
        raise RuntimeError(f"mn2mvarb: root-find did not converge: {solution.message}")
    v = _v_from_tri(solution.x)
    return mc, v.T @ v


def mse_ema(nobs: np.ndarray, tl: np.ndarray, tu: np.ndarray, mc: np.ndarray, kmom: int) -> float:
    """MSE of the kmom-th at-site EMA moment (1=mean, 2=variance, 3=skew).

    ``emafit.f`` ``mse_ema`` (line 1659): ``var_mom -> m2mn -> mn2mvarb``,
    the ``(kmom, kmom)`` diagonal entry. This is what ADJE's censoring bias
    adjustment needs (TODO.md P3): ``bias_adj = mse_ema(censored) /
    mse_ema(uncensored)``, feeding ``as_G_mse = bias_adj * mseg(...)``.

    Near skew = 0, blends the ``skew = -_SKEWMIN``/``+_SKEWMIN`` results
    linearly rather than evaluating at the true (clamped) skew directly --
    a second, smoother blend on top of ``var_mom``'s own internal clamp,
    which is a no-op at exactly ``+/-_SKEWMIN``.
    """
    if kmom < 1 or kmom > 3:
        raise ValueError("mse_ema: kmom must be between 1 and 3")

    nobs = np.asarray(nobs, dtype=float)
    tl2 = np.asarray(tl, dtype=float) - mc[0]
    tu2 = np.asarray(tu, dtype=float) - mc[0]
    k = kmom - 1

    def _diag(skew: float) -> float:
        mc2 = np.array([0.0, mc[1], skew])
        s_mn = var_mom(nobs, tl2, tu2, mc2)
        mn = m2mn(mc2)
        _, s_mc = mn2mvarb(mn, s_mn)
        return float(s_mc[k, k])

    if abs(mc[2]) > _SKEWMIN:
        return _diag(mc[2])

    tneg = _diag(-_SKEWMIN)
    tpos = _diag(_SKEWMIN)
    w = (mc[2] - _SKEWMIN) / (2.0 * _SKEWMIN)
    return w * tneg + (1.0 - w) * tpos
