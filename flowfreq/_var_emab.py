"""Native port of peakfq's VAR_EMAB/regmoms/ci_ema_m3b: the CI-shape fix.

TODO.md P3's last open item. ``compute_confidence_limits()`` forms
``log_Q +/- z*se``, symmetric by construction; peakfq skews right with
return period via Cohn's inverse-Gaussian-quadrature method (``VAR_EMAB``,
``emafit.f:1972``).

Nearly everything this needs was already ported by earlier phases: ``var_mom``,
``m2mn``, ``mn2mvarb`` (``flowfreq._var_mom``, ``flowfreq._p3_moments``,
``flowfreq._mse_ema``), the Fortran-verified ``_gridmake``/``_covw``
(``flowfreq._mse_ema`` -- already indirectly oracle-tested through
``mc2mnvb``, which is exactly ``GRIDMAKE`` + ``M2MN`` + ``COVW`` composed,
``emafit.f:2594``), ``q_p3`` (``flowfreq._p3_moments``), and the ADJE
skew-MSE machinery (``ExpectedMomentsAlgorithm._adje_bias_adjustment``,
``flowfreq.bulletin17c._b17b_skew_mse``). ``regmoms`` and ``ci_ema_m3b`` are
mostly glue code composing those pieces; ``var_emab`` itself is the one
genuinely new piece -- a nested quadrature that estimates not just the
quantile but the sampling variability of its own standard error.

Everything here is a direct, function-for-function transcription of
``vendor/peakfqr/src/emafit.f``'s ``REGMOMS``, ``VAR_EMAB`` and
``CI_EMA_M3B``, checked against three Fortran oracles
``build_fortran/_emafort.pyf`` exposes (``regmoms``, ``var_emab``,
``ci_ema_m3b``) -- see ``tests/fortran_parity/test_fortran_oracles.py``.

Call convention, worth stating explicitly since it cost time to pin down
(a wrong argument order gave a plausible-looking but silently wrong
answer -- point estimates matched exactly, confidence widths did not):

* ``pq`` is *non-exceedance* probability (``P(X <= x)``), not AEP --
  ``q_p3`` uses ``ndtri(q)`` directly. Callers pass ``pq = 1 - aep``.
* ``mc`` is the *final*, regionally-weighted (mean, variance, skew) -- the
  same one ``Bulletin17C.run_analysis()`` already produces -- not the
  at-site-only fit.
* ``r_m_mse``/``r_s2_mse`` are "no regional info" sentinels here always:
  flowfreq has no regional mean/variance inputs, only regional skew
  (``r_g_mse``, already ``Bulletin17C``'s existing MSE encoding --
  0 = generalized/no error, <0 = generalized with error, >0 = weighted,
  >=1e10 = station-only).
* ``regmoms`` and ``var_emab`` disagree on their own
  ``(r_G_mse, r_M_mse, r_S2, r_S2_mse)`` argument order -- ``emafit.f``'s
  own two signatures do -- and ``VAR_EMAB`` reorders when it calls
  ``REGMOMS`` internally (``emafit.f:2010-2011``). Mirrored here as two
  different call sites with two different orders, not "fixed" to agree.

Nothing in this module is wired into ``ExpectedMomentsAlgorithm``/
``Bulletin17C`` yet.
"""

from __future__ import annotations

import math
from typing import Tuple

import numpy as np
from scipy.stats import t as _student_t

from flowfreq._mse_ema import _covw, _gridmake, mn2mvarb
from flowfreq._p3_moments import m2mn, q_p3
from flowfreq._var_mom import var_mom

__all__ = ["regmoms", "var_emab", "ci_ema_m3b"]

#: emafit.f:2198 -- regmoms's own skew clamp before calling anything that
#: needs alpha = 4/skew**2. Looser than var_mom's internal [0.0632, 1.41]
#: (applied again, redundantly but harmlessly, on top of this one).
_REGMOMS_SKEW_CLAMP = 1.5

#: "no regional info" sentinel, matching Bulletin17C's own regional-skew
#: MSE convention (see flowfreq.bulletin17c's module docstring): >= 1e10
#: means station-only, no blending.
NO_REGIONAL_INFO = 1e10

#: emafit.f:1869 -- ci_ema_m3b's own floors on nu and the denominator.
_NU_MIN = 5.0
_C_MIN = 0.5


def _mseg_all(nobs: np.ndarray, tl: np.ndarray, tu: np.ndarray, mc: np.ndarray) -> float:
    """ADJE's at-site skew MSE. ``emafit.f`` ``mseg_all`` (line 1569), 'ADJE' branch.

    Reuses ``ExpectedMomentsAlgorithm._adje_bias_adjustment`` (already
    Fortran-verified, TODO.md P3's "Skew weighting" item) rather than a
    second transcription of the same ``bias_adj`` computation. ``tl``/``tu``
    here are already mean-shifted (``regmoms``'s own convention), which
    matches what ``_adje_bias_adjustment`` expects: it calls ``mse_ema``,
    which re-derives its own shift from ``mc[0]`` -- exactly 0 here.
    """
    from flowfreq.bulletin17c import ExpectedMomentsAlgorithm, _b17b_skew_mse

    n = int(round(float(np.sum(nobs))))
    n_adj = min(n, 150)
    bias_adj = ExpectedMomentsAlgorithm._adje_bias_adjustment(
        tuple(np.asarray(nobs, dtype=float)),
        tuple(np.asarray(tl, dtype=float)),
        tuple(np.asarray(tu, dtype=float)),
        float(mc[0]),
        float(mc[1]),
        float(mc[2]),
        n,
    )
    return bias_adj * _b17b_skew_mse(n_adj, float(mc[2]))


def regmoms(
    nobs: np.ndarray,
    tl_in: np.ndarray,
    tu_in: np.ndarray,
    mc_in: np.ndarray,
    r_g_mse: float,
    r_m_mse: float = NO_REGIONAL_INFO,
    r_s2: float = 0.0,
    r_s2_mse: float = NO_REGIONAL_INFO,
) -> np.ndarray:
    """Covariance of the fitted central moments (mean, variance, skew).

    ``emafit.f`` ``regmoms`` (line 2173): ``var_mom`` -> ``m2mn`` ->
    ``mn2mvarb`` for the at-site covariance, then blended against regional
    mean/variance/skew information if supplied. Argument order is the
    Fortran's own: ``(r_g_mse, r_m_mse, r_s2, r_s2_mse)``.
    """
    nobs = np.asarray(nobs, dtype=float)
    mean0 = float(mc_in[0])
    tl = np.asarray(tl_in, dtype=float) - mean0
    tu = np.asarray(tu_in, dtype=float) - mean0
    mc = np.array(
        [0.0, float(mc_in[1]), max(-_REGMOMS_SKEW_CLAMP, min(_REGMOMS_SKEW_CLAMP, float(mc_in[2])))]
    )

    s_mn = var_mom(nobs, tl, tu, mc)
    mn = m2mn(mc)
    _, s_mc = mn2mvarb(mn, s_mn)
    s_mc = s_mc.copy()

    if r_g_mse > 0.0:
        w_g = r_g_mse / (_mseg_all(nobs, tl, tu, mc) + r_g_mse)
    elif -98.0 <= r_g_mse <= 0.0:
        w_g = 0.0
    else:  # r_g_mse < -98, or the station-only sentinel r_g_mse >= 1e10
        w_g = 1.0

    s_mc[0, 2] = s_mc[2, 0] = w_g * s_mc[0, 2]
    s_mc[1, 2] = s_mc[2, 1] = w_g * s_mc[1, 2]
    s_mc[2, 2] = w_g**2 * s_mc[2, 2] + (1.0 - w_g) ** 2 * abs(r_g_mse)

    if 0.0 <= r_m_mse < NO_REGIONAL_INFO:
        w_m = 0.0 if r_m_mse == 0.0 else r_m_mse / (s_mc[0, 0] + r_m_mse)
        s_mc[0, 1] = s_mc[1, 0] = w_m * s_mc[0, 1]
        s_mc[0, 2] = s_mc[2, 0] = w_m * s_mc[0, 2]
        s_mc[0, 0] = w_m**2 * s_mc[0, 0] + (1.0 - w_m) ** 2 * abs(r_m_mse)

    if r_s2_mse > NO_REGIONAL_INFO:
        return s_mc

    if r_s2_mse >= NO_REGIONAL_INFO or r_s2_mse < 0.0:
        w_s2, rs2mse = 1.0, NO_REGIONAL_INFO
    elif r_s2_mse == 0.0:
        w_s2, rs2mse = 0.0, r_s2_mse
    else:
        rs2mse = r_s2_mse * (mc[1] / r_s2) ** 2
        w_s2 = rs2mse / (s_mc[1, 1] + rs2mse)

    s_mc[0, 1] = s_mc[1, 0] = w_s2 * s_mc[0, 1]
    s_mc[1, 2] = s_mc[2, 1] = w_s2 * s_mc[1, 2]
    s_mc[1, 1] = w_s2**2 * s_mc[1, 1] + (1.0 - w_s2) ** 2 * abs(rs2mse)
    return s_mc


def ci_ema_m3b(yp: float, cv_yp_syp: np.ndarray, eps: float) -> Tuple[float, float, float]:
    """Asymmetric Student-t confidence bound. ``emafit.f`` ``ci_ema_m3b`` (line 1853).

    ``beta1``, the regression of the quantile's own standard error on the
    quantile itself, is the whole source of the asymmetry: it shrinks the
    denominator on one side and grows it on the other.

    Returns ``(ci_low, ci_high, nu)`` -- ``nu`` here is ``var(yp)``
    (``cv_yp_syp[0, 0]``), not the Student-t degrees of freedom used to get
    there. The Fortran overwrites its own ``nu`` output with that value
    right before returning (``emafit.f:1902``); transcribed faithfully
    rather than "fixed", since callers may rely on it.
    """
    var_yp = float(cv_yp_syp[0, 0])
    cov_yp_syp = float(cv_yp_syp[0, 1])
    var_syp = float(cv_yp_syp[1, 1])

    beta1 = cov_yp_syp / var_yp
    var_xsi_d = var_syp - cov_yp_syp**2 / var_yp
    nu = max(0.5 * var_yp / var_xsi_d, _NU_MIN)

    p_high = (1.0 + eps) / 2.0
    t = float(_student_t.ppf(p_high, nu))
    ci_high = yp + math.sqrt(var_yp) * t / max(_C_MIN, 1.0 - beta1 * t)
    t = -t
    # Same denominator formula as ci_high, ``1 - beta1*t`` -- the whole
    # asymmetry comes from t having just been negated, not from a
    # different formula (emafit.f:1897-1899 uses "1.d0-beta1*t" both times).
    ci_low = yp + math.sqrt(var_yp) * t / max(_C_MIN, 1.0 - beta1 * t)

    return ci_low, ci_high, var_yp


def var_emab(
    nobs: np.ndarray,
    tl: np.ndarray,
    tu: np.ndarray,
    mc: np.ndarray,
    pq: np.ndarray,
    eps: float,
    r_s2: float = 0.0,
    r_m_mse: float = NO_REGIONAL_INFO,
    r_s2_mse: float = NO_REGIONAL_INFO,
    r_g_mse: float = NO_REGIONAL_INFO,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Quantile estimates, their covariance with their own SE, and CIs.

    ``emafit.f`` ``VAR_EMAB`` (line 1972): Cohn's inverse-Gaussian-
    quadrature method. Argument order here is the Fortran's own --
    ``(r_s2, r_m_mse, r_s2_mse, r_g_mse)`` -- not ``regmoms``'s; ``regmoms``
    is called internally with its own order (``emafit.f:2010-2011``).

    Builds one outer 8-point quadrature grid around ``mc`` (from
    ``regmoms(mc)``'s covariance), then -- the genuinely new part -- one
    fresh ``regmoms``/inner-grid pair *at each of those 8 outer points*, so
    the resulting ``cv_yp_syp`` captures how both the quantile estimate and
    its own standard error co-vary under parameter uncertainty. This is
    what ``ci_ema_m3b`` needs to produce an asymmetric interval; it is also
    the expensive part -- 9 ``regmoms`` calls total (1 + 8), each a full
    ``var_mom``/``mn2mvarb`` solve -- but that cost does not depend on how
    many quantiles are requested: the grids are quantile-independent, built
    once and reused for every ``pq[k]``.

    Returns
    -------
    yp, cv_yp_syp, cil, cih, nu : np.ndarray
        ``cv_yp_syp`` has shape ``(len(pq), 2, 2)`` -- transposed from the
        Fortran's own ``(2, 2, nq)``, a Pythonic-indexing choice made here,
        not a mismatch to correct for when comparing against the oracle.
    """
    nobs = np.asarray(nobs, dtype=float)
    tl = np.asarray(tl, dtype=float)
    tu = np.asarray(tu, dtype=float)
    mc = np.asarray(mc, dtype=float)
    pq = np.asarray(pq, dtype=float)

    s_mc = regmoms(nobs, tl, tu, mc, r_g_mse, r_m_mse, r_s2, r_s2_mse)
    w1, gr_mc1 = _gridmake(mc, s_mc)
    n_outer = len(w1)

    w2 = np.empty((n_outer, n_outer))
    gr_mc2 = np.empty((n_outer, n_outer, 3))
    for i in range(n_outer):
        s_mc_i = regmoms(nobs, tl, tu, gr_mc1[i], r_g_mse, r_m_mse, r_s2, r_s2_mse)
        w2[i], gr_mc2[i] = _gridmake(gr_mc1[i], s_mc_i)

    nq = len(pq)
    yp = np.empty(nq)
    cv_yp_syp = np.empty((nq, 2, 2))
    cil = np.empty(nq)
    cih = np.empty(nq)
    nu = np.empty(nq)

    for k in range(nq):
        p = float(pq[k])
        qp1 = np.array([q_p3(p, gr_mc1[i]) for i in range(n_outer)])
        sp1 = np.empty(n_outer)
        for i in range(n_outer):
            qp2_i = np.array([q_p3(p, gr_mc2[i, j]) for j in range(n_outer)])
            sp1[i] = math.sqrt(_covw(qp2_i, qp2_i, w2[i]))

        yp[k] = q_p3(p, mc)
        cv_yp_syp[k, 0, 0] = _covw(qp1, qp1, w1)
        cv_yp_syp[k, 0, 1] = cv_yp_syp[k, 1, 0] = _covw(qp1, sp1, w1)
        cv_yp_syp[k, 1, 1] = _covw(sp1, sp1, w1)
        cil[k], cih[k], nu[k] = ci_ema_m3b(yp[k], cv_yp_syp[k], eps)

    return yp, cv_yp_syp, cil, cih, nu
