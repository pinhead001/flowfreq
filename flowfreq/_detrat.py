"""Native port of peakfq's detrat: the Halloween determinant ratio, Wd.

TODO.md P3's second open item, independent of the ``var_mom``/``mse_ema``
port (``flowfreq._var_mom``/``flowfreq._mse_ema``). ``Wd`` is HWN's other
weighting input, alongside ``as_G_mse``: ``emafit.f``::

    nG = n * Wd * as_G_mse / r_G_mse

``Wd`` corrects for something ``as_G_mse`` alone does not capture: on a
censored record, the at-site mean, variance, and skew estimates are
correlated with each other, so the skew estimate carries less *independent*
information than its own marginal MSE suggests. ``detrat`` measures that
via a determinant ratio -- ``det(I - F)`` for the full (mean, variance,
skew) system versus the (mean, variance)-only subsystem -- rather than
treating skew's uncertainty in isolation. Below an at-site skew magnitude of
0.04 the Fortran short-circuits to ``Wd = 1`` (``emafit.f:3654``); flowfreq
matches that, so this module is only ever exercised above that floor.

Everything here is a direct, function-for-function transcription of
``vendor/peakfqr/src/emafit.f``'s ``detrat`` and ``vendor/peakfqr/src/probfun.f``'s
``EXPMOMCDERIV`` (Greg Schwarz's derivation, cited in-line as "EQ n" the way
the Fortran comments do), checked against two Fortran oracles
``build_fortran/_emafort.pyf`` exposes (``expmomcderiv``, and Phase 1's
``detratsub``) -- see ``tests/fortran_parity/test_fortran_oracles.py``.

``EXPMOMCDERIV`` reuses ``flowfreq._var_mom._dexpect`` for the open-tail
expected moments and their Jacobian -- the same building block
``expmomderiv``/``d_est`` use -- rather than a second implementation of the
same truncated-gamma machinery.

Nothing in this module is wired into ``Bulletin17C``/``ExpectedMomentsAlgorithm``
yet.
"""

from __future__ import annotations

import math

import numpy as np

from flowfreq._p3_moments import _fp_g3_cdf, m2p
from flowfreq._var_mom import _dexpect

__all__ = ["detrat"]

#: emafit.f:3654 -- below this at-site skew magnitude, Wd = 1 outright.
#: Matches bulletin17c.py's _HWN_SKEW_FLOOR.
_SKEW_FLOOR = 0.04

#: probfun.f:916 -- EXPMOMCDERIV's own "infinity", log10(1e20). A different
#: sentinel from flowfreq._var_mom's _D_EST_INF (1e19): this one is in the
#: same log10(flow) units as tl/tu themselves, not an arbitrarily large
#: number, and must match exactly for the open-tail Jacobian to line up
#: with what the Fortran computes.
_INF = 20.0

#: probfun.f:940 -- below this censoring probability, EXPMOMCDERIV falls
#: back to the unconditional (-inf, inf) moments/Jacobian rather than
#: dividing by a near-zero PICK.
_PICK_FLOOR = 1e-10


def _p2m(parms: np.ndarray) -> np.ndarray:
    """Pearson III params (tau, alpha, beta) -> central moments (mean, var, skew).

    ``probfun.f`` ``P2M`` (line 752): the direct gamma-distribution moment
    formulas (variance ``alpha*beta**2``, skew ``2/sqrt(alpha)`` signed by
    beta), not ``m2p``'s inverse via the noncentral-moment route -- though
    the two agree exactly, since both describe the same distribution.
    """
    tau, alpha, beta = float(parms[0]), float(parms[1]), float(parms[2])
    mean = alpha * beta + tau
    var = alpha * beta**2
    skew = 2.0 / math.sqrt(alpha)
    if beta < 0.0:
        skew = -skew
    return np.array([mean, var, skew])


def _expmomcderiv(parms: np.ndarray, tl: float, tu: float):
    """E[X^j | X outside [tl, tu]] and the Jacobian of the resulting expected

    central moments w.r.t. the fit's own central moments (mean, variance,
    skew). ``probfun.f`` ``EXPMOMCDERIV`` (line 872), Schwarz's EQ 27.

    "Outside [tl, tu]" -- the *censored* region for this threshold group,
    not the perceived one -- is deliberate: this asks what the moments
    would look like restricted to exactly the data ``detrat`` never gets to
    see, which is the whole point of a censoring bias correction.
    """
    mc = _p2m(parms)
    mc1, mc2 = float(mc[0]), float(mc[1])
    s = math.sqrt(mc2)
    tau, alpha, beta = float(parms[0]), float(parms[1]), float(parms[2])

    pil = _fp_g3_cdf(tl, parms)
    pig = 1.0 - _fp_g3_cdf(tu, parms)
    pick = pil + pig

    if pick > _PICK_FLOOR:
        mnel, jacl = _dexpect(tau, alpha, beta, -_INF, tl)
        mneg, jacg = _dexpect(tau, alpha, beta, tu, _INF)
        mne = (pil * mnel + pig * mneg) / pick
        jac = (pil * jacl + pig * jacg) / pick
    else:
        mne, jac = _dexpect(tau, alpha, beta, -_INF, _INF)

    db = np.array(
        [
            [0.0, 0.0, 0.0],
            [2.0 * mc1, 0.0, 0.0],
            [-3.0 * mc1**2 / s**3, 3.0 * mc1**3 / (2.0 * s**5), 0.0],
        ]
    )
    b = np.array(
        [
            [1.0, 0.0, 0.0],
            [-2.0 * mc1, 1.0, 0.0],
            [3.0 * mc1**2 / s**3, -3.0 * mc1 / s**3, 1.0 / s**3],
        ]
    )
    # d(tau, alpha, beta)/d(mean, variance, skew) directly -- not the same
    # basis as flowfreq._var_mom._dpdm, which differentiates w.r.t.
    # noncentral moments instead.
    skew = float(mc[2])
    tp = np.array(
        [
            [1.0, -1.0 / (s * skew), 2.0 * s / skew**2],
            [0.0, 0.0, -8.0 / skew**3],
            [0.0, skew / (4.0 * s), s / 2.0],
        ]
    )
    d3 = np.array(
        [
            [0.0, 0.0, 0.0],
            [-2.0 * mne[0], 0.0, 0.0],
            [
                6.0 * mc1 / s**3 * mne[0] - 3.0 / s**3 * mne[1],
                (-4.5 * mc1**2 * mne[0] + 4.5 * mc1 * mne[1] - 1.5 * mne[2]) / s**5,
                0.0,
            ],
        ]
    )

    dedmc = db + b @ jac @ tp + d3
    return mc, mne, dedmc


def detrat(mc: np.ndarray, n: int, nobs: np.ndarray, tl: np.ndarray, tu: np.ndarray) -> float:
    """The Halloween determinant ratio, Wd. ``emafit.f`` ``detrat`` (line 3615).

    Parameters
    ----------
    mc : central moments (mean, variance, skew) of the at-site fit.
    n : total record length (not the number of threshold groups).
    nobs, tl, tu : perception-threshold groups, ``var_mom``'s convention --
        ``nobs[i]`` observations share the pair ``(tl[i], tu[i])``.

    Returns
    -------
    float
        ``1.0`` below the 0.04 at-site skew floor; otherwise
        ``det(I - F) / det(I - F[:2, :2])``, the 3x3-vs-2x2 determinant
        ratio.
    """
    skew = float(mc[2])
    if abs(skew) < _SKEW_FLOOR:
        return 1.0

    parms = m2p(mc)
    mc1, mc2 = float(mc[0]), float(mc[1])
    s = math.sqrt(mc2)

    f0phi = np.zeros((3, 3))
    for nobs_k, tl_k, tu_k in zip(nobs, tl, tu):
        pick = 1.0 - (_fp_g3_cdf(tu_k, parms) - _fp_g3_cdf(tl_k, parms))
        etak = float(nobs_k) / n

        _, mne, dck = _expmomcderiv(parms, tl_k, tu_k)

        # Expected central moments given censored, about the *population*
        # mean/variance (not the conditional mean) -- Schwarz's EQ 23.
        mcexpect = np.array(
            [
                mne[0],
                mc1**2 - 2.0 * mc1 * mne[0] + mne[1],
                (-(mc1**3) + 3.0 * mc1**2 * mne[0] - 3.0 * mc1 * mne[1] + mne[2]) / s**3,
            ]
        )

        f1 = np.zeros((3, 3))
        f1[1, 0] = 2.0 * pick * (mcexpect[0] - mc1)
        f1[2, 0] = 3.0 * (pick * mcexpect[1] - mc2) / s**3
        f1[2, 1] = 3.0 * (pick * mcexpect[2] - skew) / (2.0 * mc2)

        f0phi += (f1 + pick * dck) * etak

    imf = np.eye(3) - f0phi
    imful = imf[:2, :2]
    det2 = np.linalg.det(imful)
    det3 = np.linalg.det(imf)
    return det3 / det2
