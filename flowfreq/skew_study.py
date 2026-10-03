"""Regional skew development: the B-WLS/B-GLS procedure of Bulletin 17C.

Bulletin 17C (England and others, 2019, p. 31 and appendix 3) recommends a
regional skew estimated with the Bayesian weighted/generalized least squares
(B-WLS/B-GLS) procedure of Veilleux and others (2011), built on Reis and
others (2005). This module implements it as the USGS applies it, following the
equations of two reports that print them in full:

- SIR 2016-5083 appendix B (Veilleux; Pacific Northwest): pseudo record
  length (eq. B2), unbiasing (B5-B6), the Griffis-Stedinger variance (B7), the
  cross-correlation of skews (B8-B10), the Fisher-Z correlation-distance model
  (B11), EVR (B12) and MBV* (B13);
- SIR 2012-5130 appendix 3 (Lamontagne and others): the B-WLS weights
  (eqs. 3-6 to 3-8), the B-GLS posterior of the model error variance with
  its exponential prior, lambda = 10 (3-9, 3-10), the GLS precision of the
  WLS estimator (3-11, 3-12) and the variance of prediction (3-13a).

The procedure:

1. Unbias each station skew with its pseudo record length ``PRL`` (B5).
2. Fit OLS; the OLS regional skew at each site sets that site's sampling
   variance through the Griffis-Stedinger formula (B6-B7), which keeps the
   weights independent of the station skews themselves.
3. B-WLS: posterior mean of the model error variance with a diagonal
   covariance (beta integrated out under a flat prior), then WLS estimates of
   the model parameters.
4. B-GLS: posterior of the model error variance given the WLS parameters,
   with the full covariance -- sampling variances plus cross-correlations of
   concurrent skews -- and the GLS precision of the WLS parameters.

Diagnostics: model error variance (posterior mean and SD), parameter SDs,
ASEV, AVPnew (the Bulletin 17 MSE of the regional skew), its effective record
length, pseudo-R^2, EVR and MBV*.

Validated against the published Pacific Northwest CONSTANT model in
``tests/test_skew_study.py``; see ``docs/MONTANA_REGIONAL_SKEW_PROVISIONAL.md``.
Roadmap: ``docs/MASTER_ROADMAP.md`` section 1.3.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import brentq, curve_fit

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd

logger = logging.getLogger(__name__)

FloatArray = NDArray[np.float64]

#: Exponential prior rate for the model error variance (SIR 2012-5130 eq. 3-10:
#: lambda = 10, a prior mean of 0.1).
PRIOR_LAMBDA = 10.0

#: Martins and Stedinger (2002) exponent in rho(G_i, G_j) = sign * cf * |rho|^kappa;
#: SIR 2016-5083 eq. B8 gives it as "between 2.8 and 3.3".
DEFAULT_KAPPA = 2.8

EARTH_RADIUS_MILES = 3958.8


# ----------------------------------------------------------------------------
# Station skew sampling variance
# ----------------------------------------------------------------------------


def gs_skew_variance(n: ArrayLike, skew: ArrayLike) -> FloatArray:
    """Griffis and Stedinger (2009) variance of the sample skew (SIR 2016-5083 eq. B7).

    Parameters
    ----------
    n : array_like
        Record length (the pseudo record length with censored or historic data).
    skew : array_like
        Skew at which to evaluate it -- in the regional procedure, the OLS
        regional estimate, not the station skew.

    Returns
    -------
    numpy.ndarray
    """
    n = np.asarray(n, dtype=float)
    g = np.asarray(skew, dtype=float)
    if np.any(n <= 0):
        raise ValueError("record length must be positive")
    a = -17.75 / n**2 + 50.06 / n**3
    b = 3.92 / n**0.3 - 31.10 / n**0.6 + 34.86 / n**0.9
    c = -7.31 / n**0.59 + 45.90 / n**1.18 - 86.50 / n**1.77
    return np.asarray((6.0 / n + a) * (1.0 + (9.0 / 6.0 + b) * g**2 + (15.0 / 48.0 + c) * g**4))


def unbiased_skew_variance(n: ArrayLike, skew: ArrayLike) -> FloatArray:
    """Variance of the unbiased station skew, ``(1 + 6/n)^2 Var[G]`` (eq. B6)."""
    n = np.asarray(n, dtype=float)
    return np.asarray((1.0 + 6.0 / n) ** 2 * gs_skew_variance(n, skew))


def unbias_skew(skew: ArrayLike, n: ArrayLike) -> FloatArray:
    """Tasker-Stedinger unbiased skew, ``(1 + 6/n) G`` (eq. B5)."""
    return np.asarray((1.0 + 6.0 / np.asarray(n, dtype=float)) * np.asarray(skew, dtype=float))


def effective_record_length(mse: float, skew: float) -> float:
    """Years of record whose unbiased station-skew variance equals ``mse``.

    With the unbiased form this reproduces the Pacific Northwest study's
    "AVPnew 0.18 ... effective record length of 41 years" (SIR 2016-5083
    p. 52) exactly. (The same report's "17 years" for the Bulletin 17B map's
    0.302 comes from the biased form, :func:`gs_skew_variance`.)

    Parameters
    ----------
    mse : float
        Variance of prediction of the regional skew.
    skew : float
        Regional skew.

    Returns
    -------
    float
    """
    if not mse > 0:
        raise ValueError(f"mse must be positive, got {mse}")
    f = lambda n: float(unbiased_skew_variance(n, skew)) - mse  # noqa: E731
    if f(1e5) > 0:
        return float("inf")
    return float(brentq(f, 3.0, 1e5))


def pseudo_record_length(
    n_systematic: float,
    mse_systematic: float,
    mse_all: float,
    historical_period: Optional[float] = None,
) -> float:
    """Pseudo record length of the skew, ``Ps * MSE(G_S) / MSE(G_C)`` (eq. B2).

    Clipped as the report requires: at least ``Ps`` and at most the historical
    period ``PH``.

    Parameters
    ----------
    n_systematic : float
        Number of systematic peaks, ``Ps``.
    mse_systematic : float
        MSE of the skew from the systematic record alone.
    mse_all : float
        MSE of the skew from all data, historic and censored included.
    historical_period : float, optional
        ``PH``, years.

    Returns
    -------
    float
    """
    if n_systematic <= 0 or mse_systematic <= 0 or mse_all <= 0:
        raise ValueError("n_systematic and both MSEs must be positive")
    prl = n_systematic * mse_systematic / mse_all
    if historical_period is not None:
        prl = min(prl, float(historical_period))
    return float(max(prl, n_systematic))


# ----------------------------------------------------------------------------
# Cross-correlation of concurrent skews
# ----------------------------------------------------------------------------


def distance_miles(lat: ArrayLike, lon: ArrayLike) -> FloatArray:
    """Great-circle distance matrix between points, miles."""
    la = np.radians(np.asarray(lat, dtype=float))
    lo = np.radians(np.asarray(lon, dtype=float))
    dla = la[:, None] - la[None, :]
    dlo = lo[:, None] - lo[None, :]
    h = np.sin(dla / 2) ** 2 + np.cos(la[:, None]) * np.cos(la[None, :]) * np.sin(dlo / 2) ** 2
    return np.asarray(2 * EARTH_RADIUS_MILES * np.arcsin(np.sqrt(np.clip(h, 0, 1))))


@dataclass(frozen=True)
class CorrelationDistanceModel:
    """Fisher-Z cross-correlation of concurrent log peaks vs. distance (eq. B11).

    ``Z = a + exp(b + c * D)``, ``rho = tanh(Z)``, D in miles. The Pacific
    Northwest study's is ``a = 0.21, b = -0.17, c = -0.0058``.
    """

    a: float
    b: float
    c: float
    n_pairs: int = 0
    min_concurrent: int = 0

    def rho(self, d: ArrayLike) -> FloatArray:
        """Modelled correlation of concurrent annual peaks at distance ``d`` miles."""
        d = np.asarray(d, dtype=float)
        return np.asarray(np.tanh(self.a + np.exp(self.b + self.c * d)))

    @classmethod
    def fit(
        cls, distance: ArrayLike, correlation: ArrayLike, min_concurrent: int = 0
    ) -> "CorrelationDistanceModel":
        """Least-squares fit of ``Z = a + exp(b + c D)`` to sample correlations.

        Parameters
        ----------
        distance : array_like
            Pair distances, miles.
        correlation : array_like
            Sample correlations of concurrent log peaks.
        min_concurrent : int
            Recorded with the model: the concurrent-years screen used.

        Returns
        -------
        CorrelationDistanceModel

        Raises
        ------
        ValueError
            Fewer than 3 pairs, or the fit fails.
        """
        d = np.asarray(distance, dtype=float)
        r = np.clip(np.asarray(correlation, dtype=float), -0.999, 0.999)
        if len(d) < 3 or len(d) != len(r):
            raise ValueError("need at least 3 aligned (distance, correlation) pairs")
        z = np.arctanh(r)
        try:
            p, _ = curve_fit(
                lambda x, a, b, c: a + np.exp(b + c * x),
                d,
                z,
                p0=(0.2, -0.2, -0.005),
                maxfev=20000,
            )
        except RuntimeError as exc:
            raise ValueError(f"correlation-distance fit failed: {exc}") from exc
        return cls(
            float(p[0]), float(p[1]), float(p[2]), n_pairs=len(d), min_concurrent=min_concurrent
        )


#: SIR 2016-5083 eq. B11 (Pacific Northwest, 33 sites with >= 85 concurrent years).
PNW_CORRELATION_MODEL = CorrelationDistanceModel(
    0.21, -0.17, -0.0058, n_pairs=411, min_concurrent=85
)


def concurrent_pseudo_record_length(
    yb: ArrayLike, ye: ArrayLike, prl: ArrayLike, ph: ArrayLike
) -> FloatArray:
    """Pseudo concurrent record length matrix (eq. B10).

    ``CY_ij = (YE_ij - YB_ij + 1) (PRL_i / PH_i) (PRL_j / PH_j)`` over the
    years common to both historical periods; zero where they do not overlap.
    ``PRL/PH`` is capped at 1.
    """
    yb = np.asarray(yb, dtype=float)
    ye = np.asarray(ye, dtype=float)
    ratio = np.minimum(np.asarray(prl, dtype=float) / np.asarray(ph, dtype=float), 1.0)
    overlap = np.minimum(ye[:, None], ye[None, :]) - np.maximum(yb[:, None], yb[None, :]) + 1
    return np.asarray(np.clip(overlap, 0, None) * ratio[:, None] * ratio[None, :])


def skew_correlation_matrix(
    lat: ArrayLike,
    lon: ArrayLike,
    prl: ArrayLike,
    yb: ArrayLike,
    ye: ArrayLike,
    ph: ArrayLike,
    model: CorrelationDistanceModel,
    kappa: float = DEFAULT_KAPPA,
) -> FloatArray:
    """Correlation of the station skew estimators (eqs. B8-B11).

    ``rho(G_i, G_j) = sign(rho_ij) * cf_ij * |rho_ij|^kappa`` with
    ``cf_ij = CY_ij / sqrt(PRL_i PRL_j)`` and ``rho_ij`` from the
    correlation-distance model. Unit diagonal.
    """
    prl = np.asarray(prl, dtype=float)
    rho_q = model.rho(distance_miles(lat, lon))
    cf = concurrent_pseudo_record_length(yb, ye, prl, ph) / np.sqrt(np.outer(prl, prl))
    out = np.sign(rho_q) * np.clip(cf, 0, 1) * np.abs(rho_q) ** kappa
    np.fill_diagonal(out, 1.0)
    return np.asarray(out)


# ----------------------------------------------------------------------------
# B-WLS / B-GLS
# ----------------------------------------------------------------------------


@dataclass
class RegionalSkewModel:
    """A fitted B-WLS/B-GLS regional skew model and its diagnostics.

    Attributes
    ----------
    names : list of str
        Parameter names (``"constant"`` first).
    beta : ndarray
        B-WLS parameter estimates.
    beta_sd : ndarray
        Their standard deviations from the B-GLS precision analysis.
    sigma2 : float
        Posterior mean of the model error variance (B-GLS).
    sigma2_sd : float
        Its posterior standard deviation.
    sigma2_wls : float
        Posterior mean of the model error variance under B-WLS (used for the
        weights).
    asev : float
        Average sampling error variance of the regional estimator at the sites.
    avp_new : float
        Average variance of prediction at a new site: the Bulletin 17 MSE.
    effective_record_length : float
        Years of record equivalent to ``avp_new`` (:func:`effective_record_length`).
    pseudo_r2 : float
        ``1 - sigma2(k) / sigma2(0)``; 0 for the constant model.
    evr : float
        Error variance ratio, sum Var(gamma_i) / (n sigma2) (eq. B12).
    mbv : float
        Misrepresentation of the beta variance, MBV* (eq. B13).
    n_sites : int
    sampling_variance : ndarray
        ``Var[gamma_i]`` per site (eq. B6).
    residuals : ndarray
        Unbiased skew minus the model's regional skew.
    leverage, influence : ndarray
        Diagnostics (Veilleux 2011): GLS hat-matrix diagonal and Cook's-D style
        influence; used to flag unusual sites.
    kappa : float, optional
        Exponent used for the skew cross-correlations, if any.
    """

    names: List[str]
    beta: FloatArray
    beta_sd: FloatArray
    sigma2: float
    sigma2_sd: float
    sigma2_wls: float
    asev: float
    avp_new: float
    effective_record_length: float
    pseudo_r2: float
    evr: float
    mbv: float
    n_sites: int
    sampling_variance: FloatArray
    residuals: FloatArray
    leverage: FloatArray
    influence: FloatArray
    kappa: Optional[float] = None
    extra: Dict[str, Any] = field(default_factory=dict)

    def regional_skew(self, x: Optional[Sequence[float]] = None) -> float:
        """Regional skew at a site with explanatory values ``x`` (excluding the constant)."""
        xv = np.r_[1.0, np.asarray(x if x is not None else [], dtype=float)]
        if len(xv) != len(self.beta):
            raise ValueError(f"expected {len(self.beta) - 1} explanatory values")
        return float(xv @ self.beta)

    def summary(self) -> Dict[str, float]:
        """The headline numbers, as a flat dict."""
        out = {f"beta_{n}": float(b) for n, b in zip(self.names, self.beta)}
        out.update({f"sd_{n}": float(s) for n, s in zip(self.names, self.beta_sd)})
        out.update(
            sigma2=self.sigma2,
            sigma2_sd=self.sigma2_sd,
            asev=self.asev,
            avp_new=self.avp_new,
            effective_record_length=self.effective_record_length,
            pseudo_r2=self.pseudo_r2,
            evr=self.evr,
            mbv=self.mbv,
            n_sites=float(self.n_sites),
        )
        return out


def _grid(v: FloatArray, n: int = 4000) -> FloatArray:
    upper = max(3.0, 20.0 * float(np.median(v)))
    return np.linspace(0.0, upper, n)


def _integrate(y: FloatArray, x: FloatArray) -> float:
    """Trapezoid rule (``np.trapezoid`` is numpy >= 2 only)."""
    return float(np.sum((y[1:] + y[:-1]) * np.diff(x)) / 2.0)


def _posterior_moments(log_post: FloatArray, grid: FloatArray) -> Tuple[float, float]:
    w = np.exp(log_post - log_post.max())
    z = _integrate(w, grid)
    mean = float(_integrate(w * grid, grid) / z)
    var = float(_integrate(w * (grid - mean) ** 2, grid) / z)
    if w[-1] > 1e-6:
        logger.warning("model error variance posterior truncated at %.2f", grid[-1])
    return mean, float(np.sqrt(max(var, 0.0)))


def fit_regional_skew(
    skew: ArrayLike,
    prl: ArrayLike,
    explanatory: Optional[ArrayLike] = None,
    names: Optional[Sequence[str]] = None,
    *,
    correlation: Optional[ArrayLike] = None,
    prior_lambda: float = PRIOR_LAMBDA,
    kappa: Optional[float] = None,
    _constant_sigma2: Optional[float] = None,
) -> RegionalSkewModel:
    """Fit a B-WLS/B-GLS regional skew model.

    Parameters
    ----------
    skew : array_like
        Station skews ``G_i`` (biased, as EMA reports them).
    prl : array_like
        Their pseudo record lengths (:func:`pseudo_record_length`).
    explanatory : array_like, optional
        ``(n, k)`` explanatory variables, without the constant column. ``None``
        fits the CONSTANT model.
    names : sequence of str, optional
        Names of the explanatory columns.
    correlation : array_like, optional
        ``(n, n)`` correlation of the skew estimators
        (:func:`skew_correlation_matrix`). ``None`` treats them as independent,
        so the B-GLS step reduces to B-WLS precision (MBV* = 1).
    prior_lambda : float
        Rate of the exponential prior on the model error variance.
    kappa : float, optional
        Recorded on the result only.

    Returns
    -------
    RegionalSkewModel

    Raises
    ------
    ValueError
        Misaligned inputs, non-positive PRL, fewer sites than parameters + 2,
        or a correlation matrix that is not positive definite.
    """
    g = np.asarray(skew, dtype=float)
    p = np.asarray(prl, dtype=float)
    n = len(g)
    if p.shape != g.shape or g.ndim != 1:
        raise ValueError("skew and prl must be 1-D and aligned")
    if not np.all(np.isfinite(g)) or not np.all(p > 0):
        raise ValueError("skews must be finite and PRL positive")
    xe = np.empty((n, 0)) if explanatory is None else np.asarray(explanatory, dtype=float)
    if xe.ndim == 1:
        xe = xe[:, None]
    if xe.shape[0] != n:
        raise ValueError("explanatory must have one row per site")
    X = np.column_stack([np.ones(n), xe])
    k = X.shape[1]
    if n < k + 2:
        raise ValueError(f"{n} sites cannot fit {k} parameters")
    col_names = ["constant"] + list(names or [f"x{i + 1}" for i in range(xe.shape[1])])
    if len(col_names) != k:
        raise ValueError("names must name every explanatory column")

    gamma = unbias_skew(g, p)  # B5
    beta_ols = np.linalg.lstsq(X, gamma, rcond=None)[0]
    var = unbiased_skew_variance(p, X @ beta_ols)  # B6-B7 at the OLS regional skew

    # --- B-WLS: diagonal covariance, beta integrated out -------------------
    grid = _grid(var)
    logp = np.empty_like(grid)
    for i, s2 in enumerate(grid):
        lam = s2 + var
        wi = 1.0 / lam
        xtwx = X.T @ (X * wi[:, None])
        b = np.linalg.solve(xtwx, X.T @ (wi * gamma))
        r = gamma - X @ b
        logp[i] = (
            -prior_lambda * s2
            - 0.5 * np.sum(np.log(lam))
            - 0.5 * np.linalg.slogdet(xtwx)[1]
            - 0.5 * np.sum(r * r * wi)
        )
    s2_wls, _ = _posterior_moments(logp, grid)
    lam_wls = s2_wls + var
    xtwx = X.T @ (X / lam_wls[:, None])
    W = np.linalg.solve(xtwx, (X / lam_wls[:, None]).T)  # (k, n), eq. 3-7
    beta = W @ gamma  # eq. 3-8

    # --- B-GLS: full covariance, conditional on beta_WLS -------------------
    sd = np.sqrt(var)
    if correlation is None:
        sigma = np.diag(var)
    else:
        c = np.asarray(correlation, dtype=float)
        if c.shape != (n, n):
            raise ValueError("correlation must be (n, n)")
        sigma = c * np.outer(sd, sd)
        np.fill_diagonal(sigma, var)
    evals, evecs = np.linalg.eigh(sigma)
    if evals.min() <= -1e-10 * evals.max():
        raise ValueError(
            f"skew covariance is not positive semidefinite (min eigenvalue {evals.min():.3g})"
        )
    evals = np.clip(evals, 0.0, None)
    r = gamma - X @ beta
    rq2 = (evecs.T @ r) ** 2
    lam_grid = grid[:, None] + evals[None, :]
    logp_gls = (
        -prior_lambda * grid
        - 0.5 * np.sum(np.log(np.maximum(lam_grid, 1e-300)), axis=1)
        - 0.5 * np.sum(rq2[None, :] / np.maximum(lam_grid, 1e-300), axis=1)
    )
    if grid[0] == 0.0 and evals.min() == 0.0:
        logp_gls[0] = -np.inf
    s2, s2_sd = _posterior_moments(logp_gls, grid)

    lam_gls = s2 * np.eye(n) + sigma  # eq. 3-12 at the posterior mean
    var_beta = W @ lam_gls @ W.T  # eq. 3-11
    beta_sd = np.sqrt(np.diag(var_beta))
    asev = float(np.mean(np.einsum("ij,jk,ik->i", X, var_beta, X)))
    avp = s2 + asev  # eq. 3-13a averaged over the sites

    # Diagnostics (B12, B13).
    evr = float(np.sum(var) / (n * s2)) if s2 > 0 else float("inf")
    w0 = 1.0 / lam_wls
    mbv = float(w0 @ lam_gls @ w0 / np.sum(w0))

    # Leverage / influence (Veilleux 2011): hat matrix of the WLS fit under
    # the GLS covariance, and a Cook's-D analogue.
    H = X @ W
    lev = np.diag(H).copy()
    resid_var = np.diag(lam_gls) * (1 - lev)
    infl = (r**2 / np.maximum(resid_var, 1e-12)) * lev / (k * np.maximum(1 - lev, 1e-12))

    if k == 1:
        pr2 = 0.0
    else:
        s2_0 = _constant_sigma2
        if s2_0 is None:
            s2_0 = fit_regional_skew(
                g, p, correlation=correlation, prior_lambda=prior_lambda
            ).sigma2
        pr2 = max(0.0, 1.0 - s2 / s2_0) if s2_0 > 0 else 0.0

    mean_skew = float(np.mean(X @ beta))
    return RegionalSkewModel(
        names=col_names,
        beta=beta,
        beta_sd=beta_sd,
        sigma2=s2,
        sigma2_sd=s2_sd,
        sigma2_wls=s2_wls,
        asev=asev,
        avp_new=avp,
        effective_record_length=effective_record_length(avp, mean_skew),
        pseudo_r2=pr2,
        evr=evr,
        mbv=mbv,
        n_sites=n,
        sampling_variance=np.asarray(var),
        residuals=np.asarray(r),
        leverage=np.asarray(lev),
        influence=np.asarray(infl),
        kappa=kappa,
    )


# ----------------------------------------------------------------------------
# Station inputs from flowfreq's own EMA
# ----------------------------------------------------------------------------


@dataclass(frozen=True)
class StationSkew:
    """One station's inputs to a regional skew study, from an EMA/MGB fit."""

    site_no: str
    skew: float
    skew_mse: float
    prl: float
    n_systematic: int
    n_peaks: int
    n_low_outliers: int
    yb: int
    ye: int
    ph: int
    at_site_option: str


def station_skew(site_no: str, peaks: "pd.DataFrame") -> StationSkew:
    """Station skew, its MSE and pseudo record length from flowfreq's EMA/MGB.

    The record is fitted station-only (no regional skew) exactly as
    :func:`flowfreq.workflow.run_ffa` would: peakfq's ``siteQT`` treatment of
    the qualification codes (historic, censored, regulated peaks removed)
    and the Multiple Grubbs-Beck low-outlier test. The skew MSE is the fit's
    own ``as_G_mse`` -- ADJE, or plain B17B ``mseg`` when MGBT finds low
    outliers (``emafit.f:707``) -- and the PRL is the fit's ``as_G_PRL_o``
    (``Ps * MSE(G_S) / MSE(G_C)``, eq. B2), clipped to ``[Ps, PH]``.

    Parameters
    ----------
    site_no : str
    peaks : pandas.DataFrame
        ``water_year``, ``peak_flow_cfs``, ``qualification_code``.

    Returns
    -------
    StationSkew

    Raises
    ------
    ValueError
        The fit fails or yields no pseudo record length.
    """
    from flowfreq.bulletin17c import Bulletin17C
    from flowfreq.workflow import peak_code_kwargs

    df = peaks.dropna(subset=["peak_flow_cfs"]).sort_values("water_year")
    years = df["water_year"].astype(int).to_numpy()
    flows = df["peak_flow_cfs"].astype(float).to_numpy()
    codes = df["qualification_code"].fillna("").astype(str).tolist()
    kwargs = peak_code_kwargs(flows.tolist(), years.tolist(), codes, site_name=site_no) or {
        "peak_flows": flows,
        "water_years": years,
    }
    b = Bulletin17C(**kwargs)
    res = b.run_analysis(method="ema")
    analyzer: Any = b._analyzer
    mse = float(
        analyzer._at_site_skew_mse(res.mean_log, res.std_log, res.skew_station, res.n_peaks)
    )
    if res.pseudo_record_length is None:
        raise ValueError(f"{site_no}: no pseudo record length")
    yb, ye = int(years.min()), int(years.max())
    ph = ye - yb + 1
    n_sys = int(res.n_systematic)
    prl = float(min(max(res.pseudo_record_length, n_sys), ph))
    return StationSkew(
        site_no=site_no,
        skew=float(res.skew_station),
        skew_mse=mse,
        prl=prl,
        n_systematic=n_sys,
        n_peaks=int(res.n_peaks),
        n_low_outliers=int(res.n_low_outliers),
        yb=yb,
        ye=ye,
        ph=ph,
        at_site_option=str(analyzer._at_site_option),
    )
