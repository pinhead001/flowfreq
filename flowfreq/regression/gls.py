"""Generalized least squares (GLS) regression for flood quantiles.

The regression-development half of roadmap §3.3: the Stedinger and Tasker
(1985) GLS procedure as the USGS WREG package implements it (Eng, Chen and
Kiang, 2009, USGS Techniques and Methods 4-A8). Station log-quantiles ``y``
are regressed on basin characteristics ``X`` with an error covariance

    Lambda = gamma^2 I + Sigma

where ``gamma^2`` is the model error variance and ``Sigma`` the sampling
covariance of the station estimates (sampling variance on the diagonal,
cross-correlated concurrent records off it).

``gamma^2`` is the method-of-moments root of

    (y - X b)' Lambda^-1 (y - X b) = n - p,

with ``b`` the GLS estimate at that ``gamma^2``, as in WREG. When the residuals
are smaller than the stated sampling error alone explains, there is no
positive root and ``gamma^2 = 0``.

This module fits; it does not choose variables or regions. Distances are
statute miles (:func:`flowfreq.skew_study.distance_miles`).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np
import numpy.typing as npt
from scipy.optimize import brentq, least_squares
from scipy.stats import pearson3

from flowfreq.skew_study import distance_miles

logger = logging.getLogger(__name__)

FloatArray = npt.NDArray[np.float64]

#: Step for the central-difference derivative dK/dG of the Pearson III frequency factor.
_DG = 1e-4


def _frequency_factor(skew: FloatArray, aep: float) -> FloatArray:
    return np.asarray(pearson3.ppf(1.0 - aep, skew), dtype=float)


def lp3_quantile_sampling_variance(
    sd: npt.ArrayLike,
    skew: npt.ArrayLike,
    n: npt.ArrayLike,
    aep: float,
    skew_mse: float,
) -> FloatArray:
    """Sampling variance of a station log10 LP3 quantile estimator.

    The WREG form (Eng and others, 2009; after Stedinger and Tasker, 1986 and
    Griffis and Stedinger, 2007), with the skew treated as a regional
    (weighted) value of known MSE::

        Var = sd^2/n * [1 + K G + K^2/2 (1 + 3/4 G^2)] + sd^2 (dK/dG)^2 MSE_G

    Parameters
    ----------
    sd : array_like
        Standard deviation of the log10 annual peaks (ideally regionally
        smoothed, so the weights are not correlated with the residuals).
    skew : array_like
        Skew used for the station quantile (weighted or regional).
    n : array_like
        Record length, years (systematic or pseudo record length). Must be > 0.
    aep : float
        Annual exceedance probability, in (0, 1).
    skew_mse : float
        MSE of the skew; 0 ignores skew uncertainty.

    Returns
    -------
    numpy.ndarray
        Variance in log10 units squared, one per site.

    Raises
    ------
    ValueError
        ``aep`` outside (0, 1), a non-positive record length, or a negative MSE.
    """
    if not 0.0 < aep < 1.0:
        raise ValueError(f"aep must be in (0, 1), got {aep}")
    if skew_mse < 0:
        raise ValueError(f"skew_mse must be >= 0, got {skew_mse}")
    sd_a = np.asarray(sd, dtype=float)
    g = np.asarray(skew, dtype=float)
    n_a = np.asarray(n, dtype=float)
    if np.any(n_a <= 0):
        raise ValueError("record length must be positive")
    k = _frequency_factor(g, aep)
    dk = (_frequency_factor(g + _DG, aep) - _frequency_factor(g - _DG, aep)) / (2 * _DG)
    base = sd_a**2 / n_a * (1.0 + k * g + k**2 / 2.0 * (1.0 + 0.75 * g**2))
    return np.asarray(base + sd_a**2 * dk**2 * skew_mse, dtype=float)


@dataclass(frozen=True)
class WregCorrelation:
    """WREG's correlation-distance model: ``rho = theta ** (d / (alpha d + 1))``.

    Parameters
    ----------
    alpha : float
        Shape, >= 0 (per mile).
    theta : float
        Correlation base, in (0, 1).
    """

    alpha: float
    theta: float

    def __post_init__(self) -> None:
        if not 0.0 < self.theta < 1.0:
            raise ValueError(f"theta must be in (0, 1), got {self.theta}")
        if self.alpha < 0:
            raise ValueError(f"alpha must be >= 0, got {self.alpha}")

    def rho(self, d_miles: npt.ArrayLike) -> FloatArray:
        """Correlation of concurrent log peaks at distance ``d_miles``."""
        d = np.asarray(d_miles, dtype=float)
        return np.asarray(self.theta ** (d / (self.alpha * d + 1.0)), dtype=float)


def fit_correlation(d_miles: npt.ArrayLike, r: npt.ArrayLike) -> WregCorrelation:
    """Fit :class:`WregCorrelation` to station-pair correlations by least squares.

    Parameters
    ----------
    d_miles : array_like
        Inter-site distances, miles.
    r : array_like
        Sample correlations of concurrent log10 annual peaks (pairs already
        screened for a minimum overlap by the caller).

    Returns
    -------
    WregCorrelation

    Raises
    ------
    ValueError
        Fewer than 3 pairs, or mismatched lengths.
    """
    d = np.asarray(d_miles, dtype=float)
    rr = np.asarray(r, dtype=float)
    if d.shape != rr.shape or d.size < 3:
        raise ValueError("need at least 3 (distance, correlation) pairs of equal length")

    def resid(p: FloatArray) -> FloatArray:
        alpha, theta = p
        return np.asarray(theta ** (d / (alpha * d + 1.0)) - rr, dtype=float)

    sol = least_squares(resid, x0=[0.002, 0.98], bounds=([0.0, 1e-6], [1.0, 1.0 - 1e-9]))
    return WregCorrelation(alpha=float(sol.x[0]), theta=float(sol.x[1]))


def sampling_covariance(
    sd: npt.ArrayLike,
    skew: npt.ArrayLike,
    n: npt.ArrayLike,
    aep: float,
    skew_mse: float,
    lat: npt.ArrayLike,
    lon: npt.ArrayLike,
    concurrent: npt.ArrayLike,
    corr: WregCorrelation,
) -> FloatArray:
    """Sampling covariance matrix Sigma of the station log quantiles.

    Diagonal: :func:`lp3_quantile_sampling_variance`. Off-diagonal::

        Sigma_ij = rho_ij sqrt(Var_i Var_j) c_ij / sqrt(n_i n_j)

    with ``c_ij`` the concurrent years and ``rho_ij`` from ``corr`` at the
    sites' distance; zero where the records do not overlap.

    Parameters
    ----------
    sd, skew, n : array_like
        As :func:`lp3_quantile_sampling_variance`, one per site.
    aep, skew_mse : float
        As :func:`lp3_quantile_sampling_variance`.
    lat, lon : array_like
        Site (or basin centroid) coordinates, decimal degrees.
    concurrent : array_like
        Square matrix of concurrent years of record.
    corr : WregCorrelation
        Correlation-distance model.

    Returns
    -------
    numpy.ndarray
        Symmetric (m, m) matrix.

    Raises
    ------
    ValueError
        Inconsistent sizes.
    """
    var = lp3_quantile_sampling_variance(sd, skew, n, aep, skew_mse)
    m = var.size
    c = np.asarray(concurrent, dtype=float)
    n_a = np.asarray(n, dtype=float)
    if c.shape != (m, m) or n_a.size != m:
        raise ValueError(f"concurrent must be ({m}, {m}) and n length {m}")
    d = distance_miles(np.asarray(lat, dtype=float), np.asarray(lon, dtype=float))
    rho = corr.rho(d)
    s = rho * np.sqrt(np.outer(var, var)) * c / np.sqrt(np.outer(n_a, n_a))
    s = np.where(c > 0, s, 0.0)
    np.fill_diagonal(s, var)
    return np.asarray((s + s.T) / 2.0, dtype=float)


@dataclass
class GLSResult:
    """A fitted GLS regression (log10 space)."""

    beta: FloatArray
    cov_beta: FloatArray
    model_error_variance: float
    avp: float
    sep_percent: float
    pseudo_r2: float
    leverage: FloatArray
    influence: FloatArray
    residuals: FloatArray
    names: List[str] = field(default_factory=list)

    def predict(self, x: npt.ArrayLike) -> Tuple[float, float]:
        """Prediction and its variance at one site's explanatory row ``x``.

        Returns
        -------
        (float, float)
            ``x'b`` and ``gamma^2 + x' cov_beta x`` (variance of prediction).
        """
        xv = np.asarray(x, dtype=float)
        return float(xv @ self.beta), float(self.model_error_variance + xv @ self.cov_beta @ xv)


def _gls_at(
    y: FloatArray, X: FloatArray, sigma: FloatArray, gamma2: float
) -> Tuple[FloatArray, FloatArray, FloatArray]:
    lam = sigma + gamma2 * np.eye(len(y))
    # np.asarray(..., float64): newer numpy stubs type linalg.inv as floating[Any].
    lam_inv = np.asarray(np.linalg.inv(lam), dtype=np.float64)
    cov = np.asarray(np.linalg.inv(X.T @ lam_inv @ X), dtype=np.float64)
    beta = np.asarray(cov @ X.T @ lam_inv @ y, dtype=np.float64)
    return beta, cov, lam_inv


def _mev(y: FloatArray, X: FloatArray, sigma: FloatArray) -> float:
    n, p = X.shape

    def f(g2: float) -> float:
        beta, _, lam_inv = _gls_at(y, X, sigma, g2)
        e = y - X @ beta
        return float(e @ lam_inv @ e) - (n - p)

    if not np.any(sigma):
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        return float(np.sum((y - X @ beta) ** 2) / (n - p))
    if f(0.0) <= 0.0:
        return 0.0
    hi = max(float(np.var(y)), 1e-6) * 10.0
    while f(hi) > 0.0:
        hi *= 10.0
    return float(brentq(f, 0.0, hi, xtol=1e-12))


def fit_gls(
    y: npt.ArrayLike,
    X: npt.ArrayLike,
    sigma: npt.ArrayLike,
    names: Optional[Sequence[str]] = None,
) -> GLSResult:
    """Fit a GLS regression of station log10 quantiles on explanatory variables.

    Parameters
    ----------
    y : array_like
        Station log10 quantiles, length m.
    X : array_like
        Design matrix (m, p), including the intercept column.
    sigma : array_like
        Sampling covariance (m, m), e.g. :func:`sampling_covariance`. All
        zeros gives ordinary least squares.
    names : sequence of str, optional
        Column names of ``X``.

    Returns
    -------
    GLSResult

    Raises
    ------
    ValueError
        Mismatched shapes, fewer sites than parameters plus one, or a design
        matrix that is not of full column rank.
    """
    y_a = np.asarray(y, dtype=float)
    X_a = np.asarray(X, dtype=float)
    S = np.asarray(sigma, dtype=float)
    if X_a.ndim != 2 or X_a.shape[0] != y_a.size or S.shape != (y_a.size, y_a.size):
        raise ValueError(f"shape mismatch: y {y_a.shape}, X {X_a.shape}, sigma {S.shape}")
    m, p = X_a.shape
    if m <= p:
        raise ValueError(f"need more sites than parameters: {m} sites, {p} parameters")
    if np.linalg.matrix_rank(X_a) < p:
        raise ValueError("design matrix is not of full column rank")

    gamma2 = _mev(y_a, X_a, S)
    beta, cov, lam_inv = _gls_at(y_a, X_a, S, gamma2)
    resid = y_a - X_a @ beta
    h = np.einsum("ij,jk,kl,li->i", X_a, cov, X_a.T, lam_inv)
    lam_diag = np.diag(S) + gamma2
    influence = resid**2 / lam_diag * h / np.maximum((1.0 - h) ** 2, 1e-12) / p
    avp = gamma2 + float(np.mean(np.einsum("ij,jk,ik->i", X_a, cov, X_a)))
    sep = 100.0 * float(np.sqrt(np.exp(np.log(10.0) ** 2 * avp) - 1.0))
    if p > 1:
        gamma0 = _mev(y_a, np.ones((m, 1)), S)
        pseudo = 1.0 - gamma2 / gamma0 if gamma0 > 0 else 0.0
    else:
        pseudo = 0.0
    return GLSResult(
        beta=beta,
        cov_beta=cov,
        model_error_variance=gamma2,
        avp=avp,
        sep_percent=sep,
        pseudo_r2=float(pseudo),
        leverage=h,
        influence=influence,
        residuals=resid,
        names=list(names) if names is not None else [f"x{i}" for i in range(p)],
    )
