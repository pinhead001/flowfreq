"""
hydrolib.core - Core data structures and utility functions
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from functools import cached_property, lru_cache
from typing import ClassVar, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import gammaln, ndtri


class SkewMethod(Enum):
    """Skew coefficient computation method."""

    STATION = auto()
    WEIGHTED = auto()
    REGIONAL = auto()


class AnalysisMethod(Enum):
    """Flood frequency analysis method."""

    MOM = auto()  # Method of Moments (traditional)
    EMA = auto()  # Expected Moments Algorithm (B17C)


@dataclass
class PeakRecord:
    """
    Represents a single peak flow record (PeakFQ-style).

    Attributes
    ----------
    year : int
        Water year of the observation
    flow : float, optional
        Peak flow in cfs (None if censored/missing)
    perception_threshold : float, optional
        Minimum flow that would have been recorded
    is_historical : bool
        Whether this is a historical (non-systematic) observation
    source : str, optional
        Data source (e.g., "USGS", "Historical")
    """

    year: int
    flow: Optional[float] = None
    perception_threshold: Optional[float] = None
    is_historical: bool = False
    source: Optional[str] = None


@dataclass
class FlowInterval:
    """
    Represents a flow interval for EMA analysis.

    For systematic peaks: lower = upper = observed value
    For censored data: lower and upper define the interval
    For historical peaks: perception_threshold defines detectability
    """

    lower: float
    upper: float
    year: int
    is_historical: bool = False
    perception_threshold: float = 0.0

    @property
    def is_censored(self) -> bool:
        return self.lower != self.upper

    @property
    def is_systematic(self) -> bool:
        return not self.is_historical

    @classmethod
    def from_peak(cls, flow: float, year: int) -> FlowInterval:
        return cls(lower=flow, upper=flow, year=year)

    @classmethod
    def from_censored(
        cls, lower: float, upper: float, year: int, perception_threshold: float = 0.0
    ) -> FlowInterval:
        return cls(lower=lower, upper=upper, year=year, perception_threshold=perception_threshold)

    @classmethod
    def from_historical(cls, flow: float, year: int, perception_threshold: float) -> FlowInterval:
        return cls(
            lower=flow,
            upper=flow,
            year=year,
            is_historical=True,
            perception_threshold=perception_threshold,
        )


@dataclass
class EMAParameters:
    """Parameters for EMA analysis."""

    systematic_start: int
    systematic_end: int
    historical_start: Optional[int] = None
    historical_end: Optional[int] = None
    historical_threshold: Optional[float] = None
    low_outlier_threshold: Optional[float] = None
    max_iterations: int = 100
    tolerance: float = 1e-6

    @property
    def systematic_years(self) -> int:
        return self.systematic_end - self.systematic_start + 1

    @property
    def historical_years(self) -> int:
        if self.historical_start and self.historical_end:
            return self.historical_end - self.historical_start + 1
        return 0


@dataclass
class FrequencyResults:
    """Results from flood frequency analysis."""

    n_peaks: int
    n_systematic: int
    n_historical: int
    n_censored: int
    n_low_outliers: int
    mean_log: float
    std_log: float
    skew_station: float
    skew_regional: Optional[float]
    skew_weighted: Optional[float]
    skew_used: float
    low_outlier_threshold: float
    mgb_critical_value: float
    method: AnalysisMethod
    quantiles: pd.DataFrame = field(default_factory=pd.DataFrame)
    confidence_limits: pd.DataFrame = field(default_factory=pd.DataFrame)
    ema_iterations: Optional[int] = None
    ema_converged: Optional[bool] = None


@dataclass
class LowFlowResults:
    """
    Results from low-flow frequency analysis.

    Mirrors :class:`FrequencyResults`' shape so the two analyses read alike.
    Fields with no high-flow analogue: `n_zero_years` and `p_zero`, from the
    conditional-probability adjustment for zero-flow years (see
    :mod:`hydrolib.lowflow`). There is no regional-skew map or MGBT-style
    outlier test for the low-flow tail, so those fields have no counterpart
    here.

    Attributes
    ----------
    n_years : int
        Number of years used in the fit (zero-flow years included).
    n_zero_years : int
        Number of those years whose annual minimum n-day mean flow is zero.
    p_zero : float
        ``n_zero_years / n_years``, the fraction of years with no flow.
    n_day : int
        Duration of the annual minimum flow (e.g. 7 for 7Q10/7Q2).
    year_type : str
        Year definition the annual minima were computed over: "climatic",
        "water", or "calendar".
    distribution : str
        Distribution fit to the positive-flow years: "lp3" or "lognormal".
    """

    n_years: int
    n_zero_years: int
    p_zero: float
    n_day: int
    year_type: str
    distribution: str
    mean_log: float
    std_log: float
    skew_station: float
    skew_used: float
    quantiles: pd.DataFrame = field(default_factory=pd.DataFrame)
    confidence_limits: pd.DataFrame = field(default_factory=pd.DataFrame)


# =============================================================================
# UTILITY FUNCTIONS WITH CACHING
# =============================================================================


#: Year definitions shared by hydrolib.lowflow and hydrolib.regime for
#: grouping a daily series into annual periods.
YEAR_TYPES = ("climatic", "water", "calendar")


def assign_year_label(dates: pd.DatetimeIndex, year_type: str) -> np.ndarray:
    """Assign each date the annual period it belongs to.

    Parameters
    ----------
    dates : pd.DatetimeIndex
    year_type : str
        One of:

        - "water": Oct 1 (Y-1) - Sep 30 Y, labeled Y (matches
          :meth:`hydrolib.usgs.USGSgage.download_peak_flow`'s inline
          water-year convention).
        - "climatic": Apr 1 Y - Mar 31 (Y+1), labeled Y. Follows Riggs
          (1972), USGS Techniques of Water-Resources Investigations Book 4
          Chapter B1, chosen for low-flow analysis to avoid splitting a
          single connected low-flow event across a year boundary.
        - "calendar": Jan 1 - Dec 31 Y, labeled Y.

        Both "water" and "climatic" label a year by the calendar year
        containing the majority of its months.

    Returns
    -------
    np.ndarray of int
    """
    if year_type == "calendar":
        return dates.year.to_numpy()
    if year_type == "water":
        return np.where(dates.month >= 10, dates.year + 1, dates.year)
    if year_type == "climatic":
        return np.where(dates.month >= 4, dates.year, dates.year - 1)
    raise ValueError(f"year_type must be one of {YEAR_TYPES}, got {year_type!r}")


#: Maximum absolute skew coefficient considered physically valid for LP3
#: fitting, matching the domain of the Bulletin 17B/17C Appendix 3
#: frequency-factor tables. Skew estimates are clipped to this range
#: because both the Wilson-Hilferty K-factor approximation and the
#: gamma-moment machinery used by EMA become numerically degenerate
#: (shape parameter alpha = 4/skew**2 collapses toward 0) well before
#: this bound, and no real annual-peak record legitimately produces a
#: station skew this extreme.
MAX_ABS_SKEW: float = 3.0


@lru_cache(maxsize=256)
def kfactor(skew: float, aep: float) -> float:
    """
    Calculate K factor for Log-Pearson Type III distribution.
    Uses Wilson-Hilferty approximation. Cached for performance.
    """
    skew = max(-MAX_ABS_SKEW, min(MAX_ABS_SKEW, skew))
    z = ndtri(1 - aep)

    if abs(skew) < 0.001:
        return z

    k = skew / 6
    return (2 / skew) * ((1 + k * z - k * k) ** 3 - 1)


def kfactor_array(skew: float, aep: np.ndarray) -> np.ndarray:
    """Vectorized K factor calculation."""
    return np.array([kfactor(skew, float(p)) for p in aep])


@lru_cache(maxsize=64)
def grubbs_beck_critical_value(n: int, alpha: float = 0.10) -> float:
    """Compute Grubbs-Beck critical value for low outlier detection."""
    exact_values = {
        10: 2.036,
        15: 2.247,
        20: 2.385,
        25: 2.486,
        30: 2.563,
        40: 2.682,
        50: 2.768,
        60: 2.837,
        70: 2.893,
        80: 2.940,
        90: 2.981,
        100: 3.017,
    }

    keys = sorted(exact_values.keys())
    if n <= keys[0]:
        return exact_values[keys[0]]
    if n >= keys[-1]:
        return -0.9043 + 3.345 * np.sqrt(np.log10(n)) - 0.4046 * np.log10(n)

    for i in range(len(keys) - 1):
        if keys[i] <= n < keys[i + 1]:
            t = (n - keys[i]) / (keys[i + 1] - keys[i])
            return exact_values[keys[i]] * (1 - t) + exact_values[keys[i + 1]] * t

    return -0.9043 + 3.345 * np.sqrt(np.log10(n)) - 0.4046 * np.log10(n)


def log_pearson3_cdf(x: float, mean: float, std: float, skew: float) -> float:
    """Compute CDF of Log-Pearson Type III distribution."""
    if x <= 0:
        return 0.0

    log_x = np.log10(x)

    if abs(skew) < 0.001:
        z = (log_x - mean) / std
        return stats.norm.cdf(z)

    alpha = 4 / (skew**2)
    beta = std * abs(skew) / 2
    xi = mean - 2 * std / skew

    if skew > 0:
        y = (log_x - xi) / beta
        if y <= 0:
            return 0.0
        return stats.gamma.cdf(y, alpha)
    else:
        y = (xi - log_x) / beta
        if y <= 0:
            return 1.0
        return 1 - stats.gamma.cdf(y, alpha)


def log_pearson3_ppf(p: float, mean: float, std: float, skew: float) -> float:
    """Compute quantile of Log-Pearson Type III distribution."""
    K = kfactor(skew, 1 - p)
    log_q = mean + K * std
    return 10**log_q


def log_pearson3_pdf(log_q: float, mean: float, std: float, skew: float) -> float:
    """Compute PDF of Log-Pearson Type III distribution in log space."""
    if abs(skew) < 0.001:
        z = (log_q - mean) / std
        return np.exp(-(z**2) / 2) / (std * np.sqrt(2 * np.pi))
    else:
        alpha = 4 / skew**2
        beta = std * abs(skew) / 2
        xi = mean - 2 * std / skew
        if skew > 0:
            y = (log_q - xi) / beta
            if y <= 0:
                return 0.0
            return y ** (alpha - 1) * np.exp(-y) / (beta * np.exp(gammaln(alpha)))
        else:
            y = (xi - log_q) / beta
            if y <= 0:
                return 0.0
            return y ** (alpha - 1) * np.exp(-y) / (beta * np.exp(gammaln(alpha)))


# =============================================================================
# PEAKFQ-STYLE LP3 FUNCTIONS
# =============================================================================


def lp3_frequency_factor_peakfq(p: float, skew: float) -> float:
    """
    Calculate LP3 frequency factor using PeakFQ methodology.

    Uses the exact gamma-distribution transformation (the same inverse-CDF
    approach as the peakfqr/PeakFQ Fortran reference, per ``qP3sub``), not the
    Wilson-Hilferty approximation used by :func:`kfactor` elsewhere in this
    module. Falls back to the normal distribution for zero skew, where the
    two are identical.

    Parameters
    ----------
    p : float
        Non-exceedance probability (0 < p < 1)
    skew : float
        Skew coefficient

    Returns
    -------
    float
        Frequency factor K

    Notes
    -----
    For negative skew the gamma distribution is reflected: the Pearson III
    variate's lower tail (small ``log_Q``) corresponds to the gamma
    distribution's *upper* tail, so the quantile must be taken at ``1 - p``,
    not ``p``. This mirrors the reflection already used by
    :func:`log_pearson3_cdf` for the CDF direction. Getting this backwards
    makes K decrease as p increases, which silently inverts every quantile
    for any negative-skew fit -- the national B17C regional skew (-0.302) is
    negative, so this is not an edge case.
    """
    if abs(skew) < 1e-8:
        return stats.norm.ppf(p)
    alpha = 4 / (skew**2)
    beta = skew / 2
    q = p if skew > 0 else 1 - p
    return beta * (stats.gamma.ppf(q, a=alpha) - alpha)


def lp3_quantile_peakfq(mu: float, sigma: float, skew: float, T: float) -> float:
    """
    Calculate LP3 quantile for a given return period using PeakFQ methodology.

    Parameters
    ----------
    mu : float
        Mean of log10(flow)
    sigma : float
        Standard deviation of log10(flow)
    skew : float
        Skew coefficient
    T : float
        Return period in years

    Returns
    -------
    float
        Flow quantile in cfs
    """
    p = 1 - 1 / T
    K = lp3_frequency_factor_peakfq(p, skew)
    return 10 ** (mu + K * sigma)


def compute_ci_lp3(
    mu: float, sigma: float, skew: float, n: int, T: float, alpha: float = 0.05
) -> Tuple[float, float, float]:
    """
    Compute confidence intervals for LP3 quantile estimates.

    Parameters
    ----------
    mu : float
        Mean of log10(flow)
    sigma : float
        Standard deviation of log10(flow)
    skew : float
        Skew coefficient
    n : int
        Sample size
    T : float
        Return period in years
    alpha : float
        Significance level (default 0.05 for 95% CI)

    Returns
    -------
    tuple
        (estimate, lower_bound, upper_bound) in cfs
    """
    p = 1 - 1 / T
    K = lp3_frequency_factor_peakfq(p, skew)
    x = mu + K * sigma
    tcrit = stats.t.ppf(1 - alpha / 2, df=n - 1)
    var = (sigma**2 / n) * (1 + (K**2) / 2)
    se = np.sqrt(var)
    return 10**x, 10 ** (x - tcrit * se), 10 ** (x + tcrit * se)
