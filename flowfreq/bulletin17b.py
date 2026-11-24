"""Bulletin 17B (1982) flood frequency analysis, superseded by Bulletin 17C.

Bulletin 17B, "Guidelines for Determining Flood Flow Frequency", is the method of
moments procedure that Bulletin 17C (2019) superseded. It is kept here to reproduce
older studies and to compare them with the Bulletin 17C fit. It is superseded, and
nothing in flowfreq runs it unless the caller asks for it.

This module holds the bulletin's equations, one function each, with the equation
numbers in the docstrings, and :class:`Bulletin17B`, which applies them in the order of
the bulletin's flow charts. All logs are base 10. Page labels and the places where the
bulletin is ambiguous are in ``tests/fixtures/b17b/procedure.md``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from decimal import ROUND_HALF_UP, Decimal
from typing import Any, Mapping, NamedTuple, Optional, Sequence

import numpy as np
import pandas as pd
from scipy.special import ndtri
from scipy.stats import t as student_t

from flowfreq.aep import same_aep
from flowfreq.core import kfactor
from flowfreq.peak_codes import (
    GREATER_THAN_CODES,
    LESS_THAN_CODES,
    PeakInterval,
    PeakTreatment,
    peak_interval,
)

logger = logging.getLogger(__name__)

__all__ = [
    "Bulletin17B",
    "Bulletin17BError",
    "Bulletin17BResults",
    "ConditionalAdjustment",
    "HistoricAdjustment",
    "HistoricMoments",
    "KN_TABLE",
    "STANDARD_AEP",
    "bulletin17b_kwargs",
    "confidence_deviates",
    "expected_probability",
    "historic_moments",
    "log_moments",
    "outlier_k",
    "removed_peaks",
    "station_skew_mse",
    "synthetic_statistics",
    "weighted_skew",
]


class Bulletin17BError(ValueError):
    """Bulletin 17B cannot be applied to this input, or the bulletin does not cover it."""


#: Appendix 4: the outlier test K_N, by sample size N = 10 to 149, as printed.
KN_TABLE: Mapping[int, float] = {
    10: 2.036,
    11: 2.088,
    12: 2.134,
    13: 2.175,
    14: 2.213,
    15: 2.247,
    16: 2.279,
    17: 2.309,
    18: 2.335,
    19: 2.361,
    20: 2.385,
    21: 2.408,
    22: 2.429,
    23: 2.448,
    24: 2.467,
    25: 2.486,
    26: 2.502,
    27: 2.519,
    28: 2.534,
    29: 2.549,
    30: 2.563,
    31: 2.577,
    32: 2.591,
    33: 2.604,
    34: 2.616,
    35: 2.628,
    36: 2.639,
    37: 2.650,
    38: 2.661,
    39: 2.671,
    40: 2.682,
    41: 2.692,
    42: 2.700,
    43: 2.710,
    44: 2.719,
    45: 2.727,
    46: 2.736,
    47: 2.744,
    48: 2.753,
    49: 2.760,
    50: 2.768,
    51: 2.775,
    52: 2.783,
    53: 2.790,
    54: 2.798,
    55: 2.804,
    56: 2.811,
    57: 2.818,
    58: 2.824,
    59: 2.831,
    60: 2.837,
    61: 2.842,
    62: 2.849,
    63: 2.854,
    64: 2.860,
    65: 2.866,
    66: 2.871,
    67: 2.877,
    68: 2.883,
    69: 2.888,
    70: 2.893,
    71: 2.897,
    72: 2.903,
    73: 2.908,
    74: 2.912,
    75: 2.917,
    76: 2.922,
    77: 2.927,
    78: 2.931,
    79: 2.935,
    80: 2.940,
    81: 2.945,
    82: 2.949,
    83: 2.953,
    84: 2.957,
    85: 2.961,
    86: 2.966,
    87: 2.970,
    88: 2.973,
    89: 2.977,
    90: 2.981,
    91: 2.984,
    92: 2.989,
    93: 2.993,
    94: 2.996,
    95: 3.000,
    96: 3.003,
    97: 3.006,
    98: 3.011,
    99: 3.014,
    100: 3.017,
    101: 3.021,
    102: 3.024,
    103: 3.027,
    104: 3.030,
    105: 3.033,
    106: 3.037,
    107: 3.040,
    108: 3.043,
    109: 3.046,
    110: 3.049,
    111: 3.052,
    112: 3.055,
    113: 3.058,
    114: 3.061,
    115: 3.064,
    116: 3.067,
    117: 3.070,
    118: 3.073,
    119: 3.075,
    120: 3.078,
    121: 3.081,
    122: 3.083,
    123: 3.086,
    124: 3.089,
    125: 3.092,
    126: 3.095,
    127: 3.097,
    128: 3.100,
    129: 3.102,
    130: 3.104,
    131: 3.107,
    132: 3.109,
    133: 3.112,
    134: 3.114,
    135: 3.116,
    136: 3.119,
    137: 3.122,
    138: 3.124,
    139: 3.126,
    140: 3.129,
    141: 3.131,
    142: 3.133,
    143: 3.135,
    144: 3.138,
    145: 3.140,
    146: 3.142,
    147: 3.144,
    148: 3.146,
    149: 3.148,
}


def outlier_k(n: int) -> float:
    """Return the Appendix 4 outlier-test value K_N for sample size ``n``.

    Raises
    ------
    Bulletin17BError
        If ``n`` is outside 10 to 149, the range the bulletin prints.
    """
    try:
        return KN_TABLE[int(n)]
    except KeyError:
        raise Bulletin17BError(
            f"Appendix 4 prints K_N for N = 10 to 149 only; got N = {n}"
        ) from None


def log_moments(x: np.ndarray) -> tuple[float, float, float]:
    """Mean, standard deviation and skew of log flows (eqs. 2, 3a and 4a).

    Parameters
    ----------
    x : numpy.ndarray
        Base-10 logs of the peaks.

    Returns
    -------
    tuple of float
        ``(mean, S, G)``, with S using N - 1 and G the bias-adjusted skew of eq. 4a.

    Raises
    ------
    Bulletin17BError
        If there are fewer than three values or the standard deviation is zero.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 3:
        raise Bulletin17BError(f"moments need at least 3 values; got {n}")
    if not np.all(np.isfinite(x)):
        raise Bulletin17BError("the log flows contain a value that is not finite (a zero flow?)")
    if np.ptp(x) == 0.0:
        raise Bulletin17BError("the log flows are all equal, so S is zero and G undefined")
    mean = float(x.mean())
    s = float(x.std(ddof=1))
    g = float(n * np.sum((x - mean) ** 3) / ((n - 1) * (n - 2) * s**3))
    return mean, s, g


def station_skew_mse(skew: float, n: float) -> float:
    """Mean square error of the station skew (eq. 6, as printed).

    ``MSE = 10 ** (A - B * log10(N / 10))``, with A = -0.33 + 0.08|G| for |G| <= 0.90
    and -0.52 + 0.30|G| above, and B = 0.94 - 0.26|G| for |G| <= 1.50 and 0.55 above.
    Table 1 of the bulletin is computed with the upper branch of A at |G| = 0.9; the
    equation is used as printed here.

    Parameters
    ----------
    skew : float
        Station skew, taken as the population skew.
    n : float
        Record length in years (H for a historically adjusted skew).
    """
    g = abs(float(skew))
    a = -0.33 + 0.08 * g if g <= 0.90 else -0.52 + 0.30 * g
    b = 0.94 - 0.26 * g if g <= 1.50 else 0.55
    return float(10.0 ** (a - b * np.log10(float(n) / 10.0)))


def weighted_skew(
    station: float, station_mse: float, generalized: float, generalized_mse: float
) -> float:
    """Weight the station and generalized skews by inverse MSE (eq. 5)."""
    total = generalized_mse + station_mse
    if total <= 0.0:
        raise Bulletin17BError("the skew mean square errors must sum to a positive number")
    return float((generalized_mse * station + station_mse * generalized) / total)


class HistoricMoments(NamedTuple):
    """Historically adjusted moments of Appendix 6."""

    w: float
    mean: float
    std: float
    skew: float


def historic_moments(
    x_systematic: np.ndarray, x_historic: np.ndarray, h: int, l: int
) -> HistoricMoments:
    """Historically weighted mean, standard deviation and skew (eqs. 6-1 to 6-4a).

    Parameters
    ----------
    x_systematic : numpy.ndarray
        Logs of the systematic peaks (N values).
    x_historic : numpy.ndarray
        Logs of the historic peaks (Z values).
    h : int
        Length of the historic period in years.
    l : int
        Number of low values (zeros, below-base years, low outliers) left out of
        ``x_systematic`` but counted in its record.

    Raises
    ------
    Bulletin17BError
        If the weights leave fewer than three effective years or the spread is zero.
    """
    xs = np.asarray(x_systematic, dtype=float)
    xz = np.asarray(x_historic, dtype=float)
    z = xz.size
    n = xs.size
    if n + l == 0:
        raise Bulletin17BError("the systematic record is empty")
    if not (np.all(np.isfinite(xs)) and np.all(np.isfinite(xz))):
        raise Bulletin17BError("the log flows contain a value that is not finite (a zero flow?)")
    w = (h - z) / (n + l)  # eq. 6-1
    if w < 1.0:
        raise Bulletin17BError(
            f"historic weight W = {w:.4f} is below 1: H = {h} is shorter than the "
            f"{n + l + z} years of record"
        )
    m = h - w * l
    if m - 2 <= 0:
        raise Bulletin17BError(f"historic period H = {h} is too short for the record")
    if np.ptp(np.concatenate([xs, xz])) == 0.0:
        raise Bulletin17BError("the log flows are all equal, so S is zero and G undefined")
    mean = float((w * xs.sum() + xz.sum()) / m)  # eq. 6-2a
    d_s = xs - mean
    d_z = xz - mean
    var = float((w * np.sum(d_s**2) + np.sum(d_z**2)) / (m - 1))  # eq. 6-3a
    if not var > 0.0:
        raise Bulletin17BError("the historically adjusted standard deviation is not positive")
    std = var**0.5
    skew = float(
        m / ((m - 1) * (m - 2)) * (w * np.sum(d_s**3) + np.sum(d_z**3)) / std**3
    )  # eq. 6-4a
    return HistoricMoments(float(w), mean, float(std), skew)


def synthetic_statistics(
    q01: float, q10: float, q50: float, *, k_skew: Optional[float] = None
) -> tuple[float, float, float]:
    """Synthetic moments of a conditional curve (eqs. 5-3, 5-4 and 5-5).

    Parameters
    ----------
    q01, q10, q50 : float
        Flows at exceedance probabilities 0.01, 0.10 and 0.50 read from the
        conditional curve, in cfs.
    k_skew : float, optional
        Skew at which K.01 and K.50 are taken. By default the synthetic skew
        itself. A hand computation passes the skew rounded to a tenth.

    Returns
    -------
    tuple of float
        ``(mean_s, std_s, skew_s)`` in log units.

    Raises
    ------
    Bulletin17BError
        If a flow is not positive, the flows are not in decreasing order, or the
        synthetic skew is outside -2.0 to 2.5, where eq. 5-3 does not apply.
    """
    if not (q01 > q10 > q50 > 0):
        raise Bulletin17BError(f"need Q.01 > Q.10 > Q.50 > 0; got {q01}, {q10}, {q50}")
    l01, l10, l50 = np.log10(q01), np.log10(q10), np.log10(q50)
    skew_s = float(-2.50 + 3.12 * (l01 - l10) / (l10 - l50))  # eq. 5-3
    if not -2.0 <= skew_s <= 2.5:
        raise Bulletin17BError(
            f"synthetic skew {skew_s:.3f} is outside -2.0 to 2.5, where eq. 5-3 applies"
        )
    kg = skew_s if k_skew is None else float(k_skew)
    k01 = kfactor(kg, 0.01)
    k50 = kfactor(kg, 0.50)
    std_s = float((l01 - l50) / (k01 - k50))  # eq. 5-4
    mean_s = float(l50 - k50 * std_s)  # eq. 5-5
    return mean_s, std_s, skew_s


def confidence_deviates(k: np.ndarray, n: int, z: float) -> tuple[np.ndarray, np.ndarray]:
    """Upper and lower confidence-limit deviates (eqs. 9-4a, 9-4b, 9-5 and 9-6).

    Parameters
    ----------
    k : numpy.ndarray
        Pearson III deviates K at the weighted skew.
    n : int
        Systematic record length.
    z : float
        Standard normal deviate of the one-sided confidence level, for example
        ``ndtri(0.95)`` for a two-sided 90 percent interval.

    Returns
    -------
    tuple of numpy.ndarray
        ``(k_upper, k_lower)``.
    """
    k = np.asarray(k, dtype=float)
    if n < 2:
        raise Bulletin17BError(f"confidence limits need N >= 2; got N = {n}")
    a = 1.0 - z**2 / (2.0 * (n - 1))  # eq. 9-5
    if a <= 0.0:
        raise Bulletin17BError(
            f"eq. 9-5 gives a = {a:.4f} <= 0 for N = {n} and z = {z:.4f}: "
            "the record is too short for this confidence level"
        )
    b = k**2 - z**2 / n  # eq. 9-6
    disc = k**2 - a * b
    if not np.all(disc >= 0.0):
        raise Bulletin17BError("eqs. 9-4a and 9-4b have a negative discriminant")
    root = np.sqrt(disc)
    return (k + root) / a, (k - root) / a  # eqs. 9-4a, 9-4b


def expected_probability(aep: np.ndarray, n: int) -> np.ndarray:
    """Expected probability P_N of a quantile with exceedance probability ``aep`` (eq. 11-1).

    ``P_N = Prob[t(N-1) > K_n * sqrt(N / (N + 1))]``, K_n the standard normal
    deviate of ``aep``.
    """
    aep = np.asarray(aep, dtype=float)
    return np.asarray(student_t.sf(ndtri(1.0 - aep) * np.sqrt(n / (n + 1.0)), df=n - 1))


#: Fewest systematic peaks the moments may rest on, once peaks have moved to the historic
#: ones and low outliers have been dropped. Appendix 4 starts at 10 as well.
MIN_SYSTEMATIC = 10

#: Largest share of the systematic record that may be truncated (p. 5-1).
MAX_TRUNCATED = 0.25

#: The bulletin's frequency-table rows (Table 12-3), as annual exceedance probabilities.
STANDARD_AEP: tuple[float, ...] = (0.99, 0.90, 0.50, 0.10, 0.05, 0.02, 0.01, 0.005, 0.002)


@dataclass
class HistoricAdjustment:
    """Historic adjustment of Appendix 6 (filled in when the adjustment is made).

    Attributes
    ----------
    h : int
        Years in the historic period.
    z : int
        Historic peaks: those given, the high outliers, and the systematic peaks above
        the smallest given historic peak.
    l : int
        Low values left out of the moments but counted in the systematic record.
    w : float
        Weight of a systematic peak, (H - Z) / (N + L) (eq. 6-1).
    mean, std, skew : float
        Adjusted moments of the logs (eqs. 6-2a to 6-4a).
    historic_peaks : list of (year, flow)
        The Z peaks, by year.
    """

    h: int
    z: int
    l: int
    w: float
    mean: float
    std: float
    skew: float
    historic_peaks: list[tuple[int, float]]


@dataclass
class ConditionalAdjustment:
    """Conditional probability adjustment of Appendix 5 (filled in when it is made).

    Attributes
    ----------
    p_tilde : float
        Probability that an annual peak is above the truncation level: N / n (eq. 5-1a),
        or (H - W L) / H after a historic adjustment (eq. 5-1b).
    n_truncated : int
        Truncated years: zeros, years below the gage base and low outliers (L).
    mean, std, skew : float
        Moments of the conditional curve, from the peaks kept (historically adjusted when
        that adjustment was made). Never rounded; hand-computation mode evaluates the
        curve at the skew rounded to a tenth.
    q01, q10, q50 : float
        Adjusted curve in cfs at exceedance probabilities 0.01, 0.10 and 0.50, which is
        the conditional curve at P / ``p_tilde`` (eq. 5-2).
    mean_s, std_s, skew_s : float
        Synthetic statistics (eqs. 5-3 to 5-5). The skew is not yet weighted.
    """

    p_tilde: float
    n_truncated: int
    mean: float
    std: float
    skew: float
    q01: float
    q10: float
    q50: float
    mean_s: float
    std_s: float
    skew_s: float


@dataclass
class Bulletin17BResults:
    """Result of a Bulletin 17B fit.

    Attributes
    ----------
    n_years : int
        Years in the systematic record, including zero, below-base and low-outlier years.
    n_systematic : int
        Peaks that enter the moments (the bulletin's N).
    n_zeros, n_below_base : int
        Years left out of the moments as zeros or below the gage base. A below-base year
        whose recorded flow is zero counts as below base only.
    outlier_order : str
        ``"high first"``, ``"low first"`` or ``"both on unadjusted statistics"``.
    high_outlier_threshold, low_outlier_threshold : float
        Outlier thresholds in cfs, from the statistics each test used (eqs. 7, 8a, 8b),
        or the low-outlier threshold the caller gave.
    high_outliers, low_outliers : list of (year, flow)
        Peaks above the high threshold and strictly below the low one.
    historic, conditional : optional
        Appendix 6 and Appendix 5 adjustments, ``None`` when not made.
    mean_log, std_log : float
        Mean and standard deviation of the log10 flows of the final curve.
    skew_station : float
        The station skew that was weighted (historically adjusted or synthetic when
        those adjustments were made).
    skew_station_mse : float
        Its mean square error from eq. 6.
    skew_generalized, skew_generalized_mse, skew_weighted : float, optional
        Generalized skew, its MSE and the eq. 5 weighted skew; ``None`` without a
        generalized skew.
    skew_used : float
        The skew the quantiles use. In hand-computation mode a weighted skew is rounded
        to a tenth; a station skew used alone is not.
    quantiles : pandas.DataFrame
        Columns ``aep, flow_cfs, log_flow, K_factor, lower, upper,
        expected_probability``. ``lower`` and ``upper`` are the confidence limits in cfs.
        After a conditional adjustment, rows with an AEP above ``conditional.p_tilde``
        extrapolate the synthetic curve below the truncation level; ``steps`` names them.
    plotting_positions : pandas.DataFrame
        Columns ``year, flow_cfs, weighted_order, exceedance``: Weibull positions, and
        Appendix 6's weighted ones (eqs. 6-6 to 6-8, a = 0) after a historic adjustment.
        Every year of the systematic record has a row, so truncated years take ranks too:
        low outliers at their flow, then below-base years (``flow_cfs`` NaN), then zeros.
    steps : list of str
        One plain sentence per decision of the bulletin flow chart.
    """

    n_years: int
    n_systematic: int
    n_zeros: int
    n_below_base: int
    outlier_order: str
    high_outlier_threshold: float
    low_outlier_threshold: float
    high_outliers: list[tuple[int, float]]
    low_outliers: list[tuple[int, float]]
    historic: Optional[HistoricAdjustment]
    conditional: Optional[ConditionalAdjustment]
    mean_log: float
    std_log: float
    skew_station: float
    skew_station_mse: float
    skew_generalized: Optional[float]
    skew_generalized_mse: Optional[float]
    skew_weighted: Optional[float]
    skew_used: float
    quantiles: pd.DataFrame
    plotting_positions: pd.DataFrame
    steps: list[str] = field(default_factory=list)


def _finite_number(value: object, name: str) -> float:
    try:
        x = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        raise Bulletin17BError(f"{name} must be a number; got {value!r}") from None
    if not np.isfinite(x):
        raise Bulletin17BError(f"{name} must be finite; got {value!r}")
    return x


def _integer_year(value: object, name: str) -> int:
    if isinstance(value, (str, bytes)):
        raise Bulletin17BError(f"{name}s must be whole numbers, not text; got {value!r}")
    year = _finite_number(value, name)
    if isinstance(value, (bool, np.bool_)) or year != int(year):
        raise Bulletin17BError(f"{name}s must be whole numbers (integers); got {value!r}")
    return int(year)


def _require_unique_years(years: Sequence[int], name: str) -> None:
    if len(set(years)) != len(years):
        dup = sorted({y for y in years if years.count(y) > 1})
        raise Bulletin17BError(f"{name} appears more than once: {dup}")


def _nearest_tenth(x: float) -> float:
    """Round to the nearest tenth with ties away from zero, on the decimal value.

    Ten decimals are kept first, so that a tie carrying floating-point noise is still
    a tie. Never returns a negative zero.
    """
    return float(Decimal(f"{x:.10f}").quantize(Decimal("0.1"), ROUND_HALF_UP)) + 0.0


def _outlier_order(skew: float) -> str:
    """Which outlier test comes first for this station skew (p. 17, p. 12-3)."""
    if skew > 0.4:
        return "high first"
    if skew < -0.4:
        return "low first"
    return "both on unadjusted statistics"


class _Screened(NamedTuple):
    """The outlier tests and the moments they and the historic adjustment lead to."""

    order: str
    high_threshold_log: float
    low_threshold: float
    high_outliers: list[tuple[int, float]]
    low_outliers: list[tuple[int, float]]
    historic: Optional[HistoricAdjustment]
    n_systematic: int
    mean: float
    std: float
    skew: float


def _year_list(years: Sequence[int]) -> str:
    return ", ".join(str(y) for y in years)


class Bulletin17B:
    """Bulletin 17B (1982) flood frequency analysis, superseded by Bulletin 17C.

    The steps are the moments of the log peaks, the outlier tests, the historic adjustment
    of Appendix 6, the conditional probability adjustment of Appendix 5, skew weighting,
    the log-Pearson III curve, confidence limits and expected probability.

    The outlier tests follow the flow chart on p. 12-3. With a station skew above +0.4
    the high test and the historic adjustment come first, and the low test after an
    adjustment uses the adjusted statistics (eq. 8b) with K read at the number of
    systematic peaks tested. Below -0.4 the low test comes first. In between, both use
    the unadjusted statistics. Low outliers are then dropped and the historic
    adjustment is made, with L counting them.

    Zero flows, years below the gage base and low outliers are truncated: they stay in the
    record length and in L but not in the moments. A record with any of them gets the
    conditional probability adjustment, and its curve comes from the synthetic statistics
    (eqs. 5-3 to 5-5). A record with more than 25 percent of its years truncated is
    refused.

    Parameters
    ----------
    peak_flows : sequence of float
        Annual peak flows in cfs. Zero is a zero flood year. Negative flows are refused.
    water_years : sequence of int
        Water year of each peak. A year may appear only once.
    generalized_skew : float, optional
        Generalized (regional) skew. Without it the station skew is used.
    generalized_skew_mse : float, default 0.302
        Mean square error of the generalized skew (Plate I value).
    historic_peaks : sequence of (year, flow), optional
        Known peaks in cfs outside the systematic record. They need a
        ``historic_period``. Systematic peaks above the smallest of them are treated as
        historic peaks too. A systematic peak equal to it stays systematic.
    historic_period : (first_year, last_year), optional
        The period, of H = last - first + 1 years, in which the historic peaks and any
        high outliers are known to be the largest. It has to cover the systematic
        record. Without it, high outliers are kept in the systematic record. At least
        10 systematic peaks have to be left once peaks have moved to the historic ones
        and truncated years are out.
    below_base_years : sequence of int, optional
        Water years of the systematic record whose peak was below the gage base. A year
        may have a flow in ``peak_flows`` or not; a flow given for it is not used. It may
        not be a historic peak as well.
    low_outlier_threshold : float, optional
        Threshold in cfs that replaces the one of eq. 8a or 8b. Peaks strictly below it
        are truncated as low outliers; a peak equal to it is kept. It is applied where the
        flow chart has the low-outlier test, so a high test that the bulletin makes
        without the low outliers is made without these.
    confidence : float, default 0.90
        Two-sided confidence level of the limits.
    aeps : sequence of float, optional
        Annual exceedance probabilities to report; ``STANDARD_AEP`` by default. Two
        that are one flood written differently (0.667 and 0.67) are refused.
    hand_computation : bool, default False
        Round the weighted skew to the nearest tenth (ties away from zero) before
        reading K, as the worked examples in the bulletin do. The same rounding applies to
        the station skew that K is read at for the conditional curve of Appendix 5, and to
        the synthetic skew that K is read at in eq. 5-4. A station skew used alone,
        without a generalized skew, stays unrounded. Exact mode never rounds.

    Raises
    ------
    Bulletin17BError
        On invalid input, or a record the bulletin or this class does not cover.
    """

    def __init__(
        self,
        peak_flows: Sequence[float],
        water_years: Sequence[int],
        *,
        generalized_skew: Optional[float] = None,
        generalized_skew_mse: float = 0.302,
        historic_peaks: Optional[Sequence[tuple[int, float]]] = None,
        historic_period: Optional[tuple[int, int]] = None,
        below_base_years: Optional[Sequence[int]] = None,
        low_outlier_threshold: Optional[float] = None,
        confidence: float = 0.90,
        aeps: Optional[Sequence[float]] = None,
        hand_computation: bool = False,
    ) -> None:
        flows = [_finite_number(q, "a peak flow") for q in peak_flows]
        if isinstance(water_years, (str, bytes)):
            raise Bulletin17BError(f"water_years must be a sequence of years; got {water_years!r}")
        raw_years = list(water_years)
        if len(flows) != len(raw_years):
            raise Bulletin17BError(
                f"{len(flows)} peak flows but {len(raw_years)} water years: lengths must match"
            )
        if any(q < 0 for q in flows):
            raise Bulletin17BError("peak flows must not be negative")
        years = [_integer_year(y, "water year") for y in raw_years]
        _require_unique_years(years, "a water year")
        self._below_base_years = self._checked_below_base(below_base_years)
        below = set(self._below_base_years)
        # Below-base years that came with a flow; the flow is not used.
        self._flows_ignored = sorted(below & set(years))
        measured = sorted((y, q) for y, q in zip(years, flows) if y not in below)
        # The peaks above zero and above the gage base, which the statistics are made of.
        self._years = [y for y, q in measured if q > 0.0]
        self._flows = [q for _, q in measured if q > 0.0]
        self._zero_years = [y for y, q in measured if q == 0.0]
        # Every year of the systematic record, whether its peak was measured or not.
        self._all_years = sorted(self._years + self._zero_years + self._below_base_years)

        self._gen_skew: Optional[float] = (
            None
            if generalized_skew is None
            else _finite_number(generalized_skew, "generalized_skew")
        )
        mse = _finite_number(generalized_skew_mse, "generalized_skew_mse")
        if mse <= 0.0:
            raise Bulletin17BError(f"generalized_skew_mse must be positive; got {mse}")
        self._gen_mse = mse
        conf = _finite_number(confidence, "confidence")
        if not 0.0 < conf < 1.0:
            raise Bulletin17BError(f"confidence must be between 0 and 1; got {conf}")
        self._confidence = conf
        aep_list = [_finite_number(a, "an AEP") for a in (STANDARD_AEP if aeps is None else aeps)]
        if not aep_list or any(not 0.0 < a < 1.0 for a in aep_list):
            raise Bulletin17BError("AEPs must be a non-empty list of values between 0 and 1")
        for i, a in enumerate(aep_list):
            for b in aep_list[:i]:
                if same_aep(a, b):
                    raise Bulletin17BError(f"{b:g} and {a:g} are the same AEP given twice")
        self._aeps = aep_list
        self._hand = bool(hand_computation)
        self._historic_period = self._checked_period(historic_period)
        self._historic_peaks = self._checked_historic_peaks(historic_peaks)
        self._low_threshold: Optional[float] = None
        if low_outlier_threshold is not None:
            threshold = _finite_number(low_outlier_threshold, "low_outlier_threshold")
            if threshold <= 0.0:
                raise Bulletin17BError(f"low_outlier_threshold must be positive; got {threshold:g}")
            self._low_threshold = threshold

    @staticmethod
    def _checked_below_base(years: Optional[Sequence[int]]) -> list[int]:
        """The below-base years as sorted, distinct integers."""
        if years is None:
            return []
        try:
            if isinstance(years, (str, bytes)):
                raise TypeError
            raw = list(years)
        except TypeError:
            raise Bulletin17BError(
                f"below_base_years must be a sequence of years; got {years!r}"
            ) from None
        checked = [_integer_year(y, "below-base year") for y in raw]
        _require_unique_years(checked, "a below-base year")
        return sorted(checked)

    def _checked_period(self, period: Optional[tuple[int, int]]) -> Optional[tuple[int, int]]:
        """The historic period as two years that cover the systematic record."""
        if period is None:
            return None
        try:
            first_raw, last_raw = period
        except (TypeError, ValueError):
            raise Bulletin17BError(
                f"historic_period must be (first_year, last_year); got {period!r}"
            ) from None
        first = _integer_year(first_raw, "historic_period year")
        last = _integer_year(last_raw, "historic_period year")
        if last < first:
            raise Bulletin17BError(f"historic_period ends before it begins: {first} to {last}")
        record = self._all_years
        if record and not (first <= record[0] and record[-1] <= last):
            raise Bulletin17BError(
                f"historic_period {first} to {last} does not cover the systematic record, "
                f"{record[0]} to {record[-1]}"
            )
        return first, last

    def _checked_historic_peaks(
        self, peaks: Optional[Sequence[tuple[int, float]]]
    ) -> list[tuple[int, float]]:
        """The historic peaks by year, each inside the period and outside the record."""
        checked: list[tuple[int, float]] = []
        for pair in [] if peaks is None else peaks:
            try:
                year_raw, flow_raw = pair
            except (TypeError, ValueError):
                raise Bulletin17BError(
                    f"historic_peaks must be (year, flow) pairs; got {pair!r}"
                ) from None
            year = _integer_year(year_raw, "historic peak year")
            flow = _finite_number(flow_raw, "a historic peak flow")
            if flow <= 0.0:
                raise Bulletin17BError(f"historic peak flows must be positive; got {flow:g}")
            checked.append((year, flow))
        if not checked:
            return []
        years = [y for y, _ in checked]
        _require_unique_years(years, "a historic peak year")
        if self._historic_period is None:
            raise Bulletin17BError("historic_peaks need a historic_period, which sets H")
        first, last = self._historic_period
        outside = sorted(y for y in years if not first <= y <= last)
        if outside:
            raise Bulletin17BError(
                f"historic peak(s) in {outside} are outside the historic period {first} to {last}"
            )
        both = sorted(set(years) & set(self._all_years))
        if both:
            raise Bulletin17BError(
                f"historic peak(s) in {both} are in the systematic record as well"
            )
        return sorted(checked)

    def run_analysis(self) -> Bulletin17BResults:
        """Fit the record and return the frequency table, limits and steps.

        Raises
        ------
        Bulletin17BError
            If the record is constant, has fewer than 10 or more than 149 peaks above zero
            and the gage base, has more than 25 percent of its years truncated, or gives a
            synthetic skew outside the range of eq. 5-3.
        """
        steps: list[str] = []
        n = len(self._all_years)
        s = self._screen_and_adjust(steps)
        mean, std, skew = s.mean, s.std, s.skew
        conditional = self._adjust_for_truncation(s, steps)
        if conditional is not None:
            mean, std, skew = conditional.mean_s, conditional.std_s, conditional.skew_s
        skew_mse = station_skew_mse(skew, n if s.historic is None else s.historic.h)
        label = "station skew" if conditional is None else "synthetic skew"
        skew_w, skew_used = self._weight_skew(skew, skew_mse, label, steps)
        quantiles = self._frequency_curve(mean, std, skew_used, n, n, steps)
        if conditional is not None:
            self._note_extrapolation(conditional.p_tilde, steps)
        return Bulletin17BResults(
            n_years=n,
            n_systematic=s.n_systematic,
            n_zeros=len(self._zero_years),
            n_below_base=len(self._below_base_years),
            outlier_order=s.order,
            high_outlier_threshold=float(10.0**s.high_threshold_log),
            low_outlier_threshold=s.low_threshold,
            high_outliers=s.high_outliers,
            low_outliers=s.low_outliers,
            historic=s.historic,
            conditional=conditional,
            mean_log=mean,
            std_log=std,
            skew_station=skew,
            skew_station_mse=skew_mse,
            skew_generalized=self._gen_skew,
            skew_generalized_mse=None if self._gen_skew is None else self._gen_mse,
            skew_weighted=skew_w,
            skew_used=skew_used,
            quantiles=quantiles,
            plotting_positions=self._plotting_positions(s.historic),
            steps=steps,
        )

    @property
    def _n_unmeasured(self) -> int:
        """Zero and below-base years, the part of L known before any outlier test."""
        return len(self._zero_years) + len(self._below_base_years)

    def _truncated_parts(self, n_low: int) -> list[str]:
        """What is truncated, in words, for messages."""
        parts = []
        if self._zero_years:
            parts.append(f"{len(self._zero_years)} zero year(s)")
        if self._below_base_years:
            parts.append(f"{len(self._below_base_years)} year(s) below the gage base")
        if n_low:
            parts.append(f"{n_low} low outlier(s)")
        return parts

    def _require_mostly_kept(self, n_low: int) -> None:
        """Refuse when more than 25 percent of the record is truncated (p. 5-1)."""
        n = len(self._all_years)
        truncated = self._n_unmeasured + n_low
        if truncated > MAX_TRUNCATED * n:
            raise Bulletin17BError(
                f"{truncated} of the {n} years are truncated "
                f"({', '.join(self._truncated_parts(n_low))}). Bulletin 17B allows the "
                "conditional probability adjustment only when not over 25 percent of the "
                "record is truncated."
            )

    def _adjust_for_truncation(
        self, s: _Screened, steps: list[str]
    ) -> Optional[ConditionalAdjustment]:
        """Make the Appendix 5 adjustment, or return None when nothing is truncated."""
        n_low = len(s.low_outliers)
        truncated = self._n_unmeasured + n_low
        if truncated == 0:
            return None
        n = len(self._all_years)
        if s.historic is None:
            p_tilde = s.n_systematic / n  # eq. 5-1a
            equation = f"eq. 5-1a, N / n = {s.n_systematic} / {n}"
        else:
            h, w = s.historic.h, s.historic.w
            p_tilde = (h - w * s.historic.l) / h  # eq. 5-1b
            equation = f"eq. 5-1b, (H - W L) / H with H = {h}, W = {w:.5f}, L = {s.historic.l}"
        if p_tilde <= 0.5:
            # Only after a historic adjustment that moves much of the record to Z. Without
            # one, the 25 percent limit keeps this at 0.75 or more.
            raise Bulletin17BError(
                f"the probability of a peak above the truncation level is {p_tilde:.4f} "
                f"({equation}). The synthetic statistics need the curve at exceedance "
                "probability 0.50, which the adjusted curve does not reach."
            )
        low_years = {y for y, _ in s.low_outliers}
        kept = min(q for y, q in zip(self._years, self._flows) if y not in low_years)
        known = [q for _, q in s.low_outliers] + [0.0 for _ in self._zero_years]
        level = f"at or below {kept:g} cfs, the smallest peak kept"
        if known:
            level = f"above {max(known):g} cfs, the largest known flow truncated, and {level}"
        steps.append(
            f"Made the conditional probability adjustment of Appendix 5: {truncated} of {n} "
            f"years truncated ({', '.join(self._truncated_parts(n_low))}). The truncation "
            f"level is {level}. P~ = {p_tilde:.4f} ({equation})."
        )
        skew_c = _nearest_tenth(s.skew) if self._hand else s.skew
        q01, q10, q50 = (
            float(10.0 ** (s.mean + kfactor(skew_c, p / p_tilde) * s.std))
            for p in (0.01, 0.10, 0.50)
        )
        mean_s, std_s, skew_s = synthetic_statistics(q01, q10, q50)
        rounding = ""
        if self._hand:
            k_skew = _nearest_tenth(skew_s)
            mean_s, std_s, skew_s = synthetic_statistics(q01, q10, q50, k_skew=k_skew)
            rounding = (
                f" K for the conditional curve is read at skew {skew_c:.1f}, and at "
                f"{k_skew:.1f} in eq. 5-4."
            )
        steps.append(
            f"The conditional curve (mean {s.mean:.4f}, standard deviation {s.std:.4f}, skew "
            f"{s.skew:.4f}) gives {q01:,.0f}, {q10:,.0f} and {q50:,.0f} cfs at exceedance "
            f"probabilities 0.01, 0.10 and 0.50 over all years.{rounding} The synthetic "
            f"statistics (eqs. 5-3 to 5-5) are mean {mean_s:.4f}, standard deviation "
            f"{std_s:.4f} and skew {skew_s:.4f}. They replace the statistics above."
        )
        return ConditionalAdjustment(
            p_tilde=p_tilde,
            n_truncated=truncated,
            mean=s.mean,
            std=s.std,
            skew=s.skew,
            q01=q01,
            q10=q10,
            q50=q50,
            mean_s=mean_s,
            std_s=std_s,
            skew_s=skew_s,
        )

    def _note_extrapolation(self, p_tilde: float, steps: list[str]) -> None:
        """Name the AEPs that lie below the truncation level."""
        beyond = [a for a in self._aeps if a > p_tilde]
        if beyond:
            steps.append(
                f"AEP {', '.join(f'{a:g}' for a in beyond)}: above P~ = {p_tilde:.4f}, so "
                "the flow extrapolates the synthetic curve below the truncation level. The "
                "curve is fitted to exceedance probabilities 0.01 to 0.50."
            )

    def _screen_and_adjust(self, steps: list[str]) -> _Screened:
        """Outlier tests and historic adjustment in the order of the p. 12-3 flow chart.

        The moments returned leave the low outliers out. The conditional probability
        adjustment that has to follow truncated years is not made here.
        """
        self._leave_out_unmeasured(steps)
        x = np.log10(np.asarray(self._flows, dtype=float))
        mean, std, g = log_moments(x)
        steps.append(
            f"Computed the mean ({mean:.4f}), standard deviation ({std:.4f}) and station "
            f"skew ({g:.4f}) of the base-10 logs of the {x.size} annual peaks."
        )
        order = _outlier_order(g)
        k_n = outlier_k(x.size)
        steps.append(
            f"The station skew is {g:.4f}, so the outlier tests run in the order: {order}. "
            f"K_N is {k_n:.3f} for N = {x.size}."
        )
        historic: Optional[HistoricAdjustment] = None
        no_lows = np.zeros(x.size, dtype=bool)
        if order == "low first":
            lo_limit, low = self._low_test(x, mean - k_n * std, "eq. 8a", "", steps)
            if low.any():
                mean, std, g = self._without_lows(x, low)
                k_n = outlier_k(int((~low).sum()))
            hi_log = mean + k_n * std
            high = x > hi_log
            if low.any():
                self._note_dropped(low, (mean, std, g), self._adjusts(high), steps)
                steps[-1] += f" The high outlier test uses them, with K_N = {k_n:.3f}."
            self._note_test("high", "eq. 7", float(10.0**hi_log), high, steps)
            historic = self._adjust_for_history(x, high, low, steps)
        elif order == "high first":
            hi_log = mean + k_n * std
            high = x > hi_log
            self._note_test("high", "eq. 7", float(10.0**hi_log), high, steps)
            historic = self._adjust_for_history(x, high, no_lows, steps)
            if historic is None:
                lo_limit, low = self._low_test(x, mean - k_n * std, "eq. 8a", "", steps)
            else:
                lo_limit, low = self._low_test(
                    x,
                    historic.mean - k_n * historic.std,
                    "eq. 8b",
                    f" on the adjusted statistics, K = {k_n:.3f} for the {x.size} "
                    "systematic peaks tested",
                    steps,
                )
            if low.any():
                mean, std, g = self._without_lows(x, low)
                self._note_dropped(low, (mean, std, g), historic is not None, steps)
                if historic is not None:
                    historic = self._adjust_for_history(x, high, low, steps, again=True)
        else:
            hi_log = mean + k_n * std
            high = x > hi_log
            self._note_test("high", "eq. 7", float(10.0**hi_log), high, steps)
            lo_limit, low = self._low_test(x, mean - k_n * std, "eq. 8a", "", steps)
            if low.any():
                mean, std, g = self._without_lows(x, low)
                self._note_dropped(low, (mean, std, g), self._adjusts(high), steps)
            historic = self._adjust_for_history(x, high, low, steps)
        n_systematic = int((~low).sum())
        if historic is not None:
            mean, std, g = historic.mean, historic.std, historic.skew
            n_systematic -= historic.z - len(self._historic_peaks)
        return _Screened(
            order=order,
            high_threshold_log=hi_log,
            low_threshold=lo_limit,
            high_outliers=self._peaks_where(high),
            low_outliers=self._peaks_where(low),
            historic=historic,
            n_systematic=n_systematic,
            mean=mean,
            std=std,
            skew=g,
        )

    def _leave_out_unmeasured(self, steps: list[str]) -> None:
        """Record the zero and below-base years, and refuse when they leave too little."""
        if not self._n_unmeasured:
            return
        parts = []
        if self._zero_years:
            parts.append(f"{len(self._zero_years)} zero year(s) ({_year_list(self._zero_years)})")
        if self._below_base_years:
            parts.append(
                f"{len(self._below_base_years)} year(s) below the gage base "
                f"({_year_list(self._below_base_years)})"
            )
        steps.append(f"Left out of the statistics: {' and '.join(parts)}.")
        if self._flows_ignored:
            steps.append(
                f"The flow given for below-base year(s) {_year_list(self._flows_ignored)} is "
                "not used: the year is below the gage base, whatever the flow."
            )
        self._require_mostly_kept(0)
        none = np.zeros(len(self._flows), dtype=bool)
        self._require_enough(none, none)

    def _peaks_where(self, mask: np.ndarray) -> list[tuple[int, float]]:
        """The systematic peaks a mask selects, by year."""
        return [(y, q) for y, q, m in zip(self._years, self._flows, mask) if m]

    def _low_test(
        self, x: np.ndarray, limit_log: float, equation: str, detail: str, steps: list[str]
    ) -> tuple[float, np.ndarray]:
        """The low-outlier threshold in cfs and the peaks below it.

        A threshold the caller gave replaces that of ``equation``, and is compared with
        the flows themselves.
        """
        if self._low_threshold is not None:
            limit = self._low_threshold
            low = np.asarray(self._flows) < limit
            label = f"the threshold given, which replaces {equation}"
        else:
            limit = float(10.0**limit_log)
            low = x < limit_log
            label = equation + detail
        self._note_test("low", label, limit, low, steps)
        return limit, low

    def _note_test(
        self, kind: str, equation: str, limit_cfs: float, flagged: np.ndarray, steps: list[str]
    ) -> None:
        """Record one outlier test and what it flagged."""
        side = "above" if kind == "high" else "below"
        shown = f"{limit_cfs:,.0f}" if limit_cfs >= 100 else f"{limit_cfs:.3g}"
        limit = f"{side} {shown} cfs ({equation})"
        if not flagged.any():
            steps.append(f"No {kind} outliers {limit}.")
            return
        listed = ", ".join(f"{y} ({q:g} cfs)" for y, q in self._peaks_where(flagged))
        text = f"{kind.capitalize()} outlier(s) {limit}: {listed}."
        if kind == "high" and self._historic_period is None:
            text += (
                " With no historic period given they are treated as part of the systematic "
                "record and kept."
            )
        steps.append(text)

    def _without_lows(self, x: np.ndarray, low: np.ndarray) -> tuple[float, float, float]:
        """Moments of the peaks left when the low outliers are dropped."""
        self._require_mostly_kept(int(low.sum()))
        self._require_enough(np.zeros(x.size, dtype=bool), low)
        return log_moments(x[~low])

    @staticmethod
    def _note_dropped(
        low: np.ndarray, moments: tuple[float, float, float], adjusted: bool, steps: list[str]
    ) -> None:
        """Record the moments without the low outliers, and whether they are final."""
        mean, std, g = moments
        text = (
            f"Dropped the {int(low.sum())} low outlier(s) and recomputed the mean "
            f"({mean:.4f}), standard deviation ({std:.4f}) and skew ({g:.4f})."
        )
        if adjusted:
            text += (
                " These still include any high outlier and are not used for the curve, which "
                "takes the historically adjusted statistics that follow."
            )
        steps.append(text)

    def _require_enough(self, in_z: np.ndarray, low: np.ndarray) -> None:
        """Refuse when fewer than ``MIN_SYSTEMATIC`` peaks are left for the moments."""
        left = int((~in_z & ~low).sum())
        if left >= MIN_SYSTEMATIC:
            return
        n_moved, n_low = int(in_z.sum()), int(low.sum())
        why = [f"{n_moved} moved to the historic peaks"] if n_moved else []
        why += self._truncated_parts(n_low)
        text = (
            f"only {left} of the {len(self._all_years)} years of the systematic record are "
            f"left for the moments{': ' if why else ''}{', '.join(why)}. "
            f"Bulletin 17B needs at least {MIN_SYSTEMATIC}."
        )
        if n_moved and self._historic_peaks:
            year, flow = min(self._historic_peaks, key=lambda e: e[1])
            n = len(self._flows)
            above = sum(q > flow for q in self._flows)
            if above > n / 2:
                text += (
                    f" The historic peak of {flow:g} cfs in {year} is below {above} of the "
                    f"{n} systematic peaks, so it may be in the wrong units."
                )
        raise Bulletin17BError(text)

    def _adjusts(self, high: np.ndarray) -> bool:
        """True when there is a historic period and something in it to adjust for."""
        return self._historic_period is not None and bool(self._historic_peaks or high.any())

    def _adjust_for_history(
        self,
        x: np.ndarray,
        high: np.ndarray,
        low: np.ndarray,
        steps: list[str],
        *,
        again: bool = False,
    ) -> Optional[HistoricAdjustment]:
        """Make the Appendix 6 adjustment, or return None when there is nothing to adjust.

        Z takes the given historic peaks, the high outliers, and every systematic peak
        above the smallest given historic peak. A systematic peak equal to it stays in
        the systematic record. L is the number of zeros, below-base years and low
        outliers. ``again`` marks a second adjustment, made because the first one brought
        out low outliers.
        """
        if self._historic_period is None:
            return None
        first, last = self._historic_period
        if not self._adjusts(high):
            steps.append(
                f"The historic period {first} to {last} has no historic peaks and no high "
                "outliers, so no historic adjustment was made."
            )
            return None
        in_z = high.copy()
        if self._historic_peaks:
            smallest = min(q for _, q in self._historic_peaks)
            in_z |= np.asarray(self._flows) > smallest
        in_z &= ~low
        self._require_enough(in_z, low)
        peaks = sorted(self._historic_peaks + self._peaks_where(in_z))
        h = last - first + 1
        n_truncated = self._n_unmeasured + int(low.sum())
        hm = historic_moments(
            x[~in_z & ~low], np.log10(np.array([q for _, q in peaks])), h=h, l=n_truncated
        )
        moved = self._peaks_where(in_z)
        listed = ", ".join(f"{y} ({q:g} cfs)" for y, q in moved)
        steps.append(
            f"Made the historic adjustment of Appendix 6 for {first} to {last}: H = {h}, "
            f"Z = {len(peaks)}, L = {n_truncated}, W = {hm.w:.5f}."
            + (" It replaces the first adjustment." if again else "")
            + (f" Systematic peaks treated as historic: {listed}." if moved else "")
            + f" Adjusted mean {hm.mean:.4f}, standard deviation {hm.std:.4f} and skew "
            f"{hm.skew:.4f}."
        )
        return HistoricAdjustment(
            h=h,
            z=len(peaks),
            l=n_truncated,
            w=hm.w,
            mean=hm.mean,
            std=hm.std,
            skew=hm.skew,
            historic_peaks=peaks,
        )

    def _weight_skew(
        self, skew_station: float, skew_mse: float, label: str, steps: list[str]
    ) -> tuple[Optional[float], float]:
        """Weight the station and generalized skews (eq. 5), or use the station skew."""
        if self._gen_skew is None:
            steps.append(
                f"No generalized skew was given, so the {label} {skew_station:.4f} "
                "is used for the curve."
            )
            return None, skew_station
        w = weighted_skew(skew_station, skew_mse, self._gen_skew, self._gen_mse)
        used = _nearest_tenth(w) if self._hand else w + 0.0  # no -0.0
        text = (
            f"Weighted the {label} {skew_station:.4f} (MSE {skew_mse:.3f}) with the "
            f"generalized skew {self._gen_skew:.4f} (MSE {self._gen_mse:.3f}) by inverse "
            f"MSE to get {w:.4f}."
        )
        if self._hand:
            text += f" Rounded to {used:.1f}, as the worked examples do."
        steps.append(text)
        return w, used

    def _frequency_curve(
        self,
        mean: float,
        std: float,
        skew: float,
        n_limits: int,
        n_prob: int,
        steps: list[str],
    ) -> pd.DataFrame:
        """Quantiles, confidence limits and expected probability at the skew used."""
        aep = np.asarray(self._aeps, dtype=float)
        k = np.array([kfactor(skew, float(a)) for a in aep])
        log_q = mean + k * std
        z = float(ndtri(1.0 - (1.0 - self._confidence) / 2.0))
        k_up, k_lo = confidence_deviates(k, n_limits, z)
        steps.append(
            f"Computed the log-Pearson III quantiles at skew {skew:.4f}, with "
            f"{100 * self._confidence:g} percent two-sided confidence limits for N = "
            f"{n_limits} and the expected probability for N = {n_prob}."
        )
        return pd.DataFrame(
            {
                "aep": aep,
                "flow_cfs": 10.0**log_q,
                "log_flow": log_q,
                "K_factor": k,
                "lower": 10.0 ** (mean + k_lo * std),
                "upper": 10.0 ** (mean + k_up * std),
                "expected_probability": expected_probability(aep, n_prob),
            }
        )

    def _plotting_positions(self, historic: Optional[HistoricAdjustment]) -> pd.DataFrame:
        """Weibull positions, weighted by Appendix 6 when ``historic`` is given.

        Events rank from the largest, the Z historic peaks first (eqs. 6-6 to 6-8 with
        a = 0). Without an adjustment, W = 1, Z = 0 and H is the number of years in the
        systematic record, n, and this is rank / (n + 1). Truncated years take ranks:
        below-base years come after every measured peak and zeros come last.
        """
        top = [] if historic is None else historic.historic_peaks
        in_top = {y for y, _ in top}
        rest = [(y, q) for y, q in zip(self._years, self._flows) if y not in in_top]
        events = sorted(top, key=lambda e: -e[1]) + sorted(rest, key=lambda e: -e[1])
        events += [(y, float("nan")) for y in self._below_base_years]
        events += [(y, 0.0) for y in self._zero_years]
        z = len(top)
        w = 1.0 if historic is None else historic.w
        h = len(events) if historic is None else historic.h
        rank = np.arange(1, len(events) + 1, dtype=float)
        order = np.where(rank <= z, rank, w * rank - (w - 1.0) * (z + 0.5))
        table = pd.DataFrame(
            {
                "year": [y for y, _ in events],
                "flow_cfs": [q for _, q in events],
                "weighted_order": order,
                "exceedance": order / (h + 1.0),
            }
        )
        return table.sort_values("year", ignore_index=True)


def _threshold_row(row: object) -> tuple[int, int, float]:
    """One perception-threshold row as (start year, end year, threshold in cfs).

    The threshold may be infinite (a period that was not observed); it is not negative.
    """
    if not isinstance(row, Mapping) or not {"start_year", "end_year", "threshold_cfs"} <= set(row):
        raise Bulletin17BError(
            "each perception threshold needs start_year, end_year and threshold_cfs; "
            f"got {row!r}"
        )
    start = _integer_year(row["start_year"], "threshold start year")
    end = _integer_year(row["end_year"], "threshold end year")
    try:
        threshold = float(row["threshold_cfs"])
    except (TypeError, ValueError):
        raise Bulletin17BError(
            f"a threshold must be a number; got {row['threshold_cfs']!r}"
        ) from None
    if threshold != threshold:
        raise Bulletin17BError("a threshold must not be NaN")
    if end < start:
        raise Bulletin17BError(f"a perception threshold ends before it begins: {start} to {end}")
    if threshold < 0.0:
        raise Bulletin17BError(f"a perception threshold must not be negative; got {threshold:g}")
    return start, end, threshold


def _as_list(value: object, name: str) -> list[Any]:
    """``value`` as a list, refusing text, mappings, frames and anything not a sequence."""
    if isinstance(value, (str, bytes, Mapping, pd.DataFrame)):
        raise Bulletin17BError(f"{name} must be a sequence; got {type(value).__name__}")
    try:
        return list(value)  # type: ignore[call-overload]
    except TypeError:
        raise Bulletin17BError(f"{name} must be a sequence; got {value!r}") from None


def _span(start: int, end: int) -> str:
    return f"{start} to {end}"


def _historic_row_start(
    rows: list[tuple[int, int, float]],
    measured: Mapping[int, float],
    first: int,
    last: int,
    record_start: int,
) -> Optional[int]:
    """The start of the one historic row, after checking the rest of the finite rows.

    Rows of zero are dropped and a repeated row counts once. A row that touches the
    systematic record, ``first`` to ``last``, may only be one that changes nothing there:
    every year it covers there has a measured peak, and its threshold is at or below
    all of them. ``record_start`` is the first year of the record before peaks were
    removed, where the historic row has to end.
    """
    historic: list[tuple[int, int, float]] = []
    for start, end, threshold in sorted(set(r for r in rows if r[2] > 0.0)):
        if start > last:
            raise Bulletin17BError(
                f"the perception threshold {_span(start, end)} begins after the last year of "
                f"the record, {last}; Bulletin 17B has no use for it"
            )
        lo, hi = max(start, first), min(end, last)
        if lo <= hi:
            unmeasured = [y for y in range(lo, hi + 1) if y not in measured]
            covered = [measured[y] for y in range(lo, hi + 1) if y in measured]
            if unmeasured:
                raise Bulletin17BError(
                    f"the perception threshold of {threshold:g} cfs for {_span(start, end)} "
                    f"covers {len(unmeasured)} year(s) of the systematic record with no "
                    f"measured peak ({_span(unmeasured[0], unmeasured[-1])}"
                    f"{'' if len(unmeasured) == hi - lo + 1 else ' and others'}); "
                    "dropping the row would lose that information, and Bulletin 17B cannot "
                    "take it"
                )
            if threshold > min(covered):
                raise Bulletin17BError(
                    f"the perception threshold of {threshold:g} cfs for {_span(start, end)} is "
                    "a change of threshold inside the systematic record, which Bulletin 17B "
                    "cannot take"
                )
        if start < first:
            historic.append((start, end, threshold))
    if len(historic) > 1:
        spans = ", ".join(_span(s, e) for s, e, _ in historic)
        raise Bulletin17BError(
            f"the perception thresholds give more than one historic period ({spans}); "
            "Bulletin 17B takes one historic period"
        )
    if not historic:
        return None
    start, end, _ = historic[0]
    if end < record_start - 1:
        raise Bulletin17BError(
            f"the historic threshold ends in {end} but the record starts in {record_start}: "
            f"{_span(end + 1, record_start - 1)} is covered by no period"
        )
    return start


class _Coded(NamedTuple):
    """One peak with its codes read."""

    year: int
    flow: float
    given: object
    interval: PeakInterval


def _read_peaks(
    peak_flows: Sequence[float],
    water_years: Sequence[int],
    peak_codes: Optional[Sequence[object]],
) -> list[_Coded]:
    """Check the three lists and classify every peak by its codes."""
    flows_raw = _as_list(peak_flows, "peak_flows")
    years = [_integer_year(y, "water year") for y in _as_list(water_years, "water_years")]
    if len(flows_raw) != len(years):
        raise Bulletin17BError(
            f"{len(flows_raw)} peak flows but {len(years)} water years: lengths must match"
        )
    codes = [None] * len(years) if peak_codes is None else _as_list(peak_codes, "peak_codes")
    if len(codes) != len(years):
        raise Bulletin17BError(
            f"{len(codes)} peak codes but {len(years)} peaks: lengths must match"
        )
    out = []
    for year, raw, code in zip(years, flows_raw, codes):
        try:
            flow = float(raw)
        except (TypeError, ValueError):
            raise Bulletin17BError(
                f"the peak flow of {year} must be a number; got {raw!r}"
            ) from None
        if flow < 0.0:
            raise Bulletin17BError(f"the peak flow of {year} is negative ({flow:g})")
        as_code = (
            str(int(code))
            if isinstance(code, (int, np.integer)) and not isinstance(code, (bool, np.bool_))
            else code
        )
        if isinstance(code, (Mapping, pd.DataFrame, pd.Series, bytes)):
            raise Bulletin17BError(f"the peak code of {year} is not a code: {code!r}")
        try:
            interval = peak_interval(0.0 if flow != flow else flow, as_code)
        except (TypeError, ValueError):
            raise Bulletin17BError(f"the peak code of {year} is not a code: {code!r}") from None
        out.append(_Coded(year, flow, code, interval))
    return out


def removed_peaks(
    peak_flows: Sequence[float],
    water_years: Sequence[int],
    peak_codes: Optional[Sequence[object]],
) -> list[tuple[int, float, str]]:
    """The peaks that :func:`bulletin17b_kwargs` removes for their codes.

    Parameters
    ----------
    peak_flows, water_years, peak_codes
        As for :func:`bulletin17b_kwargs`.

    Returns
    -------
    list of (int, float, str)
        Year, flow in cfs (NaN if none was given) and the codes as given, for each peak
        with a code 3, O, 6 or C, in the order given.

    Raises
    ------
    Bulletin17BError
        On lists that are not sequences, differ in length or hold bad values.
    """
    return [
        (
            p.year,
            p.flow,
            p.given if isinstance(p.given, str) else ",".join(sorted(p.interval.codes)),
        )
        for p in _read_peaks(peak_flows, water_years, peak_codes)
        if p.interval.treatment is PeakTreatment.REMOVED
    ]


def _years_text(years: Sequence[int]) -> str:
    return ", ".join(str(y) for y in sorted(years))


def bulletin17b_kwargs(
    peak_flows: Sequence[float],
    water_years: Sequence[int],
    peak_codes: Optional[Sequence[object]],
    perception_thresholds: Optional[Sequence[Mapping[str, Any]]] = None,
) -> dict[str, Any]:
    """:class:`Bulletin17B` arguments from an NWIS peak record and its qualification codes.

    The codes are classified as :func:`flowfreq.peak_codes.peak_interval` classifies them
    for Bulletin 17C, so a combination of codes resolves the same way:

    - code 7 peaks are historic peaks;
    - code 4 years are below the gage base: they go to ``below_base_years`` and the flow
      printed for them is not used;
    - codes 3 and O, and 6 and C, remove the peak, as ``peak_code_kwargs`` does by
      default. Removal wins over every other code on the same peak. The removal is logged
      with the years, and :func:`removed_peaks` lists them. A regulated record is refused
      elsewhere, in the report section, not here;
    - code 8, and code 7 with 4 or 8, have no Bulletin 17B treatment and are refused;
    - a peak with no flow (NaN), apart from code 4 and the removed ones, is dropped with a
      logged warning that names its year. A historic peak with no flow is dropped as
      well, and a below-base year never needs a flow.

    The historic period is the one perception-threshold row that starts before the
    systematic record and ends in the year before it, through the last year of the
    record, or the last historic peak if that is later (a historic peak after the gage
    closed). The threshold's value is not checked against the historic peaks: one above
    some of them, or a historic peak missing from the period, goes unnoticed. Without a
    row, historic peaks give a period from the first historic year (or the first year of
    the record, if that is earlier). A row with no historic peaks still gives the period,
    so that high outliers can be adjusted for.

    Rows of zero are dropped, as ``run_ffa`` does, and a repeated row counts once. A
    positive row that touches the systematic record is accepted only if it changes
    nothing there: each year it covers has a measured peak, and its threshold is at or
    below all of them. A row that covers a year of the record with no measured peak
    (a gap, a below-base year, a historic peak or a removed peak) is refused, because
    dropping it would lose information; so are a change of threshold inside the record
    and two historic rows. An infinite threshold is an unobserved period, which Bulletin
    17B has no treatment for. It is refused unless no year in it has a peak and none of
    it is in the historic period, in which case it is ignored and the choice is logged.
    The systematic record runs from the first to the last year of the peaks that are not
    historic, counting below-base years. A historic peak may fall in a gap in it, but
    not in a year that has a systematic peak.

    Parameters
    ----------
    peak_flows : sequence of float
        Annual peak flows in cfs.
    water_years : sequence of int
        Water year of each peak.
    peak_codes : sequence, optional
        One NWIS qualification code field per peak (``"7"``, ``"2,6"``, an integer code,
        NaN for none).
    perception_thresholds : sequence of mapping, optional
        Rows with ``start_year``, ``end_year`` and ``threshold_cfs``, the shape the
        ``gaged_ffa`` report section takes.

    Returns
    -------
    dict
        ``peak_flows``, ``water_years``, ``historic_peaks``, ``historic_period`` and
        ``below_base_years``, as :class:`Bulletin17B` takes them.

    Raises
    ------
    Bulletin17BError
        On an input that is not a sequence, a code with no Bulletin 17B treatment, a year
        that occurs twice, thresholds Bulletin 17B cannot take, or no peaks left.
    """
    peaks = _read_peaks(peak_flows, water_years, peak_codes)
    rows = [
        _threshold_row(r)
        for r in _as_list(
            [] if perception_thresholds is None else perception_thresholds, "perception_thresholds"
        )
    ]

    kept: list[tuple[int, float]] = []
    historic: list[tuple[int, float]] = []
    below: list[int] = []
    removed: list[int] = []
    no_flow: list[int] = []
    no_flow_historic: list[int] = []
    for p in peaks:
        interval = p.interval
        if interval.treatment is PeakTreatment.REMOVED:
            removed.append(p.year)
            continue
        if interval.treatment is PeakTreatment.GREATER_THAN or (
            interval.is_historic and interval.codes & (LESS_THAN_CODES | GREATER_THAN_CODES)
        ):
            raise Bulletin17BError(
                f"the {p.year} peak has code(s) {','.join(sorted(interval.codes))}, which have "
                "no Bulletin 17B treatment"
            )
        if interval.treatment is PeakTreatment.LESS_THAN:
            below.append(p.year)
        elif p.flow != p.flow:
            (no_flow_historic if interval.is_historic else no_flow).append(p.year)
        else:
            (historic if interval.is_historic else kept).append((p.year, p.flow))
    if removed:
        logger.info(
            "Peaks removed for their codes (3, O, 6 or C), %d: %s",
            len(removed),
            _years_text(removed),
        )
    for label, years_dropped in (("peak", no_flow), ("historic peak", no_flow_historic)):
        if years_dropped:
            logger.warning(
                "Dropped %d %s(s) with no flow: %s",
                len(years_dropped),
                label,
                _years_text(years_dropped),
            )

    systematic_years = [y for y, _ in kept] + below
    if not systematic_years:
        what = "no peaks are left" if not historic else "no systematic peaks are left"
        raise Bulletin17BError(f"{what} after the peak codes are applied")
    _require_unique_years(systematic_years + [y for y, _ in historic], "a water year")
    first, last = min(systematic_years), max(systematic_years)
    # The record before peaks were removed: a removed first year does not open a gap.
    record_start = min(p.year for p in peaks if not p.interval.is_historic)

    finite = [r for r in rows if np.isfinite(r[2])]
    start = _historic_row_start(finite, dict(kept), first, last, record_start)
    period: Optional[tuple[int, int]] = None
    last_historic = max((y for y, _ in historic), default=last)
    if start is not None:
        period = (start, max(last, last_historic))
    elif historic:
        period = (min(first, min(y for y, _ in historic)), max(last, last_historic))
    if period is not None:
        outside = sorted(y for y, _ in historic if not period[0] <= y <= period[1])
        if outside:
            raise Bulletin17BError(
                f"historic peak(s) in {outside} are outside the historic period "
                f"{_span(*period)}"
            )
    with_peaks = set(systematic_years) | {y for y, _ in historic}
    for s_, e_, _ in sorted(set(r for r in rows if not np.isfinite(r[2]))):
        in_period = period is not None and s_ <= period[1] and e_ >= period[0]
        if in_period or any(s_ <= y <= e_ for y in with_peaks):
            raise Bulletin17BError(
                f"the perception threshold {_span(s_, e_)} is infinite, an unobserved period, "
                "and it covers years with peaks or the historic period; Bulletin 17B has no "
                "treatment for an unobserved period"
            )
        logger.info(
            "Perception threshold %s is infinite and no year in it has a peak or is in the "
            "historic period: ignored, those years are not part of the record",
            _span(s_, e_),
        )
    kept.sort()
    return {
        "peak_flows": [q for _, q in kept],
        "water_years": [y for y, _ in kept],
        "historic_peaks": sorted(historic),
        "historic_period": period,
        "below_base_years": sorted(below),
    }
