"""Bulletin 17B (1982) formulas as pure functions.

Bulletin 17B, "Guidelines for Determining Flood Flow Frequency", is the method of
moments procedure that Bulletin 17C (2019) superseded. It is kept here to reproduce
older studies and to compare them with the Bulletin 17C fit. It is superseded, and
nothing in flowfreq runs it unless the caller asks for it.

This module holds the bulletin's equations, one function each, with the equation
numbers in the docstrings. All logs are base 10. The procedure that strings the
equations together lives elsewhere. Page labels and the places where the bulletin is
ambiguous are in ``tests/fixtures/b17b/procedure.md``.
"""

from __future__ import annotations

from typing import Mapping, NamedTuple, Optional

import numpy as np
from scipy.special import ndtri
from scipy.stats import t as student_t

from flowfreq.core import kfactor

__all__ = [
    "Bulletin17BError",
    "HistoricMoments",
    "KN_TABLE",
    "confidence_deviates",
    "expected_probability",
    "historic_moments",
    "log_moments",
    "outlier_k",
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
