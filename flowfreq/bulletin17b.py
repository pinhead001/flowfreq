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

import logging
from dataclasses import dataclass, field
from typing import Mapping, NamedTuple, Optional, Sequence

import numpy as np
import pandas as pd
from scipy.special import ndtri
from scipy.stats import t as student_t

from flowfreq.core import kfactor

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


#: The bulletin's frequency-table rows (Table 12-3), as annual exceedance probabilities.
STANDARD_AEP: tuple[float, ...] = (0.99, 0.90, 0.50, 0.10, 0.05, 0.02, 0.01, 0.005, 0.002)


@dataclass
class HistoricAdjustment:
    """Historic adjustment of Appendix 6 (filled in when the adjustment is made)."""

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
    """Conditional probability adjustment of Appendix 5 (filled in when it is made)."""

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
        Years left out of the moments as zeros or below the gage base.
    outlier_order : str
        ``"high first"``, ``"low first"`` or ``"both on unadjusted statistics"``.
    high_outlier_threshold, low_outlier_threshold : float
        Outlier thresholds in cfs.
    high_outliers, low_outliers : list of (year, flow)
        Peaks the tests flagged.
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
        The skew the quantiles use (rounded to a tenth in hand-computation mode).
    quantiles : pandas.DataFrame
        Columns ``aep, flow_cfs, log_flow, K_factor, lower, upper,
        expected_probability``. ``lower`` and ``upper`` are the confidence limits in cfs.
    plotting_positions : pandas.DataFrame
        Columns ``year, flow_cfs, weighted_order, exceedance`` (Weibull here).
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


def _not_implemented(what: str) -> Bulletin17BError:
    return Bulletin17BError(f"{what} is not implemented yet in Bulletin17B")


class Bulletin17B:
    """Bulletin 17B (1982) flood frequency analysis, superseded by Bulletin 17C.

    This version handles a complete systematic record: moments of the log peaks,
    the outlier tests, skew weighting, the log-Pearson III curve, confidence limits
    and expected probability. Records that need the historic adjustment, zero or
    low-outlier treatment are refused until those steps exist.

    Parameters
    ----------
    peak_flows : sequence of float
        Annual peak flows in cfs.
    water_years : sequence of int
        Water year of each peak. A year may appear only once.
    generalized_skew : float, optional
        Generalized (regional) skew. Without it the station skew is used.
    generalized_skew_mse : float, default 0.302
        Mean square error of the generalized skew (Plate I value).
    historic_peaks, historic_period, below_base_years, low_outlier_threshold
        Inputs of later steps; giving one raises :class:`Bulletin17BError`.
    confidence : float, default 0.90
        Two-sided confidence level of the limits.
    aeps : sequence of float, optional
        Annual exceedance probabilities to report; ``STANDARD_AEP`` by default.
    hand_computation : bool, default False
        Round the weighted skew to the nearest tenth before reading K, as the
        worked examples in the bulletin do. Exact mode never rounds.

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
        raw_years = list(water_years)
        if len(flows) != len(raw_years):
            raise Bulletin17BError(
                f"{len(flows)} peak flows but {len(raw_years)} water years: lengths must match"
            )
        if any(q < 0 for q in flows):
            raise Bulletin17BError("peak flows must not be negative")
        years: list[int] = []
        for y in raw_years:
            fy = _finite_number(y, "a water year")
            if isinstance(y, bool) or fy != int(fy):
                raise Bulletin17BError(f"water years must be integers; got {y!r}")
            years.append(int(fy))
        if len(set(years)) != len(years):
            dup = sorted({y for y in years if years.count(y) > 1})
            raise Bulletin17BError(f"a water year appears more than once: {dup}")
        order = np.argsort(years, kind="stable")
        self._years = [years[i] for i in order]
        self._flows = [flows[i] for i in order]

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
        self._aeps = aep_list
        self._hand = bool(hand_computation)
        self._historic_peaks = historic_peaks
        self._historic_period = historic_period
        self._below_base_years = below_base_years
        self._low_threshold = low_outlier_threshold

    def run_analysis(self) -> Bulletin17BResults:
        """Fit the record and return the frequency table, limits and steps.

        Raises
        ------
        Bulletin17BError
            If the record is constant, has fewer than 10 or more than 149 peaks, or
            needs a step that is not implemented yet.
        """
        steps: list[str] = []
        self._require_supported()
        x = np.log10(np.asarray(self._flows, dtype=float))
        n = x.size
        mean, std, g = log_moments(x)
        steps.append(
            f"Computed the mean ({mean:.4f}), standard deviation ({std:.4f}) and station "
            f"skew ({g:.4f}) of the base-10 logs of the {n} annual peaks."
        )
        order, high, low = self._screen_outliers(x, n, mean, std, g, steps)
        # Later tasks add their steps here, each as a separate method that returns the
        # adjusted moments: the historic adjustment (Appendix 6), then the low-outlier,
        # zero and conditional probability adjustments (Appendix 5).
        skew_station = g
        skew_mse = station_skew_mse(skew_station, n)
        skew_w, skew_used = self._weight_skew(skew_station, skew_mse, steps)
        quantiles = self._frequency_curve(mean, std, skew_used, n, n, steps)
        k_n = outlier_k(n)
        return Bulletin17BResults(
            n_years=n,
            n_systematic=n,
            n_zeros=0,
            n_below_base=0,
            outlier_order=order,
            high_outlier_threshold=float(10.0 ** (mean + k_n * std)),
            low_outlier_threshold=float(10.0 ** (mean - k_n * std)),
            high_outliers=high,
            low_outliers=low,
            historic=None,
            conditional=None,
            mean_log=mean,
            std_log=std,
            skew_station=skew_station,
            skew_station_mse=skew_mse,
            skew_generalized=self._gen_skew,
            skew_generalized_mse=None if self._gen_skew is None else self._gen_mse,
            skew_weighted=skew_w,
            skew_used=skew_used,
            quantiles=quantiles,
            plotting_positions=self._plotting_positions(),
            steps=steps,
        )

    def _require_supported(self) -> None:
        """Refuse inputs only the later steps can treat, rather than answer wrongly."""
        if self._historic_peaks:
            raise _not_implemented("The historic adjustment (historic_peaks)")
        if self._historic_period is not None:
            raise _not_implemented("The historic adjustment (historic_period)")
        if self._below_base_years:
            raise _not_implemented("The conditional probability adjustment (below_base_years)")
        if self._low_threshold is not None:
            raise _not_implemented("A user-supplied low outlier threshold")
        if any(q == 0.0 for q in self._flows):
            raise _not_implemented("The treatment of zero flows")

    def _screen_outliers(
        self,
        x: np.ndarray,
        n: int,
        mean: float,
        std: float,
        g: float,
        steps: list[str],
    ) -> tuple[str, list[tuple[int, float]], list[tuple[int, float]]]:
        """Apply the Appendix 4 tests in the order the station skew sets (p. 17)."""
        k_n = outlier_k(n)
        if g > 0.4:
            order = "high first"
        elif g < -0.4:
            order = "low first"
        else:
            order = "both on unadjusted statistics"
        steps.append(
            f"The station skew is {g:.4f}, so the outlier tests run in the order: {order}. "
            f"K_N is {k_n:.3f} for N = {n}."
        )
        hi_log = mean + k_n * std
        lo_log = mean - k_n * std
        high = [(y, q) for y, q, lx in zip(self._years, self._flows, x) if lx > hi_log]
        low = [(y, q) for y, q, lx in zip(self._years, self._flows, x) if lx < lo_log]
        if high:
            listed = ", ".join(f"{y} ({q:g} cfs)" for y, q in high)
            steps.append(
                f"High outlier(s) above {10 ** hi_log:,.0f} cfs: {listed}. With no historic "
                "period given they are treated as part of the systematic record and kept."
            )
        else:
            steps.append(f"No high outliers above {10 ** hi_log:,.0f} cfs.")
        if low:
            listed = ", ".join(f"{y} ({q:g} cfs)" for y, q in low)
            raise _not_implemented(
                f"The low outlier treatment (the test flags {listed} below "
                f"{10 ** lo_log:,.0f} cfs)"
            )
        steps.append(f"No low outliers below {10 ** lo_log:,.0f} cfs.")
        return order, high, low

    def _weight_skew(
        self, skew_station: float, skew_mse: float, steps: list[str]
    ) -> tuple[Optional[float], float]:
        """Weight the station and generalized skews (eq. 5), or use the station skew."""
        if self._gen_skew is None:
            steps.append(
                f"No generalized skew was given, so the station skew {skew_station:.4f} "
                "is used for the curve."
            )
            return None, skew_station
        w = weighted_skew(skew_station, skew_mse, self._gen_skew, self._gen_mse)
        used = (round(w, 1) if self._hand else w) + 0.0  # + 0.0 turns -0.0 into 0.0
        text = (
            f"Weighted the station skew {skew_station:.4f} (MSE {skew_mse:.3f}) with the "
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

    def _plotting_positions(self) -> pd.DataFrame:
        """Weibull plotting positions, rank / (N + 1) from the largest peak."""
        flows = np.asarray(self._flows, dtype=float)
        n = flows.size
        rank = np.empty(n)
        rank[np.argsort(-flows, kind="stable")] = np.arange(1, n + 1)
        return pd.DataFrame(
            {
                "year": self._years,
                "flow_cfs": flows,
                "weighted_order": rank,
                "exceedance": rank / (n + 1.0),
            }
        )
