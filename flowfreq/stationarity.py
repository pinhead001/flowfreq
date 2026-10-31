"""
Stationarity screening for annual series (roadmap section 6.1), and a mixed-population
seasonality screen.

Bulletin 17C assumes a stationary record and asks the analyst to evaluate that
assumption; it does not prescribe a method. These are the detection tests in common use
(and in USACE's Nonstationarity Detection Tool), for screening only: a detected trend or
change point is a reason to look at the record, not a correction.

- :func:`mann_kendall`: the Mann-Kendall trend test with the tie correction, optionally
  with the Hamed and Rao (1998) variance correction for serial correlation.
- :func:`sens_slope`: Sen's (1968) slope and intercept.
- :func:`spearman_trend`: Spearman's rho against time.
- :func:`pettitt`: the Pettitt (1979) change-point test, with its approximate p-value.
- :func:`stationarity_screen`: all of them on one series (in log10 for flows), with
  plain-language flags.
- :func:`mixed_population_screen`: how annual peaks split between the cool season
  (rain, rain-on-snow) and the warm season (snowmelt), and the concentration of their
  dates (Burn, 1997). A split record is a warning that one LP3 distribution may be
  fitted to two flood populations (Bulletin 17C section 9.3); it is not a test.

Not here yet (roadmap section 6.1): CUSUM and Lombard change points, trends in variance
and in peaks-over-threshold counts, field significance across a region, and the
attribution overlay on regulation and land-use dates.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)

ArrayLike = Union[Sequence[float], npt.ArrayLike]


def _series(x: ArrayLike, min_n: int = 3) -> np.ndarray:
    a = np.asarray(x, dtype=float).ravel()
    if a.size < min_n:
        raise ValueError(f"need at least {min_n} values, got {a.size}")
    if not np.all(np.isfinite(a)):
        raise ValueError("the series has missing or non-finite values; drop them first")
    return a


@dataclass(frozen=True)
class MannKendallResult:
    """Mann-Kendall test of a monotonic trend.

    Attributes
    ----------
    n : int
    s : float
        ``sum_{i<j} sign(x_j - x_i)``.
    var_s : float
        Variance of ``s`` (tie-corrected; Hamed-Rao corrected when asked).
    z : float
        Continuity-corrected normal score.
    p_value : float
        Two-sided.
    tau : float
        Kendall's tau, ``s / (n (n - 1) / 2)``.
    trend : str
        ``"increasing"``, ``"decreasing"`` or ``"no trend"`` at ``alpha``.
    method : str
        ``"original"`` or ``"hamed_rao"``.
    """

    n: int
    s: float
    var_s: float
    z: float
    p_value: float
    tau: float
    trend: str
    method: str


def _hamed_rao_factor(x: np.ndarray) -> float:
    """Hamed and Rao (1998) n / n* from the significant autocorrelations of the ranks.

    The ranks are detrended with Sen's slope first, and only lags whose autocorrelation is
    significant at the 5 percent level (two-sided, ``1.96 / sqrt(n)``) enter the sum.
    """
    n = x.size
    t = np.arange(n, dtype=float)
    slope, _ = sens_slope(x, t)
    detrended = x - slope * t
    ranks = stats.rankdata(detrended)
    r = ranks - ranks.mean()
    denom = float(np.sum(r * r))
    total = 0.0
    bound = 1.96 / math.sqrt(n)
    for k in range(1, n - 1):
        rho = float(np.sum(r[:-k] * r[k:]) / denom)
        if abs(rho) > bound:
            total += (n - k) * (n - k - 1) * (n - k - 2) * rho
    return 1.0 + 2.0 / (n * (n - 1) * (n - 2)) * total


def mann_kendall(
    x: ArrayLike, *, method: str = "original", alpha: float = 0.05
) -> MannKendallResult:
    """Mann-Kendall test for a monotonic trend in a series in time order.

    Parameters
    ----------
    x : array_like
        At least three finite values, in time order.
    method : {"original", "hamed_rao"}, default "original"
        ``"hamed_rao"`` scales ``var_s`` by the Hamed and Rao (1998) factor for serial
        correlation.
    alpha : float, default 0.05
        Level for ``trend``.

    Returns
    -------
    MannKendallResult

    Raises
    ------
    ValueError
        Fewer than three values, a non-finite value, or an unknown ``method``.
    """
    if method not in ("original", "hamed_rao"):
        raise ValueError(f"method must be 'original' or 'hamed_rao', got {method!r}")
    a = _series(x)
    n = a.size
    diff = a[None, :] - a[:, None]
    s = float(np.sum(np.sign(diff[np.triu_indices(n, k=1)])))
    _, counts = np.unique(a, return_counts=True)
    ties = float(np.sum(counts * (counts - 1) * (2 * counts + 5)))
    var_s = (n * (n - 1) * (2 * n + 5) - ties) / 18.0
    if method == "hamed_rao":
        var_s *= max(_hamed_rao_factor(a), 1e-12)
    if s > 0:
        z = (s - 1.0) / math.sqrt(var_s)
    elif s < 0:
        z = (s + 1.0) / math.sqrt(var_s)
    else:
        z = 0.0
    p = float(2.0 * stats.norm.sf(abs(z)))
    trend = "no trend" if p >= alpha else ("increasing" if z > 0 else "decreasing")
    return MannKendallResult(
        n=n,
        s=s,
        var_s=float(var_s),
        z=float(z),
        p_value=p,
        tau=s / (n * (n - 1) / 2.0),
        trend=trend,
        method=method,
    )


def sens_slope(x: ArrayLike, t: Optional[ArrayLike] = None) -> Tuple[float, float]:
    """Sen's slope (median pairwise slope) and the intercept ``median(x - slope t)``.

    Parameters
    ----------
    x : array_like
    t : array_like, optional
        Times (e.g. water years); ``0, 1, ...`` when omitted.

    Returns
    -------
    (float, float)
        Slope per unit of ``t``, and intercept.
    """
    a = _series(x, min_n=2)
    times = np.arange(a.size, dtype=float) if t is None else np.asarray(t, dtype=float)
    if times.shape != a.shape:
        raise ValueError("x and t differ in length")
    i, j = np.triu_indices(a.size, k=1)
    dt = times[j] - times[i]
    keep = dt != 0
    slope = float(np.median((a[j] - a[i])[keep] / dt[keep]))
    return slope, float(np.median(a - slope * times))


def spearman_trend(x: ArrayLike, t: Optional[ArrayLike] = None) -> Tuple[float, float]:
    """Spearman's rho of ``x`` against time and its two-sided p-value (scipy)."""
    a = _series(x)
    times = np.arange(a.size, dtype=float) if t is None else np.asarray(t, dtype=float)
    res = stats.spearmanr(times, a)
    return float(res.statistic), float(res.pvalue)


@dataclass(frozen=True)
class PettittResult:
    """Pettitt change-point test.

    Attributes
    ----------
    index : int
        Position of the first value after the change (the second segment starts here).
    year : int or None
        ``years[index]`` when years were given.
    k : float
        ``max_t |U_t|``.
    p_value : float
        Approximate two-sided p-value, ``2 exp(-6 k^2 / (n^3 + n^2))``, capped at 1.
    shift : str
        ``"up"`` or ``"down"``: the second segment's median against the first's.
    """

    index: int
    year: Optional[int]
    k: float
    p_value: float
    shift: str


def pettitt(x: ArrayLike, years: Optional[Iterable[int]] = None) -> PettittResult:
    """Pettitt (1979) test for a single change point in a series in time order.

    Parameters
    ----------
    x : array_like
        At least three finite values, in time order.
    years : iterable of int, optional
        The years of ``x``, to report the change point as a year.

    Returns
    -------
    PettittResult
    """
    a = _series(x)
    n = a.size
    sign = np.sign(a[None, :] - a[:, None])  # sign[i, j] = sign(x_j - x_i)
    # U_t = sum_{i < t} sum_{j >= t} sign(x_i - x_j), for split after position t-1.
    u = np.array([-float(np.sum(sign[:t, t:])) for t in range(1, n)])
    t_star = int(np.argmax(np.abs(u)))
    k = float(abs(u[t_star]))
    index = t_star + 1
    p = float(min(1.0, 2.0 * math.exp(-6.0 * k**2 / (n**3 + n**2))))
    year_list = list(years) if years is not None else None
    if year_list is not None and len(year_list) != n:
        raise ValueError("x and years differ in length")
    shift = "up" if np.median(a[index:]) > np.median(a[:index]) else "down"
    return PettittResult(
        index=index,
        year=int(year_list[index]) if year_list is not None else None,
        k=k,
        p_value=p,
        shift=shift,
    )


def stationarity_screen(
    values: pd.Series, *, log: bool = True, alpha: float = 0.05
) -> Dict[str, Any]:
    """Trend and change-point tests on one annual series, with plain-language flags.

    Parameters
    ----------
    values : pd.Series
        Annual values (e.g. peak flows) indexed by water year. Missing values are dropped.
    log : bool, default True
        Test ``log10(values)`` (flows); values must then be positive.
    alpha : float, default 0.05

    Returns
    -------
    dict
        ``n``, ``years`` (first, last), ``mann_kendall`` and ``mann_kendall_hamed_rao``
        (``s``, ``z``, ``p_value``, ``tau``, ``trend``), ``sens_slope_log10_per_year``
        (or ``sens_slope_per_year`` when ``log=False``), ``spearman`` (``rho``,
        ``p_value``), ``pettitt`` (``year``, ``k``, ``p_value``, ``shift``) and ``flags``:
        a list of messages, empty when nothing is significant at ``alpha``.

    Raises
    ------
    ValueError
        Fewer than ten values, or a non-positive value with ``log=True``.
    """
    s = values.dropna().sort_index()
    if len(s) < 10:
        raise ValueError(f"need at least 10 years to screen, got {len(s)}")
    x = np.asarray(s, dtype=float)
    if log:
        if np.any(x <= 0):
            raise ValueError("log=True needs positive values (zero flows: use log=False)")
        x = np.log10(x)
    years = np.asarray(s.index, dtype=int)
    mk = mann_kendall(x, alpha=alpha)
    mk_hr = mann_kendall(x, method="hamed_rao", alpha=alpha)
    slope, _ = sens_slope(x, years)
    rho, p_rho = spearman_trend(x, years)
    pt = pettitt(x, years=years)
    flags: List[str] = []
    if mk_hr.trend != "no trend":
        flags.append(
            f"{mk_hr.trend} trend (Mann-Kendall with Hamed-Rao correction, p = {mk_hr.p_value:.3f})"
        )
    elif mk.trend != "no trend":
        flags.append(
            f"{mk.trend} trend in the plain Mann-Kendall test (p = {mk.p_value:.3f}), not "
            "significant once serial correlation is allowed for"
        )
    if pt.p_value < alpha:
        flags.append(
            f"change point in {pt.year} (Pettitt, p = {pt.p_value:.3f}, shift {pt.shift}); check "
            "for regulation, urbanization, a gage move or a rating change about then"
        )
    for f in flags:
        logger.info("stationarity screen: %s", f)
    slope_key = "sens_slope_log10_per_year" if log else "sens_slope_per_year"
    return {
        "n": int(len(s)),
        "years": (int(years[0]), int(years[-1])),
        "mann_kendall": {k: getattr(mk, k) for k in ("s", "z", "p_value", "tau", "trend")},
        "mann_kendall_hamed_rao": {
            k: getattr(mk_hr, k) for k in ("s", "z", "p_value", "tau", "trend")
        },
        slope_key: slope,
        "spearman": {"rho": rho, "p_value": p_rho},
        "pettitt": {"year": pt.year, "k": pt.k, "p_value": pt.p_value, "shift": pt.shift},
        "flags": flags,
    }


#: Months whose annual peaks are rain or rain-on-snow driven in snowmelt regions.
COOL_SEASON = (10, 11, 12, 1, 2, 3)


def mixed_population_screen(
    peak_dates: Iterable[Any],
    *,
    cool_months: Sequence[int] = COOL_SEASON,
    mixed_band: Tuple[float, float] = (0.2, 0.8),
) -> Dict[str, Any]:
    """How annual peaks split between the cool and warm seasons (warn only).

    Parameters
    ----------
    peak_dates : iterable of date-like
        The date of each annual peak.
    cool_months : sequence of int, default October-March
        Months counted as the cool (rain, rain-on-snow) season.
    mixed_band : (float, float), default (0.2, 0.8)
        The record is flagged ``mixed`` when the cool-season share lies inside this band.

    Returns
    -------
    dict
        ``n``, ``cool_season_share``, ``mean_day_of_year`` (circular mean),
        ``concentration`` (mean resultant length, 0 for dates spread evenly through the
        year, 1 for a single date; Burn 1997), ``mixed`` and ``message``.
    """
    dates = pd.to_datetime(pd.Series(list(peak_dates))).dropna()
    if dates.empty:
        raise ValueError("no peak dates")
    n = int(len(dates))
    cool = float(dates.dt.month.isin(list(cool_months)).mean())
    doy = dates.dt.dayofyear.to_numpy(dtype=float)
    days_in_year = np.where(dates.dt.is_leap_year.to_numpy(), 366.0, 365.0)
    theta = 2.0 * np.pi * doy / days_in_year
    c, s = float(np.mean(np.cos(theta))), float(np.mean(np.sin(theta)))
    concentration = math.hypot(c, s)
    mean_doy = (math.degrees(math.atan2(s, c)) % 360.0) / 360.0 * 365.25
    mixed = mixed_band[0] < cool < mixed_band[1]
    message = (
        f"possibly mixed flood populations: {cool:.0%} of {n} annual peaks fall in the cool "
        "season and the rest in the warm season. Bulletin 17C section 9.3: consider whether "
        "one distribution fits both."
        if mixed
        else f"{cool:.0%} of {n} annual peaks fall in the cool season"
    )
    if mixed:
        logger.warning(message)
    return {
        "n": n,
        "cool_season_share": cool,
        "mean_day_of_year": mean_doy,
        "concentration": concentration,
        "mixed": bool(mixed),
        "message": message,
    }
