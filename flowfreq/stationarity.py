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
- :func:`cusum`: the CUSUM change-point statistic, with a permutation p-value.
- :func:`lombard`: Lombard's rank tests for an abrupt or a smooth change, in level
  (Wilcoxon score) or in variability (Mood score), as USACE's Nonstationarity Detection
  Tool runs them.
- :func:`field_significance`: Benjamini-Hochberg FDR and Walker's test across sites.
- :func:`field_significance_bootstrap`: the same question answered by resampling years,
  which keeps the correlation between nearby sites (and, in blocks, within a site).
- :func:`attribute_change_point`: candidate non-climate causes near a change point, from
  the peak qualification codes (regulation, urbanization) and caller-supplied events.
- :func:`basin_change_events`: dated events for that overlay from the National Inventory
  of Dams and the NLCD impervious epochs StreamStats reports.
- :func:`stationarity_screen`: all of them on one series (in log10 for flows), with
  plain-language flags.
- :func:`mixed_population_screen`: how annual peaks split between the cool season
  (rain, rain-on-snow) and the warm season (snowmelt), and the concentration of their
  dates (Burn, 1997). A split record is a warning that one LP3 distribution may be
  fitted to two flood populations (Bulletin 17C section 9.3); it is not a test.

- :func:`pot_counts` and :func:`pot_count_trend`: annual counts of independent
  peaks over a threshold in a daily series, and a Poisson-regression trend in them.
- :func:`stationarity_report`: the screen as a one-page figure and a Markdown summary.
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


def _pettitt_u(a: np.ndarray) -> np.ndarray:
    """Pettitt's ``U_t`` for every split, ``t = 1 .. n-1`` (split after position t-1)."""
    sign = np.sign(a[None, :] - a[:, None])  # sign[i, j] = sign(x_j - x_i)
    # U_t = sum_{i < t} sum_{j >= t} sign(x_i - x_j)
    return np.array([-float(np.sum(sign[:t, t:])) for t in range(1, a.size)])


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
    u = _pettitt_u(a)
    t_star = int(np.argmax(np.abs(u)))
    k = float(abs(u[t_star]))
    index = t_star + 1
    p = float(min(1.0, 2.0 * math.exp(-6.0 * k**2 / (n**3 + n**2))))
    year_list = list(years) if years is not None else None
    if year_list is not None and len(year_list) != n:
        raise ValueError("x and years differ in length")
    # U_t sums sign(x_i - x_j) with i before the change: negative when later values are
    # larger.
    shift = "up" if u[t_star] < 0 else "down"
    return PettittResult(
        index=index,
        year=int(year_list[index]) if year_list is not None else None,
        k=k,
        p_value=p,
        shift=shift,
    )


#: Peak-code evidence the attribution overlay looks for, with its label.
_CODE_EVIDENCE = (
    (frozenset({"5", "6"}), "regulation or diversion (peak code 5/6)"),
    (frozenset({"C"}), "urbanization, mining, agricultural or channel change (peak code C)"),
)


def attribute_change_point(
    change_year: int,
    water_years: Sequence[int],
    peak_codes: Sequence[Any],
    *,
    events: Iterable[Tuple[str, int]] = (),
    window: int = 5,
) -> Dict[str, Any]:
    """Candidate non-climate causes of a change point (roadmap section 6.1).

    A change point that lines up with the start of regulation, urbanization, a dam, a gage
    move or a rating change is not a climate signal. This lists, within ``window`` years
    of ``change_year``: the first year the record carries regulation evidence (peak codes
    5 or 6) and alteration evidence (code C), and every caller-supplied event (a dam's
    completion year from NID, a gage relocation, a datum or rating change).

    Parameters
    ----------
    change_year : int
        First year of the second segment (e.g. :func:`pettitt`'s ``year``).
    water_years, peak_codes : sequence
        One entry per annual peak; codes as NWIS reports them (``"6"``, ``"6C"``, ...).
    events : iterable of (str, int), optional
        ``(label, year)`` of other known changes.
    window : int, default 5

    Returns
    -------
    dict
        ``change_year``, ``candidates`` (``label``, ``year``, ``offset`` = year minus
        ``change_year``, ``source``; nearest first), ``explained`` (any candidate) and
        ``message``.

    Raises
    ------
    ValueError
        ``water_years`` and ``peak_codes`` differ in length.
    """
    from flowfreq.peak_codes import parse_codes

    years = [int(y) for y in water_years]
    if len(years) != len(peak_codes):
        raise ValueError("water_years and peak_codes differ in length")
    parsed = [parse_codes(c) for c in peak_codes]
    found: List[Tuple[str, int, str]] = []
    for codes, label in _CODE_EVIDENCE:
        first = [y for y, pc in sorted(zip(years, parsed)) if pc & codes]
        if first:
            found.append((label, first[0], "peak codes"))
    found += [(str(label), int(year), "event") for label, year in events]
    near = sorted(
        (
            (abs(year - int(change_year)), year, label, src)
            for label, year, src in found
            if abs(year - int(change_year)) <= window
        )
    )
    candidates: List[Dict[str, Any]] = [
        {"label": label, "year": year, "offset": year - int(change_year), "source": src}
        for _, year, label, src in near
    ]
    if candidates:
        message = "; ".join(f"{c['label']} in {c['year']}" for c in candidates)
    else:
        message = f"no recorded regulation, alteration or event within {window} years"
    return {
        "change_year": int(change_year),
        "candidates": candidates,
        "explained": bool(candidates),
        "message": message,
    }


def stationarity_screen(
    values: pd.Series,
    *,
    log: bool = True,
    alpha: float = 0.05,
    peak_codes: Optional[Sequence[Any]] = None,
    events: Iterable[Tuple[str, int]] = (),
    n_permutations: int = 999,
    seed: Optional[int] = 0,
) -> Dict[str, Any]:
    """Trend and change-point tests on one annual series, with plain-language flags.

    Parameters
    ----------
    values : pd.Series
        Annual values (e.g. peak flows) indexed by water year. Missing values are dropped.
    log : bool, default True
        Test ``log10(values)`` (flows); values must then be positive.
    alpha : float, default 0.05
    peak_codes : sequence, optional
        Peak qualification codes aligned with ``values`` (before missing values are
        dropped). With them, or with ``events``, a significant change point is passed to
        :func:`attribute_change_point` and its flag names the candidate causes.
    events : iterable of (str, int), optional
        Other known changes, as :func:`attribute_change_point` takes them (for example
        :func:`basin_change_events`'s ``events``).
    n_permutations : int, default 999
        Permutations for the Lombard p-values.
    seed : int, optional
        Seed of those permutations; the default 0 makes the screen reproducible.

    Returns
    -------
    dict
        ``n``, ``years`` (first, last), ``mann_kendall`` and ``mann_kendall_hamed_rao``
        (``s``, ``z``, ``p_value``, ``tau``, ``trend``), ``sens_slope_log10_per_year``
        (or ``sens_slope_per_year`` when ``log=False``), ``spearman`` (``rho``,
        ``p_value``), ``pettitt`` (``year``, ``k``, ``p_value``, ``shift``),
        ``lombard`` (``wilcoxon_abrupt``, ``wilcoxon_smooth``, ``mood_abrupt``,
        ``mood_smooth``, each ``statistic``, ``p_value``, ``start_year``, ``end_year``),
        ``attribution`` and ``flags``: a list of messages, empty when nothing is
        significant at ``alpha``.

    Raises
    ------
    ValueError
        Fewer than ten values, or a non-positive value with ``log=True``.
    """
    if peak_codes is not None and len(peak_codes) != len(values):
        raise ValueError("peak_codes must align with values")
    code_by_year = (
        dict(zip(np.asarray(values.index, dtype=int).tolist(), list(peak_codes)))
        if peak_codes is not None
        else {}
    )
    events = list(events)
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
    lomb = {
        f"{score}_{kind}": lombard(
            x, years, kind=kind, score=score, n_permutations=n_permutations, seed=seed
        )
        for score in ("wilcoxon", "mood")
        for kind in ("abrupt", "smooth")
    }
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
    attribution: Optional[Dict[str, Any]] = None
    if pt.p_value < alpha:
        if (code_by_year or events) and pt.year is not None:
            yrs = sorted(code_by_year)
            attribution = attribute_change_point(
                int(pt.year), yrs, [code_by_year[y] for y in yrs], events=events
            )
            why = (
                f"candidate causes: {attribution['message']}"
                if attribution["explained"]
                else f"{attribution['message']}; check for a gage move or a rating change"
            )
        else:
            why = "check for regulation, urbanization, a gage move or a rating change about then"
        flags.append(
            f"change point in {pt.year} (Pettitt, p = {pt.p_value:.3f}, shift {pt.shift}); {why}"
        )
    gradual = lomb["wilcoxon_smooth"]
    if gradual.p_value < alpha and pt.p_value >= alpha:
        flags.append(
            f"gradual change in level between {gradual.start_year} and {gradual.end_year} "
            f"(Lombard smooth, p = {gradual.p_value:.3f})"
        )
    spread = min((lomb["mood_abrupt"], lomb["mood_smooth"]), key=lambda r: r.p_value)
    if spread.p_value < alpha:
        when = (
            f"about {spread.start_year}"
            if spread.kind == "abrupt"
            else f"between {spread.start_year} and {spread.end_year}"
        )
        flags.append(
            f"change in variability {when} (Lombard Mood {spread.kind}, "
            f"p = {spread.p_value:.3f}); one LP3 fit may not describe the whole record"
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
        "lombard": {
            key: {
                "statistic": r.statistic,
                "p_value": r.p_value,
                "start_year": r.start_year,
                "end_year": r.end_year,
            }
            for key, r in lomb.items()
        },
        "attribution": attribution,
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


@dataclass(frozen=True)
class CusumResult:
    """CUSUM change-point statistic.

    Attributes
    ----------
    index : int
        Position of the first value after the change.
    year : int or None
        ``years[index]`` when years were given.
    max_abs : float
        ``max_k |sum_{i<=k} (x_i - mean(x))|``.
    p_value : float
        Permutation p-value, ``(1 + #{perm >= observed}) / (1 + n_permutations)``.
    """

    index: int
    year: Optional[int]
    max_abs: float
    p_value: float


def cusum(
    x: ArrayLike,
    years: Optional[Iterable[int]] = None,
    *,
    n_permutations: int = 999,
    seed: Optional[int] = None,
) -> CusumResult:
    """CUSUM test for a shift in the mean, with a permutation p-value.

    The cumulative sum of deviations from the mean peaks (in magnitude) at the last value
    before a shift. Its significance comes from shuffling the series, which keeps the
    values and breaks their order.

    Parameters
    ----------
    x : array_like
        At least three finite values, in time order.
    years : iterable of int, optional
    n_permutations : int, default 999
    seed : int, optional

    Returns
    -------
    CusumResult
    """
    a = _series(x)
    if n_permutations < 1:
        raise ValueError(f"n_permutations must be >= 1, got {n_permutations}")
    dev = np.cumsum(a - a.mean())
    k = int(np.argmax(np.abs(dev[:-1])))
    observed = float(abs(dev[k]))
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(n_permutations):
        perm = rng.permutation(a)
        if float(np.max(np.abs(np.cumsum(perm - perm.mean())[:-1]))) >= observed:
            exceed += 1
    year_list = list(years) if years is not None else None
    if year_list is not None and len(year_list) != a.size:
        raise ValueError("x and years differ in length")
    index = k + 1
    return CusumResult(
        index=index,
        year=int(year_list[index]) if year_list is not None else None,
        max_abs=observed,
        p_value=(1.0 + exceed) / (1.0 + n_permutations),
    )


def field_significance(p_values: ArrayLike, *, alpha: float = 0.05) -> Dict[str, Any]:
    """Field significance of single-site tests across a region.

    In a statewide run, about ``alpha`` of the sites test significant by chance alone.
    This reports the Benjamini and Hochberg (1995) false-discovery-rate rejections, and
    Walker's test that the smallest p-value is significant for ``m`` sites,
    ``p_walker = 1 - (1 - p_min)^m``. Neither accounts for spatial correlation between
    sites, which makes both somewhat liberal.

    Parameters
    ----------
    p_values : array_like
        One p-value per site, each in [0, 1].
    alpha : float, default 0.05

    Returns
    -------
    dict
        ``n_sites``, ``n_locally_significant`` (``p < alpha``), ``fdr_reject`` (list of
        bool in input order), ``fdr_threshold`` (the largest rejected p-value, NaN when
        none), ``walker_p`` and ``walker_significant``.
    """
    p = np.asarray(p_values, dtype=float).ravel()
    if p.size == 0 or np.any(~np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("p_values must be a non-empty set of values in [0, 1]")
    m = p.size
    order = np.argsort(p)
    ranked = p[order]
    below = ranked <= alpha * np.arange(1, m + 1) / m
    reject = np.zeros(m, dtype=bool)
    threshold = float("nan")
    if below.any():
        k = int(np.max(np.flatnonzero(below)))
        reject[order[: k + 1]] = True
        threshold = float(ranked[k])
    walker = float(1.0 - (1.0 - float(ranked[0])) ** m)
    return {
        "n_sites": int(m),
        "n_locally_significant": int(np.sum(p < alpha)),
        "fdr_reject": [bool(r) for r in reject],
        "fdr_threshold": threshold,
        "walker_p": walker,
        "walker_significant": bool(walker < alpha),
    }


# ------------------------------------------------------------------ Lombard --

#: Score functions of Lombard's rank tests: Wilcoxon for a change in location, Mood for a
#: change in scale (variability).
_LOMBARD_SCORES = {
    "wilcoxon": lambda u: 2.0 * u - 1.0,
    "mood": lambda u: (u - 0.5) ** 2,
}

#: Asymptotic critical values of the abrupt statistic, the Cramer-von Mises integral of a
#: squared Brownian bridge (Anderson and Darling 1952), by level.
LOMBARD_ABRUPT_CRITICAL = {0.10: 0.347, 0.05: 0.461, 0.01: 0.743}


@dataclass(frozen=True)
class LombardResult:
    """Lombard's rank test for an abrupt or a smooth (gradual) change.

    Attributes
    ----------
    kind : str
        ``"abrupt"`` or ``"smooth"``.
    score : str
        ``"wilcoxon"`` (location) or ``"mood"`` (scale).
    statistic : float
        ``S_n`` as Rouge and others (2013) normalize it.
    p_value : float
        Permutation p-value, ``(1 + #{perm >= observed}) / (1 + n_permutations)``.
    start_index, end_index : int
        Abrupt: both are the first value after the change. Smooth: the first value off
        the old level and the first at the new one. Either is the best-fitting candidate
        (:func:`_lombard_locate`), a pointer rather than an estimate with a stated
        uncertainty.
    start_year, end_year : int or None
        The years at those positions, when years were given.
    """

    kind: str
    score: str
    statistic: float
    p_value: float
    start_index: int
    end_index: int
    start_year: Optional[int]
    end_year: Optional[int]


def _lombard_scores(a: np.ndarray, score: str) -> np.ndarray:
    phi = _LOMBARD_SCORES[score](stats.rankdata(a) / (a.size + 1.0))
    sd = float(np.std(phi))
    if sd == 0.0:
        raise ValueError("the series is constant; Lombard's scores are undefined")
    return np.asarray((phi - phi.mean()) / sd)


def _lombard_stat(z: np.ndarray, kind: str) -> np.ndarray:
    """Statistic for each row of ``z`` (standardized scores, one series per row)."""
    n = z.shape[-1]
    partial = np.cumsum(z, axis=-1)  # L_k, k = 1..n
    if kind == "abrupt":
        return np.asarray(np.sum(partial**2, axis=-1) / n**2)
    s = np.cumsum(partial, axis=-1)  # S_tau = sum_{j <= tau} L_j
    # sum_{1 <= a < b <= n} (S_b - S_a)^2 = n sum S^2 - (sum S)^2
    pairs = n * np.sum(s**2, axis=-1) - np.sum(s, axis=-1) ** 2
    return np.asarray(pairs / float(n) ** 5)


def _lombard_locate(z: np.ndarray, kind: str) -> Tuple[int, int]:
    """Where the change is: the candidate whose contrast is largest relative to its null
    standard deviation (the score statistic of that one candidate).

    A candidate's contrast is ``sum_i w_i Z_i``: for a step after position ``k``,
    ``w_i = 1`` for ``i <= k``; for a ramp over ``(tau1, tau2]``,
    ``w_i = #{j in (tau1, tau2] : j >= i}``, the partial-sum weights of
    ``C_{tau1 tau2}``. Under shuffling its variance is proportional to
    ``sum (w - mean w)^2``. The raw contrast, unnormalized, favors the middle of the
    record for a step and the whole record for a ramp.
    """
    n = z.size
    if kind == "abrupt":
        k = np.arange(1, n)
        contrast = np.cumsum(z)[:-1]
        score = contrast**2 / (k * (n - k) / n)
        k_best = int(np.argmax(score))
        return k_best + 1, k_best + 1
    i = np.arange(1, n + 1)
    best_val, best = -1.0, (0, n - 1)
    for tau1 in range(0, n - 1):
        tau2 = np.arange(tau1 + 1, n + 1)[:, None]
        w = np.clip(tau2 - np.maximum(i[None, :] - 1, tau1), 0, None).astype(float)
        var = np.sum((w - w.mean(axis=1, keepdims=True)) ** 2, axis=1)
        ok = var > 1e-12
        if not ok.any():
            continue
        score = np.where(ok, (w @ z) ** 2 / np.where(ok, var, 1.0), -1.0)
        j = int(np.argmax(score))
        if score[j] > best_val:
            best_val, best = float(score[j]), (tau1, int(tau2[j, 0]) - 1)
    return best


def lombard(
    x: ArrayLike,
    years: Optional[Iterable[int]] = None,
    *,
    kind: str = "abrupt",
    score: str = "wilcoxon",
    n_permutations: int = 999,
    seed: Optional[int] = None,
) -> LombardResult:
    """Lombard's (1987) rank test for a change in location or scale.

    The set USACE's Nonstationarity Detection Tool runs (ECB 2018-14): ``kind="abrupt"``
    tests for a step and ``kind="smooth"`` for a gradual change between two unknown times;
    the Wilcoxon score detects a change in level, the Mood score a change in variability.
    With standardized rank scores ``Z_i`` and partial sums ``L_k = sum_{i<=k} Z_i``,

    - abrupt: ``S_n = n^-2 sum_k L_k^2``;
    - smooth: ``S_n = n^-5 sum_{tau1 < tau2} (sum_{tau1 < j <= tau2} L_j)^2``

    (Rouge, Ge and Cai 2013). The abrupt statistic converges to the Cramer-von Mises
    distribution (:data:`LOMBARD_ABRUPT_CRITICAL`), and the smooth one has a 5 percent
    point near 0.040 (simulated; both are checked in the tests). The p-value here comes
    from shuffling the series, which needs no asymptotics and is the same for both kinds.

    Parameters
    ----------
    x : array_like
        At least five finite values, in time order.
    years : iterable of int, optional
    kind : {"abrupt", "smooth"}, default "abrupt"
    score : {"wilcoxon", "mood"}, default "wilcoxon"
    n_permutations : int, default 999
    seed : int, optional

    Returns
    -------
    LombardResult

    Raises
    ------
    ValueError
        Fewer than five values, a non-finite or constant series, an unknown ``kind`` or
        ``score``, or ``years`` of another length.
    """
    if kind not in ("abrupt", "smooth"):
        raise ValueError(f"kind must be 'abrupt' or 'smooth', got {kind!r}")
    if score not in _LOMBARD_SCORES:
        raise ValueError(f"score must be one of {sorted(_LOMBARD_SCORES)}, got {score!r}")
    if n_permutations < 1:
        raise ValueError(f"n_permutations must be >= 1, got {n_permutations}")
    a = _series(x, min_n=5)
    year_list = list(years) if years is not None else None
    if year_list is not None and len(year_list) != a.size:
        raise ValueError("x and years differ in length")
    z = _lombard_scores(a, score)
    observed = float(_lombard_stat(z, kind))
    rng = np.random.default_rng(seed)
    exceed = 0
    for start in range(0, n_permutations, 500):
        batch = rng.permuted(np.tile(z, (min(500, n_permutations - start), 1)), axis=1)
        exceed += int(np.sum(_lombard_stat(batch, kind) >= observed * (1.0 - 1e-12)))
    start_index, end_index = _lombard_locate(z, kind)
    return LombardResult(
        kind=kind,
        score=score,
        statistic=observed,
        p_value=(1.0 + exceed) / (1.0 + n_permutations),
        start_index=start_index,
        end_index=end_index,
        start_year=int(year_list[start_index]) if year_list is not None else None,
        end_year=int(year_list[end_index]) if year_list is not None else None,
    )


# ------------------------------------------- bootstrap field significance --


def _mk_p(a: np.ndarray) -> float:
    """Two-sided plain Mann-Kendall p-value, as :func:`mann_kendall` computes it."""
    n = a.size
    i, j = np.triu_indices(n, k=1)
    s = float(np.sum(np.sign(a[j] - a[i])))
    _, counts = np.unique(a, return_counts=True)
    ties = float(np.sum(counts * (counts - 1) * (2 * counts + 5)))
    var_s = (n * (n - 1) * (2 * n + 5) - ties) / 18.0
    if var_s <= 0.0 or s == 0.0:
        return 1.0
    z = (abs(s) - 1.0) / math.sqrt(var_s)
    return float(math.erfc(z / math.sqrt(2.0)))


def _site_p_values(x: np.ndarray, min_years: int) -> np.ndarray:
    out = np.full(x.shape[1], np.nan)
    for c in range(x.shape[1]):
        col = x[:, c]
        col = col[np.isfinite(col)]
        if col.size >= min_years:
            out[c] = _mk_p(col)
    return out


def field_significance_bootstrap(
    values: pd.DataFrame,
    *,
    alpha: float = 0.05,
    n_bootstrap: int = 1000,
    block_length: int = 1,
    min_years: int = 10,
    seed: Optional[int] = None,
) -> Dict[str, Any]:
    """Field significance of Mann-Kendall trends across correlated sites, by bootstrap.

    :func:`field_significance` treats sites as independent, which makes it liberal when
    nearby gages see the same storms. This resamples whole years (Douglas, Vogel and
    Kroll 2000): every site takes the same years, so the cross-site correlation is kept
    while the time order, and with it any trend, is broken. The count of locally
    significant sites in the record is compared with its distribution over the
    resamples. ``block_length > 1`` resamples moving blocks of consecutive years, which
    also keeps serial correlation (Khaliq and others 2009).

    Parameters
    ----------
    values : pandas.DataFrame
        One row per year (in time order), one column per site; missing values allowed.
        Mann-Kendall is rank-based, so flows need no log transform.
    alpha : float, default 0.05
        Level of the local tests and of the field test.
    n_bootstrap : int, default 1000
    block_length : int, default 1
    min_years : int, default 10
        Sites with fewer values are left out of the record's count; in a resample, a
        site with fewer values than that is not counted.
    seed : int, optional

    Returns
    -------
    dict
        ``n_sites``, ``n_years``, ``site_p_values`` (site to p), ``n_significant``,
        ``critical_count`` (the ``1 - alpha`` quantile of the resampled counts),
        ``p_value``, ``field_significant`` and ``block_length``.

    Raises
    ------
    ValueError
        No site with ``min_years`` values, or a bad ``n_bootstrap``/``block_length``.
    """
    if n_bootstrap < 1:
        raise ValueError(f"n_bootstrap must be >= 1, got {n_bootstrap}")
    x = values.to_numpy(dtype=float)
    n_years = x.shape[0]
    if not 1 <= block_length <= n_years:
        raise ValueError(f"block_length must be in [1, {n_years}], got {block_length}")
    min_years = max(min_years, 3)
    keep = np.sum(np.isfinite(x), axis=0) >= min_years
    if not keep.any():
        raise ValueError(f"no site has {min_years} or more values")
    x = x[:, keep]
    sites = [str(c) for c, k in zip(values.columns, keep) if k]
    p_obs = _site_p_values(x, min_years)
    observed = int(np.sum(p_obs < alpha))
    rng = np.random.default_rng(seed)
    n_blocks = -(-n_years // block_length)
    counts = np.empty(n_bootstrap, dtype=int)
    for b in range(n_bootstrap):
        starts = rng.integers(0, n_years - block_length + 1, size=n_blocks)
        rows = (starts[:, None] + np.arange(block_length)[None, :]).ravel()[:n_years]
        counts[b] = int(np.sum(_site_p_values(x[rows], min_years) < alpha))
    p_value = float((1 + np.sum(counts >= observed)) / (1 + n_bootstrap))
    return {
        "n_sites": len(sites),
        "n_years": int(n_years),
        "site_p_values": dict(zip(sites, (float(p) for p in p_obs))),
        "n_significant": observed,
        "critical_count": float(np.quantile(counts, 1.0 - alpha)),
        "p_value": p_value,
        "field_significant": bool(p_value < alpha),
        "block_length": int(block_length),
    }


# ------------------------------------------------- basin change events --


def basin_change_events(
    *,
    geometry: Optional[Dict[str, Any]] = None,
    nid: Any = None,
    impervious: Optional[Dict[int, float]] = None,
    region: Optional[str] = None,
    min_storage_share: float = 0.1,
    urban_impervious_pct: Optional[float] = None,
) -> Dict[str, Any]:
    """Dated basin changes for :func:`attribute_change_point`, from NID and NLCD.

    - Each NID dam inside the watershed whose storage is at least ``min_storage_share``
      of the basin's total NID storage, at its ``year_completed``.
    - The year half of today's basin NID storage was in place.
    - The NLCD epoch in which basin impervious cover first reached
      ``urban_impervious_pct`` (the regulation screen's urban threshold, 5 percent).
      NLCD starts in 2001, so a basin already above it then gets a note, not a date.

    Parameters
    ----------
    geometry : dict, optional
        GeoJSON ``Polygon``/``MultiPolygon`` of the watershed (e.g. the ``geometry`` of
        ``WatershedCharacteristics.polygon_geojson``). Needed, with ``nid``, for dams.
    nid : flowfreq.regulation.NIDInventory, optional
        From :func:`flowfreq.regulation.load_nid`.
    impervious : dict of int to float, optional
        Year to percent impervious, from :func:`flowfreq.regulation.impervious_series`.
    region : str, optional
        StreamStats region; a region in
        :data:`flowfreq.regulation.IMPERVIOUS_UNVERIFIED_REGIONS` has its impervious
        series skipped, as the regulation screen skips it.
    min_storage_share : float, default 0.1
    urban_impervious_pct : float, optional
        Defaults to :data:`flowfreq.regulation.URBAN_IMPERVIOUS_PCT`.

    Returns
    -------
    dict
        ``events`` (list of ``(label, year)``, ready for ``events=``), ``notes`` (what
        could not be dated) and ``n_dams`` (dams in the basin, or None without NID).

    Raises
    ------
    ValueError
        Only one of ``geometry`` and ``nid`` given.
    """
    from flowfreq import regulation

    threshold = (
        regulation.URBAN_IMPERVIOUS_PCT if urban_impervious_pct is None else urban_impervious_pct
    )
    events: List[Tuple[str, int]] = []
    notes: List[str] = []
    n_dams: Optional[int] = None
    if (geometry is None) != (nid is None):
        raise ValueError("dam events need both geometry and nid")
    if geometry is not None and nid is not None:
        dams = regulation.nid_dams_in_basin(geometry, nid)
        n_dams = len(dams)
        total = float(dams["nid_storage_acft"].sum())
        dated = dams[dams["year_completed"].notna() & (dams["year_completed"] > 0)]
        if len(dated) < len(dams):
            notes.append(
                f"{len(dams) - len(dated)} of {len(dams)} NID dams have no completion year"
            )
        if total > 0:
            for _, d in dated.iterrows():
                share = float(d["nid_storage_acft"]) / total
                if share >= min_storage_share:
                    events.append(
                        (
                            f"dam {d['nid_id']} completed ({d['nid_storage_acft']:,.0f} acre-ft, "
                            f"{share:.0%} of basin NID storage)",
                            int(d["year_completed"]),
                        )
                    )
            by_year = dated.groupby(dated["year_completed"].astype(int))["nid_storage_acft"].sum()
            cumulative = by_year.cumsum()
            reached = cumulative[cumulative >= 0.5 * total]
            if len(reached):
                events.append(("half of today's basin NID storage in place", int(reached.index[0])))
            else:
                notes.append("undated dams hold more than half of the basin NID storage")
    if impervious:
        if region is not None and region.upper() in regulation.IMPERVIOUS_UNVERIFIED_REGIONS:
            notes.append(f"impervious cover not used in {region.upper()} (unverified there)")
        else:
            yrs = sorted(impervious)
            if impervious[yrs[0]] >= threshold:
                notes.append(
                    f"impervious cover was already {impervious[yrs[0]]:.1f}% in NLCD {yrs[0]}, "
                    "the first land-cover epoch"
                )
            else:
                for prev, cur in zip(yrs, yrs[1:]):
                    if impervious[cur] >= threshold:
                        events.append(
                            (
                                f"basin impervious cover reached {threshold:g}% "
                                f"(NLCD {prev} to {cur})",
                                int(cur),
                            )
                        )
                        break
    return {"events": events, "notes": notes, "n_dams": n_dams}


# ------------------------------------------------ peaks over a threshold --


@dataclass(frozen=True)
class PotCounts:
    """Independent peaks over a threshold in a daily series, counted by water year.

    Attributes
    ----------
    threshold : float
    counts : pandas.Series
        Events per complete water year (zeros included), indexed by water year.
    events : pandas.DataFrame
        One row per event: ``peak_date``, ``peak``, ``water_year``.
    min_separation_days : int
    quantile : float
        Fraction of the record's days at or below the threshold. Below
        :data:`POT_LOW_QUANTILE` the "events" are ordinary flows rather than floods.
    """

    threshold: float
    counts: pd.Series
    events: pd.DataFrame
    min_separation_days: int
    quantile: float = float("nan")


#: A threshold exceeded on more than a quarter of the days is not a flood threshold.
POT_LOW_QUANTILE = 0.75


def _water_year(index: pd.DatetimeIndex) -> np.ndarray:
    return np.asarray(index.year + (index.month >= 10).astype(int), dtype=int)


def _pot_events(s: pd.Series, threshold: float, min_sep: int) -> pd.DataFrame:
    above = s[s > threshold]
    if above.empty:
        return pd.DataFrame({"peak_date": [], "peak": [], "water_year": []})
    days = (above.index - above.index[0]).days.to_numpy()
    # A new event starts when the gap since the previous exceedance day is at least
    # min_sep days, so exceedances closer than that are one event.
    new = np.concatenate(([True], np.diff(days) >= min_sep))
    cluster = np.cumsum(new)
    frame = pd.DataFrame({"value": above.to_numpy(), "cluster": cluster}, index=above.index)
    idx = frame.groupby("cluster")["value"].idxmax()
    peaks = above.loc[idx.to_numpy()]
    return pd.DataFrame(
        {
            "peak_date": peaks.index,
            "peak": peaks.to_numpy(dtype=float),
            "water_year": _water_year(pd.DatetimeIndex(peaks.index)),
        }
    )


def pot_counts(
    daily: pd.Series,
    *,
    threshold: Optional[float] = None,
    events_per_year: float = 2.0,
    min_separation_days: int = 7,
    min_days: int = 330,
) -> PotCounts:
    """Count independent peaks over a threshold in each complete water year.

    Exceedance days closer together than ``min_separation_days`` form one event, whose
    peak is the largest daily value in it. Without a ``threshold``, it is the highest
    threshold at which the record averages ``events_per_year`` events per complete year
    (two is common in trend studies of rain-driven floods; Mallakpour and Villarini
    2015). A snowmelt river may reach two only near its median flow, since its freshet is
    one event; that is logged as a warning (see ``quantile``), and one event a year or an
    explicit threshold suits it better.

    Parameters
    ----------
    daily : pandas.Series
        Daily values indexed by date; missing values are dropped.
    threshold : float, optional
    events_per_year : float, default 2.0
        Used only when ``threshold`` is None.
    min_separation_days : int, default 7
    min_days : int, default 330
        Days a water year needs to count as complete.

    Returns
    -------
    PotCounts

    Raises
    ------
    ValueError
        Not a date index, duplicate dates, fewer than five complete water years, or a
        bad ``events_per_year`` / ``min_separation_days``.
    """
    if not isinstance(daily.index, pd.DatetimeIndex):
        raise ValueError("daily must be indexed by date")
    if min_separation_days < 1:
        raise ValueError(f"min_separation_days must be >= 1, got {min_separation_days}")
    s = daily.dropna().sort_index()
    s.index = pd.DatetimeIndex(s.index).tz_localize(None).normalize()
    if s.index.duplicated().any():
        raise ValueError("daily has more than one value for some dates")
    wy = pd.Series(_water_year(pd.DatetimeIndex(s.index)), index=s.index)
    years_, n_days = np.unique(wy.to_numpy(), return_counts=True)
    complete = [int(y) for y, n in zip(years_, n_days) if n >= min_days]
    if len(complete) < 5:
        raise ValueError(f"need at least 5 complete water years, got {len(complete)}")
    s = s[wy.isin(complete).to_numpy()]
    if threshold is None:
        if not events_per_year > 0:
            raise ValueError(f"events_per_year must be > 0, got {events_per_year}")
        target = events_per_year * len(complete)
        # The count is not monotone in the threshold (a snowmelt freshet is one event
        # over any moderate threshold), so scan down from the top and keep the highest
        # threshold that reaches the target.
        candidates = np.unique(np.quantile(s.to_numpy(), np.linspace(0.5, 0.9995, 300)))[::-1]
        threshold = float(candidates[-1])
        miss: Optional[float] = None
        for cand in candidates:
            if len(_pot_events(s, float(cand), min_separation_days)) >= target:
                threshold = float(cand)
                break
            miss = float(cand)
        if miss is not None and threshold < miss:
            # Refine between the last miss and the first hit.
            lo, hi = threshold, miss
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if len(_pot_events(s, mid, min_separation_days)) >= target:
                    lo = mid
                else:
                    hi = mid
            threshold = lo
    quantile = float(np.mean(s.to_numpy() <= threshold))
    if quantile < POT_LOW_QUANTILE:
        logger.warning(
            "POT threshold %.4g is exceeded on %.0f%% of days: the record does not average "
            "%g independent events a year over a high threshold (a snowmelt regime?); "
            "consider events_per_year=1 or an explicit threshold",
            threshold,
            100 * (1 - quantile),
            events_per_year,
        )
    events = _pot_events(s, float(threshold), min_separation_days)
    counts = (
        events.groupby("water_year").size().reindex(complete, fill_value=0).astype(int)
        if len(events)
        else pd.Series(0, index=complete, dtype=int)
    )
    counts.index.name = "water_year"
    return PotCounts(float(threshold), counts, events, int(min_separation_days), quantile)


def pot_count_trend(counts: pd.Series, *, alpha: float = 0.05) -> Dict[str, Any]:
    """Trend in annual event counts: Poisson regression on year, and Mann-Kendall.

    Fits ``log E[N_t] = a + b (t - mean t)`` by maximum likelihood and tests ``b = 0``
    with the likelihood-ratio statistic. When the counts are overdispersed (Pearson
    dispersion above 1), the statistic is divided by the dispersion (quasi-Poisson), so
    clustering of wet years does not pass for a trend.

    Parameters
    ----------
    counts : pandas.Series
        Events per year, indexed by year (e.g. :class:`PotCounts`'s ``counts``).
    alpha : float, default 0.05

    Returns
    -------
    dict
        ``n_years``, ``total_events``, ``mean_per_year``, ``slope_log_per_year`` (``b``),
        ``change_per_decade_pct`` (``100 (exp(10 b) - 1)``), ``rate_at_center`` (``exp(a)``)
        and ``center_year`` (mean t), ``dispersion``, ``p_value``,
        ``mann_kendall_p`` and ``trend`` (``"increasing"``, ``"decreasing"`` or
        ``"no trend"`` at ``alpha``, from the Poisson test).

    Raises
    ------
    ValueError
        Fewer than five years, negative or non-integer counts, or no events at all.
    """
    y = np.asarray(counts, dtype=float)
    t = np.asarray(counts.index, dtype=float)
    if y.size < 5:
        raise ValueError(f"need at least 5 years, got {y.size}")
    if np.any(y < 0) or np.any(y != np.round(y)):
        raise ValueError("counts must be non-negative integers")
    if y.sum() == 0:
        raise ValueError("no events in the record")
    center = float(t.mean())
    t = t - center
    x = np.column_stack([np.ones_like(t), t])
    beta = np.array([math.log(y.mean()), 0.0])
    for _ in range(100):
        mu = np.exp(x @ beta)
        step = np.linalg.solve(x.T @ (mu[:, None] * x), x.T @ (y - mu))
        beta = beta + step
        if np.max(np.abs(step)) < 1e-12:
            break
    mu = np.exp(x @ beta)
    pos = y > 0
    ll1 = float(np.sum(y[pos] * np.log(mu[pos])) - np.sum(mu))
    ll0 = float(np.sum(y[pos]) * math.log(y.mean()) - y.sum())
    deviance = max(0.0, 2.0 * (ll1 - ll0))
    dispersion = float(np.sum((y - mu) ** 2 / mu) / (y.size - 2))
    stat = deviance / max(dispersion, 1.0)
    p = float(stats.chi2.sf(stat, 1))
    b = float(beta[1])
    trend = "no trend" if p >= alpha else ("increasing" if b > 0 else "decreasing")
    return {
        "n_years": int(y.size),
        "total_events": int(y.sum()),
        "mean_per_year": float(y.mean()),
        "slope_log_per_year": b,
        "change_per_decade_pct": 100.0 * (math.exp(10.0 * b) - 1.0),
        "rate_at_center": float(math.exp(beta[0])),
        "center_year": center,
        "dispersion": dispersion,
        "p_value": p,
        "mann_kendall_p": mann_kendall(y).p_value,
        "trend": trend,
    }


# ------------------------------------------------------------- report --


@dataclass
class StationarityReport:
    """:func:`stationarity_report`'s output.

    Attributes
    ----------
    screen : dict
        :func:`stationarity_screen`'s result.
    pot : dict or None
        :func:`pot_count_trend`'s result, when counts were given.
    markdown : str
        A one-page summary: the record, a table of the tests, the flags, the candidate
        causes of a change point and the event-count trend.
    figure : matplotlib.figure.Figure
        The series with Sen's trend, the change point and the segment medians; Pettitt's
        statistic against its critical value; and the event counts when given.
    """

    screen: Dict[str, Any]
    pot: Optional[Dict[str, Any]]
    markdown: str
    figure: Any

    def save(self, path: Any) -> None:
        """Write the figure; the format follows the extension (``.png``, ``.pdf``, ...)."""
        self.figure.savefig(path, dpi=150, bbox_inches="tight")


def _yes(p: float, alpha: float) -> str:
    return "significant" if p < alpha else "not significant"


def stationarity_report(
    values: pd.Series,
    *,
    log: bool = True,
    alpha: float = 0.05,
    peak_codes: Optional[Sequence[Any]] = None,
    events: Iterable[Tuple[str, int]] = (),
    pot: Optional[PotCounts] = None,
    title: Optional[str] = None,
    seed: Optional[int] = 0,
) -> StationarityReport:
    """A one-page stationarity summary of an annual series (roadmap section 6.1).

    Runs :func:`stationarity_screen` (and :func:`pot_count_trend` on ``pot``) and lays
    the result out as a figure and a Markdown page, for attaching to a Bulletin 17C
    analysis. Screening only: nothing here changes an estimate.

    Parameters
    ----------
    values : pandas.Series
        Annual values indexed by water year (e.g. annual peaks).
    log, alpha, peak_codes, events, seed
        As :func:`stationarity_screen` takes them.
    pot : PotCounts, optional
        Event counts from :func:`pot_counts` on the gage's daily record.
    title : str, optional
        Heading, e.g. the gage number and name.

    Returns
    -------
    StationarityReport

    Raises
    ------
    ValueError
        As :func:`stationarity_screen` raises.
    """
    from matplotlib.figure import Figure
    from matplotlib.ticker import MaxNLocator

    events = list(events)
    screen = stationarity_screen(
        values, log=log, alpha=alpha, peak_codes=peak_codes, events=events, seed=seed
    )
    pot_result = pot_count_trend(pot.counts, alpha=alpha) if pot is not None else None
    s = values.dropna().sort_index()
    years = np.asarray(s.index, dtype=int)
    x = np.asarray(s, dtype=float)
    xt = np.log10(x) if log else x
    slope, intercept = sens_slope(xt, years)
    pt = screen["pettitt"]

    n_panels = 3 if pot is not None else 2
    fig = Figure(figsize=(8.5, 3.0 * n_panels))
    axes = fig.subplots(n_panels, 1)
    ax = axes[0]
    ax.plot(years, x, "o-", ms=3, lw=0.8, color="0.25", label="annual value")
    line = intercept + slope * years
    ax.plot(years, 10**line if log else line, color="C0", label="Sen's slope")
    if pt["year"] is not None:
        before, after = x[years < pt["year"]], x[years >= pt["year"]]
        if before.size and after.size:
            ax.hlines(np.median(before), years[0], pt["year"], colors="C3", linestyles=":")
            ax.hlines(np.median(after), pt["year"], years[-1], colors="C3", linestyles=":")
        ax.axvline(
            pt["year"],
            color="C3",
            ls="--",
            lw=1,
            label=f"Pettitt change {pt['year']} (p = {pt['p_value']:.3f})",
        )
    for label, year in events:
        if years[0] <= year <= years[-1]:
            ax.axvline(year, color="C2", lw=0.8, alpha=0.6)
    if log:
        ax.set_yscale("log")
    ax.set_ylabel("annual value")
    ax.set_title(title or "Stationarity screen")
    ax.legend(fontsize=7, loc="best")

    ax = axes[1]
    u = np.abs(_pettitt_u(xt))
    n = xt.size
    k_crit = math.sqrt(-math.log(alpha / 2.0) * (n**3 + n**2) / 6.0)
    ax.plot(years[1:], u, color="0.25")
    ax.axhline(k_crit, color="C3", ls="--", lw=1, label=f"critical K at {alpha:g}")
    ax.set_ylabel("|Pettitt U_t|")
    ax.legend(fontsize=7, loc="best")

    if pot is not None and pot_result is not None:
        ax = axes[2]
        cy = np.asarray(pot.counts.index, dtype=int)
        ax.bar(cy, pot.counts.to_numpy(), color="0.6", width=0.8)
        rate = pot_result["rate_at_center"] * np.exp(
            pot_result["slope_log_per_year"] * (cy - pot_result["center_year"])
        )
        ax.plot(cy, rate, color="C0", label="Poisson trend")
        ax.set_ylabel(f"events > {pot.threshold:,.0f}")
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(fontsize=7, loc="best")
    axes[-1].set_xlabel("water year")
    fig.tight_layout()

    mk = screen["mann_kendall_hamed_rao"]
    lomb = screen["lombard"]
    slope_key = "sens_slope_log10_per_year" if log else "sens_slope_per_year"
    rows = [
        ("Mann-Kendall (Hamed-Rao)", f"tau = {mk['tau']:.3f}", mk["p_value"]),
        (
            "Mann-Kendall (original)",
            f"tau = {screen['mann_kendall']['tau']:.3f}",
            screen["mann_kendall"]["p_value"],
        ),
        ("Spearman", f"rho = {screen['spearman']['rho']:.3f}", screen["spearman"]["p_value"]),
        ("Pettitt", f"{pt['year']}, shift {pt['shift']}", pt["p_value"]),
    ]
    for key, label in (
        ("wilcoxon_abrupt", "Lombard Wilcoxon, abrupt"),
        ("wilcoxon_smooth", "Lombard Wilcoxon, smooth"),
        ("mood_abrupt", "Lombard Mood, abrupt"),
        ("mood_smooth", "Lombard Mood, smooth"),
    ):
        r = lomb[key]
        where = (
            f"{r['start_year']}"
            if r["start_year"] == r["end_year"]
            else f"{r['start_year']}-{r['end_year']}"
        )
        rows.append((label, where, r["p_value"]))
    lines = [
        f"# {title or 'Stationarity screen'}",
        "",
        f"Water years {screen['years'][0]}-{screen['years'][1]} ({screen['n']} values"
        f"{', log10' if log else ''}). Sen's slope "
        f"{screen[slope_key]:+.4f} per year"
        + (f" ({100 * (10 ** (10 * screen[slope_key]) - 1):+.1f}% per decade)." if log else "."),
        "",
        "| Test | Result | p | At " + f"{alpha:g}" + " |",
        "|---|---|---|---|",
        *(f"| {name} | {res} | {p:.3f} | {_yes(p, alpha)} |" for name, res, p in rows),
        "",
        "## Flags",
        "",
        *([f"- {f}" for f in screen["flags"]] or ["- none at this level"]),
    ]
    att = screen.get("attribution")
    if att is not None:
        lines += ["", "## Change-point attribution", "", f"- {att['message']}"]
    if pot is not None and pot_result is not None:
        lines += [
            "",
            "## Peaks over a threshold",
            "",
            f"- threshold {pot.threshold:,.0f} (exceeded on {100 * (1 - pot.quantile):.1f}% "
            f"of days), events at least {pot.min_separation_days} days "
            f"apart: {pot_result['total_events']} in {pot_result['n_years']} years "
            f"({pot_result['mean_per_year']:.2f} per year)"
            + (
                "; a low threshold for floods (snowmelt regime?)"
                if pot.quantile < POT_LOW_QUANTILE
                else ""
            ),
            f"- Poisson trend {pot_result['change_per_decade_pct']:+.1f}% per decade, "
            f"p = {pot_result['p_value']:.3f} (dispersion {pot_result['dispersion']:.2f}); "
            f"Mann-Kendall p = {pot_result['mann_kendall_p']:.3f}",
        ]
    lines += [
        "",
        "Screening only (Bulletin 17C asks the analyst to evaluate stationarity; it does not "
        "prescribe a test). A flag is a reason to look at the record, not a correction.",
        "",
    ]
    return StationarityReport(screen, pot_result, "\n".join(lines), fig)
