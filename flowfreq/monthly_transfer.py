"""
Monthly transfer-variance inputs: per-month statistics and gage pairs.

Extends :mod:`flowfreq.transfer_variance` from flood quantiles to monthly flows. For
each calendar month, a gage's yearly monthly means over a water-year window (October to
September, labelled by the September year) give one statistic, the **mean** or the
**median** of those yearly means. Each comes with its log10 sampling variance and, from
concurrent water years, the sampling covariance between gages. Pairs are transferred by
the runoff-volume ratio with exponent 1, so
:func:`flowfreq.transfer_variance.fit_dissimilarity` fits ``s2_d`` per month unchanged.

Formulas
--------
With ``x`` the ``n`` yearly means of one gage and month:

- ``"mean"``: ``value = mean(x)``, ``log_q = log10(value)`` and, by the delta method,
  ``variance = var(x, ddof=1) / (n value^2 ln(10)^2)``.
- ``"median"``: ``value = median(x)``, ``log_q = median(log10 x)`` and
  ``variance = (pi/2) var(log10 x, ddof=1) / n``.

For two gages with concurrent years ``C`` (``nc = |C|``) and their own ``na``, ``nb``:

- ``"mean"``: ``nc cov(x_a[C], x_b[C]) / (na nb mean(x_a) mean(x_b) ln(10)^2)``.
- ``"median"``: ``nc sd(log10 x_a[C]) sd(log10 x_b[C]) asin(rho) / (na nb)``, with
  ``rho`` the correlation of the logs over ``C`` clipped to [-1, 1]. These are the
  bivariate-normal asymptotics ``n Var = (pi/2) s^2`` and ``n Cov = s_a s_b asin(rho)``.
- ``nc < 2`` gives 0.

A pair's transfer is ``log_q[donor] + log10(At/Ag) + log10(Rt/Rg)``, where ``R`` is the
runoff depth, and its known variance is ``V_d + V_r - 2 C_dr``.

What this module cannot do
--------------------------
- **One window per call.** Each call covers one water-year window.
- **The donor correlation is the caller's.** :func:`monthly_leave_one_out` and
  :func:`blend_month` take it as a fixed value (default 0); the flood estimator
  (:func:`flowfreq.transfer_variance.estimate_donor_correlation`) is not carried over.
- **The runoff depth is the caller's.** flowfreq does not estimate it.
- **Zero-flow months are excluded, not modelled.** They are reported in
  :attr:`MonthlyStatistics.excluded`.
"""

from __future__ import annotations

import calendar
import logging
import math
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import norm

from flowfreq.transfer_variance import BlendResult, DissimilarityModel, combine, fit_dissimilarity

logger = logging.getLogger(__name__)

_LN10 = math.log(10.0)
_STATISTICS = ("mean", "median")
_YEARLY_COLUMNS = ["site_no", "water_year", "month", "flow_cfs"]


@dataclass(frozen=True)
class MonthlyStatistics:
    """Monthly statistics of a set of gages over one water-year window.

    Attributes
    ----------
    statistic : str
        ``"mean"`` or ``"median"``.
    start_water_year, end_water_year : int
        The window, inclusive.
    yearly : pd.DataFrame
        Columns ``site_no``, ``water_year``, ``month``, ``flow_cfs``: the mean daily
        flow of each complete gage-month.
    table : pd.DataFrame
        Columns ``site_no``, ``month``, ``value_cfs``, ``log_q``, ``n_years``,
        ``variance``; a row for every site and month 1-12. Excluded rows are NaN.
    covariance : dict of int to pd.DataFrame
        Per month, the sampling covariance of ``log_q``, indexed and columned by
        ``site_no`` in input order. The diagonal is ``variance``.
    excluded : pd.DataFrame
        Columns ``site_no``, ``month``, ``reason``.
    """

    statistic: str
    start_water_year: int
    end_water_year: int
    yearly: pd.DataFrame
    table: pd.DataFrame
    covariance: Dict[int, pd.DataFrame]
    excluded: pd.DataFrame


def _check_daily(site: str, series: pd.Series) -> None:
    """Refuse a series with more than one value on a date (sub-daily or duplicated)."""
    dates = pd.DatetimeIndex(series.dropna().index).normalize()
    if dates.has_duplicates:
        raise ValueError(
            f"daily[{site!r}] must have one value per day; it repeats "
            f"{dates[dates.duplicated()][0].date()} (sub-daily or duplicated data)"
        )


def _yearly_means(
    site: str, daily: pd.Series, start: int, end: int
) -> List[Tuple[str, int, int, float]]:
    """Complete-month means of one site: (site, water_year, month, mean)."""
    values = daily.dropna()
    index = pd.DatetimeIndex(values.index)
    water_year = np.asarray(index.year + (index.month >= 10), dtype=int)
    months = np.asarray(index.month, dtype=int)
    keep = (water_year >= start) & (water_year <= end)
    frame = pd.DataFrame(
        {
            "wy": water_year[keep],
            "month": months[keep],
            "v": np.asarray(values, dtype=float)[keep],
        }
    )
    rows: List[Tuple[str, int, int, float]] = []
    if frame.empty:
        return rows
    agg = frame.groupby(["wy", "month"])["v"].agg(["mean", "count"]).reset_index()
    for wy, month, mean, count in zip(
        agg["wy"].tolist(), agg["month"].tolist(), agg["mean"].tolist(), agg["count"].tolist()
    ):
        year = int(wy) - 1 if int(month) >= 10 else int(wy)
        if int(count) == calendar.monthrange(year, int(month))[1]:
            rows.append((site, int(wy), int(month), float(mean)))
    return rows


def _own_statistic(x: np.ndarray, statistic: str) -> Tuple[float, float, float]:
    """``(value, log_q, variance)`` of the yearly means ``x`` (all positive)."""
    n = len(x)
    if statistic == "mean":
        value = float(np.mean(x))
        variance = float(np.var(x, ddof=1)) / (n * value**2 * _LN10**2)
        return value, math.log10(value), variance
    lx = np.log10(x)
    variance = (math.pi / 2.0) * float(np.var(lx, ddof=1)) / n
    return float(np.median(x)), float(np.median(lx)), variance


def _cross_covariance(xa: pd.Series, xb: pd.Series, statistic: str) -> float:
    """Sampling covariance of two sites' statistics from their concurrent years."""
    common = xa.index.intersection(xb.index)
    nc = len(common)
    if nc < 2:
        return 0.0
    a = np.asarray(xa.loc[common], dtype=float)
    b = np.asarray(xb.loc[common], dtype=float)
    na, nb = len(xa), len(xb)
    if statistic == "mean":
        cov = float(np.cov(a, b, ddof=1)[0, 1])
        return nc * cov / (na * nb * float(xa.mean()) * float(xb.mean()) * _LN10**2)
    la, lb = np.log10(a), np.log10(b)
    sa, sb = float(np.std(la, ddof=1)), float(np.std(lb, ddof=1))
    if sa == 0.0 or sb == 0.0:
        return 0.0
    rho = float(np.clip(np.corrcoef(la, lb)[0, 1], -1.0, 1.0))
    return nc * sa * sb * math.asin(rho) / (na * nb)


def monthly_statistics(
    daily: Mapping[str, pd.Series],
    *,
    start_water_year: int,
    end_water_year: int,
    statistic: str = "mean",
    min_years: int = 10,
) -> MonthlyStatistics:
    """Monthly mean or median flow of each gage, with sampling variances.

    Parameters
    ----------
    daily : mapping of str to pd.Series
        Daily mean flow (cfs) per site, indexed by a ``DatetimeIndex``.
    start_water_year, end_water_year : int
        Inclusive water-year window (October-September, labelled by the September year).
    statistic : {"mean", "median"}, default "mean"
        Statistic of the yearly monthly means.
    min_years : int, default 10
        Fewest complete years a gage-month needs. Must be at least 2.

    Returns
    -------
    MonthlyStatistics
        The yearly means, the per-site-month table, per-month covariance matrices
        and the exclusions.

    Raises
    ------
    ValueError
        ``statistic`` is not ``"mean"`` or ``"median"``, the window is reversed, or
        ``min_years`` is below 2.
    TypeError
        A series is not indexed by a ``DatetimeIndex``.
    """
    if statistic not in _STATISTICS:
        raise ValueError(f"statistic must be one of {_STATISTICS}, got {statistic!r}")
    if start_water_year > end_water_year:
        raise ValueError(
            f"start water year {start_water_year} is after end water year {end_water_year}"
        )
    if min_years < 2:
        raise ValueError(f"min_years must be >= 2, got {min_years}")
    for site, series in daily.items():
        if not isinstance(series.index, pd.DatetimeIndex):
            raise TypeError(f"daily[{site!r}] must be indexed by a DatetimeIndex")
        _check_daily(str(site), series)

    sites = [str(s) for s in daily]
    yearly_rows: List[Tuple[str, int, int, float]] = []
    for key, series in daily.items():
        yearly_rows.extend(_yearly_means(str(key), series, start_water_year, end_water_year))
    yearly = pd.DataFrame(yearly_rows, columns=_YEARLY_COLUMNS)

    # Yearly means per (site, month), indexed by water year.
    groups: Dict[Tuple[str, int], pd.Series] = {}
    for gkey, sub in yearly.groupby(["site_no", "month"]):
        gsite = str(sub["site_no"].iloc[0])
        gmonth = int(sub["month"].iloc[0])
        groups[(gsite, gmonth)] = pd.Series(
            sub["flow_cfs"].to_numpy(dtype=float), index=sub["water_year"].to_numpy(dtype=int)
        )

    table, covariance, excluded = _tabulate(
        groups, sites, range(1, 13), statistic, min_years, "month", yearly_name="yearly mean"
    )

    return MonthlyStatistics(
        statistic=statistic,
        start_water_year=start_water_year,
        end_water_year=end_water_year,
        yearly=yearly,
        table=table,
        covariance=covariance,
        excluded=excluded,
    )


def _cell(frame: pd.DataFrame, row: str, col: str) -> float:
    """Return one numeric cell of a covariance frame as a float."""
    return float(np.asarray(frame.loc[row, col], dtype=float))


def build_monthly_pairs(
    stats: MonthlyStatistics,
    month: int,
    gages: pd.DataFrame,
    connections: Sequence[Tuple[str, str]],
    *,
    attributes: Sequence[str] = (),
) -> pd.DataFrame:
    """Transfer error of every connected gage pair for one month.

    Parameters
    ----------
    stats : MonthlyStatistics
        Output of :func:`monthly_statistics`.
    month : int
        Calendar month, 1-12.
    gages : pd.DataFrame
        Columns ``site_no``, ``drainage_area_sqmi``, ``runoff_depth`` (any consistent
        units; 1.0 everywhere gives the plain area ratio) and every ``attributes`` column.
    connections : sequence of (str, str)
        Flow-connected pairs as ``(donor, receiver)``.
    attributes : sequence of str
        Gage columns carried as ``d_<attribute>`` (receiver minus donor).

    Returns
    -------
    pd.DataFrame
        Columns ``donor``, ``receiver``, ``month``, ``statistic``, ``area_ratio``,
        ``log_ar``, ``log_rr``, ``exponent``, ``transfer``, ``observed``, ``error``,
        ``known_variance`` and the ``d_`` columns. Connections with an excluded
        donor or receiver are dropped (logged at INFO). The result feeds
        :func:`flowfreq.transfer_variance.fit_dissimilarity`.

    Raises
    ------
    ValueError
        ``month`` is outside 1-12 or ``gages`` has no ``runoff_depth`` column.
    KeyError
        A connection names a site that is not in ``gages`` or not in ``stats``.
    """
    if not 1 <= int(month) <= 12:
        raise ValueError(f"month must be in 1-12, got {month}")
    if "runoff_depth" not in gages.columns:
        raise ValueError("gages needs a 'runoff_depth' column")
    log_q, cov = _slice(stats, "month", int(month))
    return _pairs_core(
        log_q,
        cov,
        stats.statistic,
        "month",
        int(month),
        gages,
        connections,
        attributes,
        "monthly statistics",
    )


# ------------------------------------------------------------------- part 2: blend and LOO


def monthly_transfers(
    stats: MonthlyStatistics,
    month: int,
    gages: pd.DataFrame,
    target: str,
    donors: Sequence[str],
    *,
    attributes: Sequence[str] = (),
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """One target's transfers from its donors for one month, as ``combine`` takes them.

    Parameters
    ----------
    stats : MonthlyStatistics
    month : int
        Calendar month, 1-12.
    gages : pd.DataFrame
        As :func:`build_monthly_pairs`.
    target : str
        The receiving site.
    donors : sequence of str
        The donor sites.
    attributes : sequence of str
        Gage columns carried as ``d_<attribute>`` (target minus donor).

    Returns
    -------
    (pd.DataFrame, pd.DataFrame)
        The transfers (``donor``, ``transfer``, ``donor_variance``, ``log_ar``,
        ``log_rr`` and the ``d_`` columns) and the donors' sampling covariance,
        labelled by donor.

    Raises
    ------
    ValueError
        ``month`` is outside 1-12, ``gages`` has no ``runoff_depth`` column, or a donor
        is excluded that month.
    KeyError
        A site is not in ``gages`` or not in ``stats``.
    """
    if not 1 <= int(month) <= 12:
        raise ValueError(f"month must be in 1-12, got {month}")
    if "runoff_depth" not in gages.columns:
        raise ValueError("gages needs a 'runoff_depth' column")
    log_q, cov = _slice(stats, "month", int(month))
    return _transfers_core(
        log_q, cov, gages, target, donors, attributes, f"month {month}", "monthly statistics"
    )


def blend_month(
    model: DissimilarityModel,
    stats: MonthlyStatistics,
    month: int,
    gages: pd.DataFrame,
    target: str,
    donors: Sequence[str],
    *,
    regression: Optional[Tuple[float, float]] = None,
    donor_correlation: float = 0.0,
    regression_correlation: float = 0.0,
    attributes: Sequence[str] = (),
    level: float = 0.90,
) -> BlendResult:
    """GLS blend of one target's monthly transfers, and optionally a regression.

    :func:`monthly_transfers`, then :func:`flowfreq.transfer_variance.combine` with the
    donors' sampling covariance from ``stats`` and no exponent variance (the exponent on
    the runoff-volume ratio is fixed at 1).

    Parameters
    ----------
    model : DissimilarityModel
        Fitted to the month's pairs (:func:`build_monthly_pairs`).
    stats, month, gages, target, donors, attributes
        As :func:`monthly_transfers`.
    regression : (float, float), optional
        A log10 estimate at the target and its prediction variance.
    donor_correlation, regression_correlation, level
        As :func:`flowfreq.transfer_variance.combine`.

    Returns
    -------
    BlendResult
    """
    transfers, dcov = monthly_transfers(stats, month, gages, target, donors, attributes=attributes)
    return combine(
        model,
        transfers,
        regression,
        donor_cov=dcov,
        donor_correlation=donor_correlation,
        regression_correlation=regression_correlation,
        level=level,
    )


def monthly_leave_one_out(
    stats: MonthlyStatistics,
    month: int,
    gages: pd.DataFrame,
    connections: Sequence[Tuple[str, str]],
    terms: Sequence[str] = ("log_ar",),
    *,
    attributes: Sequence[str] = (),
    donor_correlation: float = 0.0,
    level: float = 0.90,
    powers: Optional[Mapping[str, float]] = None,
) -> pd.DataFrame:
    """Score the monthly blend by holding out each connected gage in turn.

    For each gage ``j`` with a connection whose other end is not excluded that month,
    ``s2_d`` is refitted from the month's pairs that do not involve ``j``, its
    neighbours' transfers are blended (:func:`blend_month`), and the result is compared
    with ``j``'s own statistic. Coverage is checked against the variance of
    ``blend - y_j``, ``V_blend + V_j - 2 sum_k w_k C_kj``, since the observed value is
    itself an estimate; ``blend_covered_naive`` uses ``V_blend`` alone. This mirrors
    :func:`flowfreq.transfer_variance.leave_one_out`, without a regression and with a
    fixed donor correlation.

    Parameters
    ----------
    stats, month, gages, connections, attributes
        As :func:`build_monthly_pairs`.
    terms : sequence of str
        Passed to :func:`flowfreq.transfer_variance.fit_dissimilarity`.
    donor_correlation : float, default 0.0
        Correlation of two donors' dissimilarity errors.
    level : float, default 0.90
    powers : mapping of str to float, optional
        Passed to :func:`flowfreq.transfer_variance.fit_dissimilarity`.

    Returns
    -------
    pd.DataFrame
        One row per scored gage: ``site_no``, ``month``, ``statistic``, ``n_donors``,
        ``observed``, ``blend``, ``blend_variance``, ``blend_error``, ``blend_covered``,
        ``blend_covered_naive`` and ``note``. Errors are estimate minus observed, log10.
        A gage whose fold cannot be refitted or blended keeps NaN estimates, and its
        ``note`` says why.
    """
    if not 1 <= int(month) <= 12:
        raise ValueError(f"month must be in 1-12, got {month}")
    log_q, cov = _slice(stats, "month", int(month))
    return _loo_core(
        log_q, cov, stats.statistic, "month", int(month), gages, connections, terms,
        attributes, donor_correlation, level, powers, "monthly statistics",
    )  # fmt: skip


# ------------------------------------------------------------------ key-generic cores
#
# The monthly functions above, and flowfreq.duration_transfer, key their statistics by
# month or by exceedance percentage; these cores take the per-key log10 statistics and
# sampling covariance directly.


def _tabulate(
    groups: Mapping[Tuple[str, int], pd.Series],
    sites: Sequence[str],
    keys: Sequence[int],
    statistic: str,
    min_years: int,
    key_col: str,
    yearly_name: str = "yearly value",
) -> Tuple[pd.DataFrame, Dict[int, pd.DataFrame], pd.DataFrame]:
    """Per site and key: the statistic of its yearly series, with variance and covariance.

    ``groups[(site, key)]`` is the yearly series indexed by water year. Returns the table
    (``site_no``, ``key_col``, ``value_cfs``, ``log_q``, ``n_years``, ``variance``), the
    per-key covariance frames and the exclusions (``site_no``, ``key_col``, ``reason``).
    """
    table_rows: List[Dict[str, object]] = []
    excluded_rows: List[Dict[str, object]] = []
    valid: Dict[Tuple[str, int], pd.Series] = {}
    for site in sites:
        for key in keys:
            x = groups.get((site, key), pd.Series(dtype=float))
            n = len(x)
            value = log_q = variance = float("nan")
            reason = ""
            if n < min_years:
                reason = f"fewer than {min_years} complete years ({n})"
            elif statistic == "mean" and float(x.mean()) <= 0.0:
                reason = "zero period mean flow"
            elif statistic == "median" and bool((x <= 0.0).any()):
                reason = f"zero flow in a {yearly_name} (median works in log space)"
            else:
                value, log_q, variance = _own_statistic(np.asarray(x, dtype=float), statistic)
                valid[(site, key)] = x
            if reason:
                excluded_rows.append({"site_no": site, key_col: key, "reason": reason})
            table_rows.append(
                {
                    "site_no": site,
                    key_col: key,
                    "value_cfs": value,
                    "log_q": log_q,
                    "n_years": n,
                    "variance": variance,
                }
            )
    table = pd.DataFrame(
        table_rows, columns=["site_no", key_col, "value_cfs", "log_q", "n_years", "variance"]
    )
    excluded = pd.DataFrame(excluded_rows, columns=["site_no", key_col, "reason"])
    covariance: Dict[int, pd.DataFrame] = {}
    for key in keys:
        k = len(sites)
        matrix = np.full((k, k), np.nan)
        sel = table[key_col] == key
        var_by_site = dict(
            zip(table.loc[sel, "site_no"].tolist(), table.loc[sel, "variance"].tolist())
        )
        for i, a in enumerate(sites):
            if (a, key) not in valid:
                continue
            matrix[i, i] = float(var_by_site[a])
            for j in range(i + 1, k):
                b = sites[j]
                if (b, key) not in valid:
                    continue
                c = _cross_covariance(valid[(a, key)], valid[(b, key)], statistic)
                matrix[i, j] = matrix[j, i] = c
        covariance[int(key)] = pd.DataFrame(matrix, index=list(sites), columns=list(sites))
    return table, covariance, excluded


def _slice(stats: object, key_col: str, key: int) -> Tuple[Dict[str, float], pd.DataFrame]:
    """``(log_q by site, covariance frame)`` of one key of a statistics object."""
    table = getattr(stats, "table")
    sel = table[table[key_col] == key]
    log_q = dict(zip(sel["site_no"].astype(str).tolist(), sel["log_q"].astype(float).tolist()))
    return log_q, getattr(stats, "covariance")[int(key)]


def _gage_columns(
    gages: pd.DataFrame, attributes: Sequence[str]
) -> Tuple[Dict[str, float], Dict[str, float], Dict[str, Dict[str, float]]]:
    if "runoff_depth" not in gages.columns:
        raise ValueError("gages needs a 'runoff_depth' column")
    site_list = gages["site_no"].astype(str).tolist()
    area = dict(zip(site_list, gages["drainage_area_sqmi"].to_numpy(dtype=float).tolist()))
    depth = dict(zip(site_list, gages["runoff_depth"].to_numpy(dtype=float).tolist()))
    attr = {a: dict(zip(site_list, gages[a].to_numpy(dtype=float).tolist())) for a in attributes}
    return area, depth, attr


def _check_pairs(connections: Sequence[Tuple[str, str]]) -> None:
    """Refuse a gage connected to itself and a pair listed twice (either direction)."""
    seen = set()
    for a_raw, b_raw in connections:
        a, b = str(a_raw), str(b_raw)
        if a == b:
            raise ValueError(f"gage {a} is connected to itself")
        pair = frozenset((a, b))
        if pair in seen:
            raise ValueError(f"the pair {a}-{b} is listed twice")
        seen.add(pair)


def _pairs_core(
    log_q: Mapping[str, float],
    cov: pd.DataFrame,
    statistic: str,
    key_col: str,
    key: int,
    gages: pd.DataFrame,
    connections: Sequence[Tuple[str, str]],
    attributes: Sequence[str],
    source: str = "statistics",
) -> pd.DataFrame:
    """The pair table of :func:`build_monthly_pairs` for one key."""
    _check_pairs(connections)
    area, depth, attr = _gage_columns(gages, attributes)
    rows: List[Dict[str, object]] = []
    for donor_raw, receiver_raw in connections:
        donor, receiver = str(donor_raw), str(receiver_raw)
        for site in (donor, receiver):
            if site not in area:
                raise KeyError(f"connection names {site}, which is not in the gage table")
            if site not in log_q:
                raise KeyError(f"connection names {site}, which is not in the {source}")
        if math.isnan(log_q[donor]) or math.isnan(log_q[receiver]):
            logger.info(
                "%s %d: dropping connection %s -> %s (a gage is excluded)",
                key_col,
                key,
                donor,
                receiver,
            )
            continue
        log_ar = math.log10(area[receiver] / area[donor])
        log_rr = math.log10(depth[receiver] / depth[donor])
        transfer = log_q[donor] + log_ar + log_rr
        observed = log_q[receiver]
        known = (
            _cell(cov, donor, donor)
            + _cell(cov, receiver, receiver)
            - 2.0 * _cell(cov, donor, receiver)
        )
        row: Dict[str, object] = {
            "donor": donor,
            "receiver": receiver,
            key_col: int(key),
            "statistic": statistic,
            "area_ratio": 10.0**log_ar,
            "log_ar": log_ar,
            "log_rr": log_rr,
            "exponent": 1.0,
            "transfer": transfer,
            "observed": observed,
            "error": transfer - observed,
            "known_variance": known,
        }
        for a in attributes:
            row[f"d_{a}"] = attr[a][receiver] - attr[a][donor]
        rows.append(row)
    columns = [
        "donor",
        "receiver",
        key_col,
        "statistic",
        "area_ratio",
        "log_ar",
        "log_rr",
        "exponent",
        "transfer",
        "observed",
        "error",
        "known_variance",
    ] + [f"d_{a}" for a in attributes]
    return pd.DataFrame(rows, columns=columns)


def _transfers_core(
    log_q: Mapping[str, float],
    cov: pd.DataFrame,
    gages: pd.DataFrame,
    target: str,
    donors: Sequence[str],
    attributes: Sequence[str],
    where: str,
    source: str = "statistics",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """The transfers and donor covariance of :func:`monthly_transfers` for one key."""
    area, depth, attr = _gage_columns(gages, attributes)
    target = str(target)
    donor_list = [str(d) for d in donors]
    for site in (target, *donor_list):
        if site not in area:
            raise KeyError(f"{site} is not in the gage table")
        if site not in log_q:
            raise KeyError(f"{site} is not in the {source}")
    rows: List[Dict[str, object]] = []
    for d in donor_list:
        if math.isnan(log_q[d]):
            raise ValueError(f"{where}: donor {d} is excluded (see stats.excluded)")
        log_ar = math.log10(area[target] / area[d])
        log_rr = math.log10(depth[target] / depth[d])
        row: Dict[str, object] = {
            "donor": d,
            "transfer": log_q[d] + log_ar + log_rr,
            "donor_variance": _cell(cov, d, d),
            "log_ar": log_ar,
            "log_rr": log_rr,
        }
        for a in attributes:
            row[f"d_{a}"] = attr[a][target] - attr[a][d]
        rows.append(row)
    dcov = pd.DataFrame(
        [[_cell(cov, a, b) for b in donor_list] for a in donor_list],
        index=donor_list,
        columns=donor_list,
    )
    return pd.DataFrame(rows), dcov


def _loo_core(
    log_q: Mapping[str, float],
    cov: pd.DataFrame,
    statistic: str,
    key_col: str,
    key: int,
    gages: pd.DataFrame,
    connections: Sequence[Tuple[str, str]],
    terms: Sequence[str],
    attributes: Sequence[str],
    donor_correlation: float,
    level: float,
    powers: Optional[Mapping[str, float]],
    source: str = "statistics",
) -> pd.DataFrame:
    """The leave-one-out table of :func:`monthly_leave_one_out` for one key."""
    pairs = _pairs_core(log_q, cov, statistic, key_col, key, gages, connections, attributes, source)
    neighbours: Dict[str, List[str]] = {}
    for a_raw, b_raw in connections:
        a, b = str(a_raw), str(b_raw)
        if math.isnan(log_q.get(a, math.nan)) or math.isnan(log_q.get(b, math.nan)):
            continue
        neighbours.setdefault(a, []).append(b)
        neighbours.setdefault(b, []).append(a)
    zcrit = float(norm.ppf(0.5 + level / 2.0))
    rows: List[Dict[str, object]] = []
    for j in [str(site) for site in gages["site_no"]]:
        nbrs = neighbours.get(j, [])
        if not nbrs:
            continue
        row: Dict[str, object] = {
            "site_no": j,
            key_col: int(key),
            "statistic": statistic,
            "n_donors": len(nbrs),
            "observed": log_q[j],
            "blend": math.nan,
            "blend_variance": math.nan,
            "blend_error": math.nan,
            "blend_covered": math.nan,
            "blend_covered_naive": math.nan,
            "note": "",
        }
        fold = pairs[(pairs["donor"] != j) & (pairs["receiver"] != j)]
        try:
            if fold.empty:
                raise ValueError("no pairs left without this gage")
            model = fit_dissimilarity(fold, terms, powers=powers)
        except ValueError as exc:
            row["note"] = f"s2_d refit failed: {exc}"
            rows.append(row)
            continue
        try:
            transfers, dcov = _transfers_core(
                log_q, cov, gages, j, nbrs, attributes, f"{key_col} {key}", source
            )
            res = combine(
                model,
                transfers,
                None,
                donor_cov=dcov,
                donor_correlation=donor_correlation,
                level=level,
            )
        except ValueError as exc:
            row["note"] = f"blend failed: {exc}"
            rows.append(row)
            continue
        w = res.weights.reindex(nbrs).to_numpy(dtype=float)
        c_kj = np.array([_cell(cov, i, j) for i in nbrs])
        v_err = res.variance + _cell(cov, j, j) - 2.0 * float(w @ c_kj)
        err = res.log_q - log_q[j]
        row.update(
            blend=res.log_q,
            blend_variance=res.variance,
            blend_error=err,
            blend_covered=bool(abs(err) <= zcrit * math.sqrt(max(v_err, 0.0))),
            blend_covered_naive=bool(abs(err) <= zcrit * math.sqrt(res.variance)),
        )
        rows.append(row)
    return pd.DataFrame(rows)
