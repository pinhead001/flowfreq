"""
Duration-quantile transfer variance: flow-duration statistics and gage pairs.

Methow item 7, part 3. Extends :mod:`flowfreq.monthly_transfer` from monthly means to
flow-duration quantiles. A duration quantile is a lowercase q: q10, q50 and q90 are the
daily flows equaled or exceeded on 10, 50 and 90 percent of days.

For each complete water year (October to September, every day present), a gage's daily
flows give that year's q_p, with the Weibull plotting position
(``numpy.quantile(..., 1 - p/100, method="weibull")``). Over a water-year window the
statistic is the **median** of those yearly values, the median annual flow-duration curve
(Vogel and Fennessey, 1994), or their **mean**. Treating each year's q_p as one observation
keeps years independent, as the monthly statistics do; the period-of-record curve, whose
daily values are strongly autocorrelated, has no such simple sampling variance.

The log10 sampling variance and the concurrent-year covariance between gages are the
monthly module's formulas (see :mod:`flowfreq.monthly_transfer`), applied to the yearly
q_p series. Pairs are transferred by the runoff-volume ratio with exponent 1, and the GLS
blend and leave-one-out are the monthly ones, keyed by exceedance percentage.

What this module cannot do
--------------------------
- **One window per call**, and **complete water years only**: a year with a missing day
  is dropped.
- **Zero low flows are excluded, not modelled.** A median statistic needs every yearly
  value positive (it works in log space), so an intermittent gage's q90 is reported in
  :attr:`DurationStatistics.excluded`.
- **The donor correlation and the runoff depth are the caller's**, as in the monthly
  module.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from flowfreq.monthly_transfer import (
    _STATISTICS,
    _check_daily,
    _loo_core,
    _pairs_core,
    _slice,
    _tabulate,
    _transfers_core,
)
from flowfreq.transfer_variance import BlendResult, DissimilarityModel, combine

logger = logging.getLogger(__name__)

_KEY = "exceedance"
_YEARLY_COLUMNS = ["site_no", "water_year", _KEY, "flow_cfs"]


@dataclass(frozen=True)
class DurationStatistics:
    """Flow-duration statistics of a set of gages over one water-year window.

    Attributes
    ----------
    statistic : str
        ``"median"`` or ``"mean"`` of the yearly duration quantiles.
    start_water_year, end_water_year : int
        The window, inclusive.
    exceedances : tuple of int
        The exceedance percentages, e.g. ``(10, 50, 90)``.
    yearly : pd.DataFrame
        Columns ``site_no``, ``water_year``, ``exceedance``, ``flow_cfs``: each complete
        water year's q_p.
    table : pd.DataFrame
        Columns ``site_no``, ``exceedance``, ``value_cfs``, ``log_q``, ``n_years``,
        ``variance``; a row for every site and exceedance. Excluded rows are NaN.
    covariance : dict of int to pd.DataFrame
        Per exceedance, the sampling covariance of ``log_q`` by ``site_no``.
    excluded : pd.DataFrame
        Columns ``site_no``, ``exceedance``, ``reason``.
    """

    statistic: str
    start_water_year: int
    end_water_year: int
    exceedances: Tuple[int, ...]
    yearly: pd.DataFrame
    table: pd.DataFrame
    covariance: Dict[int, pd.DataFrame]
    excluded: pd.DataFrame


def _yearly_quantiles(
    site: str, daily: pd.Series, start: int, end: int, exceedances: Sequence[int]
) -> List[Tuple[str, int, int, float]]:
    """``(site, water_year, p, q_p)`` for each complete water year in the window."""
    values = daily.dropna()
    dates = pd.DatetimeIndex(values.index).normalize()
    water_year = np.asarray(dates.year + (dates.month >= 10), dtype=int)
    flows = np.asarray(values, dtype=float)
    rows: List[Tuple[str, int, int, float]] = []
    for wy in range(start, end + 1):
        sel = water_year == wy
        x = flows[sel]
        days = (pd.Timestamp(f"{wy}-09-30") - pd.Timestamp(f"{wy - 1}-10-01")).days + 1
        if len(dates[sel].unique()) != days or len(x) != days:
            continue
        for p in exceedances:
            rows.append(
                (site, wy, int(p), float(np.quantile(x, 1.0 - p / 100.0, method="weibull")))
            )
    return rows


def duration_statistics(
    daily: Mapping[str, pd.Series],
    *,
    start_water_year: int,
    end_water_year: int,
    exceedances: Sequence[int] = (10, 50, 90),
    statistic: str = "median",
    min_years: int = 10,
) -> DurationStatistics:
    """Median (or mean) annual duration quantiles of each gage, with sampling variances.

    Parameters
    ----------
    daily : mapping of str to pd.Series
        Daily mean flow (cfs) per site, indexed by a ``DatetimeIndex``.
    start_water_year, end_water_year : int
        Inclusive water-year window.
    exceedances : sequence of int, default (10, 50, 90)
        Exceedance percentages, each strictly between 0 and 100.
    statistic : {"median", "mean"}, default "median"
        Statistic of the yearly duration quantiles.
    min_years : int, default 10
        Fewest complete water years a gage needs. Must be at least 2.

    Returns
    -------
    DurationStatistics

    Raises
    ------
    ValueError
        An exceedance outside (0, 100) or repeated, an unknown ``statistic``, a reversed
        window, or ``min_years`` below 2.
    TypeError
        A series is not indexed by a ``DatetimeIndex``.
    """
    keys = tuple(int(p) for p in exceedances)
    if not keys or any(not 0 < p < 100 for p in keys) or len(set(keys)) != len(keys):
        raise ValueError(f"exceedances must be distinct integers in (0, 100), got {exceedances}")
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
    rows: List[Tuple[str, int, int, float]] = []
    for key, series in daily.items():
        rows.extend(_yearly_quantiles(str(key), series, start_water_year, end_water_year, keys))
    yearly = pd.DataFrame(rows, columns=_YEARLY_COLUMNS)
    groups: Dict[Tuple[str, int], pd.Series] = {}
    for _, sub in yearly.groupby(["site_no", _KEY]):
        site, p = str(sub["site_no"].iloc[0]), int(sub[_KEY].iloc[0])
        groups[(site, p)] = pd.Series(
            sub["flow_cfs"].to_numpy(dtype=float), index=sub["water_year"].to_numpy(dtype=int)
        )
    table, covariance, excluded = _tabulate(groups, sites, keys, statistic, min_years, _KEY)
    return DurationStatistics(
        statistic=statistic,
        start_water_year=start_water_year,
        end_water_year=end_water_year,
        exceedances=keys,
        yearly=yearly,
        table=table,
        covariance=covariance,
        excluded=excluded,
    )


def _key(stats: DurationStatistics, exceedance: int) -> int:
    p = int(exceedance)
    if p not in stats.exceedances:
        raise KeyError(f"exceedance {exceedance} is not in the statistics {stats.exceedances}")
    return p


def build_duration_pairs(
    stats: DurationStatistics,
    exceedance: int,
    gages: pd.DataFrame,
    connections: Sequence[Tuple[str, str]],
    *,
    attributes: Sequence[str] = (),
) -> pd.DataFrame:
    """Transfer error of every connected gage pair for one duration quantile.

    As :func:`flowfreq.monthly_transfer.build_monthly_pairs`, keyed by ``exceedance``
    instead of ``month``; the result feeds
    :func:`flowfreq.transfer_variance.fit_dissimilarity`.

    Raises
    ------
    KeyError
        ``exceedance`` is not in ``stats``, or a connection names an unknown site.
    ValueError
        ``gages`` has no ``runoff_depth`` column.
    """
    p = _key(stats, exceedance)
    log_q, cov = _slice(stats, _KEY, p)
    return _pairs_core(log_q, cov, stats.statistic, _KEY, p, gages, connections, attributes)


def duration_transfers(
    stats: DurationStatistics,
    exceedance: int,
    gages: pd.DataFrame,
    target: str,
    donors: Sequence[str],
    *,
    attributes: Sequence[str] = (),
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """One target's transfers and its donors' covariance for one duration quantile.

    As :func:`flowfreq.monthly_transfer.monthly_transfers`, keyed by ``exceedance``.
    """
    p = _key(stats, exceedance)
    log_q, cov = _slice(stats, _KEY, p)
    return _transfers_core(log_q, cov, gages, target, donors, attributes, f"q{p}")


def blend_duration(
    model: DissimilarityModel,
    stats: DurationStatistics,
    exceedance: int,
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
    """GLS blend of one target's duration-quantile transfers, and optionally a regression.

    As :func:`flowfreq.monthly_transfer.blend_month`, keyed by ``exceedance``.
    """
    transfers, dcov = duration_transfers(
        stats, exceedance, gages, target, donors, attributes=attributes
    )
    return combine(
        model,
        transfers,
        regression,
        donor_cov=dcov,
        donor_correlation=donor_correlation,
        regression_correlation=regression_correlation,
        level=level,
    )


def duration_leave_one_out(
    stats: DurationStatistics,
    exceedance: int,
    gages: pd.DataFrame,
    connections: Sequence[Tuple[str, str]],
    terms: Sequence[str] = ("log_ar",),
    *,
    attributes: Sequence[str] = (),
    donor_correlation: float = 0.0,
    level: float = 0.90,
    powers: Optional[Mapping[str, float]] = None,
) -> pd.DataFrame:
    """Score the duration-quantile blend by holding out each connected gage in turn.

    As :func:`flowfreq.monthly_transfer.monthly_leave_one_out`, keyed by ``exceedance``:
    one row per scored gage with ``site_no``, ``exceedance``, ``statistic``,
    ``n_donors``, ``observed``, ``blend``, ``blend_variance``, ``blend_error``,
    ``blend_covered``, ``blend_covered_naive`` and ``note``.
    """
    p = _key(stats, exceedance)
    log_q, cov = _slice(stats, _KEY, p)
    return _loo_core(
        log_q,
        cov,
        stats.statistic,
        _KEY,
        p,
        gages,
        connections,
        terms,
        attributes,
        donor_correlation,
        level,
        powers,
    )
