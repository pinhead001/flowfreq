"""
Naturalization add-back: gaged flow plus the consumptive part of upstream diversions.

A regulated donor carries its regulation into every transfer (:mod:`flowfreq.qppq` does
not naturalize). For a basin whose regulation is irrigation diversion, the natural flow at
the gage is approximately

    Q_natural(t) = Q_gaged(t) + sum_d  rate_d(t) * (1 - return_fraction_d) * in_effect_d(t)

where ``rate_d`` is diversion ``d``'s withdrawal, ``return_fraction_d`` the share that
comes back above the gage (canal seepage, tailwater), and ``in_effect_d`` its season and
period. This is the methow-flow low-flow template (its 2002 canal inventory).

What this module does not do
----------------------------
- **It does not estimate the inputs.** Diversion rates, return fractions and periods are
  the caller's, from water-right records, canal inventories or metering.
- **No routing or lag.** A diversion is added back on the day it is withdrawn, which suits
  daily or monthly series for diversions close upstream of the gage.
- **No storage.** Reservoir regulation needs a mass balance with storage change, not an
  add-back.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

Rate = Union[float, pd.Series]


@dataclass(frozen=True)
class Diversion:
    """One upstream diversion.

    Parameters
    ----------
    name : str
        A label, unique within one :func:`naturalize` call.
    rate_cfs : float or pd.Series
        Withdrawal in cfs: a constant, or a series indexed by date (missing dates add
        nothing and are counted in the summary). Must be >= 0.
    return_fraction : float, default 0.0
        Share of the withdrawal that returns to the stream above the gage, in [0, 1].
    months : tuple of int, default every month
        Months in which the diversion operates (its season).
    start, end : date-like, optional
        First and last dates the diversion is in effect, inclusive.
    """

    name: str
    rate_cfs: Rate
    return_fraction: float = 0.0
    months: Tuple[int, ...] = tuple(range(1, 13))
    start: Optional[Union[str, pd.Timestamp]] = None
    end: Optional[Union[str, pd.Timestamp]] = None

    def __post_init__(self) -> None:
        if isinstance(self.rate_cfs, pd.Series):
            if not isinstance(self.rate_cfs.index, pd.DatetimeIndex):
                raise ValueError(f"{self.name}: a rate series must be indexed by date")
            if bool((self.rate_cfs.dropna() < 0).any()):
                raise ValueError(f"{self.name}: rate must be >= 0")
        elif not (np.isfinite(float(self.rate_cfs)) and float(self.rate_cfs) >= 0):
            raise ValueError(f"{self.name}: rate must be finite and >= 0, got {self.rate_cfs}")
        if not 0.0 <= float(self.return_fraction) <= 1.0:
            raise ValueError(
                f"{self.name}: return_fraction must be in [0, 1], got {self.return_fraction}"
            )
        if not self.months or any(not 1 <= int(m) <= 12 for m in self.months):
            raise ValueError(f"{self.name}: months must be in 1-12, got {self.months}")
        if self.start is not None and self.end is not None:
            if pd.Timestamp(self.start) > pd.Timestamp(self.end):
                raise ValueError(f"{self.name}: start {self.start} is after end {self.end}")


@dataclass(frozen=True)
class NaturalizedFlow:
    """Result of :func:`naturalize`.

    Attributes
    ----------
    flow : pd.Series
        Naturalized flow, on the gaged series' index.
    added : pd.DataFrame
        The add-back of each diversion (cfs), one column per diversion.
    summary : pd.DataFrame
        Per diversion: ``days_in_effect``, ``days_without_rate`` (in effect but no rate
        for that date) and ``mean_added_cfs`` (over the whole series).
    """

    flow: pd.Series
    added: pd.DataFrame
    summary: pd.DataFrame = field(repr=False)


def _dates(index: pd.DatetimeIndex) -> pd.DatetimeIndex:
    """Calendar dates of an index: local dates for a tz-aware one, time of day dropped."""
    if index.tz is not None:
        index = index.tz_localize(None)
    return index.normalize()


def _in_effect(dates: pd.DatetimeIndex, d: Diversion) -> np.ndarray:
    on = np.asarray(dates.month.isin([int(m) for m in d.months]))
    if d.start is not None:
        on &= np.asarray(dates >= _dates(pd.DatetimeIndex([d.start]))[0])
    if d.end is not None:
        on &= np.asarray(dates <= _dates(pd.DatetimeIndex([d.end]))[0])
    return on


#: Share of the days in effect a rate series may miss before :func:`naturalize` warns.
MISSING_RATE_WARN = 0.25


def naturalize(gaged: pd.Series, diversions: Sequence[Diversion]) -> NaturalizedFlow:
    """Add the consumptive part of each upstream diversion back to a gaged series.

    Parameters
    ----------
    gaged : pd.Series
        Gaged flow (cfs), indexed by date (daily or monthly).
    diversions : sequence of Diversion

    Returns
    -------
    NaturalizedFlow

    Raises
    ------
    TypeError
        ``gaged`` is not indexed by a ``DatetimeIndex``.
    ValueError
        Two diversions share a name.
    """
    if not isinstance(gaged.index, pd.DatetimeIndex):
        raise TypeError("gaged must be indexed by a DatetimeIndex")
    names = [d.name for d in diversions]
    if len(set(names)) != len(names):
        raise ValueError(f"duplicate diversion names: {names}")
    index = pd.DatetimeIndex(gaged.index)
    dates = _dates(index)
    added = pd.DataFrame(index=index)
    rows: List[dict] = []
    for d in diversions:
        on = _in_effect(dates, d)
        if isinstance(d.rate_cfs, pd.Series):
            r = d.rate_cfs.dropna()
            r_dates = _dates(pd.DatetimeIndex(r.index))
            if r_dates.has_duplicates:
                raise ValueError(f"{d.name}: the rate series has duplicate dates")
            by_date = pd.Series(r.to_numpy(dtype=float), index=r_dates)
            rate = by_date.reindex(dates).to_numpy(dtype=float)
        else:
            rate = np.full(len(index), float(d.rate_cfs))
        missing = on & np.isnan(rate)
        add = np.where(on & ~np.isnan(rate), rate * (1.0 - float(d.return_fraction)), 0.0)
        added[d.name] = add
        rows.append(
            {
                "diversion": d.name,
                "days_in_effect": int(on.sum()),
                "days_without_rate": int(missing.sum()),
                "mean_added_cfs": float(np.mean(add)) if len(add) else 0.0,
            }
        )
        if missing.any():
            share = float(missing.sum()) / max(int(on.sum()), 1)
            log = logger.warning if share > MISSING_RATE_WARN else logger.info
            log(
                "%s: no rate on %d of %d day(s) in effect; nothing added there (a rate series "
                "is read by calendar date)",
                d.name,
                int(missing.sum()),
                int(on.sum()),
            )
    flow = gaged.astype(float) + (added.sum(axis=1) if len(diversions) else 0.0)
    summary = pd.DataFrame(
        rows, columns=["diversion", "days_in_effect", "days_without_rate", "mean_added_cfs"]
    ).set_index("diversion")
    flow.name = gaged.name
    return NaturalizedFlow(flow=flow, added=added, summary=summary)
