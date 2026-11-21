"""Ice-event counting: episodes of flagged days, cold air and the air index.

This module counts episodes of flagged days and rates nothing. It groups flagged
dates into events, flags days that follow cold air, and counts cold spells per
winter. Every threshold is a required argument; none has a default.

A winter is labelled by its water year (``year + (month >= 10)``), so November
2014 through March 2015 is winter 2015.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

__all__ = [
    "winter_of",
    "count_events",
    "cold_air_days",
    "cold_spells_per_winter",
    "events_in_months",
]


def winter_of(dates: pd.Series) -> pd.Series:
    """Water-year label of each date.

    Parameters
    ----------
    dates : pd.Series
        Datetime values.

    Returns
    -------
    pd.Series
        int64 water year, ``year + (month >= 10)``, on the index of ``dates``.
    """
    d = pd.to_datetime(dates)
    return (d.dt.year + (d.dt.month >= 10).astype(int)).astype("int64")


def count_events(dates: pd.Series, merge_gap_days: int) -> pd.DataFrame:
    """Group flagged dates into events.

    Two flagged days are in one event when the number of days between them is at
    most ``merge_gap_days + 1``, that is, when at most ``merge_gap_days`` days are
    missing between them.

    Parameters
    ----------
    dates : pd.Series
        Datetimes of flagged days, in any order; duplicates are ignored.
    merge_gap_days : int
        Largest run of unflagged days that still joins two flagged days.

    Returns
    -------
    pd.DataFrame
        One row per event, sorted by start, with columns ``start`` and ``end``
        (datetime64) and ``n_days`` (int64), the number of flagged days.
    """
    d = pd.Series(pd.to_datetime(dates)).dropna().drop_duplicates().sort_values()
    d = d.reset_index(drop=True)
    if d.empty:
        return pd.DataFrame(
            {
                "start": pd.Series([], dtype="datetime64[ns]"),
                "end": pd.Series([], dtype="datetime64[ns]"),
                "n_days": pd.Series([], dtype="int64"),
            }
        )
    gap_days = d.diff().dt.days
    new_event = (gap_days > merge_gap_days + 1).fillna(True)
    event_id = new_event.cumsum()
    grouped = d.groupby(event_id).agg(start="min", end="max", n_days="count")
    out = pd.DataFrame(grouped).reset_index(drop=True)
    out["n_days"] = out["n_days"].astype("int64")
    return out


def cold_air_days(tmin: pd.Series, threshold_c: float, lookback_days: int) -> pd.Series:
    """Flag days with cold air on the day or in the days before it.

    Parameters
    ----------
    tmin : pd.Series
        Daily minimum air temperature (deg C) on a daily date index.
    threshold_c : float
        A day is cold when ``tmin <= threshold_c``.
    lookback_days : int
        Number of preceding days that also count; a day is flagged when it or
        any of the ``lookback_days`` days before it is cold.

    Returns
    -------
    pd.Series
        bool, on the index of ``tmin``. Missing TMIN counts as not cold.
    """
    cold = (tmin <= threshold_c).fillna(False).astype(float)
    flagged = cold.rolling(lookback_days + 1, min_periods=1).max() > 0
    return flagged.astype(bool)


def cold_spells_per_winter(
    tmin: pd.Series,
    threshold_c: float,
    merge_gap_days: int,
    months: Sequence[int],
) -> pd.DataFrame:
    """Count cold spells and valid air-temperature days per winter.

    A cold day has ``tmin <= threshold_c``. Days outside ``months`` are dropped
    before spells are counted, and a spell is attributed to the winter of its
    first day.

    Parameters
    ----------
    tmin : pd.Series
        Daily minimum air temperature (deg C) indexed by date.
    threshold_c : float
        Cold-day threshold, deg C.
    merge_gap_days : int
        Largest run of non-cold days that still joins two cold days.
    months : Sequence[int]
        Calendar months that make up a winter.

    Returns
    -------
    pd.DataFrame
        Indexed by winter (water year), with ``n_spells`` (int64) and
        ``n_air_days`` (int64, non-missing TMIN days in ``months``).
    """
    idx = pd.DatetimeIndex(tmin.index)
    in_months = idx.month.isin(list(months))
    t = pd.Series(tmin.to_numpy(dtype=float), index=idx)[in_months]
    valid = t.dropna()
    winters = winter_of(pd.Series(valid.index, index=valid.index))
    n_air = valid.groupby(winters.to_numpy()).size().astype("int64")
    n_air.index.name = "winter"

    cold_dates = pd.Series(valid.index[(valid <= threshold_c).to_numpy()])
    ev = count_events(cold_dates, merge_gap_days)
    n_spells: pd.Series
    if ev.empty:
        n_spells = pd.Series([], dtype="int64")
    else:
        n_spells = ev.groupby(winter_of(ev["start"]).to_numpy()).size().astype("int64")
    out = pd.DataFrame({"n_spells": n_spells, "n_air_days": n_air})
    out["n_spells"] = out["n_spells"].fillna(0).astype("int64")
    out["n_air_days"] = out["n_air_days"].fillna(0).astype("int64")
    out.index = pd.Index(np.asarray(out.index, dtype="int64"), name="winter")
    return out.sort_index()


def events_in_months(dates: pd.Series, months: Sequence[int], merge_gap_days: int) -> pd.DataFrame:
    """Count events among flagged dates that fall in the given months.

    Dates outside ``months`` are dropped first, so an event cannot borrow a day
    from outside the winter. Each event is attributed to the winter and month of
    its first day.

    Parameters
    ----------
    dates : pd.Series
        Datetimes of flagged days.
    months : Sequence[int]
        Calendar months to keep.
    merge_gap_days : int
        Largest run of unflagged days that still joins two flagged days.

    Returns
    -------
    pd.DataFrame
        The :func:`count_events` columns plus ``winter`` (water year of the
        start, int64) and ``month`` (calendar month of the start, int64).
    """
    d = pd.Series(pd.to_datetime(dates))
    kept = d[d.dt.month.isin(list(months))]
    ev = count_events(kept, merge_gap_days)
    ev["winter"] = winter_of(ev["start"]).astype("int64")
    ev["month"] = ev["start"].dt.month.astype("int64")
    return ev
