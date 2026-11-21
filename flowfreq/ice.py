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
    "winter_event_counts",
    "decade_scaled",
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


def winter_event_counts(events: pd.DataFrame, counted_winters: Sequence[int]) -> pd.DataFrame:
    """Count events per winter, with a zero row for every counted winter.

    Parameters
    ----------
    events : pd.DataFrame
        Events with a ``winter`` column (as from :func:`events_in_months`).
    counted_winters : Sequence[int]
        Winters (water years) in which the site was observed.

    Returns
    -------
    pd.DataFrame
        Columns ``winter`` and ``n_events`` (both int64), one row per counted
        winter in the order given; 0 where the winter has no event. Events in
        winters not listed are dropped.
    """
    winters = pd.Series(list(counted_winters), dtype="int64")
    counts = events["winter"].astype("int64").value_counts()
    n = winters.map(counts).fillna(0).astype("int64")
    return pd.DataFrame({"winter": winters, "n_events": n})


def decade_scaled(
    n_events: pd.Series,
    air_index: pd.Series,
    reference_winters: Sequence[int],
) -> dict[str, float]:
    """Scale a site's event count to a reference period by an air index.

    The rate is events per cold spell over the site's counted winters. It is
    applied to the air index of each reference winter: the expected events in
    a winter are ``rate * spells``, and the chance of at least one event is
    ``1 - exp(-rate * spells)`` (a Poisson count).

    Parameters
    ----------
    n_events : pd.Series
        Events per counted winter, indexed by winter. Its index is the set of
        winters the site was observed.
    air_index : pd.Series
        Cold spells (``n_spells``) by winter. A winter missing from it counts
        as 0 spells.
    reference_winters : Sequence[int]
        Winters of the reference period.

    Returns
    -------
    dict[str, float]
        ``n_winters``, ``events_per_winter_obs``, ``share_winters_obs`` (share
        of counted winters with an event), ``rate_per_spell``,
        ``events_per_winter_ref``, ``events_ref_decade``
        (``events_per_winter_ref * 10``, the ten-year period of D&L 38.4),
        ``share_winters_ref`` (mean Poisson chance of an event over the
        reference winters), ``spells_obs`` (spells in the counted winters) and
        ``spells_ref_mean``. When there are no spells in the counted winters
        the rate cannot be formed and the four scaled values are NaN. When
        there are spells but no events they are 0.0.
    """
    site_winters = list(n_events.index)
    events = n_events.to_numpy(dtype=float)
    n_winters = len(site_winters)
    spells_site = air_index.reindex(site_winters).fillna(0).to_numpy(dtype=float)
    spells_ref = air_index.reindex(list(reference_winters)).fillna(0).to_numpy(dtype=float)

    spells_obs = float(spells_site.sum())
    spells_ref_mean = float(spells_ref.mean()) if spells_ref.size else float("nan")
    nan = float("nan")
    if n_winters:
        events_per_winter_obs = float(events.mean())
        share_winters_obs = float((events > 0).mean())
    else:
        events_per_winter_obs = share_winters_obs = nan

    if spells_obs > 0:
        rate = float(events.sum()) / spells_obs
        events_ref = rate * spells_ref_mean
        share_ref = float(np.mean(1.0 - np.exp(-rate * spells_ref))) if spells_ref.size else nan
    else:
        logger.info("no cold spell in the site's winters: rate per spell is undefined")
        rate = events_ref = share_ref = nan

    return {
        "n_winters": n_winters,
        "events_per_winter_obs": events_per_winter_obs,
        "share_winters_obs": share_winters_obs,
        "rate_per_spell": rate,
        "events_per_winter_ref": events_ref,
        "events_ref_decade": events_ref * 10,
        "share_winters_ref": share_ref,
        "spells_obs": spells_obs,
        "spells_ref_mean": spells_ref_mean,
    }


def ice_days_from_qualifiers(daily: pd.DataFrame) -> pd.DataFrame:
    """Flag gauge days that carry an ice or an estimated qualifier.

    The ``qualification_code`` is split on ``":"`` into tokens. ``ICE`` marks
    an ice day. ``e`` or ``ESTIMATED`` marks an estimated day, which is not
    ice: an estimate alone can be an equipment fault. A missing code gives
    False for both.

    Parameters
    ----------
    daily : pd.DataFrame
        Indexed by day (a DatetimeIndex, tz-naive or tz-aware), with a
        ``qualification_code`` column such as ``A:e:ICE`` or ``P:ICE``.

    Returns
    -------
    pd.DataFrame
        Columns ``date`` (datetime64, tz-naive, no time of day), ``ice`` and
        ``estimated`` (bool), one row per input row.
    """
    idx = pd.DatetimeIndex(daily.index)
    if idx.tz is not None:
        idx = idx.tz_localize(None)
    tokens = (
        daily["qualification_code"]
        .astype("object")
        .map(lambda c: set(c.split(":")) if isinstance(c, str) else set())
    )
    return pd.DataFrame(
        {
            "date": pd.Series(idx.normalize(), dtype="datetime64[ns]"),
            "ice": tokens.map(lambda t: "ICE" in t).to_numpy(dtype=bool),
            "estimated": tokens.map(lambda t: "e" in t or "ESTIMATED" in t).to_numpy(dtype=bool),
        }
    )


def gauge_ice_qaqc(
    ice: pd.Series, cold: pd.Series, min_cold_spell_days: int, merge_gap_days: int
) -> pd.DataFrame:
    """Compare gauge ice days with cold air.

    ``ice_without_cold`` is an event of ice days that are not cold-air days,
    of any length. ``cold_without_ice`` is an event of cold-air days that are
    not ice days with at least ``min_cold_spell_days`` days. Events are those
    of :func:`count_events`.

    Parameters
    ----------
    ice, cold : pd.Series
        Bool series indexed by day. They are aligned on the union of their
        indexes, and a day missing from one is False there.
    min_cold_spell_days : int
        Shortest ``cold_without_ice`` event reported.
    merge_gap_days : int
        As in :func:`count_events`.

    Returns
    -------
    pd.DataFrame
        Columns ``kind``, ``start``, ``end`` and ``n_days``, ``ice_without_cold``
        rows first, each kind sorted by start. Empty, with those columns, when
        there is nothing to flag.
    """
    union = ice.index.union(cold.index)
    i = ice.reindex(union).fillna(False).astype(bool)
    c = cold.reindex(union).fillna(False).astype(bool)
    parts = []
    ev = count_events(pd.Series(union[(i & ~c).to_numpy()]), merge_gap_days)
    parts.append(ev.assign(kind="ice_without_cold"))
    ev = count_events(pd.Series(union[(c & ~i).to_numpy()]), merge_gap_days)
    parts.append(ev[ev["n_days"] >= min_cold_spell_days].assign(kind="cold_without_ice"))
    out = pd.concat(parts, ignore_index=True)[["kind", "start", "end", "n_days"]]
    out["n_days"] = out["n_days"].astype("int64")
    return out.reset_index(drop=True)
