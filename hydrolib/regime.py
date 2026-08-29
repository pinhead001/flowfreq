"""
hydrolib.regime - Flow regime metrics

Intra-annual flow-pattern metrics computed from a daily flow series:
Richards-Baker flashiness, TQmean, baseflow separation and index, and
monthly/seasonal summary statistics -- plus (see the diel-variation
functions further down this module) metrics computed from an instantaneous
series. Grouped together because both describe the *shape* of a flow
record rather than the frequency of its extremes, which is what
hydrolib.bulletin17c and hydrolib.lowflow cover.

Method notes
------------
**Year definition.** Metrics here default to the *water year* (Oct 1 - Sep
30), not the climatic year hydrolib.lowflow uses. The climatic year (Apr 1
- Mar 31) exists specifically to avoid splitting a single connected
low-flow event; flashiness, TQmean, and baseflow index describe a whole
year's hydrograph shape rather than the timing of one extreme, and water
year is both the USGS convention for these indices (e.g. Richards-Baker
flashiness is conventionally reported by water year) and what published
regional statistics use, so it is the more useful default for
comparability. Pass ``year_type="climatic"`` or ``"calendar"`` to use
those instead; see :func:`hydrolib.core.assign_year_label`.

**Completeness.** Every per-year metric in this module uses the same rule:
a year needs at least ``min_days`` valid daily values (negative values
count as missing, matching hydrolib.lowflow) to get a computed value; below
that, the metric is NaN and the row's ``complete`` flag is False rather
than the year disappearing or a value being silently computed from a
partial record. One metric is only approximately exact even in a
"complete" year: the flashiness index sums day-to-day differences, and a
handful of missing days each spoil the (at most two) differences adjacent
to them, which are dropped via ``nansum`` rather than fabricated. With
``min_days=350`` (the default, matching hydrolib.lowflow) this is at most
roughly a dozen dropped difference terms out of ~364, a small but nonzero
approximation worth knowing is there.

**Baseflow separation methods.** Five are implemented, named explicitly
because a BFI value is not comparable across methods without knowing which
one produced it -- every function and results object reports the method
used alongside the number:

- ``"ih_smoothed_minima"`` (default): the UKIH/Institute of Hydrology
  smoothed-minima method (Institute of Hydrology, 1980; as codified by Wahl
  & Wahl, 1988, the most widely cited "Base Flow Index" definition and the
  one most reviewers will recognize by that name). Non-overlapping 5-day
  blocks; a block's minimum is a baseflow ordinate ("turning point") if 0.9
  times it does not exceed either adjacent block's minimum; ordinates are
  connected by straight-line interpolation, capped at actual flow (baseflow
  cannot exceed total streamflow). Needs no drainage area and is
  parameter-light, which is why it is the default here.
- ``"hysep_fixed"``, ``"hysep_sliding"``, ``"hysep_local_minimum"``: the
  three USGS HYSEP methods (Sloto & Crouse, 1996), all built on the same
  interval width ``N* = round_to_odd(2 * drainage_area_sqmi ** 0.2)`` days
  (Pettyjohn & Henning, 1979). Fixed-interval takes the minimum in each
  non-overlapping ``2N*``-day block as a flat baseflow value for that whole
  block (a step function -- the crudest of the three). Sliding-interval is
  a centered rolling minimum of window ``2N*+1`` (smoother, still a direct
  minimum, no interpolation). Local-minimum flags a day as a turning point
  when it equals the minimum of its own centered ``2N*+1`` window, then
  interpolates between turning points exactly as the IH method does. All
  three need ``drainage_area_sqmi`` -- there is no way around this; it is
  intrinsic to how HYSEP defines its window width, so a site missing that
  attribute cannot use these methods, and the error says so rather than
  guessing a drainage area or silently falling back to another method.
- ``"lyne_hollick"``: the Lyne & Hollick (1979) recursive digital filter,
  3 passes (forward, backward, forward) per Nathan & McMahon (1990),
  alpha=0.925. A single missing day would otherwise poison every
  subsequent value through the filter's own recursion (each day depends on
  the previous filtered value), so this implementation runs the filter
  separately on each gap-free run of the input and leaves gaps as NaN
  rather than bridging or resetting across them silently.

**Period-of-record summary.** :meth:`FlowRegime.summary` aggregates the
annual series two different ways, because they are not the same kind of
quantity: baseflow index and flashiness index are each a ratio of two
summable quantities (baseflow volume over total volume; total day-to-day
change over total volume), so their period-of-record value is computed by
pooling the numerator and denominator across every complete year and
taking one ratio -- volume-weighted, so a high-flow year is not given the
same influence as a low-flow year with far less water in it, which a naive
mean of per-year ratios would do. TQmean has no such pooled form, because
its threshold (that year's own mean daily flow) is itself year-relative; its
period-of-record value is the arithmetic mean of the per-year TQmean
values, which is what the term conventionally means.
"""

from __future__ import annotations

import calendar
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

from .core import YEAR_TYPES, assign_year_label

#: Baseflow separation methods accepted by separate_baseflow, baseflow_index,
#: and FlowRegime.
BASEFLOW_METHODS = (
    "ih_smoothed_minima",
    "hysep_fixed",
    "hysep_sliding",
    "hysep_local_minimum",
    "lyne_hollick",
)

#: Methods in BASEFLOW_METHODS that require drainage_area_sqmi.
_HYSEP_METHODS = ("hysep_fixed", "hysep_sliding", "hysep_local_minimum")

#: Standard meteorological seasons: month -> (season name, year offset).
#: The year offset handles December, which belongs to the following
#: meteorological winter (e.g. Dec 2019 is part of "winter 2020").
_SEASON_BY_MONTH: Dict[int, Tuple[str, int]] = {
    12: ("winter", 1),
    1: ("winter", 0),
    2: ("winter", 0),
    3: ("spring", 0),
    4: ("spring", 0),
    5: ("spring", 0),
    6: ("summer", 0),
    7: ("summer", 0),
    8: ("summer", 0),
    9: ("fall", 0),
    10: ("fall", 0),
    11: ("fall", 0),
}

#: Calendar months making up each standard meteorological season.
_SEASON_MONTHS: Dict[str, Tuple[int, int, int]] = {
    "winter": (12, 1, 2),
    "spring": (3, 4, 5),
    "summer": (6, 7, 8),
    "fall": (9, 10, 11),
}


def _days_in_season(year: int, season: str) -> int:
    """True calendar day count for a (year, season) pair, e.g. "winter 2020"
    (Dec 2019 + Jan 2020 + Feb 2020) = 31 + 31 + 29 = 91.

    Computed independent of any particular dataset's date range -- unlike
    counting rows in a reindexed frame, which only reflects the overlap
    between the data actually provided and the season, and silently
    undercounts a season's true length whenever the input doesn't happen to
    span that season's full calendar extent (always true for the first and
    last season touched by any real, finite dataset).
    """
    total = 0
    for month in _SEASON_MONTHS[season]:
        # December belongs to the following winter's year label.
        cal_year = year - 1 if month == 12 else year
        total += calendar.monthrange(cal_year, month)[1]
    return total


# =============================================================================
# Shared helpers
# =============================================================================


def _daily_with_year(daily_data: pd.DataFrame, year_type: str) -> pd.DataFrame:
    """Reindex onto a complete daily calendar and attach a year label.

    Negative values are set to NaN (a data artifact, not a legitimate low --
    matching hydrolib.lowflow.annual_minimum_flow's convention), and a real
    gap becomes an explicit NaN row rather than a missing one, since
    :meth:`hydrolib.usgs.USGSgage.download_daily_flow` drops rows with no
    value.
    """
    if year_type not in YEAR_TYPES:
        raise ValueError(f"year_type must be one of {YEAR_TYPES}, got {year_type!r}")

    full_idx = pd.date_range(daily_data.index.min(), daily_data.index.max(), freq="D")
    flows = daily_data["flow_cfs"].reindex(full_idx)
    flows = flows.where(flows >= 0)
    years = assign_year_label(full_idx, year_type)
    return pd.DataFrame({"year": years, "flow_cfs": flows.to_numpy()}, index=full_idx)


def _year_completeness(df: pd.DataFrame, min_days: int) -> Tuple[pd.Series, pd.Series]:
    """Per-year valid-day counts and the min_days completeness flag, as a
    (n_days, complete) pair sharing one groupby pass."""
    n_days = df.groupby("year")["flow_cfs"].apply(lambda s: int(s.notna().sum()))
    return n_days, n_days >= min_days


# =============================================================================
# Richards-Baker Flashiness Index
# =============================================================================


def richards_baker_flashiness(
    daily_data: pd.DataFrame, year_type: str = "water", min_days: int = 350
) -> pd.DataFrame:
    """Richards-Baker Flashiness Index (Baker et al., 2004) per year.

    RBI = sum(abs(q[i] - q[i-1])) / sum(q), using the n-1 within-year
    day-to-day differences over n days (the day-to-day change across a year
    boundary is not included -- each year is scored independently).

    Parameters
    ----------
    daily_data : pd.DataFrame
        Daily mean flow indexed by date, with a ``flow_cfs`` column.
    year_type : str
        Year definition; default "water". See the module docstring.
    min_days : int
        Per-year completeness threshold. Default 350.

    Returns
    -------
    pd.DataFrame
        Columns: ``year``, ``flashiness_index`` (NaN if incomplete),
        ``n_days``, ``complete``.
    """
    df = _daily_with_year(daily_data, year_type)

    def _rbi(group: pd.DataFrame) -> float:
        vals = group["flow_cfs"].to_numpy()
        total = np.nansum(vals)
        if total <= 0 or np.all(np.isnan(vals)):
            return np.nan
        diffs = np.abs(np.diff(vals))
        return float(np.nansum(diffs) / total)

    n_days, complete = _year_completeness(df, min_days)
    rbi = df.groupby("year").apply(_rbi, include_groups=False)
    rbi = rbi.where(complete, np.nan)

    result = pd.DataFrame(
        {"flashiness_index": rbi, "n_days": n_days, "complete": complete}
    ).reset_index(names="year")
    return result.sort_values("year").reset_index(drop=True)


# =============================================================================
# TQmean
# =============================================================================


def tqmean(daily_data: pd.DataFrame, year_type: str = "water", min_days: int = 350) -> pd.DataFrame:
    """TQmean: fraction of the year that flow exceeds that year's own mean.

    Parameters
    ----------
    daily_data : pd.DataFrame
        Daily mean flow indexed by date, with a ``flow_cfs`` column.
    year_type : str
        Year definition; default "water".
    min_days : int
        Per-year completeness threshold. Default 350.

    Returns
    -------
    pd.DataFrame
        Columns: ``year``, ``tqmean`` (NaN if incomplete), ``n_days``,
        ``complete``.
    """
    df = _daily_with_year(daily_data, year_type)

    def _tqmean(group: pd.DataFrame) -> float:
        vals = group["flow_cfs"].to_numpy()
        valid = vals[~np.isnan(vals)]
        if len(valid) == 0:
            return np.nan
        return float(np.mean(valid > valid.mean()))

    n_days, complete = _year_completeness(df, min_days)
    tq = df.groupby("year").apply(_tqmean, include_groups=False)
    tq = tq.where(complete, np.nan)

    result = pd.DataFrame({"tqmean": tq, "n_days": n_days, "complete": complete}).reset_index(
        names="year"
    )
    return result.sort_values("year").reset_index(drop=True)


# =============================================================================
# Baseflow separation
# =============================================================================


def _interpolate_baseflow(
    flows: np.ndarray, positions: np.ndarray, values: np.ndarray
) -> np.ndarray:
    """Straight-line interpolation between baseflow turning points, capped at
    actual flow (baseflow cannot exceed total streamflow). NaN outside the
    span bounded by the first and last turning point, and wherever `flows`
    itself is NaN. Shared by the IH and HYSEP-local-minimum methods, which
    differ only in how a turning point is identified.
    """
    n = len(flows)
    baseflow = np.full(n, np.nan)
    if len(positions) < 2:
        return baseflow

    lo, hi = int(positions[0]), int(positions[-1])
    idx = np.arange(lo, hi + 1)
    interpolated = np.interp(idx, positions, values)
    segment = flows[lo : hi + 1]
    valid = ~np.isnan(segment)
    baseflow[lo : hi + 1] = np.where(valid, np.minimum(interpolated, segment), np.nan)
    return baseflow


def _ih_smoothed_minima(flows: np.ndarray, block_days: int = 5, factor: float = 0.9) -> np.ndarray:
    """UKIH/Institute of Hydrology smoothed-minima baseflow separation.

    Parameters
    ----------
    flows : np.ndarray
    block_days : int
        Non-overlapping block length, default 5 (the standard UKIH value).
    factor : float
        Turning-point ratio, default 0.9 (the standard UKIH value): a
        block's minimum is a turning point if ``factor * block_min <=
        min(previous block's min, next block's min)``.
    """
    n = len(flows)
    n_blocks = n // block_days
    block_mins = np.full(n_blocks, np.nan)
    block_pos = np.full(n_blocks, -1, dtype=int)

    for b in range(n_blocks):
        seg = flows[b * block_days : (b + 1) * block_days]
        if np.isnan(seg).any():
            continue
        i = int(np.argmin(seg))
        block_mins[b] = seg[i]
        block_pos[b] = b * block_days + i

    turning_pos, turning_val = [], []
    for i in range(1, n_blocks - 1):
        if np.isnan(block_mins[i - 1]) or np.isnan(block_mins[i]) or np.isnan(block_mins[i + 1]):
            continue
        if factor * block_mins[i] <= min(block_mins[i - 1], block_mins[i + 1]):
            turning_pos.append(block_pos[i])
            turning_val.append(block_mins[i])

    return _interpolate_baseflow(flows, np.array(turning_pos), np.array(turning_val))


def _hysep_interval_days(drainage_area_sqmi: float) -> int:
    """HYSEP's interval half-width N*, in days (Pettyjohn & Henning, 1979).

    N* = 2 * A^0.2, rounded to the nearest odd integer, A in square miles.
    """
    if drainage_area_sqmi <= 0:
        raise ValueError(f"drainage_area_sqmi must be positive, got {drainage_area_sqmi}")
    n_star = 2.0 * drainage_area_sqmi**0.2
    n_odd = int(2 * round((n_star - 1) / 2) + 1)
    return max(n_odd, 1)


def _hysep_fixed(flows: np.ndarray, n_star: int) -> np.ndarray:
    """HYSEP fixed-interval method: minimum of each non-overlapping 2*N*-day
    block, held constant across the block (a step function)."""
    width = 2 * n_star
    n = len(flows)
    baseflow = np.full(n, np.nan)
    for start in range(0, n, width):
        end = min(start + width, n)
        seg = flows[start:end]
        if np.isnan(seg).any():
            continue
        baseflow[start:end] = seg.min()
    return baseflow


def _hysep_sliding(flows: np.ndarray, n_star: int) -> np.ndarray:
    """HYSEP sliding-interval method: centered rolling minimum, window
    2*N*+1 days."""
    window = 2 * n_star + 1
    return pd.Series(flows).rolling(window=window, center=True, min_periods=window).min().to_numpy()


def _hysep_local_minimum(flows: np.ndarray, n_star: int) -> np.ndarray:
    """HYSEP local-minimum method: a day is a turning point if it equals the
    minimum of its own centered 2*N*+1-day window; turning points are then
    interpolated exactly as the IH method does."""
    window = 2 * n_star + 1
    s = pd.Series(flows)
    roll_min = s.rolling(window=window, center=True, min_periods=window).min()
    is_turning = (s == roll_min) & roll_min.notna()
    positions = np.where(is_turning.to_numpy())[0]
    values = flows[positions]
    return _interpolate_baseflow(flows, positions, values)


def _gap_free_runs(values: np.ndarray) -> list:
    """Maximal contiguous non-NaN index ranges, as half-open (start, end) pairs."""
    is_valid = ~np.isnan(values)
    runs = []
    start = None
    for i, valid in enumerate(is_valid):
        if valid and start is None:
            start = i
        elif not valid and start is not None:
            runs.append((start, i))
            start = None
    if start is not None:
        runs.append((start, len(values)))
    return runs


def _lyne_hollick_pass(x: np.ndarray, alpha: float, forward: bool) -> np.ndarray:
    """One forward or backward pass of the Lyne-Hollick recursive filter."""
    y = x if forward else x[::-1]
    n = len(y)
    quickflow = np.zeros(n)
    baseflow = np.empty(n)
    baseflow[0] = y[0]
    for i in range(1, n):
        quickflow[i] = alpha * quickflow[i - 1] + ((1 + alpha) / 2) * (y[i] - y[i - 1])
        b = y[i] - quickflow[i]
        b = min(b, y[i])
        b = max(b, 0.0)
        baseflow[i] = b
    return baseflow if forward else baseflow[::-1]


def _lyne_hollick(flows: np.ndarray, alpha: float = 0.925, passes: int = 3) -> np.ndarray:
    """Lyne-Hollick recursive digital filter, run separately on each gap-free
    run of the input.

    A missing day would otherwise poison every value after it: each day's
    filtered value depends on the previous day's, so one NaN propagates
    through the entire remainder of the recursion (and, with backward
    passes, potentially the portion before it too). Splitting on gap-free
    runs contains that to exactly the missing day(s) themselves.
    """
    result = np.full(len(flows), np.nan)
    for start, end in _gap_free_runs(flows):
        segment = flows[start:end].astype(float)
        b = segment.copy()
        forward = True
        for _ in range(passes):
            b = _lyne_hollick_pass(b, alpha, forward)
            forward = not forward
        result[start:end] = b
    return result


def separate_baseflow(
    daily_data: pd.DataFrame,
    method: str = "ih_smoothed_minima",
    drainage_area_sqmi: Optional[float] = None,
    block_days: int = 5,
    factor: float = 0.9,
    alpha: float = 0.925,
    passes: int = 3,
) -> pd.Series:
    """Separate a daily flow series into baseflow and quickflow.

    See the module docstring for what each method in :data:`BASEFLOW_METHODS`
    does and why ``"ih_smoothed_minima"`` is the default.

    Parameters
    ----------
    daily_data : pd.DataFrame
        Daily mean flow indexed by date, with a ``flow_cfs`` column.
    method : str
        One of :data:`BASEFLOW_METHODS`.
    drainage_area_sqmi : float, optional
        Required for the three ``"hysep_*"`` methods; unused otherwise.
    block_days, factor : optional
        Only used by ``"ih_smoothed_minima"``; see :func:`_ih_smoothed_minima`.
    alpha, passes : optional
        Only used by ``"lyne_hollick"``; see :func:`_lyne_hollick`.

    Returns
    -------
    pd.Series
        Baseflow in cfs, indexed on the complete daily calendar spanning
        `daily_data` (gaps filled as NaN -- see
        :func:`hydrolib.lowflow.annual_minimum_flow`'s Notes on why a real
        gap must not be silently bridged). NaN outside the span any given
        method can support (e.g. before the first identified turning
        point).

    Raises
    ------
    ValueError
        `method` is not recognized, or a ``"hysep_*"`` method was requested
        without `drainage_area_sqmi`.
    """
    if method not in BASEFLOW_METHODS:
        raise ValueError(f"method must be one of {BASEFLOW_METHODS}, got {method!r}")
    if method in _HYSEP_METHODS and drainage_area_sqmi is None:
        raise ValueError(
            f"method={method!r} requires drainage_area_sqmi -- HYSEP's window width is "
            f"defined in terms of it (N* = round_to_odd(2 * A**0.2)), so there is no "
            f"drainage-area-free way to run it. Pass it explicitly, or use "
            f"method='ih_smoothed_minima' (the default) or 'lyne_hollick', neither of "
            f"which needs it."
        )

    full_idx = pd.date_range(daily_data.index.min(), daily_data.index.max(), freq="D")
    flows = daily_data["flow_cfs"].reindex(full_idx).where(lambda s: s >= 0).to_numpy(dtype=float)

    if method == "ih_smoothed_minima":
        baseflow = _ih_smoothed_minima(flows, block_days=block_days, factor=factor)
    elif method == "lyne_hollick":
        baseflow = _lyne_hollick(flows, alpha=alpha, passes=passes)
    else:
        n_star = _hysep_interval_days(drainage_area_sqmi)
        if method == "hysep_fixed":
            baseflow = _hysep_fixed(flows, n_star)
        elif method == "hysep_sliding":
            baseflow = _hysep_sliding(flows, n_star)
        else:
            baseflow = _hysep_local_minimum(flows, n_star)

    return pd.Series(baseflow, index=full_idx, name="baseflow_cfs")


def baseflow_index(
    daily_data: pd.DataFrame,
    year_type: str = "water",
    min_days: int = 350,
    method: str = "ih_smoothed_minima",
    drainage_area_sqmi: Optional[float] = None,
    **method_kwargs,
) -> pd.DataFrame:
    """Baseflow index (BFI) per year: sum(baseflow) / sum(total flow).

    Parameters
    ----------
    daily_data : pd.DataFrame
        Daily mean flow indexed by date, with a ``flow_cfs`` column.
    year_type : str
        Year definition; default "water".
    min_days : int
        Per-year completeness threshold. Default 350.
    method : str
        Baseflow separation method; see :func:`separate_baseflow` and the
        module docstring. Default "ih_smoothed_minima".
    drainage_area_sqmi : float, optional
        Required for the ``"hysep_*"`` methods.
    **method_kwargs
        Passed through to :func:`separate_baseflow` (``block_days``,
        ``factor``, ``alpha``, ``passes``).

    Returns
    -------
    pd.DataFrame
        Columns: ``year``, ``baseflow_index`` (NaN if incomplete or if the
        method found no valid baseflow ordinate in that year), ``n_days``,
        ``complete``, ``method``.
    """
    baseflow = separate_baseflow(
        daily_data, method=method, drainage_area_sqmi=drainage_area_sqmi, **method_kwargs
    )
    df = _daily_with_year(daily_data.reindex(baseflow.index), year_type)
    df["baseflow_cfs"] = baseflow.to_numpy()

    def _bfi(group: pd.DataFrame) -> float:
        total = np.nansum(group["flow_cfs"].to_numpy())
        base = np.nansum(group["baseflow_cfs"].to_numpy())
        if total <= 0 or group["baseflow_cfs"].notna().sum() == 0:
            return np.nan
        return float(base / total)

    n_days, complete = _year_completeness(df, min_days)
    bfi = df.groupby("year").apply(_bfi, include_groups=False)
    bfi = bfi.where(complete, np.nan)

    result = pd.DataFrame(
        {"baseflow_index": bfi, "n_days": n_days, "complete": complete}
    ).reset_index(names="year")
    result["method"] = method
    return result.sort_values("year").reset_index(drop=True)


# =============================================================================
# Monthly and seasonal summary statistics
# =============================================================================


def monthly_flow_summary(
    daily_data: pd.DataFrame, year_type: str = "water", min_completeness_frac: float = 0.9
) -> pd.DataFrame:
    """Per-year, per-calendar-month flow summary statistics.

    Distinct from :meth:`hydrolib.hydrograph.Hydrograph.get_summary_stats`,
    which pools every year together into one climatological day-of-water-year
    shape; this returns one row per (year, month) so month-to-month and
    year-to-year variation both stay visible for trend analysis.

    Parameters
    ----------
    daily_data : pd.DataFrame
        Daily mean flow indexed by date, with a ``flow_cfs`` column.
    year_type : str
        Year definition used for the ``year`` column; default "water".
        Calendar month is unaffected by this choice.
    min_completeness_frac : float
        Fraction of a month's calendar days that must have a valid value
        for that (year, month) row to be marked ``complete``. Default 0.9
        (at most ~3 missing days in a 31-day month). A day count rather
        than a fraction is not used here, unlike the annual ``min_days``
        gate elsewhere in this module, because months are not all the same
        length.

    Returns
    -------
    pd.DataFrame
        Columns: ``year``, ``month`` (1-12), ``mean_flow_cfs``,
        ``min_flow_cfs``, ``max_flow_cfs``, ``median_flow_cfs``, ``n_days``,
        ``days_in_month``, ``complete``.
    """
    full_idx: pd.DatetimeIndex = pd.date_range(
        daily_data.index.min(), daily_data.index.max(), freq="D"
    )
    flows = daily_data["flow_cfs"].reindex(full_idx).where(lambda s: s >= 0)
    years = assign_year_label(full_idx, year_type)
    df = pd.DataFrame(
        {
            "year": years,
            "month": full_idx.month,
            "calendar_year": full_idx.year,
            "flow_cfs": flows.to_numpy(),
        },
        index=full_idx,
    )

    grouped = df.groupby(["year", "month"])["flow_cfs"]
    summary = grouped.agg(mean_flow_cfs="mean", min_flow_cfs="min", max_flow_cfs="max")
    summary["median_flow_cfs"] = grouped.median()
    summary["n_days"] = grouped.apply(lambda s: int(s.notna().sum()))
    # A (year_label, month) group's calendar year is NOT always year_label:
    # under year_type="climatic", a January-March row is labeled by the
    # PREVIOUS calendar year (climatic year Y spans Apr Y - Mar Y+1), so
    # using year_label directly here would ask calendar.monthrange for the
    # wrong year and silently get February's leap-year day count wrong.
    # Every date within one (year_label, month) group shares the same true
    # calendar year, so read it from the data instead of the label.
    calendar_year_by_group = df.groupby(["year", "month"])["calendar_year"].first()
    summary["days_in_month"] = [
        calendar.monthrange(int(cy), int(month))[1]
        for (_, month), cy in zip(summary.index, calendar_year_by_group)
    ]
    summary["complete"] = summary["n_days"] >= (summary["days_in_month"] * min_completeness_frac)

    return summary.reset_index().sort_values(["year", "month"]).reset_index(drop=True)


def seasonal_flow_summary(
    daily_data: pd.DataFrame, min_completeness_frac: float = 0.9
) -> pd.DataFrame:
    """Per standard meteorological season flow summary statistics.

    Uses the fixed meteorological definition (winter=DJF, spring=MAM,
    summer=JJA, fall=SON), labeled by the season's own conventional year
    (e.g. December belongs to the following winter: Dec 2019 is part of
    "winter 2020"). This is intentionally independent of any ``year_type``
    argument -- unlike the annual and monthly functions in this module,
    which accept water/climatic/calendar year because those describe a full
    year's worth of data, meteorological seasons are a fixed, universally
    recognized definition, and a "water-year-aligned" reinterpretation of
    what "winter" means would need its own explanation and risks being
    misread as the standard definition. A given season's year label can
    therefore differ from the water-year label attached to the same days
    elsewhere in this module (October, for instance, is "fall" of its own
    calendar year here, but falls inside the water year numbered one
    higher) -- this is expected, not an inconsistency to reconcile.

    Parameters
    ----------
    daily_data : pd.DataFrame
        Daily mean flow indexed by date, with a ``flow_cfs`` column.
    min_completeness_frac : float
        Fraction of a season's calendar days that must have a valid value
        for that (year, season) row to be marked ``complete``. Default 0.9.

    Returns
    -------
    pd.DataFrame
        Columns: ``year``, ``season`` ("winter"/"spring"/"summer"/"fall"),
        ``mean_flow_cfs``, ``min_flow_cfs``, ``max_flow_cfs``,
        ``median_flow_cfs``, ``n_days``, ``days_in_season``, ``complete``.
    """
    full_idx: pd.DatetimeIndex = pd.date_range(
        daily_data.index.min(), daily_data.index.max(), freq="D"
    )
    flows = daily_data["flow_cfs"].reindex(full_idx).where(lambda s: s >= 0)

    season_names = np.empty(len(full_idx), dtype=object)
    season_years = np.empty(len(full_idx), dtype=int)
    for month, (name, offset) in _SEASON_BY_MONTH.items():
        mask = full_idx.month == month
        season_names[mask] = name
        season_years[mask] = full_idx.year[mask] + offset

    df = pd.DataFrame(
        {"year": season_years, "season": season_names, "flow_cfs": flows.to_numpy()},
        index=full_idx,
    )

    grouped = df.groupby(["year", "season"])["flow_cfs"]
    summary = grouped.agg(mean_flow_cfs="mean", min_flow_cfs="min", max_flow_cfs="max")
    summary["median_flow_cfs"] = grouped.median()
    summary["n_days"] = grouped.apply(lambda s: int(s.notna().sum()))
    # True calendar length of the season, e.g. winter 2020 (Dec2019+Jan2020+
    # Feb2020) = 91 days -- NOT df.groupby(...).size(), which only reflects
    # the overlap between the season and whatever date range this dataset
    # happens to span, and silently undercounts for the first/last season
    # touched by any real, finite record (exactly the practically common
    # case). Same class of bug as monthly_flow_summary's days_in_month.
    summary["days_in_season"] = [_days_in_season(year, season) for year, season in summary.index]
    summary["complete"] = summary["n_days"] >= (summary["days_in_season"] * min_completeness_frac)

    season_order = {"winter": 0, "spring": 1, "summer": 2, "fall": 3}
    result = summary.reset_index()
    result["_order"] = result["season"].map(season_order)
    return result.sort_values(["year", "_order"]).drop(columns="_order").reset_index(drop=True)


# =============================================================================
# Facade
# =============================================================================


class FlowRegime:
    """
    Flow regime metrics computed from a daily flow series: Richards-Baker
    flashiness, TQmean, baseflow index, and monthly/seasonal summary
    statistics -- each as a per-year (or per-year-month/season) table, plus
    a period-of-record summary. See the module docstring for method
    definitions and the year-type/completeness conventions shared across
    all of them.

    Parameters
    ----------
    daily_data : pd.DataFrame
        Daily mean flow indexed by date, with a ``flow_cfs`` column.
    year_type : str
        Year definition for the annual and monthly tables; default "water".
        Does not affect :attr:`seasonal`, which always uses standard
        meteorological seasons -- see :func:`seasonal_flow_summary`.
    min_days : int
        Per-year completeness threshold for the annual table. Default 350.
    min_completeness_frac : float
        Per-month/season completeness threshold for the monthly and
        seasonal tables. Default 0.9.
    baseflow_method : str
        One of :data:`BASEFLOW_METHODS`. Default "ih_smoothed_minima".
    drainage_area_sqmi : float, optional
        Required if `baseflow_method` is one of the ``"hysep_*"`` methods.
    baseflow_block_days, baseflow_factor : optional
        Passed through to :func:`separate_baseflow` for the
        "ih_smoothed_minima" method.
    baseflow_alpha, baseflow_passes : optional
        Passed through to :func:`separate_baseflow` for the "lyne_hollick"
        method.

    Examples
    --------
    >>> regime = FlowRegime(daily_df, year_type="water")
    >>> regime.annual  # doctest: +SKIP
    >>> regime.summary()  # doctest: +SKIP
    """

    def __init__(
        self,
        daily_data: pd.DataFrame,
        year_type: str = "water",
        min_days: int = 350,
        min_completeness_frac: float = 0.9,
        baseflow_method: str = "ih_smoothed_minima",
        drainage_area_sqmi: Optional[float] = None,
        baseflow_block_days: int = 5,
        baseflow_factor: float = 0.9,
        baseflow_alpha: float = 0.925,
        baseflow_passes: int = 3,
    ):
        self._year_type = year_type
        self._min_days = min_days
        self._baseflow_method = baseflow_method
        # Reindexed daily flow with year labels -- the same shape every
        # per-year function in this module builds internally, kept here so
        # summary() can pool across years without re-deriving it or, worse,
        # re-deriving it slightly differently.
        self._daily = _daily_with_year(daily_data, year_type)

        rbi = richards_baker_flashiness(daily_data, year_type=year_type, min_days=min_days)
        tq = tqmean(daily_data, year_type=year_type, min_days=min_days)
        bfi = baseflow_index(
            daily_data,
            year_type=year_type,
            min_days=min_days,
            method=baseflow_method,
            drainage_area_sqmi=drainage_area_sqmi,
            block_days=baseflow_block_days,
            factor=baseflow_factor,
            alpha=baseflow_alpha,
            passes=baseflow_passes,
        )

        annual = rbi[["year", "flashiness_index", "n_days", "complete"]].merge(
            tq[["year", "tqmean"]], on="year"
        )
        annual = annual.merge(bfi[["year", "baseflow_index"]], on="year")
        self._annual = annual[
            ["year", "flashiness_index", "tqmean", "baseflow_index", "n_days", "complete"]
        ]

        self._monthly = monthly_flow_summary(
            daily_data, year_type=year_type, min_completeness_frac=min_completeness_frac
        )
        self._seasonal = seasonal_flow_summary(
            daily_data, min_completeness_frac=min_completeness_frac
        )
        self._baseflow_series = separate_baseflow(
            daily_data,
            method=baseflow_method,
            drainage_area_sqmi=drainage_area_sqmi,
            block_days=baseflow_block_days,
            factor=baseflow_factor,
            alpha=baseflow_alpha,
            passes=baseflow_passes,
        )

    @property
    def annual(self) -> pd.DataFrame:
        """Per-year table: flashiness_index, tqmean, baseflow_index, n_days, complete."""
        return self._annual.copy()

    @property
    def monthly(self) -> pd.DataFrame:
        """Per (year, month) summary statistics table."""
        return self._monthly.copy()

    @property
    def seasonal(self) -> pd.DataFrame:
        """Per (year, season) summary statistics table (standard meteorological seasons)."""
        return self._seasonal.copy()

    @property
    def baseflow_series(self) -> pd.Series:
        """The continuous daily baseflow series from the separation method used."""
        return self._baseflow_series.copy()

    def summary(self) -> pd.Series:
        """Period-of-record summary across complete years.

        See the module docstring for why baseflow_index and flashiness_index
        are pooled (volume-weighted) while tqmean is averaged.

        Returns
        -------
        pd.Series
            Index: ``n_years``, ``flashiness_index``, ``tqmean``,
            ``baseflow_index``.
        """
        complete = self._annual[self._annual["complete"]]
        n_years = len(complete)
        if n_years == 0:
            return pd.Series(
                {
                    "n_years": 0,
                    "flashiness_index": np.nan,
                    "tqmean": np.nan,
                    "baseflow_index": np.nan,
                }
            )

        complete_years = set(complete["year"])
        rbi_pooled = self._pooled_flashiness(complete_years)
        bfi_pooled = self._pooled_baseflow_index(complete_years)
        tqmean_mean = float(self._annual.loc[self._annual["complete"], "tqmean"].mean())

        return pd.Series(
            {
                "n_years": n_years,
                "flashiness_index": rbi_pooled,
                "tqmean": tqmean_mean,
                "baseflow_index": bfi_pooled,
            }
        )

    def _pooled_flashiness(self, complete_years: set) -> float:
        """Volume-weighted flashiness across complete years: pool every
        complete year's own diff-sum and flow-sum, then take one ratio."""
        total_diff = 0.0
        total_flow = 0.0
        for year, group in self._daily.groupby("year"):
            if year not in complete_years:
                continue
            vals = group["flow_cfs"].to_numpy()
            total_diff += np.nansum(np.abs(np.diff(vals)))
            total_flow += np.nansum(vals)
        return float(total_diff / total_flow) if total_flow > 0 else np.nan

    def _pooled_baseflow_index(self, complete_years: set) -> float:
        """Volume-weighted BFI across complete years: pool baseflow and
        total-flow volume, then take one ratio."""
        baseflow = self._baseflow_series.reindex(self._daily.index).to_numpy()
        total_base = 0.0
        total_flow = 0.0
        for year, group in self._daily.assign(baseflow_cfs=baseflow).groupby("year"):
            if year not in complete_years:
                continue
            total_base += np.nansum(group["baseflow_cfs"].to_numpy())
            total_flow += np.nansum(group["flow_cfs"].to_numpy())
        return float(total_base / total_flow) if total_flow > 0 else np.nan


# =============================================================================
# Diel variation
# =============================================================================


def diel_variation(
    iv_data: pd.DataFrame, tz: str, min_completeness_frac: float = 0.9
) -> pd.DataFrame:
    """Per-day diel (within-day) variation statistics from an instantaneous series.

    Computes daily range, coefficient of variation, and related summary
    statistics for each local calendar day. This is descriptive only: it
    quantifies how much a stream's flow swings within a day, not why. A
    snowmelt-driven diel signal (afternoon melt pulse, overnight recession)
    and a hydropeaking or diversion-driven signal can produce a similar
    range or CV, and telling them apart requires knowledge of upstream
    operations (reservoirs, diversions, hydropower schedules) that this
    function has no way to know and does not attempt to infer. Treat every
    number this returns as "how much flow varied that day," not as an
    automatic natural/regulated classification.

    Parameters
    ----------
    iv_data : pd.DataFrame
        Instantaneous flow with a tz-aware datetime index and a ``flow_cfs``
        column -- the shape returned by
        :meth:`hydrolib.usgs.USGSgage.download_instantaneous_flow`.
    tz : str
        IANA time zone the data should be grouped by calendar day in, e.g.
        ``"America/Los_Angeles"``. Required, with no default: grouping by
        the index's own zone (typically UTC, per
        ``download_instantaneous_flow``'s default) instead of the gage's
        local zone silently fractures each local day across two UTC-labeled
        buckets, corrupting the range and CV of both -- this is not a
        hypothetical edge case, it happens on every single day for any
        site west or east of the UTC meridian. If `iv_data` is already in
        the zone you want, pass that same zone here; the conversion is then
        a no-op.
    min_completeness_frac : float
        Fraction of the day's *expected* observation count (inferred from
        the data's own median sampling interval; see Notes) that must be
        present for a day to be marked ``complete``. Default 0.9.

    Returns
    -------
    pd.DataFrame
        One row per local calendar day with at least one observation.
        Columns:

        - ``date`` : the local calendar date
        - ``min_flow_cfs``, ``max_flow_cfs``, ``mean_flow_cfs``,
          ``std_flow_cfs`` : float
        - ``range_cfs`` : ``max_flow_cfs - min_flow_cfs``, the day's diel
          amplitude
        - ``cv`` : ``std_flow_cfs / mean_flow_cfs``, NaN if the day's mean
          flow is not positive (a day that includes a zero or negative-
          treated-as-missing mean is not a meaningful denominator; see
          hydrolib.lowflow's module docstring on the same issue for annual
          minima)
        - ``n_obs``, ``expected_obs``, ``complete``

    Notes
    -----
    **Range and CV are reported regardless of completeness.** This is a
    deliberate departure from the annual/monthly/seasonal metrics earlier
    in this module, which report NaN for an incomplete period. A day
    missing some fraction of its readings (a logger gap, a brief outage)
    still usually captures most of the diurnal cycle, so its range and CV
    remain informative, just less precise, at partial coverage -- unlike an
    annual metric computed from a year missing an entire season. `complete`
    is there to tell you how much to trust the day's precision, not to
    gate whether a value is reported at all. A day with exactly one
    observation still gets ``range_cfs = 0`` (mathematically correct: a
    single point has no spread) but ``cv = NaN`` (a standard deviation needs
    at least two points) -- both are reported as computed, with `complete`
    correctly False.

    **Expected observations per day** are inferred from the *median* time
    step across the entire input, not a hardcoded assumption like 15
    minutes -- NWIS's most common instantaneous interval, but not the only
    one a logger might report at. If the sampling interval genuinely
    changes partway through the record (e.g. a logger upgrade from hourly
    to 15-minute reporting), this single global median will misjudge
    completeness for whichever era doesn't match it; split the record at
    the change and call this function on each piece separately in that
    case.

    **Negative values** are treated as missing, matching every other daily
    function in this library -- a data artifact, not a legitimate reading.

    Raises
    ------
    ValueError
        Fewer than two timestamps are available to infer a sampling
        interval from.
    """
    if len(iv_data) < 2:
        raise ValueError(
            "diel_variation needs at least 2 timestamps to infer a sampling interval; "
            f"got {len(iv_data)}"
        )

    local_index = iv_data.index.tz_convert(tz)
    flows = iv_data["flow_cfs"].where(iv_data["flow_cfs"] >= 0)

    step_minutes = pd.Series(local_index).diff().dropna().dt.total_seconds().median() / 60.0
    if not np.isfinite(step_minutes) or step_minutes <= 0:
        raise ValueError(
            "Could not infer a sampling interval from iv_data's index (timestamps must "
            "be strictly increasing)."
        )
    expected_obs = 1440.0 / step_minutes

    df = pd.DataFrame({"date": local_index.date, "flow_cfs": flows.to_numpy()})
    grouped = df.groupby("date")["flow_cfs"]

    result = grouped.agg(min_flow_cfs="min", max_flow_cfs="max", mean_flow_cfs="mean")
    result["std_flow_cfs"] = grouped.std(ddof=1)
    result["range_cfs"] = result["max_flow_cfs"] - result["min_flow_cfs"]
    result["cv"] = np.where(
        result["mean_flow_cfs"] > 0, result["std_flow_cfs"] / result["mean_flow_cfs"], np.nan
    )
    result["n_obs"] = grouped.apply(lambda s: int(s.notna().sum()))
    result["expected_obs"] = expected_obs
    result["complete"] = result["n_obs"] >= (expected_obs * min_completeness_frac)

    return result.reset_index().sort_values("date").reset_index(drop=True)


def diel_variation_summary(daily_diel: pd.DataFrame) -> pd.Series:
    """Period-of-record summary of a :func:`diel_variation` table.

    Computed only over rows marked ``complete``, so a period dominated by
    gappy days doesn't quietly average in unreliable single-observation
    ranges. Pass any slice of a `diel_variation` table -- filtered by
    month, by season, by year, or by any other criterion -- to get a
    summary over just that period; this function does not do any grouping
    of its own, by design, since "over a period" can mean whatever period
    the caller is working with.

    Parameters
    ----------
    daily_diel : pd.DataFrame
        A table from :func:`diel_variation`, or any row subset of one.

    Returns
    -------
    pd.Series
        Index: ``n_days`` (complete days used), ``mean_diel_amplitude_cfs``
        (mean of ``range_cfs``), ``mean_diel_cv`` (mean of ``cv``).

    Examples
    --------
    >>> daily = diel_variation(iv_data, tz="America/Los_Angeles")  # doctest: +SKIP
    >>> diel_variation_summary(daily)  # period of record  # doctest: +SKIP
    >>> july = daily[pd.DatetimeIndex(daily["date"]).month == 7]  # doctest: +SKIP
    >>> diel_variation_summary(july)  # July only  # doctest: +SKIP
    """
    complete = daily_diel[daily_diel["complete"]]
    return pd.Series(
        {
            "n_days": len(complete),
            "mean_diel_amplitude_cfs": complete["range_cfs"].mean() if len(complete) else np.nan,
            "mean_diel_cv": complete["cv"].mean() if len(complete) else np.nan,
        }
    )
