"""
flowfreq.subdaily - Sub-daily (instantaneous-series) flow metrics

Metrics computed from an instantaneous (unit-value) series rather than a
daily mean one: diel variation, the timing of each day's extremes, and
ramping rates. Grouped together because they all share the same two
problems, which no daily metric has -- every one of them must be grouped on
*local* calendar days rather than on the UTC axis the data is stored on, and
every one of them has to decide what an irregular or missing timestamp
means.

See ``docs/SUBDAILY_METRICS_DESIGN.md`` for the design, including the
reasoning behind each convention below.

Method notes
------------
**Descriptive, not diagnostic.** Nothing in this module classifies a regime
as natural or regulated, and none of it should be read that way. An
afternoon snowmelt peak and an afternoon load-following peak produce the
same hour-of-maximum; a diversion shutting off and a thunderstorm recession
produce the same down-ramp. Telling them apart needs knowledge of upstream
operations that this library has no way to obtain and does not guess at.
These numbers are inputs to that judgement.

**Local days, always.** Every per-day function here takes an explicit `tz`
with no default. Grouping a Pacific gage's record by its UTC index splits
each local day across two UTC-labeled buckets seven or eight hours out of
phase, which corrupts the range, the timing and the completeness of both.
There is no sensible default because the correct zone is a property of the
gage, not of the data.

**Timing is circular.** The mean of an hour-of-day column is *not* its
arithmetic mean: hours 23.5 and 0.5 average arithmetically to 12.0, the one
time of day such a site never peaks, and the error is worst for the sites
whose timing is most consistent. Use :func:`extreme_timing_summary` or
:func:`circular_hour_statistics`, which report a mean hour alongside the
concentration `r` that says whether that mean means anything.

**Time differences are measured on UTC, values are labelled in local
time.** A local wall-clock difference across a daylight-saving fall-back
transition is negative for a real forward hour, which turns a modest ramp
into a large one of the wrong sign. Rates therefore use the UTC axis for
``dt`` while days and hours are read from local time.

**Negative values are treated as missing**, matching every daily function in
:mod:`flowfreq.regime` -- a data artifact, not a legitimate reading.
"""

from __future__ import annotations

from typing import Any, Dict, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

#: A one-dimensional column of numbers as this module accepts it. Spelled out
#: rather than left as ``Sequence[float]`` because these functions hand each
#: other numpy arrays and pandas Series, neither of which is a Sequence --
#: and rather than as ``npt.ArrayLike``, which would also admit a bare scalar
#: that the 1-D indexing here cannot handle.
NumericColumn = Union[Sequence[float], np.ndarray, "pd.Series[Any]"]

#: The same, for a column of dates.
DateColumn = Union[Sequence[object], np.ndarray, "pd.Series[Any]"]

#: Value columns whose units are known, mapped to the unit string used in
#: ramping-rate column names. Extend this rather than passing ``units=`` at
#: every call site if a new series type becomes common.
KNOWN_VALUE_UNITS: Dict[str, str] = {
    "flow_cfs": "cfs",
    "gage_height_ft": "ft",
}

#: Units that measure a true-zero physical quantity, so that a percent change
#: in them is meaningful. Gage height is deliberately absent: it is measured
#: from an arbitrary local datum, so a "5% rise in stage" has no referent --
#: see :func:`ramping_rates` and the design doc S4.3.
TRUE_ZERO_UNITS: Tuple[str, ...] = ("cfs",)

#: How many multiples of the record's own median time step an interval may
#: span before it is treated as a gap rather than a ramp.
DEFAULT_MAX_GAP_MULTIPLE: float = 3.0

#: Which timestamp of an interval decides the local day it is counted in.
INTERVAL_LABELS: Tuple[str, ...] = ("end", "start")


def _local_index(iv_data: pd.DataFrame, tz: str) -> pd.DatetimeIndex:
    """Convert an instantaneous frame's index to `tz`, checking it is usable.

    Raises
    ------
    TypeError
        The index is not tz-aware, so it cannot be converted to a named zone
        and the local calendar day it belongs to is unknowable.
    """
    index = pd.DatetimeIndex(iv_data.index)
    if index.tz is None:
        raise TypeError(
            "iv_data must have a timezone-aware index to be grouped by local calendar "
            "day; got a naive index. download_instantaneous_flow returns a UTC index."
        )
    return index.tz_convert(tz)


def _median_step_minutes(local_index: pd.DatetimeIndex) -> float:
    """Median spacing of a timestamp index, in minutes.

    Raises
    ------
    ValueError
        The spacing is not a positive finite number, which means the
        timestamps are not strictly increasing.
    """
    step = pd.Series(local_index).diff().dropna().dt.total_seconds().median() / 60.0
    if not np.isfinite(step) or step <= 0:
        raise ValueError(
            "Could not infer a sampling interval from iv_data's index (timestamps must "
            "be strictly increasing)."
        )
    return float(step)


def _valid_values(iv_data: pd.DataFrame, value_col: str) -> pd.Series:
    """The value column with negatives masked to NaN as missing.

    Raises
    ------
    KeyError
        `value_col` is not in the frame. Named explicitly because the default
        differs between a discharge and a gage-height series.
    """
    if value_col not in iv_data.columns:
        raise KeyError(
            f"iv_data has no {value_col!r} column; columns are {list(iv_data.columns)}. "
            f"Pass value_col= to name the series to use."
        )
    values = iv_data[value_col]
    return values.where(values >= 0)


def _local_day_hours(dates: DateColumn, tz: str) -> np.ndarray:
    """Length in hours of each given local calendar day.

    Almost always 24, but 23 or 25 on a daylight-saving transition day. Used
    for `expected_obs` so that the two transition days of each year are not
    mismarked as incomplete (spring) or over-complete (autumn).
    """
    naive = pd.to_datetime(pd.Series(list(dates)))
    starts = naive.dt.tz_localize(tz, ambiguous=True, nonexistent="shift_forward")
    ends = (naive + pd.Timedelta(days=1)).dt.tz_localize(
        tz, ambiguous=True, nonexistent="shift_forward"
    )
    # Subtracted as Series, not via .to_numpy(): a tz-aware .to_numpy() is an
    # object array of Timestamps, which made expected_obs an object column.
    return np.asarray((ends - starts).dt.total_seconds() / 3600.0, dtype=float)


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
        :meth:`flowfreq.usgs.USGSgage.download_instantaneous_flow`.
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
          flowfreq.lowflow's module docstring on the same issue for annual
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

    **Expected observations per day** are the day's actual local length
    divided by the *median* time step across the entire input, not a
    hardcoded assumption like 15 minutes -- NWIS's most common
    instantaneous interval, but not the only one a logger might report at.
    The length is 24 hours on almost every day, but 23 on a daylight-saving
    spring-forward day and 25 on a fall-back day, so a full record of either
    reads as exactly complete. (Before v0.10.0 this used a fixed 24 hours,
    which marked a full spring-forward day incomplete at any
    `min_completeness_frac` above 23/24 and a full fall-back day as 25/24
    complete.) If the sampling interval genuinely
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

    local_index = _local_index(iv_data, tz)
    flows = _valid_values(iv_data, "flow_cfs")

    step_minutes = _median_step_minutes(local_index)

    df = pd.DataFrame({"date": local_index.date, "flow_cfs": flows.to_numpy()})
    grouped = df.groupby("date")["flow_cfs"]

    result = grouped.agg(min_flow_cfs="min", max_flow_cfs="max", mean_flow_cfs="mean")
    result["std_flow_cfs"] = grouped.std(ddof=1)
    result["range_cfs"] = result["max_flow_cfs"] - result["min_flow_cfs"]
    result["cv"] = np.where(
        result["mean_flow_cfs"] > 0, result["std_flow_cfs"] / result["mean_flow_cfs"], np.nan
    )
    result["n_obs"] = grouped.apply(lambda s: int(s.notna().sum()))
    # Each local day's actual length, as in daily_extreme_timing and
    # ramping_rates: 23 hours on a spring-forward day, 25 on a fall-back day.
    daily = pd.DataFrame(result).reset_index()
    daily["expected_obs"] = _local_day_hours(daily["date"], tz) * 60.0 / step_minutes
    daily["complete"] = daily["n_obs"] >= (daily["expected_obs"] * min_completeness_frac)

    return daily.sort_values("date").reset_index(drop=True)


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


def circular_hour_statistics(hours: NumericColumn) -> pd.Series:
    """Circular (directional) statistics for a collection of hour-of-day values.

    An hour-of-day is an angle on a 24-hour circle, not a point on a line,
    and must be summarized as one. The arithmetic mean of hours 23.5 and 0.5
    is 12.0 -- noon, the one time of day such a site never peaks. The error
    is largest precisely for the sites whose timing is most consistent,
    which is what makes it dangerous: it does not look like an error in the
    output, it looks like a site with a midday peak.

    Parameters
    ----------
    hours : sequence, ndarray or Series of float
        Hour-of-day values in ``[0, 24)``. NaNs are dropped. Values outside
        the range are wrapped modulo 24 rather than rejected, so an hour
        expressed as 25.0 is read as 1.0.

    Returns
    -------
    pd.Series
        - ``n`` : number of non-NaN values used
        - ``mean_hour`` : circular mean, in ``[0, 24)``. NaN if `n` is 0, or
          if the resultant is exactly zero, where no mean direction exists at
          all. Note that perfectly opposed values usually do *not* land
          exactly on zero in floating point -- hours 0 and 12 give a
          resultant around 1e-17 and so a finite but entirely arbitrary mean
          of 6.0. That is not treated as degenerate, because there is no
          principled cutoff separating it from a merely uninformative sample,
          and NaN-ing real data at an invented threshold would be worse. It
          is instead what ``concentration`` is for: read it first.
        - ``concentration`` : the resultant length `r` in ``[0, 1]``. 1 means
          every value is identical; 0 means they are spread evenly round the
          clock. **Read this before reading ``mean_hour``** -- a mean hour
          with `r` near 0 is an arbitrary direction, not a peak time.
        - ``circular_std_hours`` : ``sqrt(-2 ln r)`` in hours. Infinite as
          `r` approaches 0, correctly: scattered timing has no spread in
          hours to quote.
        - ``rayleigh_p`` : p-value of the Rayleigh test of uniformity, using
          Zar's approximation. A large value says the data are consistent
          with no preferred hour at all. NaN if `n` is 0.

    Notes
    -----
    The Rayleigh approximation is
    ``p = exp(sqrt(1 + 4n + 4(n^2 - R^2)) - (1 + 2n))`` with ``R = n*r``
    (Zar, *Biostatistical Analysis*). It is a test against the uniform
    distribution only: rejecting it says timing is not uniform, not that it
    is unimodal. A reach with a morning and an evening peak can be strongly
    non-uniform and still have a meaningless mean direction, which is why
    `concentration` is reported rather than left to the caller to derive.

    Examples
    --------
    >>> stats = circular_hour_statistics([23.5, 0.5])
    >>> round(float(stats["mean_hour"]), 6)
    0.0
    >>> round(float(stats["concentration"]), 6)
    0.991445
    """
    values = np.asarray(hours, dtype=float)
    values = values[np.isfinite(values)]
    n = int(values.size)

    if n == 0:
        return pd.Series(
            {
                "n": 0,
                "mean_hour": np.nan,
                "concentration": np.nan,
                "circular_std_hours": np.nan,
                "rayleigh_p": np.nan,
            }
        )

    theta = 2.0 * np.pi * (values % 24.0) / 24.0
    cos_mean = float(np.mean(np.cos(theta)))
    sin_mean = float(np.mean(np.sin(theta)))
    # Clamped because hypot of two float means can land a hair above 1.0 for
    # identical inputs, which would make log(r) positive and the spread below
    # the square root of a negative number -- NaN for the most consistent
    # timing possible, the opposite of the truth.
    r = float(min(np.hypot(cos_mean, sin_mean), 1.0))

    if r == 0.0:
        # Perfectly opposed values: no mean direction exists, and -2*ln(0)
        # is not a finite spread. Both are reported as such rather than as a
        # number that happens to fall out of atan2(0, 0).
        mean_hour = np.nan
        circular_std_hours = np.inf
    else:
        # The outer modulo matters: a mean direction at the wrap point comes
        # back from atan2 as a hair under zero, and `% (2*pi)` then lifts it to
        # almost exactly 2*pi -- reported as hour 24.0, which is not in [0, 24)
        # and reads as a midnight peak having been missed by a full day.
        mean_hour = (
            float(np.arctan2(sin_mean, cos_mean) % (2.0 * np.pi)) * 24.0 / (2.0 * np.pi)
        ) % 24.0
        # Branched rather than clamped, for the same reason as the clamp
        # above plus one more: at r == 1, -2*log(r) is -0.0, and sqrt(-0.0)
        # is -0.0, so perfectly consistent timing would report a *negative*
        # standard deviation.
        spread = -2.0 * np.log(r)
        circular_std_hours = 0.0 if spread <= 0.0 else float(np.sqrt(spread)) * 24.0 / (2.0 * np.pi)

    resultant = n * r
    rayleigh_p = float(
        np.exp(np.sqrt(1.0 + 4.0 * n + 4.0 * (n**2 - resultant**2)) - (1.0 + 2.0 * n))
    )

    return pd.Series(
        {
            "n": n,
            "mean_hour": mean_hour,
            "concentration": r,
            "circular_std_hours": circular_std_hours,
            "rayleigh_p": min(rayleigh_p, 1.0),
        }
    )


def daily_extreme_timing(
    iv_data: pd.DataFrame,
    tz: str,
    *,
    value_col: str = "flow_cfs",
    min_completeness_frac: float = 0.9,
) -> pd.DataFrame:
    """Hour of the daily maximum and minimum, per local calendar day.

    *When* in the day flow peaks and bottoms out, as a fractional local hour.
    The pair is the point: a morning minimum with a late-afternoon maximum is
    the snowmelt signature, while either half alone is ambiguous. Like every
    function in this module this is descriptive -- an afternoon melt peak and
    an afternoon load-following peak give the same hour, and separating them
    needs knowledge of upstream operations this cannot infer.

    Parameters
    ----------
    iv_data : pd.DataFrame
        Instantaneous series with a tz-aware datetime index and a `value_col`
        column -- the shape returned by
        :meth:`flowfreq.usgs.USGSgage.download_instantaneous_flow`.
    tz : str
        IANA time zone to group local calendar days in, e.g.
        ``"America/Los_Angeles"``. Required, with no default: hour-of-peak in
        UTC at a Pacific gage is seven or eight hours out of phase with the
        diel cycle being measured, which makes every number meaningless
        rather than merely shifted.
    value_col : str
        Column to take extremes of. Default ``"flow_cfs"``; pass
        ``"gage_height_ft"`` for a stage series.
    min_completeness_frac : float
        Fraction of the day's expected observation count that must be present
        for a day to be marked ``complete``. Default 0.9.

    Returns
    -------
    pd.DataFrame
        One row per local calendar day with at least one valid observation.
        Columns:

        - ``date`` : the local calendar date
        - ``hour_of_max``, ``hour_of_min`` : fractional local hour in
          ``[0, 24)`` of the day's extreme -- ``14.25`` is 14:15 local. NaN
          on a flat day; see Notes.
        - ``max_value``, ``min_value`` : the extremes themselves, in
          `value_col`'s units
        - ``n_tied_max``, ``n_tied_min`` : how many of the day's observations
          equal that extreme
        - ``n_obs``, ``expected_obs``, ``complete``

    Notes
    -----
    **A flat day has no timing, and gets NaN rather than 0.0.** If every
    valid observation in a day is the same value, the "hour of the maximum"
    is just whichever observation came first. Reporting that as ``0.0`` would
    be a fabricated midnight peak, and a regulated reach held at a constant
    release produces runs of exactly these days -- each one would drag the
    circular mean in :func:`extreme_timing_summary` toward midnight. Both
    hour columns are NaN when ``max_value == min_value``.

    **Ties are reported, not hidden.** When the extreme is not unique but the
    day is not flat -- two equally high afternoon crests, or a sensor whose
    resolution reports the same value twice at the peak -- the *first*
    occurrence is returned and ``n_tied_max`` says how many there were. So a
    reader can tell an identified moment from one of four equally high
    readings. Returning the first silently, with nothing recording that a
    choice was made, is the failure this avoids.

    **Hours are local wall-clock, including across daylight saving.** On a
    fall-back day two observations genuinely share a fractional hour and the
    local day is 25 hours long; on a spring-forward day the 02:00-03:00 hour
    does not exist and the day is 23 hours. Both are properties of local
    time, not errors, and local time is what a diel question asks about.
    ``expected_obs`` accounts for it: it is computed from each day's actual
    length, not from a fixed 1440 minutes, so the two transition days of a
    year are not mismarked. :func:`diel_variation` and :func:`ramping_rates`
    use the same convention.

    **Completeness does not gate the values**, matching
    :func:`diel_variation`: a day missing some readings still usually
    captures its peak, so the timing is reported and ``complete`` says how
    much to trust it. :func:`extreme_timing_summary` pools only complete
    days.

    Raises
    ------
    ValueError
        Fewer than two timestamps, so no sampling interval can be inferred.
    KeyError
        `value_col` is not a column of `iv_data`.
    TypeError
        `iv_data`'s index is not tz-aware.

    Examples
    --------
    >>> timing = daily_extreme_timing(iv, tz="America/Los_Angeles")  # doctest: +SKIP
    >>> extreme_timing_summary(timing, which="max")  # doctest: +SKIP
    """
    if len(iv_data) < 2:
        raise ValueError(
            "daily_extreme_timing needs at least 2 timestamps to infer a sampling "
            f"interval; got {len(iv_data)}"
        )

    local_index = _local_index(iv_data, tz)
    values = _valid_values(iv_data, value_col)
    step_minutes = _median_step_minutes(local_index)

    hour = local_index.hour + local_index.minute / 60.0 + local_index.second / 3600.0
    frame = pd.DataFrame(
        {"date": local_index.date, "value": values.to_numpy(), "hour": np.asarray(hour)}
    ).dropna(subset=["value"])

    if frame.empty:
        return pd.DataFrame(
            columns=[
                "date",
                "hour_of_max",
                "hour_of_min",
                "max_value",
                "min_value",
                "n_tied_max",
                "n_tied_min",
                "n_obs",
                "expected_obs",
                "complete",
            ]
        )

    grouped = frame.groupby("date", sort=True)
    result = grouped.agg(max_value=("value", "max"), min_value=("value", "min"))
    result["n_obs"] = grouped["value"].size()

    # idxmax/idxmin return the first occurrence of a tied extreme, which is
    # the documented convention; n_tied_* below records that a tie happened.
    result["hour_of_max"] = frame.loc[grouped["value"].idxmax(), "hour"].to_numpy()
    result["hour_of_min"] = frame.loc[grouped["value"].idxmin(), "hour"].to_numpy()

    extremes = frame.join(result[["max_value", "min_value"]], on="date")
    result["n_tied_max"] = (
        extremes["value"].eq(extremes["max_value"]).groupby(extremes["date"]).sum().astype(int)
    )
    result["n_tied_min"] = (
        extremes["value"].eq(extremes["min_value"]).groupby(extremes["date"]).sum().astype(int)
    )

    flat = result["max_value"] == result["min_value"]
    result.loc[flat, ["hour_of_max", "hour_of_min"]] = np.nan

    result = result.reset_index()
    result["expected_obs"] = _local_day_hours(result["date"], tz) * 60.0 / step_minutes
    result["complete"] = result["n_obs"] >= (result["expected_obs"] * min_completeness_frac)

    return result[
        [
            "date",
            "hour_of_max",
            "hour_of_min",
            "max_value",
            "min_value",
            "n_tied_max",
            "n_tied_min",
            "n_obs",
            "expected_obs",
            "complete",
        ]
    ].reset_index(drop=True)


def extreme_timing_summary(daily_timing: pd.DataFrame, which: str = "max") -> pd.Series:
    """Circular summary of a :func:`daily_extreme_timing` table.

    Computed over rows marked ``complete`` only, matching
    :func:`diel_variation_summary`, so a period dominated by gappy days does
    not average in unreliable timing. Pass any slice of a table -- one month,
    one season, one year -- to summarize just that period; this function does
    no grouping of its own, since "over a period" means whatever period the
    caller is working with.

    Parameters
    ----------
    daily_timing : pd.DataFrame
        A table from :func:`daily_extreme_timing`, or any row subset of one.
    which : str
        ``"max"`` to summarize ``hour_of_max``, ``"min"`` for ``hour_of_min``.

    Returns
    -------
    pd.Series
        ``n_days`` (complete days with a defined timing), then the fields of
        :func:`circular_hour_statistics` prefixed with ``which``:
        ``mean_hour_of_max``, ``concentration_of_max``,
        ``circular_std_hours_of_max``, ``rayleigh_p_of_max``.

    Notes
    -----
    ``n_days`` counts complete days whose timing is *defined*, so it is
    smaller than the number of complete rows when the period contains flat
    days (whose timing is NaN -- see :func:`daily_extreme_timing`'s Notes).

    **Always read ``concentration`` next to ``mean_hour``.** A mean hour with
    a concentration near 0 is an arbitrary direction on the clock, not a peak
    time, and nothing about the number itself reveals that.

    Raises
    ------
    ValueError
        `which` is not ``"max"`` or ``"min"``.
    """
    if which not in ("max", "min"):
        raise ValueError(f"which must be 'max' or 'min', got {which!r}")

    column = f"hour_of_{which}"
    complete = daily_timing[daily_timing["complete"]]
    stats = circular_hour_statistics(complete[column].to_numpy())

    return pd.Series(
        {
            "n_days": int(stats["n"]),
            f"mean_hour_of_{which}": stats["mean_hour"],
            f"concentration_of_{which}": stats["concentration"],
            f"circular_std_hours_of_{which}": stats["circular_std_hours"],
            f"rayleigh_p_of_{which}": stats["rayleigh_p"],
        }
    )


def _resolve_units(value_col: str, units: Optional[str]) -> str:
    """Unit string for a value column's rate-column names.

    Raises
    ------
    ValueError
        The column is not in :data:`KNOWN_VALUE_UNITS` and no `units` was
        given. A rate column whose units are unknown would be named
        ``max_up_ramp_per_hr``, which reads as dimensionless and gets
        misinterpreted; naming the unit is required rather than guessed.
    """
    if units is not None:
        return units
    if value_col in KNOWN_VALUE_UNITS:
        return KNOWN_VALUE_UNITS[value_col]
    raise ValueError(
        f"Cannot infer units for value_col={value_col!r} (known: "
        f"{sorted(KNOWN_VALUE_UNITS)}). Pass units= so the rate columns are named "
        f"with their unit rather than left dimensionless."
    )


def ramping_rates(
    iv_data: pd.DataFrame,
    tz: str,
    *,
    value_col: str = "flow_cfs",
    units: Optional[str] = None,
    relative: Optional[bool] = None,
    max_gap_hours: Optional[float] = None,
    interval_label: str = "end",
    min_completeness_frac: float = 0.9,
) -> pd.DataFrame:
    """Per-day rate of change of an instantaneous series -- ramping rates.

    How fast flow (or stage) changes between consecutive observations,
    summarized per local calendar day as the largest rise, the largest fall,
    and the mean absolute rate. Down-ramping is the metric salmon-habitat
    work consumes: redd dewatering and fry stranding are driven by how fast
    the water's edge retreats, not by the day's total range, which
    :func:`diel_variation` already covers.

    Parameters
    ----------
    iv_data : pd.DataFrame
        Instantaneous series with a tz-aware datetime index and a `value_col`
        column, as returned by
        :meth:`flowfreq.usgs.USGSgage.download_instantaneous_flow` or
        :meth:`flowfreq.usgs.USGSgage.download_instantaneous_stage`.
    tz : str
        IANA time zone the local calendar days are defined in. Required, with
        no default, for the reason given in this module's docstring. Note
        that this labels the *days*; the time differences themselves are
        measured on the UTC axis -- see Notes.
    value_col : str
        Column to differentiate. Default ``"flow_cfs"``; pass
        ``"gage_height_ft"`` for stage.
    units : str, optional
        Unit string used in the output column names, e.g. ``max_up_ramp_``
        + ``cfs`` + ``_per_hr``. Inferred from `value_col` via
        :data:`KNOWN_VALUE_UNITS` when omitted.
    relative : bool, optional
        Also report percent-per-hour rates. Default ``None``, meaning True
        for a true-zero quantity (discharge) and False otherwise. Passing
        True for a series whose units are not in :data:`TRUE_ZERO_UNITS`
        raises -- see Notes on why a percentage of stage is not a number.
    max_gap_hours : float, optional
        Intervals longer than this are treated as gaps rather than ramps and
        excluded. Default ``None``, meaning
        :data:`DEFAULT_MAX_GAP_MULTIPLE` times the record's own median time
        step, so the default is right for a 15-minute record and for an
        hourly one without the caller having to know which they have.
    interval_label : str
        Which end of an interval decides the local day it counts in:
        ``"end"`` (default) or ``"start"``. See Notes.
    min_completeness_frac : float
        Fraction of the day's expected observation count required for
        ``complete``. Default 0.9.

    Returns
    -------
    pd.DataFrame
        One row per local calendar day that has at least one interval
        assigned to it. With ``u`` standing for the resolved unit:

        - ``date`` : the local calendar date
        - ``max_up_ramp_<u>_per_hr`` : largest positive rate. ``0.0`` when
          intervals were measured but none rose; NaN when no valid interval
          was available at all -- see Notes.
        - ``max_down_ramp_<u>_per_hr`` : the **signed** minimum, so always
          ``<= 0``
        - ``mean_abs_ramp_<u>_per_hr`` : mean of the absolute rates
        - ``max_up_ramp_pct_per_hr``, ``max_down_ramp_pct_per_hr``,
          ``mean_abs_ramp_pct_per_hr`` : the same three in percent per hour,
          present only when `relative` resolves True
        - ``n_reversals`` : sign changes in the day's rate sequence
        - ``n_intervals`` : valid intervals used
        - ``n_intervals_gapped`` : intervals dropped for exceeding
          `max_gap_hours`
        - ``n_obs``, ``expected_obs``, ``complete``

    Notes
    -----
    **Time differences are measured on the UTC axis, never on local
    wall-clock time.** This is not stylistic. Across a daylight-saving
    fall-back transition the local clock runs 01:59 -> 01:00, so a local-time
    difference is *negative* across a real forward hour: a 15-minute interval
    spanning the transition would compute as ``-45 minutes``, turning a
    modest real ramp into a large one of the wrong sign.

    **Down-ramp is signed and always ``<= 0``.** "A maximum down-ramp of 50"
    is ambiguous about sign in exactly the context where sign decides whether
    a limit was exceeded. A value that cannot be positive cannot be misread.

    **``0.0`` and NaN mean different things for ``max_up_ramp``.** ``0.0``
    means intervals were measured and the largest rise among them was zero (a
    monotonically falling day, or a flat one); NaN means no valid interval was
    available at all (every one dropped as a gap). Collapsing them would make
    a fully-gapped day indistinguishable from a recession.

    **A gap is not a ramp.** A change computed across a six-hour hole in a
    15-minute record is the net change over six hours: it reads as a gentle
    ramp however violently flow actually moved inside the hole, and can
    equally conceal a large real ramp. Averaged in, these bias the mean
    absolute rate *downward*. Intervals over `max_gap_hours` are therefore
    excluded and counted in ``n_intervals_gapped``, so the exclusion is
    visible rather than silent.

    **Which day an interval belongs to is a real choice.** An interval spans
    two instants and can straddle local midnight. By default it counts in the
    day of its *ending* timestamp, so an overnight recession lands on the day
    whose low flow it produced -- the day a stranding question is about.
    ``interval_label="start"`` gives the other convention. The parameter is
    explicit because the choice shifts every midnight-straddling interval by
    a day, and left implicit it would make two runs disagree with no visible
    cause.

    **Percent-per-hour on a stage series is meaningless, and refused.** Gage
    height is measured from an arbitrary local datum: a gage reading 4.0 ft
    rising to 4.2 ft has not "risen 5%" in any physical sense, since the same
    water-surface change at a gage whose datum sits two feet lower reads as
    2.9%, and neither figure divides by the actual depth of water. The result
    would look exactly like a regulatory ramping rate while having no
    referent, so `relative=True` on a non-true-zero series raises instead.

    **Percent rates are NaN where the preceding value is not positive.** A
    percentage change from zero is not a number, and an intermittent reach
    produces these. The absolute rate for the same interval is still
    reported.

    Raises
    ------
    ValueError
        Fewer than two timestamps; `interval_label` not in
        :data:`INTERVAL_LABELS`; `max_gap_hours` not positive; `relative` is
        True for a series whose units are not in :data:`TRUE_ZERO_UNITS`; or
        units could not be inferred for `value_col`.
    KeyError
        `value_col` is not a column of `iv_data`.
    TypeError
        `iv_data`'s index is not tz-aware.

    Examples
    --------
    >>> ramps = ramping_rates(iv, tz="America/Los_Angeles")  # doctest: +SKIP
    >>> stage_ramps = ramping_rates(  # doctest: +SKIP
    ...     stage_iv, tz="America/Los_Angeles", value_col="gage_height_ft"
    ... )
    """
    if len(iv_data) < 2:
        raise ValueError(
            f"ramping_rates needs at least 2 timestamps to form an interval; " f"got {len(iv_data)}"
        )
    if interval_label not in INTERVAL_LABELS:
        raise ValueError(
            f"interval_label must be one of {list(INTERVAL_LABELS)}, got {interval_label!r}"
        )

    unit = _resolve_units(value_col, units)
    if relative is None:
        relative = unit in TRUE_ZERO_UNITS
    elif relative and unit not in TRUE_ZERO_UNITS:
        raise ValueError(
            f"relative=True was requested for a series in {unit!r}, which is not a "
            f"true-zero quantity (known true-zero units: {list(TRUE_ZERO_UNITS)}). A "
            f"percent-per-hour rate on gage height divides by a stage measured from an "
            f"arbitrary local datum, so it has no physical referent while looking "
            f"exactly like a regulatory ramping rate. Use the absolute ft/hr rate."
        )

    local_index = _local_index(iv_data, tz)
    values = _valid_values(iv_data, value_col)
    step_minutes = _median_step_minutes(local_index)

    if max_gap_hours is None:
        max_gap_hours = DEFAULT_MAX_GAP_MULTIPLE * step_minutes / 60.0
    elif max_gap_hours <= 0:
        raise ValueError(f"max_gap_hours must be positive, got {max_gap_hours}")

    intervals = _build_intervals(
        local_index, values, tz, max_gap_hours=max_gap_hours, interval_label=interval_label
    )
    daily = _aggregate_intervals(intervals, unit=unit, relative=relative)

    # Observation counts come from the observations themselves, not from the
    # interval table: a day's coverage is a property of its readings, and
    # interval_label would otherwise shift it.
    obs = pd.DataFrame({"date": local_index.date, "value": values.to_numpy()}).dropna(
        subset=["value"]
    )
    counts = obs.groupby("date", sort=True).size().rename("n_obs")
    daily = daily.join(counts, how="left")
    daily["n_obs"] = daily["n_obs"].fillna(0).astype(int)

    daily = daily.reset_index()
    daily["expected_obs"] = _local_day_hours(daily["date"], tz) * 60.0 / step_minutes
    daily["complete"] = daily["n_obs"] >= (daily["expected_obs"] * min_completeness_frac)

    return daily.reset_index(drop=True)


def _build_intervals(
    local_index: pd.DatetimeIndex,
    values: pd.Series,
    tz: str,
    *,
    max_gap_hours: float,
    interval_label: str,
) -> pd.DataFrame:
    """One row per consecutive pair of valid observations, with its rate.

    Pairs are formed between *valid* observations, so a masked negative or a
    NaN does not create an interval of its own -- the pair spans across it and
    is then judged against `max_gap_hours` like any other, which is the right
    treatment: a reading dropped as an artifact leaves a hole no different
    from a logger outage.
    """
    valid = pd.DataFrame(
        {"utc": local_index.tz_convert("UTC"), "local": local_index, "value": values.to_numpy()}
    ).dropna(subset=["value"])

    if len(valid) < 2:
        return pd.DataFrame(columns=["date", "dt_hours", "delta", "prev_value", "rate", "gapped"])

    dt_hours = valid["utc"].diff().dt.total_seconds().to_numpy()[1:] / 3600.0
    delta = valid["value"].diff().to_numpy()[1:]
    prev_value = valid["value"].to_numpy()[:-1]
    label_local = (
        valid["local"].to_numpy()[1:]
        if interval_label == "end"
        else (valid["local"].to_numpy()[:-1])
    )

    gapped = dt_hours > max_gap_hours
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.where(gapped, np.nan, delta / dt_hours)

    return pd.DataFrame(
        {
            "date": pd.DatetimeIndex(label_local).date,
            "dt_hours": dt_hours,
            "delta": delta,
            "prev_value": prev_value,
            "rate": rate,
            "gapped": gapped,
        }
    )


def _count_reversals(rates: pd.Series) -> int:
    """Sign changes in a day's sequence of rates.

    The standard hydropeaking "number of reversals" count. Zero-rate
    intervals are not sign changes and are skipped rather than counted as a
    crossing, so a flat spell between a rise and a fall contributes one
    reversal, not two.
    """
    signs = np.sign(rates.dropna().to_numpy())
    signs = signs[signs != 0]
    if signs.size < 2:
        return 0
    return int(np.count_nonzero(np.diff(signs) != 0))


def _aggregate_intervals(intervals: pd.DataFrame, *, unit: str, relative: bool) -> pd.DataFrame:
    """Collapse an interval table to one row per local day."""
    columns = [
        f"max_up_ramp_{unit}_per_hr",
        f"max_down_ramp_{unit}_per_hr",
        f"mean_abs_ramp_{unit}_per_hr",
    ]
    if relative:
        columns += [
            "max_up_ramp_pct_per_hr",
            "max_down_ramp_pct_per_hr",
            "mean_abs_ramp_pct_per_hr",
        ]
    columns += ["n_reversals", "n_intervals", "n_intervals_gapped"]

    if intervals.empty:
        empty = pd.DataFrame(columns=columns)
        empty.index.name = "date"
        return empty

    work = intervals.copy()
    if relative:
        with np.errstate(divide="ignore", invalid="ignore"):
            pct = np.where(
                work["prev_value"].to_numpy() > 0,
                100.0 * work["rate"].to_numpy() / work["prev_value"].to_numpy(),
                np.nan,
            )
        work["pct_rate"] = pct

    grouped = work.groupby("date", sort=True)
    out = pd.DataFrame(index=grouped.size().index)
    out.index.name = "date"

    # clip(upper=0)/clip(lower=0) rather than filtering to positive or
    # negative intervals: on a day whose every interval falls, the largest
    # rise really is 0.0 ("no rise was observed"), which is different from the
    # NaN a day with no valid interval at all gets. See the docstring.
    out[f"max_up_ramp_{unit}_per_hr"] = grouped["rate"].apply(lambda s: s.clip(lower=0).max())
    out[f"max_down_ramp_{unit}_per_hr"] = grouped["rate"].apply(lambda s: s.clip(upper=0).min())
    out[f"mean_abs_ramp_{unit}_per_hr"] = grouped["rate"].apply(lambda s: s.abs().mean())

    if relative:
        out["max_up_ramp_pct_per_hr"] = grouped["pct_rate"].apply(lambda s: s.clip(lower=0).max())
        out["max_down_ramp_pct_per_hr"] = grouped["pct_rate"].apply(lambda s: s.clip(upper=0).min())
        out["mean_abs_ramp_pct_per_hr"] = grouped["pct_rate"].apply(lambda s: s.abs().mean())

    out["n_reversals"] = grouped["rate"].apply(_count_reversals).astype(int)
    out["n_intervals"] = grouped["rate"].apply(lambda s: int(s.notna().sum())).astype(int)
    out["n_intervals_gapped"] = grouped["gapped"].sum().astype(int)

    return out[columns]


def ramping_rate_summary(
    daily_ramps: pd.DataFrame, limits: Optional[Mapping[str, float]] = None
) -> pd.Series:
    """Period summary of a :func:`ramping_rates` table, with limit exceedance.

    Computed over rows marked ``complete`` only, matching the other summary
    functions in this module. Does no grouping of its own: pass any row
    subset -- one month, one season, the whole record -- to summarize just
    that period.

    Parameters
    ----------
    daily_ramps : pd.DataFrame
        A table from :func:`ramping_rates`, or any row subset of one.
    limits : mapping of str to float, optional
        Ramping limits to count exceedances against, keyed by column name,
        e.g. ``{"max_down_ramp_cfs_per_hr": -100.0}``. For each entry,
        ``n_days_exceeding_<column>`` and ``frac_days_exceeding_<column>`` are
        added. Comparison direction follows the **sign of the limit**: a
        negative limit counts days at or below it (a down-ramp limit), a
        positive or zero limit counts days at or above it (an up-ramp limit).

    Returns
    -------
    pd.Series
        ``n_days``, then for every rate column present in `daily_ramps` its
        ``mean_`` and ``extreme_`` values, plus ``mean_n_reversals`` and any
        exceedance fields from `limits`. ``extreme_`` is the maximum for an
        up-ramp column and the minimum for a down-ramp one, so in both cases
        it is the most severe value in the period rather than the numerically
        largest.

    Raises
    ------
    KeyError
        A key of `limits` is not a column of `daily_ramps`. Silently skipping
        it would report "no exceedances" for a limit that was never checked.

    Examples
    --------
    >>> summary = ramping_rate_summary(  # doctest: +SKIP
    ...     ramps, limits={"max_down_ramp_cfs_per_hr": -100.0}
    ... )
    """
    complete = daily_ramps[daily_ramps["complete"]] if len(daily_ramps) else daily_ramps
    result: Dict[str, float] = {"n_days": float(len(complete))}

    rate_columns = [c for c in daily_ramps.columns if c.endswith("_per_hr")]
    for column in rate_columns:
        series = complete[column] if len(complete) else pd.Series(dtype=float)
        result[f"mean_{column}"] = float(series.mean()) if len(series) else np.nan
        if len(series) == 0 or series.isna().all():
            result[f"extreme_{column}"] = np.nan
        elif column.startswith("max_down_ramp"):
            result[f"extreme_{column}"] = float(series.min())
        else:
            result[f"extreme_{column}"] = float(series.max())

    if "n_reversals" in daily_ramps.columns:
        result["mean_n_reversals"] = (
            float(complete["n_reversals"].mean()) if len(complete) else np.nan
        )

    for column, limit in (limits or {}).items():
        if column not in daily_ramps.columns:
            raise KeyError(
                f"limits names column {column!r}, which is not in this table "
                f"(columns: {list(daily_ramps.columns)}). Reporting zero exceedances "
                f"for a limit that was never checked would be worse than raising."
            )
        series = complete[column] if len(complete) else pd.Series(dtype=float)
        exceeding = (series <= limit) if limit < 0 else (series >= limit)
        n_exceeding = int(exceeding.sum())
        result[f"n_days_exceeding_{column}"] = float(n_exceeding)
        result[f"frac_days_exceeding_{column}"] = (
            n_exceeding / len(series) if len(series) else np.nan
        )

    return pd.Series(result)
