"""
Streamflow record extension by Maintenance Of Variance Extension, type 1 (MOVE.1).

A short-record gage ``y`` and a concurrent long-record index gage ``x`` are
related, in log space, by the line that reproduces the short gage's mean and
standard deviation rather than minimizing squared error:

    log y = mean(log y) + sign(r) * (sd(log y) / sd(log x)) * (log x - mean(log x))

with the moments taken over the concurrent period. Applied to the index gage's
full record, it estimates the short gage's flows on every day the index gage
has one.

Why not ordinary least squares: an OLS estimate shrinks toward the mean by the
factor ``r``, so an extended record built from it has too little variance --
its high flows are too low and its low flows too high, and a duration curve or
frequency curve computed from it is too flat. MOVE.1 keeps the variance and
is the extension Hirsch (1982) recommends when the extended record will be
used for its distribution rather than for individual values. Bulletin 17C
Appendix 8 uses the same idea (MOVE.3) for annual peaks.

References
----------
Hirsch, R.M., 1982, A comparison of four streamflow record extension
techniques: Water Resources Research, v. 18, no. 4, p. 1081-1088.

Helsel, D.R., and Hirsch, R.M., 2002, Statistical methods in water
resources: U.S. Geological Survey Techniques of Water-Resources
Investigations, book 4, chap. A3 (the line of organic correlation, of which
MOVE.1 is the record-extension use).

What this module cannot do
--------------------------
- **It cannot say the extension is right.** The concurrent period's moments
  are taken as the short gage's moments; a concurrent period that is
  unusually dry or wet at the index gage carries its relation, not its
  climate, into the extended record, but a relation that changes with
  climate (a losing reach, a glacier) is not seen. Validate on held-out
  concurrent data.
- **Serial correlation is ignored.** Daily flows are far from independent, so
  ``n_concurrent`` days of overlap carry much less information than the count
  suggests, and no variance or effective record length is reported. For
  annual peaks, Bulletin 17C Appendix 8's MOVE.3 and effective-record
  accounting is the method; this module does not implement it.
- **Zero and negative flows have no logarithm.** They are excluded from the
  fit and counted, and an index-gage zero produces no estimate (``NaN``),
  not a zero. An intermittent index gage is a poor choice for that reason.
- **No lag.** The relation is same-day. If the two gages' timing differs,
  shift the index series first (``flowfreq.qppq.apply_lag``).
- **One relation for the year.** Fit by season and extend each season
  separately if the relation differs between melt and baseflow.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

#: Fewest concurrent positive pairs a fit accepts. A floor against nonsense,
#: not a recommendation: Bulletin 17C Appendix 8 asks for at least ten
#: concurrent *years* for peaks, and daily data need far more than ten days.
MIN_CONCURRENT = 10


@dataclass(frozen=True)
class Move1:
    """A fitted MOVE.1 relation, in base-10 logarithms.

    Attributes
    ----------
    mean_short, sd_short : float
        Mean and sample standard deviation of ``log10`` flow at the short
        gage over the concurrent period.
    mean_index, sd_index : float
        The same at the index gage.
    correlation : float
        Pearson correlation of the concurrent ``log10`` flows.
    n_concurrent : int
        Concurrent pairs used, both flows positive.
    n_excluded_nonpositive : int
        Concurrent pairs dropped because either flow was zero or negative.
    concurrent_start, concurrent_end : pd.Timestamp
        First and last concurrent date used.
    """

    mean_short: float
    sd_short: float
    mean_index: float
    sd_index: float
    correlation: float
    n_concurrent: int
    n_excluded_nonpositive: int
    concurrent_start: pd.Timestamp
    concurrent_end: pd.Timestamp

    @property
    def slope(self) -> float:
        """``sign(r) * sd_short / sd_index``, the log-space slope."""
        return float(np.sign(self.correlation)) * self.sd_short / self.sd_index

    def predict(self, index_flows: np.ndarray | pd.Series | float) -> np.ndarray:
        """Short-gage flows estimated from index-gage flows.

        Parameters
        ----------
        index_flows : array-like
            Index-gage flows, in the units the fit used.

        Returns
        -------
        np.ndarray
            Estimated flows. ``NaN`` where the index flow is zero, negative or
            missing: a logarithmic relation says nothing there.
        """
        x = np.asarray(index_flows, dtype=float)
        out = np.full(x.shape, np.nan)
        positive = np.isfinite(x) & (x > 0)
        out[positive] = 10.0 ** (
            self.mean_short + self.slope * (np.log10(x[positive]) - self.mean_index)
        )
        return out


def _as_series(flows: pd.Series, name: str) -> pd.Series:
    if not isinstance(flows, pd.Series):
        raise TypeError(f"{name} must be a pandas Series indexed by date")
    if not isinstance(flows.index, pd.DatetimeIndex):
        raise TypeError(f"{name} must have a DatetimeIndex")
    if flows.index.has_duplicates:
        raise ValueError(f"{name} has duplicate dates")
    return flows.astype(float)


def fit_move1(short: pd.Series, index: pd.Series, *, min_concurrent: int = MIN_CONCURRENT) -> Move1:
    """Fit MOVE.1 to the concurrent record of a short gage and an index gage.

    Parameters
    ----------
    short, index : pd.Series
        Flows indexed by date (daily values, or annual peaks keyed by water
        year's date -- the function only pairs equal index labels).
    min_concurrent : int, optional
        Fewest concurrent pairs, both positive, to accept.

    Returns
    -------
    Move1

    Raises
    ------
    TypeError
        If either input is not a date-indexed Series.
    ValueError
        With duplicate dates, fewer than ``min_concurrent`` usable pairs, no
        variation at either gage, or a correlation that is not positive. A
        negative correlation between two streamflow records is a data problem,
        not a relation to extend with.
    """
    y = _as_series(short, "short")
    x = _as_series(index, "index")
    pairs = pd.concat([y.rename("y"), x.rename("x")], axis=1, join="inner", sort=True).dropna()
    usable = (pairs["y"] > 0) & (pairs["x"] > 0)
    excluded = int((~usable).sum())
    pairs = pairs[usable]
    if len(pairs) < min_concurrent:
        raise ValueError(
            f"MOVE.1 needs at least {min_concurrent} concurrent positive pairs, got "
            f"{len(pairs)} ({excluded} excluded as zero or negative)"
        )
    ly = np.log10(pairs["y"].to_numpy())
    lx = np.log10(pairs["x"].to_numpy())
    sd_y = float(np.std(ly, ddof=1))
    sd_x = float(np.std(lx, ddof=1))
    if not (sd_y > 0 and sd_x > 0):
        raise ValueError("MOVE.1 needs variation at both gages over the concurrent period")
    r = float(np.corrcoef(ly, lx)[0, 1])
    if not r > 0:
        raise ValueError(f"concurrent log flows are not positively correlated (r = {r:.3f})")
    if excluded:
        logger.info("MOVE.1: %d concurrent pair(s) with a zero or negative flow excluded", excluded)
    return Move1(
        mean_short=float(ly.mean()),
        sd_short=sd_y,
        mean_index=float(lx.mean()),
        sd_index=sd_x,
        correlation=r,
        n_concurrent=int(len(pairs)),
        n_excluded_nonpositive=excluded,
        concurrent_start=pd.Timestamp(pairs.index.min()),
        concurrent_end=pd.Timestamp(pairs.index.max()),
    )


def extend_record(short: pd.Series, index: pd.Series, fit: Move1 | None = None) -> pd.DataFrame:
    """The short gage's record, extended over the index gage's by MOVE.1.

    Observed values are kept where they exist; every other date the index gage
    covers gets the MOVE.1 estimate.

    Parameters
    ----------
    short, index : pd.Series
        Flows indexed by date.
    fit : Move1, optional
        A relation already fitted, for example on one season or a subset of
        years. Fitted from the full concurrent record when omitted.

    Returns
    -------
    pd.DataFrame
        Indexed by the union of both records' dates, sorted, with ``flow``
        and ``source`` (``"observed"``, ``"move1"``, or ``"none"`` where the
        index flow was zero, negative or missing and so no estimate exists).

    Raises
    ------
    TypeError, ValueError
        As :func:`fit_move1`.
    """
    y = _as_series(short, "short").dropna()
    x = _as_series(index, "index")
    if fit is None:
        fit = fit_move1(y, x)
    dates = y.index.union(x.index).sort_values()
    estimated = pd.Series(fit.predict(x.reindex(dates).to_numpy()), index=dates)
    observed = y.reindex(dates)
    flow = observed.where(observed.notna(), estimated)
    source = np.where(observed.notna(), "observed", np.where(estimated.notna(), "move1", "none"))
    return pd.DataFrame({"flow": flow, "source": source}, index=dates)
