"""
hydrolib.lowflow - Low-flow frequency analysis

Fits a distribution to the annual minimum n-day mean flow and computes
quantiles at specified non-exceedance probabilities -- the low-flow
counterpart to hydrolib.bulletin17c's high-flow analysis. The standard
products (7Q10, 7Q2) are the n=7, non-exceedance-probability=0.10 and 0.50
cases of the same general computation; n and the probabilities are always
arguments, never hardcoded.

Method notes
------------
**Distribution.** Log-Pearson III fit by station skew, reusing the same
log-space machinery as the high-flow side. Quantiles are computed with
:func:`hydrolib.core.lp3_frequency_factor_peakfq`, the exact gamma-quantile
method (matching peakfqr/PeakFQ's own ``qP3sub``), rather than the
Wilson-Hilferty approximation :func:`hydrolib.core.kfactor` uses on the
high-flow side -- there is no low-flow generalized-skew map analogous to
B17C Plate 1 to justify carrying the approximation's error into a tail
where sample skew is already the least stable parameter. Pass
``distribution="lognormal"`` to force skew to zero, which is more stable on
short or noisy records at the cost of not fitting an asymmetric tail.

**Year definition.** Annual minima are computed over a *climatic year*
(April 1 - March 31) by default, following Riggs (1972), USGS Techniques of
Water-Resources Investigations Book 4 Chapter B1 -- neither the calendar
year nor the water year (Oct 1 - Sep 30) is used, because both can split a
single connected low-flow event across two years. A climatic year is
labeled by the calendar year containing the majority of its months (April 1
Y - March 31 Y+1 is climatic year Y), the same logic by which a water year
is labeled by its ending calendar year. Pass ``year_type="water"`` or
``"calendar"`` to use those instead.

**Zero-flow years.** A year whose annual minimum n-day mean is exactly zero
cannot contribute to a log-space fit -- log10(0) is undefined, and letting
it become -inf would corrupt every subsequent moment. This module never
takes that log: distribution parameters are fit to the positive-flow years
only, and zero years are folded back in through a conditional-probability
adjustment analogous to Bulletin 17B Appendix 5's treatment of censored low
outliers on the high-flow side. With ``p0 = n_zero_years / n_years``, a
requested non-exceedance probability ``p <= p0`` has quantile exactly 0 (that
fraction of years or more had no flow at all, so no positive flow value has
that low a non-exceedance probability); for ``p > p0`` the fit is queried at
the conditional probability ``(p - p0) / (1 - p0)``. ``n_zero_years`` and
``p0`` are always reported on :class:`hydrolib.core.LowFlowResults` so a
zero-inflated record is never silently indistinguishable from a well-behaved
one.

**Minimum record length.** Fewer than 10 years raises rather than returning
a fit: 10 systematic years is Bulletin 17C's own minimum for a flood
frequency analysis, and a low-flow quantile at the standard 10% or 2%
non-exceedance probability is no better supported by a shorter low-flow
record than a high-flow one is. Below 20 years a warning is logged rather
than raised -- 20 years gives roughly 2-4 observations at or below a 10%-90%
non-exceedance probability, which is little enough that the fitted quantile
should be read with real caution, but not so little the analysis is
meaningless.
"""

from __future__ import annotations

import logging
from typing import ClassVar, Optional

import numpy as np
import pandas as pd
from scipy.special import ndtri

from .core import (
    MAX_ABS_SKEW,
    YEAR_TYPES,
    LowFlowResults,
    assign_year_label,
    lp3_frequency_factor_peakfq,
)

logger = logging.getLogger(__name__)

#: Year definitions accepted by annual_minimum_flow and LowFlowFrequency.
#: Alias of hydrolib.core.YEAR_TYPES, kept under this name for backward
#: compatibility with the public name this module already exports.
LOW_FLOW_YEAR_TYPES = YEAR_TYPES


def annual_minimum_flow(
    daily_data: pd.DataFrame,
    n_day: int = 7,
    year_type: str = "climatic",
    min_days: int = 350,
) -> pd.DataFrame:
    """Compute the annual minimum n-day mean flow for each low-flow year.

    Parameters
    ----------
    daily_data : pd.DataFrame
        Daily mean flow indexed by date, with a ``flow_cfs`` column -- the
        shape returned by :meth:`hydrolib.usgs.USGSgage.download_daily_flow`.
        Negative values are treated as missing (see Notes), not as extreme
        lows.
    n_day : int
        Averaging window length in days, e.g. 7 for 7Q10/7Q2. Must be >= 1.
    year_type : str
        Year definition; see :data:`LOW_FLOW_YEAR_TYPES` and the module
        docstring. Default "climatic".
    min_days : int
        Minimum number of valid daily values a year must have to be marked
        ``complete``. Default 350, i.e. at most ~15 missing days.

    Returns
    -------
    pd.DataFrame
        One row per low-flow year that has at least one daily value,
        columns:

        - ``year`` : int, the low-flow year label
        - ``flow_cfs`` : float, minimum n-day mean flow in that year, or NaN
          if no full n-day window of valid data falls within the year
        - ``window_end_date`` : the last date of the window achieving that
          minimum (NaT if ``flow_cfs`` is NaN) -- useful for a reviewer
          checking when a reported low occurred
        - ``n_days`` : int, count of valid (non-missing, non-negative) daily
          values within that year's calendar dates
        - ``complete`` : bool, ``n_days >= min_days``
        - ``n_day`` : int, the window length used (repeated on every row so
          an exported CSV is self-describing without external context)

    Notes
    -----
    **Gaps.** :meth:`hydrolib.usgs.USGSgage.download_daily_flow` drops rows
    with no reported value, so a real gap in the record is a missing date,
    not a NaN-valued one. A plain ``rolling(window=n_day)`` would then treat
    ``n_day`` available *rows* as the window regardless of whether they are
    actually consecutive calendar days, silently bridging across the gap.
    This function reindexes onto a complete daily calendar first so that a
    gap becomes an explicit NaN and ``min_periods=n_day`` correctly refuses
    any window that touches it, rather than averaging non-adjacent days.

    **Negative values.** True streamflow cannot be negative; a negative
    daily value reflects a data artifact (commonly ice-affected or
    backwater conditions) rather than an extreme low, so it is treated as
    missing rather than as the year's candidate minimum. A true zero is not
    affected by this and is retained -- see the module docstring on
    zero-flow handling.

    **Window alignment.** Each n-day mean is labeled by the last day of its
    window (pandas' default rolling alignment). Since only the minimum value
    within a year is of interest, alignment does not change which n-day
    windows exist -- it only changes which date is attached to one -- except
    at a year boundary, where a window is credited to the year containing
    its last day.
    """
    if n_day < 1:
        raise ValueError(f"n_day must be >= 1, got {n_day}")
    if year_type not in LOW_FLOW_YEAR_TYPES:
        raise ValueError(f"year_type must be one of {LOW_FLOW_YEAR_TYPES}, got {year_type!r}")

    full_idx = pd.date_range(daily_data.index.min(), daily_data.index.max(), freq="D")
    flows = daily_data["flow_cfs"].reindex(full_idx)
    flows = flows.where(flows >= 0)

    n_day_mean = flows.rolling(window=n_day, min_periods=n_day).mean()
    window_years = assign_year_label(pd.DatetimeIndex(n_day_mean.index), year_type)
    daily_years = assign_year_label(full_idx, year_type)

    windows = pd.DataFrame(
        {"year": window_years, "flow_cfs": n_day_mean.to_numpy()}, index=n_day_mean.index
    ).dropna(subset=["flow_cfs"])

    if windows.empty:
        minima = pd.DataFrame(columns=["year", "flow_cfs", "window_end_date"])
    else:
        min_idx = windows.groupby("year")["flow_cfs"].idxmin()
        minima = windows.loc[min_idx].reset_index(names="window_end_date")
        minima = minima[["year", "flow_cfs", "window_end_date"]]

    completeness = (
        pd.DataFrame({"year": daily_years, "has_value": flows.notna().to_numpy()})
        .groupby("year")["has_value"]
        .sum()
        .rename("n_days")
    )

    result = minima.set_index("year").join(completeness, how="right")
    result["window_end_date"] = pd.to_datetime(result["window_end_date"])
    result["n_days"] = result["n_days"].astype(int)
    result["complete"] = result["n_days"] >= min_days
    result["n_day"] = n_day

    return result.reset_index(names="year").sort_values("year").reset_index(drop=True)


class LowFlowFrequency:
    """
    Low-flow frequency analysis: fits Log-Pearson III (or lognormal) to the
    annual minimum n-day mean flow and computes quantiles at specified
    non-exceedance probabilities.

    See the module docstring for the distribution, year-definition, and
    zero-flow method notes.

    Parameters
    ----------
    daily_data : pd.DataFrame
        Daily mean flow indexed by date, with a ``flow_cfs`` column.
    n_day : int
        Averaging window length in days. Default 7 (7Q10/7Q2).
    year_type : str
        Year definition; see :data:`LOW_FLOW_YEAR_TYPES`. Default
        "climatic".
    min_days : int
        Per-year completeness threshold passed to :func:`annual_minimum_flow`.
        Default 350.
    distribution : str
        "lp3" (station skew) or "lognormal" (skew forced to zero).
    require_complete_years : bool
        If True (default), years failing the ``min_days`` completeness
        check are excluded before fitting. A year missing exactly its
        lowest-flow period would otherwise report an annual minimum higher
        than the truth, biasing the fit toward understating how low flows
        actually get -- excluding it is the safer default.

    Raises
    ------
    ValueError
        ``distribution`` is not recognized, or fewer than
        :attr:`MIN_YEARS` usable years are available.

    Examples
    --------
    >>> lff = LowFlowFrequency(daily_df, n_day=7)
    >>> results = lff.run_analysis()
    >>> results.quantiles  # doctest: +SKIP
    >>> q7_10_and_7_2 = lff.compute_quantiles(np.array([0.10, 0.50]))  # doctest: +SKIP
    """

    STANDARD_NONEXCEEDANCE: ClassVar[np.ndarray] = np.array([0.50, 0.20, 0.10, 0.05, 0.02, 0.01])

    #: Bulletin 17C's own minimum systematic record length; a low-flow
    #: quantile is no better supported by a shorter record than a high-flow
    #: one is.
    MIN_YEARS: ClassVar[int] = 10
    #: Below this, a fit is produced but a warning is logged: ~2-4
    #: observations at a 10%-90% non-exceedance probability is thin support.
    RECOMMENDED_MIN_YEARS: ClassVar[int] = 20
    #: Absolute floor to compute a sample skew at all (matches
    #: MethodOfMoments' implicit requirement in bulletin17c.py).
    MIN_POSITIVE_YEARS: ClassVar[int] = 3

    def __init__(
        self,
        daily_data: pd.DataFrame,
        n_day: int = 7,
        year_type: str = "climatic",
        min_days: int = 350,
        distribution: str = "lp3",
        require_complete_years: bool = True,
    ):
        if distribution not in ("lp3", "lognormal"):
            raise ValueError(f"distribution must be 'lp3' or 'lognormal', got {distribution!r}")

        self._n_day = n_day
        self._year_type = year_type
        self._distribution = distribution

        annual = annual_minimum_flow(
            daily_data, n_day=n_day, year_type=year_type, min_days=min_days
        )
        usable = annual[annual["complete"]] if require_complete_years else annual
        usable = usable.dropna(subset=["flow_cfs"]).reset_index(drop=True)

        if len(usable) < self.MIN_YEARS:
            raise ValueError(
                f"Only {len(usable)} usable {year_type} year(s) of {n_day}-day low flow "
                f"available; at least {self.MIN_YEARS} are required for a low-flow "
                f"frequency fit. Pass require_complete_years=False to include incomplete "
                f"years if that is enough to reach the minimum (biases toward "
                f"understating how low flows get, since a year missing its lowest period "
                f"would report a falsely high minimum)."
            )
        if len(usable) < self.RECOMMENDED_MIN_YEARS:
            logger.warning(
                "Only %d years of %d-day low flow available; %d or more is recommended "
                "for a low-flow frequency fit. Quantiles at extreme non-exceedance "
                "probabilities (e.g. 0.01, 0.02) are supported by very few observations "
                "at this record length.",
                len(usable),
                n_day,
                self.RECOMMENDED_MIN_YEARS,
            )

        self._annual_minimums_all = annual
        self._annual_minimums = usable
        self._results: Optional[LowFlowResults] = None
        self._n_positive_years: Optional[int] = None

    @property
    def annual_minimums(self) -> pd.DataFrame:
        """Per-year table actually used for fitting (see require_complete_years)."""
        return self._annual_minimums.copy()

    @property
    def annual_minimums_all_years(self) -> pd.DataFrame:
        """Every year with any daily data, including incomplete ones excluded from the fit."""
        return self._annual_minimums_all.copy()

    @property
    def n(self) -> int:
        """Number of years used in the fit (zero-flow years included)."""
        return len(self._annual_minimums)

    @property
    def results(self) -> Optional[LowFlowResults]:
        return self._results

    def run_analysis(self) -> LowFlowResults:
        """Fit the distribution and compute quantiles and confidence limits
        at :attr:`STANDARD_NONEXCEEDANCE`.

        Returns
        -------
        LowFlowResults
        """
        flows = self._annual_minimums["flow_cfs"].to_numpy(dtype=float)
        n_years = len(flows)

        is_zero = flows <= 0.0
        n_zero = int(np.sum(is_zero))
        p_zero = n_zero / n_years
        positive = flows[~is_zero]
        n_pos = len(positive)

        if n_pos < self.MIN_POSITIVE_YEARS:
            raise ValueError(
                f"Only {n_pos} nonzero year(s) out of {n_years}; at least "
                f"{self.MIN_POSITIVE_YEARS} are required to fit a distribution to the "
                f"positive-flow years. A record this zero-dominated may not be suited to "
                f"this method at all -- consider whether the site is effectively "
                f"intermittent rather than perennial with occasional zero years."
            )

        log_flows = np.log10(positive)
        mean_log = float(np.mean(log_flows))
        std_log = float(np.std(log_flows, ddof=1))

        if std_log > 0:
            skew_station = (
                n_pos
                * np.sum((log_flows - mean_log) ** 3)
                / ((n_pos - 1) * (n_pos - 2) * std_log**3)
            )
            skew_station = float(np.clip(skew_station, -MAX_ABS_SKEW, MAX_ABS_SKEW))
        else:
            # Zero variance (every positive year has the identical annual minimum --
            # plausible on a reach governed by a fixed minimum-flow release) makes
            # the sample-skew formula a 0/0 division. np.clip does not clip NaN, so
            # leaving this unguarded would silently poison skew_used, and from there
            # every K-factor and quantile downstream. Zero variance implies zero
            # skew, and every quantile correctly collapses to the single observed
            # value regardless of skew once std_log is 0, so this is the correct
            # value here, not just a safe placeholder.
            skew_station = 0.0
        skew_used = 0.0 if self._distribution == "lognormal" else skew_station

        self._n_positive_years = n_pos
        self._results = LowFlowResults(
            n_years=n_years,
            n_zero_years=n_zero,
            p_zero=p_zero,
            n_day=self._n_day,
            year_type=self._year_type,
            distribution=self._distribution,
            mean_log=mean_log,
            std_log=std_log,
            skew_station=skew_station,
            skew_used=skew_used,
        )
        self._results.quantiles = self.compute_quantiles()
        self._results.confidence_limits = self.compute_confidence_limits()
        return self._results

    def _validate_probabilities(self, non_exceedance: Optional[np.ndarray]) -> np.ndarray:
        p = (
            self.STANDARD_NONEXCEEDANCE
            if non_exceedance is None
            else np.asarray(non_exceedance, dtype=float)
        )
        if np.any((p <= 0) | (p >= 1)):
            raise ValueError("non_exceedance probabilities must be strictly between 0 and 1")
        return p

    def compute_quantiles(self, non_exceedance: np.ndarray = None) -> pd.DataFrame:
        """Compute low-flow quantiles at the given non-exceedance probabilities.

        Parameters
        ----------
        non_exceedance : np.ndarray, optional
            Non-exceedance probabilities, e.g. 0.10 for 7Q10, 0.50 for 7Q2.
            Defaults to :attr:`STANDARD_NONEXCEEDANCE`.

        Returns
        -------
        pd.DataFrame
            Columns: ``non_exceedance_prob``, ``return_period`` (``1/p``),
            ``flow_cfs``, ``log_flow`` (NaN, never -inf, where ``flow_cfs``
            is exactly 0), ``K_factor``, ``conditional_prob`` (NaN where
            ``p <= p_zero``; see the module docstring on zero-flow years).
        """
        if self._results is None:
            self.run_analysis()

        p = self._validate_probabilities(non_exceedance)
        p0 = self._results.p_zero
        below = p <= p0
        p_cond = np.where(below, np.nan, (p - p0) / (1 - p0))

        K = np.full_like(p, np.nan)
        K[~below] = [
            lp3_frequency_factor_peakfq(pc, self._results.skew_used) for pc in p_cond[~below]
        ]

        log_Q = self._results.mean_log + K * self._results.std_log
        flow_cfs = np.where(below, 0.0, 10**log_Q)
        log_flow = np.where(below, np.nan, log_Q)

        return pd.DataFrame(
            {
                "non_exceedance_prob": p,
                "return_period": 1 / p,
                "flow_cfs": flow_cfs,
                "log_flow": log_flow,
                "K_factor": K,
                "conditional_prob": p_cond,
            }
        )

    def compute_confidence_limits(
        self, non_exceedance: np.ndarray = None, confidence: float = 0.90
    ) -> pd.DataFrame:
        """Compute confidence limits for the quantile estimates.

        Uses the same asymptotic LP3 variance formula as
        :meth:`hydrolib.bulletin17c.FloodFrequencyAnalysis.compute_confidence_limits`,
        applied to the positive-year fit at the conditional probability. This
        does not account for the added uncertainty in ``p_zero`` itself, so
        it understates true uncertainty for a record with any zero years;
        treat it as a lower bound on the actual interval width. Limits are
        not reported (NaN) for ``p <= p_zero``, where the point estimate is
        the deterministic value 0 rather than a value with sampling
        uncertainty in the usual sense.

        Parameters
        ----------
        non_exceedance : np.ndarray, optional
            Defaults to :attr:`STANDARD_NONEXCEEDANCE`.
        confidence : float
            Confidence level, default 0.90 (matching the high-flow side's
            default).

        Returns
        -------
        pd.DataFrame
            Columns: ``non_exceedance_prob``, ``return_period``,
            ``flow_cfs``, ``lower_{pct}pct``, ``upper_{pct}pct``.
        """
        if self._results is None:
            self.run_analysis()

        p = self._validate_probabilities(non_exceedance)
        p0 = self._results.p_zero
        below = p <= p0
        p_cond = np.where(below, np.nan, (p - p0) / (1 - p0))

        G = self._results.skew_used
        n_pos = self._n_positive_years
        alpha = 1 - confidence
        z_alpha = ndtri(1 - alpha / 2)

        K = np.full_like(p, np.nan)
        K[~below] = [lp3_frequency_factor_peakfq(pc, G) for pc in p_cond[~below]]

        var_factor = 1 / n_pos + K**2 * (1 + 0.75 * G**2) / (2 * (n_pos - 1))
        se_log = self._results.std_log * np.sqrt(var_factor)
        log_Q = self._results.mean_log + K * self._results.std_log

        flow_cfs = np.where(below, 0.0, 10**log_Q)
        lower = np.where(below, np.nan, 10 ** (log_Q - z_alpha * se_log))
        upper = np.where(below, np.nan, 10 ** (log_Q + z_alpha * se_log))
        pct = int(confidence * 100)

        return pd.DataFrame(
            {
                "non_exceedance_prob": p,
                "return_period": 1 / p,
                "flow_cfs": flow_cfs,
                f"lower_{pct}pct": lower,
                f"upper_{pct}pct": upper,
            }
        )
