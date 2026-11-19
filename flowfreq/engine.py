"""
flowfreq.engine - a plain moments fit kept from early versions. Deprecated.

``B17CEngine`` fits log-Pearson Type III by the sample moments of the peaks it is given
and nothing more: no Expected Moments Algorithm, no low-outlier test, no historical or
censored data and no regional skew. Its name says Bulletin 17C and its answer is not a
Bulletin 17C answer except on a complete record with no low outliers and station skew.
Use :class:`flowfreq.Bulletin17C` or :func:`flowfreq.run_ffa`. ``flowfreq.batch`` and
``flowfreq.plots`` are built on this class and are deprecated with it.
"""

from __future__ import annotations

import warnings
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from .core import PeakRecord, compute_ci_lp3, lp3_frequency_factor_peakfq, lp3_quantile_peakfq

# Standard return periods used in flood frequency analysis
STANDARD_RETURN_PERIODS = (1.5, 2, 5, 10, 25, 50, 100, 200, 500)


#: What the deprecation warning says, here and in ``flowfreq.batch`` and ``flowfreq.plots``.
DEPRECATION = (
    "flowfreq.engine.B17CEngine (and flowfreq.batch and flowfreq.plots, which are built on "
    "it) is deprecated and will be removed in a later release. It is a plain moments fit: "
    "no EMA, no low-outlier test, no regional skew. Use flowfreq.Bulletin17C or "
    "flowfreq.run_ffa, and flowfreq.freq_plot.plot_frequency_curve for the figure."
)


class B17CEngine:
    """
    A log-Pearson Type III fit by sample moments. Deprecated.

    Not a Bulletin 17C analysis: see the module docstring. Creating one warns.

    Examples
    --------
    >>> from flowfreq.core import PeakRecord
    >>> from flowfreq.engine import B17CEngine
    >>> records = [PeakRecord(year=y, flow=f) for y, f in zip(years, flows)]
    >>> engine = B17CEngine()
    >>> engine.fit(records)
    >>> quantiles = engine.quantiles([10, 50, 100])
    >>> ci = engine.quantiles_with_ci([10, 50, 100])
    """

    def __init__(self) -> None:
        warnings.warn(DEPRECATION, DeprecationWarning, stacklevel=2)
        self.params: Optional[Tuple[float, float, float]] = None
        self.n: int = 0
        self._records: List[PeakRecord] = []

    @property
    def mu(self) -> Optional[float]:
        """Mean of log10(flow)."""
        return self.params[0] if self.params else None

    @property
    def sigma(self) -> Optional[float]:
        """Standard deviation of log10(flow)."""
        return self.params[1] if self.params else None

    @property
    def skew(self) -> Optional[float]:
        """Skew coefficient."""
        return self.params[2] if self.params else None

    def fit(self, records: Sequence[PeakRecord]) -> Tuple[float, float, float]:
        """
        Fit LP3 distribution to peak flow records.

        Parameters
        ----------
        records : sequence of PeakRecord
            Peak flow records to analyze

        Returns
        -------
        tuple
            (mu, sigma, skew) - fitted LP3 parameters
        """
        self._records = list(records)
        obs = [r.flow for r in records if r.flow is not None]

        if len(obs) < 3:
            raise ValueError("At least 3 non-null flow values required for fitting")

        log_flows = np.log10(obs)
        n = len(obs)
        mu = float(log_flows.mean())
        sigma = float(log_flows.std(ddof=1))
        # Bulletin 17C Eq. 7-2, the unbiased sample estimator -- the same form
        # Bulletin17C.run_analysis uses, so the two paths in this library agree
        # on the station skew. Before 0.4.0 this was the biased population
        # coefficient, sum((x-mu)**3)/n / sigma**3, whose magnitude is smaller
        # by a factor of n**2 / ((n-1)(n-2)): 7.2% at n=44, 39% at n=10.
        skew = float(n * np.sum((log_flows - mu) ** 3) / ((n - 1) * (n - 2) * sigma**3))

        self.params = (mu, sigma, skew)
        self.n = n

        return self.params

    def fit_from_flows(
        self, flows: Sequence[float], years: Optional[Sequence[int]] = None
    ) -> Tuple[float, float, float]:
        """
        Fit LP3 distribution directly from flow values.

        Parameters
        ----------
        flows : sequence of float
            Peak flow values in cfs
        years : sequence of int, optional
            Corresponding water years

        Returns
        -------
        tuple
            (mu, sigma, skew) - fitted LP3 parameters
        """
        if years is None:
            years = range(len(flows))

        records = [
            PeakRecord(year=int(y), flow=float(f))
            for y, f in zip(years, flows)
            if f is not None and not np.isnan(f)
        ]
        return self.fit(records)

    def quantiles(
        self, return_periods: Sequence[float] = STANDARD_RETURN_PERIODS
    ) -> Dict[float, float]:
        """
        Compute flood quantiles for specified return periods.

        Parameters
        ----------
        return_periods : sequence of float
            Return periods in years (default: standard set)

        Returns
        -------
        dict
            Mapping of return period to flow quantile in cfs
        """
        if self.params is None:
            raise RuntimeError("Must call fit() before computing quantiles")

        mu, sigma, skew = self.params
        return {T: lp3_quantile_peakfq(mu, sigma, skew, T) for T in return_periods}

    def quantiles_with_ci(
        self,
        return_periods: Sequence[float] = STANDARD_RETURN_PERIODS,
        alpha: float = 0.05,
    ) -> Dict[float, Dict[str, float]]:
        """
        Compute flood quantiles with confidence intervals.

        Parameters
        ----------
        return_periods : sequence of float
            Return periods in years
        alpha : float
            Significance level (default 0.05 for 95% CI)

        Returns
        -------
        dict
            Mapping of return period to dict with 'estimate', 'lower', 'upper'
        """
        if self.params is None:
            raise RuntimeError("Must call fit() before computing quantiles")

        mu, sigma, skew = self.params
        result = {}

        for T in return_periods:
            estimate, lower, upper = compute_ci_lp3(mu, sigma, skew, self.n, T, alpha)
            result[T] = {"estimate": estimate, "lower": lower, "upper": upper}

        return result

    def frequency_table(
        self,
        return_periods: Sequence[float] = STANDARD_RETURN_PERIODS,
        alpha: float = 0.05,
    ) -> "pd.DataFrame":
        """
        Generate a frequency table with quantiles and confidence intervals.

        Parameters
        ----------
        return_periods : sequence of float
            Return periods in years
        alpha : float
            Significance level for CI

        Returns
        -------
        pd.DataFrame
            Frequency table with columns: Return Period, AEP, Flow (cfs),
            Lower CI, Upper CI
        """

        ci = self.quantiles_with_ci(return_periods, alpha)

        data = []
        for T in return_periods:
            data.append(
                {
                    "Return Period (yr)": T,
                    "AEP (%)": 100 / T,
                    "Flow (cfs)": ci[T]["estimate"],
                    f"Lower {int((1-alpha)*100)}%": ci[T]["lower"],
                    f"Upper {int((1-alpha)*100)}%": ci[T]["upper"],
                }
            )

        return pd.DataFrame(data)

    def summary(self) -> str:
        """Return a summary of the fitted model."""
        if self.params is None:
            return "Model not fitted. Call fit() first."

        mu, sigma, skew = self.params
        lines = [
            "B17C Engine - LP3 Fit Summary",
            "=" * 35,
            f"Sample size (n):     {self.n}",
            f"Mean (log10 Q):      {mu:.4f}",
            f"Std Dev (log10 Q):   {sigma:.4f}",
            f"Skew coefficient:    {skew:.4f}",
            "",
            "Key Quantiles:",
            f"  Q10  (10-yr):  {lp3_quantile_peakfq(mu, sigma, skew, 10):,.0f} cfs",
            f"  Q50  (50-yr):  {lp3_quantile_peakfq(mu, sigma, skew, 50):,.0f} cfs",
            f"  Q100 (100-yr): {lp3_quantile_peakfq(mu, sigma, skew, 100):,.0f} cfs",
        ]
        return "\n".join(lines)
