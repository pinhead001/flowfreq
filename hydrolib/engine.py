"""
hydrolib.engine - B17C Engine for flood frequency analysis

Provides a simplified interface for Bulletin 17C analysis with PeakFQ-style output.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .core import (
    PeakRecord,
    compute_ci_lp3,
    lp3_frequency_factor_peakfq,
    lp3_quantile_peakfq,
)

# Standard return periods used in flood frequency analysis
STANDARD_RETURN_PERIODS = (1.5, 2, 5, 10, 25, 50, 100, 200, 500)


class B17CEngine:
    """
    Simplified Bulletin 17C engine for flood frequency analysis.

    Provides PeakFQ-style quantile estimation with confidence intervals.

    Examples
    --------
    >>> from hydrolib import B17CEngine, PeakRecord
    >>> records = [PeakRecord(year=y, flow=f) for y, f in zip(years, flows)]
    >>> engine = B17CEngine()
    >>> engine.fit(records)
    >>> quantiles = engine.quantiles([10, 50, 100])
    >>> ci = engine.quantiles_with_ci([10, 50, 100])
    """

    def __init__(self):
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
        mu = float(log_flows.mean())
        sigma = float(log_flows.std(ddof=1))
        skew = float(((log_flows - mu) ** 3).mean() / (sigma**3))

        self.params = (mu, sigma, skew)
        self.n = len(obs)

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
        import pandas as pd

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
