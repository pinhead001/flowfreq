"""
hydrolib.bulletin17c - Bulletin 17C Flood Frequency Analysis

Implements both Method of Moments (MOM) and Expected Moments Algorithm (EMA).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from datetime import datetime
from functools import cached_property, lru_cache
from typing import ClassVar, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.integrate import quad as _quad
from scipy.special import gammainc, gammaincc, gammaln, ndtri

from .core import (
    MAX_ABS_SKEW,
    AnalysisMethod,
    EMAParameters,
    FlowInterval,
    FrequencyResults,
    grubbs_beck_critical_value,
    kfactor,
    kfactor_array,
    kfactor_skew_derivative,
    log_pearson3_cdf,
    log_pearson3_pdf,
)


class FloodFrequencyAnalysis(ABC):
    """Abstract base class for flood frequency analysis methods."""

    STANDARD_AEP: ClassVar[np.ndarray] = np.array(
        [0.995, 0.99, 0.95, 0.90, 0.80, 0.67, 0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002]
    )

    def __init__(
        self, peak_flows: np.ndarray, regional_skew: float = None, regional_skew_mse: float = None
    ):
        self._peak_flows = np.array(peak_flows)
        self._peak_flows = self._peak_flows[~np.isnan(self._peak_flows)]
        self._n_zeros = int(np.sum(self._peak_flows == 0))
        self._peak_flows = self._peak_flows[self._peak_flows > 0]

        self._regional_skew = regional_skew
        self._regional_skew_mse = regional_skew_mse
        self._results: Optional[FrequencyResults] = None

    @property
    def n(self) -> int:
        return len(self._peak_flows)

    @property
    def peak_flows(self) -> np.ndarray:
        return self._peak_flows.copy()

    @cached_property
    def log_flows(self) -> np.ndarray:
        return np.log10(self._peak_flows)

    @property
    def results(self) -> Optional[FrequencyResults]:
        return self._results

    @abstractmethod
    def run_analysis(self) -> FrequencyResults:
        pass

    def _compute_weighted_skew(
        self, station_skew: float, n_effective: Optional[int] = None
    ) -> Optional[float]:
        """Compute weighted skew from station and regional values.

        See `_compute_skew_weighting` for the underlying formula and MSE.
        Retained for backward compatibility; prefer `_compute_skew_weighting`
        when the MSE of the resulting skew is also needed (e.g. for CI width).
        """
        skew_weighted, _ = self._compute_skew_weighting(station_skew, n_effective)
        return skew_weighted

    def _compute_skew_weighting(
        self, station_skew: float, n_effective: Optional[int] = None
    ) -> Tuple[Optional[float], float]:
        """Compute weighted skew and its mean-square error.

        Uses the Bulletin 17C Appendix 4 variance formula (eq. A4-2):
        V(G) = [6N(N-1) / ((N-2)(N+1)(N+3))] * [1 + (6/N)G^2 + (15/N^2)G^4 + ...]

        When a regional skew is supplied, the weighted skew's MSE is the
        harmonic combination of the station and regional MSEs:
        V(G_weighted) = mse_station * mse_regional / (mse_station + mse_regional)
        (the variance-minimizing combination of two independent estimators).
        When no regional skew is supplied, the MSE returned is simply the
        station skew's own sampling variance.

        Parameters
        ----------
        station_skew : float
            Station skew coefficient.
        n_effective : int, optional
            Effective sample size for variance calculation. If None, uses self.n.
            For EMA with historical data, this should be the count of point
            observations (systematic + historical peaks).

        Returns
        -------
        tuple[float or None, float]
            (weighted_skew, skew_mse). weighted_skew is None if no regional
            skew was supplied; skew_mse is always the MSE of whichever skew
            (station or weighted) is ultimately used.
        """
        n = n_effective if n_effective is not None else self.n
        # B17C exact variance of sample skew (Appendix 4, eq. A4-2)
        base_var = (6 * n * (n - 1)) / ((n - 2) * (n + 1) * (n + 3))
        mse_station = base_var * (1 + (6 / n) * station_skew**2 + (15 / (n**2)) * station_skew**4)

        if self._regional_skew is None or self._regional_skew_mse is None:
            return None, mse_station

        w_regional = mse_station / (mse_station + self._regional_skew_mse)
        w_station = self._regional_skew_mse / (mse_station + self._regional_skew_mse)

        skew_weighted = w_station * station_skew + w_regional * self._regional_skew
        mse_weighted = (mse_station * self._regional_skew_mse) / (
            mse_station + self._regional_skew_mse
        )
        return skew_weighted, mse_weighted

    def compute_quantiles(self, aep: np.ndarray = None) -> pd.DataFrame:
        """Compute flood frequency quantiles."""
        if self._results is None:
            self.run_analysis()

        if aep is None:
            aep = self.STANDARD_AEP

        K = kfactor_array(self._results.skew_used, aep)
        log_Q = self._results.mean_log + K * self._results.std_log
        Q = 10**log_Q

        return pd.DataFrame(
            {"aep": aep, "return_period": 1 / aep, "flow_cfs": Q, "log_flow": log_Q, "K_factor": K}
        )

    def compute_confidence_limits(
        self, aep: np.ndarray = None, confidence: float = 0.90
    ) -> pd.DataFrame:
        """Compute confidence limits for quantile estimates."""
        if self._results is None:
            self.run_analysis()

        if aep is None:
            aep = self.STANDARD_AEP

        alpha = 1 - confidence
        z_alpha = ndtri(1 - alpha / 2)

        n = self._results.n_systematic or self.n
        K = kfactor_array(self._results.skew_used, aep)
        G = self._results.skew_used
        std_log = self._results.std_log

        # Classic Bulletin 17B/17C approximate variance of the LP3 quantile
        # estimate, in two parts:
        #  1. Sampling variance of mean/std, propagated through K (eq. B17B
        #     App. 9): Var/S^2 = (1/n) * [1 + K*G + (K^2/2)*(1 + 0.75*G^2)]
        #  2. Sampling variance of the skew itself (skew_used_mse), propagated
        #     through the sensitivity of K to G via (dK/dG)^2. This term is
        #     small for near-median AEPs but dominates at extreme return
        #     periods, where K is large and most sensitive to skew.
        moment_var_factor = (1 / n) * (1 + K * G + (K**2 / 2) * (1 + 0.75 * G**2))
        moment_var_factor = np.maximum(moment_var_factor, 0.0)
        var_log = std_log**2 * moment_var_factor

        skew_mse = self._results.skew_used_mse
        if skew_mse is not None and skew_mse > 0:
            dK_dG = np.array([kfactor_skew_derivative(G, float(p)) for p in aep])
            var_log = var_log + (dK_dG**2) * (std_log**2) * skew_mse

        se_log = np.sqrt(var_log)

        log_Q = self._results.mean_log + K * self._results.std_log
        log_lower = log_Q - z_alpha * se_log
        log_upper = log_Q + z_alpha * se_log

        return pd.DataFrame(
            {
                "aep": aep,
                "return_period": 1 / aep,
                "flow_cfs": 10**log_Q,
                "lower_5pct": 10**log_lower,
                "upper_5pct": 10**log_upper,
            }
        )

    def plot_frequency_curve(
        self,
        site_name: str = None,
        site_no: str = None,
        save_path: str = None,
        figsize: Tuple[int, int] = (10, 7),
        show_confidence: bool = True,
        yscale: str = "log",
    ) -> plt.Figure:
        """Plot professional flood frequency curve."""
        if self._results is None:
            self.run_analysis()

        fig, ax = plt.subplots(figsize=figsize)

        def prob_to_x(p):
            return ndtri(1 - p)

        # Plot observed data
        sorted_flows = np.sort(self._peak_flows)[::-1]
        n = len(sorted_flows)
        plotting_position = np.arange(1, n + 1) / (n + 1)

        x_data = [prob_to_x(p) for p in plotting_position]
        ax.scatter(
            x_data,
            sorted_flows,
            c="blue",
            s=40,
            zorder=5,
            label="Observed Annual Peaks",
            edgecolors="darkblue",
            linewidth=0.5,
        )

        # Plot fitted curve
        aep_curve = np.array(
            [
                0.999,
                0.995,
                0.99,
                0.95,
                0.90,
                0.80,
                0.50,
                0.20,
                0.10,
                0.04,
                0.02,
                0.01,
                0.005,
                0.002,
                0.001,
            ]
        )
        K_curve = kfactor_array(self._results.skew_used, aep_curve)
        Q_curve = 10 ** (self._results.mean_log + K_curve * self._results.std_log)
        x_curve = [prob_to_x(p) for p in aep_curve]

        method_label = "EMA" if self._results.method == AnalysisMethod.EMA else "MOM"
        ax.plot(
            x_curve,
            Q_curve,
            "b-",
            linewidth=2,
            label=f"LP3 Fitted Curve ({method_label})",
            zorder=4,
        )

        # Confidence limits
        if show_confidence:
            cl = self.compute_confidence_limits(aep_curve)
            x_cl = [prob_to_x(p) for p in cl["aep"]]
            ax.fill_between(
                x_cl,
                cl["lower_5pct"],
                cl["upper_5pct"],
                alpha=0.2,
                color="blue",
                label="90% Confidence Interval",
            )
            ax.plot(x_cl, cl["lower_5pct"], "b--", linewidth=1, alpha=0.7)
            ax.plot(x_cl, cl["upper_5pct"], "b--", linewidth=1, alpha=0.7)

        # Low outliers
        if self._results.low_outlier_threshold > 0 and self._results.n_low_outliers > 0:
            threshold = self._results.low_outlier_threshold
            low_outliers = self._peak_flows[self._peak_flows < threshold]
            for lo in low_outliers:
                idx = np.where(sorted_flows == lo)[0]
                if len(idx) > 0:
                    ax.scatter(
                        [x_data[i] for i in idx],
                        [lo] * len(idx),
                        c="red",
                        s=60,
                        marker="x",
                        zorder=6,
                    )
            ax.axhline(
                threshold,
                color="red",
                linestyle=":",
                alpha=0.7,
                label=f"Low Outlier Threshold ({self._results.n_low_outliers})",
            )

        ax.set_yscale(yscale)
        ax.set_ylabel("Peak Discharge (cfs)", fontsize=11)
        ax.set_xlabel("Annual Exceedance Probability", fontsize=11)

        # X-axis ticks
        prob_ticks = [0.99, 0.95, 0.90, 0.80, 0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002]
        x_ticks = [prob_to_x(p) for p in prob_ticks]
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(
            [f"{p*100:.1f}%" if p >= 0.01 else f"{p*100:.2f}%" for p in prob_ticks],
            fontsize=9,
            rotation=45,
            ha="right",
        )

        # Secondary axis for return period
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        rp_ticks = [1.01, 2, 5, 10, 25, 50, 100, 200, 500]
        rp_x = [prob_to_x(1 / rp) for rp in rp_ticks]
        ax2.set_xticks(rp_x)
        ax2.set_xticklabels(
            [str(int(rp)) if rp >= 1 else f"{rp:.2f}" for rp in rp_ticks], fontsize=9
        )
        ax2.set_xlabel("Return Period (years)", fontsize=11)

        # Title
        title = f"Flood Frequency Analysis (Bulletin 17C - {method_label})"
        if site_name and site_no:
            title = f"Flood Frequency Analysis (Bulletin 17C - {method_label})\nUSGS {site_no} - {site_name}"
        ax.set_title(title, fontsize=12, fontweight="bold", pad=35)

        # Stats annotation
        r = self._results
        stats_text = f"n = {r.n_peaks}"
        if r.n_historical > 0:
            stats_text += f" ({r.n_systematic} sys + {r.n_historical} hist)"
        stats_text += f"\nMean(log Q) = {r.mean_log:.4f}"
        stats_text += f"\nStd(log Q) = {r.std_log:.4f}"
        stats_text += f"\nStation Skew = {r.skew_station:.3f}"
        if r.skew_weighted is not None:
            stats_text += f"\nWeighted Skew = {r.skew_weighted:.3f}"

        ax.annotate(
            stats_text,
            xy=(0.02, 0.98),
            xycoords="axes fraction",
            fontsize=9,
            ha="left",
            va="top",
            family="monospace",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.9),
        )

        ax.legend(loc="lower right", fontsize=9)
        ax.grid(True, which="both", alpha=0.3)
        ax.set_xlim(prob_to_x(0.999), prob_to_x(0.001))

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches="tight")

        return fig


class MethodOfMoments(FloodFrequencyAnalysis):
    """Traditional Method of Moments flood frequency analysis."""

    def run_analysis(self) -> FrequencyResults:
        """Run Method of Moments analysis."""
        log_flows = self.log_flows
        n = self.n

        mean_log = np.mean(log_flows)
        std_log = np.std(log_flows, ddof=1)

        skew_station = n * np.sum((log_flows - mean_log) ** 3) / ((n - 1) * (n - 2) * std_log**3)
        skew_station = float(np.clip(skew_station, -MAX_ABS_SKEW, MAX_ABS_SKEW))

        skew_weighted, skew_used_mse = self._compute_skew_weighting(skew_station)
        skew_used = skew_weighted if skew_weighted is not None else skew_station

        k_n = grubbs_beck_critical_value(n)
        threshold_log = mean_log - k_n * std_log
        low_outlier_threshold = 10**threshold_log
        n_low_outliers = int(np.sum(self._peak_flows < low_outlier_threshold))

        self._results = FrequencyResults(
            n_peaks=n,
            n_systematic=n,
            n_historical=0,
            n_censored=0,
            n_low_outliers=n_low_outliers,
            mean_log=mean_log,
            std_log=std_log,
            skew_station=skew_station,
            skew_regional=self._regional_skew,
            skew_weighted=skew_weighted,
            skew_used=skew_used,
            skew_used_mse=skew_used_mse,
            low_outlier_threshold=low_outlier_threshold,
            mgb_critical_value=k_n,
            method=AnalysisMethod.MOM,
        )

        self._results.quantiles = self.compute_quantiles()
        self._results.confidence_limits = self.compute_confidence_limits()

        return self._results


class ExpectedMomentsAlgorithm(FloodFrequencyAnalysis):
    """
    Bulletin 17C Expected Moments Algorithm (EMA) implementation.

    Handles systematic record, historical floods, censored observations,
    and low outliers identified by Multiple Grubbs-Beck test.
    """

    def __init__(
        self,
        peak_flows: np.ndarray,
        water_years: np.ndarray = None,
        regional_skew: float = None,
        regional_skew_mse: float = None,
        ema_params: EMAParameters = None,
        historical_peaks: List[Tuple[int, float]] = None,
        perception_thresholds: Dict[Tuple[int, int], float] = None,
        user_low_outlier_threshold: Optional[float] = None,
    ):
        super().__init__(peak_flows, regional_skew, regional_skew_mse)

        if water_years is not None:
            self._water_years = np.array(water_years)
        else:
            end_year = datetime.now().year
            self._water_years = np.arange(end_year - len(peak_flows) + 1, end_year + 1)

        self._historical_peaks = historical_peaks or []
        self._perception_thresholds = perception_thresholds or {}
        self._user_low_outlier_threshold = user_low_outlier_threshold
        self._ema_params = ema_params or self._auto_configure_ema_params()
        self._intervals: List[FlowInterval] = []

    def _auto_configure_ema_params(self) -> EMAParameters:
        """Auto-configure EMA parameters from data."""
        sys_start = int(self._water_years.min())
        sys_end = int(self._water_years.max())

        all_years = set(range(sys_start, sys_end + 1))
        recorded_years = set(self._water_years.astype(int))
        gaps = sorted(all_years - recorded_years)

        hist_start = None
        hist_end = None
        hist_threshold = None

        # An explicitly provided perception threshold for the historical period is
        # authoritative and takes precedence over the derivation below, which
        # conflates two different quantities: "the biggest flood we happen to have a
        # record of" (max of historical peak values) is not the same as "the smallest
        # flood that would have been noticed at all" (the actual perception
        # threshold). Confusing the two understates how informative the
        # pre-systematic record really is -- every unlisted year in the historical
        # period gets censored against the wrong (larger, weaker) bound -- and
        # materially biases every downstream moment. self._perception_thresholds was
        # accepted and stored by __init__, but nothing previously read it back out
        # here; it only reached get_perception_thresholds_table(), a display-only
        # method that never touched the fit itself.
        #
        # EMAParameters carries a single historical threshold, so multiple historical
        # perception periods must be collapsed; min() is used because the lowest
        # threshold is the one that binds.
        for (start, end), threshold in self._perception_thresholds.items():
            if end < sys_start and threshold > 0:
                # This is a historical perception period
                if hist_start is None or start < hist_start:
                    hist_start = start
                if hist_end is None or end > hist_end:
                    hist_end = end
                # Use the perception threshold (not max flow)
                if hist_threshold is None or threshold < hist_threshold:
                    hist_threshold = threshold

        # If historical peaks exist but no explicit thresholds, derive from peaks
        if self._historical_peaks and hist_threshold is None:
            hist_years = [h[0] for h in self._historical_peaks]
            if hist_start is None:
                hist_start = min(hist_years)
            if hist_end is None:
                hist_end = max(max(hist_years), sys_start - 1)
            hist_threshold = max(h[1] for h in self._historical_peaks)
        elif gaps and hist_threshold is None:
            first_gap = gaps[0]
            pre_gap_mask = self._water_years < first_gap
            if np.any(pre_gap_mask):
                pre_gap_flows = self._peak_flows[pre_gap_mask]
                hist_threshold = float(np.max(pre_gap_flows))
                hist_start = sys_start
                hist_end = first_gap - 1

        return EMAParameters(
            systematic_start=sys_start,
            systematic_end=sys_end,
            historical_start=hist_start,
            historical_end=hist_end,
            historical_threshold=hist_threshold,
        )

    def _build_flow_intervals(self, low_threshold: float = 0.0) -> List[FlowInterval]:
        """Build flow intervals for EMA analysis."""
        intervals = []

        for flow, year in zip(self._peak_flows, self._water_years):
            year = int(year)

            if flow < low_threshold:
                intervals.append(
                    FlowInterval.from_censored(
                        lower=0, upper=low_threshold, year=year, perception_threshold=0.0
                    )
                )
            else:
                intervals.append(FlowInterval.from_peak(flow, year))

        for year, flow in self._historical_peaks:
            threshold = self._ema_params.historical_threshold or flow
            intervals.append(
                FlowInterval.from_historical(flow, year, perception_threshold=threshold)
            )

        if self._ema_params.historical_start and self._ema_params.historical_threshold:
            hist_recorded_years = {h[0] for h in self._historical_peaks}
            sys_years = set(self._water_years.astype(int))

            for year in range(
                self._ema_params.historical_start, self._ema_params.historical_end + 1
            ):
                if year not in hist_recorded_years and year not in sys_years:
                    intervals.append(
                        FlowInterval.from_censored(
                            lower=0,
                            upper=self._ema_params.historical_threshold,
                            year=year,
                            perception_threshold=self._ema_params.historical_threshold,
                        )
                    )

        self._intervals = sorted(intervals, key=lambda x: x.year)
        return self._intervals

    @staticmethod
    def _truncated_gamma_moment(alpha: float, lower: float, upper: float, k: int) -> float:
        """Compute E[X^k] for a standardized Gamma(alpha,1) truncated to [lower, upper].

        Uses the identity E[X^k | lower<X<upper] = Gamma(alpha+k)/Gamma(alpha)
        * [P(alpha+k, upper) - P(alpha+k, lower)] / [P(alpha, upper) - P(alpha, lower)]
        where P is the regularized incomplete gamma function.

        Parameters
        ----------
        alpha : float
            Shape parameter of the gamma distribution.
        lower : float
            Lower truncation bound (>=0).
        upper : float
            Upper truncation bound (can be inf).
        k : int
            Moment order (1, 2, or 3).

        Returns
        -------
        float
            The k-th raw moment of the truncated gamma.
        """
        # Regularized incomplete gamma: gammainc(a, x) = P(a, x)
        if np.isinf(upper):
            p_num = gammaincc(alpha + k, lower) if lower > 0 else 1.0
            p_den = gammaincc(alpha, lower) if lower > 0 else 1.0
        else:
            p_num = gammainc(alpha + k, upper) - (gammainc(alpha + k, lower) if lower > 0 else 0.0)
            p_den = gammainc(alpha, upper) - (gammainc(alpha, lower) if lower > 0 else 0.0)

        if p_den < 1e-30:
            return 0.0

        # ADJ = Gamma(alpha+k)/Gamma(alpha) = alpha*(alpha+1)*...*(alpha+k-1)
        adj = 1.0
        for j in range(k):
            adj *= alpha + j

        return adj * p_num / p_den

    @staticmethod
    def _truncated_normal_moments(zl: float, zu: float) -> Tuple[float, float, float]:
        """Compute E[Z^k | zl<Z<zu] for k=1,2,3 where Z ~ N(0,1).

        Uses the recurrence: E[Z^k] = (k-1)*E[Z^{k-2}] - (zu^{k-1}*phi(zu) - zl^{k-1}*phi(zl))/(Phi(zu)-Phi(zl))
        """
        phi_u = stats.norm.pdf(zu)
        phi_l = stats.norm.pdf(zl)
        cdf_u = stats.norm.cdf(zu)
        cdf_l = stats.norm.cdf(zl)
        dp = cdf_u - cdf_l
        if dp < 1e-30:
            mid = (zl + zu) / 2.0
            return mid, mid**2, mid**3

        ez1 = -(phi_u - phi_l) / dp
        ez2 = 1.0 + (-(zu * phi_u - zl * phi_l) / dp)
        ez3 = 2.0 * ez1 + (-(zu**2 * phi_u - zl**2 * phi_l) / dp)

        return ez1, ez2, ez3

    def _compute_ema_moments(
        self, mean_log: float, std_log: float, skew: float
    ) -> Tuple[float, float, float]:
        """Compute expected moments given current parameter estimates.

        Uses analytical truncated distribution moments following the Fortran
        mP3 approach: incomplete gamma moments for |skew| > ~0.001, and
        truncated normal (Wilson-Hilferty) for near-zero skew.
        """
        sum_x = 0.0
        sum_x2 = 0.0
        sum_x3 = 0.0

        for interval in self._intervals:
            if not interval.is_censored:
                x = np.log10(interval.lower)
                sum_x += x
                sum_x2 += x**2
                sum_x3 += x**3
            else:
                lower = interval.lower if interval.lower > 0 else 1e-10
                upper = interval.upper if np.isfinite(interval.upper) else np.inf

                log_lower = np.log10(lower)
                log_upper = np.log10(upper) if np.isfinite(upper) else np.inf

                if abs(skew) < 0.001:
                    # Near-zero skew: use truncated normal moments
                    zl = (log_lower - mean_log) / std_log if std_log > 0 else -20.0
                    zu = (
                        (log_upper - mean_log) / std_log
                        if std_log > 0 and np.isfinite(log_upper)
                        else 20.0
                    )
                    zl = np.clip(zl, -20.0, 20.0)
                    zu = np.clip(zu, -20.0, 20.0)

                    ez1, ez2, ez3 = self._truncated_normal_moments(zl, zu)
                    # Transform back: x = mean + std*z
                    ex = mean_log + std_log * ez1
                    ex2 = mean_log**2 + 2 * mean_log * std_log * ez1 + std_log**2 * ez2
                    ex3 = (
                        mean_log**3
                        + 3 * mean_log**2 * std_log * ez1
                        + 3 * mean_log * std_log**2 * ez2
                        + std_log**3 * ez3
                    )
                else:
                    # Convert LP3 moments (mean, variance, skew) to gamma params
                    # m = (mu, sigma^2, gamma) -> (tau, alpha, beta)
                    alpha = 4.0 / skew**2
                    beta = std_log * abs(skew) / 2.0
                    if skew < 0:
                        beta = -beta
                    tau = mean_log - alpha * beta

                    # Standardize bounds for Gamma(alpha,1)
                    if beta > 0:
                        s_lower = max(0.0, (log_lower - tau) / beta)
                        s_upper = (log_upper - tau) / beta if np.isfinite(log_upper) else np.inf
                    else:
                        # Negative beta: flip and work with |beta|
                        abs_beta = -beta
                        s_lower = (
                            max(0.0, (tau - log_upper) / abs_beta)
                            if np.isfinite(log_upper)
                            else 0.0
                        )
                        s_upper = (tau - log_lower) / abs_beta

                    if s_upper <= s_lower:
                        # Interval outside distribution support
                        mid = (
                            log_lower + (log_upper if np.isfinite(log_upper) else log_lower)
                        ) / 2.0
                        sum_x += mid
                        sum_x2 += mid**2
                        sum_x3 += mid**3
                        continue

                    # Compute truncated gamma moments E[Y^k] for Y~Gamma(alpha,1)
                    gy1 = self._truncated_gamma_moment(alpha, s_lower, s_upper, 1)
                    gy2 = self._truncated_gamma_moment(alpha, s_lower, s_upper, 2)
                    gy3 = self._truncated_gamma_moment(alpha, s_lower, s_upper, 3)

                    # Transform back: x = tau + beta*Y (positive skew)
                    # or x = tau - |beta|*Y (negative skew)
                    if beta > 0:
                        ex = tau + beta * gy1
                        ex2 = tau**2 + 2 * tau * beta * gy1 + beta**2 * gy2
                        ex3 = (
                            tau**3
                            + 3 * tau**2 * beta * gy1
                            + 3 * tau * beta**2 * gy2
                            + beta**3 * gy3
                        )
                    else:
                        abs_beta = -beta
                        ex = tau - abs_beta * gy1
                        ex2 = tau**2 - 2 * tau * abs_beta * gy1 + abs_beta**2 * gy2
                        ex3 = (
                            tau**3
                            - 3 * tau**2 * abs_beta * gy1
                            + 3 * tau * abs_beta**2 * gy2
                            - abs_beta**3 * gy3
                        )

                sum_x += ex
                sum_x2 += ex2
                sum_x3 += ex3

        return sum_x, sum_x2, sum_x3

    def _ema_iteration(
        self, mean_log: float, std_log: float, skew: float
    ) -> Tuple[float, float, float]:
        """Perform one EMA iteration."""
        n = len(self._intervals)

        sum_x, sum_x2, sum_x3 = self._compute_ema_moments(mean_log, std_log, skew)

        new_mean = sum_x / n
        var = (sum_x2 - n * new_mean**2) / (n - 1)
        new_std = np.sqrt(max(var, 1e-10))

        if n > 2:
            m3 = sum_x3 - 3 * new_mean * sum_x2 + 2 * n * new_mean**3
            new_skew = (n * m3) / ((n - 1) * (n - 2) * new_std**3)
            new_skew = float(np.clip(new_skew, -MAX_ABS_SKEW, MAX_ABS_SKEW))
        else:
            new_skew = skew

        return new_mean, new_std, new_skew

    # Memoized because it is pure and expensive: a scipy.integrate.quad over the
    # beta-order-statistic integrand, ~85 ms a call, and MGBT is ~99% of the
    # wall time of a run_analysis() on a typical record.
    #
    # It buys nothing inside a single analysis. All 22 keys one run_analysis()
    # generates are distinct -- measured, not assumed -- so the hit rate there
    # is 0%. It pays only where the same fit repeats: 66 calls over 22 distinct
    # keys across three identical analyses, a 67% hit rate. That is the test
    # suite's access pattern, not a user's. Treat this as a suite-speed lever,
    # not a performance fix.
    #
    # Decorator order is load-bearing: staticmethod must be outermost. Before
    # Python 3.10 a staticmethod object is not callable, so lru_cache above
    # staticmethod would wrap something it cannot call and break collection on
    # the 3.9 matrix job.
    @staticmethod
    @lru_cache(maxsize=4096)
    def _mgbt_pvalue(n: int, r: int, w: float) -> float:
        """Exact MGBT p-value via GGBCRITP algorithm (peakfqr/probfun.f).

        Direct Python translation of the GGBCRITP/FGGB Fortran functions.
        Integrates over the distribution of the r-th normal order statistic
        using the non-central t-distribution to compute
        P(W(r) ≤ w | H0: data are iid standard normal).

        Reference: Cohn et al. (2013), probfun.f lines 1629–1750.

        Parameters
        ----------
        n : int
            Total sample size passed to MGBT (including zero-coded values).
        r : int
            1-based rank of the candidate outlier (r=1 = smallest).
        w : float
            Test statistic W(r) = (ZT(r) − mean(upper)) / sqrt(var(upper)).

        Returns
        -------
        float
            One-sided p-value in [0, 1].  Small values indicate a low outlier.
        """
        nc = n - r  # size of upper (conditioning) sample
        if nc < 2:
            return 1.0

        def _fggb(pzr: float) -> float:
            """Integrand: conditional p-value given the order-statistic draw."""
            if pzr <= 0.0 or pzr >= 1.0:
                return 0.0
            # r-th order statistic value via beta → normal transform
            pr = float(stats.beta.ppf(pzr, r, n + 1 - r))
            if pr <= 0.0 or pr >= 1.0:
                return 0.0
            zr = float(stats.norm.ppf(pr))
            if not np.isfinite(zr):
                return 0.0

            # Truncated-normal moments of the upper (nc) observations
            h = float(stats.norm.pdf(zr)) / max(1e-10, 1.0 - pr)
            ex1 = h
            ex2 = 1.0 + h * zr
            ex3 = 2.0 * ex1 + h * zr**2
            ex4 = 3.0 * ex2 + h * zr**3

            mu_m = ex1
            mu_s2 = ex2 - ex1**2
            if mu_s2 <= 0.0:
                return 0.0

            var_m = mu_s2 / nc
            var_s2 = (
                ex4 - 4 * ex3 * ex1 + 6 * ex2 * ex1**2 - 3 * ex1**4 - mu_s2**2
            ) / nc + 2.0 / ((nc - 1.0) * nc) * mu_s2**2
            if var_s2 <= 0.0:
                return 0.0

            # Gamma approximation to the distribution of S (upper std dev)
            alpha_g = mu_s2**2 / var_s2
            beta_g = mu_s2 / alpha_g
            cov_ms2 = (ex3 - 3 * ex2 * ex1 + 2 * ex1**3) / np.sqrt(nc * (nc - 1.0))
            mu_s = np.sqrt(beta_g) * np.exp(gammaln(alpha_g + 0.5) - gammaln(alpha_g))
            cov_ms = cov_ms2 / (2.0 * mu_s)
            var_s = mu_s2 - mu_s**2
            if var_s <= 0.0:
                return 0.0

            # Non-central t parameters (Fortran FGGB lines 1741-1748)
            lambda_ = cov_ms / var_s
            eta_p = w + lambda_
            mu_mp = mu_m - lambda_ * mu_s
            var_mp = var_m - cov_ms**2 / var_s
            if var_mp <= 0.0:
                return 0.0

            df = 2.0 * alpha_g
            ncp = (mu_mp - zr) / np.sqrt(var_mp)
            q = -np.sqrt(mu_s2 / var_mp) * eta_p

            # Match Fortran FP_TNC_CDF (probfun.f lines 1341-1350):
            # when NU > 20, use Abramowitz & Stegun p.949 normal approximation
            if df > 20.0:
                z_approx = (q * (1.0 - 1.0 / (4.0 * df)) - ncp) / np.sqrt(1.0 + q**2 / (2.0 * df))
                tnc_cdf = float(stats.norm.cdf(z_approx))
            else:
                tnc_cdf = float(stats.nct.cdf(q, df, ncp))
            return 1.0 - tnc_cdf

        try:
            result, _ = _quad(_fggb, 0.0, 1.0, limit=50, epsabs=1e-4, epsrel=1e-4)
            return float(np.clip(result, 0.0, 1.0))
        except Exception:
            return 1.0

    def _multiple_grubbs_beck(
        self, user_threshold: Optional[float] = None
    ) -> Tuple[float, int, List[float]]:
        """Perform Multiple Grubbs-Beck test for low outliers.

        Implements the three-sweep MGBT algorithm following Cohn et al. (2013)
        and the Fortran MGBTP routine in peakfqr/probfun.f.

        Test statistic for rank i (1-based, ascending):
            W(i) = (ZT(i) - mean(ZT[i+1..N])) / sqrt(var(ZT[i+1..N]))

        Three sweeps determine the number of low outliers:
          1. Outward sweep from median, α = 0.005 → J1
          2. Inward sweep from J1, α = 0.0   → J2 (= J1 in practice)
          3. Zeroin sweep from minimum, α = 0.1  → J3
        Result = MAX(J1, J2, J3).

        Parameters
        ----------
        user_threshold : float, optional
            User-supplied override.  When provided and > 0, skips the
            statistical test and uses this value directly.

        Returns
        -------
        tuple
            (threshold_cfs, n_low_outliers, pilf_flow_list)
        """
        # User override takes precedence
        if user_threshold is not None and user_threshold > 0:
            # Zeros are always below any positive threshold
            zero_pilf: List[float] = [0.0] * self._n_zeros
            nonzero_pilf = sorted([f for f in self._peak_flows if f < user_threshold])
            pilf = zero_pilf + nonzero_pilf
            return float(user_threshold), len(pilf), pilf

        # Include zero-flow years as log10(1e-88) = -88, matching Fortran MGBTP
        # which uses MAX(1D-88,X) before LOG10 (probfun.f line 1547).  Zeros must
        # be present in the test so that their extreme low values push the MGBT
        # sweep into the non-zero range and produce the correct threshold.
        log_zeros = np.full(self._n_zeros, np.log10(1e-88))
        log_flows_for_mgbt = np.concatenate([self.log_flows, log_zeros])
        n = len(log_flows_for_mgbt)

        # Need at least 5 observations for a meaningful test
        if n < 5:
            return 0.0, 0, []

        # Sort ascending (zt[0] = smallest, zt[n-1] = largest)
        zt = np.sort(log_flows_for_mgbt)

        # Median position (Fortran N/2, integer division)
        n2 = n // 2

        # ── Compute W(i) and p-values for i = N2 down to 1 (1-based) ──────
        # Initial accumulator: ZT(N2+2) to ZT(N) in Fortran 1-based
        # = zt[n2+1] to zt[n-1] in Python 0-based
        s1 = float(np.sum(zt[n2 + 1 :]))
        s2 = float(np.sum(zt[n2 + 1 :] ** 2))

        pvalues = np.ones(n2 + 1)  # 1-indexed; pvalues[0] unused

        for i_f in range(n2, 0, -1):  # Fortran I = N2 down to 1
            # Add ZT(I+1) = zt[i_f] (0-based) to the upper accumulator
            s1 += zt[i_f]
            s2 += zt[i_f] ** 2
            nc = n - i_f  # number of upper observations (I+1 to N)

            if nc < 2:
                pvalues[i_f] = 1.0
                continue

            xm = s1 / nc
            xv = (s2 - nc * xm**2) / (nc - 1)

            if xv <= 0.0:
                pvalues[i_f] = 1.0
                continue

            # ZT(I) = zt[i_f - 1] (0-based)
            w_i = (zt[i_f - 1] - xm) / np.sqrt(xv)
            # float() keeps the cache key a plain float rather than a
            # numpy scalar, so repeats of the same fit land on one entry.
            pvalues[i_f] = self._mgbt_pvalue(int(n), i_f, float(w_i))

        # ── Three sweeps ────────────────────────────────────────────────────
        alpha_out = 0.005  # outward sweep
        alpha_in = 0.0  # inward sweep  (always = J1 with Alphain=0)
        alpha_zero_in = 0.1  # zeroin sweep

        # Step 1: Outward sweep – largest i from N2 to 1 where p < alpha_out
        j1 = 0
        for i_f in range(n2, 0, -1):
            if pvalues[i_f] < alpha_out:
                j1 = i_f
                break

        # Step 2: Inward sweep – first i from j1+1 to N2 where p >= alpha_in
        j2 = j1
        for i_f in range(j1 + 1, n2 + 1):
            if pvalues[i_f] >= alpha_in:
                j2 = i_f - 1
                break

        # Step 3: Zeroin sweep – first i from 1 to N2 where p >= alpha_zero_in
        j3 = n2  # fallback: all below-median positions would be outliers
        for i_f in range(1, n2 + 1):
            if pvalues[i_f] >= alpha_zero_in:
                j3 = i_f - 1
                break

        n_low_outliers = max(j1, j2, j3)

        if n_low_outliers == 0:
            return 0.0, 0, []

        # Threshold = flow of the first NON-outlier (Fortran: qs(gbnlow+1))
        # zt is 0-based; the n_low_outliers smallest are outliers
        threshold = 10.0 ** zt[n_low_outliers]
        # Zeros are coded as log10(1e-88) = -88 in zt; restore them as 0.0
        pilf = [0.0 if zt[k] < -80.0 else 10.0 ** zt[k] for k in range(n_low_outliers)]

        return threshold, n_low_outliers, pilf

    def run_analysis(self) -> FrequencyResults:
        """Run EMA flood frequency analysis."""
        log_flows = self.log_flows
        n = self.n

        mean_log = np.mean(log_flows)
        std_log = np.std(log_flows, ddof=1)
        skew_station = n * np.sum((log_flows - mean_log) ** 3) / ((n - 1) * (n - 2) * std_log**3)
        skew_station = float(np.clip(skew_station, -MAX_ABS_SKEW, MAX_ABS_SKEW))

        low_threshold, n_low_outliers, pilf_flows = self._multiple_grubbs_beck(
            user_threshold=self._user_low_outlier_threshold
        )
        self._ema_params.low_outlier_threshold = low_threshold

        self._build_flow_intervals(low_threshold)

        converged = False
        iteration = 0

        for iteration in range(self._ema_params.max_iterations):
            new_mean, new_std, new_skew = self._ema_iteration(mean_log, std_log, skew_station)

            if (
                abs(new_mean - mean_log) < self._ema_params.tolerance
                and abs(new_std - std_log) < self._ema_params.tolerance
                and abs(new_skew - skew_station) < self._ema_params.tolerance
            ):
                converged = True
                mean_log = new_mean
                std_log = new_std
                skew_station = new_skew
                break

            mean_log = new_mean
            std_log = new_std
            skew_station = new_skew

        n_systematic = sum(1 for i in self._intervals if i.is_systematic and not i.is_censored)
        n_historical = sum(1 for i in self._intervals if i.is_historical)
        n_censored = sum(1 for i in self._intervals if i.is_censored)

        # For the weighted-skew *point estimate*, use the total number of
        # intervals as the effective sample size. This accounts for information
        # from both observed peaks and censored intervals and best matches
        # PeakfqSA's weighted skew value.
        n_intervals = len(self._intervals)
        skew_weighted, _ = self._compute_skew_weighting(skew_station, n_effective=n_intervals)
        skew_used = skew_weighted if skew_weighted is not None else skew_station

        # For the skew's *variance* (used to widen confidence intervals), use
        # the count of actual observations (systematic + historical) instead.
        # Using n_intervals here would understate skew uncertainty, since
        # censored intervals contribute much less information about skew than
        # true observations do.
        n_observed = n_systematic + n_historical
        _, skew_used_mse = self._compute_skew_weighting(skew_station, n_effective=n_observed)

        self._results = FrequencyResults(
            n_peaks=len(self._intervals),
            n_systematic=n_systematic,
            n_historical=n_historical,
            n_censored=n_censored,
            n_low_outliers=n_low_outliers,
            mean_log=mean_log,
            std_log=std_log,
            skew_station=skew_station,
            skew_regional=self._regional_skew,
            skew_weighted=skew_weighted,
            skew_used=skew_used,
            skew_used_mse=skew_used_mse,
            low_outlier_threshold=low_threshold,
            mgb_critical_value=grubbs_beck_critical_value(n),
            method=AnalysisMethod.EMA,
            ema_iterations=iteration + 1,
            ema_converged=converged,
            n_zeros=self._n_zeros,
            pilf_flows=pilf_flows,
        )

        self._results.quantiles = self.compute_quantiles()
        self._results.confidence_limits = self.compute_confidence_limits()

        return self._results

    @property
    def intervals(self) -> List[FlowInterval]:
        return self._intervals.copy()

    def get_perception_thresholds_table(self) -> pd.DataFrame:
        """
        Generate a perception thresholds table similar to HEC-SSP.

        Returns a DataFrame with columns:
        - Start Year: Beginning of the perception period
        - End Year: End of the perception period
        - Low Threshold (cfs): Minimum detectable flow (0 = no lower limit)
        - High Threshold (cfs): Maximum detectable flow (inf = no upper limit)
        - Comments: Description of the record type

        Returns
        -------
        pd.DataFrame
            Perception thresholds table in HEC-SSP format
        """
        if not self._intervals:
            self._build_flow_intervals(self._ema_params.low_outlier_threshold or 0.0)

        # Build perception periods from intervals and parameters
        periods = []

        # Get systematic record period
        sys_start = self._ema_params.systematic_start
        sys_end = self._ema_params.systematic_end
        low_outlier_thresh = self._ema_params.low_outlier_threshold or 0.0

        # Historical period (if exists)
        if self._ema_params.historical_start and self._ema_params.historical_threshold:
            hist_start = self._ema_params.historical_start
            hist_end = self._ema_params.historical_end or (sys_start - 1)
            hist_threshold = self._ema_params.historical_threshold

            periods.append(
                {
                    "Start Year": hist_start,
                    "End Year": hist_end,
                    "Low Threshold (cfs)": hist_threshold,
                    "High Threshold (cfs)": np.inf,
                    "Comments": "Historical Record",
                }
            )

        # Systematic record - check for gaps and different thresholds
        # Group consecutive years with same perception thresholds
        systematic_intervals = [
            i for i in self._intervals if i.is_systematic and sys_start <= i.year <= sys_end
        ]

        if systematic_intervals:
            # Sort by year
            systematic_intervals.sort(key=lambda x: x.year)

            # For systematic record, perception is typically 0 to infinity
            # unless there are low outlier thresholds
            periods.append(
                {
                    "Start Year": sys_start,
                    "End Year": sys_end,
                    "Low Threshold (cfs)": 0.0,
                    "High Threshold (cfs)": np.inf,
                    "Comments": "Systematic Record",
                }
            )

        # Check for explicitly defined perception thresholds
        for (start, end), threshold in self._perception_thresholds.items():
            periods.append(
                {
                    "Start Year": start,
                    "End Year": end,
                    "Low Threshold (cfs)": threshold,
                    "High Threshold (cfs)": np.inf,
                    "Comments": "User-Defined Threshold",
                }
            )

        # Add low outlier censoring info if applicable
        if low_outlier_thresh > 0 and self._results and self._results.n_low_outliers > 0:
            # Find years with low outliers
            low_outlier_years = [
                i.year for i in self._intervals if i.is_censored and i.upper <= low_outlier_thresh
            ]
            if low_outlier_years:
                periods.append(
                    {
                        "Start Year": min(low_outlier_years),
                        "End Year": max(low_outlier_years),
                        "Low Threshold (cfs)": low_outlier_thresh,
                        "High Threshold (cfs)": np.inf,
                        "Comments": f"Low Outlier Censoring ({len(low_outlier_years)} years)",
                    }
                )

        # Sort periods by start year
        periods.sort(key=lambda x: x["Start Year"])

        # Create DataFrame
        df = pd.DataFrame(periods)

        # Format High Threshold for display
        if not df.empty:
            df["High Threshold (cfs)"] = df["High Threshold (cfs)"].apply(
                lambda x: "Infinity" if np.isinf(x) else f"{x:,.0f}"
            )
            df["Low Threshold (cfs)"] = df["Low Threshold (cfs)"].apply(
                lambda x: f"{x:,.0f}" if x > 0 else "0"
            )

        return df


class Bulletin17C:
    """
    Unified interface for Bulletin 17C flood frequency analysis.

    Allows selection between Method of Moments (MOM) and Expected
    Moments Algorithm (EMA) methods.

    Examples
    --------
    >>> b17c = Bulletin17C(peak_flows)
    >>> results = b17c.run_analysis(method='mom')

    >>> b17c = Bulletin17C(peak_flows, regional_skew=-0.05, regional_skew_mse=0.12)
    >>> results = b17c.run_analysis(method='ema')
    """

    def __init__(
        self,
        peak_flows: np.ndarray,
        water_years: np.ndarray = None,
        regional_skew: float = None,
        regional_skew_mse: float = None,
        historical_peaks: List[Tuple[int, float]] = None,
        perception_thresholds: Dict[Tuple[int, int], float] = None,
        ema_params: EMAParameters = None,
        user_low_outlier_threshold: Optional[float] = None,
    ):
        self._peak_flows = np.array(peak_flows)
        self._water_years = water_years
        self._regional_skew = regional_skew
        self._regional_skew_mse = regional_skew_mse
        self._historical_peaks = historical_peaks
        self._perception_thresholds = perception_thresholds
        self._ema_params = ema_params
        self._user_low_outlier_threshold = user_low_outlier_threshold

        self._analyzer: Optional[FloodFrequencyAnalysis] = None
        self._results: Optional[FrequencyResults] = None

    @property
    def results(self) -> Optional[FrequencyResults]:
        return self._results

    @property
    def quantiles(self) -> Optional[pd.DataFrame]:
        return self._results.quantiles if self._results else None

    @property
    def confidence_limits(self) -> Optional[pd.DataFrame]:
        return self._results.confidence_limits if self._results else None

    @property
    def n(self) -> int:
        return len(self._peak_flows)

    @property
    def mean_log(self) -> Optional[float]:
        return self._results.mean_log if self._results else None

    @property
    def std_log(self) -> Optional[float]:
        return self._results.std_log if self._results else None

    @property
    def skew_station(self) -> Optional[float]:
        return self._results.skew_station if self._results else None

    @property
    def skew_weighted(self) -> Optional[float]:
        return self._results.skew_weighted if self._results else None

    @property
    def skew_used(self) -> Optional[float]:
        return self._results.skew_used if self._results else None

    @property
    def low_outlier_threshold(self) -> Optional[float]:
        return self._results.low_outlier_threshold if self._results else None

    @property
    def n_low_outliers(self) -> Optional[int]:
        return self._results.n_low_outliers if self._results else None

    def run_analysis(self, method: Union[str, AnalysisMethod] = "ema") -> FrequencyResults:
        """
        Run flood frequency analysis.

        Parameters
        ----------
        method : str or AnalysisMethod
            'mom' for Method of Moments, 'ema' for Expected Moments Algorithm
        """
        if isinstance(method, str):
            method = AnalysisMethod[method.upper()]

        if method == AnalysisMethod.MOM:
            self._analyzer = MethodOfMoments(
                self._peak_flows,
                regional_skew=self._regional_skew,
                regional_skew_mse=self._regional_skew_mse,
            )
        else:
            self._analyzer = ExpectedMomentsAlgorithm(
                self._peak_flows,
                water_years=self._water_years,
                regional_skew=self._regional_skew,
                regional_skew_mse=self._regional_skew_mse,
                ema_params=self._ema_params,
                historical_peaks=self._historical_peaks,
                perception_thresholds=self._perception_thresholds,
                user_low_outlier_threshold=self._user_low_outlier_threshold,
            )

        self._results = self._analyzer.run_analysis()
        return self._results

    def compute_quantiles(self, aep: np.ndarray = None) -> pd.DataFrame:
        if self._analyzer is None:
            self.run_analysis()
        return self._analyzer.compute_quantiles(aep)

    def compute_confidence_limits(
        self, aep: np.ndarray = None, confidence: float = 0.90
    ) -> pd.DataFrame:
        if self._analyzer is None:
            self.run_analysis()
        return self._analyzer.compute_confidence_limits(aep, confidence)

    def plot_frequency_curve(
        self,
        site_name: str = None,
        site_no: str = None,
        save_path: str = None,
        figsize: Tuple[int, int] = (10, 7),
        show_confidence: bool = True,
    ) -> plt.Figure:
        if self._analyzer is None:
            self.run_analysis()
        return self._analyzer.plot_frequency_curve(
            site_name, site_no, save_path, figsize, show_confidence
        )

    def get_perception_thresholds_table(self) -> pd.DataFrame:
        """
        Generate a perception thresholds table similar to HEC-SSP.

        This table shows the flow ranges that could be detected/recorded
        during different time periods, following the HEC-SSP format.

        Returns a DataFrame with columns:
        - Start Year: Beginning of the perception period
        - End Year: End of the perception period
        - Low Threshold (cfs): Minimum detectable flow (0 = no lower limit)
        - High Threshold (cfs): Maximum detectable flow (Infinity = no upper limit)
        - Comments: Description of the record type (e.g., "Systematic Record",
          "Historical Record", "Low Outlier Censoring")

        Returns
        -------
        pd.DataFrame
            Perception thresholds table in HEC-SSP format

        Raises
        ------
        ValueError
            If analysis method is MOM (perception thresholds only apply to EMA)

        Examples
        --------
        >>> b17c = Bulletin17C(peak_flows, water_years=years)
        >>> b17c.run_analysis(method='ema')
        >>> table = b17c.get_perception_thresholds_table()
        >>> print(table)
           Start Year  End Year Low Threshold (cfs) High Threshold (cfs)          Comments
        0        1900      1930               5000             Infinity  Historical Record
        1        1931      2020                  0             Infinity  Systematic Record
        """
        if self._analyzer is None:
            self.run_analysis()

        if not isinstance(self._analyzer, ExpectedMomentsAlgorithm):
            raise ValueError(
                "Perception thresholds table is only available for EMA analysis. "
                "Run analysis with method='ema' first."
            )

        return self._analyzer.get_perception_thresholds_table()

    def to_comparison_dict(self) -> dict:
        """Convert native results to a dict compatible with FrequencyComparator.

        Returns
        -------
        dict
            Keys: 'parameters', 'quantiles', 'confidence_intervals'.

        Raises
        ------
        ValueError
            If analysis has not been run yet.
        """
        if self._results is None:
            raise ValueError("Run analysis before converting to comparison dict")

        r = self._results

        parameters = {
            "mean_log": r.mean_log,
            "std_log": r.std_log,
            "skew_at_site": r.skew_station,
        }
        if r.skew_weighted is not None:
            parameters["skew_weighted"] = r.skew_weighted

        # Extract quantiles from DataFrame
        quantiles: dict[float, float] = {}
        if r.quantiles is not None and not r.quantiles.empty:
            for _, row in r.quantiles.iterrows():
                quantiles[float(row["aep"])] = float(row["flow_cfs"])

        # Extract confidence intervals from DataFrame
        confidence_intervals: dict[float, tuple[float, float]] = {}
        if r.confidence_limits is not None and not r.confidence_limits.empty:
            for _, row in r.confidence_limits.iterrows():
                confidence_intervals[float(row["aep"])] = (
                    float(row["lower_5pct"]),
                    float(row["upper_5pct"]),
                )

        return {
            "parameters": parameters,
            "quantiles": quantiles,
            "confidence_intervals": confidence_intervals,
        }

    def validate(
        self,
        reference: "PeakfqSAResult",
        tolerance_pct: float = 1.0,
        parameter_tolerance_pct: float = 0.5,
        ci_tolerance_pct: float = 2.0,
    ) -> "ComparisonResult":
        """Validate native results against a PeakfqSA reference.

        Parameters
        ----------
        reference : PeakfqSAResult
            Reference result from PeakfqSA.
        tolerance_pct : float
            Tolerance for quantile comparisons.
        parameter_tolerance_pct : float
            Tolerance for LP3 parameter comparisons.
        ci_tolerance_pct : float
            Tolerance for confidence interval comparisons.

        Returns
        -------
        ComparisonResult
            Detailed comparison result.
        """
        from hydrolib.validation.comparisons import FrequencyComparator

        if self._results is None:
            self.run_analysis()

        native_dict = self.to_comparison_dict()
        comparator = FrequencyComparator(
            tolerance_pct=tolerance_pct,
            parameter_tolerance_pct=parameter_tolerance_pct,
            ci_tolerance_pct=ci_tolerance_pct,
        )
        return comparator.compare(native_dict, reference)
