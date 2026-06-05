"""
hydrolib.bulletin17c - Bulletin 17C Flood Frequency Analysis

Implements both Method of Moments (MOM) and Expected Moments Algorithm (EMA).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from datetime import datetime
from functools import cached_property
from typing import ClassVar, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import gammainc, gammaincc, gammaln, ndtri

from .core import (
    AnalysisMethod,
    EMAParameters,
    FlowInterval,
    FrequencyResults,
    grubbs_beck_critical_value,
    kfactor,
    kfactor_array,
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

    def _compute_weighted_skew(self, station_skew: float) -> Optional[float]:
        """Compute weighted skew from station and regional values."""
        if self._regional_skew is None or self._regional_skew_mse is None:
            return None

        n = self.n
        mse_station = (6 / n) * (1 + (6 / n) * station_skew**2 + (15 / (n**2)) * station_skew**4)

        w_regional = mse_station / (mse_station + self._regional_skew_mse)
        w_station = self._regional_skew_mse / (mse_station + self._regional_skew_mse)

        return w_station * station_skew + w_regional * self._regional_skew

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

        var_factor = 1 / n + K**2 * (1 + 0.75 * G**2) / (2 * (n - 1))
        se_log = self._results.std_log * np.sqrt(var_factor)

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

        ax.set_yscale("log")
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

        skew_weighted = self._compute_weighted_skew(skew_station)
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
    ):
        super().__init__(peak_flows, regional_skew, regional_skew_mse)

        if water_years is not None:
            self._water_years = np.array(water_years)
        else:
            end_year = datetime.now().year
            self._water_years = np.arange(end_year - len(peak_flows) + 1, end_year + 1)

        self._historical_peaks = historical_peaks or []
        self._perception_thresholds = perception_thresholds or {}
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

        if self._historical_peaks:
            hist_years = [h[0] for h in self._historical_peaks]
            hist_start = min(hist_years)
            hist_end = max(max(hist_years), sys_start - 1)
            hist_threshold = max(h[1] for h in self._historical_peaks)
        elif gaps:
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
        else:
            new_skew = skew

        return new_mean, new_std, new_skew

    def _multiple_grubbs_beck(self) -> Tuple[float, int]:
        """Perform Multiple Grubbs-Beck test for low outliers."""
        log_flows = self.log_flows
        current_flows = np.sort(log_flows)
        n_outliers = 0
        threshold_log = -np.inf

        while len(current_flows) > 10:
            mean_log = np.mean(current_flows)
            std_log = np.std(current_flows, ddof=1)

            k_n = grubbs_beck_critical_value(len(current_flows))
            test_threshold = mean_log - k_n * std_log

            if current_flows[0] < test_threshold:
                threshold_log = test_threshold
                n_outliers += 1
                current_flows = current_flows[1:]
            else:
                break

        threshold = 10**threshold_log if threshold_log > -np.inf else 0.0
        return threshold, n_outliers

    def run_analysis(self) -> FrequencyResults:
        """Run EMA flood frequency analysis."""
        log_flows = self.log_flows
        n = self.n

        mean_log = np.mean(log_flows)
        std_log = np.std(log_flows, ddof=1)
        skew_station = n * np.sum((log_flows - mean_log) ** 3) / ((n - 1) * (n - 2) * std_log**3)

        low_threshold, n_low_outliers = self._multiple_grubbs_beck()
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

        skew_weighted = self._compute_weighted_skew(skew_station)
        skew_used = skew_weighted if skew_weighted is not None else skew_station

        n_systematic = sum(1 for i in self._intervals if i.is_systematic and not i.is_censored)
        n_historical = sum(1 for i in self._intervals if i.is_historical)
        n_censored = sum(1 for i in self._intervals if i.is_censored)

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
            low_outlier_threshold=low_threshold,
            mgb_critical_value=grubbs_beck_critical_value(n),
            method=AnalysisMethod.EMA,
            ema_iterations=iteration + 1,
            ema_converged=converged,
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
    ):
        self._peak_flows = np.array(peak_flows)
        self._water_years = water_years
        self._regional_skew = regional_skew
        self._regional_skew_mse = regional_skew_mse
        self._historical_peaks = historical_peaks
        self._perception_thresholds = perception_thresholds
        self._ema_params = ema_params

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
