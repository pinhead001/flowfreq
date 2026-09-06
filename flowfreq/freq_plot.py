"""
flowfreq.freq_plot - Frequency curve plotting.

Returns matplotlib ``Figure`` objects. Nothing here imports Streamlit; the
module and its main function were named for the app that first consumed them,
which now lives in its own repository.
"""

from __future__ import annotations

import logging
import math
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm

from .core import kfactor_array, log_pearson3_cdf

if TYPE_CHECKING:
    from .bulletin17c import Bulletin17C

logger = logging.getLogger(__name__)

_FONT_SIZE = 12
_ANNOT_FONT_SIZE = 9
_BOX_PADDING = 0.03

# Per-skew-option display style: color, linestyle
_SKEW_STYLE: Dict[str, Tuple[str, str]] = {
    "Station Skew": ("steelblue", "-"),
    "Weighted Skew": ("k", "-"),
    "Regional Skew": ("darkorange", "-"),
}
# Fallback palette when label not in _SKEW_STYLE
_FALLBACK_COLORS = ["steelblue", "k", "darkorange", "forestgreen"]
# Style for overlay curves fitted from a different moment set: dashed purple,
# deliberately unlike anything in _SKEW_STYLE so the two cannot be confused.
_EXTRA_CURVE_COLOR = "mediumorchid"


def _lp3_quantiles(
    mean_log: float,
    std_log: float,
    skew: float,
    aep: np.ndarray,
) -> np.ndarray:
    """Return LP3 flow estimates (cfs) for the given AEP array."""
    K = kfactor_array(skew, aep)
    return 10.0 ** (mean_log + K * std_log)


def _lp3_ci(
    mean_log: float,
    std_log: float,
    skew: float,
    n: int,
    aep: np.ndarray,
    z_alpha: float = 1.6449,  # 90% two-sided (norm.ppf(0.95))
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (lower, upper) 90% confidence bounds for LP3 quantiles."""
    K = kfactor_array(skew, aep)
    log_Q = mean_log + K * std_log
    var_factor = 1 / n + K**2 * (1 + 0.75 * skew**2) / (2 * (n - 1))
    se_log = std_log * np.sqrt(var_factor)
    return 10.0 ** (log_Q - z_alpha * se_log), 10.0 ** (log_Q + z_alpha * se_log)


def plot_frequency_curve(
    b17c: "Bulletin17C",
    site_name: Optional[str] = None,
    site_no: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 6),
    skew_curves: Optional[Dict[str, float]] = None,
    yscale: str = "log",
    max_flow_label: Optional[Dict] = None,
    extra_curves: Optional[Dict[str, Tuple[float, float, float, int]]] = None,
) -> plt.Figure:
    """Plot a Bulletin 17C frequency curve suitable for Streamlit display.

    Parameters
    ----------
    b17c : Bulletin17C
        A fitted Bulletin17C instance.
    site_name : str, optional
        Station name for the title.
    site_no : str, optional
        USGS station number for the title.
    figsize : tuple
        Figure size in inches.
    skew_curves : dict[str, float], optional
        Mapping of ``{label: skew_value}`` for each LP3 curve to draw.
        When None or a single entry, one black curve (like the summary
        hydrograph median line) is drawn.  With multiple entries each curve
        gets a distinct color with its own 90% CI band.
        Example: ``{"Station Skew": -0.15, "Weighted Skew": -0.22}``.
    extra_curves : dict[str, tuple[float, float, float, int]], optional
        Additional LP3 parameter sets to overlay, each as
        ``{label: (mean_log, std_log, skew, n)}``. Unlike ``skew_curves``,
        which varies only the skew of the fit already on the axes, these carry
        their own moments -- which is what a second *fit* needs. Use it to show
        the effect of a refit: a PILF-threshold override, or a
        perception-threshold-aware EMA beside the base analysis. Drawn in
        dashed purple with its own 90% band.

    Returns
    -------
    matplotlib.figure.Figure
        The frequency curve figure.
    """
    if b17c.results is None:
        b17c.run_analysis()

    r = b17c.results
    mean_log = r.mean_log
    std_log = r.std_log
    n_sys = r.n_systematic or r.n_peaks

    # Default: single curve using the skew already selected by EMA/MOM
    if not skew_curves:
        skew_curves = {"Fitted Frequency Curve": r.skew_used}

    multi = len(skew_curves) > 1

    def aep_to_x(p: float) -> float:
        return norm.ppf(1 - p)

    fig, ax = plt.subplots(figsize=figsize)

    # --- Observed peaks (Weibull plotting positions) ---
    # Separate PILF (below low-outlier threshold) from normal peaks for markers
    peak_flows = b17c._peak_flows  # noqa: SLF001
    pilf_set = set(r.pilf_flows) if r.pilf_flows else set()
    pilf_threshold = r.low_outlier_threshold if r.n_low_outliers > 0 else None

    sorted_flows = np.sort(peak_flows)[::-1]
    n_obs = len(sorted_flows)
    weibull_aep = np.arange(1, n_obs + 1) / (n_obs + 1)
    x_obs = [aep_to_x(p) for p in weibull_aep]

    # Normal peaks (above PILF threshold) → filled circles
    normal_mask = np.array([pilf_threshold is None or f >= pilf_threshold for f in sorted_flows])
    if np.any(normal_mask):
        ax.scatter(
            [x for x, m in zip(x_obs, normal_mask) if m],
            sorted_flows[normal_mask],
            c="steelblue",
            s=20,
            alpha=0.5,
            zorder=5,
            label="Observed Annual Peaks",
            edgecolors="navy",
            linewidth=0.5,
        )
    # PILF peaks → red × markers
    if np.any(~normal_mask):
        ax.scatter(
            [x for x, m in zip(x_obs, normal_mask) if not m],
            sorted_flows[~normal_mask],
            c="red",
            s=40,
            alpha=0.8,
            zorder=6,
            label="PILF (below MGBT threshold)",
            marker="x",
            linewidths=1.2,
        )

    # --- Max flow label annotation ---
    if max_flow_label:
        max_flow = max_flow_label.get("flow")
        max_year = max_flow_label.get("year")
        max_ri = max_flow_label.get("ri")
        if max_flow is not None and max_year is not None:
            # Convert to Python native types
            max_flow = float(max_flow)
            max_year = int(max_year)
            # Find the x position for the max flow (it's at rank 1, so AEP = 1/(n+1))
            max_aep = 1 / (n_obs + 1)
            max_x = aep_to_x(max_aep)
            # Build label text
            if max_ri is not None:
                max_ri = float(max_ri)
                ri_str = f"{max_ri:,.0f}" if max_ri >= 10 else f"{max_ri:.1f}"
                label_text = f"{max_flow:,.0f} cfs ({max_year})\n≈ {ri_str}-yr"
            else:
                label_text = f"{max_flow:,.0f} cfs ({max_year})"
            ax.annotate(
                label_text,
                xy=(max_x, max_flow),
                xytext=(-15, 0),
                textcoords="offset points",
                fontsize=_ANNOT_FONT_SIZE,
                ha="right",
                va="center",
                bbox=dict(
                    boxstyle="round,pad=0.2", facecolor="white", edgecolor="lightgray", alpha=0.9
                ),
                arrowprops=dict(arrowstyle="-", color="gray", alpha=0.7),
            )

    # --- LP3 curves (one per skew option) ---
    aep_fine = np.linspace(0.001, 0.999, 300)
    aep_ci = np.array(
        [0.999, 0.995, 0.99, 0.95, 0.90, 0.80, 0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002]
    )
    x_curve = [aep_to_x(p) for p in aep_fine]
    x_cl = [aep_to_x(p) for p in aep_ci]

    fallback_idx = 0
    for label, skew_val in skew_curves.items():
        # Map bare skew-option key → "Fitted Frequency Curve (Skew Type)" label
        if label in _SKEW_STYLE:
            color, ls = _SKEW_STYLE[label]
            curve_label = f"Fitted Frequency Curve ({label})"
        elif multi:
            color = _FALLBACK_COLORS[fallback_idx % len(_FALLBACK_COLORS)]
            ls = "-"
            fallback_idx += 1
            curve_label = label
        else:
            color, ls = "k", "-"
            curve_label = label  # keep whatever was passed for single-curve

        # Single curve always black
        if not multi:
            color, ls = "k", "-"

        Q_curve = _lp3_quantiles(mean_log, std_log, skew_val, aep_fine)
        ax.plot(
            x_curve, Q_curve, color=color, linestyle=ls, linewidth=1, label=curve_label, zorder=4
        )

        lower, upper = _lp3_ci(mean_log, std_log, skew_val, n_sys, aep_ci)
        ax.fill_between(x_cl, lower, upper, alpha=0.12, color=color)
        ax.plot(x_cl, lower, color=color, linestyle="--", linewidth=0.7, alpha=0.55)
        ax.plot(x_cl, upper, color=color, linestyle="--", linewidth=0.7, alpha=0.55)

    # --- Overlay curves from an independent fit (own mean/std, not just skew) ---
    for label, (ec_mean, ec_std, ec_skew, ec_n) in (extra_curves or {}).items():
        ax.plot(
            x_curve,
            _lp3_quantiles(ec_mean, ec_std, ec_skew, aep_fine),
            color=_EXTRA_CURVE_COLOR,
            linestyle="--",
            linewidth=1.5,
            label=label,
            zorder=4,
        )
        lower_ec, upper_ec = _lp3_ci(ec_mean, ec_std, ec_skew, ec_n, aep_ci)
        ax.fill_between(x_cl, lower_ec, upper_ec, alpha=0.10, color=_EXTRA_CURVE_COLOR)
        ax.plot(x_cl, lower_ec, color=_EXTRA_CURVE_COLOR, linestyle=":", linewidth=0.7, alpha=0.55)
        ax.plot(x_cl, upper_ec, color=_EXTRA_CURVE_COLOR, linestyle=":", linewidth=0.7, alpha=0.55)

    # --- Low outlier threshold ---
    if r.low_outlier_threshold > 0 and r.n_low_outliers > 0:
        ax.axhline(
            r.low_outlier_threshold,
            color="red",
            linestyle="--",
            alpha=0.7,
            label=f"MGBT Threshold ({r.n_low_outliers} PILF)",
        )

    # --- Axes ---
    ax.set_yscale(yscale)
    ax.set_ylabel("Discharge (cfs)", fontsize=_FONT_SIZE)
    ax.set_xlabel("Annual Exceedance Probability", fontsize=_FONT_SIZE)

    # Y-axis ticks: powers of 10, comma-formatted (matches hydrograph plots)
    min_flow = sorted_flows[-1] if n_obs else 1.0
    max_flow = sorted_flows[0] if n_obs else 1e6
    if yscale == "log" and min_flow > 0:
        min_exp = math.floor(math.log10(min_flow))
        max_exp = math.ceil(math.log10(max_flow))
        yticks = [10**i for i in range(min_exp, max_exp + 1)]
        ax.set_yticks(yticks)
        ax.set_yticklabels([f"{int(t):,}" for t in yticks])

    # X-axis probability ticks
    prob_ticks = [0.999, 0.99, 0.95, 0.90, 0.80, 0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002]
    x_ticks = [aep_to_x(p) for p in prob_ticks]
    ax.set_xticks(x_ticks)
    tick_labels = []
    for p in prob_ticks:
        pct = p * 100
        tick_labels.append(f"{pct:g}%" if pct >= 1 else f"{pct:.1f}%")
    ax.set_xticklabels(tick_labels, fontsize=_ANNOT_FONT_SIZE, rotation=45, ha="right")
    ax.set_xlim(aep_to_x(0.999), aep_to_x(0.002))

    # --- Secondary axis: return periods ---
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    rp_ticks = [1.5, 2, 5, 10, 25, 50, 100, 200, 500]
    rp_x = [aep_to_x(1 / rp) for rp in rp_ticks]
    ax2.set_xticks(rp_x)
    ax2.set_xticklabels(
        [str(int(rp)) if rp == int(rp) else str(rp) for rp in rp_ticks],
        fontsize=_ANNOT_FONT_SIZE,
    )
    ax2.set_xlabel("Return Period (years)", fontsize=_FONT_SIZE)

    # --- Title ---
    if site_name and site_no:
        title = f"Flood Frequency Curve\nUSGS {site_no} - {site_name}"
    elif site_no:
        title = f"Flood Frequency Curve - USGS {site_no}"
    else:
        title = "Flood Frequency Curve (Bulletin 17C)"
    ax.set_title(title, fontsize=_FONT_SIZE, fontweight="bold", pad=35)

    # --- LP3 parameters annotation (top-left, matches hydrograph style) ---
    skew_lines = "\n".join(
        f"\u03b3 {lbl.replace(' Skew', '').lower():<8} = {val:.3f}"
        for lbl, val in skew_curves.items()
    )
    stats_text = (
        f"n = {n_obs}\n"
        f"\u03bc(log Q) = {mean_log:.4f}\n"
        f"\u03c3(log Q) = {std_log:.4f}\n"
        f"{skew_lines}"
    )
    if r.n_low_outliers > 0:
        stats_text += f"\n{r.n_low_outliers} peak(s) below PILF threshold"
    if r.n_zeros > 0:
        stats_text += f"\n{r.n_zeros} zero flow(s) not displayed"
    ax.annotate(
        stats_text,
        xy=(_BOX_PADDING, 1 - _BOX_PADDING),
        xycoords="axes fraction",
        fontsize=_ANNOT_FONT_SIZE,
        ha="left",
        va="top",
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="lightgray", alpha=1.0),
    )

    ax.legend(loc="lower left", fontsize=_ANNOT_FONT_SIZE)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()

    return fig


_THRESHOLD_COLORS = ["darkorange", "forestgreen", "purple", "crimson", "teal"]


def plot_peak_flows_with_thresholds(
    peak_df: pd.DataFrame,
    site_name: Optional[str] = None,
    site_no: Optional[str] = None,
    thresholds: Optional[List[dict]] = None,
    figsize: Tuple[int, int] = (12, 5),
    ffa_year_range: Optional[Tuple[int, int]] = None,
    mgbt_threshold: Optional[float] = None,
    lp3_params: Optional[Tuple[float, float, float]] = None,
    return_periods: Tuple[float, ...] = (2, 5, 10, 25, 50, 100),
    annotate_max_peak: bool = True,
) -> plt.Figure:
    """Plot annual peak flows as a bar chart with optional perception threshold lines.

    Parameters
    ----------
    peak_df : pd.DataFrame
        DataFrame with columns ``water_year`` and ``peak_flow_cfs``.
    site_name : str, optional
        Station name for the title.
    site_no : str, optional
        USGS station number for the title.
    thresholds : list of dict, optional
        Each dict has keys ``start_year`` (int), ``end_year`` (int),
        ``lower_cfs`` (float, default 0), ``upper_cfs`` (float, default inf).
        Lower bound drawn as dashed line; upper bound drawn as dashed line or
        labelled at top of axis when upper_cfs is inf.
    figsize : tuple
        Figure size in inches.
    ffa_year_range : tuple of (int, int), optional
        ``(start_year, end_year)`` of the peak-flow record used in the FFA.
        Years outside this range are drawn as hollow outline bars so the user
        can see which peaks were excluded from the analysis.
    mgbt_threshold : float, optional
        MGBT low outlier threshold.  When provided, drawn as a solid red
        horizontal line across the full record.
    lp3_params : tuple of (mean_log, std_log, skew), optional
        The fitted LP3 moments, e.g. ``(r.mean_log, r.std_log, r.skew_used)``
        from a :class:`~flowfreq.core.FrequencyResults`. When given, draws a
        reference line at the flow for each entry in ``return_periods`` and,
        when ``annotate_max_peak`` is True, labels the largest peak in
        ``peak_df`` with its approximate recurrence interval under this fit.
        These are the two features the app's ``plot_peak_timeseries`` carried
        that this function previously did not (TODO.md).
    return_periods : tuple of float, optional
        Return periods (years) to draw reference lines for, used only when
        ``lp3_params`` is given. Default ``(2, 5, 10, 25, 50, 100)``.
    annotate_max_peak : bool, optional
        When True (default) and ``lp3_params`` is given, annotate the largest
        peak in the record with its estimated recurrence interval, read off
        the LP3 CDF at that flow.

    Returns
    -------
    matplotlib.figure.Figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    years = peak_df["water_year"].values.astype(int)
    flows = peak_df["peak_flow_cfs"].values.astype(float)

    # Determine which years are inside the FFA analysis range
    if ffa_year_range is not None:
        ffa_start, ffa_end = ffa_year_range
        in_range = (years >= ffa_start) & (years <= ffa_end)
        has_excluded = bool((~in_range).any())
    else:
        in_range = np.ones(len(years), dtype=bool)
        ffa_start = ffa_end = None
        has_excluded = False

    # Bar chart — years outside FFA range drawn as hollow outline bars
    if has_excluded:
        ax.bar(
            years[in_range],
            flows[in_range],
            color="steelblue",
            alpha=0.75,
            width=0.8,
            label="Annual Peak Flow (in analysis)",
        )
        ax.bar(
            years[~in_range],
            flows[~in_range],
            facecolor="none",
            edgecolor="steelblue",
            linewidth=0.8,
            alpha=0.7,
            width=0.8,
            label="Excluded from FFA",
        )
    else:
        ax.bar(years, flows, color="steelblue", alpha=0.75, width=0.8, label="Annual Peak Flow")

    # MGBT threshold line (solid red, full record width)
    if mgbt_threshold is not None and mgbt_threshold > 0:
        ax.axhline(
            mgbt_threshold,
            color="red",
            linestyle="-",
            linewidth=1.2,
            alpha=0.8,
            label=f"MGBT threshold: {mgbt_threshold:,.0f} cfs",
            zorder=5,
        )

    # Perception threshold lines (lower and/or upper per period)
    if thresholds:
        ax_ylim_top = ax.get_ylim()[1] if ax.get_yscale() == "linear" else None

        for i, thr in enumerate(thresholds):
            start = int(thr.get("start_year", years.min()))
            end = int(thr.get("end_year", years.max()))
            lo = float(thr.get("lower_cfs", thr.get("threshold_cfs", 0)))
            hi_raw = thr.get("upper_cfs", thr.get("threshold_cfs", None))
            hi = float("inf") if hi_raw in (None, float("inf"), "inf", "") else float(hi_raw)

            color = _THRESHOLD_COLORS[i % len(_THRESHOLD_COLORS)]
            ax.axvspan(start - 0.5, end + 0.5, alpha=0.07, color=color)

            # Lower threshold line
            if lo > 0:
                lo_label = f"T\u2081\u208b: {lo:,.0f} cfs ({start}\u2013{end})"
                ax.hlines(
                    lo,
                    start - 0.5,
                    end + 0.5,
                    colors=color,
                    linestyles="--",
                    linewidth=1.5,
                    label=lo_label,
                    zorder=5,
                )

            # Upper threshold line
            if not np.isinf(hi) and hi > 0:
                hi_label = f"T\u2081\u207a: {hi:,.0f} cfs ({start}\u2013{end})"
                ax.hlines(
                    hi,
                    start - 0.5,
                    end + 0.5,
                    colors=color,
                    linestyles=":",
                    linewidth=1.5,
                    label=hi_label,
                    zorder=5,
                )
            elif np.isinf(hi):
                # Mark "T_upper = ∞" as text at top of the span
                ax.annotate(
                    f"T\u1d64\u209a\u209a\u2091\u1d63 = \u221e",
                    xy=(start + (end - start) / 2, 1.0),
                    xycoords=("data", "axes fraction"),
                    fontsize=7,
                    ha="center",
                    va="bottom",
                    color=color,
                )

    # --- Return-period reference lines + max-peak recurrence annotation ---
    if lp3_params is not None:
        mean_log, std_log, skew = lp3_params
        rp_arr = np.array(sorted(set(float(rp) for rp in return_periods)))
        rp_flows = _lp3_quantiles(mean_log, std_log, skew, 1.0 / rp_arr)

        for i, (rp, flow) in enumerate(zip(rp_arr, rp_flows)):
            # Only the first line gets a label, so the whole family collapses
            # to one legend entry instead of one per return period.
            ax.axhline(
                flow,
                color="dimgray",
                linestyle=":",
                linewidth=1.0,
                alpha=0.6,
                zorder=2,
                label="Return Period (LP3 fit)" if i == 0 else None,
            )
            rp_label = f"{rp:g}-yr"
            ax.annotate(
                rp_label,
                xy=(1.0, flow),
                xycoords=("axes fraction", "data"),
                xytext=(4, 0),
                textcoords="offset points",
                fontsize=7,
                color="dimgray",
                ha="left",
                va="center",
                annotation_clip=False,
            )

        if annotate_max_peak and len(flows) > 0:
            max_idx = int(np.argmax(flows))
            max_flow = float(flows[max_idx])
            max_year = int(years[max_idx])
            if max_flow > 0:
                aep_at_max = 1.0 - log_pearson3_cdf(max_flow, mean_log, std_log, skew)
                # An AEP near (or, at the CDF's saturation, exactly) zero makes
                # 1/aep either a meaningless multi-billion-year figure or an
                # outright division by zero. Report a bound past the largest
                # tabulated return period instead, in both cases.
                ri = 1.0 / aep_at_max if aep_at_max > 0 else float("inf")
                if ri > 10 * rp_arr[-1]:
                    ri_line = f"> {rp_arr[-1]:g}-yr"
                else:
                    ri_str = f"{ri:,.0f}" if ri >= 10 else f"{ri:.1f}"
                    ri_line = f"≈ {ri_str}-yr"
                ax.annotate(
                    f"{max_flow:,.0f} cfs ({max_year})\n{ri_line}",
                    xy=(max_year, max_flow),
                    xytext=(0, 12),
                    textcoords="offset points",
                    fontsize=_ANNOT_FONT_SIZE,
                    ha="center",
                    va="bottom",
                    bbox=dict(
                        boxstyle="round,pad=0.2",
                        facecolor="white",
                        edgecolor="lightgray",
                        alpha=0.9,
                    ),
                    arrowprops=dict(arrowstyle="-", color="gray", alpha=0.7),
                    zorder=6,
                )

    ax.set_yscale("log")
    ax.set_xlabel("Water Year", fontsize=_FONT_SIZE)
    ax.set_ylabel("Peak Flow (cfs)", fontsize=_FONT_SIZE)

    # Y-axis: powers of 10, comma-formatted (matches hydrograph style)
    valid_flows = flows[flows > 0]
    if len(valid_flows) > 0:
        min_exp = math.floor(math.log10(valid_flows.min()))
        max_exp = math.ceil(math.log10(valid_flows.max()))
        yticks = [10**i for i in range(min_exp, max_exp + 1)]
        ax.set_yticks(yticks)
        ax.set_yticklabels([f"{int(t):,}" for t in yticks])

    # Title
    if site_name and site_no:
        title = f"Annual Peak Flows\nUSGS {site_no} - {site_name}"
    elif site_no:
        title = f"Annual Peak Flows - USGS {site_no}"
    else:
        title = "Annual Peak Flows"
    ax.set_title(title, fontsize=_FONT_SIZE, fontweight="bold")

    # Record summary annotation (top-left, matches hydrograph style)
    n_years = len(years)
    if has_excluded:
        n_in = int(in_range.sum())
        stats_text = f"n = {n_years} total peaks\n" f"{n_in} in analysis ({ffa_start}–{ffa_end})"
    else:
        year_range = f"{years.min()}–{years.max()}" if n_years > 0 else ""
        stats_text = f"n = {n_years} years\n{year_range}"
    ax.annotate(
        stats_text,
        xy=(_BOX_PADDING, 1 - _BOX_PADDING),
        xycoords="axes fraction",
        fontsize=_ANNOT_FONT_SIZE,
        ha="left",
        va="top",
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="lightgray", alpha=1.0),
    )

    has_valid_thresholds = thresholds and any(
        float(t.get("lower_cfs", t.get("threshold_cfs", 0)) or 0) > 0
        or (
            t.get("upper_cfs") not in (None, "", "inf")
            and not np.isinf(float(t.get("upper_cfs", float("inf"))))
        )
        for t in thresholds
    )
    if has_valid_thresholds or mgbt_threshold or has_excluded or lp3_params is not None:
        ax.legend(loc="upper right", fontsize=_ANNOT_FONT_SIZE)

    ax.grid(True, which="both", alpha=0.3, axis="y")
    fig.tight_layout()

    return fig


#: Former name of :func:`plot_frequency_curve`, kept so the split-out app and
#: any other existing caller keep working across the rename. The "streamlit"
#: suffix was always a misnomer -- this module imports matplotlib and returns a
#: Figure -- and it became actively misleading once the app moved out. Prefer
#: the new name; this alias can go once no pinned consumer uses it.
plot_frequency_curve_streamlit = plot_frequency_curve
