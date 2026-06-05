"""
hydrolib.freq_plot - Frequency curve plotting for Streamlit display.
"""

from __future__ import annotations

import logging
import math
from typing import TYPE_CHECKING, Dict, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

from .core import kfactor_array

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


def plot_frequency_curve_streamlit(
    b17c: "Bulletin17C",
    site_name: Optional[str] = None,
    site_no: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 6),
    skew_curves: Optional[Dict[str, float]] = None,
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
        skew_curves = {"LP3 Fitted Curve": r.skew_used}

    multi = len(skew_curves) > 1

    def aep_to_x(p: float) -> float:
        return norm.ppf(1 - p)

    fig, ax = plt.subplots(figsize=figsize)

    # --- Observed peaks (Weibull plotting positions) ---
    peak_flows = b17c._peak_flows  # noqa: SLF001
    sorted_flows = np.sort(peak_flows)[::-1]
    n_obs = len(sorted_flows)
    weibull_aep = np.arange(1, n_obs + 1) / (n_obs + 1)
    x_obs = [aep_to_x(p) for p in weibull_aep]

    ax.scatter(
        x_obs,
        sorted_flows,
        c="steelblue",
        s=20,
        alpha=0.5,
        zorder=5,
        label="Observed Annual Peaks",
        edgecolors="navy",
        linewidth=0.5,
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
        if label in _SKEW_STYLE:
            color, ls = _SKEW_STYLE[label]
        else:
            color = _FALLBACK_COLORS[fallback_idx % len(_FALLBACK_COLORS)]
            ls = "-"
            fallback_idx += 1

        # Single curve → black "k-" linewidth=1 (matches summary hydrograph median)
        if not multi:
            color, ls = "k", "-"

        Q_curve = _lp3_quantiles(mean_log, std_log, skew_val, aep_fine)
        ax.plot(x_curve, Q_curve, color=color, linestyle=ls, linewidth=1, label=label, zorder=4)

        lower, upper = _lp3_ci(mean_log, std_log, skew_val, n_sys, aep_ci)
        ax.fill_between(x_cl, lower, upper, alpha=0.12, color=color)
        ax.plot(x_cl, lower, color=color, linestyle="--", linewidth=0.7, alpha=0.55)
        ax.plot(x_cl, upper, color=color, linestyle="--", linewidth=0.7, alpha=0.55)

    # --- Low outlier threshold ---
    if r.low_outlier_threshold > 0 and r.n_low_outliers > 0:
        ax.axhline(
            r.low_outlier_threshold,
            color="red",
            linestyle="--",
            alpha=0.7,
            label=f"MGBT Threshold ({r.n_low_outliers} low outliers)",
        )

    # --- Axes ---
    ax.set_yscale("log")
    ax.set_ylabel("Discharge (cfs)", fontsize=_FONT_SIZE)
    ax.set_xlabel("Annual Exceedance Probability", fontsize=_FONT_SIZE)

    # Y-axis ticks: powers of 10, comma-formatted (matches hydrograph plots)
    min_flow = sorted_flows[-1] if n_obs else 1.0
    max_flow = sorted_flows[0] if n_obs else 1e6
    if min_flow > 0:
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
