"""
hydrolib.plots - Plotting utilities for flood frequency analysis
"""

from __future__ import annotations

from typing import TYPE_CHECKING, List, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np

if TYPE_CHECKING:
    from .core import PeakRecord
    from .engine import B17CEngine


def apply_b17c_style():
    """Apply standard B17C plotting style."""
    plt.rcParams.update(
        {
            "figure.dpi": 140,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.grid": True,
            "axes.grid.which": "both",
            "grid.alpha": 0.3,
            "font.size": 10,
        }
    )


def plotting_positions(flows: Sequence[float]) -> np.ndarray:
    """
    Calculate Weibull plotting positions as return periods.

    Parameters
    ----------
    flows : sequence of float
        Flow values (will be sorted in descending order)

    Returns
    -------
    np.ndarray
        Return periods for plotting
    """
    n = len(flows)
    ranks = np.arange(1, n + 1)
    return (n + 1) / ranks


def plot_frequency_curve(
    engine: "B17CEngine",
    records: List["PeakRecord"],
    return_periods: Sequence[float] = (1.5, 2, 5, 10, 25, 50, 100, 200, 500),
    show_bankfull: bool = True,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    figsize: tuple = (10, 6),
) -> plt.Figure:
    """
    Plot flood frequency curve with observed data and confidence intervals.

    Parameters
    ----------
    engine : B17CEngine
        Fitted B17C engine
    records : list of PeakRecord
        Peak flow records used for fitting
    return_periods : sequence of float
        Return periods for fitted curve
    show_bankfull : bool
        Whether to show bankfull range annotation (1.5-2 yr)
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save figure
    figsize : tuple
        Figure size

    Returns
    -------
    plt.Figure
        Matplotlib figure
    """
    apply_b17c_style()

    # Get observed flows
    obs = sorted([r.flow for r in records if r.flow is not None], reverse=True)
    rp_obs = plotting_positions(obs)

    # Get fitted quantiles with CI
    rps = np.array(list(return_periods))
    q = engine.quantiles_with_ci(return_periods)

    est = np.array([q[T]["estimate"] for T in rps])
    lo = np.array([q[T]["lower"] for T in rps])
    hi = np.array([q[T]["upper"] for T in rps])

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot observed data
    ax.scatter(
        rp_obs,
        obs,
        s=40,
        c="blue",
        edgecolors="darkblue",
        linewidth=0.5,
        zorder=5,
        label="Observed Peaks",
    )

    # Plot fitted curve
    ax.plot(rps, est, "b-", linewidth=2, label="LP3 Fit", zorder=4)

    # Plot confidence interval
    ax.fill_between(rps, lo, hi, alpha=0.25, color="blue", label="95% CI")

    # Bankfull annotation
    if show_bankfull and 1.5 in q and 2 in q:
        q15 = q[1.5]["estimate"]
        q2 = q[2]["estimate"]
        ax.axhspan(q15, q2, alpha=0.15, color="green", zorder=1)
        ax.text(
            1.6, (q15 + q2) / 2, "Bankfull", fontsize=9, va="center", ha="left", color="darkgreen"
        )

    # Formatting
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Return Period (years)", fontsize=11)
    ax.set_ylabel("Peak Discharge (cfs)", fontsize=11)

    if title:
        ax.set_title(title, fontsize=12, fontweight="bold")
    else:
        ax.set_title("Flood Frequency Curve (Bulletin 17C)", fontsize=12, fontweight="bold")

    # Set x-axis ticks
    ax.set_xticks([1.5, 2, 5, 10, 25, 50, 100, 200, 500])
    ax.set_xticklabels(["1.5", "2", "5", "10", "25", "50", "100", "200", "500"])

    # Add stats annotation
    if engine.params:
        mu, sigma, skew = engine.params
        stats_text = (
            f"n = {engine.n}\n"
            f"Mean(log Q) = {mu:.4f}\n"
            f"Std(log Q) = {sigma:.4f}\n"
            f"Skew = {skew:.3f}"
        )
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
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig


def plot_multi_site_comparison(
    results: dict,
    return_period: float = 100,
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    figsize: tuple = (10, 6),
) -> plt.Figure:
    """
    Plot comparison of quantiles across multiple sites.

    Parameters
    ----------
    results : dict
        Output from run_multi_site or analyze_sites
    return_period : float
        Return period to compare
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save figure
    figsize : tuple
        Figure size

    Returns
    -------
    plt.Figure
        Matplotlib figure
    """
    apply_b17c_style()

    # Extract data for sites without errors
    sites = []
    estimates = []
    lowers = []
    uppers = []

    for site, result in results.items():
        if "error" not in result and return_period in result["ci"]:
            sites.append(site)
            ci = result["ci"][return_period]
            estimates.append(ci["estimate"])
            lowers.append(ci["lower"])
            uppers.append(ci["upper"])

    if not sites:
        raise ValueError("No valid results to plot")

    # Sort by estimate
    sort_idx = np.argsort(estimates)
    sites = [sites[i] for i in sort_idx]
    estimates = [estimates[i] for i in sort_idx]
    lowers = [lowers[i] for i in sort_idx]
    uppers = [uppers[i] for i in sort_idx]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    y_pos = np.arange(len(sites))
    errors = [
        [e - l for e, l in zip(estimates, lowers)],
        [u - e for e, u in zip(estimates, uppers)],
    ]

    ax.barh(y_pos, estimates, xerr=errors, align="center", color="steelblue", alpha=0.7, capsize=3)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(sites)
    ax.set_xlabel(f"Q{int(return_period)} Flow (cfs)", fontsize=11)
    ax.set_ylabel("Site", fontsize=11)

    if title:
        ax.set_title(title, fontsize=12, fontweight="bold")
    else:
        ax.set_title(f"{int(return_period)}-Year Flood Comparison", fontsize=12, fontweight="bold")

    ax.set_xscale("log")
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig
