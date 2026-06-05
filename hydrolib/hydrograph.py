"""
hydrolib.hydrograph - Hydrograph analysis and plotting
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import ClassVar, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.dates import DateFormatter

plt_color: str = "steelblue"


class Hydrograph:
    """Class for hydrograph analysis and plotting."""

    FONT_SIZE: ClassVar[int] = 12
    ANNOT_FONT_SIZE: ClassVar[int] = 9  # Smaller font for annotations

    # Consistent spacing for text boxes from axes edges (axes fraction coordinates)
    BOX_PADDING: ClassVar[float] = 0.03

    MONTH_STARTS: ClassVar[List[int]] = [1, 32, 62, 93, 123, 154, 184, 215, 246, 274, 305, 335]
    MONTH_LABELS: ClassVar[List[str]] = [
        "Oct",
        "Nov",
        "Dec",
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
    ]

    @staticmethod
    @lru_cache(maxsize=1024)
    def day_of_water_year(month: int, day: int) -> int:
        """Convert month/day to day of water year (1-366, starting Oct 1)."""
        days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        wy_month = (month - 10) % 12
        dowy = sum(days_in_month[(10 + i - 1) % 12] for i in range(wy_month)) + day
        return dowy

    @classmethod
    def _compute_dowy_series(cls, dates: pd.DatetimeIndex) -> np.ndarray:
        """Compute day of water year for a date series."""
        return np.array([cls.day_of_water_year(d.month, d.day) for d in dates])

    @classmethod
    def plot_daily_timeseries(
        cls,
        daily_data: pd.DataFrame,
        site_name: str = None,
        site_no: str = None,
        save_path: str = None,
        figsize: Tuple[int, int] = (12, 5),
        color: str = plt_color,
        por_start: str = None,
        por_end: str = None,
    ) -> plt.Figure:
        """Plot daily flow time series."""
        fig, ax = plt.subplots(figsize=figsize)

        ax.plot(daily_data.index, daily_data["flow_cfs"], color=color, linewidth=0.5, alpha=0.8)
        ax.set_yscale("log")
        ax.set_ylabel("Discharge (cfs)", fontsize=cls.FONT_SIZE)
        ax.set_xlabel("Date", fontsize=cls.FONT_SIZE)

        title = "Mean Daily Streamflow"
        if site_name and site_no:
            title = f"Mean Daily Streamflow\nUSGS {site_no} - {site_name}"
        elif site_no:
            title = f"Mean Daily Streamflow - USGS {site_no}"
        ax.set_title(title, fontsize=cls.FONT_SIZE, fontweight="bold")

        ax.grid(True, which="both", alpha=0.3)
        ax.xaxis.set_major_formatter(DateFormatter("%Y"))

        # Set y-ticks to powers of 10
        min_flow = daily_data["flow_cfs"].min()
        max_flow = daily_data["flow_cfs"].max()
        if min_flow > 0:
            min_exp = math.floor(math.log10(min_flow))
            max_exp = math.ceil(math.log10(max_flow))
            ticks = [10**i for i in range(min_exp, max_exp + 1)]
            ax.set_yticks(ticks)
            ax.set_yticklabels([f"{int(tick):,}" for tick in ticks])

        # Format dates for annotation (cross-platform, no leading zeros)
        d_min = daily_data.index.min()
        d_max = daily_data.index.max()
        plot_start = f"{d_min.month}/{d_min.day}/{d_min.year}"
        plot_end = f"{d_max.month}/{d_max.day}/{d_max.year}"

        # Build annotation text
        if por_start and por_end:
            annot_text = f"POR: {por_start} - {por_end}\nPlot: {plot_start} - {plot_end}"
        else:
            annot_text = f"POR: {plot_start} - {plot_end}"

        ax.annotate(
            annot_text,
            xy=(cls.BOX_PADDING, 1 - cls.BOX_PADDING),
            xycoords="axes fraction",
            fontsize=cls.ANNOT_FONT_SIZE,
            ha="left",
            va="top",
            bbox=dict(
                boxstyle="round,pad=0.3", facecolor="white", edgecolor="lightgray", alpha=1.0
            ),
        )

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches="tight")

        return fig

    @classmethod
    def plot_summary_hydrograph(
        cls,
        daily_data: pd.DataFrame,
        site_name: str = None,
        site_no: str = None,
        save_path: str = None,
        figsize: Tuple[int, int] = (10, 6),
        percentiles: List[int] = None,
        por_start: str = None,
        por_end: str = None,
    ) -> plt.Figure:
        """Plot summary hydrograph with day of water year on x-axis."""
        if percentiles is None:
            percentiles = [10, 25, 50, 75, 90]

        df = daily_data.copy()
        df["dowy"] = cls._compute_dowy_series(df.index)

        stats_df = df.groupby("dowy")["flow_cfs"].agg(
            ["mean", "min", "max"] + [lambda x, p=p: np.nanpercentile(x, p) for p in percentiles]
        )
        stats_df.columns = ["mean", "min", "max"] + [f"p{p}" for p in percentiles]

        fig, ax = plt.subplots(figsize=figsize)

        if len(percentiles) >= 4:
            ax.fill_between(
                stats_df.index,
                stats_df[f"p{percentiles[0]}"],
                stats_df[f"p{percentiles[-1]}"],
                alpha=0.3,
                color=plt_color,
                label=f"{percentiles[0]}-{percentiles[-1]}th percentile",
            )
            ax.fill_between(
                stats_df.index,
                stats_df[f"p{percentiles[1]}"],
                stats_df[f"p{percentiles[-2]}"],
                alpha=0.4,
                color=plt_color,
                label=f"{percentiles[1]}-{percentiles[-2]}th percentile",
            )

        ax.plot(stats_df.index, stats_df["p50"], "k-", linewidth=1, label="Median")
        ax.plot(
            stats_df.index,
            stats_df["min"],
            ":",
            linewidth=1,
            color=plt_color,
            alpha=0.4,
            label="Min",
        )
        ax.plot(
            stats_df.index,
            stats_df["max"],
            ":",
            linewidth=1,
            color=plt_color,
            alpha=0.4,
            label="Max",
        )

        ax.set_yscale("log")
        ax.set_ylabel("Discharge (cfs)", fontsize=cls.FONT_SIZE)
        ax.set_xlabel("Day of Water Year", fontsize=cls.FONT_SIZE)

        ax.set_xticks(cls.MONTH_STARTS)
        ax.set_xticklabels(cls.MONTH_LABELS)
        ax.set_xlim(1, 366)

        # Set y-ticks to powers of 10
        min_flow = daily_data["flow_cfs"].min()
        max_flow = daily_data["flow_cfs"].max()
        if min_flow > 0:
            min_exp = math.floor(math.log10(min_flow))
            max_exp = math.ceil(math.log10(max_flow))
            ticks = [10**i for i in range(min_exp, max_exp + 1)]
            ax.set_yticks(ticks)
            ax.set_yticklabels([f"{int(tick):,}" for tick in ticks])

        title = "Summary Hydrograph"
        if site_name and site_no:
            title = f"Summary Hydrograph\nUSGS {site_no} - {site_name}"
        elif site_no:
            title = f"Summary Hydrograph - USGS {site_no}"
        ax.set_title(title, fontsize=cls.FONT_SIZE, fontweight="bold")

        ax.legend(
            loc="upper right",
            fontsize=cls.ANNOT_FONT_SIZE,  # Smaller legend font
        )
        ax.grid(True, which="both", alpha=0.3)

        # Format dates for annotation (cross-platform, no leading zeros)
        d_min = daily_data.index.min()
        d_max = daily_data.index.max()
        plot_start = f"{d_min.month}/{d_min.day}/{d_min.year}"
        plot_end = f"{d_max.month}/{d_max.day}/{d_max.year}"
        n_years = len(df.index.year.unique())

        # Build annotation text
        if por_start and por_end:
            annot_text = f"POR: {por_start} - {por_end}\nPlot: {plot_start} - {plot_end}\n({n_years} years plotted)"
        else:
            annot_text = f"POR: {plot_start} - {plot_end}\n({n_years} years)"

        ax.annotate(
            annot_text,
            xy=(cls.BOX_PADDING, 1 - cls.BOX_PADDING),
            xycoords="axes fraction",
            fontsize=cls.ANNOT_FONT_SIZE,
            ha="left",
            va="top",
            bbox=dict(
                boxstyle="round,pad=0.3", facecolor="white", edgecolor="lightgray", alpha=1.0
            ),
        )

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches="tight")

        return fig

    @classmethod
    def get_summary_stats(
        cls,
        daily_data: pd.DataFrame,
        percentiles: List[int] = None,
    ) -> pd.DataFrame:
        """Compute summary hydrograph statistics by day of water year.

        Returns a DataFrame with columns for each statistic (mean, min, max, percentiles)
        and rows for each day of the water year (1-366).
        """
        if percentiles is None:
            percentiles = [10, 25, 50, 75, 90]

        df = daily_data.copy()
        df["dowy"] = cls._compute_dowy_series(df.index)

        stats_df = df.groupby("dowy")["flow_cfs"].agg(
            ["mean", "min", "max"] + [lambda x, p=p: np.nanpercentile(x, p) for p in percentiles]
        )
        stats_df.columns = ["mean", "min", "max"] + [f"p{p}" for p in percentiles]
        stats_df.index.name = "day_of_water_year"

        return stats_df.reset_index()

    @classmethod
    def plot_flow_duration_curve(
        cls,
        daily_data: pd.DataFrame,
        site_name: str = None,
        site_no: str = None,
        save_path: str = None,
        table_path: str = None,
        figsize: Tuple[int, int] = (10, 6),
        por_start: str = None,
        por_end: str = None,
    ) -> Tuple[plt.Figure, pd.DataFrame]:
        """Plot flow duration curve and return flow duration statistics table."""
        flows = daily_data["flow_cfs"].dropna().values
        flows_sorted = np.sort(flows)[::-1]  # descending
        n = len(flows_sorted)
        exceedance_pct = np.arange(1, n + 1) / n * 100

        fig, ax = plt.subplots(figsize=figsize)

        ax.plot(exceedance_pct, flows_sorted, color=plt_color, linewidth=1.5)
        ax.set_yscale("log")
        ax.set_xlabel("Percent of Time Exceeded (%)", fontsize=cls.FONT_SIZE)
        ax.set_ylabel("Discharge (cfs)", fontsize=cls.FONT_SIZE)

        title = "Flow Duration Curve"
        if site_name and site_no:
            title = f"Flow Duration Curve\nUSGS {site_no} - {site_name}"
        elif site_no:
            title = f"Flow Duration Curve - USGS {site_no}"
        ax.set_title(title, fontsize=cls.FONT_SIZE, fontweight="bold")

        ax.grid(True, which="both", alpha=0.3)

        # Set y-ticks to powers of 10
        min_flow = flows.min()
        max_flow = flows.max()
        if min_flow > 0:
            min_exp = math.floor(math.log10(min_flow))
            max_exp = math.ceil(math.log10(max_flow))
            ticks = [10**i for i in range(min_exp, max_exp + 1)]
            ax.set_yticks(ticks)
            ax.set_yticklabels([f"{int(tick):,}" for tick in ticks])

        # Set x-ticks
        ax.set_xscale("linear")
        ax.set_xlim(0, 100)
        ax.set_xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])

        # Format dates for annotation (cross-platform, no leading zeros)
        d_min = daily_data.index.min()
        d_max = daily_data.index.max()
        plot_start = f"{d_min.month}/{d_min.day}/{d_min.year}"
        plot_end = f"{d_max.month}/{d_max.day}/{d_max.year}"
        n_years = len(daily_data.index.year.unique())

        # Build annotation text
        if por_start and por_end:
            annot_text = f"POR: {por_start} - {por_end}\nPlot: {plot_start} - {plot_end}\n({n_years} years plotted)"
        else:
            annot_text = f"POR: {plot_start} - {plot_end}\n({n_years} years)"

        ax.annotate(
            annot_text,
            xy=(cls.BOX_PADDING, 1 - cls.BOX_PADDING),
            xycoords="axes fraction",
            fontsize=cls.ANNOT_FONT_SIZE,
            ha="left",
            va="top",
            bbox=dict(
                boxstyle="round,pad=0.3", facecolor="white", edgecolor="lightgray", alpha=1.0
            ),
        )

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches="tight")

        # Compute flow duration statistics
        percentiles_exceeded = [1, 5, 10, 20, 50, 80, 90, 95, 99]
        flow_stats = []
        for pct in percentiles_exceeded:
            flow_value = np.percentile(flows, 100 - pct)
            flow_stats.append({"Percent Exceeded": f"{pct}%", "Flow (cfs)": flow_value})

        stats_df = pd.DataFrame(flow_stats)

        if table_path:
            stats_df.to_csv(table_path, index=False)

        return fig, stats_df
