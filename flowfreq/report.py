"""
flowfreq.report - Technical report generation
"""

from __future__ import annotations

from datetime import datetime
from typing import Dict, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .bulletin17c import Bulletin17C, FloodFrequencyAnalysis
from .core import AnalysisMethod, FrequencyResults
from .hydrograph import Hydrograph
from .usgs import USGSgage

# Alias for type hints
USGSGage = USGSgage


class HydroReport:
    """Generate technical report for hydrologic analysis."""

    def __init__(self, gage: USGSGage, analysis: Union[Bulletin17C, FloodFrequencyAnalysis]):
        self._gage = gage

        if isinstance(analysis, Bulletin17C):
            self._analysis = analysis._analyzer
            self._results = analysis.results
        else:
            self._analysis = analysis
            self._results = analysis.results

        self._figures: Dict[str, str] = {}

    @property
    def figures(self) -> Dict[str, str]:
        return self._figures.copy()

    def generate_all_figures(self, output_dir: str = ".") -> Dict[str, str]:
        """Generate all figures for the report."""
        import os

        os.makedirs(output_dir, exist_ok=True)

        site_name = self._gage.site_name
        site_no = self._gage.site_no

        if self._gage.daily_data is not None:
            path = os.path.join(output_dir, "fig1_daily_timeseries.png")
            Hydrograph.plot_daily_timeseries(
                self._gage.daily_data, site_name, site_no, save_path=path
            )
            plt.close()
            self._figures["daily_timeseries"] = path

        if self._gage.daily_data is not None:
            path = os.path.join(output_dir, "fig2_summary_hydrograph.png")
            Hydrograph.plot_summary_hydrograph(
                self._gage.daily_data, site_name, site_no, save_path=path
            )
            plt.close()
            self._figures["summary_hydrograph"] = path

        if self._gage.peak_data is not None:
            path = os.path.join(output_dir, "fig3_annual_peaks.png")
            self._plot_annual_peaks(save_path=path)
            plt.close()
            self._figures["annual_peaks"] = path

        path = os.path.join(output_dir, "fig4_frequency_curve.png")
        self._analysis.plot_frequency_curve(site_name, site_no, save_path=path)
        plt.close()
        self._figures["frequency_curve"] = path

        return self._figures

    def _plot_annual_peaks(self, save_path: str = None) -> plt.Figure:
        """Plot annual peak flow bar chart."""
        fig, ax = plt.subplots(figsize=(12, 5))

        df = self._gage.peak_data
        ax.bar(
            df["water_year"],
            df["peak_flow_cfs"],
            color="steelblue",
            edgecolor="darkblue",
            linewidth=0.5,
        )

        ax.set_xlabel("Water Year", fontsize=11)
        ax.set_ylabel("Peak Discharge (cfs)", fontsize=11)

        title = "Annual Peak Streamflow"
        if self._gage.site_name and self._gage.site_no:
            title = f"Annual Peak Streamflow\nUSGS {self._gage.site_no} - {self._gage.site_name}"
        ax.set_title(title, fontsize=12, fontweight="bold")

        ax.grid(True, axis="y", alpha=0.3)

        mean_peak = df["peak_flow_cfs"].mean()
        max_peak = df["peak_flow_cfs"].max()
        max_year = df.loc[df["peak_flow_cfs"].idxmax(), "water_year"]

        ax.annotate(
            f"Mean: {mean_peak:,.0f} cfs\nMax: {max_peak:,.0f} cfs (WY {max_year})",
            xy=(0.02, 0.98),
            xycoords="axes fraction",
            fontsize=9,
            ha="left",
            va="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches="tight")

        return fig

    def generate_report_text(self) -> str:
        """Generate technical report in markdown format."""
        r = self._results
        site_no = self._gage.site_no
        site_name = self._gage.site_name or "Unknown"
        da = self._gage.drainage_area

        method_name = (
            "Expected Moments Algorithm (EMA)"
            if r.method == AnalysisMethod.EMA
            else "Method of Moments (MOM)"
        )

        report = f"""# Flood Frequency Analysis Report

## USGS {site_no} - {site_name}

**Report Date:** {datetime.now().strftime('%B %d, %Y')}

**Analysis Method:** {method_name}

---

## 1. Introduction

This report presents a flood frequency analysis for USGS streamgage {site_no} ({site_name}) following the methodology outlined in Bulletin 17C, "Guidelines for Determining Flood Flow Frequency" (USGS, 2019).

## 2. Site Description

| Parameter | Value |
|-----------|-------|
| USGS Site Number | {site_no} |
| Station Name | {site_name} |
"""
        if da:
            report += f"| Drainage Area | {da:,.1f} sq mi |\n"

        if self._gage.peak_data is not None:
            df = self._gage.peak_data
            report += (
                f"| Period of Record | WY {df['water_year'].min()} - {df['water_year'].max()} |\n"
            )
            report += f"| Number of Peak Flow Records | {len(df)} |\n"

        report += f"""
## 3. Statistical Analysis

### 3.1 Log-Pearson Type III Distribution Parameters

| Parameter | Value |
|-----------|-------|
| Total Observations | {r.n_peaks} |
| Systematic Record | {r.n_systematic} |
| Historical Observations | {r.n_historical} |
| Censored Observations | {r.n_censored} |
| Low Outliers (MGB) | {r.n_low_outliers} |
| Mean of Log Q | {r.mean_log:.4f} |
| Standard Deviation of Log Q | {r.std_log:.4f} |
| Station Skew Coefficient | {r.skew_station:.4f} |
"""

        if r.skew_regional is not None:
            report += f"| Regional Skew Coefficient | {r.skew_regional:.4f} |\n"
            report += f"| Weighted Skew Coefficient | {r.skew_weighted:.4f} |\n"

        if r.method == AnalysisMethod.EMA:
            report += f"| EMA Iterations | {r.ema_iterations} |\n"
            report += f"| EMA Converged | {'Yes' if r.ema_converged else 'No'} |\n"

        report += f"""
### 3.2 Low Outlier Analysis

| Parameter | Value |
|-----------|-------|
| MGB Critical Value (K_n) | {r.mgb_critical_value:.4f} |
| Low Outlier Threshold | {r.low_outlier_threshold:,.0f} cfs |
| Number of Low Outliers | {r.n_low_outliers} |

## 4. Flood Frequency Estimates

**Table 1. Flood Frequency Quantiles with 90% Confidence Intervals**

| AEP (%) | Return Period (yr) | Flow (cfs) | 5% Limit (cfs) | 95% Limit (cfs) |
|---------|-------------------|------------|----------------|-----------------|
"""

        cl = r.confidence_limits
        for _, row in cl.iterrows():
            aep_pct = row["aep"] * 100
            aep_str = f"{aep_pct:.1f}" if aep_pct >= 1 else f"{aep_pct:.2f}"
            rp = row["return_period"]
            rp_str = f"{rp:,.0f}" if rp >= 10 else f"{rp:.1f}"

            report += f"| {aep_str} | {rp_str} | {row['flow_cfs']:,.0f} | {row['lower_5pct']:,.0f} | {row['upper_5pct']:,.0f} |\n"

        report += """
## 5. References

England, J.F., Jr., et al., 2019, Guidelines for determining flood flow frequency—Bulletin 17C: U.S. Geological Survey Techniques and Methods, book 4, chap. B5, 148 p., https://doi.org/10.3133/tm4B5.

---

## Appendix A: Figures

- **Figure 1.** Mean Daily Streamflow Time Series
- **Figure 2.** Summary Hydrograph by Day of Water Year
- **Figure 3.** Annual Peak Streamflow
- **Figure 4.** Flood Frequency Curve with 90% Confidence Limits

## Appendix B: Annual Peak Flow Data

| Water Year | Peak Date | Peak Flow (cfs) |
|------------|-----------|-----------------|
"""

        if self._gage.peak_data is not None:
            for _, row in self._gage.peak_data.iterrows():
                date_str = (
                    row["peak_date"].strftime("%Y-%m-%d") if pd.notna(row["peak_date"]) else "N/A"
                )
                report += (
                    f"| {int(row['water_year'])} | {date_str} | {row['peak_flow_cfs']:,.0f} |\n"
                )

        return report

    def save_report(self, output_path: str):
        """Save report to markdown file."""
        with open(output_path, "w") as f:
            f.write(self.generate_report_text())
