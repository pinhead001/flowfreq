"""
Streamlit app for USGS Daily Flow Hydrograph Generator

Run with: streamlit run app/streamlit_app.py
"""

import io
import sys
import zipfile
from datetime import datetime
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
from scipy import stats

matplotlib.use("Agg")  # Use non-interactive backend

# Add parent directory to path for hydrolib import (needed for Streamlit Cloud)
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.ffa_export import export_comparison_csv, export_ffa_to_zip
from app.ffa_runner import format_parameters_df, format_quantile_df

# Import hydrolib
from hydrolib import Hydrograph, __version__
from hydrolib.freq_plot import plot_frequency_curve_streamlit
from hydrolib.usgs import USGSgage
from hydrolib.workflow import SKEW_OPTIONS, build_skew_curves_dict, compute_skew_tables, run_ffa

st.set_page_config(
    page_title="USGS Hydrograph-erator",
    page_icon=":ocean:",
    layout="wide",
)

st.title("USGS Hydrograph-erator")
st.markdown(
    "Generate daily flow plots for USGS gages using [hydrolib](https://github.com/pinhead001/hydrolib)"
)

# Sidebar inputs
st.sidebar.header("Input Parameters")

# Gage input - single or multiple
input_mode = st.sidebar.radio("Input Mode", ["Single Gage", "Multiple Gages"])

if input_mode == "Single Gage":
    gage_input = st.sidebar.text_input("USGS Gage Number", value="09355500")
    gage_list = [gage_input.strip()] if gage_input else []
else:
    gage_input = st.sidebar.text_area(
        "USGS Gage Numbers (one per line or comma-separated)",
        value="09355500\n11066460",
        height=100,
    )
    # Parse input
    gage_list = [g.strip() for g in gage_input.replace(",", "\n").split("\n") if g.strip()]

st.sidebar.markdown(f"**{len(gage_list)} gage(s) selected**")

# Plot options
st.sidebar.header("Plot Options")
show_timeseries = st.sidebar.checkbox("Daily Time Series", value=True)
show_summary = st.sidebar.checkbox("Summary Hydrograph", value=True)
show_fdc = st.sidebar.checkbox("Flow Duration Curve", value=True)

# Y-axis scale options
scale_mode = st.sidebar.radio("Y-Axis Scale", ["Linear", "Log", "Per Plot"], horizontal=True)

if scale_mode == "Per Plot":
    st.sidebar.markdown("**Scale per plot:**")
    yscale_timeseries = (
        "log" if st.sidebar.checkbox("Log: Daily Time Series", value=False) else "linear"
    )
    yscale_summary = (
        "log" if st.sidebar.checkbox("Log: Summary Hydrograph", value=False) else "linear"
    )
    yscale_fdc = "log" if st.sidebar.checkbox("Log: Flow Duration Curve", value=False) else "linear"
    yscale_freq = (
        "log" if st.sidebar.checkbox("Log: Flood Frequency Curve", value=False) else "linear"
    )
    yscale_peak = (
        "log" if st.sidebar.checkbox("Log: Peak Flow Time Series", value=False) else "linear"
    )
else:
    yscale_timeseries = yscale_summary = yscale_fdc = yscale_freq = yscale_peak = scale_mode.lower()

# Date Range Filters (shown when data is loaded)
# Initialize defaults
plot_start = None
plot_end = None
peak_start_year = None
peak_end_year = None

# Create placeholder for date range controls (will be populated after data loads)
date_range_container = st.sidebar.container()

# Flood Frequency Analysis section
st.sidebar.markdown("---")
st.sidebar.subheader("Flood Frequency Analysis")
enable_ffa = st.sidebar.checkbox("Enable Flood Frequency Analysis", value=False)
if enable_ffa:
    regional_skew = st.sidebar.number_input(
        "Regional Skew",
        value=-0.302,
        step=0.001,
        format="%.3f",
        help="Generalized skew from Bulletin 17C Appendix or regional study. Nationwide default: -0.302",
    )
    regional_skew_se = st.sidebar.number_input(
        "Regional Skew SE",
        value=0.55,
        step=0.01,
        format="%.2f",
        help="Standard error of the regional skew estimate. Nationwide default: 0.55",
    )
    show_freq_curve = st.sidebar.checkbox("Frequency Curve", value=True)
    if show_freq_curve:
        label_max_on_freq = st.sidebar.checkbox("Label Max Flow on Curve", value=False)
    else:
        label_max_on_freq = False
    show_peak_timeseries = st.sidebar.checkbox("Peak Flow Time Series", value=True)
    if show_peak_timeseries:
        quantile_options = [1.5, 2, 5, 10, 25, 50, 100, 200, 500]
        show_quantile_lines = st.sidebar.multiselect(
            "Show Return Period Lines",
            options=quantile_options,
            default=[1.5, 25, 100],
            help="Horizontal lines showing flood frequency quantiles",
        )
        show_max_ri = st.sidebar.checkbox("Show Max Peak RI Estimate", value=False)
    else:
        show_quantile_lines = []
        show_max_ri = False
    st.sidebar.markdown("**Skew Options**")
    skew_station_on = st.sidebar.checkbox("Station Skew", value=False)
    skew_weighted_on = st.sidebar.checkbox("Weighted Skew", value=True)
    skew_regional_on = st.sidebar.checkbox("Regional Skew", value=False)
    st.sidebar.markdown("**Low Outliers (PILFs)**")
    pilf_override = st.sidebar.number_input(
        "PILF Threshold Override (cfs)",
        min_value=0.0,
        value=0.0,
        step=100.0,
        format="%.0f",
        help=(
            "Leave at 0 to let the Multiple Grubbs-Beck Test pick the threshold. "
            "A positive value overrides it: every peak below this flow is treated "
            "as a potentially influential low flood and censored out of the fit."
        ),
    )
else:
    regional_skew = -0.302
    regional_skew_se = 0.55
    show_freq_curve = False
    label_max_on_freq = False
    show_peak_timeseries = False
    show_quantile_lines = []
    show_max_ri = False
    skew_station_on = False
    skew_weighted_on = True
    skew_regional_on = False
    pilf_override = 0.0

# Download button
download_data = st.sidebar.button("Download Data", type="primary")

# Initialize session state
if "gage_data" not in st.session_state:
    st.session_state.gage_data = {}  # Stores full daily data for each gage
if "gage_info" not in st.session_state:
    st.session_state.gage_info = {}  # Stores site info for each gage
if "figures" not in st.session_state:
    st.session_state.figures = {}
if "expanded_gage" not in st.session_state:
    st.session_state.expanded_gage = None
if "expanded_plot_idx" not in st.session_state:
    st.session_state.expanded_plot_idx = 0
if "peak_data" not in st.session_state:
    st.session_state.peak_data = {}
if "ffa_results" not in st.session_state:
    st.session_state.ffa_results = {}
if "prev_scale_settings" not in st.session_state:
    st.session_state.prev_scale_settings = None
if "prev_ffa_inputs" not in st.session_state:
    st.session_state.prev_ffa_inputs = None

# Clear figures cache when scale settings change
current_scale_settings = (
    scale_mode,
    yscale_timeseries,
    yscale_summary,
    yscale_fdc,
    yscale_freq,
    yscale_peak,
)
if st.session_state.prev_scale_settings != current_scale_settings:
    st.session_state.figures = {}
    st.session_state.prev_scale_settings = current_scale_settings


def get_plot_list(site_no):
    """Get list of available plots for a gage (excluding stats dataframe)."""
    if site_no not in st.session_state.figures:
        return []
    return [
        name
        for name in st.session_state.figures[site_no].keys()
        if not name.endswith("_stats") and not name.endswith("_data")
    ]


def get_plot_display_name(plot_key):
    """Convert plot key to display name."""
    names = {
        "daily_timeseries": "Daily Time Series",
        "summary_hydrograph": "Summary Hydrograph",
        "flow_duration_curve": "Flow Duration Curve",
        "frequency_curve": "Frequency Curve",
    }
    return names.get(plot_key, plot_key)


def format_date(dt):
    """Format datetime as M/D/YYYY string."""
    return f"{dt.month}/{dt.day}/{dt.year}"


def estimate_ri_from_lp3(flow, mean_log, std_log, skew):
    """Estimate return interval for a flow using LP3 parameters."""
    if flow <= 0 or std_log <= 0:
        return None
    log_flow = np.log10(flow)
    k = (log_flow - mean_log) / std_log
    try:
        non_exceed_prob = stats.pearson3.cdf(k, skew)
        exceed_prob = 1 - non_exceed_prob
        if exceed_prob > 0:
            return 1 / exceed_prob
    except:
        pass
    return None


def plot_peak_timeseries(
    peak_df,
    site_name,
    site_no,
    yscale="linear",
    quantiles=None,
    max_ri_info=None,
    pilf_threshold=None,
    pilf_source=None,
):
    """Plot annual peak flows with optional quantile reference lines.

    quantiles: dict of {return_period: flow_value} to draw as horizontal lines
    pilf_threshold: low-outlier cut in cfs; peaks below it are drawn hollow
    pilf_source: where that cut came from, for the legend ("MGBT"/"override")
    """
    fig, ax = plt.subplots(figsize=(10, 4))

    years = peak_df["water_year"].values
    flows = peak_df["peak_flow_cfs"].values
    # Peaks below the PILF cut are censored out of the fit. Drawing them hollow
    # is the whole point of exposing the override: without it the control
    # changes numbers in a table and nothing a user can see on the record.
    censored = (
        np.asarray(flows) < pilf_threshold if pilf_threshold else np.zeros(len(flows), dtype=bool)
    )
    if censored.any():
        ax.bar(
            years[~censored],
            flows[~censored],
            color="steelblue",
            alpha=0.7,
            label="Peak used in fit",
        )
        ax.bar(
            years[censored],
            flows[censored],
            facecolor="none",
            edgecolor="steelblue",
            linewidth=0.8,
            alpha=0.7,
            label=f"Low outlier, censored ({int(censored.sum())})",
        )
        ax.axhline(
            pilf_threshold,
            color="red",
            linestyle="-",
            linewidth=1.2,
            alpha=0.8,
            label=f"PILF threshold {pilf_threshold:,.0f} cfs ({pilf_source or 'MGBT'})",
        )
        ax.legend(fontsize=8, loc="upper left")
    else:
        ax.bar(years, flows, color="steelblue", alpha=0.7)

    # Add quantile lines with labels above
    if quantiles:
        x_min = peak_df["water_year"].min()
        for rp, flow in sorted(quantiles.items()):
            ax.axhline(y=flow, linestyle="--", color="#404040", alpha=0.9)
            rp_str = f"{rp:g}"  # Format without trailing zeros
            ax.text(
                x_min,
                flow,
                f" {rp_str}-yr: {flow:,.0f}",
                va="bottom",
                ha="left",
                fontsize=8,
                color="#404040",
            )

    # Add max RI annotation if provided
    if max_ri_info:
        max_flow = max_ri_info.get("flow")
        max_year = max_ri_info.get("year")
        max_ri = max_ri_info.get("ri")
        if max_flow and max_year and max_ri:
            ri_str = f"{max_ri:,.0f}" if max_ri >= 10 else f"{max_ri:.1f}"
            annot_text = f"{max_flow:,.0f} cfs\n≈ {ri_str}-yr"
            ax.annotate(
                annot_text,
                xy=(max_year, max_flow),
                xytext=(10, -15),
                textcoords="offset points",
                fontsize=8,
                ha="left",
                va="top",
                bbox=dict(
                    boxstyle="round,pad=0.2", facecolor="white", edgecolor="lightgray", alpha=0.9
                ),
            )

    ax.set_yscale(yscale)
    ax.set_xlabel("Water Year")
    ax.set_ylabel("Peak Flow (cfs)")
    title = "Annual Peak Flows"
    if site_name and site_no:
        title = f"Annual Peak Flows\nUSGS {site_no} - {site_name}"
    elif site_no:
        title = f"Annual Peak Flows - USGS {site_no}"
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    return fig


def generate_plots(
    site_no,
    plot_data,
    gage_info,
    por_start_str=None,
    por_end_str=None,
    yscale_ts="linear",
    yscale_sh="linear",
    yscale_fdc="linear",
):
    """Generate plots for a gage using the provided data."""
    site_figs = {}

    # Daily time series
    if show_timeseries:
        fig1 = Hydrograph.plot_daily_timeseries(
            plot_data,
            site_name=gage_info.get("name"),
            site_no=site_no,
            figsize=(10, 4),
            por_start=por_start_str,
            por_end=por_end_str,
            yscale=yscale_ts,
        )
        site_figs["daily_timeseries"] = fig1

    # Summary hydrograph
    if show_summary:
        fig2 = Hydrograph.plot_summary_hydrograph(
            plot_data,
            site_name=gage_info.get("name"),
            site_no=site_no,
            figsize=(10, 4),
            percentiles=[10, 25, 50, 75, 90],
            por_start=por_start_str,
            por_end=por_end_str,
            yscale=yscale_sh,
        )
        site_figs["summary_hydrograph"] = fig2
        # Store summary stats for export
        site_figs["summary_stats"] = Hydrograph.get_summary_stats(plot_data)

    # Flow duration curve
    if show_fdc:
        fig3, stats_df = Hydrograph.plot_flow_duration_curve(
            plot_data,
            site_name=gage_info.get("name"),
            site_no=site_no,
            figsize=(10, 4),
            por_start=por_start_str,
            por_end=por_end_str,
            yscale=yscale_fdc,
        )
        site_figs["flow_duration_curve"] = fig3
        site_figs["fdc_stats"] = stats_df

    return site_figs


# Download data when button clicked
if download_data and gage_list:
    st.session_state.gage_data = {}
    st.session_state.gage_info = {}
    st.session_state.figures = {}
    st.session_state.peak_data = {}
    st.session_state.ffa_results = {}

    progress_bar = st.progress(0)
    status_text = st.empty()

    for idx, site_no in enumerate(gage_list):
        status_text.text(f"Fetching site info for {site_no}...")

        try:
            # Create gage object and fetch site info (including POR)
            gage = USGSgage(site_no)
            gage.fetch_site_info()

            status_text.text(f"Downloading data for {site_no}...")

            # Download entire period of record using POR dates from site info
            # USGS API returns limited data without date range, so use POR dates
            start_dt = gage.daily_por_start if gage.daily_por_start else "1900-01-01"
            end_dt = gage.daily_por_end if gage.daily_por_end else None
            daily_data = gage.download_daily_flow(start_date=start_dt, end_date=end_dt)

            # Store gage info
            site_name_display = gage.site_name or "Unknown"
            por_start_date = daily_data.index.min()
            por_end_date = daily_data.index.max()

            st.session_state.gage_info[site_no] = {
                "name": site_name_display,
                "drainage_area": gage.drainage_area,
                "por_start": por_start_date,
                "por_end": por_end_date,
                "records": len(daily_data),
            }

            # Store full daily data
            st.session_state.gage_data[site_no] = daily_data

            # Download peak flows and run FFA if enabled
            if enable_ffa:
                with st.spinner(f"Downloading peak flows for {site_no}..."):
                    try:
                        peak_df = gage.download_peak_flow()
                        st.session_state.peak_data[site_no] = peak_df
                    except Exception as e:
                        st.warning(f"Could not download peak flows for {site_no}: {e}")
                        st.session_state.peak_data[site_no] = None

                if st.session_state.peak_data.get(site_no) is not None:
                    peak_df = st.session_state.peak_data[site_no]
                    with st.spinner(f"Running flood frequency analysis for {site_no}..."):
                        ffa_result = run_ffa(
                            peak_flows=peak_df["peak_flow_cfs"].values,
                            water_years=peak_df["water_year"].values.astype(int),
                            regional_skew=regional_skew,
                            regional_skew_se=regional_skew_se,
                            low_outlier_threshold_override=pilf_override or None,
                        )
                        st.session_state.ffa_results[site_no] = ffa_result
                    if ffa_result.get("error"):
                        st.warning(f"FFA error for {site_no}: {ffa_result['error']}")

        except Exception as e:
            st.error(f"Error processing {site_no}: {str(e)}")

        progress_bar.progress((idx + 1) / len(gage_list))

    status_text.text("Download complete!")
    plt.close("all")
    st.rerun()  # Refresh page to display plots with new data


# Show plot date range controls only if data is loaded
if st.session_state.gage_data:
    with date_range_container:
        st.header("Date Range Filters")
        st.markdown("*Filter plots without re-downloading*")

        # Daily flow date range
        st.markdown("**Daily Flow Range**")
        all_starts = [info["por_start"] for info in st.session_state.gage_info.values()]
        all_ends = [info["por_end"] for info in st.session_state.gage_info.values()]
        overall_start = min(all_starts)
        overall_end = max(all_ends)

        col1, col2 = st.columns(2)
        plot_start = col1.date_input(
            "Start",
            value=overall_start,
            min_value=overall_start,
            max_value=overall_end,
            key="plot_start",
        )
        plot_end = col2.date_input(
            "End",
            value=overall_end,
            min_value=overall_start,
            max_value=overall_end,
            key="plot_end",
        )

        # Peak flow year range (if FFA enabled and peak data exists)
        if enable_ffa and st.session_state.peak_data:
            st.markdown("**Peak Flow Years (for FFA)**")
            st.caption("Click 'Update Plots' to re-run FFA with new range")
            all_peak_years = []
            for site_no, peak_df in st.session_state.peak_data.items():
                if peak_df is not None and len(peak_df) > 0:
                    all_peak_years.extend(peak_df["water_year"].tolist())
            if all_peak_years:
                min_peak_year = int(min(all_peak_years))
                max_peak_year = int(max(all_peak_years))
                col3, col4 = st.columns(2)
                peak_start_year = col3.number_input(
                    "Start Year",
                    min_value=min_peak_year,
                    max_value=max_peak_year,
                    value=min_peak_year,
                    key="peak_start_year",
                )
                peak_end_year = col4.number_input(
                    "End Year",
                    min_value=min_peak_year,
                    max_value=max_peak_year,
                    value=max_peak_year,
                    key="peak_end_year",
                )

        # Update plots button
        if st.button("Update Plots", type="secondary"):
            st.session_state.figures = {}
            # Re-run FFA when anything it depends on changed. The PILF
            # override belongs in here: changing it and not refitting would
            # leave the displayed curve silently stale, which is worse than
            # not offering the control at all.
            current_ffa_inputs = (
                peak_start_year,
                peak_end_year,
                pilf_override,
                regional_skew,
                regional_skew_se,
            )
            if st.session_state.prev_ffa_inputs != current_ffa_inputs:
                st.session_state.prev_ffa_inputs = current_ffa_inputs
                # Re-run FFA with filtered peak data
                if enable_ffa:
                    for site_no in st.session_state.peak_data.keys():
                        peak_df = st.session_state.peak_data[site_no]
                        if peak_df is not None and len(peak_df) > 0:
                            # Filter by year range
                            if peak_start_year is not None and peak_end_year is not None:
                                filtered_peak = peak_df[
                                    (peak_df["water_year"] >= peak_start_year)
                                    & (peak_df["water_year"] <= peak_end_year)
                                ]
                            else:
                                filtered_peak = peak_df
                            if len(filtered_peak) > 0:
                                ffa_result = run_ffa(
                                    peak_flows=filtered_peak["peak_flow_cfs"].values,
                                    water_years=filtered_peak["water_year"].values.astype(int),
                                    regional_skew=regional_skew,
                                    regional_skew_se=regional_skew_se,
                                    low_outlier_threshold_override=pilf_override or None,
                                )
                                st.session_state.ffa_results[site_no] = ffa_result

    # Generate/display plots for each gage
    for site_no in st.session_state.gage_data.keys():
        gage_info = st.session_state.gage_info[site_no]
        full_data = st.session_state.gage_data[site_no]

        # Filter data for plot date range
        plot_data = full_data[
            (full_data.index >= pd.Timestamp(plot_start))
            & (full_data.index <= pd.Timestamp(plot_end))
        ]

        if len(plot_data) == 0:
            st.warning(f"No data for {site_no} in selected date range")
            continue

        # Display gage header
        st.subheader(f"USGS {site_no} - {gage_info['name']}")

        # Format POR and plot dates
        por_str = f"{format_date(gage_info['por_start'])} - {format_date(gage_info['por_end'])}"
        plot_str = f"{format_date(plot_data.index.min())} - {format_date(plot_data.index.max())}"
        da_str = f"{gage_info['drainage_area']:,.1f} sq mi" if gage_info["drainage_area"] else "N/A"

        # Display info
        info_cols = st.columns(4)
        info_cols[0].markdown(
            f"**Drainage Area**<br><small>{da_str}</small>", unsafe_allow_html=True
        )
        info_cols[1].markdown(f"**POR**<br><small>{por_str}</small>", unsafe_allow_html=True)
        info_cols[2].markdown(
            f"**Plot Range**<br><small>{plot_str}</small>", unsafe_allow_html=True
        )
        info_cols[3].markdown(
            f"**Records**<br><small>{len(plot_data):,} days</small>",
            unsafe_allow_html=True,
        )

        # Generate plots if not cached or need refresh
        if site_no not in st.session_state.figures:
            # Determine if we need to show POR annotation (when plot range differs from POR)
            por_start_str = None
            por_end_str = None
            if (
                plot_data.index.min() > gage_info["por_start"]
                or plot_data.index.max() < gage_info["por_end"]
            ):
                por_start_str = format_date(gage_info["por_start"])
                por_end_str = format_date(gage_info["por_end"])

            st.session_state.figures[site_no] = generate_plots(
                site_no,
                plot_data,
                gage_info,
                por_start_str,
                por_end_str,
                yscale_ts=yscale_timeseries,
                yscale_sh=yscale_summary,
                yscale_fdc=yscale_fdc,
            )

        site_figs = st.session_state.figures[site_no]

        # Determine columns based on selected plots
        plot_keys = [k for k in site_figs.keys() if not k.endswith("_stats")]
        if plot_keys:
            cols = st.columns(len(plot_keys))
            for i, plot_key in enumerate(plot_keys):
                with cols[i]:
                    st.pyplot(site_figs[plot_key])

        # Peak flow time series plot
        if show_peak_timeseries and site_no in st.session_state.peak_data:
            peak_df = st.session_state.peak_data[site_no]
            if peak_df is not None and len(peak_df) > 0:
                # Filter by peak year range if set
                if peak_start_year is not None and peak_end_year is not None:
                    peak_df = peak_df[
                        (peak_df["water_year"] >= peak_start_year)
                        & (peak_df["water_year"] <= peak_end_year)
                    ]
                if len(peak_df) == 0:
                    st.warning(f"No peak data for {site_no} in selected year range")
                else:
                    # Build quantiles dict from FFA results if available
                    quantiles = None
                    if show_quantile_lines and site_no in st.session_state.ffa_results:
                        ffa = st.session_state.ffa_results[site_no]
                        if not ffa.get("error") and "quantile_df" in ffa:
                            qdf = ffa["quantile_df"]
                            quantiles = {
                                float(row["Return Interval (yr)"]): row["Flow (cfs)"]
                                for _, row in qdf.iterrows()
                                if float(row["Return Interval (yr)"]) in show_quantile_lines
                            }
                    # Calculate max RI info if enabled
                    max_ri_info = None
                    if show_max_ri and site_no in st.session_state.ffa_results:
                        ffa = st.session_state.ffa_results[site_no]
                        if not ffa.get("error") and "parameters" in ffa:
                            params = ffa["parameters"]
                            max_idx = peak_df["peak_flow_cfs"].idxmax()
                            max_flow = peak_df.loc[max_idx, "peak_flow_cfs"]
                            max_year = int(peak_df.loc[max_idx, "water_year"])
                            ri = estimate_ri_from_lp3(
                                max_flow,
                                params["mean_log"],
                                params["std_log"],
                                params["skew_weighted"],
                            )
                            if ri:
                                max_ri_info = {"flow": max_flow, "year": max_year, "ri": ri}

                    # PILF cut actually applied to the fit, so the override is
                    # visible on the record and not only in the parameter table
                    pilf_threshold = None
                    pilf_source = None
                    if site_no in st.session_state.ffa_results:
                        ffa = st.session_state.ffa_results[site_no]
                        if not ffa.get("error"):
                            fp = ffa.get("parameters", {})
                            if fp.get("n_low_outliers", 0) > 0:
                                pilf_threshold = fp.get("low_outlier_threshold")
                                pilf_source = fp.get("low_outlier_source")

                    peak_fig = plot_peak_timeseries(
                        peak_df,
                        gage_info.get("name", ""),
                        site_no,
                        yscale=yscale_peak,
                        quantiles=quantiles,
                        max_ri_info=max_ri_info,
                        pilf_threshold=pilf_threshold,
                        pilf_source=pilf_source,
                    )
                    st.pyplot(peak_fig)

        # Frequency curve plot
        if show_freq_curve and site_no in st.session_state.ffa_results:
            ffa_result = st.session_state.ffa_results[site_no]
            if not ffa_result.get("error") and ffa_result.get("b17c"):
                selected_skew_labels = [
                    lbl
                    for lbl, on in [
                        ("Station Skew", skew_station_on),
                        ("Weighted Skew", skew_weighted_on),
                        ("Regional Skew", skew_regional_on),
                    ]
                    if on
                ]
                skew_curves = build_skew_curves_dict(ffa_result, selected_skew_labels) or None
                # Build max flow label if enabled
                max_flow_label = None
                if label_max_on_freq and site_no in st.session_state.peak_data:
                    peak_df_freq = st.session_state.peak_data[site_no]
                    if peak_df_freq is not None and len(peak_df_freq) > 0:
                        # Filter by peak year range if set
                        if peak_start_year is not None and peak_end_year is not None:
                            peak_df_freq = peak_df_freq[
                                (peak_df_freq["water_year"] >= peak_start_year)
                                & (peak_df_freq["water_year"] <= peak_end_year)
                            ]
                        if len(peak_df_freq) > 0:
                            params = ffa_result.get("parameters", {})
                            max_idx = peak_df_freq["peak_flow_cfs"].idxmax()
                            max_flow = float(peak_df_freq.loc[max_idx, "peak_flow_cfs"])
                            max_year = int(peak_df_freq.loc[max_idx, "water_year"])
                            ri = (
                                estimate_ri_from_lp3(
                                    max_flow,
                                    params.get("mean_log", 0),
                                    params.get("std_log", 1),
                                    params.get("skew_weighted", 0),
                                )
                                if params
                                else None
                            )
                            max_flow_label = {"flow": max_flow, "year": max_year, "ri": ri}
                freq_fig = plot_frequency_curve_streamlit(
                    ffa_result["b17c"],
                    site_name=gage_info.get("name", ""),
                    site_no=site_no,
                    skew_curves=skew_curves,
                    yscale=yscale_freq,
                    max_flow_label=max_flow_label,
                )
                st.pyplot(freq_fig)

        # FFA results expander
        if enable_ffa and site_no in st.session_state.ffa_results:
            ffa_result = st.session_state.ffa_results[site_no]
            if not ffa_result.get("error"):
                selected_skew_labels = [
                    lbl
                    for lbl, on in [
                        ("Station Skew", skew_station_on),
                        ("Weighted Skew", skew_weighted_on),
                        ("Regional Skew", skew_regional_on),
                    ]
                    if on
                ]
                skew_tables = compute_skew_tables(ffa_result, selected_skew_labels)

                with st.expander("Flood Frequency Results", expanded=False):
                    # Show analysis year range
                    if peak_start_year is not None and peak_end_year is not None:
                        st.info(f"Analysis Period: {peak_start_year} - {peak_end_year}")

                    # Convergence badge
                    if ffa_result.get("converged"):
                        st.success("EMA Converged")
                    else:
                        st.warning("MOM Fallback (EMA did not converge)")

                    st.markdown("**LP3 Parameters**")
                    st.dataframe(
                        format_parameters_df(ffa_result["parameters"]),
                        use_container_width=True,
                    )

                    if skew_tables:
                        for label, tbl in skew_tables.items():
                            st.markdown(
                                f"**Frequency Table — {label}** (Return intervals 1.5–500 years)"
                            )
                            st.dataframe(
                                format_quantile_df(tbl),
                                use_container_width=True,
                            )
                    else:
                        st.markdown("**Frequency Table** (Return intervals 1.5–500 years)")
                        st.dataframe(
                            format_quantile_df(ffa_result["quantile_df"]),
                            use_container_width=True,
                        )

        st.divider()

    # Multi-gage FFA comparison
    if enable_ffa and len(gage_list) > 1:
        sites_with_ffa = [
            s
            for s in gage_list
            if s in st.session_state.ffa_results
            and not st.session_state.ffa_results[s].get("error")
        ]
        if sites_with_ffa:
            st.subheader("Flood Frequency Comparison")
            rows = []
            for sno in sites_with_ffa:
                r = st.session_state.ffa_results[sno]
                q100 = r["quantile_df"][r["quantile_df"]["Return Interval (yr)"] == 100][
                    "Flow (cfs)"
                ].values
                q100_val = int(q100[0]) if len(q100) > 0 else None
                attrs = st.session_state.gage_info.get(sno, {})
                rows.append(
                    {
                        "Site No": sno,
                        "Site Name": attrs.get("name", sno),
                        "100-yr Flow (cfs)": f"{q100_val:,}" if q100_val else "N/A",
                        "Weighted Skew": f"{r['parameters']['skew_weighted']:.3f}",
                        "Method": r.get("method", "ema").upper(),
                        "Converged": "Yes" if r.get("converged") else "No",
                    }
                )
            st.dataframe(pd.DataFrame(rows), use_container_width=True)

    plt.close("all")


# Export section
if st.session_state.gage_data:
    st.sidebar.header("Export")

    # Select which gages to export
    selected_gages = st.sidebar.multiselect(
        "Select gages to export",
        options=list(st.session_state.gage_data.keys()),
        default=list(st.session_state.gage_data.keys()),
    )

    if selected_gages:
        # Create ZIP file in memory
        zip_buffer = io.BytesIO()

        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
            for site_no in selected_gages:
                gage_info = st.session_state.gage_info[site_no]
                full_data = st.session_state.gage_data[site_no]
                site_figs = st.session_state.figures.get(site_no, {})

                # Export daily flow data CSV
                csv_buffer = io.StringIO()
                export_df = full_data.copy()
                export_df.index.name = "date"
                export_df.to_csv(csv_buffer)
                zf.writestr(f"{site_no}/daily_flow.csv", csv_buffer.getvalue())

                # Export plots as PNG
                for name, fig in site_figs.items():
                    if name.endswith("_stats"):
                        # Export stats as CSV
                        csv_buffer = io.StringIO()
                        fig.to_csv(csv_buffer, index=False)
                        zf.writestr(f"{site_no}/{name}.csv", csv_buffer.getvalue())
                    else:
                        # Save PNG
                        img_buffer = io.BytesIO()
                        fig.savefig(img_buffer, format="png", dpi=300, bbox_inches="tight")
                        img_buffer.seek(0)
                        zf.writestr(f"{site_no}/{name}.png", img_buffer.read())

                # Export summary stats if not already in figures
                if "summary_stats" not in site_figs and show_summary:
                    summary_stats = Hydrograph.get_summary_stats(full_data)
                    csv_buffer = io.StringIO()
                    summary_stats.to_csv(csv_buffer, index=False)
                    zf.writestr(f"{site_no}/summary_stats.csv", csv_buffer.getvalue())

                # Export peak flow data if available
                if site_no in st.session_state.peak_data:
                    peak_df_export = st.session_state.peak_data[site_no]
                    if peak_df_export is not None and len(peak_df_export) > 0:
                        # Filter by peak year range if set
                        if peak_start_year is not None and peak_end_year is not None:
                            peak_df_export = peak_df_export[
                                (peak_df_export["water_year"] >= peak_start_year)
                                & (peak_df_export["water_year"] <= peak_end_year)
                            ]
                        if len(peak_df_export) > 0:
                            csv_buffer = io.StringIO()
                            peak_df_export.to_csv(csv_buffer, index=False)
                            zf.writestr(f"{site_no}/peak_flows.csv", csv_buffer.getvalue())

                # Export FFA results
                if enable_ffa and site_no in st.session_state.ffa_results:
                    ffa_result = st.session_state.ffa_results[site_no]
                    if not ffa_result.get("error"):
                        selected_skew_labels = [
                            lbl
                            for lbl, on in [
                                ("Station Skew", skew_station_on),
                                ("Weighted Skew", skew_weighted_on),
                                ("Regional Skew", skew_regional_on),
                            ]
                            if on
                        ]
                        skew_curves_export = (
                            build_skew_curves_dict(ffa_result, selected_skew_labels) or None
                        )
                        freq_fig_for_export = None
                        if show_freq_curve and ffa_result.get("b17c"):
                            # Build max flow label for export if enabled
                            max_flow_label_export = None
                            if label_max_on_freq and site_no in st.session_state.peak_data:
                                peak_df_for_label = st.session_state.peak_data[site_no]
                                if peak_df_for_label is not None and len(peak_df_for_label) > 0:
                                    # Filter by peak year range
                                    if peak_start_year is not None and peak_end_year is not None:
                                        peak_df_for_label = peak_df_for_label[
                                            (peak_df_for_label["water_year"] >= peak_start_year)
                                            & (peak_df_for_label["water_year"] <= peak_end_year)
                                        ]
                                    if len(peak_df_for_label) > 0:
                                        params = ffa_result.get("parameters", {})
                                        max_idx = peak_df_for_label["peak_flow_cfs"].idxmax()
                                        max_flow = float(
                                            peak_df_for_label.loc[max_idx, "peak_flow_cfs"]
                                        )
                                        max_year = int(peak_df_for_label.loc[max_idx, "water_year"])
                                        ri = (
                                            estimate_ri_from_lp3(
                                                max_flow,
                                                params.get("mean_log", 0),
                                                params.get("std_log", 1),
                                                params.get("skew_weighted", 0),
                                            )
                                            if params
                                            else None
                                        )
                                        max_flow_label_export = {
                                            "flow": max_flow,
                                            "year": max_year,
                                            "ri": ri,
                                        }
                            freq_fig_for_export = plot_frequency_curve_streamlit(
                                ffa_result["b17c"],
                                site_name=gage_info.get("name", ""),
                                site_no=site_no,
                                skew_curves=skew_curves_export,
                                yscale=yscale_freq,
                                max_flow_label=max_flow_label_export,
                            )
                        export_ffa_to_zip(zf, site_no, ffa_result, freq_fig_for_export)

            # Multi-gage FFA comparison CSV
            if enable_ffa and len(selected_gages) > 1:
                site_results_for_export = {}
                for sno in selected_gages:
                    if sno in st.session_state.ffa_results:
                        attrs = st.session_state.gage_info.get(sno, {})
                        site_results_for_export[sno] = {
                            "site_name": attrs.get("name", sno),
                            "drainage_area_sqmi": attrs.get("drainage_area", 0),
                            "ffa_results": st.session_state.ffa_results[sno],
                        }
                if site_results_for_export:
                    export_comparison_csv(zf, site_results_for_export)

        zip_buffer.seek(0)

        st.sidebar.download_button(
            label=f"Download {len(selected_gages)} Gage(s) as ZIP",
            data=zip_buffer,
            file_name=f"usgs_data_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip",
            mime="application/zip",
        )

# Footer
st.sidebar.markdown("---")
st.sidebar.markdown(
    f"""
    **hydrolib** v{__version__}
    [GitHub](https://github.com/pinhead001/hydrolib)
    """
)
