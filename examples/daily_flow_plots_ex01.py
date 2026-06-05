"""
Script to download and plot daily flow hydrographs using hydrolib.

This script demonstrates downloading daily streamflow data from USGS
and plotting both a daily time series and summary hydrograph.
"""

import os

import matplotlib
import matplotlib.pyplot as plt

from hydrolib import Hydrograph
from hydrolib.usgs import USGSGage

matplotlib.use("Agg")  # Use non-interactive backend for saving plots

# Output directory
output_dir = r"C:\atemp\NV5\NM_SanJuanR\hydrology\py\sta_09355500"
os.makedirs(output_dir, exist_ok=True)

# Example USGS site (Santa Ana River at MWD Crossing)
site_no = "09355500"

# Initialize USGS gage
gage = USGSGage(site_no)

# Download daily flow data (complete water years for POR)
print(f"Downloading daily flow data for USGS {site_no}...")
daily_data = gage.download_daily_flow(start_date="1970-10-01", end_date="2025-10-01")

print(f"Downloaded {len(daily_data)} days of data")
print(f"Site name: {gage.site_name}")
print(f"Drainage area: {gage.drainage_area} sq mi")

# Save daily data table
daily_data_path = os.path.join(output_dir, "daily_flow_data.csv")
daily_data.to_csv(daily_data_path)
print(f"Saved daily data to: {daily_data_path}")

# Plot daily time series
print("\nPlotting daily time series...")
fig1 = Hydrograph.plot_daily_timeseries(
    daily_data,
    site_name=gage.site_name,
    site_no=gage.site_no,
    save_path=os.path.join(output_dir, "daily_timeseries.png"),
    figsize=(12, 5),
)

# Plot summary hydrograph
print("Plotting summary hydrograph...")
fig2 = Hydrograph.plot_summary_hydrograph(
    daily_data,
    site_name=gage.site_name,
    site_no=gage.site_no,
    save_path=os.path.join(output_dir, "summary_hydrograph.png"),
    figsize=(10, 6),
    percentiles=[10, 25, 50, 75, 90],
)

# Plot flow duration curve
print("Plotting flow duration curve...")
fig3, stats_df = Hydrograph.plot_flow_duration_curve(
    daily_data,
    site_name=gage.site_name,
    site_no=gage.site_no,
    save_path=os.path.join(output_dir, "flow_duration_curve.png"),
    table_path=os.path.join(output_dir, "flow_duration_stats.csv"),
    figsize=(8, 6),
)

print("\nFiles saved to:", output_dir)
print("- daily_flow_data.csv")
print("- daily_timeseries.png")
print("- summary_hydrograph.png")
print("- flow_duration_curve.png")
print("- flow_duration_stats.csv")
