# %% Complete workflow - run this entire cell

from hydrolib import (
    USGSGage,
    Bulletin17C,
)
import os

# Set output directory
output_dir = r"C:\Users\chrisj.nelson\hydro_output\gage_11066460"
os.makedirs(output_dir, exist_ok=True)

# 1. Download data
site_no = "11066460"
gage = USGSGage(site_no)
gage.download_peak_flow()

print(f"Site: {gage.site_name}")
print(f"Period of Record: {gage.period_of_record}")
print(f"Number of peaks: {len(gage.peak_data)}")

# 2. Set up analysis (station skew only - no regional skew)
peak_flows = gage.peak_data["peak_flow_cfs"].values
water_years = gage.peak_data["water_year"].values
historical_peaks = []  # Empty - no historical data
regional_skew = None   # Station skew only
regional_skew_mse = None

# 3. Run analysis
b17c = Bulletin17C(peak_flows, water_years=water_years)
results = b17c.run_analysis(method='ema')

print(f"\nStation Skew: {results.skew_station:.4f}")
print(f"Skew Used: {results.skew_used:.4f}")

# 4. Get tables
freq_table = b17c.compute_confidence_limits()
thresh_table = b17c.get_perception_thresholds_table()

# 5. Plot
fig = b17c.plot_frequency_curve(
    site_name=gage.site_name,
    site_no=gage.site_no,
    save_path=os.path.join(output_dir, "frequency_curve.png"),
)
print(f"\nPlot saved to: {output_dir}/frequency_curve.png")

# 6. Generate report
weighted_skew_str = f"{results.skew_weighted:.4f}" if results.skew_weighted else "N/A"
drainage_str = f"{gage.drainage_area}" if gage.drainage_area else "N/A"
n_historical = len(historical_peaks) if historical_peaks else 0
regional_skew_str = f"{regional_skew:.3f}" if regional_skew else "N/A (Station skew only)"
regional_mse_str = f"{regional_skew_mse:.3f}" if regional_skew_mse else "N/A"

report_text = f"""
# Flood Frequency Analysis Report
## USGS {gage.site_no} - {gage.site_name}

### Site Information
- Site Number: {gage.site_no}
- Site Name: {gage.site_name}
- Drainage Area: {drainage_str} sq mi
- Period of Record: {gage.period_of_record[0]} - {gage.period_of_record[1]}
- Number of Annual Peaks: {len(gage.peak_data)}

### Analysis Parameters
- Method: Bulletin 17C - {results.method.name}
- Regional Skew: {regional_skew_str}
- Regional Skew MSE: {regional_mse_str}
- Historical Peaks: {n_historical}

### Log-Pearson Type III Statistics
- Mean (log10 Q): {results.mean_log:.4f}
- Std Dev (log10 Q): {results.std_log:.4f}
- Station Skew: {results.skew_station:.4f}
- Weighted Skew: {weighted_skew_str}
- Skew Used: {results.skew_used:.4f}
- Low Outlier Threshold: {results.low_outlier_threshold:,.0f} cfs
- Low Outliers Identified: {results.n_low_outliers}

### Flood Frequency Estimates
{freq_table.to_string(index=False)}

### Perception Thresholds
{thresh_table.to_string(index=False)}

---
*Analysis performed using hydrolib v0.0.3*
"""

report_path = os.path.join(output_dir, "flood_frequency_report.md")
with open(report_path, "w") as f:
    f.write(report_text)

print(f"Report saved to: {report_path}")

# %%
