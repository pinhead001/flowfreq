# %% Imports
import os

from hydrolib import B17CEngine, Bulletin17C, PeakRecord, USGSGage, plot_frequency_curve

# Set output directory (outside repo)
output_dir = r"C:\Users\chrisj.nelson\hydro_output\gage_11066460"
os.makedirs(output_dir, exist_ok=True)

# %% 1. Download peak flow data from USGS
site_no = "11066460"
gage = USGSGage(site_no)
gage.download_peak_flow()

print(f"Site: {gage.site_name}")
print(f"Site No: {gage.site_no}")
print(f"Drainage Area: {gage.drainage_area} sq mi")
print(f"Period of Record: {gage.period_of_record}")
print(f"Number of peaks: {len(gage.peak_data)}")

# View the data
print("\nPeak Flow Data:")
print(gage.peak_data)
# print(gage.peak_data.head(10))

# %% 2. Review data and add historical peaks
# Get records as PeakRecord objects for review
records = gage.get_peak_records()
print(f"\nTotal systematic records: {len(records)}")

# Review largest peaks
flows_sorted = sorted(records, key=lambda r: r.flow or 0, reverse=True)
print("\nTop 10 largest peaks:")
for i, r in enumerate(flows_sorted[:10], 1):
    print(f"  {i}. {r.year}: {r.flow:,.0f} cfs")

# Add historical peaks (example - adjust based on local knowledge)
# These are floods known to have occurred before systematic record
historical_peaks = [
    # (year, flow_cfs) - example historical floods
    # Uncomment and modify if you have historical data:
    # (1938, 45000),
    # (1916, 52000),
]

# %% 3. Run Bulletin 17C analysis
# Extract arrays for Bulletin17C
peak_flows = gage.peak_data["peak_flow_cfs"].values
water_years = gage.peak_data["water_year"].values

# Create analysis with optional regional skew
# (Get regional skew from USGS regional skew map for your area)
# regional_skew = -0.50  # Example value - adjust for your region
# regional_skew_mse = 0.14  # Example MSE

regional_skew = None  # Use None to compute station skew only
regional_skew_mse = None  # Use None to compute station skew only

b17c = Bulletin17C(
    peak_flows,
    water_years=water_years,
    regional_skew=regional_skew,
    regional_skew_mse=regional_skew_mse,
    historical_peaks=historical_peaks if historical_peaks else None,
)

# Run EMA analysis
results = b17c.run_analysis(method="ema")

print("\n" + "=" * 50)
print("BULLETIN 17C ANALYSIS RESULTS")
print("=" * 50)
print(f"Method: {results.method.name}")
print(f"Sample size: {results.n_peaks}")
print(f"Systematic records: {results.n_systematic}")
print(f"Historical records: {results.n_historical}")
print(f"Low outliers: {results.n_low_outliers}")
print(f"\nLog-space statistics:")
print(f"  Mean (log10 Q): {results.mean_log:.4f}")
print(f"  Std Dev (log10 Q): {results.std_log:.4f}")
print(f"  Station Skew: {results.skew_station:.4f}")
if results.skew_weighted:
    print(f"  Weighted Skew: {results.skew_weighted:.4f}")
print(f"  Skew Used: {results.skew_used:.4f}")
if results.ema_converged is not None:
    print(f"\nEMA converged: {results.ema_converged} ({results.ema_iterations} iterations)")

# %% 4. Get frequency table
print("\n" + "=" * 50)
print("FLOOD FREQUENCY TABLE")
print("=" * 50)

freq_table = b17c.compute_confidence_limits()
print(freq_table.to_string(index=False))

# %% 5. Get perception thresholds table (HEC-SSP style)
print("\n" + "=" * 50)
print("PERCEPTION THRESHOLDS TABLE")
print("=" * 50)

thresh_table = b17c.get_perception_thresholds_table()
print(thresh_table.to_string(index=False))

# %% 6. Plot frequency curve
fig = b17c.plot_frequency_curve(
    site_name=gage.site_name,
    site_no=gage.site_no,
    save_path=os.path.join(output_dir, "frequency_curve.png"),
    show_confidence=True,
)
print(f"\nFrequency curve saved to: {output_dir}/frequency_curve.png")

# %% 7. Alternative: Use B17CEngine for quick analysis
engine = B17CEngine()
engine.fit(records)

print("\n" + "=" * 50)
print("B17C ENGINE SUMMARY")
print("=" * 50)
print(engine.summary())

# Get formatted frequency table
print("\nFrequency Table with 95% CI:")
print(engine.frequency_table().to_string(index=False))

# %% 8. Generate short report (no tabulate required)

# Get the tables first
freq_table = b17c.compute_confidence_limits()
thresh_table = b17c.get_perception_thresholds_table()

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
- Regional Skew: {regional_skew:.3f}
- Regional Skew MSE: {regional_skew_mse:.3f}
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

print(f"\nReport saved to: {report_path}")


# %%
