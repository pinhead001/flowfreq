# Vignette: Jupyter Notebook Analysis

This vignette walks through a complete Bulletin 17C flood frequency analysis in a Jupyter notebook, using Big Sandy River at Bruceton, TN (USGS 03606500) as the demonstration gage.

See [`examples/flood_frequency_analysis.ipynb`](../examples/flood_frequency_analysis.ipynb) for the runnable notebook.

## Setup

```bash
pip install -e ".[dev]"
pip install jupyter
jupyter notebook
# or: jupyter lab
```

---

## Cell 1 — Imports

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from hydrolib import USGSgage, Bulletin17C, Hydrograph

%matplotlib inline
plt.rcParams["figure.dpi"] = 120
```

---

## Cell 2 — Download Data

```python
SITE_NO = "03606500"   # Big Sandy River at Bruceton, TN

gage = USGSgage(SITE_NO)
gage.fetch_site_info()

print(f"Site:           {gage.site_name}")
print(f"Drainage area:  {gage.drainage_area} sq mi")

# Annual peak flows
peak_df = gage.download_peak_flow()
print(f"\nPeak flow record: {len(peak_df)} years")
print(peak_df.tail())

# Daily flow (for hydrograph plots)
daily_df = gage.download_daily_flow()
print(f"\nDaily flow record: {len(daily_df):,} days")
```

---

## Cell 3 — Hydrograph Plots

```python
fig1 = Hydrograph.plot_daily_timeseries(
    daily_df,
    site_name=gage.site_name,
    site_no=SITE_NO,
)

fig2 = Hydrograph.plot_summary_hydrograph(
    daily_df,
    site_name=gage.site_name,
    site_no=SITE_NO,
    percentiles=[10, 25, 50, 75, 90],
)

fig3, fdc_stats = Hydrograph.plot_flow_duration_curve(
    daily_df,
    site_name=gage.site_name,
    site_no=SITE_NO,
)
plt.show()
print(fdc_stats)
```

---

## Cell 4 — Run Bulletin 17C EMA

```python
# Regional skew: nationwide B17C default
REGIONAL_SKEW    = -0.302
REGIONAL_SKEW_SE =  0.55      # standard error
REGIONAL_SKEW_MSE = REGIONAL_SKEW_SE ** 2

b17c = Bulletin17C(
    peak_flows=peak_df["peak_flow_cfs"].values,
    water_years=peak_df["water_year"].values.astype(int),
    regional_skew=REGIONAL_SKEW,
    regional_skew_mse=REGIONAL_SKEW_MSE,
)

results = b17c.run_analysis(method="ema")

print(f"Method:           {results.method.name}")
print(f"n peaks:          {results.n_peaks}")
print(f"EMA converged:    {results.ema_converged} ({results.ema_iterations} iterations)")
print(f"\nμ (log10 Q):      {results.mean_log:.4f}")
print(f"σ (log10 Q):      {results.std_log:.4f}")
print(f"Station skew:     {results.skew_station:.4f}")
print(f"Regional skew:    {results.skew_weighted:.4f}  (weighted)")
print(f"Skew used:        {results.skew_used:.4f}")
print(f"\nLow outlier threshold: {results.low_outlier_threshold:,.0f} cfs")
print(f"Low outliers (MGBT):   {results.n_low_outliers}")
```

---

## Cell 5 — Frequency Table

```python
AEP = np.array([0.667, 0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002])
RI  = [1.5, 2, 5, 10, 25, 50, 100, 200, 500]

q_df  = b17c.compute_quantiles(aep=AEP)
ci_df = b17c.compute_confidence_limits(aep=AEP, confidence=0.90)

freq_table = pd.DataFrame({
    "Return Interval (yr)": RI,
    "AEP (%)":              [f"{p*100:.1f}%" for p in AEP],
    "Flow (cfs)":           q_df["flow_cfs"].round(0).astype(int),
    "Lower 90% CI":         ci_df["lower_5pct"].round(0).astype(int),
    "Upper 90% CI":         ci_df["upper_5pct"].round(0).astype(int),
})
freq_table = freq_table.set_index("Return Interval (yr)")
print(freq_table.to_string())
```

Expected output (approximately):
```
                   AEP (%)  Flow (cfs)  Lower 90% CI  Upper 90% CI
Return Interval (yr)
1.5                  66.7%        8621          7241          9931
2                    50.0%       12891         11034         14750
5                    20.0%       24843         21403         28284
10                   10.0%       36018         30105         42032
25                    4.0%       54736         43920         66552
50                    2.0%       71914         55968         88860
100                   1.0%       92701         70182        115220
200                   0.5%      118111         87283        148939
500                   0.2%      160093        113953        206233
```

---

## Cell 6 — Compare Skew Options

```python
from hydrolib.core import kfactor_array

skew_options = {
    "Station Skew":  results.skew_station,
    "Weighted Skew": results.skew_weighted,
    "Regional Skew": REGIONAL_SKEW,
}

rows = []
for label, skew_val in skew_options.items():
    K = kfactor_array(skew_val, AEP)
    flows = 10 ** (results.mean_log + K * results.std_log)
    rows.append(pd.Series(
        {ri: int(f) for ri, f in zip(RI, flows)},
        name=label,
    ))

comparison = pd.DataFrame(rows)
comparison.columns = [f"RI={ri}yr" for ri in RI]
print(comparison.to_string())
```

---

## Cell 7 — Frequency Curve Plot

```python
from hydrolib.freq_plot import plot_frequency_curve_streamlit

# Single weighted skew curve (default)
fig = plot_frequency_curve_streamlit(
    b17c,
    site_name=gage.site_name,
    site_no=SITE_NO,
)
plt.show()
```

```python
# Multi-skew overlay
skew_curves = {
    "Station Skew":  results.skew_station,
    "Weighted Skew": results.skew_weighted,
    "Regional Skew": REGIONAL_SKEW,
}
fig = plot_frequency_curve_streamlit(
    b17c,
    site_name=gage.site_name,
    site_no=SITE_NO,
    skew_curves=skew_curves,
)
fig.savefig("frequency_curve_comparison.png", dpi=300, bbox_inches="tight")
plt.show()
```

---

## Cell 8 — Save Report

```python
from hydrolib import HydroReport
import os

output_dir = "./output"
os.makedirs(output_dir, exist_ok=True)

report = HydroReport(gage, b17c)
figures = report.generate_all_figures(output_dir)
report.save_report(os.path.join(output_dir, "flood_frequency_report.md"))

print("Saved files:")
for f in os.listdir(output_dir):
    print(f"  {f}")
```

---

## Cell 9 — Batch Analysis (multiple gages)

```python
from hydrolib.usgs import fetch_nwis_batch

GAGES = ["03606500", "09355500", "11066460"]

for site_no in GAGES:
    g = USGSgage(site_no)
    peak_df = g.download_peak_flow()

    b = Bulletin17C(
        peak_flows=peak_df["peak_flow_cfs"].values,
        water_years=peak_df["water_year"].values.astype(int),
        regional_skew=-0.302,
        regional_skew_mse=0.302,
    )
    r = b.run_analysis(method="ema")

    q100 = b.compute_quantiles(aep=np.array([0.01]))["flow_cfs"].values[0]
    print(f"{site_no}  n={r.n_peaks:3d}  γ_w={r.skew_weighted:.3f}  Q100={q100:>10,.0f} cfs")
```
