# HydroLib

A Python library for USGS streamflow data retrieval and Bulletin 17C flood frequency analysis, with an interactive Streamlit web application.

## Features

- **USGS Data Retrieval** — Download mean daily flow and annual peak flow from NWIS for any gage
- **Bulletin 17C Analysis**
  - Expected Moments Algorithm (EMA) — the current USGS standard method
  - Method of Moments (MOM) fallback
  - Weighted regional skew (MSE weighting per B17C Appendix 6)
  - Multiple Grubbs-Beck test (MGBT) for low outlier detection
  - 90% confidence intervals (5%/95% limits)
  - Station / weighted / regional skew comparison
- **Hydrograph Plotting** — Daily time series, summary hydrograph, flow duration curve
- **Frequency Curve** — Log-probability axis, LP3 fitted curve, CI band, multi-skew overlay
- **Streamlit App** — Interactive web app: download, analyze, plot, and export as ZIP
- **CLI** — `hydrolib validate` and `hydrolib benchmark` for numerical validation
- **Reports** — Automated Markdown technical reports

## Installation

```bash
# Clone or download
git clone https://github.com/pinhead001/hydrolib.git
cd hydrolib

# Create and activate virtual environment (recommended)
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS/Linux:
source .venv/bin/activate

# Install library + dev tools
pip install -e ".[dev]"

# Install Streamlit (for the web app)
pip install streamlit
```

**Dependencies:** `numpy`, `pandas`, `matplotlib`, `scipy`, `requests`, `click`

## Quick Start

### One-liner analysis

```python
from hydrolib import analyze_gage

result = analyze_gage(
    site_no="03606500",          # Big Sandy River at Bruceton TN
    regional_skew=-0.302,        # Nationwide mean from B17C
    regional_skew_mse=0.302,     # = SE² = 0.55²
    output_dir="./output",
)
# Saves frequency_curve.png, flood_frequency_report.md, etc.
```

### Step-by-step analysis

```python
import numpy as np
from hydrolib import USGSgage, Bulletin17C

# 1. Download data
gage = USGSgage("03606500")
peak_df = gage.download_peak_flow()

# 2. Run EMA
b17c = Bulletin17C(
    peak_flows=peak_df["peak_flow_cfs"].values,
    water_years=peak_df["water_year"].values.astype(int),
    regional_skew=-0.302,
    regional_skew_mse=0.302,      # SE² (0.55² ≈ 0.302)
)
results = b17c.run_analysis(method="ema")

# 3. Quantile table
aep = np.array([0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002])
q_df = b17c.compute_quantiles(aep=aep)
ci_df = b17c.compute_confidence_limits(aep=aep)
print(q_df)

# 4. Frequency curve
fig = b17c.plot_frequency_curve(site_name=gage.site_name, site_no=gage.site_no)
fig.savefig("frequency_curve.png", dpi=300, bbox_inches="tight")
```

## Module Overview

| Module | Purpose |
|--------|---------|
| `hydrolib.usgs` | `USGSgage` — NWIS data download |
| `hydrolib.bulletin17c` | `Bulletin17C` — EMA/MOM analysis, quantiles, CI, plots |
| `hydrolib.core` | `FrequencyResults`, `kfactor`, `grubbs_beck_critical_value` |
| `hydrolib.hydrograph` | `Hydrograph` — daily timeseries, summary hydrograph, FDC |
| `hydrolib.freq_plot` | `plot_frequency_curve_streamlit` — Streamlit-ready frequency curve |
| `hydrolib.report` | `HydroReport` — automated Markdown report |
| `hydrolib.validation` | Benchmark framework, comparison engine |
| `app/streamlit_app.py` | Interactive web app |
| `app/ffa_runner.py` | Analysis runner + display formatters |
| `app/ffa_export.py` | ZIP export (PNG, CSV, LP3 params) |

## Key APIs

### `USGSgage`

```python
gage = USGSgage("03606500")
gage.fetch_site_info()              # Populates site_name, drainage_area, POR dates

daily_df = gage.download_daily_flow(start_date="2000-01-01")
# → DataFrame indexed by date, column: flow_cfs

peak_df = gage.download_peak_flow()
# → DataFrame: water_year, peak_date, peak_flow_cfs, qualification_code
```

### `Bulletin17C`

```python
b17c = Bulletin17C(
    peak_flows,          # np.ndarray of annual peaks (cfs)
    water_years,         # np.ndarray of water years
    regional_skew,       # float or None (station-only if None)
    regional_skew_mse,   # float or None (SE²)
)

results = b17c.run_analysis(method="ema")   # or "mom"

# Results attributes
results.mean_log          # μ (log10)
results.std_log           # σ (log10)
results.skew_station      # station skew
results.skew_weighted     # weighted skew (None if no regional)
results.skew_used         # skew used in quantile calculation
results.ema_converged     # bool
results.n_low_outliers    # MGBT-flagged low outliers
results.low_outlier_threshold  # MGBT threshold (cfs)

# Quantile table
q_df = b17c.compute_quantiles(aep=np.array([0.10, 0.02, 0.01]))
# columns: aep, return_period, flow_cfs, log_flow, K_factor

# Confidence limits
ci_df = b17c.compute_confidence_limits(aep=..., confidence=0.90)
# columns: aep, return_period, flow_cfs, lower_5pct, upper_5pct
```

### Standard AEPs (return intervals 1.5–500 yr)

```python
aep = np.array([0.667, 0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002])
# RI =              1.5,  2,   5,   10,   25,   50,  100,  200,   500
```

## CLI

```bash
# Validate EMA against reference fixtures
hydrolib validate

# Benchmark report (text)
hydrolib benchmark

# Benchmark report (JSON)
hydrolib benchmark --format json
```

## Streamlit App

```bash
# Run locally
streamlit run app/streamlit_app.py
# Opens at http://localhost:8501
```

**App features:**
- Single or multi-gage mode
- Daily time series, summary hydrograph, flow duration curve
- Flood Frequency Analysis (enable in sidebar):
  - EMA with MOM fallback
  - Regional skew inputs
  - Skew option checkboxes: Station / Weighted / Regional
  - Frequency curve with multi-skew overlay
  - Per-skew frequency tables (RI 1.5–500 yr, 90% CI)
- ZIP export: plots (PNG), daily flow CSV, frequency table, LP3 parameters
- Multi-gage comparison table

## Vignettes

| Guide | Description |
|-------|-------------|
| [CLI Usage](docs/vignette_cli.md) | Command-line validation and benchmarking |
| [Jupyter Notebook](docs/vignette_jupyter.md) | Interactive analysis walkthrough |
| [Local Streamlit](docs/vignette_streamlit_local.md) | Run the web app on your machine |
| [Streamlit Cloud](docs/vignette_streamlit_web.md) | Deploy to Streamlit Community Cloud |

## Bulletin 17C Technical Notes

**Weighted skew** (per B17C §6):
```
MSE_station = (6/n) × (1 + (6/n)G² + (15/n²)G⁴)
w_regional  = MSE_station / (MSE_station + MSE_regional)
G_weighted  = w_station × G_station + w_regional × G_regional
```

**90% confidence intervals:**
```
Var(X_p) = σ² × [1/n + K² × (1 + 0.75G²) / (2(n-1))]
CI = Q̂ ± 1.645 × σ × √Var_factor
```

**Skew options** — Any combination of station, weighted, and regional skew can be selected to overlay multiple LP3 curves on the frequency plot and produce separate frequency tables for comparison.

## References

England, J.F., Jr., et al., 2019, Guidelines for determining flood flow frequency — Bulletin 17C: U.S. Geological Survey Techniques and Methods, book 4, chap. B5, 148 p. https://doi.org/10.3133/tm4B5

## Version History

| Version | Changes |
|---------|---------|
| **v0.1.0** | Streamlit app with FFA, skew comparison, ZIP export; freq_plot module; ffa_runner/ffa_export app modules; validation framework |
| **v0.0.3** | EMA algorithm, historical flood handling, CLI, validation benchmarks |
| **v0.0.1** | Initial release: MOM, USGS download, hydrograph plots |

## License

MIT License
