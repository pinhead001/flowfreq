# FlowFreq

A Python library for USGS streamflow data retrieval and Bulletin 17C flood frequency analysis.

The Streamlit web application that used to live here is now a separate project:
[pinhead001/flowfreq-app](https://github.com/pinhead001/flowfreq-app), which depends
on this library.

## Features

- **USGS Data Retrieval** — Download mean daily flow, annual peak flow, and
  instantaneous (unit-value, sub-daily) flow from NWIS for any gage
- **Bulletin 17C Analysis**
  - Expected Moments Algorithm (EMA) — the current USGS standard method
  - Method of Moments (MOM) fallback
  - Weighted regional skew (MSE weighting per B17C Appendix 6)
  - Multiple Grubbs-Beck test (MGBT) for low outlier detection
  - 90% confidence intervals (5%/95% limits)
  - Station / weighted / regional skew comparison
- **Low-Flow Frequency Analysis** — Annual minimum *n*-day mean flow (7Q10-style
  statistics), climatic/water/calendar year definitions, LP3 or lognormal
  fit, zero-flow-year handling, analytic and bootstrap confidence intervals
- **Flow Regime Metrics** — Richards-Baker flashiness index, TQmean, baseflow
  separation (UKIH, Lyne-Hollick, and three HYSEP variants), monthly/seasonal
  flow summaries
- **Diel (Sub-Daily) Variation** — Within-day flow range/CV from instantaneous
  data, with timezone-correct local-day grouping
- **Flow Series I/O** — Parquet-backed save/load for daily or instantaneous
  flow series (`flowfreq.flowio`)
- **Hydrograph Plotting** — Daily time series, summary hydrograph, flow duration curve
- **Frequency Curve** — Log-probability axis, LP3 fitted curve, CI band, multi-skew overlay
- **CLI** — `flowfreq validate` and `flowfreq benchmark` for numerical validation
- **Reports** — Automated Markdown technical reports

## Installation

```bash
# Clone or download
git clone https://github.com/pinhead001/flowfreq.git
cd flowfreq

# Create and activate virtual environment (recommended)
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS/Linux:
source .venv/bin/activate

# Install library + dev tools
pip install -e ".[dev]"
```

Or install a released version straight from GitHub:

```bash
pip install git+https://github.com/pinhead001/flowfreq@v0.5.0
```

**Dependencies:** `numpy`, `pandas`, `matplotlib`, `scipy`, `requests`, `click`, `pyarrow`

## Quick Start

### One-liner analysis

```python
from flowfreq import analyze_gage

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
from flowfreq import USGSgage, Bulletin17C

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
| `flowfreq.usgs` | `USGSgage` — NWIS daily/peak/instantaneous data download |
| `flowfreq.bulletin17c` | `Bulletin17C` — EMA/MOM analysis, quantiles, CI, plots |
| `flowfreq.core` | `FrequencyResults`, `LowFlowResults`, `kfactor`, `grubbs_beck_critical_value` |
| `flowfreq.lowflow` | `LowFlowFrequency` — annual *n*-day low-flow frequency analysis |
| `flowfreq.regime` | `FlowRegime` — flashiness, baseflow separation, diel variation, seasonal summaries |
| `flowfreq.flowio` | `save_flow_frame` / `load_flow_frame` — parquet flow-series I/O |
| `flowfreq.hydrograph` | `Hydrograph` — daily timeseries, summary hydrograph, FDC |
| `flowfreq.freq_plot` | `plot_frequency_curve` — frequency curve as a matplotlib `Figure` |
| `flowfreq.report` | `HydroReport` — automated Markdown report |
| `flowfreq.validation` | Benchmark framework, comparison engine |
| `flowfreq.workflow` | `run_ffa` — one call from annual peaks to a fitted frequency curve |

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

### `LowFlowFrequency`

```python
from flowfreq import LowFlowFrequency

lff = LowFlowFrequency(
    daily_df,             # from USGSgage.download_daily_flow()
    n_day=7,               # averaging window, e.g. 7 for 7Q10
    year_type="climatic",  # "climatic" (Apr-Mar, default), "water", or "calendar"
    distribution="lp3",    # or "lognormal"
)
results = lff.run_analysis()

q_df = lff.compute_quantiles(non_exceedance=np.array([0.5, 0.2, 0.1, 0.02]))
# columns: non_exceedance_prob, return_period, flow_cfs, log_flow, K_factor, conditional_prob

ci_df = lff.compute_confidence_limits(non_exceedance=np.array([0.5, 0.1]))
boot_df = lff.compute_bootstrap_confidence_limits(non_exceedance=np.array([0.5, 0.1]))
```

### `FlowRegime`

```python
from flowfreq import FlowRegime

regime = FlowRegime(
    daily_df,
    year_type="water",                   # default; independent of LowFlowFrequency's
    baseflow_method="ih_smoothed_minima",  # or hysep_*/lyne_hollick — see BASEFLOW_METHODS
)

regime.annual     # per-year: flashiness_index, tqmean, baseflow_index, completeness
regime.monthly    # per year/month means, gated on true calendar-day completeness
regime.seasonal   # per year/season means
regime.summary()  # period-of-record pooled summary
```

### Diel Variation

```python
from flowfreq import diel_variation, diel_variation_summary

iv_df = gage.download_instantaneous_flow(tz="America/Los_Angeles")
diel_df = diel_variation(iv_df, tz="America/Los_Angeles")   # tz required, no default
diel_variation_summary(diel_df)
```

See [`docs/vignette_lowflow_regime.md`](docs/vignette_lowflow_regime.md) for a
full walkthrough of all three, with a worked example and method-choice notes.

## CLI

```bash
# Validate EMA against reference fixtures
flowfreq validate

# Benchmark report (text)
flowfreq benchmark

# Benchmark report (JSON)
flowfreq benchmark --format json
```

## Web Application

The interactive Streamlit app lives in its own repository,
[pinhead001/flowfreq-app](https://github.com/pinhead001/flowfreq-app). It installs
this library as a pinned dependency; see that repo's README to run or deploy it.

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
| [Jupyter Notebook](docs/vignette_jupyter.md) | Interactive Bulletin 17C analysis walkthrough |
| [Low-Flow & Flow Regime](docs/vignette_lowflow_regime.md) | Low-flow frequency, flashiness, baseflow separation, diel variation |

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
| **v0.5.0** | Fortran engine as a selectable analysis engine (`engine=`, `compare_engines`, `flowfreq compare` CLI); tests for `hydrograph.py`/`plots.py`/`batch.py`/`cli.py`; return-period lines and max-peak annotation on `plot_peak_flows_with_thresholds` |
| **v0.4.0** | `B17CEngine.fit` station skew corrected to the Bulletin 17C Eq. 7-2 unbiased estimator (numbers move); `plot_frequency_curve_streamlit` renamed to `plot_frequency_curve`; mypy in CI; `make clean-verify` |
| **v0.3.0** | Split into `flowfreq` (library) and `flowfreq-app` (Streamlit app); renamed from `hydrolib`; native `var_mom` port complete (EMA confidence intervals, censored-interval moment iteration, ADJE/`detrat` regional skew weighting); MOM PILF conditional-probability adjustment |
| **v0.2.0** | Instantaneous (unit-value) flow retrieval; low-flow frequency analysis (`flowfreq.lowflow`); flow regime metrics — flashiness, baseflow separation, monthly/seasonal summaries, diel variation (`flowfreq.regime`); parquet flow-series I/O (`flowfreq.flowio`); LP3 negative-skew quantile fix; EMA historical perception-threshold fix |
| **v0.1.0** | Streamlit app with FFA, skew comparison, ZIP export; freq_plot module; ffa_runner/ffa_export app modules; validation framework |
| **v0.0.3** | EMA algorithm, historical flood handling, CLI, validation benchmarks |
| **v0.0.1** | Initial release: MOM, USGS download, hydrograph plots |

## License

MIT License
