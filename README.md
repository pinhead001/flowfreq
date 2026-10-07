# FlowFreq

A Python library for USGS streamflow data retrieval and Bulletin 17C flood frequency analysis.

The Streamlit web application that used to live here is now a separate project:
[pinhead001/flowfreq-app](https://github.com/pinhead001/flowfreq-app), which depends
on this library.

## Features

- **USGS Data Retrieval** — Download mean daily flow, annual peak flow, and
  instantaneous (unit-value, sub-daily) flow and stage for any gage, from the USGS Water
  Data OGC API by default (`backend="nwis-legacy"` for the legacy NWIS services). NWIS peak
  qualification codes are applied as peakfq 8.1.0 applies them
- **PeakFQ inputs** — Read PeakFQ `.psf` specification files and WATSTORE peak files, and
  run them as Bulletin 17C analyses
- **Bulletin 17C Analysis**
  - Expected Moments Algorithm (EMA) — the current USGS standard method
  - Method of Moments (MOM) fallback
  - Weighted regional skew (MSE weighting per B17C Appendix 6)
  - Multiple Grubbs-Beck test (MGBT) for low outlier detection
  - 90% confidence intervals (5%/95% limits)
  - Station / weighted / regional skew comparison
  - Regional skew lookup by site or location where a verified study applies
    (`flowfreq.regional_skew`; Pacific Northwest so far)
  - Native engine, or the vendored USGS peakfq 8.1.0 Fortran as a selectable engine
- **Regional regression** — Offline USGS peak-flow regression equations for WA, OR, ID,
  MT, CO, UT, WY, NM, AZ and NV with prediction intervals (`flowfreq.regression`), and
  StreamStats delineation, basin characteristics and NSS estimates (`flowfreq.streamstats`)
- **Transposition** — Drainage-area-ratio transfer to ungaged sites, QPPQ daily-series
  transfer, and donor screening by basin similarity
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
- **CLI** — `flowfreq validate` and `flowfreq benchmark` for numerical validation;
  `flowfreq compare` to run both engines on a peak CSV
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
pip install git+https://github.com/pinhead001/flowfreq@v0.10.2
```

**Dependencies:** `numpy`, `pandas`, `matplotlib`, `scipy`, `requests`, `click`, `pyarrow`

### USGS Water Data API key (optional)

Peaks, daily and instantaneous values, and site information come from the USGS Water Data
API (`api.waterdata.usgs.gov`). Without a key it allows about **1000 requests per hour per IP
address**; past that it answers HTTP 429. flowfreq retries 429s with backoff, so occasional
use needs nothing, but bulk work -- batch analyses, the gage catalog build, regional-skew
studies -- can hit the limit. A key raises it.

1. **Get a key:** https://api.waterdata.usgs.gov/signup/
2. **Set it** (pick one):

   ```powershell
   # PowerShell -- current session only
   $env:USGS_API_KEY = "your-key-here"

   # PowerShell or cmd -- persistent for your user (open a new terminal afterwards)
   setx USGS_API_KEY "your-key-here"
   ```

   ```bash
   # bash / Git Bash -- current session; add to ~/.bashrc to persist
   export USGS_API_KEY="your-key-here"
   ```

   ```python
   # or in code, overriding the environment variable for this process
   from flowfreq import waterdata
   waterdata.set_api_key("your-key-here")
   ```

3. **Check it is picked up:** `python -c "import os; print(bool(os.environ.get('USGS_API_KEY')))"`.
   flowfreq sends it as the `X-Api-Key` header on every Water Data API request; a rejected key
   raises an error that says so rather than being retried.

Never commit a key. In GitHub Actions, add it as the repository secret `USGS_API_KEY`
(Settings → Secrets and variables → Actions); the weekly live-service workflow reads it.

## Quick Start

### One-liner analysis

```python
from flowfreq import analyze_gage
from flowfreq.regional_skew import regional_skew_for_site

site = "12449500"                         # Methow River at Twisp, WA
skew = regional_skew_for_site(site)       # Pacific Northwest study: G = -0.07, MSE 0.18
# Raises RegionalSkewUnavailable where no verified study applies. B17C has no
# national default, so there you supply a published skew and its MSE yourself.

result = analyze_gage(
    site_no=site,
    regional_skew=skew.skew,
    regional_skew_mse=skew.skew_mse,
    output_dir="./output",
)
# Saves frequency_curve.png, flood_frequency_report.md, etc.
# A gage the regulation screen classes as regulated is refused unless
# allow_regulated=True.
```

### Step-by-step analysis

```python
import numpy as np
from flowfreq import USGSgage, Bulletin17C

# 1. Download data
gage = USGSgage("12449500")
peak_df = gage.download_peak_flow()

# 2. Run EMA. This fits every peak as an exact systematic value. To apply NWIS
# peak codes (historic, less-/greater-than, regulated), use flowfreq.workflow.run_ffa
# or analyze_gage, which apply them by default.
b17c = Bulletin17C(
    peak_flows=peak_df["peak_flow_cfs"].values,
    water_years=peak_df["water_year"].values.astype(int),
    regional_skew=-0.07,          # Pacific Northwest regional skew (SIR 2016-5083 app. B)
    regional_skew_mse=0.18,       # its MSE (SE 0.42)
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
| `flowfreq.usgs` | `USGSgage` — daily/peak/instantaneous flow and stage download |
| `flowfreq.waterdata` / `flowfreq.peak_sources` | Water Data OGC API client (default) and peak-data backends |
| `flowfreq.peak_codes` | NWIS peak qualification codes → Bulletin 17C treatment |
| `flowfreq.psf` / `flowfreq.psf_convert` / `flowfreq.watstore` | PeakFQ `.psf` and WATSTORE readers; `.psf` → `Bulletin17C` |
| `flowfreq.fortran_engine` | The vendored peakfq 8.1.0 Fortran as a selectable engine |
| `flowfreq.regional_skew` / `flowfreq.skew_study` | Regional skew lookup; B-WLS/B-GLS regional skew development |
| `flowfreq.regulation` / `flowfreq.catalog` | Regulation/urbanization screen; national gage catalog |
| `flowfreq.regression` | Offline regional regression equations and evaluator |
| `flowfreq.streamstats` | StreamStats delineation, basin characteristics, NSS estimates |
| `flowfreq.transpose` / `flowfreq.qppq` / `flowfreq.donor_similarity` | Transposition to ungaged sites, QPPQ, donor screening |
| `flowfreq.future_flow` | Future-condition change-factor framework |
| `flowfreq.subdaily` | Sub-daily metrics: extreme timing, ramping rates |
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
gage = USGSgage("03606500")          # Big Sandy River at Bruceton, TN
gage.fetch_site_info()              # Populates site_name, drainage_area, POR dates
                                    # (still the legacy NWIS site service)

daily_df = gage.download_daily_flow(start_date="2000-01-01")
# → DataFrame indexed by date, column: flow_cfs

peak_df = gage.download_peak_flow()
# → DataFrame: water_year, peak_date, peak_flow_cfs, qualification_code
# peak_date is a UTC date on the default Water Data API backend.
# backend="nwis-legacy" on any download_* method uses the legacy NWIS services.
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
results.n_low_outliers    # rows censored as low outliers (peakfq's gbnlow)
results.n_mgbt_outliers   # of which, peaks the MGBT (or fixed threshold) flagged
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

# Native vs Fortran engine on a peak CSV (needs the built Fortran extension)
flowfreq compare --peaks peaks.csv --regional-skew -0.07 --regional-skew-se 0.4243
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

**Weighted skew:** the station and regional skews are weighted by inverse MSE.
```
w_regional  = MSE_station / (MSE_station + MSE_regional)
G_weighted  = (1 − w_regional) × G_station + w_regional × G_regional
```
`MSE_station` depends on the method:
- **EMA** (the default) follows peakfq 8.1.0's HWN weighting. The at-site MSE is ADJE, the
  censoring-aware adjustment of Bulletin 17B's `10^(A − B·log10(n/10))`, and the weight also
  carries the Halloween determinant ratio `Wd`. When MGBT finds low outliers, the at-site
  MSE switches to the unadjusted Bulletin 17B formula, as `emafit.f` does.
- **MOM** uses B17C Appendix 4 eq. A4-2:
  `[6n(n−1) / ((n−2)(n+1)(n+3))] × (1 + (6/n)G² + (15/n²)G⁴)`.

**90% confidence intervals:**
- **EMA** uses Cohn and others' EMA variance with peakfq's asymmetric bounds
  (`flowfreq._var_emab`), matching peakfq 8.1.0.
- **MOM** uses the classical approximation:
```
Var(log Q_p) = S² × (1/n) × [1 + K·G + (K²/2)(1 + 0.75G²)] + S² × (∂K/∂G)² × MSE_G
CI = log Q̂_p ± z × √Var
```

**Skew options** — Any combination of station, weighted, and regional skew can be selected to overlay multiple LP3 curves on the frequency plot and produce separate frequency tables for comparison.

## References

England, J.F., Jr., et al., 2019, Guidelines for determining flood flow frequency — Bulletin 17C: U.S. Geological Survey Techniques and Methods, book 4, chap. B5, 148 p. https://doi.org/10.3133/tm4B5

## Version History

| Version | Changes |
|---------|---------|
| **v0.10.2** | Water Data peaks reduced to one annual peak per water year (duplicate time series, secondary peaks); clear error when peak codes remove every peak; `download_peak_flow` sets every site attribute; Wave 1/2 cross-border QA; AZ region 5 (AZ verified); Wave 2 catalog regions (CO/UT/AZ); opt-in NID/NLCD regulation screen; NCHRP 15-61 future-flow procedures; weekly live-service CI. See CHANGELOG.md. |
| **v0.10.1** | `fetch_site_info` on the Water Data API (discharge-only period of record); `pd.NA` peak codes accepted; native EMA 5-10x faster with bit-identical results. See CHANGELOG.md. |
| **v0.10.0** | Water Data OGC API default for instantaneous and daily values, API key and 429 backoff; regulation screen (refuses regulated gages unless overridden); national gage catalog; `n_low_outliers` = peakfq `gbnlow`; native EMA reproduces peakfq 8.1.0 on all 24 WY/MT stations; Wave 2 regression equations (CO/UT/WY/NM/AZ/NV); `regional_skew_at`; B-WLS/B-GLS skew tooling; StreamStats polygon and region selection; paired stage/discharge. See CHANGELOG.md. |
| **v0.9.0** | Water Data OGC API as the default peak backend; NWIS peak codes applied by default; no silent regional skew default (pass a published skew, `station_skew_only=True` or `use_default_skew=True`); `.psf` input; native-EMA fixes toward peakfq 8.1.0 (zero-flow years, per-year perception thresholds, interval peaks, exact LP3 quantiles, B17B skew-MSE switch, near-zero-skew bounds); WA/OR/ID/MT regression equations; Pacific Northwest regional skew. See CHANGELOG.md. |
| **v0.8.0** | `flowfreq.streamstats`: StreamStats delineation and basin characteristics (Phase 1) and NSS flow-statistic estimates (Phase 2); `download_daily_flow` timeout and date-range fixes |
| **v0.7.0** | Transposition of computed flows to ungaged sites (`flowfreq.transpose`: flood, flow-duration and low-flow, each with mandatory exponent provenance); QPPQ daily-series transfer with melt-timing tools (`flowfreq.qppq`); standalone `regime.flow_duration_curve` |
| **v0.6.1** | `plot_peak_flows_with_thresholds` gains a `yscale` toggle, closing the last app/library plot-dedupe gap (verified numerically equivalent to the app's old values) |
| **v0.6.0** | `plot_peak_flows_with_thresholds` gains PILF/MGBT hollow-bar censoring, closing the app/library plot dedupe; mypy fix in `fortran_engine.py`; `make clean` now also wipes `.mypy_cache` |
| **v0.5.0** | Fortran engine as a selectable analysis engine (`engine=`, `compare_engines`, `flowfreq compare` CLI); tests for `hydrograph.py`/`plots.py`/`batch.py`/`cli.py`; return-period lines and max-peak annotation on `plot_peak_flows_with_thresholds` |
| **v0.4.0** | `B17CEngine.fit` station skew corrected to the Bulletin 17C Eq. 7-2 unbiased estimator (numbers move); `plot_frequency_curve_streamlit` renamed to `plot_frequency_curve`; mypy in CI; `make clean-verify` |
| **v0.3.0** | Split into `flowfreq` (library) and `flowfreq-app` (Streamlit app); renamed from `hydrolib`; native `var_mom` port complete (EMA confidence intervals, censored-interval moment iteration, ADJE/`detrat` regional skew weighting); MOM PILF conditional-probability adjustment |
| **v0.2.0** | Instantaneous (unit-value) flow retrieval; low-flow frequency analysis (`flowfreq.lowflow`); flow regime metrics — flashiness, baseflow separation, monthly/seasonal summaries, diel variation (`flowfreq.regime`); parquet flow-series I/O (`flowfreq.flowio`); LP3 negative-skew quantile fix; EMA historical perception-threshold fix |
| **v0.1.0** | Streamlit app with FFA, skew comparison, ZIP export; freq_plot module; ffa_runner/ffa_export app modules; validation framework |
| **v0.0.3** | EMA algorithm, historical flood handling, CLI, validation benchmarks |
| **v0.0.1** | Initial release: MOM, USGS download, hydrograph plots |

## License

MIT License
