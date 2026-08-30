# Vignette: Low-Flow, Flow-Regime, and Sub-Daily Variation Analysis

This vignette walks through hydrolib's low-flow frequency analysis, flow-regime
metrics (flashiness, baseflow separation, seasonal summaries), and sub-daily
(diel) variation analysis, using the Methow River at Pateros, WA (USGS
12449950) as the example gage.

It complements [`vignette_jupyter.md`](vignette_jupyter.md), which covers
Bulletin 17C flood frequency analysis — the modules here answer a different
set of questions: not "how big is the 100-year flood," but "how low does this
stream get, how flashy or baseflow-dominated is it, and how much does it swing
within a single day."

**A note on the data used below.** The download cell uses the real NWIS site
number and the real download API. Everything after it runs on a reproducible,
seeded synthetic daily/instantaneous series built to resemble this basin's
snowmelt-driven seasonality (winter rain bumps, a spring freshet, a summer/fall
recession), so that every number printed in this vignette can be reproduced
exactly without a network call. Replace the synthetic series with the
downloaded `daily_df` / `iv_df` to analyze the actual gage — every function
call below is unchanged either way, since both are plain DataFrames with the
same shape.

## Setup

```bash
pip install -e ".[dev]"
```

---

## Cell 1 — Download Data

```python
import numpy as np
import pandas as pd
from hydrolib import USGSgage

SITE_NO = "12449950"   # Methow River at Pateros, WA

gage = USGSgage(SITE_NO)
gage.fetch_site_info()
print(f"Site:           {gage.site_name}")
print(f"Drainage area:  {gage.drainage_area} sq mi")

daily_df = gage.download_daily_flow()
# → DataFrame indexed by date, column: flow_cfs

iv_df = gage.download_instantaneous_flow(tz="America/Los_Angeles")
# → DataFrame indexed by tz-aware datetime, columns: flow_cfs,
#   datetime_local, tz_cd, qualification_code. Instantaneous records are
#   typically much shorter than the daily record (often starting ~2007),
#   and download_instantaneous_flow raises NoInstantaneousDataError if the
#   site or window has none.
```

---

## Cell 2 — Reproducible Demonstration Series

```python
np.random.seed(7)

years = range(1999, 2024)  # 25 water years
dates, flows = [], []
for wy in years:
    idx = pd.date_range(f"{wy - 1}-10-01", f"{wy}-09-30", freq="D")
    frac_year = np.arange(len(idx)) / len(idx)  # 0 = Oct 1

    winter = 250 * np.exp(-0.5 * ((frac_year - 0.28) / 0.10) ** 2)
    freshet = 1400 * np.exp(-0.5 * ((frac_year - 0.62) / 0.09) ** 2)
    base = 90 + 40 * np.cos(2 * np.pi * (frac_year - 0.05))
    year_factor = np.random.uniform(0.75, 1.3)

    seasonal = (base + winter + freshet) * year_factor
    noise = np.random.lognormal(mean=0.0, sigma=0.08, size=len(idx))
    daily_flow = np.clip(seasonal * noise, 5, None)

    for _ in range(2):  # a couple of rain-on-snow winter spikes
        spike_day = np.random.randint(30, 110)
        width = np.random.randint(2, 5)
        daily_flow = daily_flow + 600 * year_factor * np.exp(
            -0.5 * ((np.arange(len(idx)) - spike_day) / width) ** 2
        )

    dates.extend(idx)
    flows.extend(daily_flow)

daily_df = pd.DataFrame(
    {"flow_cfs": flows}, index=pd.DatetimeIndex(dates, name="date")
)
daily_df = daily_df[~daily_df.index.duplicated(keep="first")].sort_index()

print(f"Daily flow record: {len(daily_df):,} days, "
      f"{daily_df.index.min().date()} to {daily_df.index.max().date()}")
print(daily_df.describe().round(1))
```

Output:
```
Daily flow record: 9,131 days, 1998-10-01 to 2023-09-30
       flow_cfs
count    9131.0
mean      486.4
std       422.1
min        79.2
25%       176.7
50%       310.9
75%       678.4
max      2330.1
```

A short instantaneous series for the diel-variation section, five July days at
15-minute resolution with a synthetic afternoon-melt / overnight-recession
signal:

```python
iv_dates = pd.date_range("2023-07-01", "2023-07-05", freq="15min", tz="UTC")
local = iv_dates.tz_convert("America/Los_Angeles")
hour_local = local.hour + local.minute / 60
diel_signal = 45 + 10 * np.sin(2 * np.pi * (hour_local - 8) / 24)
iv_noise = np.random.normal(0, 0.8, size=len(iv_dates))
iv_df = pd.DataFrame({"flow_cfs": diel_signal + iv_noise}, index=iv_dates)
```

---

## Cell 3 — Save / Load a Flow Series

Either `daily_df` above is worth caching locally rather than re-downloading on
every run. `hydrolib.flowio` wraps `daily_df.to_parquet`/`pd.read_parquet`
with a consistent default compression and dtype/timezone handling:

```python
from hydrolib import save_flow_frame, load_flow_frame

save_flow_frame(daily_df, "methow_daily.parquet")
reloaded = load_flow_frame("methow_daily.parquet")
print(reloaded.equals(daily_df))
```

Output:
```
True
```

---

## Cell 4 — Low-Flow Frequency Analysis

`LowFlowFrequency` fits a distribution to the annual minimum *n*-day mean flow
— the low-flow analogue of Bulletin 17C's annual maximum peak analysis. The
default `n_day=7` and `year_type="climatic"` give the familiar 7Q2/7Q10-style
statistics; the climatic year (Apr 1 – Mar 31, per Riggs 1972) keeps a single
low-flow season from being split across a calendar- or water-year boundary.

```python
from hydrolib import LowFlowFrequency

lff = LowFlowFrequency(daily_df, n_day=7, year_type="climatic")
results = lff.run_analysis()

print(f"n years used:  {lff.n}")
print(f"mean_log:      {results.mean_log:.4f}")
print(f"std_log:       {results.std_log:.4f}")
print(f"skew_used:     {results.skew_used:.4f}")
print(f"p_zero:        {results.p_zero:.4f}")   # fraction of years with a zero-flow minimum
```

Output:
```
n years used:  24
mean_log:      2.0263
std_log:       0.0531
skew_used:     0.7971
p_zero:        0.0000
```

`compute_quantiles` returns the standard 7QN table:

```python
q_df = lff.compute_quantiles(non_exceedance=np.array([0.5, 0.2, 0.1, 0.02]))
print(q_df[["non_exceedance_prob", "return_period", "flow_cfs"]].round(1).to_string(index=False))
```

Output (7Q2, 7Q5, 7Q10, 7Q50):
```
 non_exceedance_prob  return_period  flow_cfs
                 0.50            2.0     104.6
                 0.20            5.0      95.7
                 0.10           10.0      92.1
                 0.02           50.0      87.3
```

Two ways to get a confidence interval — the fast analytic approximation, and a
parametric bootstrap that also propagates uncertainty in the zero-flow
probability `p_zero` (see the module docstring for when that gap matters):

```python
ci_df = lff.compute_confidence_limits(non_exceedance=np.array([0.5, 0.1]))
print(ci_df.round(1).to_string(index=False))

boot_df = lff.compute_bootstrap_confidence_limits(
    non_exceedance=np.array([0.5, 0.1]), n_resamples=500, random_state=42
)
print(boot_df.round(1).to_string(index=False))
```

Output:
```
 non_exceedance_prob  return_period  flow_cfs  lower_90pct  upper_90pct
                  0.5            2.0     104.6        100.3        109.0
                  0.1           10.0      92.1         86.9         97.7

 non_exceedance_prob  return_period  flow_cfs  lower_90pct  upper_90pct  n_resamples_used
                  0.5            2.0     104.6        100.6        109.4               500
                  0.1           10.0      92.1         88.7         96.1               500
```

**Method notes** (see the `lowflow` module docstring for the full rationale):
- `distribution="lognormal"` is available as an alternative to the default
  `"lp3"` when too few nonzero years are available to trust a sample skew
  estimate (`LowFlowFrequency` warns when this applies).
- A record with any zero-flow minima uses a conditional-probability
  formulation (fit the distribution to the nonzero years, then rescale by
  `p_zero`) rather than silently dropping or log-transforming zeros.
- `annual_minimum_flow(daily_df, n_day=7, year_type="climatic")` is available
  standalone if you just want the annual-minimum table without fitting a
  distribution to it.

---

## Cell 5 — Flow Regime Metrics

`FlowRegime` bundles flashiness, baseflow separation, and monthly/seasonal
summaries behind one object so they're computed from a single consistent
year/completeness convention (default `year_type="water"`, Oct 1 – Sep 30 —
deliberately different from `LowFlowFrequency`'s climatic-year default, since
a regime description is naturally organized by water year while a low-flow
extreme should not be split across one).

```python
from hydrolib import FlowRegime

regime = FlowRegime(daily_df, year_type="water", baseflow_method="ih_smoothed_minima")
print(regime.annual.head(3).round(3).to_string(index=False))
```

Output:
```
 year  flashiness_index  tqmean  baseflow_index  n_days  complete  baseflow_complete
 1999             0.092   0.315           0.858     365      True               True
 2000             0.094   0.322           0.837     366      True               True
 2001             0.092   0.321           0.867     365      True               True
```

`flashiness_index` is the Richards-Baker Flashiness Index (day-to-day flow
volatility relative to total flow); `tqmean` is the fraction of the year flow
exceeds the annual mean (a stream with strong sustained snowmelt runs well
below its mean most of the year, so a low `tqmean` here is expected); the two
completeness columns can differ because baseflow separation has its own
edge-effect losses (see `baseflow_index`'s docstring) — never assume one
implies the other.

A period-of-record summary pools these across complete years (volume-weighted
for flashiness and baseflow index, arithmetic mean for `tqmean` — see the
module docstring for why):

```python
print(regime.summary().round(4))
```

Output:
```
n_years             25.0000
flashiness_index     0.0962
tqmean               0.3258
baseflow_index       0.8560
baseflow_n_years    25.0000
dtype: float64
```

Monthly and seasonal summaries are also available (`regime.monthly`,
`regime.seasonal`), each gated on true calendar-day completeness rather than
however many rows happen to be present:

```python
print(regime.monthly.head(3).round(1).to_string(index=False))
```

Output:
```
 year  month  mean_flow_cfs  min_flow_cfs  max_flow_cfs  median_flow_cfs  n_days  days_in_month  complete
 1999      1          277.4         236.3         443.9            270.3      31             31      True
 1999      2          213.3         180.0         259.3            204.0      28             28      True
 1999      3          315.3         196.8         482.5            296.9      31             31      True
```

**Method notes:**
- Five baseflow-separation methods are available in `BASEFLOW_METHODS`:
  `"ih_smoothed_minima"` (the default — UKIH/IH smoothed local minima),
  `"lyne_hollick"` (recursive digital filter), and three HYSEP variants
  (`"hysep_fixed"`, `"hysep_sliding"`, `"hysep_local_minimum"`, which require
  `drainage_area_sqmi=gage.drainage_area`). They are not directly comparable
  to each other — see the module docstring's disclosure on UKIH's
  block-alignment sensitivity and the HYSEP window-width constant before
  reporting a BFI value without naming the method used to get it.
- `richards_baker_flashiness`, `tqmean`, `baseflow_index`, `separate_baseflow`,
  `monthly_flow_summary`, and `seasonal_flow_summary` are all available
  standalone if you only need one metric rather than the full bundle.

---

## Cell 6 — Diel (Sub-Daily) Variation

`diel_variation` summarizes how much flow swings within each local calendar
day from an instantaneous series — descriptive only, since a snowmelt diel
signal and a hydropeaking/diversion signal can look alike from streamflow
data alone (see the function's docstring). `tz` is required, with no default:
grouping by the wrong zone silently fractures each local day across two
UTC-labeled buckets.

```python
from hydrolib import diel_variation, diel_variation_summary

diel_df = diel_variation(iv_df, tz="America/Los_Angeles")
print(diel_df.round(2).to_string(index=False))
```

Output:
```
      date  min_flow_cfs  max_flow_cfs  mean_flow_cfs  std_flow_cfs  range_cfs   cv  n_obs  expected_obs  complete
2023-06-30         36.05         53.03          44.05          5.21      16.98 0.12     28          96.0     False
2023-07-01         34.49         55.35          45.04          7.08      20.86 0.16     96          96.0      True
2023-07-02         33.43         56.10          45.03          7.25      22.68 0.16     96          96.0      True
2023-07-03         33.37         55.95          45.02          7.25      22.58 0.16     96          96.0      True
2023-07-04         33.33         55.94          45.58          7.98      22.61 0.17     69          96.0     False
```

Note the first and last rows: they're marked `complete=False` because the UTC
series only partially covers those local calendar days at its edges — `range`
and `cv` are still reported for them (this is deliberate, see the function's
Notes), but a period summary should exclude them:

```python
print(diel_variation_summary(diel_df))
```

Output:
```
n_days                      3.000
mean_diel_amplitude_cfs    22.039
mean_diel_cv                0.160
dtype: float64
```

---

## References

Riggs, H.C., 1972, Low-flow investigations: U.S. Geological Survey Techniques
of Water-Resources Investigations, book 4, chap. B1, 18 p.

Richards, R.P., 1990, Measures of flow variability for Great Lakes tributaries.

Institute of Hydrology, 1980, Low flow studies report: Wallingford, UK,
Institute of Hydrology (UKIH baseflow separation method).

Sloto, R.A., and Crouse, M.Y., 1996, HYSEP: A computer program for streamflow
hydrograph separation and analysis: U.S. Geological Survey Water-Resources
Investigations Report 96-4040.

Lyne, V., and Hollick, M., 1979, Stochastic time-variable rainfall-runoff
modelling: Institute of Engineers Australia National Conference, p. 89-93.
