# Why the Washington Peak Region 2 equation under-predicts the Methow

Methow-flow priority item 5 (TODO.md, "Priorities for methow-flow"). The SIR 2016-5118 Peak
Region 2 equation is 1.4-2.5x low at the five Methow gages methow-flow transfers from. The
question was whether a GLS refit could replace methow-flow's residual correction.

Reproduce with `python tools/wa_region2_elevation.py`, which runs offline from
`data/regression_study/wa_2016_5118_r2_inputs.csv` and writes
`data/regression_study/wa_2016_5118_r2_elevation.json`. `--build --cache <dir>` rebuilds the
inputs: StreamStats GageStats characteristics, plus flowfreq EMA quantiles through WY2025 with
their peak years.

## Answer

**Mean basin elevation (`ELEV`) explains about half of the bias at common floods, and nothing
else does.**
- Region 2's design `[1, log A, P, CANOPY]` has no elevation term.
- At the 50% AEP, `ELEV` is the strongest candidate (t = 4.1). It cuts the GLS model error
  variance by 18% and the leave-one-out error by 13%. It takes the Methow geometric-mean ratio
  from 2.08 to 1.46.
- The effect is regional, not tuned to the Methow. Refit without any of the 10 calibration
  gages in the Methow HUC8 (17020008), `ELEV` is still significant (t = 3.4), and the
  out-of-basin equation predicts the Methow at 1.48 instead of 2.24.
- The effect fades with rarer floods: t = 1.8 at the 10% AEP and 0.4 at the 1% AEP. That fits
  snowmelt shaping the annual flood in high basins while rain drives the extremes.

**After `ELEV`, nothing basin-specific remains.**
- A Methow-HUC8 indicator alone is borderline at the 50% AEP (t = 2.2, factor 1.68).
- Next to `ELEV` it drops to t = 0.6 at every AEP.
- The leftover Methow ratios, about 1.3-1.5, are within one Region 2 model-error standard
  deviation: √MEV ≈ 0.27 log units, a factor of 1.9.

**No refit of the published form removes the bias at the 1% AEP.** `FOREST` is significant
there (t = -4.4), but it barely moves the Methow and is missing at four of the five Methow
targets. So for rare floods the published equation is as good as a refit.

## Results

Each Methow ratio comes from a fit that leaves that gage out. Order: Methow above Goat Creek,
Twisp River, Methow at Twisp, Methow at Winthrop, Chewuch at Winthrop. Two of the five are
SIR calibration gages. LOO is the leave-one-out mean |Δ log10 Q| over the calibration gages.

| model | AEP | MEV | LOO | t (added) | Methow ratios (geometric mean) |
|---|---|---|---|---|---|
| base `[1, log A, P, CAN]` | 0.5 | 0.0872 | 0.257 | | 1.72 1.80 2.34 2.26 2.38 (2.08) |
| + ELEV | 0.5 | **0.0716** | **0.225** | **4.11** | 1.13 1.43 1.72 1.49 1.57 (**1.46**) |
| + MINBELEV | 0.5 | 0.0786 | 0.245 | 3.26 | (1.95) |
| + JANMINTMP | 0.5 | 0.0848 | 0.260 | -2.03 | (1.79) |
| + I24H2Y | 0.5 | 0.0845 | 0.253 | 2.18 | (1.99) |
| + METHOW_HUC8 | 0.5 | 0.0831 | 0.247 | 2.19 | (1.32) |
| + ELEV + METHOW_HUC8 | 0.5 | 0.0722 | 0.228 | 3.40, 0.58 | (1.33) |
| base | 0.1 | 0.0814 | 0.235 | | (1.72) |
| + ELEV | 0.1 | 0.0788 | 0.239 | 1.81 | 1.11 1.41 1.67 1.39 1.82 (1.46) |
| base | 0.01 | 0.1128 | 0.301 | | (1.44) |
| + ELEV | 0.01 | 0.1127 | 0.309 | 0.44 | 0.99 1.38 1.47 1.27 1.79 (1.35) |
| + FOREST | 0.01 | 0.0879 | 0.275 | -4.39 | (FOREST missing at 4 of 5) |

Base-model residuals at the 50% AEP, by distance from the Methow:
- within 40 miles: +0.22 log units, 87% positive (15 gages);
- 40-80 miles: -0.04;
- 80-150 miles: -0.03.

Nine of the ten Methow-HUC8 calibration gages lie above the published equation. This is the
cluster that `ELEV` absorbs.

## Recommendation for methow-flow

- **50% and 20% AEPs:** use a Region 2 refit with `ELEV`. It is a physical, regionally
  supported term that explains the Methow bias without fitting to it. It replaces the residual
  line for these floods, and it carries its own standard error.
- **10% AEP and rarer:** `ELEV` is not significant. Keep the published equation, and report
  the remaining ratio (about 1.3-1.5) as within the model error, not as a correction to apply.
- A Methow-specific adjustment is not supported once `ELEV` is in the model.

The refit coefficients for each model and AEP are in the results JSON (`beta`, in the order
`[1, log10 A, P, CANOPY, added...]`, with `ELEV` in feet). They are a flowfreq analysis, not a
USGS product. Storing them as a provisional equation set would follow the Oregon pattern
(`provisional/`); that is not done here.

## A defect this found in `flowfreq.regression.gls`

The first pass computed concurrent years from each record's first-to-last span, while `n` was
the number of peaks. For a record with gaps, that gives effective correlations above one, and
a sampling covariance Σ with ten negative eigenvalues. `fit_gls` then estimated the model error
variance at γ² = 0, where the moment equation uses the inverse of an indefinite matrix. That
returned 0, and nonsense coefficients: an intercept of 4.4 instead of 1.0 in one leave-one-out
fold.

Two fixes follow:
- `sampling_covariance` caps concurrent years at `min(n_i, n_j)`.
- `fit_gls` searches γ² only where Σ + γ²I is positive definite. The WREG construction does
  not guarantee a positive semi-definite Σ even with correct inputs: the committed Washington
  refit's Σ has a smallest eigenvalue of about -7e-4.

The Washington refit (`docs/WA_REGION4_REFIT.md`) moves by at most 0.3% in any coefficient
and 0.14% in MEV under both fixes, and its Region 4 conclusion is unchanged.
