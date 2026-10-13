# Provisional Montana regional skew, v2 (B-WLS/B-GLS)

> **This is not a USGS study.** It follows published USGS procedure (SIR 2016-5083 app. B, the
> Pacific Northwest study, and SIR 2014-5211 app. 1, Arizona), using flowfreq's own tooling.
> Its table row is `provisional` and is returned only to a caller who passes
> `allow_provisional=True`. Bulletin 17C (p. 31) directs users without a B-GLS study to the
> USGS, and a published USGS study supersedes this one. USGS Montana currently applies the
> Bulletin 17B Plate I map statewide with SE 0.55 (SIR 2025-5019 pp. 9-10).

Northwest studies plan, tasks 5-7 (`docs/superpowers/plans/2026-10-11-northwest-regional-studies.md`).
v1 (`docs/MONTANA_REGIONAL_SKEW_PROVISIONAL.md`) stays as built. v2 closes the gaps listed in
v1's Limitations section.

Version `montana-provisional-v2`, data through water year 2025, built 2026-10-12.

## Result

| Model | n | β (SD) | σ²δ | AVPnew | ERL (yr) | pseudo-R² | added-term t |
|---|---|---|---|---|---|---|---|
| CONSTANT (all) | 300 | -0.019 (0.081) | 0.338 | 0.345 | 24.2 | 0 | |
| **PRECIP (selected)** | 300 | -0.532 (0.146), +0.227 (0.052) per 10 in | 0.306 | **0.317** | 25.8 | 0.096 | 4.39 |
| ELEV (kft) | 300 | -0.497, +0.099 (0.034) | 0.319 | 0.330 | 25.0 | 0.058 | 2.95 |
| FRONT (0 < Divide distance ≤ 60 mi) | 300 | -0.099, +0.330 (0.124) | 0.329 | 0.339 | 24.5 | 0.027 | 2.66 |
| DIVIDE (per 100 mi, + east) | 300 | 0.137, -0.144 (0.055) | 0.336 | 0.348 | 24.0 | 0.007 | 2.61 |
| EAST indicator | 300 | 0.059, -0.097 (0.152) | 0.343 | 0.353 | 23.7 | 0 | 0.64 |
| FRONT + ELEV | 300 | | 0.317 | 0.330 | 25.0 | 0.064 | 1.12, 2.26 |
| **CONSTANT, east of Divide (table row)** | 240 | **-0.039 (0.089)** | 0.377 | **0.385** | **22.2** | 0 | |
| PRECIP, east of Divide | 240 | -0.694 (0.165), +0.342 (0.069) per 10 in | 0.323 | 0.335 | 24.7 | 0.144 | 4.98 |

**How the model was chosen.** The selection rule was fixed in the plan beforehand: the smallest
AVPnew among models whose added terms all have |β/SD| ≥ 2. Under it, mean annual
precipitation wins. Wetter basins, mostly the Rocky Mountain Front and Glacier, have higher
skews. That is the same rain-on-snowmelt structure v1's influential gages pointed to. East of
the Divide:

    G = -0.694 + 0.342 * (PRECIP / 10),   PRECIP in inches (StreamStats PRECIP),   MSE 0.335

**What the table stores.** A covariate model cannot be represented in `regional_skew.csv`, so
the table row is the **east-of-Divide CONSTANT**:
- "Montana east of the Divide (provisional)": **G = -0.04, MSE 0.385, ERL 22 yr**;
- member HUC8s are the 82 in regions 10 and 09 that touch Montana;
- valid in MT only.

West of the Divide (1701), the verified Pacific Northwest study keeps precedence.

**Read the MSE before using it.** 0.385 is *larger* than the 0.302 Bulletin 17B claims for the
Plate I map. Montana skews vary far more than record length explains: σ²δ ≈ 0.34 against the
Pacific Northwest's 0.17. So a regional value earns little weight against a 40-60 year record.
Even the best model, PRECIP, carries about 25 years of information. This result confirms v1's
finding with a better-specified study; it does not overturn it.

## What v2 changed from v1

| v1 limitation | v2 |
|---|---|
| Montana gages only | WY and ID gages in the shared HUC4s 1002, 1007, 1008 and 1701 (20 WY, 9 ID used) |
| Regulation screened by peak codes only | Plus the packaged GAGES-II screen: a `regulated`/`urban` gage is out unless the record was cut at a coded regulation start (14 excluded) |
| Gage locations for redundancy and correlation distances | StreamStats basin centroids (area-weighted centroid of the delineated polygon) |
| Few explanatory variables (gage altitude, area, latitude) | StreamStats `ELEV` and `PRECIP`; distance to the Continental Divide; a Front indicator |
| Station skews not cross-checked against PeakFQ | `tools/check_montana_station_skews.py`, on the SIR 2025-5019 data release's 63 analyses: median \|Δ station skew\| 0.0003, all 63 within 0.01. 06163400.00 differed (3.0 against 0.216) until #135: flowfreq upper-cased lowercase `c` peak codes into code C and dropped 39 peaks, where peakfq's `siteQT` keeps them |

The Continental Divide is the WBD boundary between HUC2 region 17 and regions 10/09: 702
boundary vertices in the study box. Distance is great-circle miles from the basin centroid,
positive to the east.

## Data and screening

`tools/build_montana_skew_study_v2.py` writes the inputs to `data/skew_study/montana_v2_inputs.csv`
and the results to `montana_v2_results.json`.

Of the 1,262 peak-flow sites, these were excluded:

| Reason | Sites |
|---|---|
| Fewer than 25 pre-regulation peaks | 874 |
| Redundant (nested, DAR ≤ 5, standardized distance ≤ 0.5) | 48 |
| No StreamStats basin | 24 |
| GAGES-II regulated or urban, no coded regulation start | 14 |
| No HUC or no drainage area | 2 |

**300 used:** MT 271, WY 20, ID 9; Missouri 237, Columbia 60, Hudson Bay 3. Station skews have
a mean of -0.02, a median of -0.08 and an SD of 0.68. The median pseudo record length is
44 years.

The cross-correlation model was fitted to 2,759 pairs with at least 50 concurrent years:
Z = 0.155 + exp(-0.153 - 0.0087 D).

The 24 sites StreamStats could not delineate are:
- large mainstems: the Missouri at Fort Benton, the Yellowstone below Forsyth, and the Milk;
- basins that extend into Canada: the Kootenai in Idaho and its tributaries.

This loses the biggest basins. v1's sensitivity runs found basins over 3,000 mi² barely move
the constant.

Regenerate:

```bash
python tools/build_montana_skew_study_v2.py --cache .cache/mt_skew_v2 inputs   # ~1 h, StreamStats
python tools/build_montana_skew_study_v2.py --cache .cache/mt_skew_v2 fit
```

## Limitations

- **Not reviewed and not USGS.** No site-by-site review of the kind a USGS study gives every
  gage.
- **Mixed populations** are handled only as EMA/MGBT handles them. The Front gages still
  dominate the residuals, and a storm-type split could do better.
- **Historical periods** run from the first to the last recorded peak, with no `.psf`
  historic periods.
- **The covariate model** (PRECIP) is in this document only, not in the lookup.
