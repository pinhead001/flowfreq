# Provisional Montana regional skew (B-WLS/B-GLS)

> **This is not a USGS study.** It was developed with flowfreq's own tooling
> (`flowfreq.skew_study`) as a research product. It is **not** a default anywhere in
> flowfreq, it is **not** in `flowfreq/data/regional_skew.csv` (whose schema has only
> `pending`/`verified`, and MT stays `pending`), and nothing looks it up automatically.
> Bulletin 17C (p. 31) directs users without a B-GLS study to consult the USGS; the USGS
> Wyoming-Montana Water Science Center anticipates a study (SIR 2025-5019 pp. 9-10), and a
> published USGS study supersedes this one. Roadmap §1.3 stays open for it.

Version: `montana-provisional-v1` · data through water year 2025 · built 2026-09-30.

## Result

**CONSTANT model, 284 Montana gages: G = -0.03 (SD 0.09), MSE (AVPnew) = 0.36, SE 0.60,
effective record length 23 years.**

**Do not read this as a better skew for Montana.** Its MSE of 0.36 is *larger* than the 0.302
that Bulletin 17B claims for its Plate I map (SE 0.55, the value USGS Montana uses per SIR
2025-5019). The Plate I figure is a nominal national number, not a Montana verification, but
the comparison still stands.

The model error variance, σ²δ = 0.35 (SD 0.04), is twice the Pacific Northwest's 0.17: Montana
station skews vary far more than record length explains. None of the candidate explanatory
variables captures that variation. Gage altitude comes closest, with a slope of +0.14 per
1,000 ft (2.6 SD) that explains only 5% of it (pseudo-R²).

A constant Montana skew is therefore of little use for weighting. It adds about 23 years of
information, which is less than a typical 40-60-year at-site record provides.

| Model | n | β (SD) | σ²δ (SD) | ASEV | AVPnew | ERL (yr) | pseudo-R² | EVR | MBV* |
|---|---|---|---|---|---|---|---|---|---|
| CONSTANT (statewide) | 284 | -0.027 (0.086) | 0.354 (0.042) | 0.0074 | 0.361 | 23.3 | 0 | 0.50 | 3.7 |
| East of Divide only (HUC 09+10) | 238 | -0.035 (0.092) | 0.370 (0.048) | 0.0084 | 0.378 | 22.5 | 0 | 0.47 | 3.5 |
| + east-of-Divide indicator | 284 | 0.013 (0.16), -0.047 (0.16) | 0.356 | 0.011 | 0.367 | 23.1 | 0 | 0.50 | 3.7 |
| + log10 drainage area | 284 | -0.161 (0.12), +0.079 (0.042) | 0.357 | 0.010 | 0.367 | 23.1 | 0 | 0.50 | 3.8 |
| + gage altitude (kft) | 284 | -0.545 (0.22), +0.145 (0.056) | 0.337 | 0.011 | 0.348 | 24.0 | 0.05 | 0.53 | 3.9 |
| + latitude | 284 | -2.48 (2.40), +0.052 (0.051) | 0.355 | 0.011 | 0.366 | 23.1 | 0 | 0.50 | 3.8 |

Sensitivity of the CONSTANT model:

| Variant | n | β (SD) | σ²δ | AVPnew | ERL |
|---|---|---|---|---|---|
| PRL ≥ 35 yr (the PNW study's screen) | 204 | -0.070 (0.088) | 0.326 | 0.334 | 24.8 |
| Drainage area ≤ 3,000 mi² | 270 | -0.021 (0.088) | 0.364 | 0.371 | 22.8 |
| Record extends through 1980 or later | 262 | -0.051 (0.090) | 0.361 | 0.369 | 23.0 |
| κ = 3.3 instead of 2.8 | 284 | -0.027 (0.075) | 0.352 | 0.358 | 23.5 |

West of the Divide (HUC 1701, 46 gages), the published Pacific Northwest study already applies:
G = -0.07, MSE 0.18. This study neither replaces nor contradicts it there. The east-of-Divide
indicator is not significant (-0.05 ± 0.16).

## Method

The procedure is the one SIR 2016-5083 appendix B (Veilleux) used for the Pacific Northwest,
implemented in `flowfreq.skew_study` and validated against that study (below).

1. **Data.** Every Montana discharge peak (`peaks` collection, `state_code=30`) and every
   Montana stream monitoring location (`monitoring-locations`, `site_type_code=ST`) from the
   USGS Water Data OGC API: three paged requests. Montana gages only; neighbouring-state gages
   in shared basins were not added.
2. **Screening.**
   - Regulation, diversion and urbanization evidence: peaks coded 5, 6 or C
     (`peak_codes.REGULATION_EVIDENCE_CODES`, `ALTERATION_EVIDENCE_CODES`). A site with any
     such peak keeps only its record before the first one. This is the "pre-regulation part of
     the record" rule of SIR 2016-5083's data screening.
   - At least 25 peaks in the retained record.
   - A drainage area on record, and a HUC in basin 09, 10 or 17.
   - Redundancy (SIR 2016-5083 eqs. B3-B4): of pairs with drainage-area ratio ≤ 5 and
     standardized distance ≤ 0.5, the gage with the shorter pseudo record length is dropped.
     Gage locations stand in for basin centroids.
3. **Station skews.** `skew_study.station_skew` fits each retained record with flowfreq's EMA
   and MGBT, station-only, applying peakfq's `siteQT` peak-code rules. Code 7 historic peaks
   and code 4/8 censored peaks enter EMA, and MGBT censors PILFs. This follows SIR 2016-5083
   app. B's handling of mixed-population and PILF sites.
   - The skew MSE is the fit's own `as_G_mse`: ADJE, or plain B17B `mseg` when MGBT finds low
     outliers, as `emafit.f:707` does.
   - The pseudo record length is the fit's `as_G_PRL_o` (eq. B2), clipped to [Ps, PH].
4. **Cross-correlation.** A Fisher-Z correlation-distance model, `Z = a + exp(b + c·D)`
   (eq. B11), is fitted to the concurrent log peaks of retained gage pairs with at least 50
   common years. Skew-estimator correlations then follow eqs. B8-B10 with κ = 2.8.
5. **Regression.** B-WLS for the parameters, then B-GLS for the model error variance and its
   precision, with an exponential prior of λ = 10 (SIR 2012-5130 app. 3).
   - Candidate explanatory variables: an east-of-Divide indicator (HUC 09/10 vs 17),
     log10 drainage area, gage altitude, and latitude.

Regenerate everything with:

```bash
python tools/build_montana_skew_study.py --cache <dir>
```

It paces its OGC requests and waits out any 429, since the API allows 1000 requests per hour
without a key. Downloads and per-site fits are cached, so an interrupted run resumes.

- Inputs: `data/skew_study/montana_provisional_v1_inputs.csv`. It has one row per Montana
  peak-flow site, with the reason for any exclusion.
- Results: `data/skew_study/montana_provisional_v1_results.json`.

## Validation of the tooling: Pacific Northwest

Before running Montana, the same code was run on the Pacific Northwest study's own inputs:

- SIR 2016-5083 Table B1: the station skew, PRL and basin centroid of the 290 gages used;
- eq. B11: the study's correlation-distance model, with κ = 2.8;
- each gage's historical period from its NWIS peak record through WY2012, since eq. B10 needs
  it and Table B1 does not print it.

The fit reproduces Tables B2 and B3 (`tests/test_skew_study.py`):

| | published | reproduced |
|---|---|---|
| β (SD) | -0.07 (0.10) | -0.065 (0.098) |
| σ²δ (SD) | 0.17 (0.022) | 0.171 (0.022) |
| ASEV | 0.010 | 0.0095 |
| AVPnew | 0.18 | 0.180 |
| Effective record length | 41 yr | 40.9 yr |
| Sampling error SS / model error SS | 35 / 50 | 34.9 / 49.5 |
| EVR | 0.7 | 0.70 |
| MBV* | 10 | 9.3 |
| High-influence gages (19) | listed p. 55 | 17 of 19 recovered; top two in order |

The gap is MBV* (and, less visibly, the SD of β). Both depend on the pseudo concurrent record
lengths, which here use NWIS-derived historical periods in place of the study's own PeakFQ
historic periods. κ = 3.0 gives MBV* 8.3 and κ = 3.3 gives 7.0, so 2.8 is closest.

The only way to go further is the study's own PeakFQ `.psf` historic periods, which were not
published. The station skews themselves were not recomputed: Table B1's PeakFQ 7.1 values were
used. The published effective record length is reproduced exactly by the *unbiased*
Griffis-Stedinger variance (1 + 6/n)² Var[G]. The same report's "17 years" for the Bulletin
17B map's MSE 0.302 uses the biased form.

## Diagnostics

**Data and screening.**

- Montana peak-flow sites: 998 stream sites with discharge peaks (27,676 peaks, WY ≤ 2025).
- Excluded:
  - 661 with fewer than 25 pre-regulation peaks;
  - 52 as redundant nested gages;
  - 1 with no HUC on record.
- Used: 284 gages, by basin:
  - 235 Missouri (HUC 10);
  - 46 Columbia (HUC 17);
  - 3 Hudson Bay / St. Mary (HUC 09).
- 13 used gages were cut at their first peak coded for regulation:
  - Missouri River at Fort Benton, at 1953 (Canyon Ferry);
  - Kootenai River at Libby, at 1972 (Libby Dam).

**Station skews.**

- Mean -0.02, median -0.10, standard deviation 0.69; 20 gages above +1 and 11 below -1.
- The Glacier/St. Mary gages average +1.29.
- MGBT found low outliers (PILFs) at 154 of the 284 gages, which then use the B17B skew MSE
  (`emafit.f:707`).
- Pseudo record length: median 45 yr, range 23-118. At 193 gages it exceeds the systematic
  count because of historic or censored information.

**Cross-correlation model.**

- Fitted to 2,893 pairs with at least 50 concurrent years:
  Z = 0.186 + exp(-0.257 - 0.00868 D), with D in miles between gages.
- The Pacific Northwest's is Z = 0.21 + exp(-0.17 - 0.0058 D). Montana's decays faster with
  distance.
- MBV* = 3.7: WLS alone would understate the variance of β by a factor of about 4, so the GLS
  precision analysis matters, though less than in the Pacific Northwest (10).
- EVR = 0.50, so sampling error is not negligible and WLS/GLS is warranted.

**High influence.** The ten most influential gages, in order:

| Gage | Name | Skew | PRL (yr) |
|---|---|---|---|
| 06078500 | North Fork Sun River near Augusta | 2.77 | 44 |
| 06109530 | Little Sandy Creek tributary near Virgelle | 2.88 | 24 |
| 12356000 | Skyland Creek near Essex | 2.54 | 25 |
| 06031950 | Cataract Creek near Basin | 2.03 | 51 |
| 06102500 | Teton River below South Fork near Choteau | 2.10 | 36 |
| 05014500 | Swiftcurrent Creek at Many Glacier | 1.68 | 111 |
| 12356500 | Bear Creek near Essex | 2.01 | 25 |
| 06125520 | Swimming Woman Creek tributary near Living Springs | 1.71 | 38 |
| 06336500 | Beaver Creek at Wibaux | 1.82 | 31 |
| 06088500 | Muddy Creek at Vaughn | 1.41 | 80 |

Most of them sit along the Rocky Mountain Front and in Glacier. Their records are dominated by
the June 1964 and 1975 rain-on-snowmelt floods: a mixed population of rare rainfall floods on a
snowmelt regime. That is the clearest structure in the residuals. A Front-proximity or
storm-type variable, not tried here, is the obvious candidate for a model that beats CONSTANT.

## Limitations

- **Not reviewed and not USGS.** The screening and model choices below are this study's own.
  None has had the site-by-site review a USGS study gives every gage.
- **The regulation screen uses peak codes only.** Uncoded regulation stays in:
  - the Bighorn (Yellowtail Dam) above Yellowstone River near Sidney;
  - the Tongue River Reservoir;
  - irrigation diversions on the plains, and the St. Mary Canal feeding the Milk.

  The GAGES-II/NID screen of roadmap #32 is not built yet. Dropping basins over 3,000 mi²
  barely changes the result (G = -0.02, AVPnew 0.37).
- **Mixed populations are handled only as EMA/MGBT handles them**, as SIR 2016-5083 app. B did.
  No gage was split into snowmelt and rain populations. The Rocky Mountain Front gages above
  dominate the residuals and the model error variance.
- **Gage locations stand in for basin centroids** in the redundancy screen and the
  cross-correlation distances. SIR 2016-5083 used basin centroids, which Montana lacks
  offline. For small basins the difference is small; for mainstem gages it is not.
- **Montana gages only.** Shared basins are cut at the state line: the Missouri headwaters in
  WY/ID, the Hudson Bay drainages in Alberta, and the Kootenai and Clark Fork into ID. Only 3
  Hudson Bay gages survive the screens.
- **Historical periods run from each gage's first to last recorded peak.** PeakFQ
  specification files can set longer historic periods, which would raise some pseudo record
  lengths.
- **Explanatory variables were few:** gage altitude (not mean basin elevation), drainage area,
  latitude and an east/west indicator. StreamStats basin characteristics were not collected.
- **Record screen: 25 years**, the minimum set for this study, against the Pacific Northwest's PRL of 35
  years. With PRL ≥ 35 the result is G = -0.07, AVPnew 0.33.
- **Station skews are flowfreq's own EMA.** They match peakfq 8.1.0 on the parity sites, but
  were not cross-checked against PeakFQ for these 284 gages.

Where this leaves Montana:

- East of the Divide, users still need a USGS-recommended skew (B17C p. 31: consult the USGS).
  The study here says a statewide constant would carry about 23 years of information.
- The better route is a study that models the Rocky Mountain Front effect and uses basin
  characteristics. That is the USGS WY-MT study SIR 2025-5019 anticipates, and roadmap §1.3's
  Montana item stays open for it.
- In `regional_skew.csv`, MT remains `pending`, and this result is recorded only here and in
  the roadmap.
