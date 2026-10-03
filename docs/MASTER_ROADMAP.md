# Master Roadmap: national flood frequency coverage

**Goal.** Make flowfreq able to produce a defensible flood-frequency estimate at **any gaged or
ungaged site in the 50 states, DC, and the five inhabited territories**. It should do this
with Bulletin 17C at gages, published USGS regional regression equations (RREs) at ungaged
sites, a weighted combination of the two where both exist, and documented methods for
transposition, nonstationarity, and future-condition quantiles.

**Status legend:** `[ ]` open · `[~]` partial (scaffolded; the linked issue carries the live status) · `[x]` done in this repo today.
**Phase A epic:** [#27](https://github.com/pinhead001/flowfreq/issues/27) · **Wave 1 epic:** [#28](https://github.com/pinhead001/flowfreq/issues/28) · generated per-state status: `docs/REGRESSION_COVERAGE.md`.
**Relationship to `TODO.md`:** `TODO.md` is the day-to-day log of open items and history.
This file is the long-range plan. When an item here starts, give it an entry in `TODO.md`.

Standing rules carried over from this repo's history, which apply to every phase:

- **Every coefficient has a citation.** An exponent, a regression constant, or a regional skew
  with no report ID and table number is not accepted (see `RegressionExponents`).
- **A 200 is not an answer.** Validate inputs against the published parameter limits
  *before* evaluating, not after (see `STREAMSTATS_NSS_ADDENDUM.md`).
- **Every data source is verified live before code is written against it.** Earlier guesses
  about endpoints were wrong three times.
- **Each jurisdiction's equations are validated against at least one worked example from
  its source report** before that jurisdiction is marked supported.
- Do not edit `vendor/`. Known gaps are `xfail(strict=True)`.

---

## 0. Where the repo stands (baseline)

| Capability | State | Where |
|---|---|---|
| B17C EMA + MGBT, parity with peakfq 8.1.0 | [x] | `bulletin17c.py`, `tests/fortran_parity/` |
| Fortran engine as a selectable backend | [x] | `fortran_engine.py`, `peakfqr/` |
| NWIS peak/daily retrieval | [~] legacy NWIS services only; OGC backend is a stub | `usgs.py`, `peak_sources.py` |
| StreamStats delineation + basin characteristics | [x] Phase 1, with validated watershed polygon | `streamstats.py` |
| NSS regression evaluation via web service | [~] PFS verified live in WA only | `streamstats.py` |
| Regional skew | [~] verified PNW rows (WA/OR/ID); location lookup by hydrologic unit (`regional_skew_at`, `regional_skew_for_site`); other states pending (§1.3) | `regional_skew.py`, `data/regional_skew*.csv` |
| Offline regression-equation library | [~] schema and evaluator; no equations transcribed | `regression/`, `data/regression/` |
| Gage/regression weighting (B17C App. 9 style) | [ ] | — |
| Drainage-area-ratio transposition | [x] | `transpose.py` |
| QPPQ daily-series transfer | [x] | `qppq.py` |
| Record extension (MOVE) | [ ] | — |
| Nonstationarity tests / models | [ ] | — |
| Future-condition quantiles | [~] change-factor framework; no factor sets | `future_flow.py` |
| Low-flow frequency, regime, sub-daily metrics | [x] | `lowflow.py`, `regime.py`, `subdaily.py` |

---

## 1. Data foundation (prerequisite for everything else)

### 1.1 Gage data retrieval
- [~] [#29](https://github.com/pinhead001/flowfreq/issues/29) **Migrate off legacy NWIS web services.** USGS is retiring `waterservices.usgs.gov`
      in favor of the Water Data OGC APIs (`api.waterdata.usgs.gov`). Verify the peak,
      daily-value, and site-metadata endpoints live, then add an adapter layer so
      `USGSgage` works with both backends. Existing fixtures must stay byte-identical.
      *Substantive:* the `PeakDataBackend` protocol, `validate_peak_frame`, and
      `LegacyNwisBackend` (`peak_sources.py`). *Stub:* `WaterDataApiBackend.fetch_peaks`
      raises `NotImplementedError`, and `USGSgage` does not route through the adapter.
- [~] [#30](https://github.com/pinhead001/flowfreq/issues/30) Parse **all peak qualification codes** (historic `7`, estimated, regulated `6`/`C`,
      urbanization, dam failure, `<`/`>` values) into `PeakRecord` flags. Map each code to
      an explicit B17C treatment: include, censor, exclude, or flag for review.
- [~] [#31](https://github.com/pinhead001/flowfreq/issues/31) Retrieve and store **perception thresholds and historical-period information** when
      published. Accept a PeakFQ `.psf` specification file as input so a published analysis
      can be replicated exactly (`psf` → `EMAParameters`).
- [~] [#32](https://github.com/pinhead001/flowfreq/issues/32) **Regulation and urbanization screen.** Use the peak codes, GAGES-II/NID dam storage,
      and NLCD impervious fraction to classify each gage as reference, regulated, or urban.
      Refuse, or require an override for, B17C on a regulated record.
      *Substantive:* `flowfreq/regulation.py` and `data/regulation_screen.csv.gz` (all
      9,322 GAGES-II gages, built by `tools/build_regulation_screen.py`): GAGES-II `Ref`
      → reference; NID-2009 storage over 127.8 days of mean runoff (Dudley and others
      2018) → regulated; NLCD-2006 impervious over 5% (SIR 2016-5118 p. 23) → urban;
      peak code 6 → regulated. `analyze_gage` refuses a regulated gage, and `run_ffa`
      does too when it is given `site_no=`, unless `allow_regulated=True`, which is
      recorded in provenance. *Open:* current-NID storage (the live national CSV needs
      basin polygons to sum upstream), the NLCD impervious time series, and gages outside
      GAGES-II, which only have peak-code evidence.
- [~] [#33](https://github.com/pinhead001/flowfreq/issues/33) Offline **gage catalog**: expand `flowfreq/data/gage_attributes.csv` (3 rows today) into
      a versioned national table covering every active and inactive peak-flow site with
      ≥10 years of record. Columns: site_no, name, lat/lon, DA, state, HUC8, years of
      record, regulation class, and the regression region it falls in.
      Rebuild it with a script in `tools/`; never edit it by hand.
      *Substantive:* `data/gage_catalog.csv.gz` (national, built 2026-10 from the Water
      Data API `peaks` and `monitoring-locations` collections by
      `tools/build_gage_catalog.py`), with the schema validated on load. The regulation
      class comes from #32. *Partial:* `regression_region` is filled for Wave 1 states
      only, by point-in-polygon against the StreamStats peak-region layers. Other states
      are added in their waves.

### 1.2 Geospatial inputs
- [x] **Watershed polygon** (the open Phase 1 gap). Verify a live source, such as the
      StreamStats `delineate` geometry or the NLDI `basin` endpoint, and populate
      `WatershedCharacteristics.polygon_geojson`. *Done 2026-09-27:* it is in the
      `delineate/sshydro` response already fetched; verified live and validated against
      `DRNAREA` (`docs/STREAMSTATS_MODULE_DESIGN.md` S10).
- [ ] **NLDI integration** for upstream/downstream navigation. This is needed for
      same-stream donor search (§5) and for nested-gage checks.
- [ ] **Regression-region polygons** for every jurisdiction. Source these from the report
      data releases on ScienceBase or from the StreamStats state layers, and store them as
      compressed GeoPackage/Parquet with the citation. This removes the open "region
      selection without a polygon" item: point-in-polygon picks the region, including
      area-weighted multi-region basins. *Partly superseded 2026-09-27:* online, NSS's
      own `regressionregions/bylocation` now selects regions from the watershed polygon,
      with area weights (`streamstats.locate_regression_regions`, addendum S5). An
      offline store is still needed for the offline backend.
- [ ] **Local basin-characteristic computation (fallback).** For jurisdictions or
      characteristics StreamStats does not serve, compute the common RRE variables (DRNAREA,
      PRECIP, slope, forest %, impervious %, storage %, elevation) from NHDPlus HR, 3DEP,
      NLCD, PRISM, and SSURGO. Validate each one against StreamStats on ≥10 basins per state
      and document the tolerance.

### 1.3 Regional skew
- [~] [#34](https://github.com/pinhead001/flowfreq/issues/34) **National regional-skew table.** For each jurisdiction, record the currently adopted
      skew source: the B17C-era Bayesian WLS/GLS study value and MSE, or a documented
      fallback. Do not use the B17B Plate I map silently.
      The table records the value (or map/raster), MSE, effective record length, report
      citation, and validity region.
- [x] `regional_skew_at(lat, lon)` lookup returning `(skew, mse, citation)` and plugging
      straight into `Bulletin17C`. Raise, rather than default, where no study applies.
      (The silent `-0.302` default is gone: since #44, `run_ffa`/`compare_engines` and the
      CLI require a supplied skew, station-only, or an explicit `use_default_skew=True`.)
      Membership is by hydrologic unit (`data/regional_skew_hucs.csv`, built by
      `tools/build_regional_skew_hucs.py`): `regional_skew_at_huc`, `regional_skew_at`
      (WBD point query) and `regional_skew_for_site` (the monitoring location's HUC). Only
      the Pacific Northwest study is defined so far, including its Snake River Plain
      exclusion. The Oregon closed basins and the Klamath are left `unresolved`: the report's
      region wording omits them, although 5 of its gages are there. Every further verified
      study needs its own HUC rows.
- [x] Tooling to **develop a new regional skew** (B-WLS/B-GLS, Veilleux/Reis-Stedinger) for
      jurisdictions without one. This is lower priority and research-grade.
      `flowfreq.skew_study`: B-WLS/B-GLS with a correlation-distance cross-correlation model
      and the full diagnostic set. Validated by refitting the Pacific Northwest CONSTANT
      model from SIR 2016-5083 Table B1, which reproduces Tables B2/B3 to their printed
      precision. MBV* is the exception: 9.3 against 10, because the gage historical periods
      are approximated from NWIS.
- [ ] **Complete a Bayesian (B-WLS/B-GLS) regional skew study for Montana.** Montana is the
      one Wave 1 state with no B17C-era skew: SIR 2025-5019 (pp. 9-10) states no B-WLS/B-GLS
      study covers any part of MT/ND/SD/WY, and USGS Montana applies the Bulletin 17B map
      statewide (SE 0.55; SIR 2018-5046 used 0.64), which B17C p. 31 does not recommend. The
      Pacific Northwest study (SIR 2016-5083 app. B, G = -0.07, MSE 0.18) covers only western
      MT's Columbia basin (HUC 1701, 23 gages), and USGS Montana does not adopt it.
      `flowfreq/data/regional_skew.csv` keeps MT `pending` until this exists.
      Scope: the Missouri basin east of the Divide (and its SW-MT headwaters) and the St.
      Mary/Belly (Hudson Bay) drainages at minimum, ideally statewide for consistency with
      USGS MT practice; reuse the gage catalog (#33), peak-code screening (#30/#32), and the
      tooling item above. Coordinate with the USGS WY-MT Water Science Center, which B17C
      directs users to consult and which SIR 2025-5019 says anticipates such a study -- a
      published USGS study, once available, supersedes anything developed here.
      *Provisional, not USGS:* `docs/MONTANA_REGIONAL_SKEW_PROVISIONAL.md` gives a statewide
      B-WLS/B-GLS result built from the Water Data OGC API with flowfreq's EMA/MGB. Its
      inputs are in `data/skew_study/` and `tools/build_montana_skew_study.py` regenerates
      them. It is recorded in that report and here only, not in `regional_skew.csv`: that
      table has no provisional status, and MT stays `pending`. This item stays open for a
      published USGS study.

---

## 2. Bulletin 17C completeness at gaged sites

- [x] EMA, MGBT, weighted skew, CIs (parity with peakfq 8.1.0).
- [ ] Close the Cains Coulee `as_G_mse` discrepancy or document it as permanent (`TODO.md` P3).
- [ ] **Batch B17C for a whole state.** Run every catalog gage with its regional skew and
      emit a table matching the form of the state report's at-site appendix. Diff it against
      the published appendix values as a per-state validation (§7).
- [ ] **B17C Appendix 8 record extension** (MOVE.1/MOVE.3) from a correlated long-record
      gage, feeding the extended record into EMA with correct effective-record accounting.
- [ ] **B17C Appendix 9-style weighting** of the at-site estimate with the regression
      estimate by inverse variance: `log Qw = (log Qg·V_r + log Qr·V_g)/(V_g+V_r)`. This
      needs the regression prediction variance (§3.3) and the EMA variance, which already
      exists.
- [ ] Mixed-population analysis (e.g., snowmelt vs. rain, tropical vs. frontal). At minimum,
      detect and warn; then support separate fits combined by the probability-union rule.
- [ ] Crippen-Bue / envelope-curve sanity check on every quantile. It is a flag, not a
      correction.

---

## 3. Regional regression equations: all states and territories

### 3.1 Architecture
- [~] [#35](https://github.com/pinhead001/flowfreq/issues/35) **Equation schema** (`flowfreq/regression/`): a typed, serializable definition holding
      region ID, statistic (AEP), functional form (log-linear, power, with transforms such as
      `log10(X+1)` and rescaled forms such as `log10(X/100 + 1)`, via `Variable.scale`
      and `Variable.offset`), coefficients, variable definitions and units, calibrated
      min/max per variable, SEP / average variance of prediction, model-error variance, and
      the `(XᵀΛ⁻¹X)⁻¹` matrix where published. Carry the report citation, table number, and
      effective date.
- [ ] **Two evaluation backends** behind one interface:
      1. *Online:* the existing NSS `Scenarios/Estimate` client.
      2. *Offline:* a local evaluator over the equation library, so results are reproducible
         with no network and pinned to a library version.
      A cross-check test asserts that both backends agree to a tight tolerance for every
      equation.
      Both halves exist separately (`streamstats.estimate_flow_statistics`,
      `regression.evaluate`); the shared interface and the cross-check do not.
- [ ] **Bootstrapping the offline library from NSS.** `GET /nssservices/regions/{r}/Scenarios`
      returns templates with variable limits but **no equation field**; equation strings
      appear only in `POST /nssservices/Scenarios/Estimate` results. So
      `tools/snapshot_nss.py` (not yet written) must POST synthetic in-range inputs for every
      region, parse each returned equation string into intercept and coefficients, and dump
      every region/statistic group to versioned JSON under `flowfreq/data/regression/`, with
      a diff report on re-run so equation changes upstream are visible in review. NSS does
      not carry the report table number, covariance, model-error variance or `n_sites`;
      those come from the report (`TODO.md` #35).
- [x] Out-of-range handling: raise by default; `allow_extrapolation=True` records the
      violation (same pattern as `transpose.py`). `regression.evaluate` raises
      `OutOfRangeError`; with the flag it logs and records each violation in
      `RegressionEstimate.out_of_range`.
- [~] **Multi-region basins**: area-weighted estimates, following the method each state's
      report prescribes, since some reports weight flows and others weight logs.
      The generic mechanism exists: `regression.evaluate_weighted(..., space="log"|"linear")`,
      with the caller choosing the space. Per-state rules (which space, and variants such as
      OR's elevation interpolation between regions 2A/2B) are not recorded anywhere, and
      the weighted result carries no prediction interval.
- [x] **Prediction intervals** at ungaged sites from the equation's variance terms
      (Tasker & Driver form). Fall back to average SEP where the covariance matrix is not
      published, and label it as such. `regression.equations._prediction_sd`: model error +
      `xᵀUx` when covariance and model-error variance are present, else AVP, else SEP;
      `RegressionEstimate.interval_method` names which. Tested on synthetic equations only.
- [ ] **Region-of-influence (ROI)** evaluator for states whose published method is ROI
      rather than fixed equations. This requires the state's gage dataset and the ROI
      parameters from the report.
- [ ] Urban adjustment equations (e.g., Sauer et al. national urban equations, and state
      urban equations where published), applied only where the report says they apply.

### 3.2 Jurisdiction coverage, rolled out in regional waves
`docs/REGRESSION_COVERAGE.md` is generated by `tools/gen_regression_coverage.py` rather
than typed, and `tests/test_catalog.py` fails if it is stale. It gets one row per
jurisdiction. Target columns: NSS peak-flow equations present? · report
citation(s) and year · number of regions · ROI? · urban equations? · regional skew source ·
future-flow guidance source · validation worked example located? · offline library
validated? · status.

**Why waves by river basin instead of alphabetical order:** regional skew studies,
regression calibration datasets, and donor gages cross state lines along major basins.
Neighbors share gages, border effects, and often a single USGS report, so doing them
together saves re-work and exposes cross-border inconsistencies. Each state is assigned
to exactly one wave. Where a state drains to several basins, the assignment follows the
basin that dominates its flood hydrology.

| Wave | Region | Jurisdictions (in order) | Notes |
|---|---|---|---|
| 1 | Columbia River basin ([#28](https://github.com/pinhead001/flowfreq/issues/28)) | [#37](https://github.com/pinhead001/flowfreq/issues/37) WA (pilot) · [#38](https://github.com/pinhead001/flowfreq/issues/38) OR · [#39](https://github.com/pinhead001/flowfreq/issues/39) ID · [#40](https://github.com/pinhead001/flowfreq/issues/40) MT | WA already has a live-verified NSS peak-flow estimate, which makes it the pilot that proves the per-state pipeline end to end. Snowmelt/rain mixed populations and regulated mainstem gages are the main hazards. |
| 2 | Colorado River basin | [ ] CO · [ ] UT · [ ] WY · [ ] NM · [ ] AZ · [ ] NV | CA is a Colorado compact state but goes to wave 3. Arid regions have large SEPs, zero-flow years, and heavy use of PILF (potentially influential low flood) screening. |
| 3 | West coast and Pacific | [ ] CA · [ ] AK · [ ] HI | OR and WA were done in wave 1. AK and HI are grouped here as Pacific states with sparse networks. |
| 4a | Upper Mississippi and Great Lakes | [ ] MN · [ ] WI · [ ] IA · [ ] IL · [ ] MO · [ ] MI | MI is mostly Great Lakes drainage, placed here as the nearest wave. |
| 4b | Missouri basin | [ ] ND · [ ] SD · [ ] NE · [ ] KS | Snowmelt/ice-jam peaks and prairie non-contributing area. |
| 4c | Ohio basin | [ ] OH · [ ] IN · [ ] KY · [ ] TN · [ ] WV | TN holds the Big Sandy primary test site. |
| 4d | Lower Mississippi, Arkansas-Red, and Gulf | [ ] AR · [ ] LA · [ ] MS · [ ] OK · [ ] TX | TX is mostly Gulf drainage, placed here for its Red River border and shared skew work with OK/LA. |
| 5 | East coast | [ ] ME · [ ] NH · [ ] VT · [ ] MA · [ ] RI · [ ] CT · [ ] NY · [ ] NJ · [ ] PA · [ ] DE · [ ] MD · [ ] DC · [ ] VA · [ ] NC · [ ] SC · [ ] GA · [ ] FL · [ ] AL | AL (Mobile basin, Gulf) is grouped with GA/FL because Southeast regional skew and regression studies span those states. |
| 6 | Territories | [ ] PR · [ ] USVI · [ ] Guam · [ ] CNMI · [ ] American Samoa | Expect old, sparse, or absent equations. Document the fallback for each: nearest applicable equations, index-flood from the few gages, or "not supported" with the reason. |

All 56 jurisdictions are covered: 50 states, DC, and 5 territories. Wave boundaries are
milestones. Do not start wave *n+1* until wave *n*'s cross-border consistency check (below)
passes.

Per jurisdiction, the definition of done:
- [ ] Equations for every published peak AEP (typically 50% down to 0.2%) are in the offline
      library with citation and limits.
- [ ] One or more published worked examples are reproduced to the report's printed precision.
- [ ] Online and offline backends agree.
- [ ] A regional skew is available (§1.3).
- [ ] Region polygons are loaded and point-in-polygon selection is tested on a known site.
- [ ] The state's at-site B17C appendix is re-run and compared (§7).
- [ ] **State-specific future-flow guidance** is researched and recorded (§6.3.2). "None
      published" counts as a valid, recorded answer.

Per wave, the definition of done:
- [ ] **Cross-border consistency check:** evaluate each state's equations on border gages
      and on basins that straddle the line. Report the discontinuity at the border by AEP.
      Large jumps are documented, not smoothed.
- [ ] Transposition LOOCV (§5.2) run on the wave's gage network.
- [ ] A wave vignette with one worked ungaged-site example per state.

**Expected gaps to plan for (verify each one; these are not assumptions to code against):**
- Some states publish equations only in report PDFs that are not yet in NSS, or have newer
  reports than NSS carries. Track the "newest report" separately from the "NSS version".
- Arid Southwest and Alaska: many regions use few variables and have very large SEP. Surface
  the SEP prominently.

### 3.3 Regression development toolkit (for gaps and updates)
- [ ] OLS/WLS/**GLS regression** (Stedinger & Tasker) with cross-correlation from concurrent
      record, in Python, matching the capabilities of the USGS WREG R package.
      Validate it against WREG output on a published dataset.
- [ ] Model diagnostics: VIF, leverage/influence, pseudo-R², SEP, AVP, residual maps.
- [ ] Hydrologic-region delineation aids: residual clustering and a region-of-influence
      search.

---

## 4. Ungaged-site estimation workflow

A single entry point, `estimate_at_site(lat, lon, …)`, that decides and documents:

- [ ] **On a gage** → B17C at-site, weighted with RRE (§2).
- [ ] **Near a gage on the same stream** (area ratio within the state's band, typically
      0.5–1.5) → regression-weighted transposition (§5.1).
- [ ] **Otherwise** → RRE (offline or NSS), with prediction intervals.
- [ ] Every result carries full provenance: method chosen and why, donor(s), equations,
      citations, out-of-range flags, and the per-quantile arithmetic (as
      `TransposedResults.to_markdown` does today).
- [ ] Report output suitable for a FEMA/DOT submittal (`report.py` extension).

---

## 5. Transposition of gage data to ungaged locations

### 5.1 Methods to implement (most common in USGS state reports first)
- [x] **Drainage-area ratio** with regression exponent `b(p)` (`transpose_frequency`).
- [ ] **Regression-weighted drainage-area ratio** (the Sauer / USGS state-report method):
      `Q_u,w = Q_u,r · [R − (2|A_g−A_u|/A_g)(R−1)]`, with `R = Q_g,w/Q_g,r`, valid for
      `0.5 ≤ A_u/A_g ≤ 1.5`. This is the form most state reports and StreamStats use for
      sites near a gage; each state's exact variant must be taken from its report.
- [ ] **Interpolation between two gages** on the same stream (upstream and downstream),
      per the state's report where one is prescribed.
- [x] **QPPQ** daily-series transfer (`qppq.py`). Follow-up: **map-correlation donor
      selection** (Archfield & Vogel 2010) as an alternative to nearest/most-similar donor.
- [ ] **Multiple-donor weighting** (inverse-distance or inverse-variance) and donor ranking by
      hydrologic similarity (DA, slope, precip, BFI screen already in `transpose_low_flow`).
      [#13](https://github.com/pinhead001/flowfreq/issues/13) measures why area alone is not
      enough (Methow, orographic precipitation) and asks for a predictor-similarity screen.
- [ ] **Record extension** (MOVE.1/MOVE.3/KTRL) to lengthen a short gage before transposing
      (shares code with §2's Appendix 8 item).
- [ ] **Index-flood / regional growth curves** as an option for data-sparse regions
      (territories, Alaska), clearly labeled as non-USGS-standard where it is.
- [ ] **Geostatistical methods** (top-kriging along the network): research-grade, optional.

### 5.2 Guardrails
- [x] Area-ratio band, probability-kind checks, extrapolation flags (existing pattern).
- [ ] Same-stream / nested-basin check via NLDI. Refuse cross-basin DAR unless forced.
- [ ] Regulation check. Refuse a regulated donor for an unregulated target and the
      reverse, unless forced.
- [ ] Leave-one-out validation harness for each method on a state's gage network. This
      reports bias and RMSE in log space by AEP, so method choice is evidence-based per
      region.

---

## 6. Nonstationarity and future flood quantiles

B17C assumes stationarity and says to evaluate it. It does not prescribe a nonstationary
method. Everything here is **advisory tooling with explicit labeling**. None of it replaces
the B17C estimate in a regulatory product unless the reviewing agency accepts it.

### 6.1 Detection (implement first; widely accepted)
- [ ] **Monotonic trend tests** on annual peaks and on log peaks: Mann-Kendall (with
      Hamed-Rao / pre-whitening options for autocorrelation), Sen's slope, and Spearman's rho.
- [ ] **Change-point tests**: Pettitt, CUSUM, and Lombard (smooth and abrupt), matching the
      set in USACE's Nonstationarity Detection Tool, so results are comparable with
      USACE ECB 2018-14 practice.
- [ ] Trends in **variance and in peaks-over-threshold counts**, not only in the mean.
- [ ] **Field significance** across a region (Walker / FDR, block-bootstrap for spatial
      correlation). Single-site significance is expected by chance in a statewide run.
- [ ] **Attribution screen**: overlay detected change points on regulation dates (NID dam
      completion year), urbanization history (NLCD time series), gage relocation, and
      datum/rating changes. A change point that lines up with a dam is not a climate signal.
- [ ] `stationarity_report(peaks)` giving a one-page summary with plots. **Recommendation:**
      run this on every B17C analysis by default and attach it to the output.

### 6.2 Handling detected nonstationarity (recommended hierarchy)
1. [ ] **Remove the cause if it is known and physical.** Analyze only the homogeneous period
       (post-dam, post-urbanization) using B17C with historical information where possible,
       or adjust the record, for example with an urbanization adjustment. This is the most
       defensible option and is B17C-consistent.
2. [ ] **Time-varying-parameter LP3/GEV.** Location, and optionally scale, are linear in
       time or in a covariate such as a climate index, urban %, or CO₂. Estimate by MLE
       and Bayesian MCMC, with GEV/GLO comparison. Use profile-likelihood or posterior
       intervals, a likelihood-ratio test and AIC/BIC against the stationary model, and
       reject the nonstationary model unless it is clearly supported.
3. [ ] **Covariate models with physically meaningful drivers** (ENSO/PDO/AMO indices,
       basin precipitation, impervious fraction) are preferred over time as a covariate,
       because time trends extrapolate without bound.
4. [ ] Report **design-life risk** rather than a single "return period" under
       nonstationarity: the probability of ≥1 exceedance over an N-year service life, and
       expected waiting time (Salas & Obeysekera; Read & Vogel).

### 6.3 Future flood quantiles

Build the **general, nationally applicable framework first (6.3.1)**. Then add
**state-specific guidance (6.3.2) one state at a time, as each state is developed in the
§3.2 waves**. Never add state guidance ahead of that state's regression work.

#### 6.3.1 General framework (build first, applies everywhere)
Recommended approaches, most defensible first:
- [~] [#36](https://github.com/pinhead001/flowfreq/issues/36) **Change-factor ("delta") scaling of B17C/RRE quantiles.** Use a pluggable
      `ChangeFactorSet` with per-AEP multipliers, scenario, horizon, source citation, and
      applicable geography. The first factor sets are national/federal sources only: FHWA
      HEC-17, NCHRP 15-61 guidance, and NOAA Atlas 15 future precipitation once released,
      mapped to flow through a documented precipitation-to-flow elasticity. This is the
      simplest option to review and should be the default.
      *Substantive:* the framework (`future_flow.py`), the `data/future/*.factors.json`
      loader (`available_factor_sets`), and a review of the national sources against the
      primary documents (`data/future/national_sources.json`, `docs/FUTURE_FLOW_GUIDANCE.md`).
      HEC-17 and NCHRP 15-61 give methods and single-site examples, **not tabulated national
      factors**, so no national set ships; NOAA Atlas 15 Volume 2 is unpublished (planned
      2027). Revisit then.
- [ ] **Regression space-for-time.** In RREs whose explanatory variables include
      precipitation or temperature, substitute downscaled projected values (e.g.,
      LOCA2 / CMIP6 ensemble). Report the ensemble spread, not only the median, and flag
      it when projected values fall outside the equation's calibrated range, which is
      common.
- [ ] **Nonstationary-model projection** (§6.2 #2/#3) by extrapolating covariates to a
      horizon. Allow this only with a covariate projection from an external source, never by
      extrapolating a fitted time trend by itself.
- [ ] **Rainfall-runoff pathway** (continuous simulation or design storm with future IDF
      curves). This is out of scope for the library core; provide an interface to accept
      simulated annual-maximum series and run B17C on them.
- [ ] Uncertainty accounting that combines sampling (EMA CI), model (regression SEP), and
      climate-ensemble spread. Present these separately and never collapse them into one
      number without showing the parts.
- [~] Output labeling: every future-condition quantile carries scenario, horizon, method, and
      source, and is visually distinct from the regulatory B17C value in plots and reports.
      `apply_change_factors` keeps `regulatory_flow_cfs` and `future_flow_cfs` as separate
      columns, and `FutureQuantiles.provenance` carries scenario, horizon, citation, level and
      legal status. Nothing in `plots.py`/`report.py` renders it yet.
- [~] **Precedence rule:** a state factor set, when present and applicable, overrides the
      national default. The override is always recorded in provenance, and the national
      result is still reported alongside it for comparison.
      `future_flow.select_factor_set` implements the override (and raises if two state sets
      apply), and `FutureQuantiles.provenance["overrode"]` records it. The side-by-side
      national result is not produced automatically; the caller applies both sets.

#### 6.3.2 State-specific guidance (added per state, during its wave)
For each state, as part of its definition of done in §3.2:
- [ ] Survey the state DOT drainage/hydraulics manual, the state water-resources or
      environmental agency, the state climate office, and any USGS/state cooperative
      future-flow study.
- [ ] Record the findings in `docs/FUTURE_FLOW_GUIDANCE.md`, one section per state, with
      these fields: source and citation, legal status (required, recommended, or
      informational), method type (multiplier, projected-precip regression, or
      simulation), scenario, horizon, and applicable geography.
- [ ] Where the guidance is quantitative, encode it as a `ChangeFactorSet` under
      `flowfreq/data/future/` with a test reproducing one published example.
- [ ] Where the guidance is qualitative only, record it as text, and have outputs cite it
      without inventing numbers.
- [ ] Where no guidance exists, record "none published as of <date>". The national default
      then applies and is labeled as such.

---

## 7. Validation and QA program

- [ ] **Per-state at-site validation:** re-run B17C on each state report's gage appendix and
      compare quantiles. The expected tolerance is small but non-zero, because record
      lengths differ from the report's data cutoff, so compare against the report's cutoff
      year.
- [ ] **Per-state regression validation:** worked examples (§3.2).
- [ ] **Transposition LOOCV** per region (§5.2).
- [ ] **Nonstationarity methods** validated against published USGS/USACE case studies
      (e.g., USACE NSD tool outputs on the same record).
- [ ] Golden-file regeneration tooling analogous to `tools/gen_fortran_golden.py` for the
      regression library snapshot.
- [ ] CI: offline tests always; `requires_network` live checks (NSS, Water Data API) on a
      scheduled workflow, not on every push. This handles the egress limits recorded in
      `TODO.md`.

---

## 8. Packaging, data management, docs

- [ ] Data-size budget: regression library JSON (small), region polygons (large). Consider a
      separate `flowfreq-data` package or a lazy download with checksum pinned to the
      library version.
- [ ] Data versioning and a `flowfreq.data_versions()` report embedded in every output's
      provenance.
- [ ] Vignettes: "Ungaged site in state X", "Gage-to-ungaged transposition",
      "Stationarity screening", "Future-condition quantiles with change factors".
- [ ] API-stability policy for the new modules; mypy-enforced by default (see
      `pyproject.toml` override list).
- [ ] Coordinate pin bumps with `pinhead001/flowfreq-app` as each phase ships.

---

## 9. Suggested sequencing

| Phase | Scope | Unblocks |
|---|---|---|
| A | §1.1 Water Data API migration, gage catalog; §1.3 regional-skew table schema | Everything |
| B | §3.1 equation schema + NSS snapshot tool + offline evaluator; §3.2 matrix generated | Regression everywhere |
| C | §1.2 region polygons + watershed polygon; §3.1 multi-region, prediction intervals | Correct region selection |
| D | §2 App. 9 weighting, App. 8 record extension; §5.1 regression-weighted DAR | Full ungaged workflow |
| E | §4 `estimate_at_site`, proven end to end on the **WA pilot** | The per-state template |
| F | §6.1 stationarity detection + report; §6.3.1 general future-flow framework | Default screening; national future-flow default |
| G | **State rollout by wave** (§3.2): 1 Columbia → 2 Colorado → 3 West coast/Pacific → 4a–d Mississippi → 5 East coast → 6 Territories. Each state includes its §6.3.2 guidance. | "Supported" status per jurisdiction |
| H | §6.2 nonstationary models; §3.3 GLS toolkit, ROI, geostatistics | Gap-filling and research |

Phases A–E are the critical path to the first supported state. F can proceed in parallel
because it depends only on peak data, which already exists. G is the long tail: it is
mostly data work per state, once the machinery exists.

---

## 10. Tracking: roadmap file vs. GitHub issues vs. PRs

- **This file is the plan.** Land it on `main` through one docs-only PR, so everyone
  reviews and works from the same version.
- **GitHub issues are the execution tracker.** Create them just in time, not all up front:
  - One **milestone per phase A–F** and **one per wave**.
  - One issue per roadmap item in the current phase.
  - One issue per state, created when its wave starts, using a template that copies the
    §3.2 definition-of-done checklist. Creating 56 issues today would mostly produce
    stale tickets.
- **Each issue is closed by its own PR.** Keep each PR to roughly one item or one state, so
  it stays reviewable and CI stays green.
- **Avoid two sources of truth for status.** When an item gets an issue, replace its
  checkbox here with the issue link and let the issue carry the status. Update this file
  only when the *plan* changes.
