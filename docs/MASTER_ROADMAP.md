# Master Roadmap: national flood frequency coverage

**Goal.** Make flowfreq able to produce a defensible flood-frequency estimate at **any gaged or
ungaged site in the 50 states, DC, and the five inhabited territories**. It should do this
with Bulletin 17C at gages, published USGS regional regression equations (RREs) at ungaged
sites, a weighted combination of the two where both exist, and documented methods for
transposition, nonstationarity, and future-condition quantiles.

**Status legend:** `[ ]` open · `[~]` partial · `[x]` done in this repo today.
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
| NWIS peak/daily retrieval | [~] legacy NWIS services only | `usgs.py` |
| StreamStats delineation + basin characteristics | [x] Phase 1 (no polygon) | `streamstats.py` |
| NSS regression evaluation via web service | [~] PFS verified live in WA only | `streamstats.py` |
| Regional skew | [~] caller supplies a scalar | `Bulletin17C(regional_skew=…)` |
| Offline regression-equation library | [ ] | — |
| Gage/regression weighting (B17C App. 9 style) | [ ] | — |
| Drainage-area-ratio transposition | [x] | `transpose.py` |
| QPPQ daily-series transfer | [x] | `qppq.py` |
| Record extension (MOVE) | [ ] | — |
| Nonstationarity tests / models | [ ] | — |
| Future-condition quantiles | [ ] | — |
| Low-flow frequency, regime, sub-daily metrics | [x] | `lowflow.py`, `regime.py`, `subdaily.py` |

---

## 1. Data foundation (prerequisite for everything else)

### 1.1 Gage data retrieval
- [ ] **Migrate off legacy NWIS web services.** USGS is retiring `waterservices.usgs.gov`
      in favor of the Water Data OGC APIs (`api.waterdata.usgs.gov`). Verify the peak,
      daily-value, and site-metadata endpoints live, then add an adapter layer so
      `USGSgage` works with both backends. Existing fixtures must stay byte-identical.
- [ ] Parse **all peak qualification codes** (historic `7`, estimated, regulated `6`/`C`,
      urbanization, dam failure, `<`/`>` values) into `PeakRecord` flags. Map each code to
      an explicit B17C treatment: include, censor, exclude, or flag for review.
- [ ] Retrieve and store **perception thresholds and historical-period information** when
      published. Accept a PeakFQ `.psf` specification file as input so a published analysis
      can be replicated exactly (`psf` → `EMAParameters`).
- [ ] **Regulation and urbanization screen.** Use the peak codes, GAGES-II/NID dam storage,
      and NLCD impervious fraction to classify each gage as reference, regulated, or urban.
      Refuse, or require an override for, B17C on a regulated record.
- [ ] Offline **gage catalog**: expand `flowfreq/data/gage_attributes.csv` (3 rows today) into
      a versioned national table covering every active and inactive peak-flow site with
      ≥10 years of record. Columns: site_no, name, lat/lon, DA, state, HUC8, years of
      record, regulation class, and the regression region it falls in.
      Rebuild it with a script in `tools/`; never edit it by hand.

### 1.2 Geospatial inputs
- [ ] **Watershed polygon** (the open Phase 1 gap). Verify a live source, such as the
      StreamStats `delineate` geometry or the NLDI `basin` endpoint, and populate
      `WatershedCharacteristics.polygon_geojson`.
- [ ] **NLDI integration** for upstream/downstream navigation. This is needed for
      same-stream donor search (§5) and for nested-gage checks.
- [ ] **Regression-region polygons** for every jurisdiction. Source these from the report
      data releases on ScienceBase or from the StreamStats state layers, and store them as
      compressed GeoPackage/Parquet with the citation. This removes the open "region
      selection without a polygon" item: point-in-polygon picks the region, including
      area-weighted multi-region basins.
- [ ] **Local basin-characteristic computation (fallback).** For jurisdictions or
      characteristics StreamStats does not serve, compute the common RRE variables (DRNAREA,
      PRECIP, slope, forest %, impervious %, storage %, elevation) from NHDPlus HR, 3DEP,
      NLCD, PRISM, and SSURGO. Validate each one against StreamStats on ≥10 basins per state
      and document the tolerance.

### 1.3 Regional skew
- [ ] **National regional-skew table.** For each jurisdiction, record the currently adopted
      skew source: the B17C-era Bayesian WLS/GLS study value and MSE, or a documented
      fallback. Do not use the B17B Plate I map silently.
      The table records the value (or map/raster), MSE, effective record length, report
      citation, and validity region.
- [ ] `regional_skew_at(lat, lon)` lookup returning `(skew, mse, citation)` and plugging
      straight into `Bulletin17C`. Raise, rather than default, where no study applies.
- [ ] Tooling to **develop a new regional skew** (B-WLS/B-GLS, Veilleux/Reis-Stedinger) for
      jurisdictions without one. This is lower priority and research-grade.

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
- [ ] **Equation schema** (`flowfreq/regression/`): a typed, serializable definition holding
      region ID, statistic (AEP), functional form (log-linear, power, with transforms such as
      `log10(X+1)`), coefficients, variable definitions and units, calibrated
      min/max per variable, SEP / average variance of prediction, model-error variance, and
      the `(XᵀΛ⁻¹X)⁻¹` matrix where published. Carry the report citation, table number, and
      effective date.
- [ ] **Two evaluation backends** behind one interface:
      1. *Online:* the existing NSS `Scenarios/Estimate` client.
      2. *Offline:* a local evaluator over the equation library, so results are reproducible
         with no network and pinned to a library version.
      A cross-check test asserts that both backends agree to a tight tolerance for every
      equation.
- [ ] **Bootstrapping the offline library from NSS.** `GET /nssservices/regions/{r}/Scenarios`
      and the related equation/citation endpoints expose equations and limits. Write
      `tools/snapshot_nss.py` to dump every region/statistic group to versioned JSON under
      `flowfreq/data/regression/`, with a diff report on re-run so equation changes upstream
      are visible in review.
- [ ] Out-of-range handling: raise by default; `allow_extrapolation=True` records the
      violation (same pattern as `transpose.py`).
- [ ] **Multi-region basins**: area-weighted estimates, following the method each state's
      report prescribes, since some reports weight flows and others weight logs.
- [ ] **Prediction intervals** at ungaged sites from the equation's variance terms
      (Tasker & Driver form). Fall back to average SEP where the covariance matrix is not
      published, and label it as such.
- [ ] **Region-of-influence (ROI)** evaluator for states whose published method is ROI
      rather than fixed equations. This requires the state's gage dataset and the ROI
      parameters from the report.
- [ ] Urban adjustment equations (e.g., Sauer et al. national urban equations, and state
      urban equations where published), applied only where the report says they apply.

### 3.2 Jurisdiction coverage matrix
Build `docs/REGRESSION_COVERAGE.md`, generated by the snapshot tool rather than typed.
It gets one row per jurisdiction, with columns: NSS peak-flow equations present? · report
citation(s) and year · number of regions · ROI? · urban equations? · regional skew source ·
validation worked example located? · offline library validated? · status.

Jurisdictions (56):

- [ ] Alabama · Alaska · Arizona · Arkansas · California · Colorado · Connecticut · Delaware
- [ ] District of Columbia · Florida · Georgia · Hawaii · Idaho · Illinois · Indiana · Iowa
- [ ] Kansas · Kentucky · Louisiana · Maine · Maryland · Massachusetts · Michigan · Minnesota
- [ ] Mississippi · Missouri · Montana · Nebraska · Nevada · New Hampshire · New Jersey
- [ ] New Mexico · New York · North Carolina · North Dakota · Ohio · Oklahoma · Oregon
- [ ] Pennsylvania · Rhode Island · South Carolina · South Dakota · Tennessee · Texas · Utah
- [ ] Vermont · Virginia · Washington · West Virginia · Wisconsin · Wyoming
- [ ] Puerto Rico · U.S. Virgin Islands · Guam · Commonwealth of the Northern Mariana Islands
- [ ] American Samoa

Per jurisdiction, the definition of done:
- [ ] Equations for every published peak AEP (typically 50% down to 0.2%) are in the offline
      library with citation and limits.
- [ ] One or more published worked examples are reproduced to the report's printed precision.
- [ ] Online and offline backends agree.
- [ ] A regional skew is available (§1.3).
- [ ] Region polygons are loaded and point-in-polygon selection is tested on a known site.

**Expected gaps to plan for (verify each one; these are not assumptions to code against):**
- Territories (USVI, Guam, CNMI, American Samoa) may have old, sparse, or no NSS equations.
  Document the fallback for each: nearest applicable equations, index-flood from the few
  gages, or "not supported" with the reason.
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
Recommended approaches, most defensible first:
- [ ] **Scaling / "delta" factors applied to B17C quantiles.** Use published agency change
      factors where they exist (state DOT/DEC future-flow multipliers, NOAA Atlas 15
      future precipitation once released, NCHRP 15-61 guidance, FHWA HEC-17). Record the
      source, emission scenario, and planning horizon. This is the simplest option to
      review and should be the default.
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
- [ ] Output labeling: every future-condition quantile carries scenario, horizon, method, and
      source, and is visually distinct from the regulatory B17C value in plots and reports.

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
| A | §1.1 Water Data API migration, gage catalog; §1.3 regional-skew table | Everything |
| B | §3.1 equation schema + NSS snapshot tool + offline evaluator; §3.2 matrix generated | Regression everywhere |
| C | §1.2 region polygons + watershed polygon; §3.1 multi-region, prediction intervals | Correct region selection |
| D | §2 App. 9 weighting, App. 8 record extension; §5.1 regression-weighted DAR | Full ungaged workflow |
| E | §4 `estimate_at_site`; §7 per-state validation, rolled out state by state | "Supported" status per jurisdiction |
| F | §6.1 stationarity detection + report | Default screening on every analysis |
| G | §6.2 nonstationary models; §6.3 future quantiles (change factors first) | Future-condition products |
| H | §3.3 GLS toolkit, ROI, territories fallbacks, geostatistics | Gap-filling and research |

Phases A–C are the critical path. F can proceed in parallel with any of them because it
depends only on peak data, which already exists.
