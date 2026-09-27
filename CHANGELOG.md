# Changelog

All notable changes to FlowFreq (formerly HydroLib) are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Western Oregon peak-flow regression equations** (`flowfreq/data/regression/OR.json`,
  status `partial`, #38): 21 equations, Regions 1, 2A and 2B at the report's 7 recurrence
  intervals, from USGS SIR 2005-5116 (Cooper, 2005) Tables 10-12 and 15. Double-entered
  against the PDF text and reproducing the report's worked examples and its Appendix D
  regression estimates; Region 1 matches live NSS (GC729). No covariance is published, and
  the 2A/2B elevation transition zone (eq. 8) is not modelled by the schema.
- **Montana peak-flow regression equations** (`flowfreq/data/regression/MT.json`, status
  `verified`, #40): 80 basin-characteristics equations, all 8 hydrologic regions at all
  10 AEPs, from USGS SIR 2015-5019-F ver. 1.1 Table 1-4, Table 3 and Table 1-5
  (covariance). Double-entered against the workbook and PDF text, reproducing the report's
  worked examples and every region's MVP, and matching live NSS at its
  three-significant-figure output for every region and AEP. The Northwest region was fit by
  WLS. The channel-width equations (SIR 2020-5142) are a separate method and out of scope.
- **Washington peak-flow regression equations** (`flowfreq/data/regression/WA.json`, status
  `partial`, #37): 32 equations, regions 1-4 at all 8 AEPs, from USGS SIR 2016-5118 ver. 1.2
  Table 6 and Table 7, with Regions 1-2 covariance at full precision from the report's Flood
  Q Tools workbook. Double-entered against the PDF text, reproducing the report's worked
  example and its Table 8 regression estimates, and matching live NSS at its
  three-significant-figure output for every region and AEP. Region 4 carries no covariance
  because Table 7's Region 4 matrix is inconsistent with the equation's raw-P basis.
- **`flowfreq.donor_similarity`: screen and rank donors on basin similarity, not drainage
  area alone** (#13). Opt-in; no existing number changes and `flowfreq.transpose`'s area
  screen is untouched.
  - `screen_donors`: the area-ratio screen (`DEFAULT_AREA_RATIO_RANGE` by default) plus an
    optional per-predictor ratio band (e.g. `PRECPRIS10`, `CANOPY_PCT`), reporting each
    donor's ratios, a combined `pass`/`fail`/`excluded` status and the criterion that failed.
    No predictor thresholds are shipped; they are local, like the area band.
  - `rank_donors_by_similarity`: region-of-influence distance in standardized
    basin-characteristic space (Burn, 1990; log10 drainage area as in Tasker and others,
    1996), equal weights unless supplied, with per-attribute standardized differences.
  - `similarity_donor_chooser`: plugs the ranking into `qppq.loocv_qppq`'s `donor_chooser`.
  - Missing, non-numeric or unit-mismatched attributes exclude a donor with the reason stated,
    or raise for the target. Nothing is imputed.
- **`flowfreq.psf_convert`: a PeakFQ `.psf` station converted into a Bulletin 17C analysis**
  (#31, the converter half). `station_rows` builds the station's EMA rows exactly as peakfq
  8.1.0's `siteQT` does: codes, `Peak`/`Interval` lines, perception thresholds in file order,
  gap years, and the no-information rule. `convert_station`/`convert_psf` then map the
  settings the way `main.R` does. The skew option keeps peakfq's `r_G_mse` sign encoding
  (weighted `SkewSE**2`, generalized `-SkewSE**2`, station-only sentinel). `LOType`
  MGBT/FIXED/NONE keeps the `gbthrsh0` encoding, along with `WeightOpt`, `ConfInterval`,
  `EXTENDED` AEPs and `O EMA`. Output-only options are logged as ignored.
  `StationInputs.bulletin17c_kwargs(engine)` gives `Bulletin17C` arguments only after
  rebuilding that engine's rows from them and checking the rows equal `siteQT`'s. When
  they don't, it raises `UnsupportedSpecError` and names the years.
  `StationInputs.fortran_reference()` runs any station through `emafitpr` directly.
  Validated on the vendored WY/MT file (24 stations):
  - BegYear/EndYear/HistPeaks/GagedPeaks match peakfq for all 24.
  - Powder River's and Cains Coulee's `emafitpr` rows equal the committed goldens' inputs.
  - Powder River refits to the golden within 1e-6.
  - The `Bulletin17C` Fortran route is bit-identical to a direct `siteQT` call on 23 of 24.
- `flowfreq.watstore`: a reader for legacy WATSTORE (`I ASCI`) peak files, following
  `readWATSTORE`. It returns the standard peak-frame columns plus `site_no`.
- `flowfreq.peak_codes.peak_frame_intervals` applies the `siteQT` code rules to a peak frame,
  and `count_acted_on_codes` counts the codes that change a peak's treatment. `analyze_gage`
  gains an opt-in `apply_peak_codes` (via `psf_convert.convert_peak_frame`). By default it
  now logs which treatment-changing codes it is ignoring. It computes nothing different. (#30)
- `PsfFile.peak_file`, and `O Plot Style`/`O Plot Position` are now read as separate
  options (both were stored under one `Plot` key, so the second overwrote the first).
- **Wave 1 regional skew values** in `flowfreq/data/regional_skew.csv` (#34, #37-#40).
  Each value was transcribed from its primary source and cross-checked against a second
  report. WA, OR and ID are `verified` for the Pacific Northwest B-GLS CONSTANT model
  (A.G. Veilleux): G = -0.07, MSE = 0.18 (SE 0.4243), effective record length 41 years.
  Sources: SIR 2016-5083 app. B, Table B2 p. 53, and SIR 2016-5118 app. A, Table A2 p. 66.
  Two rows stay `pending`, with the reason in their notes:
  - Idaho's Snake River Plain gets its own row. SIR 2016-5083 p. 52 says the PNW model is
    not valid there.
  - MT: SIR 2025-5019 pp. 9-10 says no B-GLS study covers the state. USGS Montana practice
    uses the spatially varying B17B Plate I map, which this table never returns.
- **Idaho peak-flow regression equations** (`flowfreq/data/regression/ID.json`, status
  `verified`, #39): 66 equations, all 6 regions (1_2, 3, 4, 5, 6_8 and 7) at all 11 AEPs,
  from USGS SIR 2016-5083 ver. 1.1 Table 4, Table 5 and Table A5 (covariance),
  double-entered against the PDF text and matching live NSS at its
  three-significant-figure output for every region and AEP. Region 4's
  `(F/100 + 1)^b` and region 7's `(Emin/1,000)^b` terms are stored as published with
  `Variable.scale`/`offset`, so every covariance matrix is Table A5's, unrebased.
  `partial`, #39): 55 equations, regions 1_2, 3, 5, 6_8 and 7 at all 11 AEPs, from USGS
  SIR 2016-5083 ver. 1.1 Table 4, Table 5 and Table A5 (covariance), double-entered
  against the PDF text and matching live NSS at its three-significant-figure output for
  every region and AEP. Region 4 is left out because its `(F/100 + 1)^b` term has no
  schema transform yet.
- **`tools/snapshot_nss.py` and `flowfreq.regression.nss`** (#35). NSS templates carry no
  equations, and neither `apiconfig` nor `RegressionRegions` exposes them. The tool therefore
  fills each region's parameters with in-range midpoints, calls `Scenarios/Estimate`, and
  harvests each equation string, `sep` (stored verbatim as `nss_sep`, never as `sep_log`),
  ASEp, limits and citation into `data/nss_snapshots/<STATE>_<date>.json`. A re-run prints a
  diff against the previous snapshot. A region that Estimate will not answer, even after
  retrying its unbounded gating parameters, is recorded as `unresolved` with the reason.
  `parse_equation` rewrites an NSS string exactly into `RegressionEquation` form: the
  `log10`, `log10_plus1` and `10^(c*X)` identity terms, with scales folded into the
  intercept. An affine base in one variable, such as `(X+0.01)^b`, `(X-20)^b` or
  `(X/100+1)^b`, becomes `log10` with `Variable.scale`/`offset`, stored as written. It
  refuses anything else, such as `max(0, E-1)` or a base in two variables, rather than
  approximating it. `statistic_code_to_aep` maps `PK66_7AEP` to 0.667. The
  initial snapshot covers WA, ID, MT and OR peak flow: 437 equations, 420 parsed, and every
  parsed equation reproduces NSS within NSS's own 3-significant-figure rounding. Oregon
  regions 2A and 2B are unresolved.
- **Water Data OGC API backend for instantaneous values** (`flowfreq.waterdata`, #29).
  `USGSgage.download_instantaneous_flow` / `download_instantaneous_stage` take a new
  `backend=` keyword; `"waterdata-ogc"` reads 00060/00065 from
  `api.waterdata.usgs.gov/ogcapi/v1/collections/continuous` and returns the same frame as
  the legacy path (UTC index, value column, `datetime_local`, `tz_cd`,
  `qualification_code`). The default stays `"nwis-legacy"`, whose behaviour is unchanged.
  A site with several series for the parameter raises `AmbiguousTimeSeriesError` listing
  them, and `ts_id=` takes the 32-hex `time_series_id`; series are never merged, and
  conflicting duplicate timestamps within one series raise. Local time is derived from the
  monitoring location's time-zone fields (IANA zone with DST, fixed offset without). Windows
  are local calendar days, chunked (`chunk_years` at most 3 under the API's 1100-day cap),
  paged, half-open so no boundary instant is counted twice, and default to the series'
  period of record from `time-series-metadata`. `Approved`/`Provisional`/`ESTIMATED` map to
  `A`/`P`/`e`; other qualifier tokens such as `ICE` are kept verbatim. Tests run offline
  against trimmed live captures from 2026-09-25 (`tests/fixtures/waterdata_ogc/`, regenerated
  by `tools/capture_waterdata_fixtures.py`).
- **`waterdata-ogc` peak backend implemented** (`flowfreq.peak_sources.WaterDataApiBackend`),
  built against the live Water Data OGC API `peaks` collection, which was verified on
  2026-09-25. It adds the `USGS-` prefix to the site ID, filters to `parameter_code=00060` so
  gage-height rows do not duplicate water years, and follows `rel=next` paging. It takes
  `water_year` from the API's own field. `qualifiers_to_codes` translates qualifier tokens
  (`HISTORIC`, `LESSTHAN`, …) into the comma-separated NWIS code string that
  `flowfreq.peak_codes` expects, and drops date-precision and gage-height flags. Unknown tokens
  are logged at warning level and are never split into characters. An empty result raises.
  Big Sandy (03606500) and Orestimba (11274500) match their fixtures live. `nwis-legacy`
  remains the default until the parity test passes. (#29)
- **Roadmap scaffolding for Phase A (data foundation) and Wave 1 (Columbia River basin)**,
  per `docs/MASTER_ROADMAP.md`. Types, loaders, validation and tests only: no endpoint was
  live-verified and no published coefficient or skew value was transcribed.
  - `flowfreq.peak_codes`: NWIS peak qualification codes, turned into EMA intervals using
    peakfq 8.1.0's `siteQT` rules (codes 4, 8, 3/O, and 6/C), plus a regulation classifier
    based only on peak codes. It deliberately has no "reference" class. (#30, #32)
  - `flowfreq.psf`: a PeakFQ `.psf` specification reader, checked against the vendored WY/MT
    file and peakfq's own `RegSkew`/`RegMSEG` output. (#31)
  - `flowfreq.peak_sources`: a pluggable peak-data backend. `nwis-legacy` wraps
    `USGSgage`. `waterdata-ogc` raises until it is live-verified. (#29)
  - `flowfreq.catalog` and `tools/build_gage_catalog.py`: the gage catalog schema, a
    loader, and a build tool that needs network access. (#33)
  - `flowfreq.regional_skew` and `flowfreq/data/regional_skew.csv`: a regional skew lookup
    that never returns a `pending` row. Wave 1 rows are pending. (#34)
  - `flowfreq.regression`: the regression-equation schema, an offline evaluator with range
    checks and prediction intervals, and multi-region weighting. It includes the table of all
    56 jurisdictions and their waves, and per-state JSON libraries (Wave 1 files are pending,
    with no equations). (#35)
  - `flowfreq.future_flow`: `ChangeFactorSet`, `apply_change_factors`, and precedence of
    state sets over the national set. It ships no factor values. (#36)
  - `tools/gen_regression_coverage.py` generates `docs/REGRESSION_COVERAGE.md`. A test fails
    if the committed report is stale.
  - `docs/FUTURE_FLOW_GUIDANCE.md`, and a state-rollout issue template.
- **`flowfreq.subdaily`** -- a new module for instantaneous-series metrics, per
  `docs/SUBDAILY_METRICS_DESIGN.md`. Two families of metric plus the circular
  statistics they need:
  - `daily_extreme_timing()` -- hour of the daily maximum *and* minimum, per local
    calendar day, as a fractional local hour. The pair is what separates a snowmelt
    signal (morning minimum, late-afternoon peak) from an operations-driven one;
    either half alone is ambiguous. A flat day reports NaN rather than a fabricated
    midnight peak, and ties report their count rather than silently returning the
    first of several.
  - `circular_hour_statistics()` and `extreme_timing_summary()` -- pooled timing as
    **circular** statistics. An hour-of-day is an angle, not a point on a line: the
    arithmetic mean of hours 23.5 and 0.5 is noon, the one time such a site never
    peaks, and the error is worst for the sites whose timing is most consistent.
    Reports a mean hour beside the concentration `r` and a Rayleigh p-value, so a
    meaningless mean direction announces itself.
  - `ramping_rates()` and `ramping_rate_summary()` -- rate of change between
    consecutive observations, per day, as up-ramp, down-ramp and mean absolute rate
    in **cfs/hr, %/hr and ft/hr**, plus the standard hydropeaking reversal count and
    exceedance counts against a stated licence limit. Down-ramp is signed and never
    positive; `dt` is measured on the UTC axis, never on local wall-clock time (which
    runs backwards across a fall-back transition and turns a modest ramp into a large
    one of the wrong sign); an interval longer than `max_gap_hours` is counted as a
    gap rather than averaged in as a gentle ramp.
  - **Percent-per-hour on a stage series raises.** Gage height is measured from an
    arbitrary local datum, so a "5% rise in stage" has no physical referent while
    looking exactly like a regulatory ramping rate. Same stance as `transpose`'s
    `_require_kind`: refuse rather than return the plausible wrong number.
- **`USGSgage.download_instantaneous_stage()`** -- NWIS parameter 00065 (gage height),
  which is the series most FERC and HCP ramping-rate conditions are written against.
  Shares the whole retrieval path with `download_instantaneous_flow` (chunking, UTC
  axis, tz-abbreviation mapping, multi-sensor `ts_id` disambiguation), now
  parameterized over `IV_PARAMETERS` rather than duplicated. Cached on a separate
  `instantaneous_stage` property so a stage download cannot land where a caller
  expects `flow_cfs`.
- **`Variable.scale` and `Variable.offset`** in the regression-equation schema
  (`flowfreq.regression`). The value entering an equation is now
  `T(scale * x + offset)`, so published terms such as `(FOREST/100 + 1)^b`,
  `(ELEV/1000)^b` and `(GUTTER + 0.1)^b` are stored exactly as printed, with the
  published covariance matrix unchanged. Both fields are optional (defaults 1 and 0),
  so every existing file loads and evaluates identically and `schema_version` stays 1.
  `log10_plus1` remains as an alias for `log10` with offset 1. `minimum`/`maximum`
  still apply to the raw value in published units. A transform argument outside its
  domain (`<= 0` for `log10`) raises a `ValueError` naming the variable, even with
  `allow_extrapolation=True`. `Variable.to_dict()` writes the JSON form and omits
  the new fields at their defaults.

### Changed
- **`peak_sources.DEFAULT_BACKEND` is now `waterdata-ogc`** (the USGS Water Data OGC API),
  replacing `nwis-legacy`, which USGS is retiring. `get_backend()` with no name returns the new
  backend; `get_backend("nwis-legacy")` still works. Switched because the #29 parity test
  (`tests/test_peak_backend_parity.py`) passed live on 03606500, 11274500 and 01638500: same
  water years, flows and discharge codes. One visible difference: `peak_date` is the API's UTC
  date when the time of day is known, so an evening peak can read one day later than legacy's
  local date. Water years are unaffected. `USGSgage.download_peak_flow` is unchanged and still
  reads the legacy service directly.
- **Breaking: no silent regional skew.** `run_ffa`, `compare_engines` and `flowfreq compare`
  used to fall back to a regional skew of -0.302 (SE 0.55) whenever none was given. They now
  raise `ValueError` (the CLI raises a usage error) unless the caller makes exactly one choice:
  - `regional_skew` plus `regional_skew_se`, from a published study;
  - `station_skew_only=True` (`--station-skew`) for the at-site skew alone; or
  - `use_default_skew=True` (`--default-skew`) to accept -0.302 explicitly, which logs a
    warning.

  Why: Bulletin 17C (p. 31) gives no national default, and says the 17B plate 1 estimates
  "are not recommended for use in flood frequency studies". No publication giving -0.302 has
  been found; 0.302 is 0.55², plate 1's MSE. peakfq 8.1.0 likewise stops with an error when a
  weighted skew has no `GenSkew`/`SkewSE`.
  - The new `resolve_regional_skew()` implements the rule. `run_ffa` now records the choice
    under `parameters["regional_skew_source"]` (`user`, `default` or `station`).
  - `regional_skew_se` has no default either, because 0.55 was plate 1's SE.
  - **Migration:** add `use_default_skew=True` to reproduce old results exactly. Better, pass
    a published value; `flowfreq.regional_skew.regional_skew_for` covers the states that have
    one.
  - `Bulletin17C` itself is unchanged; it never defaulted, and `regional_skew=None` has always
    meant station skew there.
- **`diel_variation` and `diel_variation_summary` moved** from `flowfreq.regime` to
  `flowfreq.subdaily`, alongside the other instantaneous-series metrics they share
  their local-day and completeness conventions with. **No import path changed** --
  both are re-exported from `flowfreq.regime` and from the package top level, and
  `tests/test_regime.py` was deliberately left importing them from `regime` as the
  back-compatibility check. The move is to get them type-checked: `pyproject.toml`'s
  mypy override list exempts `flowfreq.regime` as pre-existing debt and states that
  nothing new is added to it, so a new module is the only place new public functions
  are checked at all.

### Fixed
- **Native EMA: at-site skew MSE now follows peakfq's B17B switch.** `emafit.f:707-711` uses
  the plain Bulletin 17B `mseg(n, G)` over the whole record (uncapped), not ADJE's
  censoring-adjusted MSE, whenever MGBT computes the low-outlier threshold and finds low
  outliers, and keeps that choice for the confidence bounds. The native engine always used
  ADJE. It now switches the same way (`ExpectedMomentsAlgorithm._at_site_option`); a
  user-supplied threshold keeps ADJE, as peakfq's FIXED option does.
  **Results change** for sites where MGBT finds low outliers and a regional skew is weighted
  in: the weighted skew, quantiles and confidence bounds move to peakfq's values.

  | Site | Weighted skew vs peakfq, before -> after | Worst quantile error, before -> after |
  |---|---|---|
  | Cains Coulee 06327450 (11 PILFs) | 0.058 -> 6e-6 | 9.06% -> 0.0012% (Q100: 2.06% -> 0.00013%) |
  | 06324500.00 / .01 (17 PILFs) | off -> within 3e-6 | -- |

  Big Sandy, Powder River and 12363000 (no MGBT low outliers) are unchanged. This resolves the
  long-standing Cains Coulee `skew_weighted` xfail and the three that followed from it; the
  "unexplained ~3x `as_G_mse` discrepancy" it was attributed to is this switch (a standalone
  `mseg_all` call runs ADJE). CLAUDE.md's Validation Status is corrected.
- **`USGSgage.download_peak_flow` silently dropped peaks with an unknown day or
  month.** NWIS writes these as `00` (`1897-03-00`), which `pd.to_datetime`
  coerced to NaT, so the row was lost. Those rows are typically the historic
  peaks, which a B17C analysis relies on most. Partial dates now follow peakfq
  8.1.0's reader: an unknown month becomes January, which keeps the water year
  equal to the calendar year, and an unknown day becomes the 1st. **This changes
  results**: sites with partial-date peaks now return more peaks than before.
- **`download_peak_flow` read an all-numeric `peak_cd` column as float**, so
  code `7` came back as `"7.0"`. Code columns are now read as strings.
- **The native EMA ignored perception thresholds that were not entirely before the
  systematic record.** `ExpectedMomentsAlgorithm` honoured only threshold periods ending
  before the first systematic year, and merged several of those into one period at the
  smallest threshold. A period over a gap inside the record was dropped, so those years were
  treated as having no information. It now applies every period to the years it covers, in
  order with the later period winning, as peakfq 8.1.0's `siteQT` does
  (`vendor/peakfqr/R/readInputs.R`). A year with no peak inside a nonzero-threshold period is
  censored below that threshold. An observed peak carries its year's threshold. Two
  `emafit.f` `gbtest` rules that only matter once such rows exist now apply too. The
  low-outlier cutoff applies to every row, not only systematic peaks (lines 1062-1075). A
  censored year whose threshold is no larger than the smallest peak enters MGBT as a
  systematic observation (lines 966-978). `EMAParameters` is unchanged. Its single
  historical period is still reported, but when it was summarised from `perception_thresholds`
  it no longer drives the fit. **This changes results** for any input with a threshold period
  that overlaps or lies inside the record, or with more than one historical period. Big
  Sandy, 12363000, Powder River and Cains Coulee are bit-identical, including confidence
  limits. On the WY/MT stations with in-record thresholds, native against live `emafitpr` on
  the same inputs:

  | Station | Rows (native / Fortran) | Weighted skew, before → after (Fortran) | Q(1% AEP) cfs, before → after (Fortran) |
  |---|---|---|---|
  | 06185500.10 | 74 → 81 / 81 | 1.0372 → 1.0195 (1.0195) | 83,919 → 82,918 (82,607) |
  | 06324500.00 | 143 → 145 / 145 | 0.3575 → 0.3553 (0.4894) | 38,362 → 38,209 (39,877) |
  | 06324500.01 | 145 → 145 / 145 | 0.2485 → 0.2485 (0.3100) | 41,353 → 41,353 (42,290) |
  | 06324710.00 | 17 → 47 / 47 | 0.1943 → 0.1210 (0.1210) | 33,668 → 27,938 (27,934) |
  | 06325500.00 | 25 → 32 / 32 | -0.1897 → -0.2166 (-0.2166) | 4,087 → 3,923 (3,922) |
  | 06327700.00 | 11 → 50 / 50 | 0.2507 → -0.0117 (-0.0117) | 14,256 → 7,849 (7,849) |

  Mean, standard deviation and both skews now agree with `emafitpr` to within 3e-6 on four
  of the six stations. The remaining quantile gap is in the quantile function, not the
  fit. The 06324500 pair still miss on the weighted skew because `emafitb` switches to the
  B17B skew MSE when MGBT finds low outliers (`emafit.f` lines 706-710), which the native
  engine does not do yet. That is a separate defect, marked `xfail(strict=True)` in
  `tests/fortran_parity/test_live_wymt_perception_thresholds.py`.
- **LP3 quantiles used the Wilson-Hilferty approximation instead of the exact Pearson III
  inverse.** `flowfreq.core.kfactor` computed the frequency factor K with Wilson-Hilferty at
  every skew. peakfq 8.1.0 computes its quantiles with `emafit.f`'s `qP3` (line 3266). That
  inverts the incomplete gamma function exactly and uses Wilson-Hilferty only for
  |skew| < 0.001, where the two agree. The approximation's error grows with |skew| and
  return period, so native quantiles drifted from peakfq's even where the fitted moments were
  identical. `kfactor` now calls `flowfreq._p3_moments.q_p3`, the existing port of `qP3`
  already checked against the Fortran. **This changes results**: every quantile, and every K
  factor, from `ExpectedMomentsAlgorithm`, `MethodOfMoments`, `Bulletin17C`, the frequency
  plots and anything else built on `kfactor`/`kfactor_array`. The change grows with |skew|
  and with smaller AEP. For example, at skew 0.87 and AEP 0.002 the flow is 1.3% lower. EMA
  confidence bounds are unchanged; they already came from `qP3` inside `var_emab`, so the
  table's central `flow_cfs` now matches the quantile those bounds were built around.
  Largest quantile error against peakfq 8.1.0 over each site's AEPs (0.995 to 0.002):

  | Site | Skew used | Worst error, before → after | Q(0.2% AEP) cfs, before → after (peakfq) |
  |---|---|---|---|
  | Big Sandy 03606500 | -0.156 | 0.056% → 0.0002% | 31,653 → 31,636 (31,636) |
  | 12363000 | 0.286 | 0.106% → 1.4e-6% | 154,279 → 154,116 (154,116) |
  | Powder River 06326500 | -0.184 | 0.100% → 2.7e-7% | 62,145 → 62,083 (62,083) |
  | Cains Coulee 06327450 | -0.662 | 9.74% → 9.06% | 2,853 → 2,822 (2,928) |
  | 06185500.11 (live `emafitpr`) | 0.871 | 1.35% → 3e-7% | 160,844 → 158,706 (158,706) |
  | 06329200.00 (live `emafitpr`) | -0.445 | 0.96% → 1e-5% | 10,358 → 10,260 (10,260) |

  The residual Big Sandy error is the moments' own 2e-6 gap, not the quantile function. On
  06329200.00 the 1e-5% (4e-8 in log10, at AEP 0.995) is the Fortran's own gamma inverse:
  `q_p3` and `emafit.f`'s `qp3sub` differ by exactly that on identical moments.
  Cains Coulee's error is its known weighted-skew gap, which Wilson-Hilferty had been partly
  hiding. Its Q100 error moves from 1.53% to 2.06%, so the recorded "< 2%" bound becomes a
  new `xfail(strict=True)` (`test_q100_error_under_two_percent`) with that reason, not a
  wider tolerance.
- `peak_codes.parse_codes` read the date-precision codes `Bd`/`Bm` one character at a
  time, so the legacy NWIS value `"7,Bd"` (Big Sandy 03606500, WY1897) parsed as
  {7, B, D}, and `"Bm"` (Orestimba 11274500, WY1947) as {B, M} with an "unrecognised code
  M" warning. They are now whole codes, and neither changes a peak's EMA interval. A
  float-formatted code (`"7.0"`, from a code column pandas inferred as float) is read as
  code 7 rather than {7, 0}. No numeric result changes: D and M were never acted on.

## [0.8.0]

### Added
- **`flowfreq.streamstats`** (Phase 1) -- USGS StreamStats watershed delineation and
  basin-characteristics retrieval for a pour point: `pourpoint` snap, `ss-delineate`
  `delineate/sshydro`, `ss-hydro`, per `docs/STREAMSTATS_MODULE_DESIGN.md`. Validates
  every response against what it's supposed to contain rather than trusting a 200 (an
  unsnappable point returns HTTP 200 with a `WarningMsg` string rather than an error);
  typed/indexable characteristics with provenance; an offline-capable cache; a batch
  entry point that returns results and failures side by side rather than aborting on
  one bad point. Confirmed against the live service, reproducing the design doc's
  published values exactly. Flow-statistics regression estimation (NSS) is deliberately
  out of scope for this phase, and so is the watershed polygon itself --
  `WatershedCharacteristics.polygon_geojson` is always `None`; no call in the verified
  protocol returns it, and getting it is an open question for later.
- **`flowfreq.streamstats`** (Phase 2) -- NSS (National Streamflow Statistics)
  flow-statistic estimation from real basin characteristics, per
  `docs/STREAMSTATS_NSS_ADDENDUM.md`: `list_statistic_groups()`,
  `estimate_flow_statistics()`, `batch_estimate_flow_statistics()`. Every estimate
  carries NSS's own regression-equation string and a resolved citation
  (title/author/DOI). Validates every parameter against the region's own valid range
  *before* submission -- confirmed live that NSS returns a plausible, silently-wrong
  extrapolated value for an out-of-range input with no warning at all. Region
  selection across NSS's several regressionRegions per state is not automatic (no
  available call filters by location without a watershed polygon, which Phase 1 does
  not produce) -- every geographically-plausible region is returned, unlabelled as to
  which is geographically correct; picking it is the caller's responsibility.

### Fixed
- **`USGSgage.download_daily_flow`**: the timeout was hardcoded at 30s even after an
  earlier fix made a full period-of-record request (tens of thousands of rows for a
  long-running site) the default -- now configurable, default 60s. The default end date
  used local wall-clock time rather than UTC, so a host clock behind UTC could silently
  request a narrower range than intended; now computed in UTC. `start_date`/`end_date`
  are now validated (parseable, `start_date` not after `end_date`) before any request is
  sent, rather than reaching NWIS unvalidated.

## [0.7.0]

### Added
- **`flowfreq.transpose`** -- moving a computed statistic from a gaged donor basin to a
  nearby ungaged target by drainage-area ratio, `Q_t(p) = Q_d(p) * (A_t/A_d)**b(p)`, with
  `b(p)` from the applicable published regression rather than assumed to be 1. Three
  functions, deliberately separate: `transpose_frequency` (floods),
  `transpose_duration` (flow-duration statistics) and `transpose_low_flow`.

  The failure mode this is built around is silent: a wrong exponent does not raise, it
  returns a plausible, monotone, right-order-of-magnitude discharge that is wrong by
  10-30%. So `RegressionExponents` will not construct without a citation, every quantile
  records whether its exponent was published, interpolated or extrapolated, and
  `TransposedResults` deliberately carries no LP3 moments -- no distribution was fitted at
  the target, and reporting the donor's would invent a fit nobody performed.

  `probability_kind` ("aep" / "exceedance" / "non_exceedance") makes the category error
  unreachable rather than merely documented: a flood regression's exponent describes flood
  response and says nothing about 7Q10, so handing one to `transpose_low_flow` raises.
  Separate functions were not enough on their own, because the *arguments* were still
  interchangeable.

  Transposed curves are checked for monotonicity in both directions and raise if they
  invert. That was a real defect found during development: with a groundwater-dominated
  donor (flat dry end) and an exponent set falling toward the dry end -- the combination
  the duration literature actually describes -- `transpose_duration` could return
  `Q99 > Q95` at area ratios inside the supported band.

- **`flowfreq.qppq`** -- QPPQ daily-series transfer through two flow-duration curves, with
  an invertible empirical `FlowDurationCurve`, seasonal curve construction, donor ranking,
  goodness-of-fit and a leave-one-out harness.

  Includes a measured result that contradicts the obvious fix. QPPQ assumes donor and
  target sit at the same position in their own curves on the same day; in a snowmelt basin
  a higher target melts later, breaking that systematically. Grouping the curves by season
  is the intuitive remedy and *makes it worse*: on a 25-day imposed offset, log-NSE was
  0.273 annual, 0.040 for three seasons, 0.605 monthly, and 0.895 once the donor series was
  lagged. A coarse season is longer than the offset being corrected. `estimate_donor_lag`
  recovers the offset from rank correlation; `center_of_timing` gives the route for a
  target with no record, since melt timing regresses on basin elevation.

  `performance()` reports NSE, log-NSE, KGE and dry-end bias together, because an estimate
  that is exact through the freshet and wrong by 10x in September still scores NSE > 0.99.

- **`regime.flow_duration_curve`** -- duration statistics were previously reachable only as
  a by-product of drawing a figure, at nine hardcoded percentiles.
  `Hydrograph.plot_flow_duration_curve` now delegates to it, so the table it has always
  returned is one computation rather than two that can drift.

### Fixed
- `fortran_engine.build_emafit_arrays` reused the scalar loop-local names `ql`/`qu`/`tl`
  for the arrays of the same name later in the function -- harmless at runtime, a real
  mypy error.

## [0.6.1]

### Added
- `plot_peak_flows_with_thresholds` gained a `yscale: str = "log"` parameter (same convention
  `plot_frequency_curve` already uses), the last gap against `flowfreq-app`'s
  `plot_peak_timeseries` -- a per-plot linear/log toggle the app exposes as a real sidebar
  control. Verified live that switching the app over loses nothing numerically: the
  quantile-line values this function computes analytically from `lp3_params` match
  `run_ffa`'s actual `quantile_df["Flow (cfs)"]` to 0.0% (uncensored and censored sites both
  checked, including the app's full return-period list), and `core.log_pearson3_cdf` matches
  the app's own max-peak-recurrence formula to ~1e-15.

## [0.6.0]

### Added
- `plot_peak_flows_with_thresholds`'s `mgbt_threshold` now also draws peaks below the cut as
  hollow outline bars, matching what `flowfreq-app`'s own `plot_peak_timeseries` did that this
  function previously didn't -- the third and last of the three features the app carried
  (0.5.0 moved the other two, return-period lines and the max-peak annotation). A new
  `mgbt_threshold_source` argument labels the threshold line's legend entry (`"override"` vs.
  the default `"MGBT"`), matching the app's own PILF-source distinction. This closes the plot
  dedupe: `flowfreq-app` can now switch to this function and delete its own copy entirely.

### Fixed
- `fortran_engine.build_emafit_arrays` reused the scalar loop-local names `ql`/`qu`/`tl` for
  the final `ql`/`qu`/`tl` arrays later in the same function -- harmless at runtime (Python
  doesn't care), but a real mypy violation (`error: Incompatible types in assignment`).
  Renamed the loop-locals to `row_ql`/`row_qu`/`row_tl`.
- **`make clean` never wiped `.mypy_cache`**, which is how the error above shipped in 0.5.0
  despite `clean-verify` passing repeatedly: mypy's incremental cache can mask a real error on
  a file it already has a (stale, error-free) entry for, the exact failure mode
  `clean-verify` exists to rule out for the rest of the tree. `clean` now removes it too.

## [0.5.0]

### Added
- **The vendored Fortran as a selectable analysis engine.** `Bulletin17C.run_analysis(engine=)`
  accepts `"fortran"` alongside the default `"native"` (which stays the default forever); the
  real feature is `flowfreq.workflow.compare_engines`, which runs both and returns an
  `EngineComparisonReport` with `.max_quantile_deviation_pct` and `.to_markdown()`, and the
  `flowfreq compare` CLI subcommand built on it. Both require the built f2py extension and
  raise the library's existing, actionable `ImportError` rather than falling back to committed
  golden files -- a comparison meant for a LOMR/CLOMR submittal has to mean "measured against
  PeakFQ 8.1.0 just now," not "replayed from a file in this repository" (`docs/
  FORTRAN_ENGINE_DESIGN.md` section 9). New module `flowfreq/fortran_engine.py`: a library-side
  interval builder translating a `Bulletin17C` input set into `emafitpr`'s `ql/qu/tl/tu/dtype`
  arrays, following `siteQT` (`vendor/peakfqr/R/readInputs.R`) rather than the parity suite's
  test-only builder. Handles gap years with and without a declared perception threshold (the
  worked example in the design doc: site 12363000's at-site skew is +0.435 omitting the four
  unmeasured years, +0.250 censoring them at the lower threshold -- both reproduced exactly
  through the new builder), historic-peak flagging restricted to the USGS historic flag rather
  than the whole historical period, zero flows floored at `Qmin = 1e-20` (peakfq has no special
  zero-flow case), a user-supplied PILF threshold, and overlapping perception-threshold periods
  (last-declared wins, per `siteQT`'s own documented priority rule). The adapter from
  `ReferenceResult` to `FrequencyResults` never synthesizes a field the Fortran did not report --
  `ema_iterations`/`ema_converged` are `None`, not a guess. Verified live against the built
  extension on all four parity sites: Big Sandy 0.0587%, Powder River 0.1003%, site 12363000
  0.1057% (matching the design doc's own worked example), all under tolerance; Cains Coulee
  9.741%, correctly surfaced as a failure -- it is the site carrying the pre-existing, still-open
  `skew_weighted` residual (see P3 below), not a defect in this feature.
- Tests for the four modules that had none: `hydrograph.py`, `plots.py`, `batch.py`, and
  `cli.py` (the last covering `validate`/`benchmark` and `compare`'s extension-agnostic paths;
  `compare`'s Fortran-backed happy path lives in `tests/fortran_parity/test_live_cli_compare.py`
  instead). Writing the `batch.py` tests surfaced a real, pre-existing defect, not introduced
  here and not fixed by this release: `batch.run_multi_site` passes NWIS records straight from
  `usgs.fetch_nwis_batch` (plain dicts) into `B17CEngine.fit`, which reads `.flow` as an
  attribute -- every real multi-site call has therefore been failing per-site, silently, caught
  by a broad `except Exception` and turned into `{"error": ...}`. Pinned as
  `tests/test_batch.py::TestAnalyzeSites::test_real_fetch_output_shape_is_analyzable`,
  `xfail(strict=True)`, so it stops the build the moment someone's fix makes it pass.
- `plot_peak_flows_with_thresholds` gained opt-in return-period reference lines and a max-peak
  recurrence annotation via a new `lp3_params=(mean_log, std_log, skew)` argument -- the
  library half of retiring `flowfreq-app`'s own duplicate plotting code. The app side (deleting
  its local copy) is a separate, later change gated on a release and a pin bump.

## [0.4.0]

### Added
- mypy runs in CI, enforced on the modules that already pass; the rest are
  exempted individually in `pyproject.toml` so the debt is countable and shrinks
  by deleting a stanza. The library ships `py.typed`, so downstream checkers
  trust these annotations -- nothing verified them before.
- `make cov` and a CI coverage step. `pytest-cov` had been a declared dev
  dependency that nothing invoked.
- `make clean-verify`, which wipes build artifacts and every `__pycache__`
  before running the gate, so a green result reflects the tree rather than
  whatever was left lying around.

### Changed
- **`B17CEngine.fit` now uses the Bulletin 17C Eq. 7-2 station skew. The
  numbers this public API returns have moved.** It computed
  `((x - mean)**3).mean() / std**3`, the biased population coefficient, while
  `Bulletin17C.run_analysis` in this same library used the unbiased sample
  estimator `n * sum((x - mean)**3) / ((n-1)(n-2) * std**3)`. The two differ by
  `n**2 / ((n-1)(n-2))` -- 7.2% at n=44, 39% at n=10, and short records are
  ordinary in flood frequency work.

  Measured on Big Sandy (n=44): station skew moves from -0.1748 to -0.1874,
  which is now exactly what the Bulletin 17C path reports for the same record.
  Quantiles move **Q2 +0.13%, Q10 -0.10%, Q100 -0.57%, Q500 -0.93%**. Anything
  derived from `B17CEngine` moves with it, including `batch.batch_summary_table`
  and `plots`. If you have reported a discharge from this class, it will not
  reproduce under 0.4.0 -- pin `v0.3.0` if you need the old figures, and expect
  the 0.4.0 value to be the defensible one.

  `Bulletin17C` is unaffected: it was always correct, and the release exists to
  make the two agree. Recorded in 0.3.0's test suite as a strict xfail; the
  three tests that replace it in `tests/test_engine.py` now guard against a
  revert, one of them naming the old estimator explicitly.

  One deliberate difference remains: `Bulletin17C` clips the station skew to
  ±3.0 (`MAX_ABS_SKEW`) and `B17CEngine` does not, so the two can still diverge
  on a record with extreme skew. That is a separate question from the estimator
  and was left alone.

- `flowfreq.freq_plot.plot_frequency_curve_streamlit` is now
  `plot_frequency_curve`. The old name remains as an alias, so pinned consumers
  keep working; it can go once none use it. The module imports matplotlib and
  returns a `Figure` -- the suffix was always a misnomer and became misleading
  once the app moved to its own repository.

### Fixed
- Four public signatures annotated names their modules never imported, so
  `typing.get_type_hints()` raised `NameError` on `engine.B17CEngine
  .frequency_table`, `batch.batch_summary_table`,
  `freq_plot.plot_peak_flows_with_thresholds` and `Bulletin17C.validate`.
  `from __future__ import annotations` kept it from raising at import, which is
  why it went unnoticed. The first three now resolve; `validate` keeps a
  `TYPE_CHECKING` import to avoid inverting the package's layering, and says so.

- `import flowfreq.peakfqr` without the f2py extension built raised a bare
  `ModuleNotFoundError` naming a private submodule. It now explains that the
  extension is built on demand, gives the command and the toolchain, and says
  that nothing else in the library depends on it -- the native EMA is the
  default path. Its docstring also cited `_shared/peakfqr/src/emafit.f`, a path
  that does not exist here; the sources are under `vendor/peakfqr/`.

## [0.3.0]

### Changed
- **Split into two repositories.** This repo is the analysis library; the Streamlit
  application moved to [pinhead001/flowfreq-app](https://github.com/pinhead001/flowfreq-app),
  which installs this library as a pinned dependency. `app/` and its three test
  modules are gone from here, along with the `smoke` make target and the CI job
  that ran it.
- `flowfreq.workflow` — new module holding the high-level entry points that used to
  live in the app: `run_ffa`, `compute_skew_tables`, `build_skew_curves_dict`. The
  display formatters stayed with the app.
- The gage attributes table moved into the package at `flowfreq/data/`, replacing a
  `package-data` entry that reached outside the package and only ever resolved in a
  source checkout.
- **Renamed to `flowfreq`** — package, import name, distribution and display name.
  `from hydrolib.core import kfactor` is now `from flowfreq.core import kfactor`;
  the console script is `flowfreq`. Historical entries below keep the old name,
  since they describe what shipped at the time.

  The old name was unusable. Deltares publishes `hydrolib` and `hydrolib-core` on
  PyPI, both installing a top-level `hydrolib/` package, and `hydrolib-core` ships
  `hydrolib/core/` against this project's `hydrolib/core.py`. Installed together
  one silently destroys the other -- in one order this library's entire API
  disappears, in the other neither package imports at all -- and `pip check`
  reports nothing wrong. `flowfreq` names what the library does (flood *and*
  low-flow frequency) and cannot be mistaken for HYDROLIB.

### Added
- Native Python port of `var_mom` and its dependency tree (`mn2mvarb`/`mse_ema`, `detrat`,
  `VAR_EMAB`/`regmoms`/`ci_ema_m3b`), verified routine-by-routine against the vendored
  peakfq 8.1.0 Fortran. See TODO.md's P3 section for the full account.
- `MethodOfMoments` now applies the Bulletin 17B conditional-probability adjustment: a
  low-outlier (PILF) threshold — Grubbs-Beck or user-supplied — censors the fit instead of
  only being reported.

### Changed
- `ExpectedMomentsAlgorithm` confidence intervals are now Cohn's asymmetric bounds
  (`hydrolib._var_emab.var_emab`) instead of the symmetric `log_Q ± z*se` approximation.
- `ExpectedMomentsAlgorithm`'s at-site EMA moment iteration on censored intervals now uses
  the Fortran-verified truncated-moment code (`hydrolib._p3_moments.m_p3`) and the correct
  bias-correction sample size, closing a real accuracy gap on any record with censored
  intervals (Big Sandy's historical gap years included, not just MGBT-flagged PILFs).
- Regional skew weighting now includes ADJE's censoring bias adjustment and the Halloween
  determinant ratio (`detrat`), matching peakfq's default `at_site_option`.

### Fixed
- `hydrolib/peakfqsa/` (a subprocess wrapper around a PeakfqSA binary that does not exist)
  removed; it was mock-tested only. `hydrolib/validation/reference.py` covers what it
  contributed, pointed at references that actually exist.
- Bare `except:` in `usgs.py` narrowed to the actual failure modes.
- `analyze_gage()` no longer prints unconditionally to stdout; uses `logging` like the rest
  of the library.

### Removed
- `hydrolib/peakfqsa/` and its mock-only test suite (see Fixed, above).

---

## [0.2.0] - 2026-08-31

### Added
- Instantaneous (unit-value) flow retrieval from USGS NWIS
- Low-flow frequency analysis module (`hydrolib.lowflow`)
  - Annual n-day low-flow frequency with LP3 or lognormal distribution
  - Climatic/water/calendar year definitions
  - Zero-flow-year handling
  - Analytic and bootstrap confidence intervals
- Flow regime metrics module (`hydrolib.regime`)
  - Richards-Baker flashiness index
  - TQmean metric
  - Baseflow separation (UKIH, Lyne-Hollick, HYSEP variants)
  - Monthly and seasonal flow summaries
- Diel (sub-daily) variation analysis
  - Within-day flow range and coefficient of variation
  - Timezone-correct local-day grouping
- Flow series I/O with Parquet backend (`hydrolib.flowio`)
  - Save/load for daily and instantaneous flow data

### Changed
- Improved EMA algorithm convergence handling for edge cases
- Documentation vignettes reorganized (Low-Flow & Flow Regime guide)

### Fixed
- MGBT outlier detection edge case with small sample sizes
- Flow duration curve calculation precision

---

## [0.1.0] - 2026-01-28

### Added
- **USGS Data Retrieval** — Download mean daily, annual peak, and instantaneous flow from NWIS
- **Bulletin 17C Analysis**
  - Expected Moments Algorithm (EMA) — USGS standard method
  - Method of Moments (MOM) fallback
  - Weighted regional skew (MSE weighting per B17C Appendix 6)
  - Multiple Grubbs-Beck test (MGBT) for low outlier detection
  - 90% confidence intervals (5%/95% limits)
- **Hydrograph Plotting**
  - Daily time series plots
  - Summary hydrographs (day of water year with percentile bands)
  - Flow duration curves
- **Frequency Curve Plotting**
  - Log-probability axis
  - LP3 fitted curve with confidence interval band
  - Multi-skew overlay (station / weighted / regional)
- **Streamlit Web Application**
  - Interactive single/multi-gage analysis
  - Regional skew input controls
  - ZIP export (PNG plots, CSV data, LP3 parameters)
  - Multi-gage comparison tables
- **CLI Tools**
  - `hydrolib validate` — EMA validation against reference fixtures
  - `hydrolib benchmark` — Numerical benchmarking (text/JSON output)
- **Technical Reports** — Automated Markdown report generation
- **Validation Framework** — Parity testing against USGS Fortran reference implementation

### Fixed
- Initial release

---

## Notes on Versioning

- **0.x.x** — Pre-release. API may change without warning.
- **1.0.0** — Stable API. Breaking changes require major version bump.

For upgrade guidance, see the [migration guides](docs/) directory.
