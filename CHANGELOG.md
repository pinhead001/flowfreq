# Changelog

All notable changes to FlowFreq (formerly HydroLib) are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
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

### Changed
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
