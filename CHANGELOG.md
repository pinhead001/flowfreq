# Changelog

All notable changes to FlowFreq (formerly HydroLib) are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
