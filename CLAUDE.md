# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**FlowFreq** is a Python library for hydrologic frequency analysis with USGS data retrieval and
Bulletin 17C flood frequency analysis. It has a native Python EMA implementation, validated
against the USGS Fortran reference that is vendored into this repository.

**Reference implementation:** `vendor/peakfqr/` — the USGS `peakfq` R package v8.1.0 (CC0),
containing the authoritative Fortran EMA code. Read `vendor/peakfqr/src/emafit.f` for the
computation and `vendor/peakfqr/R/fortranWrappers.R` for call conventions. It is a verbatim
reference copy: **do not edit anything under `vendor/`**, since a change there silently
invalidates every comparison made against it.

There is no PeakfqSA binary. Earlier work assumed one existed and built `flowfreq/peakfqsa/`
as a subprocess wrapper around it; that premise was wrong, and the subsystem has been removed
— it was mock-tested only and its `requires_peakfqsa` tests collected nothing. What it
contributed to validation now lives in `flowfreq/validation/reference.py`, pointed at
references that exist. The f2py bridge in `flowfreq/peakfqr/` is the only working route to
the reference implementation. `AGENT_BUILD_INSTRUCTIONS_Claude.md` is the original build
brief and explains that history.

## Repository Layout

Everything lives in this repository; there is no external workspace to consult.

```
flowfreq/            the library (flat package, NOT src/flowfreq/)
├── peakfqr/         f2py bridge to the vendored Fortran (built, gitignored)
└── validation/      comparison engine, benchmarks, reports, reference loaders
vendor/peakfqr/      USGS peakfq 8.1.0 reference — Fortran, R, test data (do not edit)
build_fortran/       f2py build script and .pyf signature file
tools/               gen_fortran_golden.py — regenerates parity golden files
tests/fixtures/      shared fixture data (Big Sandy, HU02, Wyoming/Montana, peakfqr respec)
tests/fortran_parity/  native-vs-Fortran parity suite and committed golden files
```

The Streamlit application is no longer in this repository. It lives in
[pinhead001/flowfreq-app](https://github.com/pinhead001/flowfreq-app) and consumes
this library as a pinned dependency, so a change here reaches it only when that
repo bumps its pin.

## Build & Development Commands

```bash
pip install -e ".[dev]"

# Full suite. ~70 s, measured 2026-10-03 on a Windows developer machine (CPython 3.12,
# Fortran extension absent; 2541 tests), and ~355 s with the extension present.
# That extra ~285 s is almost all live emafitpr calls in the vendored Fortran;
# test_fixed_threshold_live.py alone takes ~190 s (twelve ~13 s calls in
# test_regional_weighting_inputs_match). The native code accounts for very little of it.
# History, same machine: ~200 s at v0.9.0 (1970 tests), then ~340 s (absent) / ~880 s (present) at
# 0.10.0, because #81 made the EMA fixed point iterate as p3est_ema does (Big Sandy 13 -> 45
# iterations, a 6000 cfs FIXED threshold 100 -> 827). Each iteration costs one m_p3 per
# censored group, and m_p3 is a 50-digit mpmath incomplete gamma. The perf fix that
# followed changed no results and no iteration counts. The incomplete gamma runs its loops in
# fixed-point integers and gets P(a+k, x) by recurrence; var_emab's regmoms grids are cached
# apart from the AEPs; MGBT's integrand skips scipy.stats' wrappers. Profile before
# assuming the fixed point is the problem again: it is cheap per iteration now, and
# _lower_gamma_reg / var_emab are where a cold fit's time goes. A single cold
# run_analysis() with a regional skew supplied costs ~0.5-1.5 s.
pytest tests/ -v

# What CI actually runs -- the same thing. The marker selection lives in
# pyproject.toml's addopts, so a bare `pytest` is the CI selection.
pytest tests/

# Build the Fortran extension and check the golden files against it. CI runs
# this too, in its own job (gfortran + meson). The parity tests read committed
# golden files and pass without it.
make parity

# Regenerate golden files after changing anything under vendor/peakfqr
python tools/gen_fortran_golden.py

black flowfreq/ tests/ && isort flowfreq/ tests/

# mypy over flowfreq/. A separate CI job (`ci / lint`) from the black/isort
# check, and easy to miss: this command list ran black+isort+pytest for a
# while with no mention of it, which is exactly how v0.7.0 shipped ten mypy
# errors in transpose.py/qppq.py straight to main without anyone running it
# locally first -- pytest was, and stayed, fully green throughout. The
# pyproject.toml [tool.mypy] override list documents which modules are
# exempted (legacy debt) versus enforced (new modules by default, including
# any you add).
mypy flowfreq/

# make check == lint (black --check + isort --check-only) + typecheck (mypy) + test,
# in CI's exact order -- the one command that actually reproduces the `ci / lint`
# and `ci / test` jobs together.
make check
```

**Reproduce CI faithfully.** Two traps have both cost a red build:

- `python -m pytest` puts the working directory on `sys.path`; CI runs the `pytest` console
  script, which does not. The `test`, `cov`, `test-all` and `parity` targets set
  `PYTHONSAFEPATH=1` for you; invoking pytest by hand, set it yourself.
- CI tests Python 3.11–3.14 (see Python Compatibility below). A green run on one interpreter
  says nothing about the others -- fixture plumbing and dependency behavior can both diverge
  across the matrix, which is why all four are tested rather than assumed from one.

**Prefer the `make` targets to the bare commands above**, and on Windows treat them as the
only supported route. Four things differ there. The last three the Makefile handles; the
first is on you at every invocation:

- **Run make with the virtualenv on PATH** — activate it, or go through a runner that does
  (`uv run make check`). The Makefile invokes `python`, and a bare `python` in `cmd.exe`
  resolves to the Microsoft Store shim rather than the venv's interpreter when nothing has
  put the venv on PATH. The failure names the wrong problem — `No module named black`
  against a `WindowsApps\...` path — which reads as a missing package but is the wrong
  interpreter, so reinstalling the package cannot help. Passing `PYTHON=<venv>\python.exe`
  covers every target that only needs the interpreter, but not `fortran` or `parity`:
  `build_fortran/build.py` locates meson with `shutil.which`, a PATH lookup, and meson is
  normally pip-installed into the venv — so with the venv off PATH it exits
  `ERROR: meson not found on PATH` however correct the interpreter is.
- `make` is not installed by default (`winget install ezwinports.make`), and its recipes run
  through `cmd.exe`, which cannot parse Unix inline env-var syntax. Both variables are
  therefore set with `export`, which make applies itself rather than handing to a shell.
- `PYTHONUTF8=1` is exported globally. Eleven files under `flowfreq/` and `tests/` contain
  non-ASCII characters; without UTF-8 mode isort cannot encode them to cp1252, reports
  "Unable to parse file", and **skips the file** — `isort --check` then passes without
  having checked it. Linux CI defaults to UTF-8 and never sees this. Note that a
  `Skipped N files` line in isort's output is unrelated: that is its default skip list
  catching `.venv` and similar, and it appears with or without UTF-8 mode. The signal is
  the "Unable to parse" warning, not the count.
- Install the dev tooling with `pip install -e ".[dev]"`, never by naming packages. The
  extra pins `black>=24.0,<25`; black 26 formats differently, so an unpinned install has you
  "fix" correctly formatted files and turn CI red.

## Python Compatibility

flowfreq supports Python 3.11–3.14.

Python 3.9 and 3.10 are not supported.

Do not introduce syntax or standard-library features that prevent
support for Python 3.11 unless explicitly approved.

When changing dependencies, verify compatibility across all supported
Python versions.

CI should test the supported Python matrix.

## Architecture

- **`core.py`** — data structures (`PeakRecord`, `FlowInterval`, `EMAParameters`,
  `FrequencyResults`, `LowFlowResults`), enums, LP3 utilities (`kfactor`,
  `grubbs_beck_critical_value`, `log_pearson3_*`, `compute_ci_lp3`)
- **`bulletin17c.py`** — `FloodFrequencyAnalysis` (ABC), `MethodOfMoments`,
  `ExpectedMomentsAlgorithm`, `Bulletin17C` (facade), and the MGBT three-sweep
- **`fortran_engine.py`** — the vendored `emafitpr` as a selectable engine (`engine="fortran"`)
- **`workflow.py`** — `run_ffa` / `compare_engines` entry points, `peak_code_kwargs`
- **`usgs.py`** — `USGSgage` retrieval (peaks, daily, instantaneous, paired flow+stage),
  `GageAttributes`
- **`waterdata.py`** — USGS Water Data OGC API client: `request` (API key, 429/503 backoff),
  instantaneous and daily values
- **`peak_sources.py`** — peak-data backends (`waterdata-ogc` default, `nwis-legacy`)
- **`peak_codes.py`** — NWIS peak qualification codes → B17C treatment (peakfq's `siteQT`)
- **`psf.py`** / **`psf_convert.py`** / **`watstore.py`** — PeakFQ `.psf` reader, `.psf` →
  `Bulletin17C` arguments, legacy WATSTORE peak-file reader
- **`regulation.py`** — regulation/urbanization screen (GAGES-II + peak code 6) and the
  regulated-record refusal
- **`catalog.py`** — national gage catalog (`data/gage_catalog.csv.gz`)
- **`regional_skew.py`** — verified regional skew by state, HUC or location;
  **`skew_study.py`** — B-WLS/B-GLS regional skew development
- **`regression/`** — offline regional regression equations: `equations` (schema,
  evaluator, prediction intervals), `library` (per-state JSON, status definitions),
  `jurisdictions` (the 56 jurisdictions and their waves), `nss` (NSS equation parser),
  `oregon` / `montana` (state-specific procedures)
- **`streamstats.py`** — StreamStats delineation, basin characteristics, NSS estimates
- **`transpose.py`** / **`qppq.py`** / **`donor_similarity.py`** — transposition to ungaged
  sites, QPPQ daily transfer, donor screening and ranking by basin similarity
- **`future_flow.py`** — future-condition change-factor framework
- **`lowflow.py`** / **`regime.py`** / **`subdaily.py`** — low-flow frequency, flow-regime and
  sub-daily metrics
- **`engine.py`**, **`batch.py`**, **`report.py`**, **`hydrograph.py`**, **`plots.py`**,
  **`freq_plot.py`**, **`flowio.py`**, **`cli.py`**

## Data Sources

The USGS Water Data OGC API (`api.waterdata.usgs.gov`) is the default backend for annual
peaks (#63/#70), instantaneous values (#77) and daily values (#79), each switched behind a
live parity test against legacy NWIS. Every call goes through `waterdata.request` (#84): it
sends an `X-Api-Key` from `USGS_API_KEY` or `waterdata.set_api_key()`, and retries HTTP
429/503 with `Retry-After` or exponential backoff. Without a key the API allows 1000
requests/hour per IP. Legacy NWIS remains available with `backend="nwis-legacy"`. On the OGC
backend an IV `ts_id` is the 32-hex `time_series_id`, not the NWIS DD number, and
`peak_date` is a UTC date.

## Conventions

- Run `black` and `isort` before every commit; CI lints `flowfreq/` and `tests/` only
- Type hints on all signatures; NumPy-format docstrings on public API
- No bare `except:`; no `print()` in library code — use `logging.getLogger(__name__)`
- Tests for every public function (happy path + one error case minimum)
- All data in log10 space when interfacing with Fortran conventions
- **Known-failing tests are `xfail(strict=True)`, never skipped or hidden behind a widened
  tolerance.** Strict means the build fails the moment one starts passing, which is the alarm
  you want. Do not add one without saying so.

## Fortran Reference

`emafitpr` in `vendor/peakfqr/src/emafit.f`:

- Inputs: n, ql, qu, tl, tu (all log10), dtype, reg_M/mse, reg_SD/mse, r_G/mse, gbthrsh0,
  pq, nq, eps, wght_opt_n
- Outputs: cmoms(3,3), yp, ci_low, ci_high, var_est, as_G_PRL_o, Wdout, MGBT results
- Call pattern: `vendor/peakfqr/R/fortranWrappers.R::emafit()`; interval construction:
  `vendor/peakfqr/R/readInputs.R::siteQT()` (every year in range gets a row; a year with no
  peak is censored at the *lower* threshold; `dtype` is 1 only for the USGS historic flag)
- Skew MSE encoding: 0 = generalized (no error), <0 = generalized (MSE = −value),
  >0 = weighted, >1e10 = station-only
- MGBT encoding: `gbthrsh0 <= -6` computes MGBT, `> -6` uses it as the threshold
- Weight options: 1=HWN, 2=ERL, 3=INV
- Build note: `build_fortran/build.py` **must** pass `_emafort.pyf`. Without it f2py wraps
  every symbol including QUADPACK's `dqag`, whose callback wrapper does not compile.

## Validation Status

`tests/fortran_parity/` compares the native EMA against committed golden files generated from
peakfq 8.1.0. The `var_mom` port (TODO.md P3) is complete — ADJE's censoring bias adjustment,
`detrat` (the Halloween determinant ratio), the at-site EMA moment iteration on censored
intervals, and the confidence-interval shape (`flowfreq._var_emab.var_emab`, Cohn's
inverse-Gaussian-quadrature method) are all ported and wired into `Bulletin17C`. On the parity
sites this reproduces peakfq 8.1.0 to within measurement noise: weighted skew to 2.4e-6 on Big
Sandy, confidence bounds within 0.06% at every AEP tested, asymmetry ratio within
0.0007-0.0022 of peakfq's own.

**Beyond the golden files: all 24 WY/MT stations** in
`vendor/peakfqr/inst/testdata/wymt_ffa_2022A.psf` now run natively (through `psf_convert`)
and match a live `emafitpr` call on the same inputs. What it took, by PR:

- every perception-threshold period applied per year, as `siteQT` does (#61);
- the exact Pearson III inverse for LP3 quantiles (#62);
- the B17B at-site skew-MSE switch when MGBT finds low outliers (#64, below);
- zero-flow rows, the `gbtmin` censoring bound and the exact MGBT cutoff (#66);
- interval peaks (codes 4/8) and upper perception thresholds (#72);
- confidence bounds interpolated between ±`skewmin` at near-zero skew (#74);
- the EMA fixed point iterating as `p3est_ema` does (#81): at-site skew to 1e-10 and
  quantiles to ~1e-5 % (the near-zero-skew stations 06328100 and 06329350 to 5e-4 % and
  0.012 %).

#90 then added historic interval peaks (code 7 with 4 or 8, `historical_interval_peaks`),
the one `siteQT` construct still refused. No vendored record has one, so it is checked
against live `emafitpr` on a synthetic record
(`tests/fortran_parity/test_historic_interval_peaks_live.py`).

`FrequencyResults.n_low_outliers` is peakfq's `gbnlow` on both engines (#88): every EMA row
`gbtest` recodes below the cutoff. MGBT's own count of flagged peaks is `n_mgbt_outliers`.

The strict xfails, from `grep -rn -A1 "xfail(" tests/` (the marker spans lines, so a
one-line `grep "xfail(strict=True"` finds only docstrings):

- `tests/validation/test_big_sandy.py`: four comparisons against the 2012 PeakfqSA manual
  (quantiles at AEP 0.99 and 0.995, confidence bounds at 0.01 and 0.02). That reference is not
  reproducible by peakfq 8.1.0 (see Test Data).
- `tests/fortran_parity/test_fixed_threshold_live.py` (`requires_fortran`, from #81), two
  places where `emafitpr` itself is the problem:
  - `test_weighted_skew_matches_emafitpr`: a synthetic record with a FIXED 300 cfs threshold
    and at-site skew 0.056. `emafitpr`'s weighted skew is ill-conditioned there, because
    `mP3`'s incomplete gamma rounds at 1e-3 for skews of a few thousandths. Native -0.0037
    sits inside the Fortran's own 1e-5-perturbation band of -0.0073 to -0.0003.
  - `test_big_sandy_6000_weighted_skew_matches`: Big Sandy's systematic record with a
    6000 cfs FIXED threshold (29 of 44 censored). `emafitpr`'s `MN2MVARB` stops after 100
    Newton iterations without converging and returns that iterate unflagged. Its ADJE
    `as_G_mse` is 2.74 against 0.064 at the true root, giving weighted skew -0.281 (peakfq)
    vs -0.166 (native). The at-site fit agrees to 5e-8.

The B17B switch (#64) closed the last golden-file xfail. Cains Coulee's `skew_weighted`, 0.058
skew units off for a long time, now matches to 6e-6 (quantiles to 0.0012%). The cause was `emafit.f:707-711`: when MGBT
computes the low-outlier threshold and finds low outliers, `emafitpr` switches the at-site skew
MSE from ADJE to the plain Bulletin 17B `mseg(n, G)` over the whole record (uncapped), and
leaves `at_site_option` there for the confidence bounds too. The native engine now follows the
same switch (`ExpectedMomentsAlgorithm._at_site_option`); a user-supplied threshold keeps ADJE,
as peakfq's FIXED option does. The "unexplained ~3x `as_G_mse` discrepancy" this was once
attributed to is that switch: a standalone `mseg_all` call runs ADJE
(`tests/fortran_parity/test_fortran_oracles.py::TestCainsCouleeAsGMseDiscrepancy`).

MGBT is the one part verified line-by-line against the Fortran (`GGBCRITP` / `FP_TNC_CDF`),
validated on Orestimba Creek (USGS 11274500, B17C Appendix 10).

**A `requires_fortran` parity test failing on a local Windows/MSYS2 build is not necessarily a
regression.** The golden files were generated with gfortran 13.3.0 on Linux
(`docs/FORTRAN_UPLOAD.md` §6.0). `emafitpr`'s EMA fit is a fixed point with condition number
~1e13 -- one ulp of input moves the converged at-site skew by 3e-3 (see
`tests/fortran_parity/test_fortran_oracles.py`'s module docstring) -- so a newer MSYS2
gfortran (16.2.0 confirmed to reproduce this) computing the same deterministic Fortran call a
few ulps differently is enough to fail a `rel=1e-12` oracle check directly, and to get
amplified into a visibly different `skew_at_site`/`mse_skew` further downstream. Confirmed via
`git stash` that the failure is identical with or without unrelated changes, and that CI's
`ci / Fortran parity` job (Ubuntu, `apt-get install gfortran`, the same toolchain family as the
golden files) stays green. Treat that CI job as the authoritative parity check; a local
Windows mismatch on these two tests alone, with `ci / Fortran parity` green, is toolchain
drift, not a code defect.

## Test Data

- Primary site: Big Sandy River at Bruceton, TN (USGS 03606500). The fixture carries **two**
  sets of expectations: `PEAKFQ_810_*` for parity work, and `EXPECTED_*` from the 2012
  PeakfqSA manual kept as historical record. The 2012 values are **not reproducible** by
  peakfq 8.1.0 — its HWN skew weighting differs by design when censored data are present.
  See `docs/FORTRAN_UPLOAD.md` §6.0b.
- Additional fixtures resolve from `vendor/peakfqr/inst/testdata/`, auto-skipping if absent.
