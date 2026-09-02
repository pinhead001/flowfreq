# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**HydroLib** is a Python library for hydrologic frequency analysis with USGS data retrieval and
Bulletin 17C flood frequency analysis. It has a native Python EMA implementation, validated
against the USGS Fortran reference that is vendored into this repository.

**Reference implementation:** `vendor/peakfqr/` — the USGS `peakfq` R package v8.1.0 (CC0),
containing the authoritative Fortran EMA code. Read `vendor/peakfqr/src/emafit.f` for the
computation and `vendor/peakfqr/R/fortranWrappers.R` for call conventions. It is a verbatim
reference copy: **do not edit anything under `vendor/`**, since a change there silently
invalidates every comparison made against it.

There is no PeakfqSA binary. Earlier work assumed one existed and built `hydrolib/peakfqsa/`
as a subprocess wrapper around it; that premise was wrong, and the subsystem has been removed
— it was mock-tested only and its `requires_peakfqsa` tests collected nothing. What it
contributed to validation now lives in `hydrolib/validation/reference.py`, pointed at
references that exist. The f2py bridge in `hydrolib/peakfqr/` is the only working route to
the reference implementation. `AGENT_BUILD_INSTRUCTIONS_Claude.md` is the original build
brief and explains that history.

## Repository Layout

Everything lives in this repository; there is no external workspace to consult.

```
hydrolib/            the library (flat package, NOT src/hydrolib/)
├── peakfqr/         f2py bridge to the vendored Fortran (built, gitignored)
└── validation/      comparison engine, benchmarks, reports, reference loaders
vendor/peakfqr/      USGS peakfq 8.1.0 reference — Fortran, R, test data (do not edit)
build_fortran/       f2py build script and .pyf signature file
tools/               gen_fortran_golden.py — regenerates parity golden files
tests/fixtures/      shared fixture data (Big Sandy, HU02, Wyoming/Montana, peakfqr respec)
tests/fortran_parity/  native-vs-Fortran parity suite and committed golden files
app/                 Streamlit application (lint + import smoke test in CI)
```

## Build & Development Commands

```bash
pip install -e ".[dev]"

# Full suite. ~90-95 s. It was ~15 s until the confidence-interval shape fix
# (hydrolib._var_emab.var_emab, TODO.md P3): nine regmoms calls per analysis,
# each a full var_mom/mn2mvarb solve. @lru_cache'd like the rest of this
# port's expensive pieces, so repeated fits of the same fixture are cheap,
# but the first fit of any given fixture still pays it. A single cold
# run_analysis() with a regional skew supplied costs a few seconds.
pytest tests/ -v

# What CI actually runs -- the same thing. The marker selection lives in
# pyproject.toml's addopts, so a bare `pytest` is the CI selection.
pytest tests/

# Build the Fortran extension and check the golden files against it. CI runs
# this too, in its own job (gfortran + meson). The parity tests read committed
# golden files and pass without it.
make parity

# Lint app/ and import the Streamlit app (needs -r app/requirements.txt)
make smoke

# Regenerate golden files after changing anything under vendor/peakfqr
python tools/gen_fortran_golden.py

black hydrolib/ tests/ && isort hydrolib/ tests/
```

**Reproduce CI faithfully.** Two traps have both cost a red build:

- `python -m pytest` puts the working directory on `sys.path`; CI runs the `pytest` console
  script, which does not. Use `PYTHONSAFEPATH=1 python -m pytest` to match CI.
- CI tests Python 3.9–3.12. A green run on one interpreter says nothing about the others —
  fixture plumbing in particular diverges (a `@staticmethod` fixture works on 3.11 and breaks
  collection on 3.9).

## Architecture

- **`core.py`** — data structures (`PeakRecord`, `FlowInterval`, `EMAParameters`,
  `FrequencyResults`, `LowFlowResults`), enums, LP3 utilities (`kfactor`,
  `grubbs_beck_critical_value`, `log_pearson3_*`, `compute_ci_lp3`)
- **`bulletin17c.py`** — `FloodFrequencyAnalysis` (ABC), `MethodOfMoments`,
  `ExpectedMomentsAlgorithm`, `Bulletin17C` (facade), and the MGBT three-sweep
- **`usgs.py`** — `USGSgage` NWIS retrieval, `GageAttributes`
- **`lowflow.py`** / **`regime.py`** — low-flow frequency and flow-regime metrics
- **`engine.py`**, **`batch.py`**, **`report.py`**, **`hydrograph.py`**, **`plots.py`**,
  **`freq_plot.py`**

## Conventions

- Run `black` and `isort` before every commit; CI lints `hydrolib/` and `tests/` only
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
intervals, and the confidence-interval shape (`hydrolib._var_emab.var_emab`, Cohn's
inverse-Gaussian-quadrature method) are all ported and wired into `Bulletin17C`. On the parity
sites this reproduces peakfq 8.1.0 to within measurement noise: weighted skew to 2.4e-6 on Big
Sandy, confidence bounds within 0.06% at every AEP tested, asymmetry ratio within
0.0007-0.0022 of peakfq's own.

One `xfail(strict=True)` remains, on Cains Coulee's `skew_weighted` (off by 0.058 skew units).
It is **not** a `var_mom`/`mse_ema` precision limit — that was checked directly and ruled out
(`mse_ema` matches the Fortran oracle to 3e-8 relative at this site's real input). The actual
gap: `emafitpr`'s own internally-computed `as_G_mse` for this one site disagrees with what
calling the same `mseg_all` Fortran routine gives standalone on identical inputs, by roughly
3x, for a reason not pinned down despite investigation (see
`tests/fortran_parity/test_fortran_oracles.py::TestCainsCouleeAsGMseDiscrepancy` and TODO.md
P3 for the full account). hydrolib's own computation may be more correct here than the golden
reference, not less — there is no known defect to fix.

MGBT is the one part verified line-by-line against the Fortran (`GGBCRITP` / `FP_TNC_CDF`),
validated on Orestimba Creek (USGS 11274500, B17C Appendix 10).

## Test Data

- Primary site: Big Sandy River at Bruceton, TN (USGS 03606500). The fixture carries **two**
  sets of expectations: `PEAKFQ_810_*` for parity work, and `EXPECTED_*` from the 2012
  PeakfqSA manual kept as historical record. The 2012 values are **not reproducible** by
  peakfq 8.1.0 — its HWN skew weighting differs by design when censored data are present.
  See `docs/FORTRAN_UPLOAD.md` §6.0b.
- Additional fixtures resolve from `vendor/peakfqr/inst/testdata/`, auto-skipping if absent.
