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

# Full suite. ~95-100 s. It was ~15 s until the confidence-interval shape fix
# (flowfreq._var_emab.var_emab, TODO.md P3): nine regmoms calls per analysis,
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
- **`usgs.py`** — `USGSgage` NWIS retrieval, `GageAttributes`
- **`lowflow.py`** / **`regime.py`** — low-flow frequency and flow-regime metrics
- **`engine.py`**, **`batch.py`**, **`report.py`**, **`hydrograph.py`**, **`plots.py`**,
  **`freq_plot.py`**

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

No parity `xfail` remains; the only strict xfails left are the four 2012 PeakfqSA manual
comparisons (see Test Data). Cains Coulee's `skew_weighted`, 0.058 skew units off for a long
time, now matches to 6e-6 (quantiles to 0.0012%). The cause was `emafit.f:707-711`: when MGBT
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
