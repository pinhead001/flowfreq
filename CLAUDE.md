# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**HydroLib** is a Python library for hydrologic frequency analysis with USGS data retrieval and Bulletin 17C flood frequency analysis. This branch (`dev` in `hybrid-17c-cld`) implements a hybrid Bulletin 17C approach that wraps PeakfqSA (USGS Fortran reference) alongside the existing native EMA implementation.

**Reference implementation:** `C:\a\hal\_shared\peakfqr` — USGS R package (peakfq v8.1.0) containing the most recent working copy of the Fortran EMA code. This is the authoritative reference (not a standalone PeakfqSA binary). Read `src/` for computation logic, `R/fortranWrappers.R` for call conventions.

## Workspace Layout

```
C:\a\hal\
├── main/                    # HydroLib main branch (read-only reference)
├── hybrid-17c-cld/          # This worktree — dev branch for hybrid 17C
├── _shared/peakfqr/         # USGS R package reference (Fortran + R)
└── _shared/testdata/        # Shared test data
```

## Build & Development Commands

```bash
# Install dependencies
pip install -e ".[dev]"

# Install dev tools
pip install pytest pytest-cov pytest-mock black isort mypy ruff pandas numpy scipy click rich

# Run all tests
pytest tests/ -v

# Run a single test file
pytest tests/test_bulletin17c.py -v

# Run tests excluding PeakfqSA-dependent tests
pytest tests/ -v -m "not requires_peakfqsa"

# Coverage
pytest tests/ --cov=hydrolib --cov-report=term-missing

# Formatting
black src/ tests/ hydrolib/
isort src/ tests/ hydrolib/

# Type checking
mypy src/hydrolib/ --ignore-missing-imports --strict
```

## Existing Architecture

### Package: `hydrolib/`
- **`core.py`** — Data structures (`PeakRecord`, `FlowInterval`, `EMAParameters`, `FrequencyResults`), enums (`SkewMethod`, `AnalysisMethod`), LP3 utility functions (`kfactor`, `grubbs_beck_critical_value`, `log_pearson3_*`, `compute_ci_lp3`)
- **`bulletin17c.py`** — Main analysis: `FloodFrequencyAnalysis` (ABC), `MethodOfMoments`, `ExpectedMomentsAlgorithm`, `Bulletin17C` (facade/factory)
- **`usgs.py`** — `USGSgage` for NWIS data retrieval, `GageAttributes` singleton
- **`engine.py`** — `B17CEngine` simplified LP3 fitting interface
- **`batch.py`** — Multi-site parallel analysis
- **`report.py`** — `HydroReport` technical report generation
- **`hydrograph.py`** / **`plots.py`** — Visualization

### New Hybrid Modules (being built)
- **`src/hydrolib/peakfqsa/`** — PeakfqSA subprocess wrapper (config, wrapper, io_converters, parsers, validators)
- **`src/hydrolib/validation/`** — Comparison engine, benchmarks, reports
- **`tests/peakfqsa/`**, **`tests/validation/`**, **`tests/integration/`** — Tests

## Key Conventions

- **Commit frequently** — minimum one commit per step and substep (Step0a, Step1b, etc.). Reference step number in commit message: `feat(peakfqsa): add config module [Step 5]`
- **Lint and document as you go** — run `black` and `isort` before each commit; write NumPy-format docstrings on all public classes/functions
- Type hints on all function signatures
- No bare `except:` — always catch specific exceptions
- No `print()` in library code — use `logging.getLogger(__name__)`
- Tests for every public function (happy path + error case minimum)
- Bulletin 17C compliance is critical — use peakfqr `src/` as the authoritative specification
- All data in log10 space when interfacing with Fortran conventions (peakfqr logs base-10)

## Fortran/R Reference Quick Reference

Key Fortran routine: `emafitpr` in `_shared/peakfqr/src/emafit.f`
- Inputs: n, ql, qu, tl, tu (all log10), dtype, reg_M/mse, reg_SD/mse, r_G/mse, gbthrsh0, pq, nq, eps, wght_opt_n
- Outputs: cmoms(3,3), yp (quantiles), ci_low, ci_high, var_est, MGBT results (gbval, gbns, gbnzero, gbnlow, gbp, gbqs)
- R wrapper: `R/fortranWrappers.R` — `emafit()` function shows exact call pattern and output field extraction
- Skew MSE encoding: 0 = generalized (no error), <0 = generalized (MSE = -value), >0 = weighted, >1e10 = station-only
- MGBT encoding: gbthrsh0 <= -6 = compute MGBT, > -6 = use as threshold, ~-5.9 = no test
- Weight options: 1=HWN (Halloween), 2=ERL, 3=INV (PeakFQ 7.4.1)

## Test Data

- Primary validation: Big Sandy River at Bruceton, TN (USGS 03606500) — from PeakfqSA manual
- Additional fixtures from `_shared/peakfqr/inst/testdata/` and `_shared/peakfqr/tests/testthat/`
