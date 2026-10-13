# Oregon Peak-Flow Regression Update: Follow-on Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Bulletin 17C-based regional regression for Oregon peak flows, statewide, by the
GLS procedure USGS uses (WREG). It replaces neither of the current equations until reviewed:
- Cooper (2005, SIR 2005-5116), western Oregon, B17B era;
- Cooper (2006, OWRD OFR SW 06-001), eastern Oregon, non-USGS.

**Architecture:**
- Station quantiles from flowfreq EMA/MGBT, weighted with the verified PNW skew (the
  provisional closed-basin extension in 1712).
- Basin characteristics computed by a `tools/` script from national rasters, *not* StreamStats
  (see Finding).
- `flowfreq.regression.gls` fits per region and AEP. It was validated on Washington SIR
  2016-5118 in `docs/WA_REGION4_REFIT.md`.
- Results go in a new `provisional` equation status, which needs a schema decision (Task 1).

**Tech Stack:** numpy, scipy, pandas, `flowfreq.regression.gls`, `flowfreq.skew_study`. Raster
work stays in `tools/` with optional dependencies, never in `flowfreq/` (NFR-6, no hard GIS
dependency).

**Spec:**
- Northwest studies plan task 10 (`docs/superpowers/plans/2026-10-11-northwest-regional-studies.md`);
- roadmap §3.3;
- Finding below.

## Finding (2026-10-12) that shapes this plan

Oregon StreamStats returned only 9-13 characteristics for test delineations at 14306500 (Alsea,
western) and 14046500 (John Day, eastern):
- **Returned:** DRNAREA, DRNDENSITY, LFPLENGTH, CSL1085LFP, STRMTOT, ORREG2 and road densities;
  PRECIP, RELIEF, ELEVMAX and MINBELEV in the east only.
- **Unavailable:** everything else, including I24H2Y and JANMAXT2K, which Cooper (2005) needs.
- **Server error:** "Cannot find a 'awcrisleycpg' raster for watershed 1710".

That may be a transient server fault, so it must be re-checked before Task 3 starts. If it
persists:
- flowfreq's online Oregon workflow (`delineate_and_get_characteristics` → `OR.json`) cannot
  evaluate the western equations;
- a regression update cannot rely on StreamStats for predictors.

Report it to the StreamStats team either way.

## Global Constraints

- Python 3.11–3.14; black 24, isort and mypy on `flowfreq/`. No new runtime dependency.
  Raster libraries (rasterio or similar) are allowed only as an optional extra used by `tools/`.
- Nothing developed here is `verified`. Equations enter the library as `provisional` and are
  returned only on opt-in, mirroring `regional_skew` (PR #130).
- Known-failing comparisons are strict xfails, never widened tolerances.

## Review Focus

1. **A provisional equation set changing the default `load_state("OR")`.** It must not. Pinned
   in Task 1.
2. **Gages on regulated mainstems** (Willamette, Owyhee, Klamath) slipping through.
   GAGES-II plus peak-code screening; pinned in Task 2.
3. **Basin characteristics** computed from a different DEM or precipitation normal than the
   equations will be evaluated with. One source, recorded per characteristic; pinned in Task 3.
4. **Closed basins (1712)** using the provisional skew. Every station quantile records its
   skew source; pinned in Task 2.
5. **Region boundaries** that make an ungaged basin straddle two regions. Use area weighting,
   as `oregon.blend` does for 2A/2B; pinned in Task 5.

---

### Task 1: `provisional` status for regression equation sets (schema decision) — done

Approved by the maintainer 2026-10-12, implemented in `tests/test_regression_provisional.py`'s PR.

- [x] `provisional` added to `flowfreq/regression/library.py` `STATUSES`.
- [x] `load_state(code)` keeps returning the published set; a published file claiming
  `provisional` is refused.
- [x] `load_state(code, allow_provisional=True)` returns
  `flowfreq/data/regression/provisional/<CODE>.json` when it exists, with `advisory` set and
  a logged warning. A subdirectory rather than the proposed `OR_provisional.json`: the
  top-level `*.json` glob (`available_states`, the coverage generator, two test loops) treats
  every stem as a state code.
- [x] `provisional_states()` lists them; `docs/REGRESSION_COVERAGE.md` marks a state that has
  one ("+ provisional set").

### Task 2: Oregon station dataset

- `tools/build_oregon_regression_dataset.py`. Every OR peak-flow gage, plus border gages within
  25 mi (WA, ID, NV, CA), from bulk OGC queries with `USGS_API_KEY`.
- Screens as the Montana v2 study.
- EMA with the verified PNW skew (provisional for 1712; at-site only where neither applies,
  flagged).
- Output: station quantiles at the 8 AEPs, station SD and skew, pseudo record length, and the
  skew source.
- Validation, in tests:
  - flowfreq EMA reproduces SIR 2016-5083 Table 1 quantiles at the Idaho/Oregon border gages
    both studies used;
  - the median ratio is within 1 %, as on the Washington controls.

### Task 3: Basin characteristics

- Re-check Oregon StreamStats first.
- If it still fails, compute from national sources in `tools/`, each recorded with source and
  date:
  - DRNAREA from StreamStats (which works);
  - PRISM 1991-2020 mean annual precipitation and 2-year 24-hour intensity from NOAA Atlas 14
    rasters, over the StreamStats polygon;
  - mean elevation and slope from 3DEP;
  - forest cover from NLCD 2021.
- Candidate set: DRNAREA, PRECIP, I24H2Y, ELEV, slope, forest, JANMAXT (PRISM).

### Task 4: Regions

- Start from the current regions (western 1, 2A, 2B; OWRD's six eastern).
- Test with residual clustering from a statewide GLS, and merge where coefficients don't differ
  significantly.
- Document every boundary change.

### Task 5: GLS fits, diagnostics and comparison

- Per region and AEP: variable selection by t ≥ 2 and VIF < 5, then SEP, pseudo-R², leverage
  and influence.
- LOOCV against the current equations at the same gages.
- Store as `provisional` (Task 1) with the full covariance.
- Doc: `docs/OREGON_REGRESSION_UPDATE.md`.
