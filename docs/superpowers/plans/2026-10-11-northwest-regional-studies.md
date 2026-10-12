# Northwest Regional Studies (WA, MT, ID, OR) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fill the remaining Wave 1 gaps (WA region 4 covariance, the MT regional skew, the
Idaho Snake River Plain skew, and Oregon's unresolved skew basins) with studies that follow
published USGS procedures, since no newer USGS report covers any of them (searched
2026-10-11).

**Architecture:** First, two pieces of infrastructure: a `provisional` status for
non-USGS skews, which callers must opt into, and a WREG-equivalent GLS regression
module. Each state study is then a `tools/build_*` script that writes committed inputs and
results under `data/`, plus a `docs/*.md` report and tests. A study whose
method is validated against a published USGS result can feed a packaged table only as
`provisional`. `verified` stays reserved for transcribed USGS values.

**Tech Stack:** Python 3.11–3.14, numpy, scipy, pandas, requests (no new dependencies);
`flowfreq.skew_study`, `flowfreq.regional_skew`, `flowfreq.regulation`,
`flowfreq.streamstats`, the Water Data OGC API.

**Spec:** this plan stands in for one. The requirements come from:
- `docs/MASTER_ROADMAP.md` §1.3 and §3.3;
- TODO.md "Carried into 0.11.0";
- PR #65, the WA region 4 covariance evidence;
- `docs/MONTANA_REGIONAL_SKEW_PROVISIONAL.md` §Limitations, the MT v1 gaps.

## Research findings (2026-10-11)

- No newer USGS peak-flow or regional-skew report covers WA, MT, ID or OR.
  - SIR 2016-5118 (WA), SIR 2016-5083 (ID, PNW skew) and SIR 2005-5116 / OWRD 2006 (OR) remain current.
  - SIR 2025-5019 (MT/ND/SD/WY) says no B-WLS/B-GLS skew study covers MT yet.
  - Oregon StreamStats has no eastern Oregon peak-flow equations.
- **MT:** the SIR 2025-5019 data release (doi:10.5066/P1WHRK8H) has a PeakFQ `.psf` and
  `ffa_2022examples.xlsx` for 63 example stations. Its `.psf` confirms that USGS applies
  B17B Plate I map skews per gage, with `SkewSE 0.55`. That workbook is a PeakFQ station-skew
  cross-check for Task 6.
- **WA:** SIR 2016-5118 publishes three spreadsheets that make a Region 4 refit possible:
  - `table2.xlsx`: every gage's EMA inputs (perception thresholds, historic periods, MGBT);
  - `table5.xlsx`: basin characteristics, regression region and the "used in regression" flag;
  - `table8.xlsx`: station, regression and weighted quantiles.
- **NM (Wave 2, out of scope here):** Tillery and others (2026) supersedes Waltemeyer (2008).
  - Nine regions, data through WY2020, regional skew -0.145 (SE 0.454).
  - Data release doi:10.5066/P9QWEXCY. Logged in TODO.md by Task 1.

## Global Constraints

- Python 3.11–3.14; no syntax beyond 3.11. No new runtime or dev dependency: read `.xlsx`
  with `zipfile` + `xml.etree` (the approach `tools/build_regulation_screen.py` takes).
- `black>=24,<25` formatting, isort, and mypy are enforced on every new `flowfreq/` module.
- No `print()` in `flowfreq/`; use `logging.getLogger(__name__)`. `tools/` scripts may print.
- NumPy docstrings and type hints on public API. Tests: happy path plus at least one error case.
- Known-failing comparisons are `xfail(strict=True)`. Never widen a tolerance to pass.
- Never edit `vendor/`.
- `verified` = a transcribed USGS value only. Nothing computed here is `verified`.
  - A flowfreq-developed value is `provisional` and is returned only when the caller passes
    `allow_provisional=True`.
  - It always carries an advisory.
- Live-service tests are `@pytest.mark.requires_network`. Offline tests use committed data.
- One PR per task (or per tightly coupled pair), merged in order once CI is green
  (memory: merge-workflow).

## Review Focus

1. **A `provisional` skew leaking into a default call.** `regional_skew_at_huc(huc)`,
   `regional_skew_for(state)` and `run_ffa` without opt-in must raise exactly as today.
   Pinned in Task 1.
2. **GLS on a singular or near-singular design**, such as a region with fewer sites than
   parameters, or duplicate sites. It must raise `ValueError`, not return NaNs. Pinned in
   Task 2.
3. **A gage with no concurrent record with any other gage.** The cross-correlation term must
   be 0, not NaN. Pinned in Task 2.
4. **Washington's Region 4 refit landing between the two bases.** The test asserts which
   basis is closer by a margin. If neither is within tolerance, the result is reported and WA
   stays `partial`; it is never forced. Pinned in Task 4.
5. **Study gages in a state the row does not cover.** A provisional MT row must not apply
   in WY or ID HUCs even though WY/ID gages were used. Pinned in Task 7.

---

### Task 1: `provisional` status in the regional skew tables

**Files:**
- Modify: `flowfreq/regional_skew.py` (`STATUSES`, `load_table`, `regional_skew_for`,
  `regional_skew_at_huc`, `regional_skew_at`, `regional_skew_for_site`, the `RegionalSkew` dataclass)
- Modify: `TODO.md` (log the NM Tillery 2026 finding and this plan)
- Test: `tests/test_regional_skew.py`, `tests/test_regional_skew_at.py`

**Interfaces:**
- Produces:
  - `STATUSES = frozenset({"pending", "verified", "provisional"})`.
  - Keyword `allow_provisional: bool = False` on `regional_skew_for`, `regional_skew_at_huc`,
    `regional_skew_at` and `regional_skew_for_site`.
  - New field `RegionalSkew.provisional: bool = False`.
  - A provisional result always has a non-empty `advisory` starting with
    `"PROVISIONAL (not USGS)"`.
  - `load_table` requires provisional rows to carry the same value columns as verified rows.

- [ ] **Step 1: Write the failing tests**

```python
def test_provisional_row_needs_opt_in(tmp_path):
    p = _table(tmp_path, 'XX,Test,-0.1,0.3,25,"Cite, 2026",T1,provisional,n\n')
    with pytest.raises(RegionalSkewUnavailable):
        regional_skew_for("XX", path=p)
    rs = regional_skew_for("XX", path=p, allow_provisional=True)
    assert rs.provisional and rs.skew == -0.1
    assert rs.advisory.startswith("PROVISIONAL (not USGS)")

def test_verified_preferred_over_provisional(tmp_path):
    body = ('XX,A,-0.1,0.3,25,"C",T,verified,n\n'
            'XX,B,0.2,0.3,25,"C",T,provisional,n\n')
    rs = regional_skew_for("XX", path=_table(tmp_path, body), allow_provisional=True)
    assert rs.skew_region == "A" and not rs.provisional

def test_provisional_row_missing_value_rejected(tmp_path):
    with pytest.raises(ValueError, match="provisional row"):
        load_table(_table(tmp_path, 'XX,T,,0.3,25,"C",T,provisional,n\n'))
```

Plus a HUC-table test in `test_regional_skew_at.py`. Use a tmp skew table with a provisional
region and a tmp HUC table with one `member` row for it. Assert that `regional_skew_at_huc`
raises without opt-in and returns `provisional=True` with it.

- [ ] **Step 2: Run them; expect FAIL** (`pytest tests/test_regional_skew.py tests/test_regional_skew_at.py -k provisional -v`)
- [ ] **Step 3: Implement.** In `regional_skew_at_huc`, the `verified` filter becomes
  "verified, else provisional when `allow_provisional`". `regional_skew_for` prefers verified
  over provisional, so when both exist for a state the verified one wins.
- [ ] **Step 4: Run the full skew tests; expect PASS, and every existing test unchanged.**
- [ ] **Step 5: Commit** with the plan doc and the TODO.md entries, then open a PR.

### Task 2: GLS regression for flood quantiles (`flowfreq/regression/gls.py`)

**Files:**
- Create: `flowfreq/regression/gls.py`
- Test: `tests/test_regression_gls.py`

**Interfaces:**
- Consumes: `flowfreq.skew_study.distance_miles`.
- Produces:
  - `lp3_quantile_sampling_variance(sd, skew, n, aep, skew_mse) -> ndarray`: per-site
    variance of the log quantile estimator. WREG (Eng and others 2009, TM 4-A8) form,
    including the regional-skew MSE term.
  - `WregCorrelation(alpha: float, theta: float)` with `.rho(d_miles) -> ndarray`:
    rho = theta ** (d / (alpha * d + 1)).
  - `fit_correlation(...)`: fits that model to station-pair correlations by least squares
    on (d, r) pairs with >= `min_overlap` concurrent years.
  - `sampling_covariance(sd, skew, n, aep, skew_mse, lat, lon, concurrent, corr) -> ndarray`:
    the Σ matrix. The off-diagonal uses rho_ij * sqrt(Var_i Var_j) * c_ij / sqrt(n_i n_j),
    with c_ij the concurrent years. Zero where c_ij = 0.
  - `@dataclass GLSResult`:
    - fields `beta`, `cov_beta` (that is (X' Lambda^-1 X)^-1), `model_error_variance`, `avp`,
      `sep_percent`, `pseudo_r2`, `leverage`, `influence`, `residuals`, `names`;
    - method `predict(x) -> (yhat, var_pred)`.
  - `fit_gls(y, X, sigma, names=None) -> GLSResult`.

- [ ] **Step 1: Write the failing tests.**
  - **Identity:** with Σ = 0 and a well-posed design, `fit_gls` equals OLS. β matches
    `numpy.linalg.lstsq` to 1e-10, and MEV equals SSE/(n-p) to 1e-10.
  - **Known MEV:** simulate y = Xβ + δ + ε with γ² = 0.04, independent ε with known
    diagonal Σ, and n = 2000 (seeded). The fitted MEV is within 0.004 of 0.04.
  - **Singular X:** two identical columns raise `ValueError, match="rank"`.
  - **Fewer sites than parameters** raises `ValueError`.
  - **Zero concurrent years:** that Σ off-diagonal is exactly 0.0.
  - **Variance formula:** `lp3_quantile_sampling_variance(0.3, 0.0, 50, 0.01, 0.0)` equals
    0.3²/50 · (1 + K²/2), with K = 2.3263479 (normal z at AEP 0.01), to 1e-8.
- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Implement.** MEV is the method-of-moments root of
  `(y - Xβ̂)' Λ⁻¹ (y - Xβ̂) = n - p`, with `Λ = γ² I + Σ` (WREG's GLS). Use
  `scipy.optimize.brentq` on γ² ∈ [0, var(y)·10], and return 0 when the equation has no root
  above 0. Pseudo-R² is `1 - γ²(k)/γ²(0)`, where γ²(0) comes from a constant-only fit with
  the same Σ.
  - AVP = γ² + mean_i x_i' cov_beta x_i.
  - SEP% = 100·sqrt(exp(ln(10)²·AVP) - 1).
- [ ] **Step 4: Run; expect PASS. `mypy flowfreq/regression/gls.py` is clean.**
- [ ] **Step 5: Commit and open a PR.**

### Task 3: Washington SIR 2016-5118 regression dataset

**Files:**
- Create: `tools/build_wa_gls_dataset.py`, `tools/_xlsx.py` (a minimal sheet reader)
- Create: `data/regression_study/wa_2016_5118_sites.csv`
- Test: `tests/test_wa_gls_dataset.py`

**Interfaces:**
- Produces `wa_2016_5118_sites.csv`, one row per Table 5 gage, with columns:
  - `site_no, region, used, lat_cent, lon_cent, DRNAREA, PRECIP, CANOPY_PCT`;
  - `n_sys, hist_begin, hist_end` (from Table 2);
  - `q_<aep>` for the eight AEPs (Table 8 line 1, station EMA);
  - `mu, sd, skew_w`: log10 LP3 moments, recovered by least squares from the eight station
    quantiles, since Table 8 line 1 is an LP3 with weighted skew.
- `tools/_xlsx.py: read_sheet(path: Path, sheet_index: int = 0) -> list[list[str]]`.

- [ ] **Step 1: Write the failing tests.**
  - The CSV has one row per Table 5 data row (count it from the sheet, then pin it), unique
    `site_no`, and `region` in {1, 2, 3, 4}.
  - Region 4 `used == True` has 139 rows (the calibration count in PR #65).
  - For three spot gages (12024400, 12020800, 12398000), `q_0.01` equals the value in PR
    #65's tests.
  - The recovered LP3 reproduces all eight station quantiles within 0.5 % at every used gage.
- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Implement the builder.** It downloads the three `.xlsx` from
  `https://pubs.usgs.gov/sir/2016/5118/sir20165118_table{2,5,8}.xlsx`, caches them in
  `--cache`, and writes the CSV. Its module docstring records the provenance.
- [ ] **Step 4: Run it, then the tests; expect PASS.**
- [ ] **Step 5: Commit the data and the tool. Open a PR (it may be combined with Task 4).**

### Task 4: Washington GLS refit and the Region 4 covariance decision

**Files:**
- Create: `tools/wa_region4_refit.py`, `docs/WA_REGION4_REFIT.md`
- Create: `data/regression_study/wa_2016_5118_refit.json`
- Modify: `flowfreq/data/regression/WA.json`, but only when the outcome rule below says so
- Test: `tests/test_wa_region4_refit.py`

**Interfaces:**
- Consumes: Task 2 `fit_gls`, `sampling_covariance` and `fit_correlation`; Task 3's CSV.
- Produces a refit JSON keyed `region -> aep -> {beta, cov_beta, mev, avp, basis}`. Region 4
  has two entries, `basis: "P"` and `basis: "P/10"`.

- [ ] **Step 1: Write the failing tests.**
  - **Validation:** Regions 1-3 refit at AEP 0.01 reproduce Table 6.
    - log10 intercept within 0.03;
    - every slope within 10 % of the stored `WA.json` coefficient;
    - MEV within 20 % of Table 7.
  - **Region 4:** the test that records the outcome compares the refit (X'Λ⁻¹X)⁻¹ against
    Table 7's printed matrix, by the ratio of diagonals, in each basis.
    - It asserts that the basis whose median diagonal ratio is closest to 1 is the one
      written in the JSON as `"selected_basis"`.
    - It asserts that that ratio is within [0.67, 1.5].
  - If the validation fails, mark it `xfail(strict=True)` with the measured numbers and stop
    before changing `WA.json`.
- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Implement.** Station sampling variances use `sd` and `skew_w` from Task 3;
  `skew_mse` is 0.18 (PNW) and `n` is the pseudo record length (`n_sys` plus historic
  years). The correlation model is fitted from WA gage pairs with ≥ 20 concurrent years,
  using the Water Data API peaks of the used gages (cached).
- [ ] **Step 4: Apply the outcome rule.** Washington's status is decided here, so the rule
  is strict:
  - Store Table 7's printed Region 4 matrix with `Variable.scale` in the basis the refit
    selects, and set WA to `verified`, only if all three hold:
    - the validation passed;
    - the selected basis is within [0.67, 1.5];
    - the selected basis also reproduces Table 6's Sp within 20 %.
  - Otherwise WA stays `partial`. `docs/WA_REGION4_REFIT.md` and the `WA.json` notes record
    the numbers either way.
- [ ] **Step 5: Run the full suite, regenerate `docs/REGRESSION_COVERAGE.md`
  (`python tools/gen_regression_coverage.py`), commit, and open a PR.**

### Task 5: Montana skew study v2: inputs to USGS standard

**Files:**
- Modify: `tools/build_montana_skew_study.py` (add `--version v2`; v1 output stays reproducible)
- Create: `data/skew_study/montana_v2_inputs.csv`
- Test: `tests/test_montana_skew_study_v2.py`

**Interfaces:**
- Produces `montana_v2_inputs.csv`. It carries the v1 columns plus:
  - `state`;
  - `lat_cent, lon_cent, ELEV, PRECIP, DRNAREA_ss` (StreamStats basin characteristics);
  - `div_dist_mi`: signed distance from the basin centroid to the Continental Divide,
    positive east. The divide is the WBD HUC2 boundary between region 17 and regions 10/09.
  - `reg_class` (from `flowfreq.regulation.classify_site`);
  - `excluded_reason`.

- [ ] **Step 1: Write the failing tests** against the committed CSV.
  - Gages from WY and ID in HUC4 1002, 1007, 1008 and 1701 are present.
  - Every used gage has a basin centroid and ELEV.
  - Every `regulated`/`urban` class gage is excluded.
  - `div_dist_mi` is > 0 for 06090800 (Fort Benton) and < 0 for 12340500 (Missoula).
  - No used pair is redundant (SD ≤ 0.5 and DAR ≤ 5), now computed with centroids.
- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Implement.**
  - Use `streamstats` batch delineation with its cache, at concurrency ≤ 3.
  - Neighbor gages come from the same OGC queries with `state_code` 16 (ID) and 56 (WY),
    restricted to those HUC4s.
  - The regulation screen is `classify_site` with the packaged GAGES-II screen plus peak codes.
- [ ] **Step 4: Run the builder (long; run it in the background and resume from the cache),
  then the tests; expect PASS.**
- [ ] **Step 5: Commit and open a PR.**

### Task 6: Montana v2 station skews checked against USGS PeakFQ

**Files:**
- Create: `tools/check_montana_station_skews.py`
- Create: `data/skew_study/sir2025_5019_station_skews.csv` (MT rows of `ffa_2022examples.xlsx`)
- Test: `tests/test_montana_station_skew_check.py`

**Interfaces:**
- Consumes: Task 5's inputs and `skew_study.station_skew`.
- Produces the CSV with `site_no, usgs_station_skew, usgs_n, flowfreq_station_skew`.

- [ ] **Step 1: Write the failing test.** For the MT example stations also in the v2 data,
  the flowfreq station skew is within 0.05 of USGS's, by median absolute difference. Any
  station off by more than 0.2 is listed by name in an `xfail(strict=True)` parametrization,
  each with its reason.
- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Implement.** Parse the workbook with `tools/_xlsx.py` (Task 3). Station
  skews use the `.psf` specs (thresholds, historic periods) through `psf_convert`, so both
  sides see the same inputs.
- [ ] **Step 4: Run; expect PASS (the documented xfails apart).**
- [ ] **Step 5: Commit and open a PR (it may be combined with Task 7).**

### Task 7: Montana v2 B-WLS/B-GLS fit and provisional table rows

**Files:**
- Modify: `tools/build_montana_skew_study.py` (the model set)
- Create: `data/skew_study/montana_v2_results.json`, `docs/MONTANA_REGIONAL_SKEW_V2.md`
- Modify: `flowfreq/data/regional_skew.csv`, `tools/build_regional_skew_hucs.py`,
  `flowfreq/data/regional_skew_hucs.csv`, `docs/MONTANA_REGIONAL_SKEW_PROVISIONAL.md`
  (a pointer to v2)
- Test: `tests/test_montana_skew_v2.py`

**Interfaces:**
- Consumes: Task 1's provisional status and Task 5's inputs.
- Produces the candidate models: CONSTANT; + `div_dist_mi`; + ELEV; + PRECIP; + a Front
  indicator (0 < `div_dist_mi` ≤ 60 mi); + an east/west indicator.

- [ ] **Step 1: Write the failing tests.**
  - The results JSON has every candidate.
  - The selected model is the one with the smallest AVPnew among models whose added
    parameters have |β/SD| ≥ 2.
  - With `allow_provisional=True`, `regional_skew_at_huc("100301021602", "MT")` (Fort
    Benton) returns `provisional=True`.
  - Without opt-in it still raises.
  - `regional_skew_at_huc("170401030103", "WY", allow_provisional=True)` raises "outside"
    (Review Focus 5).
- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Implement the fit and the table rows.**
  - **CONSTANT or indicator models:** write one provisional row per skew region, plus HUC
    member rows (HUC4 1002-1011 and 0904 within MT; `valid_states` MT).
  - **Continuous covariate models** cannot be represented in the table. Then the CONSTANT
    value is stored as the provisional row, and the doc gives the covariate equation.
  - West of the Divide (1701), the verified PNW row keeps precedence.
- [ ] **Step 4: Run the full suite; expect PASS.**
- [ ] **Step 5: Commit and open a PR.**

### Task 8: Snake River Plain skew feasibility (ID)

**Files:**
- Create: `tools/snake_river_plain_skew.py`, `docs/SNAKE_RIVER_PLAIN_SKEW.md`,
  `data/skew_study/srp_inputs.csv`
- Test: `tests/test_snake_river_plain_skew.py`

- [ ] **Step 1: Write the failing test.** The inputs CSV lists every peak-flow gage whose
  location is in SIR 2016-5083 region 0, with screening outcomes. The doc's decision line
  ("feasible" or "not feasible") agrees with the rule below.
- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Implement.** Apply Task 5's screening (regulation class, ≥ 25 pre-regulation
  peaks, redundancy).
  - **Feasible** only with ≥ 20 screened gages. Then fit CONSTANT B-GLS, compare it with the
    PNW value by t = (β - (-0.07)) / SD(β), and write a provisional row (with the
    `excluded` HUC12s reassigned as members of the new region).
  - **Otherwise** record "not feasible" with the counts; no table change.
- [ ] **Step 4: Run; expect PASS.**
- [ ] **Step 5: Commit and open a PR.**

### Task 9: Oregon closed basins (1712) and Klamath (1801): the PNW applicability test

**Files:**
- Create: `tools/pnw_edge_skew_test.py`, `docs/OREGON_UNRESOLVED_SKEW.md`
- Modify (conditional): `flowfreq/data/regional_skew.csv`, `tools/build_regional_skew_hucs.py`
- Test: `tests/test_pnw_edge_skew.py`

- [ ] **Step 1: Write the failing test.** For OR gages in HUC4 1712 and 1801 (screened as in
  Task 5), the GLS mean residual against G = -0.07 uses
  σδ² = 0.17 plus each gage's sampling variance and Task 2's correlation model.
  - If |t| < 2, the doc says "consistent". A provisional region named
    `"Pacific Northwest (OR closed basins, Klamath)"` then carries -0.07/0.18, and the HUC
    rows switch from `unresolved` to `member` for it, valid in OR only.
  - Otherwise the doc says "inconsistent", and no table changes.
- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run; expect PASS.**
- [ ] **Step 5: Commit and open a PR.**

### Task 10: Close-out

- [ ] Update TODO.md (the carried items), roadmap §1.3/§3.3, CHANGELOG [Unreleased] and the
  `flowfreq-state` memory with each study's outcome.
- [ ] Write a follow-on plan for the Oregon regression update. It is out of scope here: it
  needs Task 2's GLS validated, and StreamStats computes no eastern-Oregon variables.
