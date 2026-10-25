# Nationwide Rollout Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Every one of the 56 jurisdictions gets an offline regression estimate where NSS serves one, a single `estimate_at_site` entry point ties delineation, region selection, the library and gage weighting together, and the app exposes both the gaged and the ungaged path, by 2026-11-15.

**Architecture:** Workstream A adds an `nss_snapshot` opt-in status to the existing library loader and a generator from `data/nss_snapshots/` files, after three parser extensions. Workstream B adds `flowfreq/ungaged.py` on top of the existing StreamStats, regression and gage-weighting modules. C is data (regional skew rows). D is the app. E is the evidence run, docs and releases. A, B, C and D1 are independent and run in parallel; D2 and E follow the 0.16.0 release.

**Tech Stack:** Python 3.11-3.14, pandas, requests, Streamlit (app), pytest, `make check` (black, isort, mypy, pytest).

**Spec:** `docs/superpowers/specs/2026-10-25-nationwide-rollout-design.md`

## Agents and models

| Workstream | Tasks | Implementer model | Reviewer model | Why |
|---|---|---|---|---|
| A parser + tier | A1-A5 | opus | opus | parser edge cases and a loader contract other code depends on |
| A data | A6-A7 | sonnet | opus | generated files; the review checks counts, not prose |
| B ungaged | B1-B4 | fable | opus | the one design-heavy module; interfaces everything else uses |
| C skew rows | C1 (4 agents, one per wave group) | opus | opus (C2 is itself the review, a different agent) | citation accuracy from PDFs |
| D app | D1 | sonnet | opus | cherry-pick and drift fixes |
| D app | D2 | sonnet | opus | Streamlit page with bare-mode tests |
| E smoke, docs, release | E1-E3 | sonnet | opus | scripts and prose |
| Final whole-branch review per workstream | | | fable | the last gate before merge |

Each workstream is its own branch and PR, in its own worktree (`superpowers:using-git-worktrees`; set `PYTHONPATH=$PWD` there, see memory). Merge order: A, B, C in any order once green; then release 0.16.0 (E3); then D2; then E1-E2.

## Schedule (working days from 2026-10-27)

| Days | Work |
|---|---|
| 1-3 | A1-A5, B1-B3, C1 (4 agents), D1, all in parallel |
| 4-5 | A6-A7, B4, C2, reviews, merges, coverage doc regenerated |
| 6 | Release 0.16.0 (E3, first half); app pin bump PR to the user |
| 7-9 | D2; E1 smoke run and the fixes it finds |
| 10-12 | E2 docs; 0.16.1 or 0.17.0; app pin bump (E3, second half) |
| 13-15 | buffer to 2026-11-14 |

## Global Constraints

- Python 3.11 syntax and stdlib only; no new dependencies.
- `make check` green before every PR: black 24, isort, mypy on `flowfreq/` (new modules are mypy-enforced by default), pytest with `PYTHONSAFEPATH=1 PYTHONUTF8=1`.
- Never edit `vendor/`. Known-failing tests are `xfail(strict=True)` and announced.
- Every coefficient has a citation: an `nss_snapshot` equation cites the NSS citation record (author, title, URL, region id, citation id); a skew row cites a report with a DOI or URL that the agent opened.
- Notation: lowercase q for duration quantiles, capital Q or AEP in percent for floods.
- No `print()` in library code; `logging.getLogger(__name__)`.
- Docs regenerated, never hand-edited: `docs/REGRESSION_COVERAGE.md` via `tools/gen_regression_coverage.py` (`tests/test_catalog.py` fails when stale).
- Commit messages end with the attribution line in effect for the session.
- The app repository: modules at the root, never an `app/` directory; the pin is a tag, never `main`; the user merges and deploys app PRs.

## Review Focus

1. A basin whose located NSS regions none of the library's regions match (a verified state whose `region_name` lacks the NSS code): `estimate_at_site` must raise `EquationsUnavailable` naming the NSS codes, not return an empty table. Pinned in B2.
2. An `nss_snapshot` equation whose `asep_pct` is 0 or missing: `sep_log` stays `None` and the estimate's `interval_method` says there is no published standard error, rather than a zero-width interval. Pinned in A4.
3. A characteristic StreamStats reports as unavailable (`-999`) that an equation needs: `regression.evaluate` raises `KeyError`; `estimate_at_site` must report the missing code per region in `provenance["missing_characteristics"]` and skip that region, not crash. Pinned in B2.
4. A jurisdiction with no library file at all (DC, VI, GU, MP, AS): the error message must say NSS serves no peak-flow regions there, so the app shows the right reason. Pinned in B2 and D2.
5. A negative `scale` on an affine base (`(-1*(BDF-13))`) evaluated at BDF above 13: the argument of `log10` is non-positive and must raise a clear `ValueError` naming the code, not return NaN. Pinned in A3.

---

## Workstream A: parser extensions and the `nss_snapshot` tier

### Task A1: Identifiers with lowercase letters

**Files:**
- Modify: `flowfreq/regression/nss.py:96` (`_VARIABLE`)
- Test: `tests/test_nss_snapshot.py`

**Interfaces:**
- Produces: `parse_equation` accepts variable codes matching `^[A-Za-z][A-Za-z0-9_]*$`, excluding the function names in `_FUNCTIONS` and `logN`/`LogN` (A2).

- [ ] **Step 1: Write the failing tests**

```python
@pytest.mark.parametrize(
    "equation, codes",
    [
        ("10^2.183*DRNAREA^0.678*CSLBlue^0.239*BSHAPELFP^(-0.324)", ("DRNAREA", "CSLBlue", "BSHAPELFP")),
        ("207*DRNAREA^0.676*CSL10_85fm^0.159*BSHAPELFP^-0.390", ("DRNAREA", "CSL10_85fm", "BSHAPELFP")),
        ("5.18*DRNAREA^0.694*(IMPERV+1)^(0.382)*(SOILCorD+1)^(0.414)", ("DRNAREA", "IMPERV", "SOILCorD")),
    ],
)
def test_mixed_case_identifiers_parse(equation, codes):
    parsed = parse_equation(equation)
    assert parsed.codes == codes
    assert parsed.transforms[1] in ("log10", "log10_plus1")
```

and the round trip: `evaluate_expression(equation, inputs)` equals `10 ** parsed.log10_value(inputs)` to `rel=1e-12` at `DRNAREA=50, CSLBlue=20, BSHAPELFP=3, CSL10_85fm=20, IMPERV=5, SOILCorD=30`.

- [ ] **Step 2: Run them, expect `NSSEquationError: unsupported identifier`**

Run: `pytest tests/test_nss_snapshot.py -k mixed_case -v`

- [ ] **Step 3: Relax `_VARIABLE` and keep function names out of it**

`_VARIABLE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")`; `_check_variable` also refuses names in `_FUNCTIONS` and `{"logN", "LogN"}` so a function used as a variable still fails loudly.

- [ ] **Step 4: Run the whole parser module, expect PASS**

Run: `pytest tests/test_nss_snapshot.py -v`

- [ ] **Step 5: Commit** `feat(nss): accept mixed-case NSS variable codes`

### Task A2: `logN(X, 10)` in the exponent of 10

**Files:**
- Modify: `flowfreq/regression/nss.py:426-466` (`_linear`), `_FUNCTIONS`
- Test: `tests/test_nss_snapshot.py`

**Interfaces:**
- Produces: inside `10^( ... )`, a term `c*logN(E, b)` (`logN` or `LogN`) becomes a `log10` term on the affine base `E` (same rules as `_base`) with coefficient `c / log10(b)`; a term `c*k*(X/d)` becomes an `identity` term with coefficient `c*k/d` and a note `"constant factors folded into the coefficient of X"`.

- [ ] **Step 1: Write the failing tests**

```python
MN = "10^(24.074+0.744*logN(DRNAREA,10)-0.297*logN((LAKEAREA+1),10)-13.720*logN(LAT_OUT,10))"
WI = "10^(0.9357+0.8398*logN(DRNAREA,10)-0.6685*logN(SSURGOKSAT,10)-0.3740*logN((LC01HERB+1),10)-0.5143*logN((LC01WATER+1),10)+0.4227*I24H10Y)"
MA = "10^(1.631+0.801*LogN(DRNAREA,10)+0.589*0.001*(ELEV/3.28084)-1.137*0.01*(LC06STOR))"

def test_logN_terms_parse_to_log10():
    p = parse_equation(MN)
    assert p.intercept == pytest.approx(24.074)
    assert p.codes == ("DRNAREA", "LAKEAREA", "LAT_OUT")
    assert p.transforms == ("log10", "log10_plus1", "log10")
    assert p.coefficients == pytest.approx((0.744, -0.297, -13.720))

def test_identity_terms_with_folded_constants():
    p = parse_equation(MA)
    i = p.codes.index("ELEV")
    assert p.transforms[i] == "identity"
    assert p.coefficients[i] == pytest.approx(0.589 * 0.001 / 3.28084)
    assert any("folded" in n for n in p.notes)
```

plus the WI round trip to `rel=1e-12` at in-range inputs, and `logN(X, 2)` giving coefficient `c / log10(2)`.

- [ ] **Step 2: Run them, expect FAIL with `exponent of 10 is not linear in the variables`**
- [ ] **Step 3: Implement `logN` handling in `_linear`**

Add `"logN"` and `"LogN"` to a `_LOG_FUNCTIONS` set; in `_linear`, a `call` node of that name with two arguments whose second is a positive constant contributes a log10 term via `_base(args[0])` scaled by `1/log10(b)`; constant products around a variable fold into its coefficient. `_Linear` grows a fourth mapping, code to `(transform, scale, offset)`, so `_loglin`'s `10^` branch builds `log10`/`log10_plus1`/`identity` terms from it.

- [ ] **Step 4: Run the parser module, expect PASS**
- [ ] **Step 5: Commit** `feat(nss): parse logN terms in the exponent of 10`

### Task A3: Negated affine base

**Files:**
- Modify: `flowfreq/regression/nss.py:513-546` (`_base`), `flowfreq/regression/equations.py` (`Variable.__post_init__` and `Variable.apply`)
- Test: `tests/test_nss_snapshot.py`, `tests/test_regression_equations.py` (create if absent)

**Interfaces:**
- Produces: `(-1*(X-13))^p` and `(13-X)^p` parse to `Variable(code="X", transform="log10", scale=-1.0, offset=13.0)`; `Variable.apply` raises `ValueError(f"{code}: log10 argument {arg:g} is not positive")` when the transformed argument is non-positive.

- [ ] **Step 1: Write the failing tests**

```python
FL = "214*(DRNAREA)^(1.13)*(-1*(BDF-13))^(-0.59)*(CSL10_85)^(0.73)"

def test_negated_affine_base():
    p = parse_equation(FL)
    i = p.codes.index("BDF")
    assert p.transforms[i] == "log10"
    assert (p.scales[i], p.offsets[i]) == (-1.0, 13.0)
    assert p.coefficients[i] == pytest.approx(-0.59)

def test_negated_base_round_trip():
    inputs = {"DRNAREA": 10.0, "BDF": 4.0, "CSL10_85": 20.0}
    assert 10 ** parse_equation(FL).log10_value(inputs) == pytest.approx(evaluate_expression(FL, inputs), rel=1e-12)

def test_negative_scale_out_of_domain_names_the_code():
    v = Variable(code="BDF", transform="log10", scale=-1.0, offset=13.0)
    with pytest.raises(ValueError, match="BDF: log10 argument"):
        v.apply(14.0)
```

- [ ] **Step 2: Run, expect FAIL with `negated factor` / `non-positive scale`**
- [ ] **Step 3: Implement**

In `_base`, treat a `neg` node and a `bin *` with a negative constant as an affine base and drop the `scale > 0` check (keep `scale != 0`). In `Variable`, allow a negative `scale` and raise the named `ValueError` from `apply` for a non-positive log argument.

- [ ] **Step 4: Run `pytest tests/test_nss_snapshot.py tests/test_regression_equations.py -v`, expect PASS**
- [ ] **Step 5: Commit** `feat(nss): parse negated affine bases such as (-1*(BDF-13))^p`

### Task A4: The `nss_snapshot` status and `snapshot_to_library`

**Files:**
- Modify: `flowfreq/regression/library.py` (STATUSES, OPT_IN_SUBDIRS, `StateLibrary`, `_load_opt_in`, `load_state`, listings), `flowfreq/regression/__init__.py` (exports)
- Create: `flowfreq/regression/nss_library.py`
- Test: `tests/test_regression_nss_snapshot.py`

**Interfaces:**
- Produces:
  - `NSS_SNAPSHOT_SUBDIR = "nss_snapshot"`; `"nss_snapshot"` in `STATUSES` and `OPT_IN_SUBDIRS`.
  - `StateLibrary.unparsed: List[Dict[str, Any]]` (JSON key `unparsed`, default empty); `StateLibrary.nss_snapshot -> bool`; `StateLibrary.region_for_nss_code(code: str) -> Optional[str]`.
  - `load_state(code, data_dir=None, *, allow_provisional=False, allow_awaiting_nss=False, allow_nss_snapshot=False)`; precedence awaiting_nss, provisional, nss_snapshot.
  - `nss_snapshot_states(data_dir=None) -> List[str]`.
  - `nss_library.snapshot_to_library(snapshot: Dict[str, Any]) -> Optional[Dict[str, Any]]` returning the JSON form of a `StateLibrary` (schema_version 1, status `nss_snapshot`), or `None` when the snapshot has no regions.
  - `nss_library.asep_to_sep_log(asep_pct: float) -> float` = `sqrt(ln(1 + (asep/100)^2)) / ln(10)`.
  - The advisory text, verbatim from the spec, as `_nss_snapshot_advisory(code, lib)`.

- [ ] **Step 1: Write the failing tests** (model on `tests/test_regression_awaiting_nss.py`)

```python
def test_status_is_defined(): assert "nss_snapshot" in STATUSES
def test_default_lookup_ignores_the_snapshot(data_dir): ...          # EquationsUnavailable when no published file
def test_opt_in_returns_it_with_advisory(data_dir, caplog): ...        # advisory startswith "NSS SNAPSHOT (not report-verified)"
def test_published_set_wins_without_the_flag(...): ...
def test_awaiting_nss_wins_over_snapshot_when_both_opted_in(...): ...
def test_snapshot_status_cannot_replace_the_published_file(tmp_path): ...  # ValueError "belongs under nss_snapshot/"
def test_unparsed_is_loaded_and_counted_in_advisory(...): ...          # "2 published equations or regions"
def test_region_for_nss_code_matches_code_then_name(): ...             # "GC2015" -> "GC2015"; "(NSS GC1239)" in name -> "1"
def test_listings_keep_the_sets_apart(data_dir): ...                   # nss_snapshot_states(data_dir) == ["TN"]
```

and for the generator, on a trimmed copy of the TN snapshot (two regions, three equations, one with `parse_error`, one region `unresolved`) saved as `tests/fixtures/nss_snapshot_TN_trimmed.json`:

```python
def test_snapshot_to_library_builds_equations_and_unparsed(trimmed):
    lib = parse_state(snapshot_to_library(trimmed))
    assert lib.status == "nss_snapshot" and lib.state == "TN"
    eq = lib.equation("GC2015", 0.5)
    assert eq.intercept == pytest.approx(0.6866362692622934)
    assert eq.citation.table == "NSS region 1062 (GC2015), citation 210"
    assert eq.citation.url.startswith("http")
    assert eq.sep_log == pytest.approx(asep_to_sep_log(46.7))
    assert [u["reason"] for u in lib.unparsed] == ["<parse_error text>", "unresolved: <reason text>"]

def test_asep_zero_or_missing_leaves_sep_log_unset(trimmed): ...  # Review Focus 2
def test_snapshot_without_regions_returns_none(): assert snapshot_to_library({"state": "DC", "regions": [], ...}) is None
def test_non_aep_statistics_are_skipped_and_counted(trimmed): ...  # notes mention the count
```

- [ ] **Step 2: Run `pytest tests/test_regression_nss_snapshot.py -v`, expect ImportError / AttributeError**
- [ ] **Step 3: Implement the status in `library.py` and the generator in `nss_library.py`**

`snapshot_to_library` uses `nss.parse_equation` on each resolved region's `equation` string (not the stored `parsed` dict, so A1-A3 apply), `nss.to_regression_equation(...)` with the Citation from the spec, then `dataclasses.replace(eq, sep_log=asep_to_sep_log(asep))` when `asep_pct` is present and positive; `source_reports` and `notes` per the spec; equations are serialised with `RegressionEquation.to_dict()`.

- [ ] **Step 4: Run the test module and `tests/test_regression_awaiting_nss.py`, expect PASS**
- [ ] **Step 5: Commit** `feat(regression): nss_snapshot opt-in status and snapshot_to_library`

### Task A5: Indiana gating and the generator tool

**Files:**
- Modify: `tools/snapshot_nss.py:199-217` (`_gating_candidates`)
- Create: `tools/build_nss_snapshot_library.py`
- Test: `tests/test_nss_snapshot.py` (gating), `tests/test_regression_nss_snapshot.py` (tool main)

**Interfaces:**
- Consumes: A4's `snapshot_to_library`.
- Produces: `python tools/build_nss_snapshot_library.py --snapshot <file> [--snapshot ...] --out flowfreq/data/regression/nss_snapshot` writes `<CODE>.json` per snapshot with regions, skips a jurisdiction that has `flowfreq/data/regression/<CODE>.json` (prints `skipped CODE: published set exists`), skips one with no regions (prints `skipped CODE: no NSS peak-flow regions`), and prints `CODE: N equations, M unparsed`. Exit code 0.

- [ ] **Step 1: Write the failing tests**

```python
def test_gating_candidates_include_a_bounded_region_id_parameter():
    region = {"code": "GC1005", "parameters": [{"code": "HIGHREG", "limits": {"min": 1005, "max": 1012}}, {"code": "CONTDA", "limits": {"min": 0.27, "max": 13706}}]}
    assert {"HIGHREG": 1005.0} in _gating_candidates(region)

def test_build_tool_writes_skips_and_reports(tmp_path, capsys): ...  # TN trimmed -> written; WY -> skipped (published); DC -> skipped (no regions)
```

- [ ] **Step 2: Run, expect FAIL**
- [ ] **Step 3: Implement**

In `_gating_candidates`, a bounded parameter whose limits contain the region's GC number is also tried at that number, first. The tool is `argparse` over `snapshot_to_library`, writing `json.dumps(..., indent=1)` with a trailing newline.

- [ ] **Step 4: Run both test modules, expect PASS; then `make check`**
- [ ] **Step 5: Commit** `feat(tools): build the nss_snapshot library from snapshots; Indiana gating`

### Task A6: Harvest and commit the snapshots and the 41 library files

**Files:**
- Create: `data/nss_snapshots/<CODE>_2026-10-2N.json` for the 41 jurisdictions with regions, plus the 5 with none
- Create: `flowfreq/data/regression/nss_snapshot/<CODE>.json` (41)
- Modify: `tools/gen_regression_coverage.py`, `docs/REGRESSION_COVERAGE.md`
- Test: `tests/test_regression_nss_snapshot.py`, `tests/test_catalog.py`

**Interfaces:**
- Consumes: A1-A5. The re-harvest must run after A1-A3 so the `parsed` fields reflect the new parser; `snapshot_to_library` re-parses anyway.
- Produces: `nss_snapshot_states()` returns the 41 codes; coverage rows read `nss snapshot (N eqs, M unparsed)`.

- [ ] **Step 1: Re-run the sweep into `data/nss_snapshots/`**

Run (venv on PATH, `PYTHONUTF8=1`): `for s in CA AK HI MN WI IA IL MO MI ND SD NE KS OH IN KY TN WV AR LA MS OK TX ME NH VT MA RI CT NY NJ PA DE MD DC VA NC SC GA FL AL PR VI GU MP AS; do python tools/snapshot_nss.py --state $s; done`
Expected: IN now resolves 8 regions; refused totals drop from 220 to 24 (GA/NC/SC only).

- [ ] **Step 2: Build the library files**

Run: `python tools/build_nss_snapshot_library.py --snapshot data/nss_snapshots/*_2026-10-*.json --out flowfreq/data/regression/nss_snapshot`
Expected: 41 written, 5 skipped for no regions, 10 skipped for a published set.

- [ ] **Step 3: Write the data test**

```python
@pytest.mark.parametrize("code", nss_snapshot_states())
def test_every_snapshot_library_loads_and_evaluates(code):
    lib = load_state(code, allow_nss_snapshot=True)
    assert lib.nss_snapshot and lib.equations
    eq = lib.equations[0]
    mid = {v.code: math.sqrt(v.minimum * v.maximum) if v.minimum and v.minimum > 0 else (v.minimum + v.maximum) / 2 for v in eq.variables}
    assert evaluate(eq, mid).flow_cfs > 0

def test_snapshot_states_are_the_41(): assert len(nss_snapshot_states()) == 41
def test_unparsed_total_is_only_the_region_indicator_forms(): ...  # sum of len(lib.unparsed) == 24 + unresolved count
```

- [ ] **Step 4: Teach `gen_regression_coverage.py` the status and regenerate**

Run: `python tools/gen_regression_coverage.py && pytest tests/test_catalog.py tests/test_regression_nss_snapshot.py -v`
Expected: PASS; the coverage doc shows 41 `nss snapshot` rows.

- [ ] **Step 5: Commit** `data(regression): NSS snapshot library for 41 jurisdictions`

### Task A7: Docs for workstream A

**Files:**
- Modify: `CHANGELOG.md` (Unreleased, Added), `docs/MASTER_ROADMAP.md` §3.2 (status legend, the wave table's 46 rows), `TODO.md` (new entry under Open Items), `README.md` (coverage paragraph)

- [ ] **Step 1: Write the entries** (CHANGELOG: the status, the flag, the generator, the parser extensions, the counts; roadmap: `nss_snapshot` defined beside `verified`/`partial`/`pending`)
- [ ] **Step 2: `make check`, expect green**
- [ ] **Step 3: Commit** `docs: nss_snapshot tier and the 41-jurisdiction sweep`

---

## Workstream B: `estimate_at_site`

### Task B1: `SiteEstimate` and the pure assembly from pre-fetched inputs

**Files:**
- Create: `flowfreq/ungaged.py`
- Test: `tests/test_ungaged.py`

**Interfaces:**
- Consumes: `regression.evaluate`, `regression.evaluate_weighted`, `regression.library.StateLibrary.region_for_nss_code` (A4), `streamstats.RegressionRegionLocation`.
- Produces:
  - `@dataclass SiteEstimate` with the fields in the spec and `to_markdown() -> str`.
  - `assemble(lib: StateLibrary, regions: Sequence[RegressionRegionLocation], characteristics: Mapping[str, float], *, aeps: Optional[Sequence[float]], allow_extrapolation: bool, weighting_space: str) -> Tuple[pd.DataFrame, Dict[float, List[Tuple[float, RegressionEstimate]]], List[RegressionRegionLocation], Dict[str, List[str]]]` returning the quantile frame, the per-AEP estimates, the unmatched regions and the missing characteristics by region code. Pure: no network.

- [ ] **Step 1: Write the failing tests** on a synthetic `StateLibrary` with two regions (`GC1`, `GC2`) at AEPs 0.5, 0.1, 0.01, one variable `DRNAREA`:

```python
def test_single_region_keeps_its_interval(): ...                     # lower/upper from the estimate; interval_method is the estimate's
def test_two_regions_area_weighted_in_log_space(): ...               # flow == 10**(0.6*log q1 + 0.4*log q2); interval None; interval_method == "none: area-weighted across 2 regions"
def test_fractions_renormalise_after_an_unmatched_region(): ...     # 60/30/10 with GC3 unmatched -> 2/3, 1/3; GC3 in unmatched
def test_no_matching_region_raises_naming_the_codes(): ...           # Review Focus 1: EquationsUnavailable, "GC9" in message
def test_missing_characteristic_skips_the_region_and_records_it(): ...  # Review Focus 3
def test_extrapolated_flag_follows_out_of_range(): ...
def test_default_aeps_are_the_intersection_largest_first(): ...
def test_to_markdown_names_regions_fractions_and_citations(): ...
```

- [ ] **Step 2: Run `pytest tests/test_ungaged.py -v`, expect ImportError**
- [ ] **Step 3: Implement `SiteEstimate`, `assemble` and `to_markdown`**
- [ ] **Step 4: Run, expect PASS; `mypy flowfreq/ungaged.py` clean**
- [ ] **Step 5: Commit** `feat(ungaged): SiteEstimate and region assembly`

### Task B2: `estimate_at_site` with StreamStats and the library

**Files:**
- Modify: `flowfreq/ungaged.py`, `flowfreq/__init__.py` (export `estimate_at_site`, `SiteEstimate`)
- Test: `tests/test_ungaged.py`

**Interfaces:**
- Consumes: `streamstats.delineate_and_get_characteristics`, `streamstats.locate_regression_regions`, `regression.library.load_state`, `gage_weighting.weight_with_regression`, B1's `assemble`.
- Produces: `estimate_at_site(lat: float, lon: float, *, state: str, aeps: Optional[Sequence[float]] = None, at_site: Optional[FrequencyResults] = None, allow_nss_snapshot: bool = False, allow_provisional: bool = False, allow_extrapolation: bool = False, weighting_space: str = "log", cache: Optional[StreamStatsCache] = None, timeout: float = 60.0) -> SiteEstimate`.

- [ ] **Step 1: Write the failing tests**, monkeypatching `flowfreq.ungaged.delineate_and_get_characteristics` and `flowfreq.ungaged.locate_regression_regions` with fakes returning a `WatershedCharacteristics` built from `docs/vignettes/wave1_wave2_sites.json`'s WA entry (DRNAREA 82.2, PRECPRIS10 21.05, CANOPY_PCT 17.698; region weights 1: 0.0092, 4: 0.9908):

```python
def test_washington_site_uses_the_verified_library_and_two_regions(fakes): ...  # library_status == "verified"; regions == ["1","4"] mapped from GC codes
def test_snapshot_state_needs_the_flag(fakes, tmp_data_dir): ...               # TN without flag -> EquationsUnavailable; with flag -> advisory in result
def test_no_library_file_names_nss_coverage(fakes): ...                        # Review Focus 4: DC -> message contains "NSS serves no peak-flow regions"
def test_at_site_results_are_weighted_in_a_single_region(fakes, big_sandy_results): ...  # method == "at-site weighted with regression"; weighted frame has variance_log
def test_at_site_in_two_regions_is_not_weighted(fakes, big_sandy_results): ...           # method == "regression"; provenance says why
def test_provenance_carries_snap_delineation_and_status(fakes): ...
```

- [ ] **Step 2: Run, expect FAIL**
- [ ] **Step 3: Implement `estimate_at_site`**

The "no library file" message: `f"{code}: no regression library file. NSS {'serves' if code in NSS_SERVED else 'serves no'} peak-flow regions for this jurisdiction"` where `NSS_SERVED` is the set of codes with a published or snapshot file plus the ten verified states; keep it as a module constant computed from `available_states()` and `nss_snapshot_states()` at call time.

- [ ] **Step 4: Run, expect PASS; `make check`**
- [ ] **Step 5: Commit** `feat(ungaged): estimate_at_site`

### Task B3: Live check and vignette

**Files:**
- Create: `tests/test_ungaged_live.py` (`requires_network`), `docs/vignettes/ungaged_nationwide.py`

- [ ] **Step 1: Write the live test**: Little Klickitat (45.844124, -120.796392, WA) through `estimate_at_site`; assert the 1 % AEP is within a factor of 2 of `b17c_q1pct` 7183.2 from the sites file, and the regions are 1 and 4 with 4 above 0.95. A second case: a Tennessee point (36.1627, -86.7816, TN; Nashville) with `allow_nss_snapshot=True` returns a positive 1 % AEP and an advisory.
- [ ] **Step 2: Run with `PYTHONSAFEPATH=1 PYTHONUTF8=1 pytest tests/test_ungaged_live.py -m requires_network -v`, expect PASS** (record the values in the test's docstring)
- [ ] **Step 3: Write the vignette**: one WA site and one TN site, printing `to_markdown()`.
- [ ] **Step 4: Commit** `test(ungaged): live check on WA and TN; vignette`

### Task B4: CLI `ungaged`

**Files:**
- Modify: `flowfreq/cli.py`
- Test: `tests/test_cli.py`

**Interfaces:**
- Produces: `flowfreq ungaged --lat FLOAT --lon FLOAT --state CODE [--aep FLOAT ...] [--allow-nss-snapshot] [--allow-extrapolation]` printing `SiteEstimate.to_markdown()`; exit 1 with the message on `EquationsUnavailable`, `UnsnappablePointError`, `UnsupportedRegionError`.

- [ ] **Step 1: Write the failing test** with `click.testing.CliRunner` and `estimate_at_site` monkeypatched to return a fixed `SiteEstimate`; and one for the exit-1 path.
- [ ] **Step 2: Run, expect FAIL**
- [ ] **Step 3: Implement the command**
- [ ] **Step 4: Run `pytest tests/test_cli.py -v`, expect PASS; `make check`**
- [ ] **Step 5: Commit** `feat(cli): ungaged command`

---

## Workstream C: regional skew rows for every jurisdiction

### Task C1: Transcribe one `pending` row per jurisdiction (four agents)

**Files:**
- Modify: `flowfreq/data/regional_skew.csv`
- Test: `tests/test_regional_skew.py` (existing `load_table` validation)

Agent split: (a) Wave 3 + 4a + 4b: CA AK HI MN WI IA IL MO MI ND SD NE KS; (b) 4c + 4d: OH IN KY TN WV AR LA MS OK TX; (c) Wave 5 north: ME NH VT MA RI CT NY NJ PA DE MD DC; (d) Wave 5 south + 6: VA NC SC GA FL AL PR VI GU MP AS. Each agent works on its own branch and the rows are merged into one PR.

**Interfaces:**
- Produces: for each code, a row `state, skew_region, "", "", "", citation, table, pending, notes`. `skew_region` is the report's region name or `"<State> (statewide)"`. `citation` is the full USGS citation with DOI or URL **that the agent opened**. `table` names where the skew is printed. `notes` records the printed G, MSE (or SE) and effective record length, whether the report used a B-GLS study, the Bulletin 17B map, or a station skew, and the year of data. If no report could be opened, `notes` says `"citation not opened: <URL> <error>"` and the row stays pending.

- [ ] **Step 1: For each jurisdiction, find the current USGS peak-flow report** (NSS citation records in `data/nss_snapshots/<CODE>_*.json` give the report each equation set cites; start there), open the PDF or HTML, find the regional skew section.
- [ ] **Step 2: Append the row**; run `python -c "from flowfreq.regional_skew import load_table; load_table()"` after each.
- [ ] **Step 3: Run `pytest tests/test_regional_skew.py -v`, expect PASS**
- [ ] **Step 4: Commit per agent** `data(skew): pending rows for <codes>`

### Task C2: Second-check and verify

**Files:**
- Modify: `flowfreq/data/regional_skew.csv`, `docs/REGRESSION_COVERAGE.md` (regenerated)
- Test: `tests/test_regional_skew.py`

A different agent from C1. For every C1 row: open the cited source, confirm the printed values in `notes`. Where the report prints one statewide (or per-region) B-GLS skew with its MSE and the values match, fill `skew`, `skew_mse`, `effective_record_years`, set `status` to `verified`, and add `"second-checked 2026-11-0N"` to `notes`. Where the report says to use the Bulletin 17B map, the row stays `pending` with that stated. Where values disagree, the row stays `pending` with both readings in `notes`.

- [ ] **Step 1: Check and edit each row**
- [ ] **Step 2: Add a test** `test_every_jurisdiction_has_a_skew_row()` asserting all 56 codes appear in `load_table()["state"]`.
- [ ] **Step 3: Run `python tools/gen_regression_coverage.py && pytest tests/test_regional_skew.py tests/test_catalog.py -v`, expect PASS**
- [ ] **Step 4: Commit** `data(skew): second-checked rows; N verified`

---

## Workstream D: the app (`D:/a/flowfreq-app`, the user merges)

### Task D1: Revive the pin bump at v0.15.0

**Files:**
- Branch `bump-flowfreq-v0.15.0` from `origin/main`; cherry-pick `b21fc5e`, `1b7a9b2`, `2d9ef33` from `bump-flowfreq-v0.10.0`
- Modify: `requirements.txt`, `CLAUDE.md`, `README.md`, `CHANGELOG.md` (pin and wording: v0.15.0)

- [ ] **Step 1: Create the branch and cherry-pick; resolve conflicts keeping the branch's app code**
- [ ] **Step 2: Set the pin to `v0.15.0`; `pip install -r requirements.txt -r requirements-dev.txt`; `pip show flowfreq` says 0.15.0**
- [ ] **Step 3: Run `make check`; fix any drift** against the 0.15 API (`run_ffa` signature in `flowfreq/workflow.py:306-320`; `RegionalSkewUnavailable`; `peak_code_kwargs`)
- [ ] **Step 4: Smoke with `streamlit.testing.v1.AppTest`** on 12449500 (published study) and 01426500 (regulated refusal), as PR #5 did; record the outputs in the PR body
- [ ] **Step 5: Open the PR for the user**; do not merge

### Task D2: Ungaged-site page (after 0.16.0)

**Files:**
- Create: `ungaged_ui.py`
- Modify: `streamlit_app.py` (sidebar mode radio; the ungaged branch calls `ungaged_ui.render()`), `requirements.txt` (pin `v0.16.0`), `README.md`, `CLAUDE.md`, `CHANGELOG.md`
- Test: `tests/test_ungaged_ui.py`

**Interfaces:**
- Consumes: `flowfreq.estimate_at_site`, `flowfreq.ungaged.SiteEstimate`, `flowfreq.regression.jurisdictions.JURISDICTIONS`, `flowfreq.regression.available_states`, `flowfreq.regression.nss_snapshot_states`, the StreamStats and regression exceptions named in the spec.
- Produces: `render() -> None`; `coverage_badge(code: str) -> str` returning `"verified"`, `"partial"`, `"nss snapshot (not report-verified)"` or `"none: NSS serves no peak-flow regions"`; `run_estimate(lat, lon, state, aeps, allow_nss_snapshot) -> Dict[str, Any]` that never raises (`{"estimate": SiteEstimate}` or `{"error": str, "kind": str}`).

- [ ] **Step 1: Write the failing tests** in bare mode (`app_module` fixture pattern): `coverage_badge` for WA, TN, DC; `run_estimate` with `estimate_at_site` monkeypatched to raise each exception kind returns the error dict; the download payload equals `to_markdown()`; the advisory is shown when `library_status == "nss_snapshot"`.
- [ ] **Step 2: Run, expect FAIL**
- [ ] **Step 3: Implement `ungaged_ui.py` and the sidebar switch**
- [ ] **Step 4: `make check`; AppTest smoke on a WA point and a TN point with the flag; record outputs**
- [ ] **Step 5: Open the PR for the user**

---

## Workstream E: evidence, docs, releases

### Task E1: Nationwide smoke run

**Files:**
- Create: `tools/nationwide_smoke.py`, `docs/NATIONWIDE_SMOKE.md` (generated)
- Test: `tests/test_nationwide_smoke.py` (the selection and the table writer, offline)

**Interfaces:**
- Produces: `select_site(catalog: pd.DataFrame, state: str) -> Optional[pd.Series]` (longest `n_years`, `regulation_class == "reference"` first, else any unregulated); `smoke_one(row, *, allow_nss_snapshot=True) -> Dict[str, Any]` with keys `state, site_no, ffa_ok, ffa_q1pct, ungaged_ok, ungaged_q1pct, library_status, error`; `write_report(rows, path)`.

- [ ] **Step 1: Write the offline tests** for `select_site` and `write_report`
- [ ] **Step 2: Implement; run offline tests**
- [ ] **Step 3: Run live** (`USGS_API_KEY` from `~/.usgs_api_key`, never printed): `python tools/nationwide_smoke.py --out docs/NATIONWIDE_SMOKE.md`
Expected: a 56-row table. Known gaps go in the doc's header: territories and DC (no library), AK/HI where StreamStats cannot delineate.
- [ ] **Step 4: File each unexpected failure as an issue or a fix task; rerun**
- [ ] **Step 5: Commit** `qa: nationwide smoke run 2026-11-NN`

### Task E2: Docs

**Files:**
- Modify: `README.md` (coverage summary with the three tiers and the smoke table link), `docs/MASTER_ROADMAP.md` §0 table and §4 (`estimate_at_site` done, transposition branch open), `TODO.md` (status, the parked methow-flow items, the GA/NC/SC interaction forms, the open region polygons), `CHANGELOG.md`
- Modify: memory file `flowfreq-state.md` at the end of the work

- [ ] **Step 1: Write; `make check`**
- [ ] **Step 2: Commit** `docs: nationwide coverage, ungaged workflow, parked items`

### Task E3: Releases

- [ ] **Step 1: 0.16.0** after A, B, C merge: branch `release/prep-0.16.0`, version, CHANGELOG section, README install line and history row (as #191); merge; tag `v0.16.0` on the merge commit; `release/0.16.x`; bump main to `0.17.0.dev0` with the CLAUDE.md lines.
- [ ] **Step 2: App pin bump PR** to `v0.16.0` (part of D2).
- [ ] **Step 3: After E1 fixes**: 0.16.1 from `release/0.16.x` (or 0.17.0 from main if a feature landed); app pin bump PR.
