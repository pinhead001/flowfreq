# flowfreq: agentops adoption and the 0.17 follow-ups — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Clear the follow-ups the 0.16.0 nationwide rollout left (issues #203–#218, TODO.md
"Nationwide rollout (2026-10-25 plan)") and ship 0.17.0, running the mechanical work through
the `agentops` lead-plus-delegate machinery (`D:/a/agentops`) so the Claude Code session
spends its quota on specs, data reading, review and release, not on typing code.

**Architecture:** flowfreq becomes an `agentops` consumer (`config/agents.toml`, hooks,
`AGENTS.md`, `run.bat`, task cards under `config/cards/`), with its own state dir so cards
never mix with the methow GU sprint. Work is split three ways: **delegate cards** (code with
a named test and no data or network), **lead items** (investigation, data transcription from
USGS reports, spec writing, review, merge, release) and **user decisions** (hydrologist
picks, app merges, blocked GitHub actions). The lead writes a card only after the spec in it
leaves the delegate exactly one reasonable implementation.

**Tech Stack:** agentops 0.1.0 (Python 3.12, uv, `agy` 1.4.10 as the delegate), flowfreq
(pip/setuptools, Python 3.11–3.14, `make check`), GitHub CLI.

**Spec:** this plan argues from `TODO.md` ("Open Items" → "Nationwide rollout"), GitHub issues
#203–#218, `docs/MASTER_ROADMAP.md` §4, and `D:/a/agentops/docs/REQUIREMENTS.md` (consumer
contract) plus `D:/a/agentops/docs/adr/0001-agy-integration.md` (delegate rules). There is no
separate design document; the decisions are recorded in the next section.

## Global Constraints

- flowfreq stays a pip/setuptools project: **no `[tool.uv]` in flowfreq's `pyproject.toml`**.
  agentops runs from its own venv (`D:/a/agentops/.venv`), which already exists. CI is unchanged.
- The runner's worktrees have no editable install. Every acceptance command must set
  `PYTHONPATH=<worktree>`, `PYTHONSAFEPATH=1`, `PYTHONUTF8=1` and call the main venv's
  interpreter `D:\a\flowfreq\.venv\Scripts\python.exe` (memory: worktree-editable-install).
- Delegates get code and synthetic fixtures only (ADR-0001): no card downloads data, reads a
  USGS report, or needs `USGS_API_KEY`. Anything that does is a lead item.
- Card prompts name every file the delegate may read and edit (`templates/AGENTS.md`). The
  delegate has no shell; it edits files and ends its turn.
- Mechanical cards use `gemini-3.8-flash-medium` (the config default); algorithmic cards set
  `model = "gemini-3.8-flash-high"`. Weekly delegate capacity is about 1.2M Gemini tokens
  (ADR-0001 addendum 1); keep cards small and name files narrowly.
- The lead's existing push-and-merge workflow stands (memory: merge-workflow). So flowfreq does
  **not** install agentops' `pre-push` block and its `[hooks] deny_bash` does **not** deny
  `git push`. The runner itself never pushes or merges regardless.
- Known-failing tests are `xfail(strict=True)`, never skipped (CLAUDE.md). A live-service
  outage is not a known failure; see Task 7.
- Python 3.11 syntax only in flowfreq; `black>=24,<25` formatting; isort; mypy on new modules.
- Release 0.17.0 comes from `main` (CHANGELOG Unreleased already has fixes). A 0.16.1 from
  `release/0.16.x` is not planned unless the user asks.

## Review Focus

Inputs the follow-ups imply but no issue spells out; each has a test pinned to its task.

1. **A `region=` list whose weights sum to less than 100 %** (one side of a boundary basin
   named, the other forgotten) must still be refused with the listing, not silently
   area-weighted over a partial basin. — Task 3, card `ff-multicode-region`.
2. **A `region=` list naming two sets of the same side** (weights sum well over 100 %) must
   be refused as overlapping, not averaged. — Task 3.
3. **The HUC8-states table with a HUC8 in `CN` or `MX` only** (border units the WBD tags with
   a non-US code) must not be read as "unknown state". — Task 4, card `ff-huc8-states`.
4. **Silencing the advisory must leave `StateLibrary.advisory` set**, so a report still
   shows the status even when the sweep did not log it. — Task 2, card `ff-advisory-logger`.
5. **A worktree acceptance run importing the main checkout's `flowfreq`** would pass against
   the wrong code. — Task 1, `run.bat accept` prints `flowfreq.__file__` and fails unless it
   is under the current directory.

---

## Next steps inventory (what this plan sequences)

| # | Item | Source | Owner | Task |
|---|---|---|---|---|
| 1 | Make flowfreq an agentops consumer | agentops REQUIREMENTS | lead | 1 |
| 2 | Advisory warning silenceable for sweeps | #214 | card | 2 |
| 3 | Multi-code `region=` for boundary basins | #203 | card (flash-high) | 3 |
| 4 | Nationwide HUC8-to-states table | #205 | card (script) + lead (run it) | 4 |
| 5 | Catalog defect: 15388960 under DC | #206 | lead (find cause) → card | 5 |
| 6 | Smoke script snap distance / manual coordinates | #210 | lead (API check) → card | 5 |
| 7 | NSS region mismatches WY, NC, NV | #207 | lead | 5 |
| 8 | GA/NC/SC region-indicator equations (24) | #211 | lead (spec) → card (flash-high) | 5 |
| 9 | SC GC1587 range transcription | #208 | lead | 6 |
| 10 | `sep_log` on 32 LA/AR/MS equations | #213 | lead | 6 |
| 11 | HUC membership rows for 15 skew studies | #204 | lead | 6 |
| 12 | Per-state default for overlapping regions; MT and NM picks | #209 | **user** | 6 |
| 13 | Live-service test failure (WA PRECIP raster) | #218 | lead | 7 |
| 14 | Dependabot #219, #220, #221 | PRs | lead | 7 |
| 15 | Local `D:/a/flowfreq-app` checkout is behind origin (pins v0.10.1; origin v0.16.0) | observed | lead | 7 |
| 16 | v0.10.4 GitHub release recreation | TODO small items | **user** (`gh release create` is classifier-blocked) | 7 |
| 17 | Release 0.17.0; app pin bump PR | flowfreq-state memory, E3 step 3 | lead; **user** merges the app PR | 8 |
| — | NY urban regions (#212), 41 `pending` skew rows, methow-flow tiers 2–3, roadmap §4 donor search and submittal report | TODO | parked | — |

---

### Task 1: Make flowfreq an agentops consumer

**Files:**
- Create: `config/agents.toml`, `run.bat`, `AGENTS.md`, `.claude/agents/reviewer.md`,
  `.claude/agents/log-digest.md`, `.githooks/pre-commit`, `config/cards/README.md`
- Modify: `.claude/settings.json` (keep its `permissions`, `extraKnownMarketplaces`,
  `enabledPlugins`; add `env`, `statusLine`, `hooks`)
- Modify: `CLAUDE.md` (a short "Delegation" paragraph pointing at `run.bat` and this plan)
- Modify: `.gitignore` (no change needed: the state dir is outside the repo)

**Interfaces:**
- Produces: `run.bat gov|cards|ledger|statusline|hook ...` (routed to
  `uv run --project D:\a\agentops agentops --config %HERE%\config\agents.toml`);
  `run.bat accept <pytest paths...>` (lint + mypy + pytest in the current checkout with the
  Global Constraints environment); `run.bat check` (`make check` with the venv on PATH).
  Every card in Tasks 2–5 uses `acceptance = "run.bat accept <paths>"`.

- [ ] **Step 1: Sync agentops and confirm it runs**

Run: `uv sync --project D:/a/agentops && uv run --project D:/a/agentops agentops --help`
Expected: usage text listing `statusline`, `gov`, `cards`, `ledger`, `hook`.

- [ ] **Step 2: Write `config/agents.toml`** from `D:/a/agentops/templates/agents.toml` with:

```toml
[paths]
state_dir = "D:/a/sbx/state/agentops-flowfreq"
worktree_root = "D:/a/sbx/state/agentops-flowfreq/worktrees"
worktree_prefix = ""
[repos]
flowfreq = "D:/a/flowfreq"
```

Keep `[governor]` and `[agy]` as the template has them except
`pre_acceptance = []` (black/isort run inside `run.bat accept`, and the runner's cwd shell
has no `uv` project here) and `print_timeout = "25m"` (flowfreq's suite is slower than the
GU repos'). `[hooks]`: the template's `deny_bash` **minus** the `git push` pattern;
`deny_substrings = ["_restricted", "USGS_API_KEY", ".usgs_api_key"]`;
`teammate_tools = ["TeamCreate"]`; `task_acceptance = "run.bat check"`.

- [ ] **Step 3: Write `run.bat`** from `templates/run.bat`, changing `:agentops` to
`uv run --project "%AGENTOPS%" agentops --config "%AGENTOPS_CONFIG%" %*` and adding:

```bat
:accept
rem run.bat accept <pytest paths>  -- lint, mypy and pytest in THIS checkout (worktree-safe)
shift
set "PYTHONPATH=%HERE%"
set "PYTHONSAFEPATH=1"
set "PYTHONUTF8=1"
set "PY=D:\a\flowfreq\.venv\Scripts\python.exe"
"%PY%" -c "import flowfreq,os,sys; p=os.path.abspath(flowfreq.__file__); print(p); sys.exit(0 if p.startswith(os.path.abspath('%HERE%')) else 3)" || exit /b 3
"%PY%" -m black --check flowfreq tests || exit /b 1
"%PY%" -m isort --check-only flowfreq tests || exit /b 1
"%PY%" -m mypy flowfreq/ || exit /b 1
"%PY%" -m pytest -q -x %1 %2 %3 %4 %5 %6 %7 %8 %9
exit /b %ERRORLEVEL%

:check
set "PATH=D:\a\flowfreq\.venv\Scripts;%PATH%"
set "PYTHONPATH=%HERE%"
make check
exit /b %ERRORLEVEL%
```

- [ ] **Step 4: Verify `run.bat accept` in the main checkout and in a throwaway worktree**

Run: `run.bat accept tests/test_version.py`
Expected: prints `D:\a\flowfreq\flowfreq\__init__.py`, lint and mypy pass, `1 passed`.

Run: `git worktree add D:/a/sbx/state/agentops-flowfreq/worktrees/probe main && cd D:/a/sbx/state/agentops-flowfreq/worktrees/probe && run.bat accept tests/test_version.py; cd D:/a/flowfreq && git worktree remove D:/a/sbx/state/agentops-flowfreq/worktrees/probe`
Expected: the printed path is under `...\worktrees\probe\flowfreq\__init__.py` (Review Focus 5) and `1 passed`.

- [ ] **Step 5: Write `AGENTS.md`** from `templates/AGENTS.md`, replacing the last bullet's
toolchain line with: "Python 3.11-compatible syntax; the runner checks black, isort, mypy
and pytest with `run.bat accept`. NumPy docstrings on public functions, type hints on every
signature, `logging.getLogger(__name__)` not `print`, no bare `except`." Add one bullet:
"Never read `vendor/`, `flowfreq/data/*.csv.gz` or anything under `.cache/`."

- [ ] **Step 6: Merge the template `settings.json` into `.claude/settings.json`**

`env.AGENTOPS_CONFIG = "D:/a/flowfreq/config/agents.toml"`; `statusLine.command` and the
three hook commands use `D:/a/agentops/.venv/Scripts/python.exe -m agentops.cli --config
D:/a/flowfreq/config/agents.toml ...` (not `${CLAUDE_PROJECT_DIR}/.venv`, which is
flowfreq's venv without agentops). Keep the existing `permissions.allow` entries and add
`"Bash(run.bat *)"`, `"Bash(./run.bat *)"`, `"Bash(uv run --project D:/a/agentops *)"`.
Do not add the template's `deny` list (it includes `git push`).

- [ ] **Step 7: Copy `templates/agents/*.md` to `.claude/agents/` and `templates/githooks/pre-commit` to `.githooks/`**

Edit `.githooks/pre-commit` to run `D:\a\flowfreq\.venv\Scripts\python.exe -m black --check flowfreq tests && ... isort --check-only flowfreq tests`
(flowfreq has no pre-commit config). Do **not** copy `pre-push`. Run
`git config core.hooksPath .githooks`.

- [ ] **Step 8: Write `config/cards/README.md`**: cards are TOML per `agentops` `queue.py`
(`id`, `repo = "flowfreq"`, `title`, `prompt`, `acceptance`, optional `base`, `model`,
`depends_on`); queue one with `run.bat cards add config/cards/<id>.toml`; run with
`run.bat cards run --workers 1`; results in `D:/a/sbx/state/agentops-flowfreq/results/<id>/`;
the lead reviews `card/<id>` with `git diff main..card/<id>`, then opens the PR.

- [ ] **Step 9: Verify the governor and hooks**

Run: `run.bat gov status` → prints a governor line (mode `normal` or `unknown`), agy pools.
Run: `run.bat cards list` → empty queue, four state dirs created.
Run: `echo {} | run.bat statusline` → prints a status line without a traceback.
Start a new Claude Code session in the repo → the SessionStart hook prints the governor status.

- [ ] **Step 10: Add the CLAUDE.md paragraph and commit**

```bash
git checkout -b ops/agentops-consumer
git add config/agents.toml config/cards/README.md run.bat AGENTS.md .claude/settings.json .claude/agents .githooks/pre-commit CLAUDE.md docs/superpowers/plans/2026-10-29-agentops-adoption-and-0.17-follow-ups.md
git commit -m "ops: make flowfreq an agentops consumer (cards, governor, hooks)"
```

Open the PR, wait for CI (lint/test/parity unchanged by these files), merge.

---

### Task 2: Card `ff-advisory-logger` — silence the snapshot advisory for sweeps (#214)

**Files:**
- Create: `config/cards/ff-advisory-logger.toml`
- Delegate modifies: `flowfreq/regression/library.py:240-256`, `tools/nationwide_smoke.py`,
  `tests/test_regression_nss_snapshot.py`, `tests/test_nationwide_smoke.py`, `CHANGELOG.md`

**Interfaces:**
- Produces: module-level `advisory_logger = logging.getLogger(__name__ + ".advisory")` in
  `flowfreq/regression/library.py`; every advisory (`provisional`, `awaiting_nss`,
  `nss_snapshot`) is logged through it and nothing else is. `tools/nationwide_smoke.py`
  gains `--quiet-advisories` (sets that logger's level to ERROR).

- [ ] **Step 1: Write the card** with this prompt (abridged here; the file carries it in full):

> Read AGENTS.md first. Files you may read: `flowfreq/regression/library.py`,
> `tools/nationwide_smoke.py`, `tests/test_regression_nss_snapshot.py`,
> `tests/test_nationwide_smoke.py`, `CHANGELOG.md`. Files you may edit: the same five.
> Task: in `library.py` add `advisory_logger = logging.getLogger(__name__ + ".advisory")`
> next to `logger` and make the single `logger.warning("%s", lib.advisory)` call use it.
> Document in the `load_state` docstring Notes: "Advisories are logged on the
> `flowfreq.regression.library.advisory` logger; a sweep over many states silences them
> with `logging.getLogger("flowfreq.regression.library.advisory").setLevel(logging.ERROR)`
> and reads `StateLibrary.advisory` instead." Add tests to
> `tests/test_regression_nss_snapshot.py`: `test_advisory_logged_on_advisory_logger`
> (caplog at WARNING on that logger name captures exactly one record whose message equals
> `lib.advisory`) and `test_advisory_silenced_keeps_library_advisory` (with the logger
> level set to ERROR, caplog has no records and `lib.advisory` is non-empty). In
> `tools/nationwide_smoke.py` add `--quiet-advisories` (store_true) to `main`, applied
> before the loop, and one test in `tests/test_nationwide_smoke.py`,
> `test_quiet_advisories_flag_parses`, asserting `main(["--help"])` output lists it or that
> the parser accepts it. Add a CHANGELOG Unreleased line under Changed.

`acceptance = "run.bat accept tests/test_regression_nss_snapshot.py tests/test_nationwide_smoke.py"`

- [ ] **Step 2: Queue and run**

Run: `run.bat cards add config/cards/ff-advisory-logger.toml && run.bat cards run --workers 1 --card ff-advisory-logger`
Expected: card moves to `done`; `results/ff-advisory-logger/` has the envelope and an
acceptance log ending in `passed`.

- [ ] **Step 3: Lead review and PR**

Run: `git diff main..card/ff-advisory-logger` — check Review Focus 4 (advisory still set),
no other logger renamed, docstring wording. Then
`git push -u origin card/ff-advisory-logger && gh pr create --fill`, CI, merge. Close #214
is the user's (classifier-blocked); note it in the PR body with `Closes #214`.

---

### Task 3: Card `ff-multicode-region` — `region=` accepts one code per side (#203)

**Files:**
- Create: `config/cards/ff-multicode-region.toml` (`model = "gemini-3.8-flash-high"`)
- Delegate modifies: `flowfreq/ungaged.py` (`OverlappingRegionsError` ~145–200, `assemble`
  285–520, `estimate_at_site` 806+), `flowfreq/cli.py` (`ungaged --region`),
  `tests/test_ungaged.py`, `CHANGELOG.md`

**Interfaces:**
- Consumes: `assemble(..., region: Optional[str])` and `estimate_at_site(..., region: Optional[str])`,
  `PERCENT_WEIGHT_TOL`, `OverlappingRegionsError(state, candidates, region)`.
- Produces: `region: Optional[Union[str, Sequence[str]]]` on both; a string may hold
  comma-separated codes; `SiteEstimate.provenance["region_choice"]` is the list of codes used.

Semantics the card fixes (the lead's design decision):

1. Normalise `region` to a list of distinct codes (NSS or library codes, case-insensitive,
   as today). One code behaves exactly as now.
2. Every code must name a located region, else `OverlappingRegionsError(..., region=<that code>)`.
3. The chosen regions' `percent_weight` must sum to within `max(PERCENT_WEIGHT_TOL, 0.5*len(chosen))`
   of 100. Below → refused with the message "cover N % of the basin: name one set per side"
   (Review Focus 1). Above → refused as overlapping (Review Focus 2).
4. Chosen regions are then `used` with fractions `percent_weight/total`, as the existing
   several-region path does (`interval_method="none: area-weighted across N regions"`);
   `alternatives` lists the located regions not chosen.
5. CLI `--region` may be repeated or comma-separated; the prompt names the exact click
   option change (`multiple=True`).

- [ ] **Step 1: Write the card**, naming the files above, the five rules verbatim, and the tests
to add to `tests/test_ungaged.py` using the module's existing fixture style for located
regions with weights: `test_region_list_covering_both_sides_area_weights` (sides 60/40,
two overlapping sets per side; `region=["X1","Y1"]` → fractions 0.6/0.4, alternatives
`{X2, Y2}`); `test_region_list_partial_cover_refused` (`region=["X1"]` → `OverlappingRegionsError`,
message contains "60 %"); `test_region_list_same_side_refused` (`region=["X1","X2"]` →
`OverlappingRegionsError`); `test_region_comma_string_equals_list`; one CLI test
`test_cli_ungaged_region_repeated` in `tests/test_cli.py` using `CliRunner` with the
`flowfreq.ungaged.estimate_at_site` monkeypatched as the existing `ungaged` CLI tests at `tests/test_cli.py:400` do.

`acceptance = "run.bat accept tests/test_ungaged.py tests/test_cli.py"`

- [ ] **Step 2: Queue and run** (`--card ff-multicode-region`); on escalation the lead finishes
in the worktree rather than re-running (ADR-0001 §5, `max_retries = 1`).

- [ ] **Step 3: Lead review**: `region` typing is `Optional[Union[str, Sequence[str]]]` on
both public functions and the docstrings say so; provenance records the list; `README.md`
"Ungaged sites" gets the one-line example `--region GC729,GC1618`. PR, CI, merge
(`Closes #203`). Update `docs/MASTER_ROADMAP.md` §4 checkbox in the same PR.

---

### Task 4: Card `ff-huc8-states` — nationwide HUC8-to-states build (#205)

**Files:**
- Create: `config/cards/ff-huc8-states.toml`
- Delegate creates: `tools/build_huc8_states.py`, `tests/test_build_huc8_states.py`,
  `tests/fixtures/wbd_huc8_page.json` (synthetic, 3 features)
- Lead then runs the tool and commits `flowfreq/data/regional_skew_huc8_states.csv`

**Interfaces:**
- Consumes: the WBD map service `tools/build_regional_skew_hucs.py:78` already uses
  (`https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer`, HU8 layer 4,
  fields `huc8`, `name`, `states`), and the CSV schema `huc8,name,states` (states as a
  comma-joined string inside quotes, e.g. `"CN,MT"`).
- Produces: `fetch_huc8_pages(session, *, page_size: int = 1000) -> Iterator[dict]`
  (paginated `query?where=1=1&outFields=huc8,name,states&resultOffset=…&f=json`),
  `rows_from_features(features: list[dict]) -> list[tuple[str, str, str]]`,
  `write_table(rows, path: Path) -> None` (sorted by `huc8`, no duplicates), and
  `main(argv) -> int` with `--out` defaulting to the packaged CSV.

- [ ] **Step 1: Write the card**: delegate reads `tools/build_regional_skew_hucs.py` (for the
service constants and the request style only), `flowfreq/regional_skew.py:360-380`
(`_states_of`, the table reader, so the CSV shape matches), writes the three new files. Tests, all
offline on the fixture: `test_rows_from_features_keeps_non_us_codes` (a feature with
`states="CN,MT"` → `("09040001","St. Marys","CN,MT")`, Review Focus 3),
`test_write_table_sorted_unique`, `test_pages_stop_when_exceeded_transfer_limit_false`
(a fake session returning two pages). No network in tests.

`acceptance = "run.bat accept tests/test_build_huc8_states.py tests/test_regional_skew.py"`

- [ ] **Step 2: Queue, run, review** as in Task 2.

- [ ] **Step 3 (lead): Run the tool live and check the table**

Run: `python tools/build_huc8_states.py` then
`python -c "import pandas as pd; d=pd.read_csv('flowfreq/data/regional_skew_huc8_states.csv',dtype=str); print(len(d), d.huc8.is_unique)"`
Expected: about 2,300 rows (the WBD has ~2,270 HUC8s plus territories), `True`. Every one
of the current 417 rows is still present with the same `states` (diff the two files).

- [ ] **Step 4 (lead): Remove the "lists only some HUC8s" caveat** in
`flowfreq/regional_skew.py:418` and `:807` comments and the TODO/roadmap lines; run
`run.bat accept tests/test_regional_skew.py tests/test_regional_skew_at.py`; PR; CI; merge
(`Closes #205`).

---

### Task 5: Lead investigations that unlock the next cards (#206, #210, #207, #211)

Each sub-item ends either in a card file under `config/cards/` or a decision written into
the issue. None is delegated until its spec leaves one reasonable implementation.

- [ ] **5a. Catalog defect #206.** Run
`python -c "from flowfreq import waterdata; print(waterdata.request('monitoring-locations', params={'monitoring_location_number':'15388960'}).json()['features'][0]['properties'])"`
(with `USGS_API_KEY` exported) and read `state_code`, `district_code`, `country_code`.
If the API itself says DC (district 11 collects foreign-sited gages), the fix in
`tools/build_gage_catalog.py` is: a gage whose `country_code != "US"` gets `state` from its
coordinates' side of the border (`CA` → use `country_code` as the state column, `"CN"`, as
the HUC8 table does), never `DC`. Write card `ff-catalog-foreign-gages`: delegate edits
`tools/build_gage_catalog.py` plus `tests/test_catalog.py` with a synthetic feature; the
lead rebuilds only `--states DC AK` and commits the catalog. Then `tools/nationwide_smoke.py`
must say "no USGS gage in the catalog" for DC instead of failing (same card or a follow-up).

- [ ] **5b. StreamStats snap distance #210.** Read the StreamStats snap endpoint that
`flowfreq/streamstats.py:512 snap_point` calls and check whether it accepts a search
distance; if it does, card `ff-snap-distance` adds `snap_distance_m: Optional[float]` to
`snap_point` and `delineate_and_get_characteristics`, and `--snap-distance` to the smoke
script, with tests on a recorded response. If it does not, the fix is manual coordinates:
card `ff-smoke-coords` adds `--coords STATE=lat,lon` overrides to `tools/nationwide_smoke.py`
(tests in `tests/test_nationwide_smoke.py`), and the lead picks downstream points for AR,
GA, ND, PA, PR, TX from the gage pages. Rerun FL, MI, NE first (`--states FL MI NE
--resume`): their 404s may be gone.

- [ ] **5c. Region mismatches #207.** WY: the verified library carries WRIR 03-4107 regions
GC1239–GC1244, while NSS now locates GC2098–GC2115 (18 codes): run
`python tools/snapshot_nss.py WY` to see which report those codes cite. Decision for the
user: keep the 2003 equations and map the new codes to them (if the new codes are the same
regions re-cut), or transcribe the newer study. NC: compare `nss_snapshot/NC.json` region
codes with what `locate_regression_regions("NC", ...)` returns at 02087500's basin; a
mismatch inside NSS itself goes to the issue as a service finding. NV: call
`locate_regression_regions("NV", polygon)` for 10316500 and record whether NSS has no
peak-flow regions there (the WSP 2433-era NV library covers region 6 only; see
`docs/NEVADA_REGION6_SCHEMA.md`).

- [ ] **5d. Region-indicator parser #211.** Collect the 24 `unparsed` entries from
`nss_snapshot/{GA,NC,SC}.json`. The GA form seen is
`10^(2.81+0.00188*PCTREG1-0.00333*PCTREG3+0.00342*PCTREG5)*DRNAREA^(0.597+0.00150*PCTREG2)`:
affine expressions in characteristic codes inside the exponent of 10 and in a predictor's
exponent. If all 24 fit "affine exponent in characteristics", write card
`ff-nss-affine-exponents` (`flash-high`): extend the parser in
`flowfreq/regression/nss.py` to an exponent that is a sum of `coef*CODE` terms, carry them in
the equation schema (`flowfreq/regression/equations.py`), evaluate them in the evaluator,
with the 24 strings as fixture cases in `tests/test_nss_snapshot.py`; rebuild the three
snapshots with `tools/build_nss_snapshot_library.py --states GA NC SC` afterwards (lead).

---

### Task 6: Lead data items from the reports (#208, #213, #204) and the #209 decision packet

- [ ] **6a. SC GC1587 #208.** Open the SC report NSS cites for GC1587 (the snapshot's
`source_reports`), find the DRNAREA applicable range for that region, and correct `range_max`
in `flowfreq/data/regression/nss_snapshot/SC.json` with a `notes` line citing page and
table. Test: `tests/test_regression_nss_snapshot.py::test_sc_gc1587_drnarea_range` asserts
the corrected max. If the report range really is 1.0 mi², record that instead and leave the
refusal (then the smoke row needs a smaller basin).

- [ ] **6b. `sep_log` on LA/AR/MS #213.** Read the three reports' standard-error tables; if
ASEp is in percent, convert with the existing helper the snapshot builder uses for the other
1441 equations; if it is already log10, store it directly. Set `sep_log` on the 32 equations
in the three JSON files with `notes`, and add one assertion per state in
`tests/test_regression_nss_snapshot.py` (`sep_log` non-null for every equation).

- [ ] **6c. HUC membership rows #204.** For each of the 15 skew studies added in #195, add
rows to `flowfreq/data/regional_skew_hucs.csv` from the study's HUC listing (the row's
`notes` names the page). Add `test_new_studies_have_huc_rows` to `tests/test_regional_skew.py`, asserting each of
the 15 study ids has at least one row in `regional_skew_hucs.csv`. Expect several
sessions; do the states the smoke run touched first.

- [ ] **6d. Decision packet for #209 (user).** Add to the issue one table: state, candidate
sets (code, report, year), the smoke run's pick, and a recommendation (newest set unless
the newest is urban-only or a sub-region). For MT and NM state what each candidate's
applicable range is at the smoke point, so the user can pick. Nothing in code until the
user answers; then card `ff-region-defaults` adds `flowfreq/data/regression/region_defaults.json`
(state → default NSS code) consulted by `assemble` when `region is None` and the located
sets overlap, with provenance `region_choice_source="state default"`.

---

### Task 7: Hygiene (#218, dependabot, stale app checkout, v0.10.4)

- [ ] **7a. Live tests #218.** The failing test is
`tests/test_streamstats.py::TestLiveNSS::test_skookumchuck_low_flow_matches_its_equations`;
ss-hydro returned "Cannot find a 'wshprecip' raster for watershed 1710", so WA PRECIP was
unavailable server-side (it recovered on its own once before, 2026-10-13).
Run: `gh workflow run live.yml && gh run watch`. If it passes, comment on #218 and leave it
for the user to close. If PRECIP is still missing after two weekly runs, the change is a
`pytest.skip` **only** when `watershed.unavailable` names PRECIP with a raster error
(a service outage, which the workflow is not there to catch), never a widened assertion;
state in the PR that a skip, not an xfail, was added and why.

- [ ] **7b. Dependabot #219 (pypi-publish 1.14.3), #220 (setup-python 7.1.0), #221
(filelock <5).** Wait for green CI on each, `gh pr merge --merge` in order 219, 220, 221.
For #221 confirm `pip install -e ".[dev]"` still resolves on 3.11 (filelock 4.x drops
nothing flowfreq uses).

- [ ] **7c. Local app checkout.** `git -C D:/a/flowfreq-app pull --ff-only origin main`
(origin/main is at #7, pins v0.16.0; the local branch is at the v0.10.1 pin).

- [ ] **7d. v0.10.4 release (user).** `gh release create v0.10.4 --verify-tag --notes-file
<CHANGELOG 0.10.4 section>` attaching PyPI's dists; the session cannot run `gh release
create`. Hand the exact command and notes file path to the user.

---

### Task 8: Release 0.17.0 and the app pin bump

Depends on Tasks 2–5 merged (Task 6 items land as they finish; none blocks the release).

- [ ] **Step 1:** `CHANGELOG.md`: rename Unreleased to `0.17.0 - <date>`; `pyproject.toml`
version `0.17.0`; `tests/test_version.py` passes after `pip install -e ".[dev]"`.
- [ ] **Step 2:** PR `release/prep-0.17.0`; CI green on 3.11–3.14 and Fortran parity; merge.
- [ ] **Step 3:** `git tag -a v0.17.0 -m "flowfreq 0.17.0" <merge sha> && git push origin v0.17.0`;
`release.yml` publishes to PyPI and creates the GitHub release. **Do not** run `gh release create`.
- [ ] **Step 4:** Branch `release/0.17.x` from the tag, push; bump `main` to `0.18.0.dev0`;
update CLAUDE.md "Branches and Releases".
- [ ] **Step 5:** flowfreq-app PR: `requirements.txt` pin `@v0.17.0`, ungaged page passes a
list for the region selector when the estimate names more than one side (Task 3). Open the
PR; **the user merges and deploys** (memory: merge-workflow).

---

## Self-review notes

- Coverage: every open issue #203–#218 and every unchecked TODO rollout line has a task or
  is named as parked in the inventory. The 41 `pending` skew rows, #212 and methow-flow tiers
  stay parked by the user's earlier decision.
- Delegation boundary: the three cards written in full (Tasks 2–4) need no data, no network
  and no report; everything that does sits with the lead (Tasks 5–7) or the user (#209, app
  merge, releases the classifier blocks).
- The agentops contract is met except two deliberate deviations recorded in Global
  Constraints: no uv dependency in flowfreq's `pyproject.toml`, and no push block.
