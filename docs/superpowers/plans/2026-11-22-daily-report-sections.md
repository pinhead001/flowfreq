# Daily Report Sections Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Correct the HYSEP windows, then add `monthly_statistics`, `summary_hydrograph` and `baseflow` report sections.

**Architecture:** Task 1 is its own PR against `main` (and a port to `release/0.28.x`). Tasks 2 to 5 are a second PR: one module, `flowfreq/reporting/regime_sections.py`, registers three kinds over `daily.py`'s helpers, then the chat tools and docs.

**Tech Stack:** Python 3.11+, pandas, matplotlib, pytest; no new dependency.

**Spec:** `docs/superpowers/specs/2026-11-22-daily-report-sections-design.md`

## Global Constraints

- Must parse on Python 3.11. `make check` (black, isort, mypy, pytest) passes before each PR; `PYTHONSAFEPATH=1 PYTHONUTF8=1`.
- NumPy docstrings, type hints, no `print`, no bare `except`. Duration quantiles are lowercase q. Do not edit `vendor/`.
- No network in tests. A spec stores inputs, and a section rebuilt from its saved options gives the same blocks.
- Commit trailer `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Never push `release/*`; never tag.

## Review Focus

1. A stream with zero-flow days: figures use a linear axis and draw; a zero-flow year has BFI "-", not a crash (Tasks 2-4).
2. A record too short for one complete month or one complete year: `ReportError` naming what is missing (Tasks 2, 4).
3. Negative daily values (tidal): treated as missing in every table and figure, not as minima (Tasks 2-4).
4. `baseflow` with a HYSEP method and no drainage area known: `ReportError` that says to give `drainage_area_sqmi` (Task 4).
5. Fewer than five years for a date, and Feb 29: the row says how many years it rests on (Task 3).

---

### Task 1: HYSEP windows (branch `fix/hysep-windows` from `origin/main`)

**Files:** Modify `flowfreq/regime.py` (module docstring, `_hysep_interval_days`, three methods, `separate_baseflow`), `tests/test_regime.py` (`TestHysep`), `docs/vignette_lowflow_regime.md`, `CHANGELOG.md` (Fixed). Create `tests/fixtures/hysep/01472157_daily.csv` (`date,flow_cfs`, 1968-12-01 to 1970-02-28 and 1992-03-01 to 1992-05-31, from the Water Data API) and a `README.md` there (source, retrieval date, what the report prints, the 11-day offset).

**Interfaces:** Produces `_hysep_interval_days(drainage_area_sqmi: float) -> int` (2N*, odd, 3 to 11) and `_hysep_fixed`, `_hysep_sliding`, `_hysep_local_minimum`, each `(flows: np.ndarray, interval_days: int) -> np.ndarray`.

- [ ] **Step 1: failing tests** in `TestHysep` (existing `n_star=2` calls become `interval_days=3`, and the fixed step test checks `bf[0:3]`):

```python
@pytest.mark.parametrize(
    "area,days", [(0.5, 3), (32.0, 3), (59.1, 5), (100, 5), (1000, 7), (1e6, 11)]
)
def test_interval_days(area, days):
    assert _hysep_interval_days(area) == days
# replaces test_interval_days_increases_with_drainage_area and the "commonly cited" case

# French Creek, area 59.1, through separate_baseflow:
# fixed, on a frame starting 1992-04-05: bf["1992-04-05":"1992-04-09"] all == 49.0
# sliding: bf["1992-04-18"] == 42.0
# local minimum: bf == flow on April 9, 15, 21, 24 (1992)
# 1969: index moved back 11 days; monthly means of local-minimum baseflow within 0.3 of the report's twelve,
# annual within 0.5 percent of 27.00, and monthly streamflow within 0.005 of 42.87 42.18 57.45 46.60 28.55 19.30
# 94.13 30.39 21.70 22.48 43.10 50.48 (the offset is right)
```

- [ ] **Step 2:** `pytest tests/test_regime.py -k Hysep` fails on those.
- [ ] **Step 3:** implement: nearest of `(3, 5, 7, 9, 11)` to `2 * A**0.2`, a tie to the smaller; fixed width, and sliding and local-minimum centered windows, all `interval_days`. Rewrite the docstring's HYSEP paragraphs: definition from WRIR 96-4040, what it was checked against, the 0.4 percent left, and drop the "not checked" warning (and the vignette's mention of it).
- [ ] **Step 4:** `pytest tests/test_regime.py` passes, then `make check`.
- [ ] **Step 5:** CHANGELOG "Fixed": windows were about twice as wide and unbounded, so every `hysep_*` result changes (French Creek 1969 local minimum: 25.1 to 27.1 cfs, published 27.00). Commit, push, PR. Cherry-pick onto `fix/hysep-windows-0.28` from `origin/release/0.28.x`, PR into `release/0.28.x`.

### Task 2: `monthly_statistics` (branch `feat/daily-report-sections`, on Task 1 merged)

**Files:** Create `flowfreq/reporting/regime_sections.py`, `tests/test_reporting_regime_sections.py`. Modify `flowfreq/reporting/__init__.py` (import to register, docstring), `session.py` (`MAIN_TABLES`), `tests/test_reporting_build.py` (kinds set).

**Interfaces:** Consumes `daily._daily_frame`, `_record_table`, `_label`, `_STATION_OPTIONS`. Produces `MONTHLY_CAPTION: str = "Statistics of the monthly mean flows."`, `_build_monthly_statistics(given) -> SectionContent`, and `_valid_frame(owner: str, options: Dict[str, Any]) -> pd.DataFrame` (`_daily_frame` with negative values dropped; `ReportError` when none is left), which Tasks 3 and 4 use.

- [ ] **Step 1: failing tests.** Main table: 12 rows, first "October"; columns month, years, mean, median, standard deviation, maximum (water year), minimum (water year), percent of annual volume; each number equals one computed in the test from `monthly_flow_summary(FRAME)` complete rows through `cfs`; percents (mean of monthly means times average days in month) sum to 100 within rounding. A month with 20 missing days is not in its year count. A one-year month shows "-" for standard deviation. `annual_table=True` adds one row per water year with 12 month cells ("-" when incomplete). Title "Monthly flow statistics at ...". One Figure. Round trip: rebuilding from `content.options` gives equal tables. Refusals: 10-day record ("no complete month"); unknown option; `annual_table` not a bool.
- [ ] **Step 2:** run, fail. **Step 3:** implement. **Step 4:** the new file and `tests/test_reporting_build.py` pass. **Step 5:** commit.

### Task 3: `summary_hydrograph`

**Files:** Modify `flowfreq/hydrograph.py`, `flowfreq/reporting/regime_sections.py`, `session.py`, `tests/test_hydrograph.py`, the two report test files.

**Interfaces:** Produces `SUMMARY_HYDROGRAPH_CAPTION: str = "Daily mean flow by date, for the 1st and 15th of each month."` and `_build_summary_hydrograph(given) -> SectionContent`.

- [ ] **Step 1: failing tests.** `plot_summary_hydrograph` and `plot_flow_duration_curve` annotate water years: a record 1990-10-01 to 1992-09-30 says "(2 years". With `yscale="log"` and a 0.1 cfs minimum no y tick label is "0". Section: 24 rows, October 1 first; years and percentiles equal `Hydrograph.get_summary_stats` at that day; header names percentiles of non-exceedance; a zero-flow record still draws; a one-month record builds, and a date with no value in any year is a row with 0 years and "-" cells.
- [ ] **Steps 2-5:** fail, implement, pass, commit.

### Task 4: `baseflow`

**Files:** Modify `flowfreq/reporting/regime_sections.py`, `session.py`, the two report test files.

**Interfaces:** Produces `BASEFLOW_CAPTION: str = "Baseflow index and mean flows for the complete water years."`, `_METHOD_WORDS: Dict[str, Tuple[str, str]]` (method to sentence and citation, keys exactly `BASEFLOW_METHODS`), `_build_baseflow(given) -> SectionContent`.

- [ ] **Step 1: failing tests.** Main table: method, complete water years, mean streamflow, mean baseflow, baseflow index (equals `FlowRegime(...).summary()["baseflow_index"]` to 3 decimals), and both means in inches per year only with an area (`cfs * 365.25 * 86400 * 12 / (area * 5280**2)`). Annual table equals `baseflow_index(FRAME, method=...)`, incomplete years left out and counted in the text. Two Figures. `references` holds only the method's citation. Text says BFIs of different methods are not comparable. HYSEP without area, unknown method, `min_days` outside 1-366, no complete year, all-zero flow: `ReportError`, each message matched.
- [ ] **Steps 2-5:** fail, implement, pass, commit. Typical year: the complete year whose mean flow is nearest the median, the earliest on a tie.

### Task 5: chat tools and docs

**Files:** Modify `flowfreq/mcp_server.py`, `tests/test_mcp_server.py` (`TOOLS`, `TestDailyAndFutureTools`), `docs/CHAT_REPORTS.md`, `README.md`, `CLAUDE.md`, `CHANGELOG.md` (Added), `TODO.md`.

**Interfaces:** tools after `flow_duration`: `monthly_statistics(site_no, annual_table=False, start_date=None, end_date=None, title=None, position=None)`, `summary_hydrograph(site_no, start_date=None, end_date=None, title=None, position=None)`, `baseflow(site_no, method="ih_smoothed_minima", drainage_area_sqmi=None, annual_table=True, start_date=None, end_date=None, title=None, position=None)`, each returning `report.add_section(...)`.

- [ ] **Step 1: failing tests.** Each tool adds its kind with the stubbed download and prints its main table; a non-digit `site_no` is refused for all five daily tools; `TOOLS` has 27 in order.
- [ ] **Step 2-3:** fail, implement, with the server's instructions naming the three; "twenty-four" becomes "twenty-seven" in `CLAUDE.md`, `README.md` (not its v0.28.0 history row) and `docs/CHAT_REPORTS.md`.
- [ ] **Step 4:** `make check`; build one real report at USGS 03606500 with the three sections in DOCX and Markdown, and look at every figure and table.
- [ ] **Step 5:** commit, push, PR.
