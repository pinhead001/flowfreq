# HydroLib Hybrid 17C Implementation — Agent Build Instructions

> **For the Agent:** Read this entire document before writing a single line of code.  
> Start by generating the TODO list in Step 0. Then work through every task in order.  
> Prompt the user only when you genuinely cannot proceed without their input.  
> Do as much work as possible autonomously.

---

## Context

**Project:** HydroLib — a Python library for hydrologic frequency analysis  
**Goal:** Implement a hybrid Bulletin 17C approach that wraps PeakfqSA (USGS Fortran reference implementation) alongside the existing native EMA implementation  
**Reference Software:** PeakfqSA by Tim Cohn (USGS) — ~45,000 lines of Fortran-95  
**Test Case:** Big Sandy River at Bruceton, TN (USGS gage 03606500)

**The hybrid approach means:**
- Native Python EMA for speed and integration
- PeakfqSA subprocess wrapper for validation and edge cases
- Side-by-side comparison tools to measure agreement
- Unified `FrequencyAnalyzer` API that supports both modes

---

## Reference Repository: peakfqr

**Path:** `C:\a\hal\_shared\peakfqr`  
**Source:** USGS R package implementing Bulletin 17C / PeakfqSA methods

### How to use peakfqr as a reference

The agent must read this repository **before writing any computation code**. It is the authoritative guide to how the underlying Fortran is called, parameterized, and interpreted.

**Read `\src\` for computation logic:**
- These are the R source files that implement EMA, MGBT, confidence intervals, and frequency curve fitting
- Treat each `.R` file in `\src\` as the specification for the equivalent Python module in HydroLib
- Mirror the function signatures, parameter names, and logic flow wherever practical
- Pay close attention to how R calls into Fortran — the same call patterns apply to the Python subprocess wrapper

**Read all other folders to understand context:**

| Folder | What to learn from it |
|--------|----------------------|
| `\R\` | Public-facing API — shows intended user workflow, input/output contracts, argument names and defaults |
| `\tests\` or `\test\` | Expected inputs and outputs — use these as additional validation fixtures |
| `\man\` | Function documentation — parameter descriptions, units, valid ranges, edge case notes |
| `\vignettes\` | End-to-end worked examples — shows how inputs flow through to frequency curve output |
| `\data\` or `\inst\` | Example datasets — extract any with gage data usable as additional test fixtures |
| `DESCRIPTION` | Package dependencies — any referenced USGS packages may contain additional Fortran or methods |
| `NAMESPACE` | Exported functions — identifies which functions are the primary public API |

### Specific things to extract from peakfqr before coding

Before writing a single HydroLib module, read peakfqr and document the following in `TODO.md`:

1. **Fortran call signatures** — what arguments does R pass to the Fortran routines in `\src\`? What are the argument order, types, and units?
2. **EMA parameter conventions** — how does peakfqr define and pass perception thresholds, historical period, and systematic record length?
3. **MGBT implementation** — how does peakfqr call and interpret the Multiple Grubbs-Beck Test results?
4. **Confidence interval method** — does peakfqr use the same CI approach as the PeakfqSA manual? Note any differences.
5. **Regional skew weighting** — exact formula and MSE combination method used
6. **Output format** — what fields does peakfqr return? Map each to the equivalent `PeakfqSAResult` field in HydroLib
7. **Edge cases handled** — zero flows, all-censored records, single threshold, ties — note how peakfqr handles each

Add a section `## peakfqr Reference Notes` to `TODO.md` with findings from steps 1-7 before proceeding to Step 3.

---

## Step 0 — Read peakfqr, Then Generate the TODO List

**Before doing anything else**, read the USGS peakfqr reference repository, then generate the TODO list.

### Step 0a — Read peakfqr

The peakfqr repository is at: `C:\a\hal\_shared\peakfqr`

Read in this order:

```
1. C:\a\hal\_shared\peakfqr\DESCRIPTION        — package overview, dependencies
2. C:\a\hal\_shared\peakfqr\NAMESPACE           — exported (public) functions
3. C:\a\hal\_shared\peakfqr\src\               — ALL files (Fortran call layer — primary reference)
4. C:\a\hal\_shared\peakfqr\R\                 — ALL files (R API layer — input/output contracts)
5. C:\a\hal\_shared\peakfqr\man\               — function docs (parameter names, units, valid ranges)
6. C:\a\hal\_shared\peakfqr\vignettes\         — worked examples (end-to-end workflows)
7. C:\a\hal\_shared\peakfqr\tests\ or \test\   — test fixtures (extract as additional HydroLib fixtures)
8. C:\a\hal\_shared\peakfqr\data\ or \inst\    — example datasets
```

While reading, extract and record in `TODO.md` under a `## peakfqr Reference Notes` section:

- Fortran routine names in `\src\` and their argument lists (name, type, order, units)
- How R passes perception thresholds, historical period length, and systematic record length to Fortran
- MGBT call signature and how results are interpreted
- Confidence interval method — does it match PeakfqSA manual? Note differences
- Regional skew weighting formula and MSE combination
- Complete output field list — map each to `PeakfqSAResult` fields in HydroLib
- Edge cases explicitly handled (zero flows, all-censored, ties, single threshold)

**Do not skip this step.** The peakfqr `\src\` files are the specification for every computation module written below.

### Step 0b — Generate TODO.md

After reading peakfqr, scan the existing HydroLib codebase and generate a comprehensive TODO list in `TODO.md`:

```
- [ ] Every file to be created (with purpose)
- [ ] Every existing file to be modified (with what changes)
- [ ] Every test to be written
- [ ] Every dependency to be installed
- [ ] Every configuration change needed
- [ ] Open questions requiring user input (mark with ❓)
```

Format as nested Markdown checkboxes grouped by phase. This is the live tracking document — check items off as they complete.

**Do not proceed past Step 0 until `TODO.md` exists and the `## peakfqr Reference Notes` section is populated.**

---

## Step 1 — Gather Information  *(Ask the user these questions once, up front)*

Before building anything, ask the user the following questions in a **single grouped message**. Do not ask questions one at a time across multiple turns.

```
Questions to ask:

1. Where is your HydroLib project root?
   (e.g., ~/projects/hydrolib)

2. Do you have PeakfqSA installed?
   - Yes — provide the path to the executable
   - No — should the agent attempt to download and compile from source?
   - No — I want wrapper stubs only (skips validation tests)

3. What Python version and virtual environment manager are you using?
   (e.g., Python 3.11, conda / venv / uv)

4. What is your regional skew and MSE for initial testing?
   (Default will be -0.302 / 0.3025 per Bulletin 17C national map)

5. Do you want the agent to push commits to Git as it works?
   - Yes — provide branch name (e.g., feature/hybrid-17c)
   - No — I will handle Git manually
```

Store all answers in memory. Use them throughout the build. Do not ask again.

---

## Step 2 — Environment Setup

**Do autonomously. Do not prompt the user.**

### 2a — Verify Project Structure

```bash
# Confirm HydroLib is structured as expected
find . -name "*.py" | head -40
cat pyproject.toml || cat setup.py || cat setup.cfg
```

If the structure differs significantly from the expected layout below, **stop and describe the difference to the user, then ask how to proceed:**

```
Expected layout:
src/hydrolib/
├── analysis/
├── data/
├── visualization/
└── __init__.py
tests/
pyproject.toml or setup.py
```

### 2b — Install Dependencies

```bash
pip install --break-system-packages pytest pytest-cov pytest-mock \
    black isort mypy ruff pandas numpy scipy click rich
```

Verify each installs cleanly. If any fail, attempt alternatives and document the issue in TODO.md.

### 2c — Confirm Test Suite Passes Baseline

```bash
pytest tests/ -q --tb=short
```

Record the baseline pass/fail count in TODO.md. If tests were already failing, note which ones and continue — do not attempt to fix pre-existing failures.

---

## Step 3 — Create Directory Structure

**Do autonomously.**

Create every directory and stub file listed below. Each stub must contain:
- Module-level docstring explaining its purpose
- All planned class/function signatures with `...` bodies
- Full type hints
- `TODO:` comments on each stub body

```
src/hydrolib/peakfqsa/
├── __init__.py
├── config.py         # PeakfqSAConfig, executable detection
├── wrapper.py        # Subprocess execution of PeakfqSA
├── io_converters.py  # .psf and .dat file generation
├── parsers.py        # .out file parsing → Python objects
└── validators.py     # Tolerance-based result comparison

src/hydrolib/validation/
├── __init__.py
├── benchmarks.py     # Bulletin 17C Appendix 10 test cases
├── comparisons.py    # Native vs PeakfqSA comparison engine
└── reports.py        # HTML/text/JSON report generation

tests/peakfqsa/
├── __init__.py
├── test_config.py
├── test_wrapper.py
├── test_io_converters.py
├── test_parsers.py
└── fixtures/
    ├── __init__.py
    ├── big_sandy.py          # Test data and expected results
    └── expected/
        ├── big_sandy.out     # Copy of expected PeakfqSA output
        └── big_sandy.csv     # Expected quantiles as CSV

tests/validation/
├── __init__.py
└── test_benchmarks.py

tests/integration/
├── __init__.py
└── test_hybrid_workflow.py
```

After creating the structure, verify it with `find src/hydrolib tests -name "*.py" | sort` and print the result.

---

## Step 4 — Big Sandy River Test Fixtures

**Do autonomously.** This data comes directly from the PeakfqSA user manual.

Also check `C:\a\hal\_shared\peakfqr\tests\`, `\data\`, and `\inst\` for any additional gage datasets or expected-output files. If found, create matching fixture files alongside `big_sandy.py` so they can be used as additional validation cases later.

Create `tests/peakfqsa/fixtures/big_sandy.py` with the following data exactly as specified:

```python
"""
Big Sandy River at Bruceton, TN (USGS gage 03606500)
Test fixture data sourced from PeakfqSA User Manual (Tim Cohn, USGS, 2012).
Used as the primary validation case for hybrid Bulletin 17C implementation.
"""

# Systematic annual peaks 1930-1973 (cfs)
SYSTEMATIC_PEAKS = {
    1930: 9100,  1931: 2060,  1932: 7820,  1933: 3220,
    1934: 5580,  1935: 17000, 1936: 6740,  1937: 13800,
    1938: 4270,  1939: 5940,  1940: 1680,  1941: 1200,
    1942: 10100, 1943: 3780,  1944: 5340,  1945: 5630,
    1946: 12000, 1947: 3980,  1948: 6130,  1949: 4740,
    1950: 9880,  1951: 5230,  1952: 4260,  1953: 5000,
    1954: 3320,  1955: 5480,  1956: 11800, 1957: 5150,
    1958: 3350,  1959: 2400,  1960: 1460,  1961: 3770,
    1962: 7480,  1963: 2740,  1964: 3100,  1965: 7180,
    1966: 1920,  1967: 9060,  1968: 3080,  1969: 2800,
    1970: 4330,  1971: 5080,  1972: 12000, 1973: 7640,
}

# Historical floods (known to exceed 18,000 cfs threshold)
HISTORICAL_PEAKS = {
    1897: 25000,
    1919: 21000,
    1927: 18500,
}

# Perception thresholds
THRESHOLDS = [
    {"start": 1890, "end": 1929, "lower": 18000.0, "upper": 1e50},
    {"start": 1930, "end": 1973, "lower": 0.0,     "upper": 1e50},
]

BEGYEAR = 1890
ENDYEAR = 1973
REGIONAL_SKEW = -0.5
REGIONAL_SKEW_SD = 0.55
STATION_NAME = "BIG SANDY RIVER AT BRUCETON, TN, 1890-1973"

# Expected results from PeakfqSA manual (page 26-27)
EXPECTED_PARAMETERS = {
    "mean_log":    3.717272,
    "std_log":     0.289200,
    "skew_weighted": -0.118702,
}

EXPECTED_QUANTILES = {
    # AEP: discharge (cfs) — from PeakfqSA manual output
    0.9950:  871.25,
    0.9900:  1045.59,
    0.9500:  1706.18,
    0.9000:  2203.77,
    0.8000:  2990.15,
    0.6667:  3957.50,
    0.5000:  5284.36,
    0.2000:  9166.15,
    0.1000:  12134.65,
    0.0400:  16276.60,
    0.0200:  19617.73,
    0.0100:  23158.65,
    0.0050:  26912.12,
    0.0020:  32217.14,
}

EXPECTED_CONFIDENCE_INTERVALS = {
    # AEP: (lower_95, upper_95)
    0.1000: (9766.00,  15218.32),
    0.0200: (15154.99, 29124.18),
    0.0100: (17388.03, 37986.08),
}

TOLERANCE_PERCENT = 1.0  # Results must match PeakfqSA within 1%
```

---

## Step 5 — Configuration Module

**Do autonomously.** Implement `src/hydrolib/peakfqsa/config.py`:

Requirements:
- `PeakfqSAConfig` dataclass with fields: `executable_path`, `timeout_seconds`, `temp_dir`, `keep_temp_files`
- `find_peakfqsa()` function — search in order: env var `PEAKFQSA_PATH`, user-supplied path, common system paths (`/usr/local/bin`, `~/.local/bin`, `~/tools/peakfqsa/`)
- `validate_peakfqsa(path)` — runs `PeakfqSA` with no args, checks it starts without crashing, returns version string
- Graceful `PeakfqSANotFoundError` with a message that tells the user exactly how to fix it
- Full type hints, NumPy-format docstrings, logging

Write tests in `tests/peakfqsa/test_config.py`:
- Test executable detection with mocked filesystem
- Test version extraction from mock output
- Test `PeakfqSANotFoundError` message quality
- Test env var override

Run tests. Fix until they pass.

---

## Step 6 — I/O Converters

**Do autonomously.** Implement `src/hydrolib/peakfqsa/io_converters.py`.

**Before writing any code**, re-read the relevant files in `C:\a\hal\_shared\peakfqr\src\` to confirm the exact argument names, field order, and formatting that the Fortran routines expect. The `.psf` and `.dat` formats must match what the Fortran layer consumes — use peakfqr's R-to-Fortran calls as the ground truth.

### SpecificationFile class

Must produce output that exactly matches PeakfqSA `.psf` format from the manual:

```
# Generated by HydroLib
I {data_filename}
STATION {station_name}
BEGYEAR {begyear}
ENDYEAR {endyear}
GAGEBASE {gagebase}
SKEWOPT {skew_option}
GENSKEW {regional_skew}
SKEWSD {skew_sd}
LOMETHOD {lo_method}
CSV YES
PCPT_THRESH {start} {end} {lower} {upper}
```

### DataFile class

Must produce output matching PeakfqSA `.dat` format:
- `Q {year} {discharge}` for point observations
- `QINT {year} {lower} {upper}` for interval observations
- Comment header with station name

### Required methods on both classes:
- `from_analysis_params(...)` — class method
- `to_string()` — returns file content as string
- `write(path: Path)` — writes to file
- `validate()` — checks internal consistency, raises `ValueError` with clear message

Write the Big Sandy `.psf` and `.dat` as expected output strings in the fixture. Compare generated output character-by-character (allowing whitespace normalization). Tests must confirm the generated files exactly reproduce the Big Sandy example from the manual.

Run tests. Fix until they pass.

---

## Step 7 — PeakfqSA Wrapper

**Do autonomously.** Implement `src/hydrolib/peakfqsa/wrapper.py`.

**Before writing any code**, read `C:\a\hal\_shared\peakfqr\src\` to identify how the R package invokes the Fortran executable — argument order, environment variables, working directory assumptions, and temp file conventions. Mirror the same invocation pattern in the Python subprocess wrapper.

```python
class PeakfqSAWrapper:
    """
    Subprocess wrapper for the PeakfqSA Fortran executable.
    Handles file I/O, execution, and result parsing.
    """

    def __init__(self, config: PeakfqSAConfig) -> None: ...

    def run(
        self,
        peaks: dict[int, float],
        historical: dict[int, float],
        thresholds: list[dict],
        begyear: int,
        endyear: int,
        regional_skew: float,
        regional_skew_sd: float,
        station_name: str = "",
        **kwargs,
    ) -> "PeakfqSAResult": ...

    def _write_input_files(self, tmp_dir: Path, ...) -> tuple[Path, Path]: ...
    def _execute(self, spec_file: Path, tmp_dir: Path) -> str: ...
    def _parse_output(self, output_file: Path) -> "PeakfqSAResult": ...
    def is_available(self) -> bool: ...
```

Error handling requirements:
- `PeakfqSANotFoundError` — executable missing
- `PeakfqSAExecutionError` — non-zero exit code, include stdout/stderr in message
- `PeakfqSATimeoutError` — process exceeded timeout
- `PeakfqSAParseError` — output file could not be parsed

Tests in `test_wrapper.py`:
- Mock `subprocess.run` to test without PeakfqSA installed
- Confirm input files are generated correctly (inspect temp directory)
- Confirm error conditions raise correct exceptions
- Add `@pytest.mark.requires_peakfqsa` marker for tests that need the real binary

Register the `requires_peakfqsa` marker in `pyproject.toml` or `conftest.py`. Tests with this marker must be **skipped** (not failed) when PeakfqSA is not installed.

Run tests. Fix until they pass.

---

## Step 8 — Output Parser

**Do autonomously.** Implement `src/hydrolib/peakfqsa/parsers.py`.

**Before writing any code**, read `C:\a\hal\_shared\peakfqr\R\` and `\man\` to get the complete list of fields returned by the R package from the Fortran output. Use the peakfqr field names and structure to ensure `PeakfqSAResult` captures every value peakfqr exposes — nothing should be silently dropped.

Parse the PeakfqSA `.out` file format. The output has these sections (from manual page 26-31):

1. **Header block** — station name, spec file, data file, period summary
2. **EMA frequency curve parameters** — M, S, G (with and without regional skew)
3. **EMA frequency estimates** — AEP, quantile, CI-low, CI-high
4. **Low-outlier summary** — count, threshold, critical value
5. **Perception thresholds** — table of T_lower, T_upper by year range
6. **Sorted observations** — year, plotting position, observed Q, fitted Q

```python
@dataclass
class PeakfqSAResult:
    station_name: str
    begyear: int
    endyear: int
    n_peaks: int
    n_systematic: int
    n_historical: int
    low_outlier_count: int
    low_outlier_threshold: float
    parameters: dict[str, float]        # mean_log, std_log, skew_at_site, skew_weighted
    quantiles: dict[float, float]       # AEP → discharge
    confidence_intervals: dict[float, tuple[float, float]]  # AEP → (lower, upper)
    plotting_positions: pd.DataFrame    # year, plot_pos, q_obs, q_fit
    raw_output: str                     # full text for debugging
```

Parser implementation notes:
- Use `re.compile` patterns defined at module level (not inside functions)
- Parse defensively — log a warning and continue if a section is missing
- Never raise on missing optional sections; raise only on missing required sections

Test with the exact output text from the manual (embed it as a fixture string). Confirm every field of `PeakfqSAResult` matches expected values. Include a test that verifies `raw_output` is stored.

Run tests. Fix until they pass.

---

## Step 9 — Comparison Engine

**Do autonomously.** Implement `src/hydrolib/validation/comparisons.py`:

```python
@dataclass
class ComparisonResult:
    passed: bool
    tolerance_pct: float
    parameter_diffs: dict[str, float]     # param → % difference
    quantile_diffs: dict[float, float]    # AEP → % difference
    ci_diffs: dict[float, float]          # AEP → % difference on bounds
    max_diff_pct: float
    summary: str                          # human-readable one-liner

class FrequencyComparator:
    def __init__(self, tolerance_pct: float = 1.0) -> None: ...

    def compare(
        self,
        native: dict,           # output from HydroLib native analysis
        reference: PeakfqSAResult,
    ) -> ComparisonResult: ...

    def compare_parameters(self, native: dict, ref: PeakfqSAResult) -> dict[str, float]: ...
    def compare_quantiles(self, native: dict, ref: PeakfqSAResult) -> dict[float, float]: ...
    def compare_ci(self, native: dict, ref: PeakfqSAResult) -> dict[float, float]: ...
```

Default tolerances:
- Parameters (M, S, G): 0.5%
- Quantiles: 1.0%
- Confidence intervals: 2.0%

Write tests that:
- Confirm two identical results produce `passed=True` with all diffs = 0
- Confirm results just outside tolerance produce `passed=False`
- Confirm `summary` is a meaningful one-liner (not empty)

Run tests. Fix until they pass.

---

## Step 10 — Update FrequencyAnalyzer API

**Prompt the user before making changes here:**

```
I'm about to modify your existing FrequencyAnalyzer class to add hybrid mode support.
Before I do, I need to understand the current API:

1. Can you show me the current FrequencyAnalyzer signature and main analyze() method?
   (Or I can read it from the file — just confirm the path)

2. Are there any planned breaking changes I should be aware of?

3. Should I maintain 100% backward compatibility, or is it okay to add
   required parameters if they have good defaults?

I will not modify any existing method signatures without your approval.
```

After receiving the user's answer, add the following without breaking existing behavior:

```python
# New parameters on __init__ (all optional with defaults):
use_peakfqsa: bool = False
validate_with_peakfqsa: bool = False
peakfqsa_config: PeakfqSAConfig | None = None

# New method:
def validate(self, data, **params) -> tuple[AnalysisResult, ComparisonResult]:
    """Run both native and PeakfqSA, return both results and comparison."""
    ...

# Modified analyze() — unchanged behavior when use_peakfqsa=False:
def analyze(self, data, **params) -> AnalysisResult | tuple[AnalysisResult, ComparisonResult]:
    ...
```

Write a backward-compatibility test:
```python
def test_existing_api_unchanged():
    """Confirm existing callers need zero changes."""
    analyzer = FrequencyAnalyzer()          # no new args
    result = analyzer.analyze(sample_data)  # same return type as before
    assert isinstance(result, AnalysisResult)
```

Run all tests. Fix until they pass.

---

## Step 11 — Integration Test: Big Sandy End-to-End

**Do autonomously.**

In `tests/integration/test_hybrid_workflow.py`:

```python
@pytest.mark.requires_peakfqsa
def test_big_sandy_end_to_end():
    """
    Full integration test: run Big Sandy through native + PeakfqSA,
    compare, assert agreement within 1%.
    """
    from tests.peakfqsa.fixtures.big_sandy import (
        SYSTEMATIC_PEAKS, HISTORICAL_PEAKS, THRESHOLDS,
        BEGYEAR, ENDYEAR, REGIONAL_SKEW, REGIONAL_SKEW_SD,
        EXPECTED_QUANTILES, TOLERANCE_PERCENT,
    )

    analyzer = FrequencyAnalyzer(validate_with_peakfqsa=True)
    result, comparison = analyzer.validate(
        peaks=SYSTEMATIC_PEAKS,
        historical=HISTORICAL_PEAKS,
        thresholds=THRESHOLDS,
        begyear=BEGYEAR,
        endyear=ENDYEAR,
        regional_skew=REGIONAL_SKEW,
        regional_skew_sd=REGIONAL_SKEW_SD,
    )

    assert comparison.passed, (
        f"Native vs PeakfqSA disagreement exceeds {TOLERANCE_PERCENT}%:\n"
        f"{comparison.summary}\n"
        f"Max diff: {comparison.max_diff_pct:.2f}%"
    )
```

Also write a test that runs **without** PeakfqSA (native only) and verifies the quantiles against the hardcoded expected values from the fixture. This test must **not** require PeakfqSA.

---

## Step 12 — Benchmark Module

**Do autonomously.** Implement `src/hydrolib/validation/benchmarks.py`:

```python
"""
Bulletin 17C Appendix 10 benchmark test cases.
Each benchmark is a named, self-contained test scenario with
known inputs and expected outputs from the official USGS examples.
"""

class Benchmark:
    name: str
    description: str
    peaks: dict[int, float]
    historical: dict[int, float]
    thresholds: list[dict]
    expected_quantiles: dict[float, float]
    tolerance_pct: float

    def run_native(self) -> AnalysisResult: ...
    def run_peakfqsa(self, config: PeakfqSAConfig) -> PeakfqSAResult: ...
    def validate(self, config: PeakfqSAConfig | None = None) -> ComparisonResult: ...

# Register benchmarks
BENCHMARKS: dict[str, Benchmark] = {
    "big_sandy": BigSandyBenchmark(),
    # Additional benchmarks added here as they are validated
}

def run_all_benchmarks(peakfqsa_config=None) -> dict[str, ComparisonResult]: ...
def print_benchmark_report(results: dict[str, ComparisonResult]) -> None: ...
```

---

## Step 13 — CLI Validation Command

**Do autonomously.** Add a CLI command:

```bash
# Usage after install:
hydrolib validate --gage 03606500 --method both
hydrolib benchmark --list
hydrolib benchmark --run big_sandy
```

Implement using `click`. The command should:
- Accept a USGS gage ID or a path to a CSV of annual peaks
- Run native analysis, optionally PeakfqSA, print comparison
- Output: colored terminal table using `rich`
- Exit code 0 = passed, 1 = failed, 2 = error

Register the command in `pyproject.toml` under `[project.scripts]`.

---

## Step 14 — Documentation

**Do autonomously.**

### 14a — Docstrings

Every public class and function must have a NumPy-format docstring with:
- One-line summary
- Parameters section with types
- Returns section
- Raises section
- Examples section (at least one runnable example)

Run `python -m pytest --doctest-modules src/hydrolib/peakfqsa/` to verify doctest examples work.

### 14b — README Section

Add a "Hybrid Bulletin 17C Validation" section to the project README:

```markdown
## Hybrid Bulletin 17C Validation

HydroLib supports two modes for Bulletin 17C flood frequency analysis:

**Native mode** (default — no external dependencies):
```python
from hydrolib.analysis import FrequencyAnalyzer
analyzer = FrequencyAnalyzer()
result = analyzer.analyze(peaks, begyear=1930, endyear=2020)
```

**Validation mode** (requires PeakfqSA — compares native vs USGS reference):
```python
analyzer = FrequencyAnalyzer(validate_with_peakfqsa=True)
result, comparison = analyzer.validate(peaks, ...)
print(comparison.summary)
```

Install PeakfqSA: [instructions link]
```

### 14c — CHANGELOG

Add an entry to `CHANGELOG.md`:
```markdown
## [Unreleased]
### Added
- PeakfqSA integration for Bulletin 17C validation (`hydrolib.peakfqsa`)
- `FrequencyAnalyzer.validate()` method for side-by-side comparison
- Benchmark framework with Big Sandy River test case
- `hydrolib validate` and `hydrolib benchmark` CLI commands
- Comparison reports with per-quantile tolerance checking
```

---

## Step 15 — Final Quality Check

**Do autonomously.** Run all of the following and fix any issues before declaring done:

```bash
# Formatting
black src/ tests/
isort src/ tests/

# Type checking
mypy src/hydrolib/ --ignore-missing-imports --strict

# Tests (without PeakfqSA)
pytest tests/ -v -m "not requires_peakfqsa" --cov=hydrolib \
    --cov-report=term-missing --cov-fail-under=80

# Doctest
python -m pytest --doctest-modules src/hydrolib/peakfqsa/ src/hydrolib/validation/

# Check for leftover TODO stubs (warn but don't fail)
grep -rn "TODO:" src/hydrolib/ | grep -v ".pyc"
```

After the run, print a summary:
```
===== FINAL QUALITY REPORT =====
Tests:        X passed, Y skipped (requires_peakfqsa), 0 failed
Coverage:     XX%
Type errors:  0
Format:       Clean
TODOs remain: N (list them)
================================
```

---

## Step 16 — Update TODO.md

**Do autonomously.**

Go back to `TODO.md` (created in Step 0). Check off every completed item. Add any newly discovered issues as unchecked items.

At the top of `TODO.md`, add a status block:
```markdown
## Status
Last updated: [timestamp]
Phases complete: 1 / 1
Tests passing: X / Y
Coverage: XX%
PeakfqSA available: Yes / No (stubs only)
Open items: N
```

---

## Agent Behavior Rules

These rules apply throughout the entire build:

### When to Ask the User

**DO ask when:**
- Modifying existing public APIs (Step 10)
- PeakfqSA binary path is unknown and auto-detection fails
- A test reveals an unexpected behavior in the existing HydroLib implementation
- The existing codebase structure differs from what was expected
- A design decision has significant tradeoffs (document the options, ask once)

**DO NOT ask when:**
- Choosing between two implementation approaches that are both valid
- Deciding what to name a private method or variable
- Adding log statements or comments
- Choosing error message wording
- Formatting code
- Writing tests for code you just wrote

### Code Quality Minimums (Non-negotiable)

Every file produced must meet these standards:
- Type hints on **all** function signatures
- NumPy docstrings on **all** public classes and functions
- No `pass` without a `# TODO:` comment explaining what's missing
- No bare `except:` — always catch specific exception types
- No `print()` in library code — use `logging.getLogger(__name__)`
- Tests for every public function (at minimum: happy path + one error case)
- `black` and `isort` compliant

### Git Commits (if user enabled this in Step 1)

Commit after each step completes using conventional commits:
```bash
git add -A
git commit -m "feat(peakfqsa): add config and executable detection [Step 5]"
```

Commit messages must reference the step number.

### When Tests Fail

If a test fails:
1. Read the full traceback
2. Fix the code (not the test, unless the test is genuinely wrong)
3. Re-run the specific failing test
4. Re-run the full suite
5. If you cannot fix it within 3 attempts, **stop, describe the problem clearly, and ask the user**

Do not comment out failing tests or mark them `xfail` without explicit user approval.

---

## Quick Reference: Key File Paths

| File | Purpose |
|------|---------|
| `TODO.md` | Live task tracker — created in Step 0 |
| `src/hydrolib/peakfqsa/config.py` | PeakfqSA detection and config |
| `src/hydrolib/peakfqsa/wrapper.py` | Subprocess execution |
| `src/hydrolib/peakfqsa/io_converters.py` | .psf / .dat file generation |
| `src/hydrolib/peakfqsa/parsers.py` | .out file parsing |
| `src/hydrolib/validation/comparisons.py` | Native vs PeakfqSA comparison |
| `src/hydrolib/validation/benchmarks.py` | Bulletin 17C test cases |
| `tests/peakfqsa/fixtures/big_sandy.py` | Primary test data |
| `tests/integration/test_hybrid_workflow.py` | End-to-end integration test |

---

## Begin

**Generate `TODO.md` now. Then ask the user the questions in Step 1. Then build.**
