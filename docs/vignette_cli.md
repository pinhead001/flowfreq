# Vignette: Command-Line Usage

HydroLib installs a `hydrolib` CLI entry point (via `pyproject.toml`) for running validation and benchmarks against known reference results.

## Setup

```bash
cd /path/to/hydrolib

# Activate your virtual environment
.venv\Scripts\activate         # Windows
# source .venv/bin/activate    # macOS/Linux

# Install (if not already done)
pip install -e ".[dev]"
```

After install, `hydrolib` is available on your PATH.

---

## Commands

### `hydrolib validate`

Runs the full numerical validation suite — comparing the native EMA implementation against fixtures extracted from the USGS `peakfqr` R package (the authoritative reference). Exits with code 1 if any benchmark fails.

```bash
hydrolib validate
```

Example output:
```
Running validation benchmarks...

Big Sandy River (03606500)
  mean_log    : 4.6841  expected 4.6841  ✓ PASS  (Δ = 0.000000)
  std_log     : 0.2103  expected 0.2103  ✓ PASS  (Δ = 0.000000)
  skew_station: 0.3152  expected 0.3152  ✓ PASS  (Δ = 0.000000)
  q_100yr     : 52432   expected 52418   ✓ PASS  (Δ = 14 cfs, 0.03%)

All 12 benchmarks passed.
```

### `hydrolib benchmark`

Runs benchmarks and prints a detailed report. Use `--format json` to pipe results into other tools.

```bash
# Human-readable text report
hydrolib benchmark

# JSON (for programmatic consumption)
hydrolib benchmark --format json

# Save JSON to file
hydrolib benchmark --format json > benchmark_results.json
```

Example JSON output structure:
```json
{
  "big_sandy_mean_log": {
    "passed": true,
    "value": 4.6841,
    "expected": 4.6841,
    "delta": 0.0,
    "tolerance": 0.001
  },
  ...
}
```

---

## Exit Codes

| Code | Meaning |
|------|---------|
| `0`  | All benchmarks passed |
| `1`  | One or more benchmarks failed (`validate` only) |

---

## Running Tests

The full pytest suite covers all public APIs including EMA, plotting, export, and the Streamlit runner:

```bash
# The whole suite. The default marker selection lives in pyproject.toml's
# addopts, so a bare pytest is already what CI runs.
pytest tests/

# Specific module
pytest tests/test_ffa_runner.py -v

# With coverage
pytest tests/ --cov=hydrolib --cov-report=term-missing
```

---

## Formatting and Linting

```bash
# Format all source files
black hydrolib/ app/ tests/
isort hydrolib/ app/ tests/

# Type checking
mypy hydrolib/ --ignore-missing-imports --strict
```
