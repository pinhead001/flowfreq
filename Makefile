# Developer entry points. CI calls these same targets, so what you run locally
# and what the build runs cannot drift apart.

PYTHON ?= python
PKGS := flowfreq/ tests/

# The marker deselection lives in pyproject.toml's addopts, not here, so a bare
# `pytest` is already correct and there is exactly one place to get it wrong.

PYTEST := $(PYTHON) -m pytest

# Both variables below are set with `export` rather than inline as
# `VAR=1 command`, which is shell syntax cmd.exe does not understand -- make
# puts them in the recipe's environment itself, so these targets work on
# Windows as well as Unix.

# UTF-8 regardless of the platform's locale. Eleven files under flowfreq/ and
# tests/ contain non-ASCII characters, and on a Windows console isort cannot
# encode them to cp1252: it reports "Unable to parse file" and *skips the
# file*, so `isort --check` passes without having checked it. Exported
# globally rather than per-target because this project is UTF-8 everywhere --
# that is already what Linux CI gets by default, so this makes Windows match
# rather than changing anything.
export PYTHONUTF8 := 1

# PYTHONSAFEPATH=1 stops Python prepending the working directory to sys.path,
# which is what the `pytest` console script does and `python -m pytest` does
# not. Without it a local run can pass while CI fails on imports. Scoped to the
# four targets that invoke pytest rather than exported globally, so it does not
# reach `fortran` or `golden`. On `parity` it does also cover
# build_fortran/build.py, which is safe: that script imports only the standard
# library and locates its sources through Path(__file__), never through the
# working directory.

.DEFAULT_GOAL := help
.PHONY: help check lint typecheck fmt test test-all cov clean-verify fortran parity golden clean

help:  ## Show this help
	@awk -F':.*?## ' '/^[a-z-]+:.*## /{printf "  %-10s %s\n", $$1, $$2}' Makefile

check: lint typecheck test  ## Everything CI checks, in CI's order

lint:  ## Formatting check (does not modify files)
	$(PYTHON) -m black --check --diff $(PKGS)
	$(PYTHON) -m isort --check-only --diff $(PKGS)

typecheck:  ## mypy over the modules that pass today (see pyproject overrides)
	$(PYTHON) -m mypy flowfreq/

fmt:  ## Apply formatting
	$(PYTHON) -m black $(PKGS)
	$(PYTHON) -m isort $(PKGS)

test: export PYTHONSAFEPATH := 1
test:  ## Run the suite as CI does
	$(PYTEST) tests/

cov: export PYTHONSAFEPATH := 1
cov:  ## Test suite with a coverage report
	$(PYTEST) tests/ --cov=flowfreq --cov-report=term-missing --cov-report=xml

# Four times during the repo split a check passed locally for a reason that did
# not hold in a clean checkout: a stale build/lib/ shipping two packages in one
# wheel, a working directory shadowing an installed package, and a leftover
# __pycache__ directory steering isort's first-party classification. Each looked
# like a passing check. This target removes everything that can do that, so a
# green result means the tree is green rather than the leftovers are helpful.
# Note it depends on `clean`, which removes the f2py extension: this reports
# the ~476 tests a fresh checkout runs, not the ~608 available once the
# extension is built. That is the point -- it answers "is the tree green", not
# "is my machine green". Run `make parity` afterwards to get the Fortran
# comparisons back.
clean-verify: clean  ## Wipe every build artifact, then run the full gate
	rm -rf build/ dist/ *.egg-info .mypy_cache .pytest_cache
	find . -path ./vendor -prune -o -name "__pycache__" -type d -print0 2>/dev/null | xargs -0 rm -rf
	$(MAKE) check

test-all: export PYTHONSAFEPATH := 1
test-all:  ## Run everything, including the network tests
	$(PYTEST) tests/ -m ""

fortran:  ## Build the f2py extension from vendor/peakfqr (needs gfortran + meson)
	$(PYTHON) build_fortran/build.py

# The import check is not redundant. test_live_vs_golden.py calls importorskip,
# so a failed build would skip every parity test and still exit 0 -- which is
# exactly the silent pass this target exists to prevent.
parity: export PYTHONSAFEPATH := 1
parity:  ## Build the extension and check the golden files against it (needs gfortran + meson)
	$(PYTHON) build_fortran/build.py
	$(PYTHON) -c "from flowfreq.peakfqr import emafitpr"
	$(PYTEST) tests/fortran_parity/

golden:  ## Regenerate Fortran parity golden files (needs the extension)
	$(PYTHON) tools/gen_fortran_golden.py

clean:  ## Remove build and test artifacts
	rm -rf build_fortran/mbuild build_fortran/native.ini build_fortran/_emafort*.so
	rm -rf flowfreq/peakfqr/_emafort*.so .pytest_cache .mypy_cache
