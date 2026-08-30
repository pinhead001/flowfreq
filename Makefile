# Developer entry points. CI calls these same targets, so what you run locally
# and what the build runs cannot drift apart.

PYTHON ?= python
PKGS := hydrolib/ tests/

# Tests needing something the environment may not have: the PeakfqSA binary
# (which does not exist -- see CLAUDE.md), network access, or the peakfqr
# reference data. Kept in one place because getting this list wrong by hand is
# how a local "green" run stops meaning anything.
DESELECT := not requires_peakfqsa and not requires_network and not requires_peakfqr_testdata

# PYTHONSAFEPATH=1 stops Python prepending the working directory to sys.path,
# which is what the `pytest` console script does and `python -m pytest` does
# not. Without it a local run can pass while CI fails on imports.
PYTEST := PYTHONSAFEPATH=1 $(PYTHON) -m pytest

.DEFAULT_GOAL := help
.PHONY: help check lint fmt test test-all fortran golden clean

help:  ## Show this help
	@awk -F':.*?## ' '/^[a-z-]+:.*## /{printf "  %-10s %s\n", $$1, $$2}' Makefile

check: lint test  ## Everything CI checks, in CI's order

lint:  ## Formatting check (does not modify files)
	$(PYTHON) -m black --check --diff $(PKGS)
	$(PYTHON) -m isort --check-only --diff $(PKGS)

fmt:  ## Apply formatting
	$(PYTHON) -m black $(PKGS)
	$(PYTHON) -m isort $(PKGS)

test:  ## Run the suite as CI does
	$(PYTEST) tests/ -m "$(DESELECT)"

test-all:  ## Run everything, including network and reference-data tests
	$(PYTEST) tests/

fortran:  ## Build the f2py extension from vendor/peakfqr (needs gfortran + meson)
	$(PYTHON) build_fortran/build.py

golden:  ## Regenerate Fortran parity golden files (needs the extension)
	$(PYTHON) tools/gen_fortran_golden.py

clean:  ## Remove build and test artifacts
	rm -rf build_fortran/mbuild build_fortran/native.ini build_fortran/_emafort*.so
	rm -rf hydrolib/peakfqr/_emafort*.so .pytest_cache
