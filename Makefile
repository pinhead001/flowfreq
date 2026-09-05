# Developer entry points. CI calls these same targets, so what you run locally
# and what the build runs cannot drift apart.

PYTHON ?= python
PKGS := flowfreq/ tests/

# The marker deselection lives in pyproject.toml's addopts, not here, so a bare
# `pytest` is already correct and there is exactly one place to get it wrong.

# PYTHONSAFEPATH=1 stops Python prepending the working directory to sys.path,
# which is what the `pytest` console script does and `python -m pytest` does
# not. Without it a local run can pass while CI fails on imports.
PYTEST := PYTHONSAFEPATH=1 $(PYTHON) -m pytest

.DEFAULT_GOAL := help
.PHONY: help check lint fmt test test-all fortran parity golden clean

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
	$(PYTEST) tests/

test-all:  ## Run everything, including the network tests
	$(PYTEST) tests/ -m ""

fortran:  ## Build the f2py extension from vendor/peakfqr (needs gfortran + meson)
	$(PYTHON) build_fortran/build.py

# The import check is not redundant. test_live_vs_golden.py calls importorskip,
# so a failed build would skip every parity test and still exit 0 -- which is
# exactly the silent pass this target exists to prevent.
parity:  ## Build the extension and check the golden files against it (needs gfortran + meson)
	$(PYTHON) build_fortran/build.py
	$(PYTHON) -c "from flowfreq.peakfqr import emafitpr"
	$(PYTEST) tests/fortran_parity/

golden:  ## Regenerate Fortran parity golden files (needs the extension)
	$(PYTHON) tools/gen_fortran_golden.py

clean:  ## Remove build and test artifacts
	rm -rf build_fortran/mbuild build_fortran/native.ini build_fortran/_emafort*.so
	rm -rf flowfreq/peakfqr/_emafort*.so .pytest_cache
