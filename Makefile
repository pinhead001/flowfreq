# Developer entry points. CI calls these same targets, so what you run locally
# and what the build runs cannot drift apart.

PYTHON ?= python
# app/ is in here because CI once linted flowfreq/ and tests/ only, and that
# blind spot is why Streamlit changes could not be reviewed with confidence.
PKGS := flowfreq/ tests/ app/

# The marker deselection lives in pyproject.toml's addopts, not here, so a bare
# `pytest` is already correct and there is exactly one place to get it wrong.

# PYTHONSAFEPATH=1 stops Python prepending the working directory to sys.path,
# which is what the `pytest` console script does and `python -m pytest` does
# not. Without it a local run can pass while CI fails on imports.
PYTEST := PYTHONSAFEPATH=1 $(PYTHON) -m pytest

.DEFAULT_GOAL := help
.PHONY: help check lint fmt test test-all smoke fortran parity golden clean

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

# Importing app/streamlit_app.py executes the whole script in Streamlit's bare
# mode, so this exercises the app end to end short of a button press. Needs
# Streamlit, which needs Python >= 3.10 and is not in the dev extra.
smoke:  ## Smoke-test the Streamlit app (pip install -r app/requirements.txt first)
	$(PYTEST) tests/test_streamlit_app.py tests/test_ffa_runner.py tests/test_ffa_export.py

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
