# Phase 2 Runbook — Create the app repository

The Streamlit half of the split. Companion to `docs/REPO_SPLIT_PLAN.md` and
`docs/PHASE1_RUNBOOK.md`.

The repository is **`pinhead001/pyhydroapp`**, pairing with `pyhydrolib` on the
`py-` prefix. Settled before Step 2 deliberately: renaming after the Streamlit
Cloud repoint in Step 12 means performing that repoint twice.

---

## Step 0 — Prerequisites

Phase 2 cannot start until all three hold:

- [ ] `pyhydrolib` exists with `main` green on three CI jobs
- [ ] **`v0.3.0` is pushed as a tag** — `git ls-remote --tags` on the new repo must show it
- [ ] The original `pinhead001/hydrolib` is still **unarchived** and still serving the live app

The tag matters more here than anywhere else. `hydrolib` has never carried a
tag — Phase 1 Step 7 creates the first one in the project's history — and this
repo's `requirements.txt` pins to it. Without the tag on the remote, the app
cannot install its own dependency and Streamlit Cloud will fail the build with a
pip resolution error rather than anything that points at the real cause.

Verify before going further:

```bash
git ls-remote --tags https://github.com/pinhead001/pyhydrolib | grep v0.3.0
```

---

## Step 1 — Layout: flatten `app/` to the repo root

The combined repo nested everything under `app/` to keep it separate from
`hydrolib/`. Once the library is gone that nesting has nothing left to separate
from, so the app modules move to the repo root:

```
streamlit_app.py        was app/streamlit_app.py
ffa_runner.py           was app/ffa_runner.py
ffa_export.py           was app/ffa_export.py
requirements.txt        was app/requirements.txt
tests/
  __init__.py
  fixtures/{__init__.py,big_sandy.py}
  test_streamlit_app.py
  test_ffa_runner.py
  test_ffa_export.py
```

This is the Streamlit convention (`streamlit_app.py` at the root is what the
Cloud UI defaults to) and it is what removes the `sys.path.insert` hacks for
real. The hacks exist because `from app.ffa_runner import ...` needs the
repository *root* on `sys.path`, and Streamlit only ever puts the **script's own
directory** there. Flatten the layout and those two become the same directory,
so the imports resolve with nothing added.

> Plan §4.3 gave a different reason for deleting the hacks — that they would
> shadow the installed package. That reasoning was wrong: the hack inserts the
> repo root, which shadows nothing once `hydrolib/` is not beside it. The hacks
> are removable because of the flattening, not because the library moved. If you
> keep the `app/` subdirectory, **keep the hack in `streamlit_app.py`** or the
> deployed app will fail on `ModuleNotFoundError: No module named 'app'`.

Verified end to end against this layout before writing this runbook: 26/26 tests
pass, `streamlit run streamlit_app.py` serves HTTP 200 with no import errors,
and `import ffa_runner, ffa_export, hydrolib` all resolve from a foreign working
directory with **only** the script's directory on `sys.path`.

**Alternative:** keep `app/` as a subdirectory. Smaller diff, no import edits,
Streamlit Cloud's main file path stays `app/streamlit_app.py`. Cost is that the
`sys.path` hack in `streamlit_app.py` stays forever and a reader has to work out
why. Choose this only if you would rather not touch the Cloud config.

---

## Step 2 — Create the empty repository

```bash
gh repo create pinhead001/pyhydroapp --public \
  --description "Streamlit web app for Bulletin 17C flood frequency analysis"
```

Empty — no README, no .gitignore, no license, for the same reason as Phase 1: an
auto-created file becomes a conflicting root commit.

---

## Step 3 — Push the full history

```bash
git clone https://github.com/pinhead001/hydrolib pyhydroapp
cd pyhydroapp
git remote set-url origin https://github.com/pinhead001/pyhydroapp
git push -u origin main
```

Only 21 of the ~110 commits touch `app/`, but the clone is still worth it:
`git log --follow` on the app files keeps resolving, and the deletion commit
below becomes a legible marker of where the split happened.

---

## Step 4 — Delete the library half

```bash
git rm -r hydrolib vendor build_fortran tools examples \
          tests/fortran_parity tests/validation tests/integration \
          docs/FORTRAN_UPLOAD.md docs/vignette_cli.md \
          docs/vignette_jupyter.md docs/vignette_lowflow_regime.md \
          .bumpversion.cfg AGENT_BUILD_INSTRUCTIONS_Claude.md \
          docs/REPO_SPLIT_PLAN.md docs/PHASE1_RUNBOOK.md docs/PHASE2_RUNBOOK.md \
          flood_frequency_curve.png

# Library test modules: everything except the three app ones
git rm tests/test_bulletin17c.py tests/test_core.py tests/test_detrat.py \
       tests/test_freq_plot.py tests/test_lowflow.py tests/test_mse_ema.py \
       tests/test_p3_moments.py tests/test_r_fixtures.py tests/test_regime.py \
       tests/test_usgs.py tests/test_var_emab.py tests/test_var_mom.py \
       tests/test_workflow.py

# Fixtures: keep only big_sandy (Step 6 explains why)
git rm tests/fixtures/fortran_respec.py tests/fixtures/hu02_stations.py \
       tests/fixtures/moments_wymt.py tests/fixtures/nwis_rdb.py \
       tests/fixtures/paths.py tests/fixtures/skew_weighting.py \
       tests/fixtures/wymt_peaks.py tests/fixtures/wyoming_montana.py
```

`.gitattributes` can go too — its rules protect column-sensitive Fortran that no
longer lives here. Harmless either way.

---

## Step 5 — Flatten and rewrite the imports

```bash
git mv app/streamlit_app.py app/ffa_runner.py app/ffa_export.py .
git mv app/requirements.txt .
rmdir app
```

Then, in the four files that reference the old package path:

| Replace | With | Files |
|---|---|---|
| `from app.ffa_export import` | `from ffa_export import` | `streamlit_app.py`, `tests/test_ffa_export.py` |
| `from app.ffa_runner import` | `from ffa_runner import` | `streamlit_app.py`, `tests/test_ffa_runner.py` |
| `"app.streamlit_app"` | `"streamlit_app"` | `tests/test_streamlit_app.py` L57, L61 |
| `app/streamlit_app.py` | `streamlit_app.py` | docstrings and messages, L1, L67 |

Delete the three `sys.path` blocks — in `streamlit_app.py` (the one commented
"needed for Streamlit Cloud"), `tests/test_ffa_runner.py`, and
`tests/test_ffa_export.py` — along with the now-unused `import sys` and
`from pathlib import Path` that only served them. Leave `pathlib` where the file
still uses it: `tests/test_streamlit_app.py` reads the app source at L106/L111.

**Do not touch any `hydrolib` import.** `from hydrolib.workflow import run_ffa`
is correct and stays correct — the distribution is named `pyhydrolib`, the
package inside it is `hydrolib`.

---

## Step 6 — Keep the Big Sandy fixture

`tests/test_streamlit_app.py` (L95) and `tests/test_ffa_runner.py` import
`tests.fixtures.big_sandy`. Keep the file; the other eight fixtures go.

Keep `tests/__init__.py` as well, and read its docstring before deciding it is
redundant — it is load-bearing. Because `tests` is a package, pytest puts the
**repository root** on `sys.path`, which is what lets both `tests.fixtures...`
and the flat `ffa_runner` resolve under the `pytest` console script CI uses.
Delete it and collection fails in CI while still passing under
`python -m pytest` locally.

Duplicating 74 lines of static USGS 03606500 peak data is the right trade here:
it is historical gage data that cannot drift, and the alternative is promoting a
test fixture to public library API to serve a downstream test.

---

## Step 7 — `requirements.txt`

At the repo root, which is where Streamlit Cloud looks:

```
streamlit>=1.28.0
pyhydrolib @ git+https://github.com/pinhead001/pyhydrolib@v0.3.0
numpy>=1.20.0
pandas>=1.3.0
matplotlib>=3.4.0
scipy>=1.7.0
```

`numpy`, `pandas`, `matplotlib` and `scipy` stay declared because
`streamlit_app.py` imports all four directly, and a direct import should not
lean on a transitive dependency. `requests` is dropped — only the library's NWIS
client used it.

Bumping the library later is a one-line edit to the tag here. That is the point
of pinning: a `pyhydrolib` regression cannot reach the deployed app on its own.

---

## Step 8 — New config files

### `pyproject.toml`

The app is a deployed script, not a distribution, so this exists only to hold
the lint config — without it `black` and `isort` drift between the two repos and
every cross-repo copy/paste reformats.

```toml
[tool.black]
line-length = 100
target-version = ['py310', 'py311', 'py312']

[tool.isort]
profile = "black"
line_length = 100
known_first_party = ["ffa_runner", "ffa_export", "streamlit_app"]

[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = ["test_*.py"]
addopts = "-v"
```

No `requires-python` key is needed with no `[project]` table, but note in the
README that Streamlit needs **≥ 3.10** — that floor is why the app was a
separate CI job in the combined repo.

### `Makefile`

```make
PYTHON ?= python
PKGS := . tests/
PYTEST := PYTHONSAFEPATH=1 $(PYTHON) -m pytest

.DEFAULT_GOAL := help
.PHONY: help check lint fmt test run

help:  ## Show this help
	@awk -F':.*?## ' '/^[a-z-]+:.*## /{printf "  %-8s %s\n", $$1, $$2}' Makefile

check: lint test  ## Everything CI checks, in CI's order

lint:  ## Formatting check (does not modify files)
	$(PYTHON) -m black --check --diff $(PKGS)
	$(PYTHON) -m isort --check-only --diff $(PKGS)

fmt:  ## Apply formatting
	$(PYTHON) -m black $(PKGS)
	$(PYTHON) -m isort $(PKGS)

test:  ## Run the suite as CI does
	$(PYTEST) tests/

run:  ## Serve the app locally
	streamlit run streamlit_app.py
```

### `.github/workflows/tests.yml`

```yaml
name: Tests

on:
  push:
    branches: [main]
  pull_request:
    branches: [main]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.12'   # Streamlit needs >= 3.10

      # Installs pyhydrolib from the pinned tag, so CI exercises the same
      # resolution Streamlit Cloud performs on deploy. A broken pin fails
      # here rather than in production.
      - name: Install
        run: |
          pip install -r requirements.txt
          pip install black isort pytest

      - name: Lint
        run: make lint

      - name: Test
        run: make test
```

One job, not the library's three: there is no Fortran toolchain to provision and
no version matrix to cover, since Streamlit sets the floor at 3.10 and the app
is deployed on one interpreter.

### `.github/CODEOWNERS`

```
* @pinhead001
```

---

## Step 9 — Docs

| File | Action |
|---|---|
| `README.md` | Rewrite. Purpose, `pip install -r requirements.txt`, `streamlit run streamlit_app.py`, the live URL, Python ≥ 3.10, and a link to `pinhead001/pyhydrolib` for the analysis code |
| `CLAUDE.md` | Rewrite. An agent instruction file describing a library that is no longer here will send the next session hunting for `vendor/peakfqr` and `make parity` |
| `docs/vignette_streamlit_local.md` | Keep; update paths (`streamlit run streamlit_app.py`, no repo-root caveat) |
| `docs/vignette_streamlit_web.md` | Keep; rewrite §"requirements.txt" and the main-file-path step for the flat layout. L79's example deploy URL is derived from the repo and entry-point names, so it changes too |
| `CHANGELOG.md` | Reset to `0.1.0` — the app's version is a display string, not a release contract |
| `SECURITY.md`, `LICENSE` | Keep as-is |
| `TODO.md`, `PRODUCTION.md` | Both are library documents. Delete, or keep only the app-relevant sections |

`docs/vignette_streamlit_web.md` currently claims the library installs via
`pip install -e .` from the repo root. That is now false and actively
misleading — it is the pinned git dependency that installs it.

---

## Step 10 — Verify

The clean-venv install is the one that matters: it is the only step that
exercises the pinned git dependency the way Streamlit Cloud will.

```bash
python -m venv /tmp/appvenv
/tmp/appvenv/bin/pip install -r requirements.txt      # must resolve the v0.3.0 tag
/tmp/appvenv/bin/python -c "import hydrolib; print(hydrolib.__version__)"   # 0.3.0
```

If that pip call fails, stop — the tag is missing or the URL is wrong, and every
later step will fail in a way that points somewhere else.

```bash
make lint                                  # exit 0
make test                                  # expect 26 passed
grep -rn "sys.path\|from app\." . --include="*.py"    # expect no hits
```

Then boot it for real and confirm no import error:

```bash
streamlit run streamlit_app.py --server.headless true --server.port 8599
curl -s -o /dev/null -w "%{http_code}\n" http://localhost:8599/    # 200
```

Load `http://localhost:8599` in a browser and run one gage end to end. HTTP 200
only proves the server started; Streamlit executes the script on first client
connection, so a runtime error inside the app body will not show up until
something actually connects.

---

## Step 11 — Commit and push

```bash
git add -A
git commit
git push origin main
```

Say in the message that this is the app half of a split, that the analysis moved
to `pinhead001/pyhydrolib` rather than being deleted, and that the layout was
flattened so the `sys.path` hacks could go.

Watch CI: `gh run watch`. It must pass **before** the repoint — a red build here
means Streamlit Cloud will deploy something broken.

---

## Step 12 — Repoint Streamlit Cloud

At [share.streamlit.io](https://share.streamlit.io), on the existing app:

| Setting | New value |
|---|---|
| Repository | `pinhead001/pyhydroapp` |
| Branch | `main` |
| Main file path | `streamlit_app.py` (was `app/streamlit_app.py`) |

The main file path changes because of the Step 1 flattening. Missing it is the
single most likely way to deploy a 404 — the build succeeds and the app then
cannot find its entry point.

Some Cloud accounts do not allow changing the repository on an existing app. If
so, delete the app and create it fresh from the new repo; the URL is
reassignable if you set the same custom subdomain.

Confirm the deployed app loads and analyses one gage before Step 13.

---

## Step 13 — Archive the original

**Only now.** `pinhead001/hydrolib` is what the live deployment pointed at until
Step 12; archiving it before the new app is confirmed serving takes the app down.

```
Settings → Archive this repository
```

Add a line to both new READMEs pointing back to it as the pre-split record.

---

## Done when

- [ ] `pyhydroapp` CI green
- [ ] `pip install -r requirements.txt` resolves `pyhydrolib@v0.3.0` in a clean venv
- [ ] Deployed app loads from the new repo and analyses a gage
- [ ] No `sys.path` manipulation left in any file
- [ ] `CLAUDE.md` describes only what is in the repo
- [ ] `pinhead001/hydrolib` archived, linked from both new READMEs
