# Repository Split Plan

Splitting `pinhead001/hydrolib` into a reusable analysis library and a separate
Streamlit application.

## 1. Decisions

| Question | Decision |
|---|---|
| Library repo | `pinhead001/flowfreq` |
| Library distribution name | `flowfreq` (on PyPI, and in `pip install`) |
| Library **import** name | `flowfreq` — `from flowfreq.bulletin17c import Bulletin17C` |
| App repo | `pinhead001/flowfreq-app` — new repo, separate from both |
| Git history | preserved in both repos; the original repo is archived, not deleted |
| App → library dependency | pinned git URL, bumped deliberately |
| `app/ffa_runner.py` | split — analysis half moves into the library, formatting stays in the app |

Repo, distribution and import name all agree. That is a change of plan, and the
reason matters:

**The `hydrolib` import name was already taken.** Deltares publishes `hydrolib`
(0.5.1) and `hydrolib-core` (1.1.0) on PyPI, both installing a top-level
`hydrolib/` package — and `hydrolib-core` ships `hydrolib/core/`, colliding
directly with what was then this project's `hydrolib/core.py`. Installed together in a clean
venv, verified in both orders:

- ours then Deltares: our entire API disappears — `hydrolib.Bulletin17C` gone,
  `from hydrolib.core import kfactor` raises ImportError
- Deltares then ours: **both** break; `import hydrolib` fails on `MAX_ABS_SKEW`

`pip check` reports no problem in either case. Nothing warns.

An earlier draft of this plan kept `import hydrolib` and renamed only the
distribution, citing `scikit-learn`/`sklearn` and `pillow`/`PIL`. Those precedents
hold only because nobody else owns the import name. Renaming the distribution
alone would have fixed PyPI publication and left the collision fully intact,
aimed squarely at this library's own audience — people doing hydrology in Python
are exactly the people likely to have Deltares tooling installed.

The project display name changed with it, `HydroLib` to `FlowFreq`, because
`HYDROLIB` is Deltares' project name: keeping it would have preserved the
confusion the rename exists to end.

`flowfreq` covers what the library actually does — flood frequency *and* low-flow
frequency — and carries no `hydro-` prefix to be mistaken for its neighbour.

## 2. Why the split is cheap

The coupling is already one-way and shallow. `app/` imports the library; nothing in
`flowfreq/` imports `app`. The only library module that even mentions Streamlit is
`flowfreq/freq_plot.py`, and that is a docstring and a function name — it imports
matplotlib and returns a `Figure`, not a Streamlit call.

Four real boundary defects have to be fixed as part of the split. They are listed in
§4; none is large, but each one silently breaks a split repo if it is skipped.

## 3. File allocation

### → `flowfreq` (library)

```
flowfreq/                 all 25 modules (renamed from hydrolib/)
  peakfqr/                f2py bridge (built, gitignored)
  validation/             comparison engine, benchmarks, reports, reference loaders
  workflow.py             NEW — see §4.1
  data/gage_attributes.csv   MOVED from top-level data/ — see §4.2
vendor/peakfqr/           USGS peakfq 8.1.0 reference (7.2 MB, do not edit)
build_fortran/            build.py + _emafort.pyf
tools/gen_fortran_golden.py
tests/                    everything except the three app test modules
  fixtures/               all 9 fixture modules
  fortran_parity/         parity suite + golden files
  validation/, integration/
  test_workflow.py        NEW — the run_ffa half of test_ffa_runner.py
examples/                 notebooks, scripts, user guide PDF
docs/FORTRAN_UPLOAD.md, vignette_cli.md, vignette_jupyter.md,
     vignette_lowflow_regime.md
.gitattributes            Fortran column-sensitivity rules — load-bearing here
AGENT_BUILD_INSTRUCTIONS_Claude.md, CHANGELOG.md, TODO.md, PRODUCTION.md
```

### → `flowfreq-app` (application)

```
app/streamlit_app.py      minus the sys.path hack (§4.3)
app/ffa_export.py         unchanged
app/ffa_runner.py         formatting/display half only (§4.1)
app/requirements.txt      rewritten (§4.4)
tests/test_streamlit_app.py
tests/test_ffa_export.py
tests/test_ffa_runner.py  formatter tests only
tests/fixtures/big_sandy.py   COPIED — 74 lines of static gage data (§4.5)
docs/vignette_streamlit_local.md, vignette_streamlit_web.md
```

### Archived with the original repo

Nothing is lost — `pinhead001/hydrolib` stays readable as the pre-split record.

## 4. Boundary work

### 4.1 Split `app/ffa_runner.py`

423 lines mixing Bulletin 17C computation with pandas display formatting. The
computation half is exactly what a second consumer would otherwise have to
reimplement, so it belongs in the library.

**Down into `flowfreq/workflow.py`:**

| Symbol | Note |
|---|---|
| `run_ffa()` | wraps `Bulletin17C`; the high-level entry point |
| `_low_outlier_source()` | helper for the above |
| `_B17C_DEFAULT_SKEW` / `B17C_DEFAULT_SKEW` | −0.302, England et al. 2019 |
| `DISPLAY_RETURN_INTERVALS`, `DISPLAY_AEP` | rename to `DEFAULT_*` — no longer display-specific |
| `_skew_values_from_result()` | pure extraction |
| `compute_skew_tables()` | computes quantiles + CIs per skew option |
| `build_skew_curves_dict()` | feeds `flowfreq.freq_plot`, which is already in the library |
| `SKEW_OPTIONS` | the label vocabulary both halves key on |

**Stays in `app/ffa_runner.py`:**
`format_parameters_df()`, `format_quantile_df()`, `build_station_summary_df()` — all
three exist to produce f-string-formatted DataFrames for `st.dataframe`, and bake in
column labels and decimal places that a library should not dictate.

`compute_skew_tables()` is the one judgement call: it returns
`Dict[str, pd.DataFrame]` keyed by display label. It moves because the DataFrame
holds computed numbers, not formatted strings — the formatting pass is a separate
function that stays up.

Re-export the moved names from `flowfreq/__init__.py` and add `workflow` to the
module list in that docstring.

**Do this refactor before the fork**, in the current repo. Then the split itself is a
pure file move, both new repos inherit the clean boundary, and the refactor is
verifiable against the existing test suite in one place.

### 4.2 Fix the `data/` packaging escape

`pyproject.toml` currently ships the CSV with:

```toml
"hydrolib" = ["../data/*.csv", "py.typed"]
```

A `..` in `package-data` reaches outside the package directory. It works with the
current setuptools but is not contractual, and it is why `data*` also has to appear in
`packages.find.include`. Move the file to `flowfreq/data/gage_attributes.csv` and drop
both hacks. `flowfreq/usgs.py:73-81` already searches
`Path(__file__).parent / "data"` as one of three candidate locations, so the move
needs no code change — only the removal of the now-dead sibling paths.

### 4.3 Drop the `sys.path` hacks

Three of them exist only because `app/` and `flowfreq/` share a checkout:

- `app/streamlit_app.py:23` — `sys.path.insert(0, parent)`, commented "needed for Streamlit Cloud"
- `tests/test_ffa_runner.py:11`
- `tests/test_ffa_export.py:13`

Once `flowfreq` is a declared pip dependency, all three are wrong rather than
merely ugly: they would shadow the installed package with whatever happens to sit
beside the app. Delete them.

### 4.4 App `requirements.txt`

```
streamlit>=1.28.0
flowfreq @ git+https://github.com/pinhead001/flowfreq@v0.3.0
numpy>=1.20.0
pandas>=1.3.0
matplotlib>=3.4.0
scipy>=1.7.0
```

`numpy`/`pandas`/`matplotlib`/`scipy` stay declared because `streamlit_app.py`
imports all four directly — a direct import should not rely on a transitive
dependency. `requests` is dropped; only the library's NWIS client uses it.

Streamlit Cloud installs git dependencies from `requirements.txt` without extra
configuration, so no `packages.txt` or secrets are needed. The pin is a tag, so a
library change cannot reach the deployed app until the pin is bumped.

### 4.5 Test fixtures across the boundary

`tests/test_streamlit_app.py:95` and `tests/test_ffa_runner.py` both import
`tests.fixtures.big_sandy`. That module is 74 lines of static USGS 03606500 peak data
and published expectations — copy it into the app repo's own `tests/fixtures/`.

Duplicating it is correct rather than lazy: it is historical gage data that cannot
drift, and the alternative is promoting a test fixture to public library API purely to
serve a downstream test. `tests/test_ffa_export.py` builds its own mock and needs
nothing.

### 4.6 CI, Makefile, lint scope

**`flowfreq`** keeps `.github/workflows/tests.yml` with the `app` job deleted, and
`Makefile` with `app/` removed from `PKGS` and the `smoke` target dropped. The
`lint`, `test` (3.9–3.12 matrix), and `fortran` jobs are unchanged — including the
deliberate non-chaining of `lint` and `test` documented in that workflow.

**`flowfreq-app`** gets a single-job workflow on Python 3.12: install
`requirements.txt`, `black --check`, `isort --check-only`, then run the three test
modules. `requires-python = ">=3.10"` — Streamlit's floor, which is why the app job is
separate from the matrix today.

The app repo needs no `pyproject.toml` for packaging (it is a deployed script, not a
distribution), but one is still worth adding to hold the `black`/`isort` config so the
line-length-100 and `profile = "black"` settings do not drift between repos.

### 4.7 Versioning

- `flowfreq` starts at **0.3.0** — a new distribution name and a changed public
  surface (`flowfreq.workflow`) both warrant the minor bump. Keep `.bumpversion.cfg`;
  it already drives `pyproject.toml` and `flowfreq/__init__.py`.
- `flowfreq-app` starts at **0.1.0** with its own independent version. It has no
  `.bumpversion.cfg` to inherit — the app's version is a display string, not a
  release contract.
- The library must **tag** releases from now on, because the app pins a tag. A PyPI
  release workflow (PRODUCTION.md §"Release Process") is not required for the git-URL
  pin and can follow later.

### 4.8 Optional: rename `plot_frequency_curve_streamlit`

`flowfreq/freq_plot.py` is named for a consumer it no longer lives with. It imports no
Streamlit and returns a plain matplotlib `Figure`. Renaming the function to
`plot_frequency_curve` (keeping the old name as an alias) and rewording the module
docstring costs one commit and removes a misleading signpost for the next consumer.
Not required for the split to work — sequence it after.

## 5. Execution sequence

### Phase 0 — Refactor in place (current repo) — **DONE**

1. Create `flowfreq/workflow.py`; move the §4.1 symbols; re-export from `__init__.py`.
2. Reduce `app/ffa_runner.py` to the formatters; update its imports.
3. Split `tests/test_ffa_runner.py` → `tests/test_workflow.py` (library) + a slimmed
   `tests/test_ffa_runner.py` (formatters).
4. Move `data/gage_attributes.csv` → `flowfreq/data/`; fix `pyproject.toml` and prune
   the dead search paths in `usgs.py`.
5. `make check && make smoke` — both green before anything is forked.

Phase 0 lands as a normal PR on `main`. It is independently valuable even if the split
stalls. Done in `abfe9fb`: 502 passed, 5 xfailed, lint clean, and the packaged gage
table verified from a clean venv install.

Two items in §4 are deliberately **not** in Phase 0 and remain for the repo split:

- The three `sys.path.insert` calls (§4.3). They are harmless while both trees share a
  checkout, and `app/` moves to its own repo root in Phase 2 where the import situation
  is different — removing them now would be a change made against conditions that are
  about to stop applying.
- `plot_frequency_curve_streamlit` (§4.8), which is a rename, not a boundary fix.

### Phase 1 — Create `flowfreq`

GitHub does not allow forking a personal repo into the account that owns it, so the
history-preserving move is a clone-and-push, not the Fork button:

```bash
gh repo create pinhead001/flowfreq --public
git clone https://github.com/pinhead001/hydrolib flowfreq && cd flowfreq
git remote set-url origin https://github.com/pinhead001/flowfreq
git push -u origin main            # full history, all 174 commits

git rm -r app tests/test_streamlit_app.py tests/test_ffa_export.py \
          tests/test_ffa_runner.py docs/vignette_streamlit_*.md
# pyproject: name = "flowfreq", version 0.3.0
# Makefile: drop app/ from PKGS, drop the smoke target
# workflow: delete the app job
# README/CLAUDE.md: rewrite for a library-only repo
git commit && git push
git tag v0.3.0 && git push --tags
```

### Phase 2 — Create `flowfreq-app`

Same clone-and-push, then the inverse deletion:

```bash
gh repo create pinhead001/flowfreq-app --public
# clone, repoint origin, push main
git rm -r hydrolib vendor build_fortran tools examples \
          tests/fortran_parity tests/validation tests/integration \
          docs/FORTRAN_UPLOAD.md docs/vignette_cli.md \
          docs/vignette_jupyter.md docs/vignette_lowflow_regime.md \
          .bumpversion.cfg AGENT_BUILD_INSTRUCTIONS_Claude.md
git rm $(ls tests/test_*.py | grep -v 'ffa_\|streamlit')
git rm tests/fixtures/*.py   # then restore big_sandy.py per §4.5
# write requirements.txt (§4.4), new workflow, new README/CLAUDE.md
```

Keeping the full history here is worth the clone even though only 35 commits touch
`app/` — `git log --follow` on the app files still resolves, and the deletion commit
is a legible marker of where the split happened.

### Phase 3 — Wire and verify

1. Push `flowfreq` and tag `v0.3.0` **before** the app's `requirements.txt` pin can
   resolve.
2. In a clean venv: `pip install -r requirements.txt && make smoke` in the app repo.
3. Repoint the Streamlit Cloud deployment at `pinhead001/flowfreq-app`, main branch,
   `app/streamlit_app.py`. Confirm the deployed app loads and runs one gage.

### Phase 4 — Archive

Archive `pinhead001/hydrolib` (Settings → Archive) once both repos are green and the
deployment is live. Add a one-line README note in each new repo pointing back to it.

## 6. Verification gates

| Gate | Command | Must show |
|---|---|---|
| Phase 0 | `make check` | full suite green, one known `xfail` on Cains Coulee |
| Phase 0 | `make smoke` | app imports and runs against the refactored boundary |
| flowfreq | `PYTHONSAFEPATH=1 python -m pytest tests/` | green on 3.9 and 3.12 |
| flowfreq | `grep -rn "app\." hydrolib/ tests/` | no hits |
| flowfreq | `make parity` | extension builds, golden files match |
| flowfreq | wheel into a clean venv, then `GageAttributes.status()` from a dir with no `data/` | `file_exists: True`, `num_gages: 3` (full command in the Phase 1 runbook §6) |
| flowfreq-app | `pip install -r requirements.txt` in a clean venv | pinned tag resolves |
| flowfreq-app | `grep -rn "sys.path\|from app\." . --include="*.py"` | no hits |
| flowfreq-app | `make test` | 26 passed against the *installed* library |
| Deployment | Streamlit Cloud | one gage analysed end to end |

The packaged-CSV check matters because it is the one failure the test suite cannot
catch from a source checkout: `usgs.py` finds the file via a sibling path that only
disappears once the package is installed elsewhere. Importing `USGSgage` is not
that check — the lookup lives on `GageAttributes`, and nothing reads the CSV until
something asks it for a gage.

## 7. Risks

**Streamlit Cloud repoint.** The deployment currently tracks
`pinhead001/hydrolib`; archiving that repo without repointing first takes the app
down. Phase 3 step 3 precedes Phase 4 for this reason.

**Fortran toolchain.** Only `flowfreq` needs gfortran + meson, and only in the
`fortran` CI job. The app repo requires no compiler — worth stating in its README,
since the vendored Fortran is the loudest thing in the current repo and its absence
will read as an omission.

**The `hydrolib` name on PyPI.** `flowfreq` is the distribution name, so nothing is
claimed under `hydrolib`. Anyone who previously installed this repo from a git URL as
`hydrolib` keeps working — the import name is unchanged — but their pin needs updating
to the new repo URL. Note it in the flowfreq CHANGELOG.

**Doc drift.** `README.md`, `CLAUDE.md`, `TODO.md` (71 KB) and `PRODUCTION.md` all
describe a combined repo. `CLAUDE.md` in particular is instructions to an agent and
will actively mislead if it still claims `app/` is present. Rewriting the first two is
Phase 1/2 work, not follow-up; `TODO.md` and `PRODUCTION.md` can be pruned afterwards.
