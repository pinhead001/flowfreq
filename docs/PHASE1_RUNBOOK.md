# Phase 1 Runbook — Create `flowfreq`

> **STATUS: COMPLETE.** `pinhead001/flowfreq` exists, `main` is at `c76a3e6` with
> all 174 commits, and `v0.3.0` is tagged and pushed. What follows is the record
> of how it was done, kept because Phase 2 references its conventions and because
> the traps called out below are real ones this run hit.

Step-by-step for the library half of the split. Companion to
`docs/REPO_SPLIT_PLAN.md`; that document holds the reasoning, this one holds the
commands. Phase 2 (the app repo) is a separate runbook and must not start until
Phase 1 is tagged.

Every command runs from your machine, not from a Claude session — Phase 1 needs
repo-creation rights on the `pinhead001` account.

---

## Step 0 — Prerequisite: get Phase 0 onto `main`

**Phase 0 is not on `main`.** It lives on `claude/repo-split-planning-m4bj30`
(commits `2f41e36`, `abfe9fb`, `653289d`). Cloning `main` today gives you a
repo with the old `app/ffa_runner.py` and the `../data/*.csv` packaging bug, and
the split boundary would have to be redone by hand.

Merge it first:

```bash
gh pr create --base main --head claude/repo-split-planning-m4bj30 \
  --title "Split analysis into flowfreq.workflow (repo-split Phase 0)" \
  --body "See docs/REPO_SPLIT_PLAN.md"
# review, then merge
gh pr merge --squash --delete-branch
```

Do not proceed until `git log origin/main` shows the refactor. Everything below
assumes `main` contains it.

> Squash-merging collapses the three commits into one. That is fine for the
> history `flowfreq` inherits — the reasoning survives in the squash message
> and in `docs/REPO_SPLIT_PLAN.md`. Use `--merge` instead if you would rather
> keep them separate.

---

## Step 1 — Create the empty repository

GitHub does not allow forking a personal repo into the account that owns it, so
this is a create-then-push, not the Fork button. The history is preserved either
way.

```bash
gh repo create pinhead001/flowfreq --public \
  --description "Bulletin 17C flood frequency analysis with USGS data retrieval"
```

Create it **empty** — no README, no .gitignore, no license. Anything GitHub adds
becomes a conflicting root commit you then have to reconcile against the pushed
history.

---

## Step 2 — Push the full history

```bash
git clone https://github.com/pinhead001/hydrolib flowfreq
cd flowfreq
git remote set-url origin https://github.com/pinhead001/flowfreq
git push -u origin main
```

Check: `git log --oneline | wc -l` should report **174** commits, and the GitHub
repo should show the same. If it shows 1, you cloned wrong or the repo was not
empty.

**If it reports far fewer, your clone is shallow** and the push will fail with
`remote unpack failed: index-pack failed` / `did not receive expected object` —
the pack references objects past the shallow boundary that the remote cannot
resolve. `git log` looks plausible, which is what makes it easy to miss; check
for `.git/shallow`. Fix it with a bounded deepen, not `--unshallow`, which some
repo policy hooks block:

```bash
[ -f .git/shallow ] && git fetch --depth=1000 origin main
```

This bit the real Phase 1 run: a shallow clone reported 110 commits when the
history is 174, and cost two failed pushes before the cause was clear.

---

## Step 3 — Delete the app half

```bash
git rm -r app
git rm tests/test_streamlit_app.py tests/test_ffa_runner.py tests/test_ffa_export.py
git rm docs/vignette_streamlit_local.md docs/vignette_streamlit_web.md
```

Do **not** delete these — they are library code that only sounds app-shaped:

| Path | Why it stays |
|---|---|
| `flowfreq/freq_plot.py` | imports matplotlib, not Streamlit; returns a `Figure` |
| `tests/test_freq_plot.py` | tests the above |
| `tests/test_workflow.py` | the analysis half that came out of `ffa_runner` in Phase 0 |
| `tests/fixtures/big_sandy.py` | used by six library test modules |

---

## Step 4 — Config edits

### `pyproject.toml`

```toml
name = "flowfreq"      # was "hydrolib"
version = "0.3.0"        # was "0.2.0"
```

Leave everything else. In particular `[tool.setuptools.packages.find]` still
says `include = ["hydrolib*"]` and `[project.scripts]` still says
`flowfreq = "flowfreq.cli:cli"`. Repo, distribution, import name and console
script now all read `flowfreq` — the rename landed before this runbook runs, so
Phase 1 edits no import statements. If you find yourself editing one, something
has gone wrong.

### `.bumpversion.cfg`

```ini
current_version = 0.3.0
```

### `flowfreq/__init__.py`

Line 156: `__version__ = "0.3.0"`.

> `.bumpversion.cfg` drives both files, so after this one manual sync all future
> bumps are `bump2version minor` and nothing else.

### `Makefile`

Drop `app/` from `PKGS`, and delete the `smoke` target with its comment block
and its `.PHONY` entry:

```make
PKGS := hydrolib/ tests/
```

Delete the two-line comment above `PKGS` that explains why `app/` is in it — it
documents a decision that no longer applies, and a stale rationale is worse than
none.

`.PHONY: help check lint fmt test test-all fortran parity golden clean`

### `.github/workflows/tests.yml`

Delete the whole `app:` job — the `# app/ had no tests...` comment block through
`run: make smoke` (currently lines 59–78). Keep `lint`, `test`, and `fortran`
exactly as they are, including the deliberate non-chaining of `lint` and `test`
that the file's own comment explains.

### `.github/CODEOWNERS`

Delete the `app/ @pinhead001` entry and its comment.

---

## Step 5 — Docs

Rewrite these to describe a library-only repo. `CLAUDE.md` matters most: it is
instructions to an agent, and one that still claims `app/` exists will send the
next session hunting for files that are not there.

| File | Edit |
|---|---|
| `CLAUDE.md` | Drop `app/` from the Repository Layout block; drop the `make smoke` entry and its comment from Build & Development Commands |
| `README.md` | Drop the `app/*` rows from Module Overview (~L115–117), the `pip install streamlit` step (~L50), the `streamlit run` section (~L239), and the two Streamlit vignette rows (~L262–263). Add an install-from-git line: `pip install git+https://github.com/pinhead001/flowfreq@v0.3.0` |
| `CHANGELOG.md` | New `0.3.0` entry: renamed to `flowfreq` (package, distribution and display name), app split to its own repo, `flowfreq.workflow` added, gage table moved into the package |
| `PRODUCTION.md` | L141 mentions `app/` is smoke-tested only — now false |

`docs/FORTRAN_UPLOAD.md` references the Streamlit vignettes at L42 and L452, but
as a historical account of what those documents once claimed. Leave it; editing
a record of past state to match present state destroys the record.

---

## Step 6 — Verify before you tag

Run all of these. The last one is the one that matters most, because it is the
only failure a source checkout cannot show you.

```bash
pip install -e ".[dev]"
make lint                              # black + isort, must exit 0
PYTHONSAFEPATH=1 python -m pytest tests/   # expect 476 passed, 4 skipped, 5 xfailed
grep -rn "app\.\|app/" hydrolib/ tests/    # expect no hits
make parity                            # needs gfortran + meson
```

476 is the combined repo's 502 minus the 26 app tests deleted in Step 3
(18 in `test_streamlit_app.py`, 4 each in `test_ffa_runner.py` and
`test_ffa_export.py`). A different total means Step 3 removed the wrong thing.
The 5 xfails are expected and documented in `CLAUDE.md`; one is
`xfail(strict=True)` on Cains Coulee's weighted skew, so a *passing* xfail
fails the build on purpose.

Packaged install — build a wheel, install it into a clean venv, and read the
gage table back from a directory that has no `data/` in it:

Clear stale build output first. `build/` and `*.egg-info/` are gitignored, so a
tree that ever built under the old name still holds `build/lib/hydrolib/`, and
setuptools will happily fold it into the new wheel — producing a distribution
that ships **both** packages and recreates the very collision the rename
removed. It bit this rename during verification.

```bash
rm -rf build *.egg-info
python -m build --wheel
python -c "
import zipfile, glob
tops = {n.split('/')[0] for n in zipfile.ZipFile(glob.glob('dist/*.whl')[0]).namelist()}
assert 'hydrolib' not in tops, f'stale hydrolib in wheel: {tops}'
print('wheel is clean:', sorted(tops))
"
python -m venv /tmp/pyhl && /tmp/pyhl/bin/pip install dist/flowfreq-0.3.0-*.whl
mkdir -p /tmp/elsewhere && cd /tmp/elsewhere
/tmp/pyhl/bin/python -c "
from flowfreq.usgs import GageAttributes
st = GageAttributes.status()
assert st['file_exists'] and st['num_gages'] == 3, st
assert GageAttributes.get_site_name('09355500') == 'SAN JUAN RIVER NEAR ARCHULETA NM'
print('PASS')
"
```

Two things this proves that `pytest` cannot: the wheel is named `flowfreq`
and the package inside it, and the gage CSV resolves from
site-packages rather than from a sibling directory that only exists in a
checkout.

---

## Step 7 — Commit, push, tag

```bash
git add -A
git commit    # see message sketch below
git push origin main
git tag -a v0.3.0 -m "First release as flowfreq; app split to its own repo"
git push origin v0.3.0
```

The tag is not optional and not cosmetic: Phase 2 pins the app to
`git+https://github.com/pinhead001/flowfreq@v0.3.0`, and that reference cannot
resolve until the tag exists on the remote.

**Push the tag from your own machine.** A Claude Code web session can push
branches to this repo but not tags — `git push origin v0.3.0` returns HTTP 403
from GitHub (not the egress proxy, which records no denial), for annotated and
lightweight tags alike. That is a credential boundary, not a transient failure,
so retrying will not help. This is the one step of Phase 1 that cannot be
automated from a session.

Commit message should say what a reader six months out needs: that this repo is
the analysis half of a split, that it is the former `hydrolib` under a new name
did not, and that the app moved rather than being deleted.

---

## Step 8 — Confirm CI

Watch the first run on `main`:

```bash
gh run watch
```

Expect three jobs — `lint`, `test` (3.9–3.12), `fortran`. If a fourth appears,
the `app:` job survived Step 4 and will fail on missing `app/requirements.txt`.

---

## Done when

- [ ] `main` has 174 commits and no `app/` directory
- [ ] `v0.3.0` is pushed and visible in `gh release list` / the tag list
- [ ] CI green on three jobs
- [ ] A clean-venv wheel install resolves the gage table
- [ ] `CLAUDE.md` describes only what is in the repo

Phase 2 starts from a green Phase 1. The original `pinhead001/hydrolib` stays
**unarchived** until the app repo is deployed — archiving it early takes the
live Streamlit deployment down with it.
