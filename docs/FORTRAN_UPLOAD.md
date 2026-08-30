# Uploading the peakfqr Fortran Source & Reference Materials

**Purpose:** get the USGS Fortran EMA source and its reference/test material out of the
local-only workspace (`C:\a\hal\_shared\peakfqr`) and into this repository, so that the
remaining numerical errors can be root-caused against the authoritative implementation and
a parity test suite can be built on direct Fortran-vs-Python comparison.

Audience: whoever is working on the machine that holds `C:\a\hal\_shared\peakfqr`.
Everything below runs there.

**If you only do one thing:** copy `_shared/peakfqr/src/emafit.f` into
`vendor/peakfqr/src/`. A previous session narrowed the last open numerical defect to three
routines inside that one file — `var_mom`, `EXPMOMCDERIV`, `DEXPECT` — and stopped rather
than guess at the formula (§1.4). Everything else here makes the work reproducible and
testable; that file makes it *possible*.

---

## 1. Why this is blocking

### 1.1 The Fortran bridge does not exist outside one machine

`hydrolib/peakfqr/` currently contains **only build output**, no source:

```
hydrolib/peakfqr/
├── __init__.py
├── _emafort.cp312-win_amd64.pyd   <- Windows-only, CPython 3.12-only binary
├── libgfortran-5.dll
├── libgcc_s_seh-1.dll
├── libquadmath-0.dll
└── libwinpthread-1.dll
```

`build_fortran/build.py` compiles it from `C:\a\hal\_shared\peakfqr\src`, an absolute path
that exists on exactly one computer. The consequences:

- `from hydrolib.peakfqr import emafitpr` raises `ModuleNotFoundError: No module named
  'hydrolib.peakfqr._emafort'` anywhere else — including this checkout (Linux, CPython 3.11).
- `.github/workflows/tests.yml` runs on `ubuntu-latest` across Python 3.9–3.12. The bridge
  can never load there, so **no CI job has ever exercised the Fortran path**.
- `docs/vignette_streamlit_web.md` claims the shipped `.pyd`/`.so` "should work on Linux
  (Streamlit Cloud runs Ubuntu)". That is not correct — a `.pyd` is a PE32+ Windows DLL.
  Fix that line as part of this work.

Without the sources in-repo, the reference implementation cannot be rebuilt, read in
context, or run by anyone but its author.

### 1.2 Three tests fail purely from missing reference data

`tests/peakfqsa/test_r_fixtures.py` resolves data through `../../../../_shared/peakfqr/...`,
i.e. it reaches *outside* the repository:

| Test | Needs |
|---|---|
| `TestSkewWeightingFixtures::test_load_whist_cases` | `inst/testdata/results_WHIST.csv` |
| `TestSkewWeightingFixtures::test_whist_case_values_in_range` | same |
| `TestMomentsWymtFixtures::test_expected_csv_files_exist` | the four `wymt_ffa_2022A_*` CSVs |

### 1.3 The remaining numerical failures need the Fortran to adjudicate

**Read the branch state before trusting any failure count.** Three branches carry different
Big Sandy results, and the two that matter are *both unmerged*:

| Branch | Big Sandy | Q100 CI lower | Q100 CI upper | Notes |
|---|---|---:|---:|---|
| `main` | 16 failed / 127 passed repo-wide | 1.19 % | 12.29 % | 11 quantiles + `std_log` also fail |
| `claude/library-overview-JbAcS` | 2 failed / 19 passed | 9.00 % | **10.53 %** | best upper CI, symmetric; degrades lower |
| `claude/hydrolib-edt-attributes-mm2aif` | **4 failed / 17 passed** | 3.36 % | 17.34 % | best quantiles; worst upper CI |

CI figures are measured at a common 5 % / 5 % tolerance. `library-overview-JbAcS` reports
"21 passed" on its own tolerances (10 % lower, 12 % upper) — it is green because its gates
are wider, **not** because it is more accurate. Compare deviations, never pass counts.

`claude/hydrolib-edt-attributes-mm2aif` is the current authoritative state: **4 failures.**

- `test_quantile[0.995]`, `test_quantile[0.99]` — the most *frequent* events, off ~2–2.7 %.
- `test_confidence_interval[0.02]`, `[0.01]` — the *rarest* events, upper bound only, off
  ~10–17 % and growing with return period. Lower bounds and the point quantiles at those
  same AEPs match.

### 1.4 What earlier sessions already established

Do not re-derive these. Recorded in `TODO.md` under "Open Questions" and in the commit
messages on the two branches above:

- **Root cause of the bulk of the old failures: found and fixed.**
  `_auto_configure_ema_params()` was guessing the historical perception threshold from
  `max(historical peak values)` — 25 000 instead of Big Sandy's actual 18 000. Honoring an
  explicitly-provided `perception_thresholds` entry resolved **9 of the original 13**
  failures (`f1d334c`). `b8fb275` on the other branch fixed the same bug independently.
- **Ruled out by independent verification.** The two lowest-level EMA primitives — truncated
  gamma moments and the LP3→gamma parameter transform — were checked against brute-force
  numerical integration and are **exact**. The residual is not there.
- **The standing diagnosis.** `FloodFrequencyAnalysis.compute_confidence_limits()` uses one
  closed-form asymptotic variance (`1/n + K²(1+0.75G²)/(2(n−1))`, the Bulletin 17B/MOM
  approximation) with `n = n_systematic`. `emafitpr` does not do this. It uses **Inverse
  Modified Cholesky Gaussian Quadrature** (added to `emafit.f` Oct 2012) and derives a
  **Pseudo Effective Record Length** (`as_G_PRL_o`) from the historical/censored portion of
  the record to widen and reshape the interval. Neither is implemented here. That gap
  matches the observed signature precisely: a defect confined to the upper CI at rare
  events, on a record whose censored portion is exactly what PRL is computed from.
- **The blocker.** Pinning down and *verifying* the replacement formula requires reading
  `emafit.f` — specifically `var_mom`, `EXPMOMCDERIV`, and `DEXPECT`. A previous session
  stopped here deliberately rather than guess at a formula, and left the 4 tests **failing
  rather than skipped**, on the grounds that this is a real open numerical question and not
  a missing-dependency one. That judgement is why this upload exists.

### 1.5 A merge decision has to happen too

`claude/library-overview-JbAcS` (4 commits: B17C Appendix 4 skew MSE, `n_intervals` skew
weighting, the skew-uncertainty CI variance term with `kfactor_skew_derivative()`, and a
systematic-only EMA regression test) is on **neither** `main` **nor**
`claude/hydrolib-edt-attributes-mm2aif`. It is stranded, and it is the only work that
attacks the upper-CI term directly — it is what took that error from ~19 % to 10.5 %.

The two branches are complementary and have never been combined: one has the better
quantiles, the other the better CI. Combining them is the obvious next experiment, and the
Fortran is what settles whether the combination is *correct* or merely closer. Land that
merge before generating any golden files, or the goldens will encode a state nobody ships.


## 2. What to upload

Vendor the material under `vendor/peakfqr/`, mirroring the upstream layout so paths in the R
package and in `TODO.md` still read correctly.

### Group 1 — Fortran sources (required)

Copy the **entire** `src/` directory; it is small and inter-dependent.

| Source | Destination | Why |
|---|---|---|
| `_shared/peakfqr/src/emafit.f` | `vendor/peakfqr/src/emafit.f` | **The single highest-value file.** `emafitpr` (authoritative EMA entry point) plus the three routines the open CI defect turns on — `var_mom`, `EXPMOMCDERIV`, `DEXPECT` — and the Inverse Modified Cholesky Gaussian Quadrature block. Also `mseg_all_sub`, `p3est_ema`, `detratsub`. |
| `_shared/peakfqr/src/dcdflib1.f90` | `vendor/peakfqr/src/dcdflib1.f90` | distribution functions |
| `_shared/peakfqr/src/imslfake.f` | `vendor/peakfqr/src/imslfake.f` | IMSL shims |
| `_shared/peakfqr/src/probfun.f` | `vendor/peakfqr/src/probfun.f` | `qP3sub`, `PLOTPOSHS`, probability functions |
| `_shared/peakfqr/src/*` (anything else) | `vendor/peakfqr/src/` | `Makevars`, `Makevars.win`, headers, any further `.f`/`.f90` |

### Group 2 — Call-convention reference (required)

| Source | Destination | Why |
|---|---|---|
| `_shared/peakfqr/R/fortranWrappers.R` | `vendor/peakfqr/R/fortranWrappers.R` | the specification `build_fortran/_emafort.pyf` was written against: exact argument order, log10 conversion, `pq = 1 − AEP`, output field extraction |
| `_shared/peakfqr/R/*.R` (the rest) | `vendor/peakfqr/R/` | surrounding logic: skew weighting, MGBT, PSF parsing |
| `_shared/peakfqr/DESCRIPTION`, `NAMESPACE` | `vendor/peakfqr/` | version provenance — records *which* peakfqr the numbers came from |

### Group 3 — Test data (unblocks §1.2, feeds §5)

| Source (under `_shared/peakfqr/inst/testdata/`) | Destination (under `vendor/peakfqr/inst/testdata/`) |
|---|---|
| `results_WHIST.csv` | same |
| `wymt_ffa_2022A.psf` | same |
| `wymt_ffa_2022A_WATSTORE.TXT` | same |
| `wymt_ffa_2022A_EXPinfo_7_4.csv` | same |
| `wymt_ffa_2022A_EXPdata_7_4.csv` | same |
| `wymt_ffa_2022A_EMPdata_7_4.csv` | same |
| `wymt_ffa_2022A_MGBT_7_5_1.csv` | same |
| `extra_tests/HU02/` (whole directory) | same — backs `tests/peakfqsa/fixtures/hu02_stations.py` |

### Group 4 — Upstream R tests (recommended)

| Source | Destination | Why |
|---|---|---|
| `_shared/peakfqr/tests/testthat/test-fortran.R` | `vendor/peakfqr/tests/testthat/` | origin of `tests/peakfqsa/fixtures/fortran_respec.py` — lets us verify the transcription |
| `_shared/peakfqr/tests/testthat/test-skewweight.R` | same | origin of `skew_weighting.py` |
| `_shared/peakfqr/tests/testthat/test-moments.R` | same | origin of `moments_wymt.py` |

### Do NOT upload

- Compiled objects: `*.o`, `*.mod`, `*.dll`, `*.so`, `*.pyd`, `src-x64/`, `mbuild/`.
- The R package's own `man/`, `.Rproj`, `.Rd` — not needed here.

---

## 3. Before you copy: three checks

**a. License / redistribution.** peakfqr is USGS-authored and normally public domain, but
confirm before vendoring: open `DESCRIPTION` and any `LICENSE`/`LICENSE.note` in the package
root and copy them to `vendor/peakfqr/` alongside the code. Add a short
`vendor/peakfqr/README.md` recording upstream name, version, retrieval date, and license.
If the license turns out to restrict redistribution, stop and raise it — do not push.

**b. Size.** Fortran sources and R are a few hundred KB. `extra_tests/HU02/` may not be.
Measure first:

```powershell
Get-ChildItem -Recurse "C:\a\hal\_shared\peakfqr\inst\testdata" |
  Measure-Object -Property Length -Sum |
  ForEach-Object { "{0:N1} MB" -f ($_.Sum / 1MB) }
```

Under ~50 MB: commit it plainly. Over that: upload the Group 3 files the tests actually
name, and subset `HU02/` to the stations listed in `tests/peakfqsa/fixtures/hu02_stations.py`.
Do not reach for Git LFS for text data of this size.

**c. Line endings — this one silently corrupts Fortran.** `.f` files are **fixed-form**:
columns 1–5 label, column 6 continuation, columns 7–72 statement. An editor that converts
tabs to spaces, strips trailing whitespace, or rewrites line endings can break compilation in
ways that are invisible in a diff. Land the `.gitattributes` in §4.1 **before** the first
`git add` of any Fortran file, and copy with a byte-preserving tool (`copy`, `robocopy`,
`cp`) — never paste through an editor.

---

## 4. Changes to land alongside the upload

### 4.1 `.gitattributes` (new file, repo root) — add this first

```gitattributes
# Fixed-form Fortran is column-sensitive: normalize to LF, never touch whitespace.
*.f      text eol=lf whitespace=-trailing-space,-tab-in-indent
*.f90    text eol=lf whitespace=-trailing-space,-tab-in-indent
*.pyf    text eol=lf
*.R      text eol=lf

# Reference fixtures are byte-exact expected output: do not normalize.
vendor/peakfqr/inst/testdata/**  -text
```

### 4.2 `.gitignore` — one trap

Line 7 is `*.so`. Once the extension builds on Linux/CI, `hydrolib/peakfqr/_emafort*.so`
will be **silently ignored**. That is the correct default (build output does not belong in
git), but it is worth making explicit so nobody loses an hour to it:

```gitignore
# Build output — the Fortran extension is compiled from vendor/peakfqr/src, never committed.
hydrolib/peakfqr/_emafort*.so
hydrolib/peakfqr/_emafort*.pyd
build_fortran/mbuild/
```

Verified: nothing in the current `.gitignore` blocks any path in the §2 manifest.

### 4.3 `build_fortran/build.py` — remove the hardcoded path

It currently pins `src = r"C:\a\hal\_shared\peakfqr\src"` and a `LOCALAPPDATA` Windows Store
Python path, so it runs on one machine only. Replace the path resolution with:

```python
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC = Path(os.environ.get("PEAKFQR_SRC", REPO_ROOT / "vendor" / "peakfqr" / "src"))

if not SRC.is_dir():
    sys.exit(
        f"Fortran sources not found at {SRC}.\n"
        "See docs/FORTRAN_UPLOAD.md — vendor/peakfqr/src must be populated."
    )

# Order matters: emafit.f references symbols from the others.
sources = [SRC / n for n in ("emafit.f", "dcdflib1.f90", "imslfake.f", "probfun.f")]
missing = [str(p) for p in sources if not p.is_file()]
if missing:
    sys.exit("Missing Fortran sources:\n  " + "\n  ".join(missing))

if sys.platform == "win32":
    mingw_bin = os.environ.get("MINGW_BIN", r"C:\msys64\mingw64\bin")
    if Path(mingw_bin).is_dir():
        os.environ["PATH"] = mingw_bin + os.pathsep + os.environ["PATH"]

for tool in ("gfortran", "meson"):
    if shutil.which(tool) is None:
        sys.exit(f"ERROR: {tool} not found on PATH")

build_dir = Path(__file__).resolve().parent
cmd = [
    sys.executable, "-m", "numpy.f2py",
    "-c", *(str(p) for p in sources),
    "-m", "_emafort",
    "--backend", "meson",
    "--build-dir", str(build_dir / "mbuild"),
]
print("Running:", " ".join(cmd))
sys.exit(subprocess.run(cmd, cwd=build_dir, env=os.environ).returncode)
```

This keeps the existing Windows/MSYS2 behaviour, adds an env override, and makes the build
work unchanged on Linux and macOS once the sources are in-repo.

### 4.4 Fixture paths — stop reaching outside the repository

Two places walk up four directory levels to `_shared/`. Point them at the vendored copy.

`tests/peakfqsa/fixtures/skew_weighting.py` — replace the `_WHIST_CSV` block:

```python
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
TESTDATA_DIR = REPO_ROOT / "vendor" / "peakfqr" / "inst" / "testdata"
_WHIST_CSV = str(TESTDATA_DIR / "results_WHIST.csv")
```

`tests/peakfqsa/test_r_fixtures.py` — replace the `testdata_dir` expression in
`TestMomentsWymtFixtures::test_expected_csv_files_exist`:

```python
from tests.peakfqsa.fixtures.skew_weighting import TESTDATA_DIR
testdata_dir = str(TESTDATA_DIR)
```

Confirm `parents[3]` resolves to the repo root from
`tests/peakfqsa/fixtures/skew_weighting.py` — `fixtures` → `peakfqsa` → `tests` → root.

### 4.5 `docs/vignette_streamlit_web.md`

Correct the claim that the bundled `.pyd`/`.so` works on Linux (§1.1). State instead that the
Fortran extension is optional, is built from `vendor/peakfqr/src` via `build_fortran/build.py`,
and that HydroLib falls back to the native EMA path when it is absent.

---

## 5. Upload procedure

Run on the machine holding `C:\a\hal\_shared\peakfqr`. Commands are Git Bash; PowerShell
equivalents noted where they differ.

**Step 1 — get onto the designated branch.**

```bash
cd /c/a/hal/hybrid-17c-cld
git fetch origin
git checkout -B claude/fortran-code-github-upload-5gj1s1 origin/claude/fortran-code-github-upload-5gj1s1
git status --short          # must be clean before you start
```

**Step 2 — land `.gitattributes` first (§4.1), so the Fortran is never normalized.**

```bash
# create .gitattributes with the §4.1 content, then:
git add .gitattributes
git commit -m "chore: add .gitattributes protecting fixed-form Fortran sources"
```

**Step 3 — copy the material.**

```bash
SHARED=/c/a/hal/_shared/peakfqr
mkdir -p vendor/peakfqr/{src,R,inst/testdata,tests/testthat}

cp -r "$SHARED"/src/*                       vendor/peakfqr/src/
cp -r "$SHARED"/R/*.R                       vendor/peakfqr/R/
cp    "$SHARED"/DESCRIPTION "$SHARED"/NAMESPACE  vendor/peakfqr/
cp    "$SHARED"/LICENSE*                    vendor/peakfqr/ 2>/dev/null || true

cd "$SHARED/inst/testdata"
cp results_WHIST.csv wymt_ffa_2022A.psf wymt_ffa_2022A_WATSTORE.TXT \
   wymt_ffa_2022A_EXPinfo_7_4.csv wymt_ffa_2022A_EXPdata_7_4.csv \
   wymt_ffa_2022A_EMPdata_7_4.csv wymt_ffa_2022A_MGBT_7_5_1.csv \
   /c/a/hal/hybrid-17c-cld/vendor/peakfqr/inst/testdata/
cp -r extra_tests /c/a/hal/hybrid-17c-cld/vendor/peakfqr/inst/testdata/

cd /c/a/hal/hybrid-17c-cld
cp "$SHARED"/tests/testthat/test-{fortran,skewweight,moments}.R vendor/peakfqr/tests/testthat/

# strip any build output that came along
find vendor/peakfqr -type f \( -name '*.o' -o -name '*.mod' -o -name '*.dll' \
     -o -name '*.so' -o -name '*.pyd' \) -delete
```

**Step 4 — verify before staging.**

```bash
# Fortran arrived intact and un-normalized
head -5 vendor/peakfqr/src/emafit.f
grep -c 'subroutine emafitpr' vendor/peakfqr/src/emafit.f     # expect >= 1
file vendor/peakfqr/src/*.f                                    # "ASCII text", not "CRLF"

# nothing in the manifest is being ignored
git status --porcelain --ignored vendor/ | grep '^!!' || echo "OK: nothing ignored"

# size sanity
du -sh vendor/peakfqr
```

**Step 5 — commit in reviewable pieces.**

```bash
git add vendor/peakfqr/src vendor/peakfqr/DESCRIPTION vendor/peakfqr/NAMESPACE vendor/peakfqr/LICENSE*
git commit -m "feat(fortran): vendor peakfqr Fortran sources (emafitpr and dependencies)"

git add vendor/peakfqr/R
git commit -m "docs(fortran): vendor peakfqr R wrappers as call-convention reference"

git add vendor/peakfqr/inst/testdata vendor/peakfqr/tests
git commit -m "test(fortran): vendor peakfqr reference test data and upstream R tests"

# then the §4.2–4.5 code changes
git add .gitignore build_fortran/build.py tests/ docs/
git commit -m "fix(fortran): resolve build and fixture paths from vendored sources"

git push -u origin claude/fortran-code-github-upload-5gj1s1
```

**Step 6 — rebuild from the vendored copy and confirm nothing regressed.**

```bash
python build_fortran/build.py
python -c "from hydrolib.peakfqr import emafitpr; print('bridge OK')"
pytest tests/ -q
```

### Definition of done

- [ ] `vendor/peakfqr/src/emafit.f` present, `subroutine emafitpr` greps clean.
- [ ] `python build_fortran/build.py` succeeds **from the vendored sources**, with no path
      outside the repo (temporarily rename `C:\a\hal\_shared\peakfqr` to prove it).
- [ ] The three `test_r_fixtures.py` failures from §1.2 pass.
- [ ] Big Sandy failure count on `claude/hydrolib-edt-attributes-mm2aif` stays at 4 —
      the upload must not perturb the numerics, only make them diagnosable.
- [ ] `git status --porcelain --ignored vendor/` shows no `!!` entries.
- [ ] `vendor/peakfqr/README.md` records upstream version, retrieval date, and license.

---

## 6. Root-causing the remaining 4 failures

Work the ladder below **in order** on Big Sandy (USGS 03606500), from the merged state
of §1.5. Each rung feeds the next; the *first* mismatch is the root cause and everything
below it is downstream noise. Do not skip ahead — a `cmoms` difference means nothing if the
two codes censored different observations.

What changed since this doc was first drafted: quantiles now agree at 10 of 12 AEPs, so
rungs 1–3 are **expected to pass**. Run them anyway — they are the cheap confirmation that
the two codes are fitting the same data, and they are what makes a rung 4–6 mismatch
trustworthy. The real work is rungs 4–6, and the headline target is rung 6.

| # | Fortran output | Question it answers | Expectation / if it differs |
|---|---|---|---|
| 1 | `gbval`, `gbnlow`, `gbnzero`, `gbns` | Same MGBT low-outlier threshold and censored count? | Expected to match. If not, the two codes are fitting **different data** — stop and fix MGBT; nothing downstream is meaningful. |
| 2 | `qlema`, `quema`, `tlema`, `tuema`, `nu` | Same EMA interval representation after thresholds? | Expected to match now that the 18 000 perception threshold is honored (§1.4). A mismatch means the fix is incomplete — check `lmissing = -80.0` and the ±1e20 conventions in `TODO.md` §2. |
| 3 | `cmoms[0,0]`, `cmoms[1,1]`, `cmoms[2,2]` | Same mean, variance, skew? | Expected to match within ~0.3 %. Compare all three columns — col 1 regional+at-site, col 2 at-site only, col 3 B17B MSE. |
| 4 | `as_G_mse_o`, `as_G_mse_Syst_o`, `Wdout` | Same skew MSE and determinant-ratio weight? | The `n_intervals`-vs-`n_observed` choice (84 vs 47) was reasoned about but never verified against Fortran. This is where that gets settled. Read `detratsub` and `mseg_all_sub`. |
| 5 | `yp` | Same quantiles given identical moments? | Two known residuals at AEP 0.99 / 0.995 (~2–2.7 %). Both are at the *frequent* end where `K` is small — suspect `qP3sub` / the K-factor, not the fit. |
| 6 | **`as_G_PRL_o`**, then `var_est`, `ci_low`, `ci_high` | Same Pseudo Effective Record Length, and same confidence limits? | **The headline defect (10–17 % upper CI at rare events).** Compare `as_G_PRL_o` *first*: it is a single scalar, it is the input the whole CI method is built on, and nothing here implements it. Read `var_mom`, `EXPMOMCDERIV`, `DEXPECT`, and the Inverse Modified Cholesky Gaussian Quadrature block (`emafit.f`, monotonicity enforcement at ~lines 485–491). |

Rung 6 is the reason for the upload. `as_G_PRL_o` is one number; getting it out of the
Fortran for Big Sandy tells you immediately whether the gap is the PRL value itself or the
quadrature built on top of it, and that single comparison is worth more than any amount of
further reasoning from the Python side.

### 6.1 Test architecture: commit Fortran outputs as golden files

The extension builds only where gfortran and the sources are present, but parity tests must
run in CI on Linux. Resolve that by **capturing Fortran output once and committing it**:

```
tools/gen_fortran_golden.py          # runs emafitpr, writes golden JSON (dev machine only)
tests/fortran_parity/
├── conftest.py                      # pytest.importorskip on the extension
├── golden/
│   ├── big_sandy_03606500.json
│   ├── wymt_<site>.json
│   └── hu02_<site>.json
├── test_native_vs_golden.py         # runs everywhere, incl. CI
└── test_live_vs_golden.py           # dev machine only; catches golden-file drift
```

Each golden file records the **full** `emafitpr` output for one input set — every field in
the §6 ladder, not just the quantiles — so a failure points at a rung rather than just
saying "the answer is wrong":

```json
{
  "meta": {"site": "03606500", "peakfqr_version": "8.1.0", "generated": "2026-..."},
  "inputs": {"n": 0, "ql": [], "qu": [], "tl": [], "tu": [], "dtype": [],
             "reg_M": 0.0, "reg_M_mse": 0.0, "reg_SD": 0.0, "reg_SD_mse": 0.0,
             "r_G": 0.0, "r_G_mse": 0.0, "gbthrsh0": -99.0, "pq": [], "eps": 0.0,
             "wght_opt_n": 1},
  "outputs": {
    "mgbt":  {"gbval": 0.0, "gbns": 0, "gbnzero": 0, "gbnlow": 0},
    "ema":   {"qlema": [], "quema": [], "tlema": [], "tuema": [], "nu": []},
    "cmoms": [[0,0,0],[0,0,0],[0,0,0]],
    "skew":  {"as_G_mse_o": 0.0, "as_G_mse_Syst_o": 0.0, "as_G_PRL_o": 0.0, "Wdout": 0.0},
    "quantiles": {"yp": [], "var_est": [], "ci_low": [], "ci_high": []}
  }
}
```

Design rules that make this pay off:

- **Store log10 space, exactly as Fortran returns it.** Convert to real space only in
  assertions. Round-tripping through `10^x` hides small discrepancies — precisely the ones
  being hunted.
- **Record inputs beside outputs.** A golden file whose inputs are implicit is unfalsifiable
  and cannot be regenerated.
- **Assert per rung, tightest first** — `gbval` and the counts exactly, `cmoms` at ~1e-9,
  `yp` at the existing 1 % / 2 % tolerances. A test that only checks Q100 tells you the
  answer is wrong; a laddered test tells you *where*.
- **Pin `peakfqr_version` in `meta`** and fail loudly if the live rebuild disagrees, so a
  future upstream bump cannot silently invalidate the goldens.
- Mark the live tests with the existing `requires_peakfqsa` marker (or add
  `requires_fortran`) so `pytest -m "not requires_fortran"` stays green in CI.

### 6.2 Once parity is established

Extend beyond Big Sandy to the multi-site data this upload brings in — `wymt_ffa_2022A`
(WY/MT, PeakFQ 7.4 expected output) and `extra_tests/HU02` (Northeast US, PeakfqSA 7.5.1).
Those cover EMA branches Big Sandy alone does not: historical-period censoring, zero flows,
and each of the three skew-weighting options. Fixtures for all three already exist under
`tests/peakfqsa/fixtures/`; they have simply never had reference data to run against.

---

## 7. Reference

- `TODO.md` §1–3 — `emafitpr` signature, `cmoms` layout, EMA interval conventions, MGBT and
  skew-MSE encodings. Written from the Fortran; use it as the map when reading `emafit.f`.
- `TODO.md` §4 "Confidence Interval Method" and the `as_G_PRL_o` row of the output table —
  the two paragraphs that describe what is missing from `compute_confidence_limits()`.
- `TODO.md` "Open Questions", last entry — the previous session's own statement of what it
  needed and why it stopped. This document is the answer to that entry.
- `build_fortran/_emafort.pyf` — the f2py signature, derived from `R/fortranWrappers.R`.
  If `emafit.f` disagrees with it, the `.pyf` is wrong and the bridge is silently mis-marshalling.
- `CLAUDE.md` — repository conventions, including that all Fortran-facing data is log10 base-10.
