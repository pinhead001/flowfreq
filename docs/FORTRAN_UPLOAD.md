# Uploading the peakfqr Fortran Source & Reference Materials

**Purpose:** get the USGS Fortran EMA source and its reference/test material out of the
local-only workspace (`C:\a\hal\_shared\peakfqr`) and into this repository, so that the
remaining numerical errors can be root-caused against the authoritative implementation and
a parity test suite can be built on direct Fortran-vs-Python comparison.

Audience: whoever is working on the machine that holds `C:\a\hal\_shared\peakfqr`.
Everything below runs there.

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

### 1.3 Thirteen tests fail on numerical divergence — and cannot be diagnosed from Python alone

Baseline on this checkout: **16 failed, 127 passed**. Three are §1.2; the other thirteen are
`tests/validation/test_big_sandy.py` (USGS 03606500, the PeakfqSA manual reference site):

| Quantity | Expected (PeakfqSA) | HydroLib native | Δ |
|---|---:|---:|---:|
| `std_log` | 0.289200 | 0.296061 | **+2.37 %** |
| Q1.005 (AEP 0.995) | 871.25 | 840.03 | −3.58 % |
| Q1.01 (AEP 0.99) | 1 045.59 | 1 015.55 | −2.87 % |
| Q2 (AEP 0.50) | 5 284.36 | 5 401.02 | +2.21 % |
| Q5 | 9 166.15 | 9 474.44 | +3.36 % |
| Q10 | 12 134.65 | 12 598.98 | +3.83 % |
| Q25 | 16 276.60 | 16 962.79 | +4.22 % |
| Q50 | 19 617.73 | 20 483.07 | +4.41 % |
| Q100 | 23 158.65 | 24 212.25 | +4.55 % |
| Q200 | 26 912.12 | 28 162.43 | +4.65 % |
| Q500 | 32 217.14 | 33 739.16 | +4.72 % |
| Q50 CI upper | 29 124.18 | 27 493.33 | −5.60 % |
| Q100 CI upper | 37 986.08 | 33 318.98 | **−12.29 %** |

Read the shape, not the individual numbers: the lower tail is **too low**, the upper tail is
**too high**, and the error grows monotonically with return period. That is a curve pivoting
about its mean — the signature of a fitted standard deviation roughly 2.4 % too large, which
`std_log` confirms directly. The confidence-interval error is a *separate* defect: the lower
limits are within ~1 %, so the variance path (`var_est` / the determinant-ratio weight) is
diverging on its own, not merely inheriting the sigma error.

At least three causes produce that same sigma signature, and **Python-side inspection cannot
tell them apart**:

1. MGBT selects a different low-outlier threshold, so the two codes are fitting *different
   censored data sets*.
2. The EMA moment iteration converges to a different second moment on identical data.
3. The skew estimate or its MSE weighting differs, changing the fitted curve.

Distinguishing them requires the Fortran's *intermediate* outputs — `gbval`, `gbnlow`,
`qlema`/`quema`, `cmoms` — on the same inputs. That is what this upload delivers. §6 gives
the ordered diagnostic procedure.

---

## 2. What to upload

Vendor the material under `vendor/peakfqr/`, mirroring the upstream layout so paths in the R
package and in `TODO.md` still read correctly.

### Group 1 — Fortran sources (required)

Copy the **entire** `src/` directory; it is small and inter-dependent.

| Source | Destination | Why |
|---|---|---|
| `_shared/peakfqr/src/emafit.f` | `vendor/peakfqr/src/emafit.f` | `emafitpr` (the authoritative EMA entry point), `mseg_all_sub`, `p3est_ema`, `var_mom`, `detratsub` |
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
- [ ] Failure count drops from 16 to 13 — the remaining 13 are the §1.3 numerical work.
- [ ] `git status --porcelain --ignored vendor/` shows no `!!` entries.
- [ ] `vendor/peakfqr/README.md` records upstream version, retrieval date, and license.

---

## 6. Root-causing the 13 numerical failures

With the Fortran runnable, work the diagnostic ladder below **in order** on Big Sandy
(USGS 03606500). Each rung feeds the next; the *first* mismatch is the root cause, and
everything below it is downstream noise. Do not skip ahead — a `cmoms` difference means
nothing if the two codes censored different observations.

| # | Fortran output | Question it answers | If it differs |
|---|---|---|---|
| 1 | `gbval`, `gbnlow`, `gbnzero`, `gbns` | Did MGBT pick the same low-outlier threshold and censor the same count? | The two codes are fitting **different data**. Fix MGBT first; nothing else is meaningful. |
| 2 | `qlema`, `quema`, `tlema`, `tuema`, `nu` | Same EMA interval representation after thresholds? | Interval/perception-threshold encoding bug. Check the `lmissing = -80.0` and ±1e20 conventions in `TODO.md` §2. |
| 3 | `cmoms[0,0]`, `cmoms[1,1]`, `cmoms[2,2]` | Same mean, **variance**, skew? | **Most likely home of the +2.37 % `std_log` error.** Compare all three columns — col 1 regional+at-site, col 2 at-site only, col 3 B17B MSE. |
| 4 | `as_G_mse_o`, `as_G_mse_Syst_o`, `Wdout` | Same skew MSE and determinant-ratio weight? | Skew weighting divergence — re-read `detratsub` and the Halloween/ERL/INV option handling. |
| 5 | `yp` | Same quantiles given identical moments? | Isolated `qP3sub` / K-factor issue rather than a fitting issue. |
| 6 | `var_est`, `ci_low`, `ci_high` | Same confidence limits? | **The separate −12.29 % upper-CI defect.** Lower limits are within ~1 %, so suspect the variance formula's asymmetry, not the fit. |

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
- `build_fortran/_emafort.pyf` — the f2py signature, derived from `R/fortranWrappers.R`.
  If `emafit.f` disagrees with it, the `.pyf` is wrong and the bridge is silently mis-marshalling.
- `CLAUDE.md` — repository conventions, including that all Fortran-facing data is log10 base-10.
