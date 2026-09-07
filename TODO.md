# TODO — FlowFreq Hybrid 17C Implementation

## Status
Last updated: 2026-09-06. Version **0.4.0**.
Tests: **595 passed, 7 skipped, 1 deselected, 5 xfailed** in ~46 s, via `make clean-verify`
with the extension genuinely absent (563 before the Fortran-as-a-selectable-engine work below
started). With the extension built, `pytest tests/fortran_parity/` (`make parity`'s exact
selection) adds 259 passed, 1 pre-existing platform-precision failure, 2 xfailed -- see "Done
-- the Fortran as a selectable engine" for the full account, including the one known failure
and why it is not this work's to fix.
CI green on main. (The 608 recorded here on 2026-09-02 was the pre-split suite, which still
carried the Streamlit app's tests; those moved to `flowfreq-app` in 0.3.0. Run the number, do
not carry it forward -- it has been wrong in a commit message and a PR body already.)
Fortran reference: **vendored** at `vendor/peakfqr/` (peakfq 8.1.0, CC0).
Fortran bridge: builds from those sources via `python build_fortran/build.py`
(gfortran + meson) and is now **built and checked in CI** by `make parity`.

### Environment constraints -- read before starting work

These bit repeatedly and are not discoverable from the code:

- **Tag pushes and branch deletes return HTTP 403** from a Claude Code session. This is a
  credential boundary, not a transient failure: retrying and re-authenticating do not help.
  Anything requiring a tag or a branch deletion has to be done from a local clone.
- **`v0.4.0` is now tagged** (confirmed via `git ls-remote --tags origin`, 2026-09-05) --
  `README.md`'s install line resolves. This item is done; no action needed.
- **NWIS (`nwis.waterdata.usgs.gov`) is blocked by the egress proxy.** Tests marked
  `requires_network` cannot run in a Claude session; they are deselected by default via
  `addopts` anyway. Use the committed fixtures under `tests/fixtures/`.
- **The Fortran extension is not built by default.** `make fortran` needs gfortran and meson.
  Everything marked `requires_fortran` auto-skips without it, which is why a local run can
  report fewer tests than CI's parity job.
- **Reproduce CI with `make clean-verify`, not a bare `pytest`.** Five separate local-pass /
  CI-fail incidents during the repo split traced to environment artifacts: a stale
  `build/lib/` shipping two packages in one wheel, the working directory shadowing an
  installed package, a `__pycache__`-only directory changing isort's first-party
  classification, and a mypy version skew. `clean-verify` wipes the tree first.
- **`vendor/` is a verbatim USGS reference copy. Do not edit anything under it** -- a change
  there silently invalidates every parity comparison made against it.

Every P1 and P2 item is done. **The entire `var_mom` port (TODO.md P3) is now done**, including
`detrat`, the at-site EMA moment-iteration fix, and the confidence-interval shape fix
(`VAR_EMAB`/`regmoms`/`ci_ema_m3b`) -- all wired into `Bulletin17C`. What's left is one small,
understood numerical residual, not a missing routine:

- Big Sandy (37 censored historical gap-year intervals, at-site skew under the HWN floor):
  weighted skew matches peakfq 8.1.0 to **2.4e-6**, every quantile to **≤0.06%**, mean/std to
  **~1e-5%** -- the same level Powder River (no censoring at all) has. Confidence intervals now
  match too: **both bounds within 0.06%** at every AEP tested, asymmetry ratio within 0.0007-0.0022
  of peakfq's own (was exactly 1.000 everywhere, symmetric by construction, off by as much as 0.4
  in ratio terms at the AEP-0.01 tail).
- Cains Coulee (11 PILFs from MGBT, the one parity case whose at-site skew clears the 0.04 HWN
  floor): at-site skew matches to **0.0002**, native `Wd` to **0.002** against peakfq's 0.184.
  One real residual remains: `skew_weighted` is still 0.058 skew units off. **Not** a
  `var_mom`/`mse_ema`/`mn2mvarb` precision limit, as earlier documentation here assumed and a
  deeper investigation has now corrected -- `mse_ema` called standalone with this site's real
  post-MGBT group matches the Fortran oracle to 3e-8 relative. The actual gap: `emafitpr`'s own
  internally-computed, reported `as_G_mse` for this site (0.2212) does not match what the same
  `mseg_all` Fortran routine gives called standalone with identical inputs (0.0749) -- a ~3x
  discrepancy, reproducible even from a from-scratch, single-case golden regeneration (ruling
  out cross-case state contamination), whose exact mechanism inside `emafitpr` was not pinned
  down despite substantial investigation. See
  `tests/fortran_parity/test_fortran_oracles.py::TestCainsCouleeAsGMseDiscrepancy` for the full,
  reproducible account. `tests/fortran_parity/test_wymt_vs_golden.py`'s `test_weighted_skew_matches`
  is the one still-open `xfail(strict=True)` this leaves in the whole `var_mom`/`detrat`/`var_emab`
  tree -- and, per its own reasoning, flowfreq's own computation may well be *more* correct than
  the golden reference here, not less.

The other four `xfail(strict=True)`s left in the suite are all the 2012 PeakfqSA manual
comparisons (`tests/validation/test_big_sandy.py`), which is a different, non-reproducible
reference (CLAUDE.md's Test Data section) -- not evidence of anything left to port.

---

## Open Items (prioritised)

### Done — the Fortran as a selectable engine

Specified in **`docs/FORTRAN_ENGINE_DESIGN.md`**; all five pieces from that doc's estimate
table are implemented, tested, and (unlike an earlier checkpoint of this same item) verified
against a **live** `emafitpr` call, not just committed goldens.

**Environment note, resolved**: an earlier checkpoint of this item recorded no
gfortran/meson/ninja available. MSYS2 (`C:\msys64\mingw64\bin`) was installed and the
extension now builds on Windows too -- `PATH` needs the mingw64 bin dir appended *after* the
existing `PATH` (MSYS2 ships its own `python.exe`, which would otherwise shadow the
numpy-having interpreter), and the built `.pyd` needs four runtime DLLs
(`libgcc_s_seh-1.dll`, `libgfortran-5.dll`, `libquadmath-0.dll`, `libwinpthread-1.dll`, all
from the same mingw64 `bin/`) copied alongside it in `flowfreq/peakfqr/` -- Windows does not
search `PATH` for an extension module's own DLL dependencies the way Linux's loader does.
Not written into `build_fortran/build.py` (a Windows-specific packaging step, out of scope for
this item), but worth remembering for whoever next builds this on Windows.

**A real, environment-specific gap found in `make clean-verify` itself, not fixed here**: the
`clean` target's `rm -rf flowfreq/peakfqr/_emafort*.so` only matches the Linux extension
suffix. On Windows the built artifact is `_emafort.cp3xx-win_amd64.pyd`, which that glob never
touches -- so `make clean-verify` on a Windows machine that has ever built the extension does
**not** actually test the "extension absent" state `clean-verify`'s own docstring promises; it
silently keeps testing with the extension present. Not a portability bug worth fixing in this
item's scope (CI only ever runs `make clean-verify`/`make parity` on Linux, where the glob is
correct), but it cost real time to notice here, and manually deleting the `.pyd` and the four
DLLs was needed to get a true baseline read. Verified both ways below.

- [x] **Library-side interval builder** -- `flowfreq/fortran_engine.py::build_emafit_arrays`.
      Translates a `Bulletin17C` input set (peaks, water years, historical peaks, perception
      thresholds, a low-outlier override, optional `EMAParameters`) into `emafitpr`'s
      `ql/qu/tl/tu/dtype`, following `siteQT` (`vendor/peakfqr/R/readInputs.R`) directly rather
      than reusing `ExpectedMomentsAlgorithm._build_flow_intervals` -- confirmed by reading that
      method (`bulletin17c.py` ~lines 763-803) that it only fills gap years over the *historical*
      range (via `self._ema_params.historical_start/end`), never over a perception-threshold
      period declared across the *systematic* range, which is exactly what site 12363000's
      gap-year case exercises. The prior research's claim about this limitation is correct.

      Verified, with actual numbers, all without needing the built extension:
      - **Byte-equality against `tests/fortran_parity/cases.py::build_emafit_inputs`** on all
        four registered parity cases (`tests/fortran_parity/test_interval_builder.py::
        TestMatchesExistingParityCases`) -- Big Sandy (44 systematic + 3 historical + 37
        censored gap rows = 84 total, matching `n_censored == 37` from the P3 table), Powder
        River, Cains Coulee (both 0 censored, contiguous), and site 12363000 (98 rows, gaps
        omitted). Compared as a sorted multiset of rows, not by iteration order -- the two
        builders append systematic/historical/fill in different block orders, and row order has
        no meaning to `emafitpr` (a sum over independent intervals).
      - **The 12363000 gap-year switch**, at the row-construction level: without a declared
        threshold the four 1924-1927 gap years are omitted (98 rows total, on a 5-year slice
        used for a fast unit test); with one declared over the span, they are censored
        (`ql=Qmin`, `qu=tl=`threshold, `tu=Qmax`) and counted in `n_censored`.
        **Now verified against a live `emafitpr` call**
        (`tests/fortran_parity/test_live_interval_builder.py::TestGapYearSwitchLive`), and it
        reproduces the design doc's published example on the nose: the 98-row (omitted)
        construction gives at-site skew **0.43504703453300786**, bit-identical to the golden
        file's own `cmoms[2][1]`, and Q100 **120,064.33** cfs (design doc: "120,064"); the
        102-row (censored) construction, with a synthetic 5,000 cfs threshold declared over the
        gap (any value from roughly 1,000-10,000 cfs gives the same result to 4 significant
        figures -- found empirically, not assumed), gives at-site skew **0.2500537568149216**
        and Q100 **119,472.61** cfs -- the design doc's own "+0.250" and "119,473" (rounded),
        matched independently, not curve-fit to hit them.
      - **Declared-but-vacuous-threshold omission** (`tl <= Qmin`, e.g. Big Sandy's own
        1930-1973 systematic threshold of literally 0.0): confirmed directly against
        `readInputs.R` lines ~1030-1051 (the `keepNoInfo` filter after the missing-year branch)
        and unit-tested (`TestVacuousThreshold`) -- a gap year under a vacuous threshold is
        dropped, same as an undeclared gap year, not turned into an exact peak at `Qmin`.
      - Zero flows, overlapping perception-threshold periods, `water_years` omitted, and the
        `gbthrsh0` encoding: all unit-tested directly (`TestZeroFlows`,
        `TestOverlappingPerceptionPeriods`, `TestWaterYearsOmitted`, `TestGbthrsh0Encoding`).

      **Open question 1 (zero flows) -- resolved and implemented.** A zero flow is an exact
      peak at `log10(Qmin) = -20`, not an omission or a special censored row -- confirmed
      against `readInputs.R` (`QT$ql <- QT$peak_va` sets it directly, later clamped by
      `QT$ql[QT$ql < Qmin] <- Qmin`). `EmafitArrays.systematic_peaks` records the *real* value
      (0.0) alongside the clamped log10 row, so the adapter's `pilf_flows` can report a true
      zero rather than `Qmin`.

      **Open question 2 (overlapping perception threshold periods) -- resolved and
      implemented.** `siteQT`'s own doc comment (readInputs.R ~lines 805-806): "the last one
      specified is given priority", implemented as a sequential year->threshold map overwrite in
      the order `perception_thresholds.items()` iterates (a plain dict, so insertion order).
      Unit-tested both directions (`TestOverlappingPerceptionPeriods`) to confirm it is really
      *last-declared*, not e.g. lowest-value or widest-range.

      **Open question 3 (`mgb_critical_value`) -- decided, with a correction to the prior
      research's justification.** The design doc's adapter table (§4) leaves this open between
      "recompute natively" and "report unavailable"; the decision here is recompute, via
      `flowfreq.core.grubbs_beck_critical_value(n)` -- same function, same policy the *native*
      EMA path already uses for this exact diagnostic field (it is not itself what determines
      MGBT's censoring in either engine). **The prior research's citation for why this is safe
      was wrong and needed correcting, not just spot-checking**: it claimed
      `grubbs_beck_critical_value` is "already verified line-by-line against the Fortran
      (GGBCRITP/FP_TNC_CDF, validated on Orestimba Creek)" -- that citation is CLAUDE.md's
      Validation Status section, but it describes a *different* function,
      `ExpectedMomentsAlgorithm._mgbt_pvalue` (the actual MGBT p-value machinery, whose docstring
      says "Direct Python translation of the GGBCRITP/FGGB Fortran functions" and whose variable
      names -- `EX1..EX4`, `MuM`, `MuS2`, `Lambda` -- match `probfun.f`'s `FGGB` one-for-one).
      `grubbs_beck_critical_value` is a wholly separate, simpler critical-value table for the
      *classical single-outlier* Grubbs-Beck test (`gbtype = 'GBT'`), not the default MGBT path.
      Checked directly against `vendor/peakfqr/src/emafit.f` lines 1020-1039: its `n >= 100`
      asymptotic branch (`-0.9043 + 3.345*sqrt(log10(n)) - 0.4046*log10(n)`) is an *exact*
      transcription of that file's line 1031-1032 (labelled "Lu formula [JRS]" in a comment),
      used only when `gbtype == 'GBT'` -- not peakfq's default, which hardcodes MGBT
      (`emafit.f:982`). The `n < 100` table of exact values is not sourced from this vendored
      Fortran at all (grepped for the literal constants, found nowhere in `vendor/`); it reads as
      the standard published Bulletin 17B table. Net effect on the decision: recomputing it is
      still the right call (it is what the native engine already does with this field, so the
      two engines report it the same, consistent way), but "already Fortran-verified" is not the
      reason -- the honest reason is "consistent with the native engine's own existing policy for
      an unreported diagnostic", which is a weaker but still sufficient justification. Documented
      in `flowfreq/fortran_engine.py::_frequency_results_from_reference`'s docstring.

- [x] `ReferenceResult` -> `FrequencyResults` adapter --
      `flowfreq/fortran_engine.py::_frequency_results_from_reference`. Follows the design doc §4
      table field by field. `ema_iterations`/`ema_converged` forced `None`; `n_censored`/
      `pilf_flows` derived from the builder's own `EmafitArrays`, never read off
      `ReferenceResult`; `skew_weighted`/`skew_regional` only populated when the caller actually
      supplied regional-skew information (station-only inputs make `emafitpr`'s own "weighted"
      column numerically equal to its at-site one, and reporting that as a weighted result would
      claim a weighting that never happened -- mirrors `ExpectedMomentsAlgorithm`'s own policy).
      The one field-semantics mismatch the design doc flagged (`ReferenceResult.n_systematic`
      counts every `dtype==0` row including gap-fill censoring; the native engine's own
      `n_systematic` counts only uncensored rows) is passed through directly, as instructed, and
      documented in the adapter's docstring rather than silently reinterpreted.

      **Verified two ways.** Unit-tested against a synthetic `ReferenceResult`/`EmafitArrays`
      pair (`tests/fortran_parity/test_fortran_adapter.py`, 9 tests) -- confirms the
      non-fabrication contract (`ema_iterations`/`ema_converged is None`), the skew-weighting
      policy in both branches, `pilf_flows` derivation, and that `n_censored` comes from the
      arrays not the reference. **Now also verified end to end against a live `emafitpr` call**:
      `Bulletin17C(...).run_analysis(method="ema", engine="fortran")` on Big Sandy (which has
      historic peaks, exercising the historic-only-dtype-1 rule through the whole stack, not
      just the builder) reproduces the golden's weighted skew (**-0.15634** live vs **-0.15631**
      golden -- the small residual is the same cross-machine EMA-fixed-point noise
      `test_live_vs_golden.py`'s own docstring already measures for this gfortran build, not an
      adapter defect), `n_censored == 37`, `n_peaks == 84`, and -- the actual non-fabrication
      check that matters, not just the synthetic-data one -- `ema_iterations is None` and
      `ema_converged is None` on real Fortran output too.

- [x] `engine=` on `Bulletin17C.run_analysis`, defaulting to `"native"` forever --
      `flowfreq/bulletin17c.py`. `engine="fortran"` requires `method="ema"` (raises `ValueError`
      naming why for `method="mom"`, checked before any Fortran import so it is testable without
      the extension built -- `tests/test_bulletin17c.py::TestEngineParameter`) and delegates to
      `fortran_engine.run_fortran_ema`, keeping the live `ReferenceResult`/`EmafitArrays` on the
      instance so `compute_quantiles`/`compute_confidence_limits` can re-invoke `emafitpr` at a
      *different* AEP list afterward -- confirmed necessary and implemented: Fortran quantiles
      come from `qP3sub`'s exact gamma quantile evaluated at the AEPs passed into `emafitpr`
      itself, not from applying a K-factor to already-fitted moments the way the native path's
      `compute_quantiles` can for any AEP after one fit, so a new AEP list needs a fresh Fortran
      call built from the same arrays (same fit, same MGBT decision, different probabilities).
      Verified live: requesting a 500-year quantile after fitting at the default AEP list
      returns a fresh, correctly-computed value (`Bulletin17C(...).compute_quantiles(aep=[1/500])`
      on Big Sandy through `engine="fortran"`, checked by hand against a direct `emafitpr` call
      at that AEP).

- [x] `compare_engines` + markdown output -- `flowfreq/workflow.py::compare_engines` /
      `EngineComparisonReport`. Reuses the existing `Bulletin17C.validate()` /
      `FrequencyComparator` machinery exactly as planned, rather than new comparison logic;
      recomputes the *native* side's quantiles/confidence limits at the requested `aeps`
      explicitly after fitting (`run_analysis` always fits at `STANDARD_AEP` internally, so a
      caller-supplied `aeps` would otherwise leave the two sides compared at different AEP sets,
      which `FrequencyComparator` would silently treat as "no keys in common" rather than raise
      -- found and fixed while implementing this, not a hypothetical). Imports
      `flowfreq.peakfqr` up front so the actionable `ImportError` (design doc §9 q4: no
      golden-file fallback, ever) raises before any native work is wasted, rather than partway
      through. `max_quantile_deviation_pct` reads `ComparisonResult.quantile_diffs` specifically,
      not the mixed `max_diff_pct` (which also folds in parameters/CIs).

      **Verified live on all four parity sites**, matching TODO.md's own P3 table and the design
      doc's own worked example almost to the decimal:

      | site | `max_quantile_deviation_pct` | overall | design doc / P3 table says |
      |---|---:|---|---|
      | Big Sandy 03606500 | **0.0587%** | PASS | "every quantile to <= 0.06%" |
      | Powder River 06326500 | **0.1003%** | PASS | "quantiles <= 0.10%" |
      | 12363000 | **0.1057%** | PASS | design doc §5's own example: "0.106 ... at the 500-yr" |
      | Cains Coulee 06327450 | **9.741%** | **FAIL** | P3 table: "quantiles 0.08% to 9.7%" |

      Cains Coulee's FAIL is the known, already-`xfail(strict=True)`'d `skew_weighted` residual
      (worst skew diff reported: **0.0580**, matching the documented 0.058) surfacing through
      this tool exactly as section 1 of the design doc says it should -- not a defect in
      `compare_engines` itself. Tested as such: `max_quantile_deviation_pct` is asserted under a
      measured bound (no xfail needed, that number is not in question), and the overall
      `.passed` is `xfail(strict=True)`'d with a reason pointing at the standing xfail in
      `test_wymt_vs_golden.py`, so the two flip together if that residual is ever resolved.
      Tests: `tests/fortran_parity/test_live_compare_engines.py`.

- [x] CLI `compare` subcommand -- `flowfreq/cli.py`. `--peaks` reads a CSV with `water_year`/
      `peak_flow_cfs` columns, the same shape `USGSgage.download_peak_flow` produces (so that
      DataFrame can be saved straight to CSV and used here); `--regional-skew`/`--regional-skew-se`/
      `--low-outlier-threshold`/`--tolerance-pct`/`--output` cover the common case. Historical
      peaks and perception thresholds are **not** exposed as CLI flags yet -- a record needing
      them goes through `compare_engines` directly; noted in the command's own `--help`, not
      silently unsupported. Exits nonzero when the comparison fails tolerance (verified: Cains
      Coulee's CSV run through the CLI exits nonzero and prints "FAIL"), and raises a
      `click.ClickException` wrapping the extension's own `ImportError` if it is not built,
      rather than a bare traceback. Verified live end to end via Click's `CliRunner`, including
      writing the markdown to `--output` and rejecting a CSV with the wrong column names.
      Tests: `tests/fortran_parity/test_live_cli_compare.py`.

**A live-Fortran testing hazard found and worked around, worth recording for whoever adds the
next test file here**: `test_fortran_oracles.py`'s `TestCainsCouleeAsGMseDiscrepancy` and
`TestSkewMseOracle` depend on being the *first* code in the whole pytest process to call any
`emafitpr`-family entry point -- documented in that file's own module docstring as a real,
already-known `mseg_all_sub` `SAVE`-state leak, not something introduced here. The first drafts
of this session's three new live-Fortran test files (`test_interval_builder_live.py`,
`test_compare_engines.py`, `test_cli_compare.py`) sorted alphabetically *before*
`test_fortran_oracles.py` and, by calling `emafitpr` on several different cases each, left that
file reading contaminated state (`mseg_all_sub` returning 0.2212 instead of 0.0749) even though
none of the new files call `mseg_all_sub` directly. Fixed by renaming all three with a
`test_live_` prefix, the same convention `test_live_vs_golden.py`/`test_native_vs_golden.py`/
`test_wymt_vs_golden.py` already use and which already sorts after `test_fortran_oracles.py`.
Confirmed fixed by running `pytest tests/fortran_parity/` (exactly `make parity`'s selection)
before and after the rename: 3 spurious failures before, the one known pre-existing 1e-12
platform-precision mismatch after (see below). **Rule for next time**: a new test file that
calls `emafitpr` on more than one case must sort after `test_fortran_oracles.py`, or must be
checked against the whole `tests/fortran_parity/` directory (not run alone) before trusting it.

**Full verification status, both with and without the extension**:

- `make clean-verify` with the extension genuinely absent (the `.pyd` and its four DLLs deleted
  by hand, since -- see the note above -- the `clean` target's own glob does not remove them on
  Windows): **595 passed, 7 skipped, 1 deselected, 5 xfailed**, lint and mypy clean. Same 5
  xfails as the pre-existing baseline; no new ones outside the live-Fortran files, which skip
  entirely here.
- `pytest tests/fortran_parity/` with the extension built (exactly `make parity`'s selection):
  **1 failed, 259 passed, 2 xfailed**. The one failure --
  `TestSkewMseOracle::test_reproduces_emafitpr_as_g_mse[big_sandy_03606500]`, off by ~1.1e-8
  relative at a `rel=1e-12` tolerance -- is the one the environment fix-up message flagged in
  advance as a near-certain MinGW-vs-Linux-gfortran platform artifact, not a correctness issue;
  left as is, not loosened, per instructions.
- A separate, **pre-existing** fragility (reproducible with none of this session's files added,
  confirmed by isolating it to only the original `tests/fortran_parity/*.py` plus
  `tests/validation/test_reference.py`) surfaces only in a *full local* `pytest tests/`
  run with the extension built: `tests/validation/test_reference.py::TestFromEmafit::
  test_live_call_matches_the_golden_file` fails by the same few-ulps-of-cross-machine-noise
  margin. Never reaches CI either way -- CI's `fortran` job runs only `tests/fortran_parity/`,
  and CI's default `test` job never builds the extension -- so it is recorded here rather than
  touched; not this item's file, not introduced by this item's changes.

All five design-doc pieces are done; nothing is left open under this item except the two
already-recorded, pre-existing platform/environment artifacts above (the 1e-12 gfortran
mismatch, and `make clean-verify`'s Windows `.pyd` cleanup gap) and possibly-worthwhile
follow-ups, not blockers:

- Historical peaks / perception thresholds are not exposed as CLI flags on `flowfreq compare`
  (noted in its `--help`); a caller needing them uses `compare_engines` directly. Small, if
  ever wanted -- the arguments already exist on the Python side.
- No binary-wheel distribution story, deliberately -- see below.

Binary wheels are deliberately **out of scope** for the first pass; source checkouts serve
the CLOMR/LOMR case this is for.

### Modules with no tests

`engine.py` and `report.py` were covered in 0.4.0. The remaining four are now covered too:

- [x] `hydrograph.py` (375 lines) -- `tests/test_hydrograph.py`, happy path + error case per public
      function.
- [x] `plots.py` (255) -- `tests/test_plots.py`.
- [x] `batch.py` (128) -- `tests/test_batch.py`, including a test pinning the 0.4.0 station-skew
      fix against `Bulletin17C`'s own unbiased estimator rather than a hardcoded number. Found (not
      fixed, out of this lane's scope) a real bug while writing these: `usgs.fetch_nwis_batch`
      returns plain dicts, but `batch.run_multi_site` reads `.flow` as an attribute, so every real
      multi-site call fails per-site with `'dict' object has no attribute 'flow'`, silently
      swallowed by a broad `except Exception` into `{"error": ...}`. Pinned as
      `tests/test_batch.py::TestAnalyzeSites::test_real_fetch_output_shape_is_analyzable`,
      `xfail(strict=True)`. The real fix belongs in `usgs.py` (should it return `PeakRecord`
      instead of a dict, or should `batch.py` adapt?) -- not done here.
- [x] `cli.py` (57, now larger with the `compare` subcommand) -- `tests/test_cli.py` covers
      `validate`/`benchmark` and `compare`'s extension-agnostic paths (CSV validation, the
      `ImportError`-to-`ClickException` mapping, exit codes) via a monkeypatched
      `compare_engines` so it runs the same with or without the built extension.
      `compare`'s real, Fortran-backed happy path is covered end-to-end in
      `tests/fortran_parity/test_live_cli_compare.py` instead (`requires_fortran`).

### Downstream

- [x] **Bump the app's pin once `v0.4.0` is tagged.** Done: `flowfreq-app/requirements.txt`
      now pins `flowfreq@v0.4.0` (local branch `bump-flowfreq-0.4.0` in that checkout, not
      yet merged to its `main` or pushed to origin -- see `flowfreq-app/TODO.md`).

      **The app's reported numbers do not change.** An earlier version of this entry said
      the bump was "a visible change to the deployed app"; that was wrong. `streamlit_app.py`
      goes through `workflow.run_ffa` to `Bulletin17C` and never constructs a `B17CEngine`,
      and the 0.4.0 station-skew fix is confined to `B17CEngine`. Checked against the whole
      `v0.3.0..v0.4.0` diff of `flowfreq/`: on the app's path, `freq_plot.py` is a rename
      plus an alias and docstrings, `workflow.py` is one docstring reference, and
      `bulletin17c.py` is a `TYPE_CHECKING` import block. No arithmetic moves.

      So the expectation when verifying the bump is **identical output**, and a difference
      is a defect rather than the change landing. What is worth smoke-testing is the
      `plot_frequency_curve` rename: the app still imports the old
      `plot_frequency_curve_streamlit`, which survives only as an alias. `make test` in the
      app repo catches that, since importing `streamlit_app.py` executes it in Streamlit's
      bare mode.

      `batch.batch_summary_table` and `plots` are the consumers that *do* move with the
      skew fix -- see the CHANGELOG for the measured deltas -- and neither is used by the
      app.

      Still open on this item: running that verification and smoke test against the bumped
      branch, and deleting the app's local `plot_peak_timeseries` copy -- the latter also
      needs the library-side move from the Follow-ups item below, which has not happened yet.

### P3 — The `var_mom` port, now complete

One `xfail(strict=True)` remains from this whole item: Cains Coulee's `skew_weighted` rung, in
`tests/fortran_parity/test_wymt_vs_golden.py`. The build fails the moment it starts passing.
Nothing else is blocked on it, and nothing else in this item is still open.

Everything here bottomed out in `var_mom` and its dependency tree, the one piece of the
reference implementation that had never been ported before this item started. What follows is
the full account, phase by phase, ending in the confidence-interval shape fix that closed it.

**The defect was censoring-specific, and that is measured rather than assumed.** Two
Wyoming/Montana parity cases were added for exactly this question, plus Big Sandy (which turns
out to have real censoring of its own -- 37 historical gap-year intervals -- that the parity
cases used for oracle-level testing don't construct, since they merge historical peaks into the
systematic record rather than configuring a real historical period):

| case | censoring | reference `Wd` | native `Wd` | native vs peakfq 8.1.0 |
|---|---|---:|---:|---|
| Powder River 06326500 | none | 1.0 | 1.0 | mean **0.0**, sd **3.7e-14**, at-site skew **4.5e-12**, weighted skew **7.5e-11**; quantiles ≤ 0.10% |
| Big Sandy 03606500 | 37 censored historical gap years | 1.0 (below HWN floor) | 1.0 | mean **3.5e-7%**, sd **3.3e-5%**, weighted skew **2.4e-6** (skew units); quantiles ≤ 0.06% |
| Cains Coulee 06327450 | 11 PILFs from MGBT | **0.184** | **0.186** | at-site skew **0.0002**, weighted skew 0.058 (skew units); quantiles 0.08% to 9.7%, 1.5% at Q100 |

Cains Coulee's remaining weighted-skew gap is now the *only* open numerical residual on the
censored path -- everything upstream of it (the at-site fit, `Wd`, ADJE, and -- confirmed by a
dedicated investigation once this looked like the last item, see the "Confidence-interval shape"
entry below -- `mse_ema`/`var_mom`/`mn2mvarb` themselves, called at this site's own real,
sensitive input) matches peakfq 8.1.0 closely. It does **not** trace to `mn2mvarb`'s own
numerical-differentiation gap the way earlier revisions of this document assumed; that
explanation did not survive being checked directly. Big Sandy's own weighted skew, by contrast,
is now correct to 2.4e-6 -- its at-site skew (0.0066) sits under the 0.04 HWN floor, so it never
exercises `detrat`/ADJE's nontrivial path the way Cains Coulee's -0.708 does.

So with nothing censored the native EMA reproduces peakfq to machine precision — the in-loop
regional-skew weighting is *right*, not merely closer — and with heavy censoring but no PILFs
(Big Sandy) it now also reproduces to machine precision. Cains Coulee, the one case with both
PILFs and an at-site skew above the HWN floor, is where the one remaining residual lives, and
it is the oracle for `detrat`.

- [x] **Port `var_mom` and the routines it needs.** ~1,100 lines of Fortran in `emafit.f`
      alone, plus `CHOL33` (`probfun.f`), `DLGINV` (`imslfake.f`), and `expmomderiv`, `m2mn`,
      `m2p`, `fp_tnc_icdf` which live outside `emafit.f`. The tree:

      ```
      mseg_all(ADJE)  -> mse_ema -> var_mom -> mP3(k=6), pP3, varb, varc,
                                               d_est -> m2p, qP3, expmomderiv
                                    m2mn
                                    mn2mvarb -> mn2m_var, mc2mnvb, jmc2mnvb,
                                                chol33, dlginv   (iterative solver)
      VAR_EMAB        -> regmoms, gridmake (the inverse-Cholesky quadrature),
                         ci_ema_m3b
      ```

      Sizes, for planning: `var_mom` 111, `mP3` 141, `regmoms` 111, `gridmake` 93,
      `mn2mvarb` 77, `qP3` 73, `VAR_EMAB` 66, `jmc2mnvb` 61, `mc2mnvb` 58, `mseg_all` 57,
      `pP3` 54, `ci_ema_m3b` 51, `mse_ema` 49, `d_est` 39, `varb`/`varc` 23 each.

      **Verify per routine, not end to end.** `build_fortran/_emafort.pyf` now exposes
      `mseg_all_sub`, `detratsub`, `var_mom` and `moms_p3` alongside `emafitpr`, and
      `tests/fortran_parity/test_fortran_oracles.py` pins what each says. Checking a ported
      routine only through `emafitpr` means checking it through a fixed point with a
      condition number around 1e13, where a correct routine can look wrong and a wrong one
      can look right. All four were already externally callable — peakfqr's own R code calls
      them via `.Fortran()`, and `vendor/peakfqr/R/fortranWrappers.R` documents the
      conventions — so no new Fortran was needed.

      Two usage notes that cost time to find:

      * `detrat` takes the **post-MGBT** thresholds. Cains Coulee's input thresholds are
        uncensored throughout, so calling it with those returns 1.0; MGBT raises the lower
        threshold to log10(332) = 2.521 inside the fit, and only then does it give the 0.184
        `emafitpr` reports. Those arrays come back as `tlema`/`tuema`.
      * Group thresholds on **exact** values. Rounding to 12 decimals to be "robust" moves
        Big Sandy's at-site skew MSE by 2.2e-4.

      What the oracles found immediately, both now covered by tests:

      * `_ema_iteration` reproduces `moms_p3` **exactly** on uncensored rows (0.0 on the
        mean, ~1e-14 variance, ~1e-12 skew) and diverges only where intervals are censored
        (Cains Coulee: 0.70% variance, 4.94% skew). So the transcribed formulas are right
        and the residual is in the truncated-P3 moment code for censored intervals — that is
        what the port has to replace, and it is a smaller target than "the EMA is off".
      * `_b17b_skew_mse` is exact against `mseg()` up to n = 150 and **31% high at n = 200**
        (0.0479 against 0.0365). `mseg_all` evaluates `mseg()` at `min(n, 150)` then lets
        ADJE's bias adjustment partially undo the cap; flowfreq applies the cap alone. On a
        record longer than 150 years that over-weights the regional skew. No parity case is
        that long, so only the oracle test detects it, and fixing it needs `mse_ema`.

      **Phase 1 done: the leaf layer.** `flowfreq/_p3_moments.py` now carries `m2p`, `m2mn`,
      `mn2m`, `p_p3`, `q_p3` and `m_p3` — everything `var_mom` calls directly except `varb`,
      `varc` and `d_est` (Phase 2). Each is a direct transcription of `emafit.f`/`probfun.f`,
      checked routine-by-routine against six new oracles `_emafort.pyf` exposes (`m2p`,
      `m2mn`, `mn2m`, `pp3`, `qp3sub`, `mp3` — lower-cased there because gfortran's symbol
      table is all-lowercase and f2py's generated wrapper does not re-case a mixed-case name
      in the `.pyf`, which cost a build before it was caught). Tests:
      `tests/fortran_parity/test_fortran_oracles.py` (Fortran-gated) and
      `tests/test_p3_moments.py` (pure Python, no build needed). Nothing here is wired into
      `ExpectedMomentsAlgorithm`/`Bulletin17C` yet — Phase 2 is composing `varb`/`varc`/
      `d_est`/`var_mom` on top of this layer and only then deciding how (or whether) it
      replaces `_truncated_gamma_moment`/`_truncated_normal_moments` in `bulletin17c.py`.

      **A real defect in the reference, found while verifying `mP3`, not in this port.**
      `mP3` blends an incomplete-gamma solution with a Wilson-Hilferty one; the 2024 upstream
      patch (`FP_G1_CDF`/`FP_G1_MOM_TRC`, `probfun.f`) promoted the incomplete-gamma call
      itself to real128 for large `alpha = 4/skew**2`, but not the surrounding
      `choose(i,j)*tau**(i-j)*fp(j)` expansion `mP3` builds on top of it (`emafit.f:3049`).
      On Big Sandy's own censoring group (at-site skew 0.0066, so alpha ≈ 9e4), that
      expansion cancels by roughly 11 orders of magnitude at k = 6 and the Fortran result
      goes **negative** — impossible for `E[X**6]` of a positive-truncated variate.
      `flowfreq._p3_moments.m_p3` keeps the whole expansion in `mpmath` at 50 decimal digits
      (`_GAMMA_MOMENT_DPS`) rather than only the CDF call, and its k = 4..6 values are
      confirmed against independent arbitrary-precision quadrature of the truncated density
      (not against the Fortran, which is the thing in question) —
      `test_fortran_itself_loses_precision_on_big_sandy_s_censored_group`.

      **Resolved by Phase 2, as it turns out: `var_mom` itself never reaches this regime.**
      `var_mom` clamps `|skew|` to `[0.0632, 1.41]` before ever calling `mP3`/`pP3`
      internally (`emafit.f:2324`, easy to miss reading linearly) — so even Big Sandy's raw
      0.0066 becomes 0.0632 inside `var_mom`, capping alpha at ~1000 rather than the ~9e4 the
      finding above used directly. `flowfreq._var_mom.var_mom` matches the Fortran oracle to
      1e-9 relative on Powder River and Cains Coulee and 1e-3 on Big Sandy (see Phase 2 below
      for where that residual is), not the orders-of-magnitude gap `mP3` alone showed. The
      finding above stands as a real, documented Fortran defect reachable by calling `mP3`
      directly (as the oracle tests do) — it just is not reachable *through* `var_mom`.

      New dependency: `mpmath` (pure Python, no compiled extension), added for exactly this.

      **Phase 2 done: `varb`, `varc`, `d_est`, `expmomderiv` and the `var_mom` composition
      itself.** `flowfreq/_var_mom.py`. `expmomderiv` differentiates `_p3_moments`'s own
      (already Fortran- and quadrature-verified) truncated-moment function numerically via
      `mpmath.diff` rather than transcribing `DEXPECT`'s closed-form chain, which needs
      `DDGAM` — a nontrivial derivative of the incomplete gamma function w.r.t. its shape
      parameter. Checking the same derivative at 50/80/120 `mpmath` working digits agrees to
      30+ stable digits, so the ~1e-5 relative gap against the Fortran on Big Sandy's one
      real censored group is on the Fortran side (`DDGAM`'s own series truncates at
      `TOL=1e-11` per term, `probfun.f:1054`), not evidence the numerical-differentiation
      shortcut is wrong — though that is not independently proven the way the `mP3` finding
      is. `DPDM` and the matrix bookkeeping around all four (`DSET`/`DMSUM`/`DMRRRR`/
      `DMXYTF`/`DMMULT`/`DLGINV` in the Fortran) are plain closed-form algebra and linear
      algebra respectively, done directly with `numpy`/`numpy.linalg` rather than transcribed
      routine by routine — there is no numerical subtlety in a 3x3 sum, product, or inverse
      the way there is in the incomplete-gamma work.

      Six more oracles in `_emafort.pyf` (`varb`, `varc`, `d_est`, `expmomderiv`, alongside
      Phase 1's four): `varb`/`varc` match to 1e-9, `d_est` to 1e-3 (inherits `expmomderiv`'s
      gap), and the full `var_mom` composition to 1e-9 on the two cases that never exercise
      `d_est`'s nonzero path (both tail probabilities stay under its 1e-6 cutoff) and 1e-3 on
      Big Sandy, the one case that does. Tests: `tests/fortran_parity/test_fortran_oracles.py`
      (Fortran-gated) and `tests/test_var_mom.py` (pure Python, including an independent
      cross-check against the classical delta-method moment covariance for the fully
      uncensored case, computed without going through `var_mom`'s threshold-group machinery
      at all). Still nothing wired into `ExpectedMomentsAlgorithm`/`Bulletin17C` at this
      point — see the Phase 3 note below the "Skew weighting" item for `mn2mvarb`/`mse_ema`,
      which is now done, and what is still open (`VAR_EMAB`/`regmoms`/`ci_ema_m3b`, the
      CI-shape formula, plus the wiring decision for both).

- [x] **Skew weighting — wired in, `detrat` included; what's left is upstream of weighting.**
      The structural half is done
      (see below): the regional skew is now folded into the EMA fixed point as `moms_p3`
      does it, which took Big Sandy from 35% to **24.1%** (−0.1187 against peakfq's −0.1563)
      and improved the mean and variance at the same time.

      All of the remainder is one input. peakfq's default `at_site_option` is `ADJE`
      (`emafit.f:3888`): `as_G_mse = bias_adj * mseg(min(n,150), G)`, where `bias_adj` is the
      censoring inflation `mse_ema(censored)/mse_ema(uncensored)` — **1.4844** on Big Sandy.
      Only `mseg()` is implemented (`_b17b_skew_mse`), so a censored record under-weights the
      regional skew. Feeding peakfq's own `as_G_mse` (0.09437) through this code gives
      **−0.1592 against −0.1563, a 1.9% gap**, which is inside the xfail's 0.02 bound. So the
      structure is right and the input is not.

      Also then-unimplemented: `detrat`, the Halloween determinant ratio — now done (below).
      `emafit.f:763` applies it only when the at-site skew is ≥ 0.04 in magnitude, and Big
      Sandy's is 0.0066, so `Wd` is 1 there either way. **Cains Coulee 06327450 covers it**:
      its at-site skew is −0.708 (peakfq's; flowfreq's own is −0.830) and its reference `Wd`
      is **0.184**, so flowfreq's old implicit 1 over-weighted the regional skew by more than
      fivefold on that site. That case was the acceptance test for `detrat`; see below for
      where the resulting `xfail(strict=True)` assertions landed — not flipped, for a reason
      unrelated to `detrat` itself.

      **Phase 3 done: `mn2mvarb`/`mse_ema` — `as_G_mse` is now computable, in Python, matching
      the Fortran.** `flowfreq/_mse_ema.py`. `mse_ema(nobs, tl, tu, mc, kmom) →`
      `var_mom → m2mn → mn2mvarb`, diagonal element `kmom`; feeding it through both a censored
      and an equivalent uncensored call reproduces the documented numbers on Big Sandy exactly:
      `bias_adj` **1.4843** (documented 1.4844) and `as_G_mse` **0.094366** (documented 0.09437,
      and matching the `mseg_all_sub` oracle's own `as_G_mse_o`, **0.094375**, to 0.01%). The
      skew-weighting gap this closes is therefore no longer an open question — the input peakfq
      uses is now reproducible — only wiring it into the fixed point remains (see below).

      `mn2mvarb` (`emafit.f:2514`) solves an inverse problem: find the central-moment covariance
      whose forward map through `mc2mnvb` (an eight-point "Inverse Modified Cholesky Gaussian
      Quadrature" — `gridmake`, `normquad`/`gammaquad` via `numpy.polynomial.hermite.hermgauss`/
      `scipy.special.roots_genlaguerre`, ported faithfully) reproduces a given noncentral-moment
      covariance. The Fortran solves this with a bespoke step-halved Newton iteration, checking
      positive-semi-definiteness via `chol33` at every trial step. Reparametrizing the unknown as
      `chol33`'s own six free entries — so the candidate covariance is `V.T @ V`, automatically
      PSD for *any* real `V`, no guard needed — turns it into unconstrained root-finding
      (`scipy.optimize.root`), started from the same linearized estimate (`mn2m_var`) the Fortran
      uses. Both land on the same root: `mn2mvarb` matches the Fortran to 1e-6 relative on
      Powder River/Cains Coulee, 1e-3 on Big Sandy (inherits `expmomderiv`'s gap — see Phase 2).

      **Performance mattered here in a way it had not before**, because `mse_ema` sits on what
      would be the fixed point's hot path: the first attempt (`expmomderiv`'s Jacobian via
      9 separate `mpmath.diff` calls, one per moment/parameter pair) took **1.14 s** for one
      `mse_ema(censored)` call on Big Sandy — untenable, since `Bulletin17C.run_analysis` would
      need two such calls (censored and uncensored) per analysis. Profiling
      (`cProfile`) pointed at the redundancy directly: each `mpmath.diff` re-evaluated the
      truncated-moment function from scratch, recomputing the incomplete-gamma "down" term
      (which does not depend on which moment k is being computed) up to 3 times per call. Batching
      the three moments into one evaluation (`_gamma_trunc_moments`/`_fp_g1_mom_trc_batch` in
      `flowfreq/_p3_moments.py`, sharing `down` across k) and replacing the 9 `mpmath.diff` calls
      with 7 manual, batched central differences (2 evaluations per parameter plus 1 for the
      center point, instead of ~2 per moment-parameter pair) cut it to **0.32 s** — a ~3.6x
      speedup, bit-identical to the un-optimized version once the finite-difference step was
      tightened to `1e-12` relative (`_DEXPECT_STEP`; the first attempt at `1e-6` was too coarse
      and cost three tests ~0.1-1% accuracy against the Fortran, since found by re-running the
      full suite — measured, not assumed). Still not fast enough to call on every EMA
      iteration, but `_regional_skew_equivalent_years` (`bulletin17c.py:903`) only needs it once
      per `run_analysis()` call, using the converged at-site skew — not once per iteration — so
      this is viable for the wiring below, at roughly 0.3-0.4 s added per analysis with a
      regional skew supplied.

      **A real, serious defect found in the reference while testing this, not in the port**:
      `mseg_all_sub` is not safe to call more than once per process. Confirmed outside pytest
      entirely — call it once (correct), call `emafitpr` for a *different, unrelated* case, call
      `mseg_all_sub` again with the *original* (unchanged) inputs: the second call silently
      returns the *uncensored* value instead, and stays wrong for every subsequent call. The
      arrays going in are unmutated; this is Fortran-side `SAVE`d state leaking across calls to
      different entry points in `emafit.f`. Documented in
      `tests/fortran_parity/test_fortran_oracles.py`'s module docstring (found while adding the
      Phase 3 oracles, which is why that test file's own assertions are careful never to call
      `mseg_all_sub` a second time in the same process). Not fixed here — `vendor/` is not
      edited — but worth knowing before calling `flowfreq.peakfqr`/`flowfreq.validation.reference`
      for more than one site in a single long-running process (a Streamlit session, a batch
      script): `emafitpr` itself appears unaffected (the existing multi-case parity tests already
      interleave it across cases and pass), but anything downstream that calls `mseg_all`'s ADJE
      path a second time should not be trusted without a fresh process.

      **Done: wired into `Bulletin17C`, and the parity xfail actually flips.**
      `ExpectedMomentsAlgorithm._perception_threshold_groups` builds the `(nobs, tl, tu)`
      groups from `self._intervals` (the `perception_threshold > 0` rule above, verified
      against `tests/fortran_parity/cases.py::build_emafit_inputs`); `_adje_skew_mse` calls
      `mse_ema` through a `@staticmethod`/`@lru_cache` method (`_adje_bias_adjustment`, same
      pattern as `_mgbt_pvalue` and for the same reason — repeated fits of the same fixture
      are common, and one call is not free) and feeds the result into
      `_regional_skew_equivalent_years`, which now takes `mean_log`/`std_log` too. Falls back
      to `bias_adj = 1` (today's behavior) with a logged warning if `mse_ema` raises, rather
      than let an ancillary correction fail the whole analysis — `mn2mvarb`'s root-find is not
      guaranteed to converge on every input.

      Confirmed end to end, not assumed: Big Sandy's `TestRung3Moments::test_weighted_skew`
      in `tests/fortran_parity/test_native_vs_golden.py` — peakfq's −0.1563 against a 0.02
      tolerance — now passes, so its `xfail(strict=True)` is removed (that test file's own
      alarm going off is what caught it). Measured weighted-skew gap against the reference,
      `tests/integration/test_hybrid_workflow.py`: **0.0026** in skew units, down from 0.0376
      — matching the ~1.9% figure predicted above almost exactly, with the small residual
      being `mn2mvarb`'s own ~1e-3 relative gap on Big Sandy (Phase 2's `expmomderiv` note),
      not `detrat`, which does not apply here (Big Sandy's at-site skew is under the 0.04 HWN
      floor). Three quantiles in `tests/validation/test_big_sandy.py` (AEP 0.002, 0.99, 0.995)
      moved 2.2–2.8% *further* from the 2012 PeakfqSA manual as a direct, expected consequence
      — that manual predates HWN/ADJE and was already documented as not reproducible by peakfq
      8.1.0 — so they are now `xfail(strict=True)` there instead of silently passing at a
      widened tolerance, with the actual `PEAKFQ_810_*`/`tests/fortran_parity/` parity checks
      (the ones that matter) unaffected.

      **Cost**: the full suite went from ~29 s to ~40 s. `_adje_bias_adjustment` costs
      ~0.3–0.4 s the first time a given (fixture, moments) combination is fit with a regional
      skew supplied; the `@lru_cache` absorbs repeats. Acceptable for now; worth another pass
      if it grows further.

      **`detrat` done too — ported, wired, and verified.** `flowfreq/_detrat.py`, a direct
      transcription of `emafit.f`'s `detrat` (the 3x3-vs-2x2 determinant ratio) and
      `probfun.f`'s `EXPMOMCDERIV` (the censored-region expected-moment Jacobian it needs,
      reusing `flowfreq._var_mom._dexpect` for the open-tail pieces rather than a second
      truncated-gamma implementation). Verified against two oracles `_emafort.pyf` exposes
      (`expmomcderiv`, and Phase 1's existing `detratsub`): matched Cains Coulee's real,
      post-MGBT `Wd` of 0.184 to ~1e-9 relative precision on the first attempt — noticeably
      tighter than `mse_ema`'s ~1e-3 on Big Sandy, since `detrat` needs no Newton-solve/
      quadrature machinery. Wired into `_regional_skew_equivalent_years` via a new
      `@staticmethod`/`@lru_cache` method, `_detrat_wd`, mirroring `_adje_bias_adjustment`'s
      pattern; falls back to `Wd = 1` with a logged warning if it raises, the same posture
      `_adje_skew_mse` already takes on `mse_ema`. Tests: `TestExpMomCDerivPort`/
      `TestDetratPort` in `tests/fortran_parity/test_fortran_oracles.py` (Fortran-gated),
      `tests/test_detrat.py` (pure Python, no build needed).

      **Wiring it surfaced a real gap in `_perception_threshold_groups`, not in `detrat`
      itself.** That method reconstructs the `(nobs, tl, tu)` groups `mse_ema`/`detrat` need
      from `self._intervals`' `perception_threshold` field — but MGBT-censored PILFs get
      `perception_threshold = 0.0` in flowfreq's model (a censored *value*, correctly, for the
      moment iteration itself), while peakfq's own `tlema`/`tuema` (confirmed via a direct
      `emafitpr` call on Cains Coulee) raise the perception threshold for the **entire**
      systematic record to the MGBT cutoff once one is found, not just the flagged PILF years.
      Left unfixed, `detrat` would only ever see Cains Coulee as one uncensored group and
      return `Wd = 1` regardless of how correct the port was — the same wrong answer as before,
      for a different reason. Fixed by having `_perception_threshold_groups` take
      `max(interval.perception_threshold, self._ema_params.low_outlier_threshold)` for every
      *systematic* interval (`self._ema_params.low_outlier_threshold` is `0.0` when MGBT finds
      no PILFs, verified by reading `_multiple_grubbs_beck`'s early-return branches, so this is
      a no-op — `max(x, 0.0) == x` — on every record without one); historical-period intervals
      are left alone, since MGBT runs on systematic peaks only and their own threshold is an
      unrelated restriction. This is what actually moved Cains Coulee's native `Wd` from 1.0 to
      0.174.

      **The result is not a clean win, and that is the honest reading of it, not a regression.**
      `Wd` and `bias_adj` are now correct (or close to it — see the table above); but
      `skew_weighted` is a blend of the (now-correct) regional weight and `skew_station`, which
      is still 0.122 off from a separate, deeper defect: flowfreq's at-site EMA moment
      iteration (`_ema_iteration`/`_compute_ema_moments`) never used the Fortran-verified
      truncated-moment code, only its own approximation, and got the bias-correction sample
      size wrong on top of that. With the old buggy `Wd = 1`, more weight landed on the
      regional skew, which happened to be closer to peakfq's answer here than the broken
      at-site fit was — so the bug was accidentally diluting a different error, not fixing
      anything. Once `Wd` dropped to its (then) correct ~0.17, that dilution shrank and
      `skew_weighted`'s own gap against peakfq grew from 0.098 to 0.172, which was the state of
      things until the fix below. **That fix is now done** — see "At-site EMA moment iteration"
      immediately below — and it is what actually flips this case's `skew_at_site` xfail; fixing
      `detrat` alone, as this bullet originally suspected, could not have.

- [x] **At-site EMA moment iteration — the actual dominant defect, now fixed.**
      `_ema_iteration`/`_compute_ema_moments` (`bulletin17c.py`) is flowfreq's own transcription
      of `moms_p3` (`emafit.f:1344`), verified against a `moms_p3` Fortran oracle
      (`tests/fortran_parity/test_fortran_oracles.py::TestMomentIterationOracle`) since Phase 1.
      That oracle test showed the transcription was **exact** on uncensored rows (0.0 mean,
      ~1e-14 variance, ~1e-12 skew) and diverged only where intervals were censored (Cains
      Coulee: 0.70% variance, 4.94% skew) — correctly pointing at the censored-interval code,
      not the surrounding formulas. Fixing it took two changes, not one:

      1. **The truncated-moment formula itself.** `_compute_ema_moments`'s censored branch had
         its own approximate truncated-gamma/truncated-normal moment code (`scipy`'s `gammainc`/
         `gammaincc`, standardized bounds, a Wilson-Hilferty blend by hand) predating this
         session's `var_mom` port — but Phase 1 had *already* ported and Fortran-verified the
         real thing, `flowfreq._p3_moments.m_p3` (`mP3`, `emafit.f:2983`), for `var_mom`'s own
         use, and never wired it back into the E-step that inspired the port in the first place.
         Swapped in directly: `m_p3(tl, tu, [mean, var, skew], 3)` returns `E[X^k]` for
         `k=1..3` in real (unstandardized) log10-flow space, replacing ~110 lines of
         standardization/branching with one call. Grouped by distinct `(lower, upper)` pairs
         across intervals — Cains Coulee's 11 PILFs all share one MGBT cutoff, so this is one
         `m_p3` call per iteration, not eleven.
      2. **The bias-correction sample size.** `moms_p3`'s closing lines apply correction
         factors `c2`/`c3` sized by `n_bcf`, whose value depends on a Fortran flag (`bcf`,
         common block `/tac002/`) with two branches: `bcf=1997` (Cohn et al.) uses `n_bcf = n_e`
         (the **exact-peak count**), `bcf=2004` (Griffis et al.) uses `n_bcf = n` (the **total**
         interval count). `emafit.f:3898` sets the vendored default to `1997`; the `2004` line
         right below it (`emafit.f:3899`) is commented out and has been since before this
         repository forked from upstream. flowfreq's `_ema_iteration` used `n` (the `2004`
         convention) unconditionally — silently the wrong default, on every fit with any
         censored interval, since before this session. Fixed by computing `c2`/`c3` from
         `sums.n_exact` instead.

      Confirmed both were needed and together sufficient: applying only fix 1 left the same
      ~0.7%/4.9% gap essentially unchanged (traced by hand, not assumed — see the git history
      for the intermediate measurement); applying both together closed
      `TestMomentIterationOracle`'s Cains Coulee case from that gap to **1.6e-10 / 3.0e-10**
      relative (mean, var, skew) — the same level the uncensored cases already had. Renamed that
      test from `test_censored_rows_are_where_it_diverges` to `test_matches_on_censored_rows_too`
      to match.

      **End-to-end impact, measured against the real peakfq 8.1.0 goldens** (not the `moms_p3`
      oracle alone): both bugs bit **Big Sandy**, not just Cains Coulee — its 37 censored
      historical gap-year intervals (missing years within the historical perception period)
      exercise the exact same code path, even though Big Sandy has no MGBT PILFs at all. See
      the Status section and the P3 table above for the before/after numbers on both sites;
      `tests/integration/test_hybrid_workflow.py::test_weighted_skew_gap_is_closed`,
      `tests/validation/test_big_sandy.py` and `tests/fortran_parity/test_wymt_vs_golden.py`
      were all updated to the new measurements, and two more `xfail(strict=True)` assertions
      flipped and were un-xfailed (Cains Coulee's `skew_at_site`, and Big Sandy's AEP-0.002
      quantile against the 2012 manual) — both real, both confirmed by rerunning against their
      respective golden references before removing the marker, not assumed from the mechanism
      alone.

      **Cost**: the full suite went from ~42 s to ~76 s. `m_p3` is `mpmath`-backed (50 decimal
      digits when the incomplete-gamma branch is live), and unlike `_adje_bias_adjustment`/
      `_detrat_wd` it cannot be `@lru_cache`d the same way — it runs inside the fixed-point
      iteration itself, with different `(mean, std, skew)` on every call. Grouping by distinct
      censoring bounds keeps it to one call per group per iteration rather than one per
      interval, which is what makes this tractable at all (Cains Coulee: 11 intervals, 1 group).
      Not revisited further this session; worth a look if the suite grows past this budget.

- [x] **Confidence-interval shape — done.** `compute_confidence_limits()` used to form
      `log_Q ± z·se`, symmetric by construction (ratio 1.000 at every AEP), where peakfq skews
      right with return period — 1.03 → 1.31 → 1.41 at AEP 0.1 / 0.02 / 0.01 on Big Sandy.

      `ci_ema_m3b` (`emafit.f:1853`) itself is short:

      ```
      beta1     = cov(yp, syp) / var(yp)            # regression of the s.e. on the quantile
      var_xsi_d = var(syp) - cov(yp, syp)**2/var(yp)
      nu        = max(5, 0.5 * var(yp)/var_xsi_d)   # Student t degrees of freedom
      t         = t_nu((1 + eps)/2)
      ci_high   = yp + sqrt(var(yp)) *  t / max(0.5, 1 - beta1*t)
      t         = -t
      ci_low    = yp + sqrt(var(yp)) *  t / max(0.5, 1 - beta1*t)   # same denominator formula
      ```

      Both lines use the *same* `1 - beta1*t` denominator formula -- the whole asymmetry comes
      from `t` being reassigned to `-t` before the second line, not from a different formula for
      the two sides. Easy to mistranscribe as `1 + beta1*t` for `ci_low` (this session did,
      once, before checking the exact source text): that version gives a plausible-looking but
      wrong answer -- point estimates matched the Fortran exactly, and the interval was still
      asymmetric, just by the wrong amount, which is a harder bug to notice than an outright
      crash.

      What made this a small port after all: `beta1` needs `cov(yp, syp)` from `VAR_EMAB`
      (`emafit.f:1972`), and `VAR_EMAB` needs `regmoms` (`emafit.f:2173`) and `GRIDMAKE`
      (`emafit.f:2039`) -- but `regmoms` is `var_mom` (Phase 2) → `m2mn` (Phase 1) → `mn2mvarb`
      (Phase 3) plus regional-info blending arithmetic, and `GRIDMAKE` is exactly
      `flowfreq._mse_ema`'s existing `_gridmake`/`_covw` (already Fortran-verified indirectly,
      through `mc2mnvb`, which is `GRIDMAKE + M2MN + COVW` composed). The only routine with no
      existing counterpart was `VAR_EMAB` itself -- a nested quadrature (one 8-point grid around
      the fit, then a fresh `regmoms`/grid pair *at each of those 8 points*, to capture how the
      quantile and its own standard error co-vary) -- and `ci_ema_m3b`, the short formula above.
      New module: `flowfreq/_var_emab.py`.

      **Call convention, worth recording since it cost real time to pin down**: `VAR_EMAB`'s own
      probability argument is *non-exceedance* probability (`q_p3` uses `ndtri(q)` directly), so
      callers pass `pq = 1 - aep`, not `aep` itself -- passing `aep` directly gives a plausible
      but backwards-ordered `yp` array. And `regmoms`/`VAR_EMAB`'s two Fortran signatures use
      *different* argument orders for `(r_G_mse, r_M_mse, r_S2, r_S2_mse)`; mixing them up (this
      session did, once) silently sends the real regional-skew MSE into the wrong slot and
      produces a materially narrower, less-asymmetric interval that still looks plausible enough
      to not obviously be wrong.

      **A fourth instance of the `SAVE`d-state leak already documented for `mseg_all_sub`**: a
      direct call to the raw `var_emab`/`regmoms` Fortran oracle drifts ~2e-3 relative once other
      tests earlier in the same process have exercised `emafitpr`/MGBT/`detrat`, even though it
      matches exactly when called first in a clean process. `TestVarEmabPort` (`tests/fortran_parity
      /test_fortran_oracles.py`) checks against the committed golden file only, never a live
      oracle call, for exactly this reason.

      **Verified end to end, not just at the oracle level.** Big Sandy's own confidence bounds
      now match peakfq 8.1.0 within **0.06%** at every AEP tested (was symmetric by construction);
      the asymmetry ratio itself matches within 0.0007–0.0022 (was off by as much as 0.4 in ratio
      terms at the tail). Two `xfail(strict=True)` rungs in
      `tests/fortran_parity/test_native_vs_golden.py::TestRung6ConfidenceIntervals` flipped and
      were un-xfailed; a third assertion pinning the old symmetric behavior was rewritten into a
      positive check of the new one (`test_native_bounds_match_peakfq_closely`). Cains Coulee
      inherits its usual larger residual here too (`skew_weighted`'s own 0.058-skew-unit gap,
      the one item this whole port still leaves open — see the table above).

      **Cost, measured**: full suite ~76 s → ~93 s. `var_emab` is nine `regmoms` calls per
      confidence-limit computation (one outer + eight inner grid points, each a full
      `var_mom`/`mn2mvarb` solve) -- the single most expensive piece in the whole `var_mom` port,
      more than `mse_ema`'s own ~0.3–0.4 s. `ExpectedMomentsAlgorithm._cohn_confidence_bounds`
      is `@lru_cache`d the same way `_adje_bias_adjustment`/`_detrat_wd` are, so repeated fits of
      the same fixture (common across this test suite) pay it once. Falls back to the base
      class's symmetric formula, logged, if `var_emab` raises for any reason -- `MethodOfMoments`
      is untouched (it has no `_perception_threshold_groups`, and this whole item was always
      about the EMA path).

      **Not done, separately**: the pseudo effective record length (`as_G_PRL_o`, 54.373 for Big
      Sandy) is `eff_n * as_G_mse_Syst / as_G_mse` (`emafit.f:758`) -- a diagnostic value peakfq
      reports that flowfreq does not currently surface anywhere in `FrequencyResults`. Not needed
      for the confidence bounds themselves (`VAR_EMAB` never asks for it), so it was out of scope
      here; a small, separate follow-up if that diagnostic is ever wanted.

- [x] **Re-investigated Cains Coulee's `skew_weighted` residual — the earlier explanation was
      wrong, and now there's a much more precise, evidence-backed one.** With every other P3
      item closed, this was the one thing left labeled "small, understood" -- Phase 2 had
      documented it as `mn2mvarb`'s own numerical-differentiation gap in `expmomderiv`, "not
      independently proven the way the `mP3` finding is." Checking that claim directly, rather
      than continuing to repeat it, was the actual completion of P3.

      **The claim did not survive being checked.** `flowfreq._mse_ema.mse_ema(kmom=3)`, called
      standalone with Cains Coulee's real post-MGBT censoring group (`nobs=32, tl=2.521, tu=20`,
      its 332 cfs MGBT cutoff) and its real at-site fit, matches the Fortran oracle to **3e-8**
      relative -- nowhere near a ~1e-3 precision limit large enough to explain a 0.058-skew-unit
      gap. `var_mom`/`mn2mvarb` are not the bottleneck at Cains Coulee's own sensitive input
      after all.

      **What is actually happening**: `emafitpr`'s own internally-computed, reported `as_G_mse`
      for Cains Coulee -- 0.2212, committed in the golden file as `skew.as_G_mse_o`, and what the
      golden `skew_weighted` (-0.604) was built from -- does not match what calling the *same*
      `mseg_all` Fortran routine standalone gives for the *identical* `(nobs, tl, tu, mc)`: 0.0749,
      a ~3x difference. That 0.0749 is exactly what `flowfreq.bulletin17c.ExpectedMomentsAlgorithm
      ._adje_bias_adjustment` computes too (it is the same formula) -- so flowfreq's native fit
      is *internally consistent* with a clean, from-first-principles composition of independently
      Fortran-verified routines; it is `emafitpr`'s own reported value that disagrees with that
      composition, not the other way around.

      **Ruled out, in order, before concluding the mechanism itself is unresolved**:
      * *Existing test coverage passing was a false signal, not confirmation.*
        `TestSkewMseOracle.test_reproduces_emafitpr_as_g_mse[cains_coulee_06327450]` (added in an
        earlier phase) "passes," but it calls `mseg_all_sub` with `golden["inputs"]`'s tl/tu --
        which, as already documented for `detrat`, are *uncensored* for this site (MGBT creates
        the real censoring inside the fit). ADJE's bias adjustment is a no-op on an uncensored
        group, so that call reduces to the plain B17B `mseg()` value, which happens to equal
        `as_G_mse_o` (0.2212) -- meaning that test never actually exercised ADJE's bias
        adjustment for Cains Coulee at all. It was quietly testing the wrong group the entire
        time.
      * *Cross-case `SAVE`d-state contamination* (the already-documented `mseg_all_sub` bug,
        findings 1/3/4 in `test_fortran_oracles.py`'s module docstring) -- ruled out by
        regenerating Cains Coulee's golden file in total isolation
        (`python tools/gen_fortran_golden.py cains_coulee_06327450`, the only case in that
        process): still 0.2212. Whatever this is, it does not require an intervening,
        *different* case's `emafitpr` call the way the other four findings do.
      * *`momsadj`'s skew floor* (`emafit.f:1487`, clamps skew `>= max(-1.41, skxmax)`,
        `skxmax` itself a no-op since `lskewXmax` defaults `.FALSE.`) -- a no-op at Cains
        Coulee's -0.6 to -0.8 magnitude, nowhere near -1.41.
      * *`nG` recomputed every EMA iteration instead of once* -- would have been a genuine
        algorithmic difference from flowfreq's own "compute `nG` once, from the converged
        at-site skew, then iterate" approach. Checked directly against `p3est_ema`
        (`emafit.f:1149`): `nG = n*Wd*as_G_mse/r_G_mse` is computed once, before the iteration
        loop (`emafit.f:1194`), from the `Wd`/`as_G_mse` values passed in as arguments -- same
        structure flowfreq uses. Not the explanation.
      * *Replicating `emafitpr`'s own internal call sequence* -- `mse_ema(kmom=1)`, then
        `(kmom=2)`, then `mseg_all` (`as_G_mse`), then `mseg_all` again with the different,
        uncensored "Syst"/ERL group (`as_G_mse_Syst`), then `mseg_all` once more -- via
        standalone oracle calls in that exact order never reproduces the drift; `as_G_mse`
        stays at 0.0749 throughout. So it is not simply "enough repeated `mse_ema`/`mseg_all`
        calls with varying arguments," the way `mseg_all_sub`'s documented bug is.

      What was **not** ruled out, for lack of a way to isolate it further without transcribing
      large parts of `emafitpr`/`p3est_ema` (out of scope -- flowfreq already has its own native
      EMA fit, verified separately): something in `emafitpr`'s *first* internal fitting pass
      (MGBT, or the initial at-site-only `p3est_ema` call under `at_site_option='B17B'`,
      `emafit.f:745-754`, which runs before the ADJE-branch `mseg_all` call this item traces)
      leaves `mseg_all`/`at_site_option`-adjacent state in a condition that a standalone
      `mseg_all_sub` call from a clean process cannot reproduce. Given `at_site_option` is a
      Fortran `COMMON` variable explicitly toggled `'B17B'` → (reset) `at_site_default`/`at_site_std`
      partway through `emafitpr` (`emafit.f:751`, `790`), and given this whole class of routine
      already carries one confirmed `SAVE`d-state bug, a *second*, subtler one in the same
      family is the leading hypothesis -- but it is a hypothesis, not a confirmed finding, and is
      recorded as such.

      **Consequence**: the `skew_weighted` xfail's reason was rewritten to this account (both in
      `tests/fortran_parity/test_wymt_vs_golden.py` and a new,
      dedicated `tests/fortran_parity/test_fortran_oracles.py::TestCainsCouleeAsGMseDiscrepancy`,
      which pins the 3e-8 `mse_ema` match and the 0.0749-vs-0.2212 `mseg_all_sub` disagreement as
      committed, reproducible assertions rather than prose). flowfreq's own computation is left
      as is -- there is no principled way to deliberately reproduce a number whose mechanism
      is not understood, and doing so would mean curve-fitting to one data point rather than
      transcribing an understood routine, which is what every other line of this port has been.
      It is entirely possible flowfreq's native fit is *more* correct here than the golden
      reference, not less; that is unresolved, not something to act on speculatively.

### Follow-ups found while clearing P1 and P2

Small, specified, none blocking. (The second-parity-case item that used to head this list is
done — see the P3 table above and the Done section.)

- [x] **`plot_peak_flows_with_thresholds` dedupe -- library side done.**
      `streamlit_app.py`'s own `plot_peak_timeseries` carried three things the library
      function didn't: return-period lines, a max-peak recurrence annotation, and hollow bars
      for peaks below a PILF/MGBT cut. All three are now on `plot_peak_flows_with_thresholds`:
      the first two via `lp3_params=(mean_log, std_log, skew)` (dotted reference lines at each
      of `return_periods`, default 2/5/10/25/50/100-yr, each labelled at the right margin, plus
      -- when `annotate_max_peak=True`, the default -- an annotation on the largest peak giving
      its recurrence interval, read exactly off `core.log_pearson3_cdf` rather than interpolated
      off the drawn lines; a recurrence beyond 10x the largest requested return period, including
      the CDF's exact 0/1 saturation, reports `"> {max return period}-yr"` instead of a
      meaningless multi-billion-year number); the third by extending the existing
      `mgbt_threshold` parameter (previously line-only) to also hollow-out bars below it, plus a
      new `mgbt_threshold_source` argument for the line's label (`"override"` vs. the default
      `"MGBT"`), matching the app's own PILF-source labelling. Tested in
      `tests/test_freq_plot.py::TestReturnPeriodLinesAndMaxPeakAnnotation` and
      `TestPlotPeakFlowsWithThresholds` (censored-hollow-bar count, no-threshold-stays-solid,
      source label). App-side switch (calling this instead of `plot_peak_timeseries`, then
      deleting the app's copy) is in progress on `flowfreq-app`'s `v0.5.0` bump branch --
      see that repo's `TODO.md` for status.

- [x] **`FrequencyComparator` compares every parameter by percent difference.** That was the
      wrong metric for skew, which legitimately crosses zero: Big Sandy's reference at-site
      skew is 0.0066, so an absolute gap of 0.016 read as 249% and dominated `max_diff_pct`. Done
      in commit `dc563c5` (skews now go to their own `skew_diffs` dict, compared in skew units
      against `skew_tolerance_abs` and excluded from `max_diff_pct`) but this checkbox was never
      updated to say so. Verified still correct and covered by
      `tests/validation/test_comparisons.py::TestSkewComparedInSkewUnits`, including the exact
      Big Sandy numbers above (`test_a_near_zero_skew_no_longer_dominates_max_diff`).

- [ ] **`origin/dev` can be deleted.** The read-and-judge pass is done — see the commit
      "Port extra_curves from dev". `extra_curves` was the one library delta worth keeping
      and is on this branch; `flowfreq/setup.py` is stale packaging that contradicts
      `pyproject.toml`; everything else is superseded or older than main. Tip is `86cb147`,
      recorded here so the branch is recoverable by SHA. Left undeleted deliberately: that is
      not reversible from a commit, so it is the owner's call.

### Blocked

- [ ] **Tag pushes and branch deletes return HTTP 403** from a Claude Code session -- a
      credential boundary, not a transient failure. `v0.4.0` is now tagged (verified via
      `git ls-remote --tags origin` on 2026-09-05), and `typecheck`/`tests-engine-report` no
      longer exist on `flowfreq`'s remote at all -- both resolved, presumably from a local
      clone, so nothing outstanding there. What remains: delete the stray `parity-12363000`
      branch -- it is on **`flowfreq`** itself, not `hydrolib` as this item previously said
      (`git ls-remote --heads origin` on both repos, 2026-09-05: present on flowfreq, absent
      on hydrolib). Historically this class of block also hit `v0.2.0` and an
      `archive/dev-2026-02` tag.

---

## Done

### P1 — Reduce time to verify and commit

- [x] **One-command verify.** The default selection lives in `pyproject.toml`'s `addopts`, so
      a bare `pytest` is the CI selection on every platform. Override from the CLI when you
      want the excluded tests (`-m ""` for everything). `make check` / `make help` still
      exist; CI runs plain `pytest tests/`.

- [x] **Un-gate the test job from lint.** ~~`needs: lint`~~ removed; the jobs are
      independent, and `fail-fast: false` means all four matrix jobs report.

- [x] **Memoize `_mgbt_pvalue`.** Measured first, because the obvious assumption is wrong:
      within one analysis all 22 keys are distinct, so an `lru_cache` buys **0%** for a
      user-facing run. It pays only where the same fit repeats — 66 calls over 22 distinct
      keys across three identical analyses, 67% — which is the suite's access pattern.
      Measured after, rather than assumed: **75.4 s → 13.6 s**, two runs each way, identical
      outcomes, and refitting with the cache bypassed gives bit-identical results (max |diff|
      0.0). Bigger than this entry originally estimated, because the suite shares fixtures
      more heavily than the three-analysis probe suggested. Decorator order is load-bearing
      and now has a test: `staticmethod` outermost, or Python 3.9 breaks and 3.11 does not.

- [x] **Build the Fortran in one CI job.** New `fortran` job runs `make parity`: build,
      assert the extension imports, then run `tests/fortran_parity/`. The import assertion is
      the point — the parity tests call `importorskip`, so without it a failed build would
      skip all twelve live-vs-golden tests and report green. Verified on Linux/CPython 3.11:
      the extension builds from the vendored sources and the committed goldens are not
      drifted.

- [x] **Lint and smoke-test `app/`.** `app/` joined `PKGS`, so the existing lint job covers it
      (it was already clean). `tests/test_streamlit_app.py` imports `app/streamlit_app.py`,
      which is a top-level script, so importing it executes the whole body in Streamlit's bare
      mode — widgets return defaults, `download_data` is False, nothing contacts NWIS.
      Confirmed non-vacuous: a `NameError` in an executed branch turns all 13 tests red. Runs
      in its own `app` job because Streamlit needs Python ≥ 3.10 and the matrix still covers
      3.9.

### Fixed after P1/P2: the MOM PILF threshold gap

- [x] **`MethodOfMoments` now censors on a PILF threshold instead of only reporting it.**
      Peaks below the threshold (Grubbs-Beck by default, or `user_low_outlier_threshold`) are
      dropped from the moments, and quantiles/confidence limits are evaluated at the Bulletin
      17B conditional probability `Pc = P * n / n_conditional` (§4.2.9-4.2.10) rather than the
      requested AEP directly — the standard treatment for a MOM fit with low outliers. An AEP
      whose `Pc` would reach or exceed 1 (the requested return period falls at or below the
      threshold itself) has no conditional-distribution answer and comes back as `NaN`, logged
      once per call rather than raising. A threshold that leaves fewer than 3 conditional peaks
      raises `ValueError` — there is nothing to fit a skew to.

      No Fortran oracle exists for this: peakfq 8.1.0 only implements EMA, so unlike the rest of
      the `var_mom` port this is verified by construction (conditional moments checked against a
      direct fit on the surviving peaks; K-factors checked against the `Pc` formula applied by
      hand) and against `app/ffa_runner`'s existing override tests, not against vendored Fortran.
      `_low_outlier_source()` no longer needs the "reported only" caveat — MOM acts on the
      number it reports now, same as EMA.

### P2 — Cleanup that stops the same confusion recurring

- [x] **Reviewed `dev`'s remaining library deltas.** See the follow-up above for the branch
      itself.

- [x] **Removed the 5.2 MB of Windows binaries** in `flowfreq/peakfqr/`. `.gitignore` now also
      excludes `flowfreq/peakfqr/*.dll`, and the Streamlit vignette no longer tells readers
      the repository ships a prebuilt extension.

- [x] **Decided what `flowfreq/peakfqsa/` is for: nothing.** Deleted — `config.py`,
      `wrapper.py`, `io_converters.py`, `validators.py` (imported by nothing at all) and the
      `.out` parser, with their ~50 mock tests and the dead `requires_peakfqsa` marker. The
      result container survives as `flowfreq/validation/reference.py::ReferenceResult`,
      pointed at references that exist: `from_golden()` reads the committed golden file,
      `from_emafit()` calls the vendored Fortran live. They agree bitwise on Big Sandy.
      `tests/peakfqsa/fixtures/` was never about the binary — it merged into the existing
      `tests/fixtures/`, and that move silently broke `paths.py`'s hardcoded `parents[3]`
      (three passing tests would have skipped rather than failed), now anchored on
      `pyproject.toml`.

- [x] **Wired the PILF override into the Streamlit UI.** Sidebar control, both `run_ffa` call
      sites, and the refit trigger. Peaks below the applied cut are drawn hollow under a red
      threshold line, so the control changes something a user can see.

### Found and fixed while clearing the above

- [x] **`flowfreq benchmark` could not work for an installed user.** Three defects, all
      pre-existing: it imported fixture data from `tests/`, which is not packaged, so the
      command raised `ModuleNotFoundError` outside a source checkout; two of its three
      registered benchmarks carried no peaks at all and errored on every run; and the third
      validated against the 2012 PeakfqSA manual, which peakfq 8.1.0 does not reproduce, so a
      correct implementation was guaranteed a FAIL. Case data now ships in
      `flowfreq/validation/data/`, tolerances are measured, and known deviations are reported
      rather than compared. 1/1 passed, from any directory.

- [x] **`Benchmark.run_native()` dropped the perception thresholds**, comparing a
      systematic-only fit against the reference's censored one and calling the modelling
      difference an error.

### Per-routine oracles for the port

- [x] **Extended `_emafort.pyf`** with `mseg_all_sub`, `detratsub`, `var_mom` and `moms_p3`.
      The signature file is still not merely a filter — without it f2py tries to wrap
      QUADPACK's `dqag` and the build fails in generated C — so each symbol is listed
      deliberately. `mseg_all_sub` reproduces `emafitpr`'s `as_G_mse_o` exactly on all three
      parity cases, which makes it the oracle the ADJE work is written against; `detratsub`
      reproduces `Wdout` to 1e-9 given the post-MGBT thresholds. See the P3 entry for the
      two defects this immediately surfaced.

### Parity beyond Big Sandy

- [x] **Two Wyoming/Montana parity cases**, from peaks the repository already vendored
      (`wymt_ffa_2022A_EMPdata_7_4.csv`). Both are contiguous systematic records with no
      historic peaks and no zero flows, so their EMA inputs are unambiguous. Powder River
      (85 years, no PILFs) and Cains Coulee (32 years, 11 PILFs under MGBT). Registered in
      `CASES`, goldens generated from the vendored 8.1.0 Fortran; the peakfq 7.4 numbers in
      those CSVs are a cross-check only, never a parity target.

      `test_live_vs_golden.py` is now parametrised over every registered case rather than
      Big Sandy alone, which is what says whether its tolerances generalise. They do, with
      room: the 1-ulp conditioning response is **0.0** on Powder River and ~1e-13 on Cains
      Coulee, against Big Sandy's 1e-5 to 1e-4. Censoring drives the conditioning, not record
      length — Big Sandy's 37 censored intervals are why it is the ill-conditioned one, and
      the tolerances calibrated on it are the loosest any case needs. Minimum headroom across
      all three: 10.7x.

      Both new modules skip cleanly when the reference tree is absent — verified by hiding
      `vendor/peakfqr/inst/testdata` and re-running: 406 passed, 44 skipped, exit 0. The
      `requires_peakfqr_testdata` marker alone does not skip anything, so it is paired with
      a `skipif` on `TESTDATA_AVAILABLE`, the same way `tests/test_r_fixtures.py` does it.

### The structural half of the skew defect

- [x] **The regional skew is weighted inside the EMA fixed point.** peakfq does not average
      two skews after fitting — the explicit average is commented out at `emafit.f:1259`, and
      `moms_p3` folds the regional skew in as `nG` pseudo-observations every iteration, where
      it also moves the mean and variance. Two bugs fixed: the bias corrections `c2`/`c3`
      apply to exact peaks only, not to censored intervals' expected moments; and the
      weighting belongs in the loop, not after it. Big Sandy, against the golden file:

      | | before | after | peakfq 8.1.0 |
      |---|---:|---:|---:|
      | `mean_log` | 0.05% | **0.01%** | 3.717508 |
      | `std_log` | 1.72% | **0.63%** | 0.291043 |
      | at-site skew | 0.0165 abs | **0.0046** | 0.006601 |
      | weighted skew | 35.4% | **24.1%** | −0.156306 |
      | Q100 | — | **0.88%** | 22 959.4 |

---

The 16 phases below were completed in February against the assumption that PeakfqSA was
a standalone binary. It is not, and `flowfreq/peakfqsa/` wraps an executable that does
not exist. Kept for the record; see `AGENT_BUILD_INSTRUCTIONS_Claude.md`.

---

## peakfqr Reference Notes

### 1. Fortran Call Signatures

**Primary routine: `emafitpr`** (`vendor/peakfqr/src/emafit.f`)

```
subroutine emafitpr(n, ql, qu, tl, tu, dtype,
    reg_M, reg_M_mse, reg_SD, reg_SD_mse, r_G, r_G_mse,
    gbthrsh0, pq, nq, eps, wght_opt_n,
    gbval, gbns, gbnzero, gbnlow, gbp,
    gbqs, as_G_mse_o, as_G_mse_Syst_o,
    as_G_PRL_o, cmoms, yp, ci_low, ci_high, var_est, Wdout,
    qlema, quema, tlema, tuema, nu)
```

**Input arguments:**
| Arg | Type | Description |
|-----|------|-------------|
| n | i*4 | Number of observations (censored/uncensored) |
| ql(n) | r*8 | Lower bounds on log10(floods) |
| qu(n) | r*8 | Upper bounds on log10(floods) |
| tl(n) | r*8 | Lower bounds on log10(flood thresholds) |
| tu(n) | r*8 | Upper bounds on log10(flood thresholds) |
| dtype(n) | i*4 | 0=systematic, 1=historic |
| reg_M | r*8 | Regional mean |
| reg_M_mse | r*8 | MSE of regional mean |
| reg_SD | r*8 | Regional standard deviation |
| reg_SD_mse | r*8 | MSE of regional SD |
| r_G | r*8 | Regional skew |
| r_G_mse | r*8 | MSE of regional skew (encoding below) |
| gbthrsh0 | r*8 | MGBT control (encoding below) |
| pq(nq) | r*8 | Quantile probabilities (1-AEP) |
| nq | i*4 | Number of quantiles |
| eps | r*8 | CI coverage (0.90 = 90%) |
| wght_opt_n | i*4 | Skew weighting: 1=HWN, 2=ERL, 3=INV |

**Output arguments:**
| Arg | Type | Description |
|-----|------|-------------|
| gbval | r*8 | MGBT low outlier critical value |
| gbns | i*4 | Number of peaks used in MGBT |
| gbnzero | i*4 | Number of zero flows |
| gbnlow | i*4 | Number of PILFs detected |
| gbp(20000) | r*8 | MGBT p-values |
| gbqs(20000) | r*8 | Peaks used in MGBT |
| as_G_mse_o | r*8 | At-site skew MSE |
| as_G_mse_Syst_o | r*8 | At-site skew MSE (gaged only) |
| as_G_PRL_o | r*8 | Pseudo effective record length |
| cmoms(3,3) | r*8 | Central moments matrix (see below) |
| yp(nq) | r*8 | Quantile estimates (log10) |
| ci_low(nq) | r*8 | Lower CI bounds (log10) |
| ci_high(nq) | r*8 | Upper CI bounds (log10) |
| var_est(nq) | r*8 | Variance of estimates |
| Wdout | r*8 | Censored data adjustment factor |
| qlema..tuema(25000) | r*8 | EMA-adjusted data representation |
| nu(nq) | r*8 | Degrees of freedom for CI |

**cmoms matrix layout:**
- Column 1: Using regional info + at-site data → `cmoms[1,1]=Mean`, `cmoms[2,1]=Variance`, `cmoms[3,1]=Skew`
- Column 2: At-site only → `cmoms[1,2]=AtSiteMean`, `cmoms[2,2]=AtSiteVariance`, `cmoms[3,2]=AtSiteSkew`
- Column 3: B17B MSE formula for at-site → `cmoms[*,3]`

**Other Fortran routines called from R:**
- `PLOTPOSHS` — Hirsch-Stedinger plotting positions
- `EXPMOMCDERIV` — Expected central moments derivative
- `DEXPECT` — Expected noncentral moments derivative
- `detratsub` — Determinant ratio for skew weighting (Halloween method)
- `moms_p3` — Pearson Type III expected moments
- `p3est_ema` — EMA moment estimation with regional weighting
- `var_mom` — Variance-covariance of moment estimators
- `qP3sub` — Pearson Type III quantile (inverse CDF)
- `mseg_all_sub` — Mean-square error of at-site skew

### 2. EMA Parameter Conventions

**Flow intervals** (ql, qu):
- Exact observation: `ql[i] = qu[i] = log10(Q)`
- Less-than (censored): `ql[i] = log10(Qmin)` (~-20), `qu[i] = log10(threshold)`
- Greater-than: `ql[i] = log10(Q)`, `qu[i] = log10(Qmax)` (~+20)
- Interval: `ql[i] = log10(lower)`, `qu[i] = log10(upper)`
- Zero flows: enter as `lmissing = -80.0`

**Perception thresholds** (tl, tu):
- Systematic period: `tl = log10(0)` (uses Qmin=1e-20), `tu = log10(1e20)`
- Historical period: `tl = log10(threshold)`, `tu = log10(1e20)`
- Missing data (no info): `tl = log10(1e20)`, `tu = log10(1e20)` (both infinity)

**R wrapper** (`fortranWrappers.R:emafit`):
- Accepts real-space data, converts to log10 before calling Fortran
- pq = 1 - AEPs (converts exceedance prob to non-exceedance)
- Converts output quantiles back: `10^yp`, `10^ci_low`, `10^ci_high`
- Default AEPs: 0.995, 0.99, 0.98, 0.975, 0.96, 0.95, 0.90, 0.80, 0.70, 0.667, 0.60, 0.5704, 0.50, 0.4292, 0.40, 0.30, 0.20, 0.10, 0.05, 0.04, 0.025, 0.02, 0.01, 0.005, 0.002

### 3. MGBT Implementation

**Encoding via `gbthrsh0`:**
- `<= -6`: Run MGBT (default; R passes `-99`)
- `> -6`: Use as fixed threshold (R passes `log10(user_value)`)
- `~-5.9`: Disable low outlier test (R passes `log10(Qmin)`)

**R-side logic** (`fortranWrappers.R` lines 93-128):
- FIXED: `LOthresh >= 1e-6` → set to `log10(LOthresh)`
- NONE: `LOthresh > 1e-99` → set to `log10(Qmin)` (effectively disables)
- MGBT: `LOthresh <= 1e-99` → set to `-99`, validates ≥5 nonzero peaks, checks upper thresholds > median

**MGBT output interpretation:**
- `gbnlow` = number of PILFs detected
- `gbqs[1:nPILFs]` = peak values identified as PILFs
- `gbp[1:nPILFs]` = associated p-values
- PILF threshold = `10^gbval`

### 4. Confidence Interval Method

- Uses **Inverse Modified Cholesky Gaussian Quadrature** (added Oct 2012, emafit.f)
- CI coverage controlled by `eps` parameter (default 0.90 = 90% CI)
- Output: `ci_low` and `ci_high` in log10 space (5th and 95th percentiles for 90% CI)
- Small-sample correction applied: monotonicity enforcement on CI bounds (lines 485-491)
- Variance of estimate (`var_est`) also returned for each quantile
- This matches peakfqr's approach — no differences noted

### 5. Regional Skew Weighting

**r_G_mse encoding (four cases):**
1. `r_G_mse = 0`: Use fixed `g = r_G` with MSE=0 ("Generalized skew, no error")
2. `-98 < r_G_mse < 0`: Use fixed `g = r_G` with MSE = `-r_G_mse` ("Generalized skew, MSE > 0")
3. `0 < r_G_mse < 1e10`: Weighted average of at-site and regional skew ("Weighted")
4. `r_G_mse > 1e10`: Use at-site skew only ("Station")

**R conversion** (`main.R` lines 415-424):
- Station: `SkewMSE = -1e99`
- Weighted: `SkewMSE = SkewSE^2`
- Regional (generalized): `SkewMSE = -(SkewSE^2)` (negative)

**Weighting formula** (Bulletin 17C):
- `skew_weighted = (skew_atsite * MSE_regional + skew_regional * MSE_atsite) / (MSE_regional + MSE_atsite)`
- Halloween method (HWN, default): Applies determinant ratio `Wd` correction for censored data
- `detrat()` computes Wd from EXPMOMCDERIV subroutine
- Wd=1.0 when no censored data present (equivalent to INV method)

### 6. Output Field Mapping

**peakfqr output → PeakfqSAResult mapping:**

| peakfqr field | Source | PeakfqSAResult field |
|--------------|--------|---------------------|
| Mean | cmoms[1,1] | parameters["mean_log"] |
| StandDev | sqrt(cmoms[2,1]) | parameters["std_log"] |
| Skew | cmoms[3,1] | parameters["skew_weighted"] |
| AtSiteMean | cmoms[1,2] | parameters["mean_log_at_site"] |
| AtSiteStandDev | sqrt(cmoms[2,2]) | parameters["std_log_at_site"] |
| AtSiteSkew | cmoms[3,2] | parameters["skew_at_site"] |
| AtSiteMSEG | EMAout[[24]] | parameters["mse_skew"] |
| AtSiteMSEG_GagedOnly | EMAout[[25]] | parameters["mse_skew_systematic"] |
| RegSkew | input rG | parameters["regional_skew"] |
| RegMSEG | input rGmse | parameters["regional_skew_mse"] |
| RecordLength | n | n_peaks |
| HistPeaks | count(dtype==1) | n_historical |
| PILF_Method | derived | low_outlier_method |
| PILF_Thresh | 10^gbval | low_outlier_threshold |
| PILFs | gbnlow | low_outlier_count |
| PILF_0s | gbnzero | (informational) |
| WeightOpt | input | (config) |
| WeightCo (Wd) | EMAout[[32]] | (informational) |
| EXC_Prob | AEPs | quantiles keys |
| Estimate | 10^yp | quantiles values |
| Variance | var_est | (per-quantile) |
| Conf_Low | 10^ci_low | confidence_intervals lower |
| Conf_Up | 10^ci_high | confidence_intervals upper |

### 7. Edge Cases Handled

- **Zero flows**: Enter as `lmissing = -80.0` in log space; counted separately as PILF_0s
- **All values positive required**: R validates `all(QT[,c("ql","qu","tl","tu")] > 0)` before log transform
- **Minimum data**: Requires ≥3 rows for skew calculation; peakfq() requires ≥10 total or ≥8 exact
- **MGBT with few peaks**: Requires >5 non-zero exactly-known peaks
- **Upper threshold < median**: Rejected when using MGBT (causes numerical issues)
- **Censored data with code 4**: Converted to interval `[0, stated_value]`
- **Greater-than peaks (code 8)**: Converted to interval `[stated_value, infinity]`
- **Historic peaks (code 7)**: Set perception thresholds based on lowest historic peak in period
- **Regulated/urbanized (codes 6, C)**: Excluded by default; included if `Urb/Reg = Yes`
- **Dam failure (code 3), Opportunistic (code O)**: Always excluded

---

## Phase 0: Setup & Reference

- [x] Step 0a: Read peakfqr reference repository
- [x] Step 0a: Document Fortran call signatures
- [x] Step 0a: Document EMA parameter conventions
- [x] Step 0a: Document MGBT implementation
- [x] Step 0a: Document CI method
- [x] Step 0a: Document regional skew weighting
- [x] Step 0a: Map output fields to PeakfqSAResult
- [x] Step 0a: Document edge cases
- [x] Step 0b: Scan existing FlowFreq codebase
- [x] Step 0b: Generate this TODO list

## Phase 1: Information Gathering

- [x] Step 1: Resolve questions (peakfqr = reference code, not PeakfqSA binary)

## Phase 2: Environment Setup

- [x] Step 2a: Verify project structure
- [x] Step 2b: Install dependencies
- [x] Step 2c: Run baseline tests and record results

## Phase 3: Directory Structure

- [x] Step 3: Create `flowfreq/peakfqsa/__init__.py`
- [x] Step 3: Create `flowfreq/peakfqsa/config.py` (stub)
- [x] Step 3: Create `flowfreq/peakfqsa/wrapper.py` (stub)
- [x] Step 3: Create `flowfreq/peakfqsa/io_converters.py` (stub)
- [x] Step 3: Create `flowfreq/peakfqsa/parsers.py` (stub)
- [x] Step 3: Create `flowfreq/peakfqsa/validators.py` (stub)
- [x] Step 3: Create `flowfreq/validation/__init__.py`
- [x] Step 3: Create `flowfreq/validation/benchmarks.py` (stub)
- [x] Step 3: Create `flowfreq/validation/comparisons.py` (stub)
- [x] Step 3: Create `flowfreq/validation/reports.py` (stub)
- [x] Step 3: Create `tests/peakfqsa/__init__.py`
- [x] Step 3: Create `tests/peakfqsa/test_config.py`
- [x] Step 3: Create `tests/peakfqsa/test_wrapper.py`
- [x] Step 3: Create `tests/peakfqsa/test_io_converters.py`
- [x] Step 3: Create `tests/peakfqsa/test_parsers.py`
- [x] Step 3: Create `tests/peakfqsa/fixtures/__init__.py`
- [x] Step 3: Create `tests/peakfqsa/fixtures/big_sandy.py`
- [x] Step 3: Create `tests/validation/__init__.py`
- [x] Step 3: Create `tests/validation/test_benchmarks.py`
- [x] Step 3: Create `tests/integration/__init__.py`
- [x] Step 3: Create `tests/integration/test_hybrid_workflow.py`

## Phase 4: Test Fixtures

- [x] Step 4: Create Big Sandy River fixture data
- [x] Step 4: Create sample PeakfqSA output fixtures

## Phase 5: Configuration Module

- [x] Step 5: Implement `PeakfqSAConfig` dataclass
- [x] Step 5: Implement `find_peakfqsa()` discovery function
- [x] Step 5: Implement `validate_peakfqsa()` validation
- [x] Step 5: Implement `PeakfqSANotFoundError`
- [x] Step 5: Write tests in `test_config.py`
- [x] Step 5: Run tests and fix

## Phase 6: I/O Converters

- [x] Step 6: Implement `SpecificationFile` class
- [x] Step 6: Implement `DataFile` class
- [x] Step 6: Implement `from_analysis_params()`, `to_string()`, `write()`, `validate()`
- [x] Step 6: Write tests with Big Sandy expected output
- [x] Step 6: Run tests and fix

## Phase 7: PeakfqSA Wrapper

- [x] Step 7: Implement `PeakfqSAWrapper` class
- [x] Step 7: Implement `run()`, `_write_input_files()`, `_execute()`, `_parse_output_text()`
- [x] Step 7: Implement error classes (NotFound, Execution, Timeout, Parse)
- [x] Step 7: Register `requires_peakfqsa` marker
- [x] Step 7: Write mock-based tests
- [x] Step 7: Run tests and fix

## Phase 8: Output Parser

- [x] Step 8: Implement `PeakfqSAResult` dataclass
- [x] Step 8: Implement `.out` file parser with regex patterns
- [x] Step 8: Write tests with fixture output text
- [x] Step 8: Run tests and fix

## Phase 9: Comparison Engine

- [x] Step 9: Implement `ComparisonResult` dataclass
- [x] Step 9: Implement `FrequencyComparator` class
- [x] Step 9: Write tests (identical results, tolerance boundary)
- [x] Step 9: Run tests and fix

## Phase 10: FrequencyAnalyzer API Update

- [x] Step 10: Add `to_comparison_dict()` to Bulletin17C
- [x] Step 10: Add `validate()` method to Bulletin17C
- [x] Step 10: Write backward-compatibility test

## Phase 11: Integration Tests

- [x] Step 11: Write Big Sandy systematic-only test
- [x] Step 11: Write Big Sandy with historical test (documents convergence limitation)
- [x] Step 11: Write validation workflow test

## Phase 12: Benchmark Module

- [x] Step 12: Implement `Benchmark` class with `run_native()`, `validate_against_expected()`
- [x] Step 12: Register Big Sandy benchmark
- [x] Step 12: Implement `run_all_benchmarks()` and `print_benchmark_report()`
- [x] Step 12: Implement text and JSON report generators
- [x] Step 12: Write tests

## Phase 13: CLI Commands

- [x] Step 13: Implement `flowfreq validate` command
- [x] Step 13: Implement `flowfreq benchmark` command
- [x] Step 13: Register in `pyproject.toml`

## Phase 14: Documentation

- [x] Step 14a: All new modules have NumPy-format docstrings
- [x] Step 14b: CLAUDE.md updated with hybrid 17C architecture

## Phase 15: Final Quality Check

- [x] Step 15: Run black + isort
- [x] Step 15: Run full test suite (96/96 passing)
- [x] Step 15: Check for remaining TODOs (0 in source code)

## Phase 16: Update TODO.md

- [x] Step 16: Check off all completed items
- [x] Step 16: Update status block

---

## Known Limitations

- **Native EMA with historical data**: The native Python EMA implementation can diverge
  (produce NaN) when processing historical/censored intervals due to numerical instability
  in the quad integration. This is a known limitation of the existing `bulletin17c.py`
  implementation. Systematic-only analyses converge reliably.

- **PeakfqSA binary not available**: The PeakfqSA Fortran executable is not available as
  a standalone binary. The peakfqr R package contains the authoritative Fortran source.
  The wrapper module is implemented and tested via mocks but cannot be used end-to-end
  without a standalone executable.

- **Native EMA confidence intervals** — superseded. The Fortran is now vendored, the
  extension builds, and the question this entry asked has been answered: see the Status
  block above and `docs/FORTRAN_UPLOAD.md` §6.0 and §6.0b. Two findings worth carrying
  forward. First, the residual is interval *shape*, not variance magnitude. Second, the
  2012 PeakfqSA manual values this was measured against are not reproducible by peakfq
  8.1.0 at all, so part of the original discrepancy was never a defect — the fixture now
  carries `PEAKFQ_810_*` values for parity work alongside the 2012 ones.

## Resolved Questions

- PeakfqSA: Not a standalone binary. Use the vendored `vendor/peakfqr/src/` Fortran.
  `flowfreq/peakfqsa/` wraps a binary that does not exist and is mock-tested only.
- Reference material: vendored into `vendor/peakfqr/` (CC0). No external workspace needed.
- Big Sandy 2012 manual values: not reproducible by peakfq 8.1.0 — the HWN skew weighting
  postdates the manual and diverges by design on censored records.
- MGBT: verified line-by-line against the Fortran (`GGBCRITP`/`FP_TNC_CDF`), validated on
  Orestimba Creek (USGS 11274500).
- Python: CI tests 3.9–3.12; `requires-python >= 3.9`.
- Regional skew defaults: -0.302 / MSE 0.3025 (Bulletin 17C national map)
- Git: Commit to dev branch, no push unless asked
- FrequencyAnalyzer API: Added validate() and to_comparison_dict() to Bulletin17C facade
