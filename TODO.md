# TODO — HydroLib Hybrid 17C Implementation

## Status
Last updated: 2026-08-31
Tests: **407 passed, 1 deselected, 5 xfailed** in ~18 s (was ~75 s; see the MGBT memoization
below). CI green on main.
Fortran reference: **vendored** at `vendor/peakfqr/` (peakfq 8.1.0, CC0).
Fortran bridge: builds from those sources via `python build_fortran/build.py`
(gfortran + meson) and is now **built and checked in CI** by `make parity`.

Every P1 and P2 item is done. What is left is P3: two numerical defects that turn out to be
the same missing machinery, specified precisely below.

---

## Open Items (prioritised)

### P3 — The open numerical defects

Both are `xfail(strict=True)` in `tests/fortran_parity/test_native_vs_golden.py`, so the build
fails the moment either starts passing. Nothing else is blocked on them.

They are **not two problems**. Both bottom out in `var_mom` and its dependency tree, which is
the one piece of the reference implementation that has never been ported. Read this section
as one item with two symptoms.

- [ ] **Port `var_mom` and the routines it needs.** ~1,100 lines of Fortran in `emafit.f`
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

      Verify against the bridge, not by eye: `emafitpr` returns `as_G_mse_o`,
      `as_G_mse_Syst_o`, `as_G_PRL_o` and `var_est`, so every intermediate has an oracle.
      `make parity` runs it.

- [ ] **Skew weighting — 24% left, and it is all `as_G_mse`.** The structural half is done
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

      Also unimplemented: `detrat`, the Halloween determinant ratio. `emafit.f:763` applies it
      only when the at-site skew is ≥ 0.04 in magnitude, and Big Sandy's is 0.0066, so `Wd`
      is 1 here either way — but a record above that floor currently gets INV weighting, not
      HWN. Big Sandy cannot detect that regression; a second parity case would.

- [ ] **Confidence-interval shape.** `compute_confidence_limits()` forms `log_Q ± z·se`,
      symmetric by construction (ratio 1.000 at every AEP). peakfq skews right with return
      period — 1.03 → 1.31 → 1.41 at AEP 0.1 / 0.02 / 0.01. **Shape, not size**: total width
      is within ~2% and point estimates within 0.5%, so the variance magnitude is essentially
      right and `var_mom`'s output is not the suspect — its absence is.

      The formula is no longer a mystery. `ci_ema_m3b` (`emafit.f:1853`) is short and fully
      readable:

      ```
      beta1     = cov(yp, syp) / var(yp)            # regression of the s.e. on the quantile
      var_xsi_d = var(syp) - cov(yp, syp)**2/var(yp)
      nu        = max(5, 0.5 * var(yp)/var_xsi_d)   # Student t degrees of freedom
      t         = t_nu((1 + eps)/2)
      ci_high   = yp + sqrt(var(yp)) *  t / max(0.5, 1 - beta1*t)
      ci_low    = yp + sqrt(var(yp)) * -t / max(0.5, 1 + beta1*t)
      ```

      The whole asymmetry is the `1 - beta1*t` denominator: it shrinks on the high side and
      grows on the low side. Implementing this is an afternoon. Getting `beta1` is not —
      `cov(yp, syp)` comes from `VAR_EMAB`, hence from `var_mom`. Implementing `ci_ema_m3b`
      alone with `beta1 = 0` reproduces exactly today's symmetric interval, so there is no
      partial credit available here: the port comes first.

      The pseudo effective record length (`as_G_PRL_o`, 54.373 for Big Sandy) is
      `eff_n * as_G_mse_Syst / as_G_mse` (`emafit.f:758`) — two more `mseg_all` calls, so it
      is blocked on the same port.

### Follow-ups found while clearing P1 and P2

Small, specified, none blocking.

- [ ] **`plot_peak_flows_with_thresholds` is still uncalled.** `app/streamlit_app.py` has its
      own `plot_peak_timeseries`, which carries return-period lines and the max-peak
      recurrence annotation that the library function does not; switching to the library one
      today would lose features. The dedupe is: move those two features into
      `hydrolib/freq_plot.py`, then delete the app copy. One peak-flow plotter, tested.

- [ ] **`MethodOfMoments` ignores a user PILF threshold.** `Bulletin17C.run_analysis` does not
      pass `user_low_outlier_threshold` down the MOM path (`bulletin17c.py:1494`), and MOM
      does not censor low outliers at all — it computes a Grubbs-Beck value for reporting
      only. So when EMA fails to converge and `app/ffa_runner.run_ffa` falls back to MOM, the
      user's override is dropped. `_low_outlier_source()` reports that honestly rather than
      claiming a threshold that had no effect, but the underlying limitation stands. Doing
      this properly means censoring plus the conditional-probability adjustment, which is a
      real piece of B17B, not a one-liner.

- [ ] **`FrequencyComparator` compares every parameter by percent difference.** That is the
      wrong metric for skew, which legitimately crosses zero: Big Sandy's reference at-site
      skew is 0.0066, so an absolute gap of 0.016 read as 249% and dominated `max_diff_pct`.
      A denominator floor (`parameter_scale_floor`, default 0.1) stops it hiding the rest of
      the report, but comparing skew in skew units with its own tolerance would be better.

- [ ] **`origin/dev` can be deleted.** The read-and-judge pass is done — see the commit
      "Port extra_curves from dev". `extra_curves` was the one library delta worth keeping
      and is on this branch; `hydrolib/setup.py` is stale packaging that contradicts
      `pyproject.toml`; everything else is superseded or older than main. Tip is `86cb147`,
      recorded here so the branch is recoverable by SHA. Left undeleted deliberately: that is
      not reversible from a commit, so it is the owner's call.

### Blocked

- [ ] **Tag pushes return HTTP 403** from the agent environment, so neither `v0.2.0` nor an
      `archive/dev-2026-02` tag could be pushed. Needs tag-push permission, or someone
      pushing tags from a local clone.

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

### P2 — Cleanup that stops the same confusion recurring

- [x] **Reviewed `dev`'s remaining library deltas.** See the follow-up above for the branch
      itself.

- [x] **Removed the 5.2 MB of Windows binaries** in `hydrolib/peakfqr/`. `.gitignore` now also
      excludes `hydrolib/peakfqr/*.dll`, and the Streamlit vignette no longer tells readers
      the repository ships a prebuilt extension.

- [x] **Decided what `hydrolib/peakfqsa/` is for: nothing.** Deleted — `config.py`,
      `wrapper.py`, `io_converters.py`, `validators.py` (imported by nothing at all) and the
      `.out` parser, with their ~50 mock tests and the dead `requires_peakfqsa` marker. The
      result container survives as `hydrolib/validation/reference.py::ReferenceResult`,
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

- [x] **`hydrolib benchmark` could not work for an installed user.** Three defects, all
      pre-existing: it imported fixture data from `tests/`, which is not packaged, so the
      command raised `ModuleNotFoundError` outside a source checkout; two of its three
      registered benchmarks carried no peaks at all and errored on every run; and the third
      validated against the 2012 PeakfqSA manual, which peakfq 8.1.0 does not reproduce, so a
      correct implementation was guaranteed a FAIL. Case data now ships in
      `hydrolib/validation/data/`, tolerances are measured, and known deviations are reported
      rather than compared. 1/1 passed, from any directory.

- [x] **`Benchmark.run_native()` dropped the perception thresholds**, comparing a
      systematic-only fit against the reference's censored one and calling the modelling
      difference an error.

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
a standalone binary. It is not, and `hydrolib/peakfqsa/` wraps an executable that does
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
- [x] Step 0b: Scan existing HydroLib codebase
- [x] Step 0b: Generate this TODO list

## Phase 1: Information Gathering

- [x] Step 1: Resolve questions (peakfqr = reference code, not PeakfqSA binary)

## Phase 2: Environment Setup

- [x] Step 2a: Verify project structure
- [x] Step 2b: Install dependencies
- [x] Step 2c: Run baseline tests and record results

## Phase 3: Directory Structure

- [x] Step 3: Create `hydrolib/peakfqsa/__init__.py`
- [x] Step 3: Create `hydrolib/peakfqsa/config.py` (stub)
- [x] Step 3: Create `hydrolib/peakfqsa/wrapper.py` (stub)
- [x] Step 3: Create `hydrolib/peakfqsa/io_converters.py` (stub)
- [x] Step 3: Create `hydrolib/peakfqsa/parsers.py` (stub)
- [x] Step 3: Create `hydrolib/peakfqsa/validators.py` (stub)
- [x] Step 3: Create `hydrolib/validation/__init__.py`
- [x] Step 3: Create `hydrolib/validation/benchmarks.py` (stub)
- [x] Step 3: Create `hydrolib/validation/comparisons.py` (stub)
- [x] Step 3: Create `hydrolib/validation/reports.py` (stub)
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

- [x] Step 13: Implement `hydrolib validate` command
- [x] Step 13: Implement `hydrolib benchmark` command
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
  `hydrolib/peakfqsa/` wraps a binary that does not exist and is mock-tested only.
- Reference material: vendored into `vendor/peakfqr/` (CC0). No external workspace needed.
- Big Sandy 2012 manual values: not reproducible by peakfq 8.1.0 — the HWN skew weighting
  postdates the manual and diverges by design on censored records.
- MGBT: verified line-by-line against the Fortran (`GGBCRITP`/`FP_TNC_CDF`), validated on
  Orestimba Creek (USGS 11274500).
- Python: CI tests 3.9–3.12; `requires-python >= 3.9`.
- Regional skew defaults: -0.302 / MSE 0.3025 (Bulletin 17C national map)
- Git: Commit to dev branch, no push unless asked
- FrequencyAnalyzer API: Added validate() and to_comparison_dict() to Bulletin17C facade
