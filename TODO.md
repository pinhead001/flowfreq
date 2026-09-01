# TODO — HydroLib Hybrid 17C Implementation

## Status
Last updated: 2026-09-01
Tests: **564 passed, 1 skipped, 1 deselected, 9 xfailed** in ~40 s (up from ~29 s -- the ADJE
bias adjustment costs real time on a censored fit; see the "Skew weighting" item). CI green
on main.
Fortran reference: **vendored** at `vendor/peakfqr/` (peakfq 8.1.0, CC0).
Fortran bridge: builds from those sources via `python build_fortran/build.py`
(gfortran + meson) and is now **built and checked in CI** by `make parity`.

Every P1 and P2 item is done. The `var_mom` port (TODO.md P3) is done through Phase 3 and
**wired into `Bulletin17C`**: `mn2mvarb`/`mse_ema`, the ADJE censoring bias adjustment, now
feeds `_regional_skew_equivalent_years`, and Big Sandy's weighted skew matches peakfq 8.1.0
closely enough that the parity suite's `xfail(strict=True)` for it now passes -- so it is no
longer marked xfail, per the whole point of `strict=True`. The "Skew weighting" item below is
effectively closed for records like Big Sandy (at-site skew under the HWN floor); what is left
of P3 is `VAR_EMAB`/`regmoms`/`ci_ema_m3b` (the confidence-interval shape fix, unstarted) and
`detrat` (the Halloween determinant ratio, needed only when at-site skew exceeds 0.04 in
magnitude, which none of the three parity cases needing regional weighting do).

---

## Open Items (prioritised)

### P3 — The open numerical defects

Both are `xfail(strict=True)` in `tests/fortran_parity/test_native_vs_golden.py`, so the build
fails the moment either starts passing. Nothing else is blocked on them.

They are **not two problems**. Both bottom out in `var_mom` and its dependency tree, which is
the one piece of the reference implementation that has never been ported. Read this section
as one item with two symptoms.

**The defect is censoring-specific, and that is now measured rather than assumed.** Two
Wyoming/Montana parity cases were added for exactly this question:

| case | censoring | reference `Wd` | native vs peakfq 8.1.0 |
|---|---|---:|---|
| Powder River 06326500 | none | 1.0 | mean **0.0**, sd **3.7e-14**, at-site skew **4.5e-12**, weighted skew **7.5e-11**; quantiles ≤ 0.10% |
| Cains Coulee 06327450 | 11 PILFs from MGBT | **0.184** | at-site skew 0.122, weighted skew 0.098; quantiles 0.30% to 18.7%, 2.7% at Q100 |

So with nothing censored the native EMA reproduces peakfq to machine precision — the in-loop
regional-skew weighting is *right*, not merely closer. Everything that remains is
censored-path work. Cains Coulee is also the first case in the repository where `detrat`
actually bites, so it is the oracle for that routine.

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
        ADJE's bias adjustment partially undo the cap; hydrolib applies the cap alone. On a
        record longer than 150 years that over-weights the regional skew. No parity case is
        that long, so only the oracle test detects it, and fixing it needs `mse_ema`.

      **Phase 1 done: the leaf layer.** `hydrolib/_p3_moments.py` now carries `m2p`, `m2mn`,
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
      `hydrolib._p3_moments.m_p3` keeps the whole expansion in `mpmath` at 50 decimal digits
      (`_GAMMA_MOMENT_DPS`) rather than only the CDF call, and its k = 4..6 values are
      confirmed against independent arbitrary-precision quadrature of the truncated density
      (not against the Fortran, which is the thing in question) —
      `test_fortran_itself_loses_precision_on_big_sandy_s_censored_group`.

      **Resolved by Phase 2, as it turns out: `var_mom` itself never reaches this regime.**
      `var_mom` clamps `|skew|` to `[0.0632, 1.41]` before ever calling `mP3`/`pP3`
      internally (`emafit.f:2324`, easy to miss reading linearly) — so even Big Sandy's raw
      0.0066 becomes 0.0632 inside `var_mom`, capping alpha at ~1000 rather than the ~9e4 the
      finding above used directly. `hydrolib._var_mom.var_mom` matches the Fortran oracle to
      1e-9 relative on Powder River and Cains Coulee and 1e-3 on Big Sandy (see Phase 2 below
      for where that residual is), not the orders-of-magnitude gap `mP3` alone showed. The
      finding above stands as a real, documented Fortran defect reachable by calling `mP3`
      directly (as the oracle tests do) — it just is not reachable *through* `var_mom`.

      New dependency: `mpmath` (pure Python, no compiled extension), added for exactly this.

      **Phase 2 done: `varb`, `varc`, `d_est`, `expmomderiv` and the `var_mom` composition
      itself.** `hydrolib/_var_mom.py`. `expmomderiv` differentiates `_p3_moments`'s own
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

- [x] **Skew weighting — wired in; `detrat` is what's left, and it doesn't bite these cases.**
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

      Also unimplemented: `detrat`, the Halloween determinant ratio. `emafit.f:763` applies it
      only when the at-site skew is ≥ 0.04 in magnitude, and Big Sandy's is 0.0066, so `Wd`
      is 1 there either way. **Cains Coulee 06327450 now covers it**: its at-site skew is
      −0.708 and its reference `Wd` is **0.184**, so hydrolib's implicit 1 over-weights the
      regional skew by more than fivefold on that site. That case is the acceptance test for
      `detrat`, and its two `xfail(strict=True)` assertions will fail the build the moment
      the routine lands and starts working.

      **Phase 3 done: `mn2mvarb`/`mse_ema` — `as_G_mse` is now computable, in Python, matching
      the Fortran.** `hydrolib/_mse_ema.py`. `mse_ema(nobs, tl, tu, mc, kmom) →`
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
      `hydrolib/_p3_moments.py`, sharing `down` across k) and replacing the 9 `mpmath.diff` calls
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
      edited — but worth knowing before calling `hydrolib.peakfqr`/`hydrolib.validation.reference`
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

      What is still open: `detrat`, the Halloween determinant ratio (below) — needed only when
      at-site skew is ≥ 0.04, which none of the three parity cases needing regional weighting
      are.

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

Small, specified, none blocking. (The second-parity-case item that used to head this list is
done — see the P3 table above and the Done section.)

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
