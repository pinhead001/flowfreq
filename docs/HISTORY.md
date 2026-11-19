# History moved from TODO.md

Moved on 2026-11-19, verbatim, so that `TODO.md` holds what is open. Nothing here is a task. Paths and test counts are as they were when each block was written; several name files that no longer exist (`flowfreq/peakfqsa/`, `tests/peakfqsa/`, `app/`).

### Done — StreamStats Phase 2 follow-ups

Phase 2 itself (NSS flow-statistics estimation) is done -- see "Done -- StreamStats
module Phase 2" below. All three things it found are closed (#80, #82):

- [x] **Region selection without a watershed polygon.** NSS defines several
      independently-calibrated `regressionRegions` per statistic group within a state
      (WA Peak-Flow: four), and confirmed live, nothing filters them by location --
      not `ByLocation` with a bare point (Phase 1 produces no polygon to try instead).
      `estimate_flow_statistics()` currently returns every region whose parameters are
      merely in numeric range, unlabelled as to which is geographically correct; the
      caller must know. Revisit once/if Phase 1's own polygon gap (below) is closed --
      a real polygon through `Scenarios/ByLocation` may resolve both at once.
      *Update 2026-09-27:* the polygon gap is closed (design doc S10), so this is now
      unblocked. **Done 2026-09-27** (addendum S5): `locate_regression_regions` POSTs
      the polygon to `regions/{region}/regressionregions/bylocation`, verified live, and
      `estimate_flow_statistics(watershed_polygon=...)` estimates only the located
      regions. Each carries NSS's area `percent_weight`, plus NSS's own area-weighted
      `areaave` result when a basin spans regions. Goat Creek is `GC1751` (3,290 cfs), not
      `GC1750` (4,370). The old all-in-range behaviour stays available without a polygon or
      with `include_unlocated=True`. Still open from this: WA low-flow regions
      `GC1434`/`GC1556`/`GC1558` have no geometry in NSS and can never be located.
- [x] **`characteristic_codes` does not filter.** Found live 2026-09-27 (design doc
      S10): ss-hydro ignores the `bcLabels` query parameter this module sends -- the
      same 22 characteristics came back for `*`, a comma list and a semicolon list.
      USGS's own notebook sets `bcLabels` inside the POSTed `bcrequest` body,
      semicolon-delimited. Harmless for correctness (a superset comes back) but the
      documented subset-for-speed never happens, and a code not in the region's default
      set (e.g. `ELEV1000`) may never be computed. **Fixed 2026-09-28** (design doc
      S11): the real query parameter is `BCs` (ss-hydro OpenAPI); `bcLabels` is
      ignored in the query and in the body alike. A code the region cannot compute
      comes back as HTTP 200 with value -999 and is now kept out of `characteristics`
      (listed in `WatershedCharacteristics.unavailable`).
- [x] Low-Flow Statistics (`LFS`, the other group WA supports) was exercised only
      through its client-side validation (correctly skipping 3 of 4 regions as
      out-of-range, 1 for a missing `ELEV1000` characteristic) -- no region actually
      returned a Low-Flow estimate live. Only Peak-Flow (`PFS`) has a confirmed,
      complete live estimate. Worth a live pass with a point that actually has
      `ELEV1000` and falls in a valid drainage-area range before trusting `LFS`
      results the way `PFS` is now trusted.
      **Done 2026-09-28** (addendum S6), with one part not possible as written: **no WA
      point has `ELEV1000`.** WA's ss-hydro computes 11 characteristics and `ELEV1000`
      is not among them. Asking for it returns the -999 "not found" sentinel. Its region,
      `GC1434` (Nooksack, SIR 2009-5170), also has no geometry in NSS. So `GC1434` cannot
      be reached through StreamStats at all; it needs a caller-supplied `ELEV1000`, and
      what elevation that means was not verified. The live pass therefore used the western
      WA regions StreamStats does serve:
      - Skookumchuck River near Vail (DRNAREA 39.9, PRECIP 71.07) gives M7D10Y of 16.2 cfs
        (`GC1556`, `0.15*DRNAREA^1.27`) and 12.0 cfs (`GC1557`,
        `0.000848*DRNAREA^1.17*PRECIP^1.23`). Both reproduce independently through
        `flowfreq.regression.nss.evaluate_expression` to NSS's 3-significant-figure
        rounding.
      - It found a real bug: LFS reports its standard error under code `SE`, not `ASEp`,
        so `standard_error_pct` was always `None` for low flows. Fixed, and the code used
        is recorded in `standard_error_code`.
      - Pinned by `TestLiveNSS::test_skookumchuck_low_flow_matches_its_equations`.

### Done — StreamStats module Phase 2: NSS flow-statistics estimation

`docs/STREAMSTATS_NSS_ADDENDUM.md` -- the live-verification pass Phase 1 required
before any Phase 2 code, same discipline as the original design doc. Every
NSS-related detail the `ff-idea02` PDF/py transcript guessed was wrong, in the same
pattern as its `ss-hydro` guess: `GET /nssservices/regions/{region}` carries no
statistic groups, `POST /nssservices/estimate` does not exist at any path tried, and
the real estimate call (`POST /nssservices/Scenarios/Estimate`) takes a **bare JSON
array** body -- wrapping it in `{"scenarioList": [...]}`, matching the API's own
machine-readable parameter name, is a confirmed-live `500` with no other diagnostic.
Found the real protocol by pulling the NSS docs viewer's Angular bundle apart to its
`assets/config.json`, which named the actual machine-readable spec
(`GET /nssservices/apiconfig`, a custom non-OpenAPI format) -- there is no published
OpenAPI/Swagger for `nssservices`, and the human-readable per-endpoint docs it
references live on an internal, sign-in-gated USGS GitLab, unreachable without
authentication.

Implemented in `flowfreq/streamstats.py`: `list_statistic_groups(region=None)`,
`estimate_flow_statistics(region, characteristics, statistic_group_codes=None)`,
`batch_estimate_flow_statistics(points, ...)`, plus `RegressionCitation`,
`FlowStatisticEstimate`, `RegionFlowEstimates` (typed, indexable, mirroring Phase 1's
own style). Every estimate carries NSS's own literal regression-equation string and a
resolved citation (title/author/DOI, via `GET /nssservices/citations`) -- closes the
design doc's FR-9 "uncited exponent" concern for real, not just in principle.

**NSS's own "a 200 is not an answer," found live and recorded as this module's
governing check for Phase 2**: submitting a wildly out-of-range parameter value
(`DRNAREA=999999` against a calibrated range of `[0.25, 3310]`) still returned a
computed, plausible-looking flow estimate with **no error and no warning at all** --
worse than Phase 1's `WarningMsg`-bearing failures, since there is no message to
detect. `estimate_flow_statistics` therefore validates every parameter against that
same request's own echoed-back `limits.min`/`limits.max` *before* submission
(`_fill_and_validate_regions`), so an invalid region is never even asked for, rather
than trusting a returned number after the fact. Confirmed live end to end
(`tests/test_streamstats.py::TestLiveNSS`): Goat Creek's real basin characteristics
reproduce a live NSS 50-percent-AEP peak-flow estimate with its real equation and
citation. 17 new tests (13 mocked + 1 live, plus dataclass serialization and
`_fill_and_validate_regions` unit tests).

### Done — StreamStats module Phase 1, and the `download_daily_flow` precedent fix

`flowfreq/streamstats.py`, per `docs/STREAMSTATS_MODULE_DESIGN.md` (merged 2026-09-10,
written from hands-on verification against the live service, not assumed from
documentation) -- watershed delineation and basin characteristics for a pour point:
`pourpoint` snap -> `ss-delineate` `delineate/sshydro` -> `ss-hydro` basin
characteristics, typed/indexable by StreamStats code with unit/description/service
message, stamped with provenance (service versions, request URLs, timestamp, server
used). Batch entry point mirrors `usgs.fetch_nwis_batch`'s `(results, errors)` shape --
one bad point never aborts the batch. Offline-capable JSON-file cache keyed on the
*requested* (not snapped) coordinate specifically so a cache hit needs no network call
at all, not even the snap -- a deliberate, documented departure from the design doc's
literal wording of the cache key, in favor of its own stronger NFR-5 requirement.
`tests/test_streamstats.py` (37 mocked tests + 3 `requires_network` live tests). PR #23.

**Confirmed live 2026-09-11, by the user running the `requires_network` tests by hand
(no network access in the Claude Code session that built this) -- and it caught two
real bugs the mocked suite couldn't, exactly as intended:**

- `snap_point`'s `output` field-name guess was wrong in a more specific way than
  missing keys: the real `pourpoint` response nests the snapped coordinate as a
  GeoJSON Point (`{"coordinates": [lon, lat]}`), not flat fields. Fixed with
  `_extract_point_lat_lon`, GeoJSON-order-aware, falling back to the old guessed keys
  defensively.
- A second bug the first fix then exposed: this implementation had added its own
  extra call to `ss-delineate/v1/delineate/features/{region}` (`x`/`y`/`crs` params,
  borrowed from the unverified ff-idea02 PDF/py transcript -- the same source whose
  `ss-hydro` endpoint guess was already known wrong) to fetch a watershed polygon for
  FR-3. Live, first a 422 (wrong param names -- fixed to `lat`/`lon`, matching the
  rest of the API), then a 200 carrying an unrelated zero-area Point feature, not a
  polygon. **Removed** rather than patched further -- it was never part of the design
  doc's own literally-verified protocol, which only ever called snap ->
  `delineate/sshydro` -> `ss-hydro`. The pipeline now matches that protocol exactly
  (two GETs, one POST, nothing more); FR-3's `WarningMsg` check is a recursive scan of
  the `sshydro` response instead of a polygon check.

**One known, honest gap, not a bug: `WatershedCharacteristics.polygon_geojson` is
always `None`.** No call in the verified, working protocol returns the watershed's
geometry. Obtaining it is an open question for a future session -- worth its own
live-verification pass (same discipline as everything else here) before attempting
again, not another guess from the PDF/py script. **Closed 2026-09-27**: the polygon
was in the `sshydro` response all along (`bcrequest.wsresp.featurecollection[0]`,
feature `globalwatershed`), found via USGS's own workflow notebook and verified live
on four basins before any code (design doc S10). Now returned and validated.

All three `requires_network` tests now pass against the real service, including both
Methow points reproducing the design doc's own published `DRNAREA`/`PRECPRIS10`/
`CANOPY_PCT` values exactly.

- [x] **The design doc's own "read this first" precedent, fixed alongside this work as
      it said to be.** `usgs.py::download_daily_flow` already sent an explicit date range
      (the earlier range-less-request defect was fixed before this session, commit
      `04f6b65`), but a code review of that fix found three more issues from combining it
      with the method's pre-existing hardcoded 30s timeout and local-clock default: the
      timeout was too tight for a full period-of-record request (now the method's own
      default) on a long-running, high-frequency site -- made configurable, default 60s;
      the default end date used local wall-clock (`date.today()`) rather than UTC, so a
      host clock behind UTC could silently request a narrower range than intended, the
      exact "silently wrong, not missing" failure mode this method's own docstring warns
      about; and `start_date`/`end_date` were sent to NWIS unvalidated -- a reversed or
      malformed range now raises before any request is made. 8 new tests in
      `tests/test_usgs.py`.
- [x] Pinned `pip-audit>=2.10.1,<3` in `ci.yml`, closing a gap a review found: it was the
      one unpinned tool in a CI file that otherwise pins everything (black<25, pytest<9,
      SHA-pinned actions) with an explicit stated rationale for why reproducibility
      matters -- an unpinned pip-audit release could silently change scan behavior between
      runs with no diff in this repo to review. (The step itself -- see the "Next --
      pip-audit in CI" item this replaces -- was already added to `ci.yml` sometime after
      that item was written and before this session found it; this session's contribution
      was closing the pin gap, not adding the step.)

### Done — transposing computed flows to an ungaged site, and QPPQ

Specified in **`docs/TRANSPOSITION_DESIGN.md`**. Read that before writing any more of it; as
with the Fortran engine, the parts with a silent failure mode are written down there rather
than left to be improvised.

Drainage-area-ratio transposition of a donor gage's fitted statistics to a nearby ungaged
site, `Q_t(p) = Q_d(p) * (A_t/A_d)**b(p)`, with `b(p)` from the applicable published USGS
regional regression rather than assumed to be 1.

- [x] **`RegressionExponents`** (exponents by AEP + a **mandatory** citation -- the class
      will not construct without one) and the interpolation/extrapolation of `b(p)` onto the
      AEPs actually computed. Interpolated against the normal deviate; `linear` by default
      because a reviewer can reproduce it with a calculator, `pchip` available and measured
      to differ by <0.005 over a normal published set. `flowfreq/transpose.py`.
- [x] **`transpose_frequency`** + guardrails (area-ratio band 0.5-1.5, raises by default
      outside it, `allow_out_of_range=True` transposes anyway and records the violation) +
      `to_markdown()` showing the per-quantile arithmetic and every caveat that applies.
- [x] **`regime.flow_duration_curve`** — duration statistics were a by-product of drawing a
      figure, at nine hardcoded percentiles. Now a function, at arbitrary percentiles;
      `Hydrograph.plot_flow_duration_curve` delegates to it, so the table it has always
      returned is one computation rather than two that can drift, pinned by a test asserting
      the plot's numbers and column names did not move.
- [x] **`transpose_duration`** — same machinery, exponents indexed by exceedance fraction,
      because `b` is not constant across a duration curve (near-linear at the wet end where
      area sets the flow, falling away at the dry end where geology does). A zero donor
      statistic comes back `NaN`, never `0`: scaling a zero would assert the target is dry
      rather than estimate it.
- [x] **`transpose_low_flow`** — separate function, tighter band (0.7-1.3), a **required**
      `hydrogeologic_setting` assertion recorded in provenance, refusal on a zero statistic
      or a donor that ever goes dry, and an optional BFI similarity screen.

**The category error is now unreachable, not just documented.** The design said a flood
exponent applied to 7Q10 is a category error and answered it with separate functions -- but
separate functions were not enough, because the *arguments* were still interchangeable and
nothing stopped a flood exponent set being passed to the low-flow one. `RegressionExponents`
now carries a `probability_kind` (`"aep"` / `"exceedance"` / `"non_exceedance"`) and each
function requires its own, so the mistake raises rather than returning a plausible wrong
number. Worth remembering as a pattern: a guardrail that lives only in which function you
call is not a guardrail if the inputs still fit both.

**The whole standard AEP set transposes, 0.995 down to 0.002** -- all fourteen points, not
just the flood range. The six more frequent than Q2 sit below anything a flood regression
publishes: in normal-deviate space AEP 0.995 is at `z = -2.58`, as far *below* the published
range as Q500 is above it. They are transposed on a clamped (endpoint-held) exponent, and
every one of them is flagged `extrapolated` in the provenance, logged once, and named in the
markdown caveats. That is the honest posture, not a fix: a flood regression's exponent
describes flood response, and at AEP 0.9 it is being asked about flows that are not floods.
Where duration-regression exponents exist for the frequent end, supplying them in the same
`RegressionExponents` removes the extrapolation entirely -- deliberately a matter of passing
more points, not calling a different function.

Two things the design doc argues and this list should not lose:

- **A flood exponent applied to 7Q10 is a category error**, not an approximation. Low flows
  are controlled by baseflow storage and geology; that is why published low-flow regressions
  carry a geology term and flood regressions do not. Hence `transpose_low_flow` as a
  separate function rather than a flag on the flood path.
- [x] **QPPQ (probability-preserving) transfer** of a whole daily series -- the right tool
      when the caller wants a time series rather than statistics. `flowfreq/qppq.py`: an
      invertible `FlowDurationCurve`, seasonal curve construction, donor ranking and
      goodness-of-fit, and a leave-one-out harness.
      `performance()` reports NSE, log-NSE, KGE and dry-end bias together, because an
      estimate that is exact through the freshet and wrong by 10x in September still scores
      NSE > 0.99. The counterintuitive result worth remembering: a coarse seasonal split of
      the curves makes snowmelt melt-timing phase error *worse*, not better, since a season
      block is longer than the offset it is meant to correct -- `estimate_donor_lag`
      recovering the offset from rank correlation is what actually works (log-NSE 0.273
      annual, 0.040 three-season, 0.605 monthly, 0.895 once lagged, on a 25-day imposed
      offset). See CHANGELOG.md's 0.7.0 entry.

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
      instead of a dict, or should `batch.py` adapt?) -- not done here. (Fixed since, in #46;
      the test is a plain passing test.)
- [x] `cli.py` (57, now larger with the `compare` subcommand) -- `tests/test_cli.py` covers
      `validate`/`benchmark` and `compare`'s extension-agnostic paths (CSV validation, the
      `ImportError`-to-`ClickException` mapping, exit codes) via a monkeypatched
      `compare_engines` so it runs the same with or without the built extension.
      `compare`'s real, Fortran-backed happy path is covered end-to-end in
      `tests/fortran_parity/test_live_cli_compare.py` instead (`requires_fortran`).

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
      more heavily than the three-analysis probe suggested. Decorator order has a test
      (`staticmethod` outermost); the reverse broke only on Python 3.9, which is no longer
      supported.

- [x] **Build the Fortran in one CI job.** New `fortran` job runs `make parity`: build,
      assert the extension imports, then run `tests/fortran_parity/`. The import assertion is
      the point — the parity tests call `importorskip`, so without it a failed build would
      skip all twelve live-vs-golden tests and report green. Verified on Linux/CPython 3.11:
      the extension builds from the vendored sources and the committed goldens are not
      drifted.

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
      hand) and against the app's `ffa_runner` override tests (now in `flowfreq-app`), not
      against vendored Fortran.
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
a standalone binary. It is not, and the `flowfreq/peakfqsa/` wrapper they built has since
been deleted (P2 above). Kept for the record; see `AGENT_BUILD_INSTRUCTIONS_Claude.md`.

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
