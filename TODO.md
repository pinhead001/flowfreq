# TODO — FlowFreq Hybrid 17C Implementation

> Long-range plan (national 17C + regional regression coverage, transposition,
> nonstationarity, future quantiles): see `docs/MASTER_ROADMAP.md`. Roadmap items
> in progress are tracked under "Roadmap Phase A and Wave 1" in Open Items below.

## Status
Last updated: 2026-09-29. Version **0.9.0** (`pyproject.toml`; tag `v0.9.0`, released
2026-09-27 in PR #75). Unreleased on `main` since: Water Data OGC API defaults for
instantaneous (#77) and daily (#79) values, the Water Data API key and 429/503 backoff
(#84), and the StreamStats watershed-polygon / NSS region-selection / Low-Flow work
(#78, #80, #82) -- see CHANGELOG.md's Unreleased section.
Python 3.11–3.14 (`requires-python >= 3.11`; CI matrix in `.github/workflows/ci.yml`).

Tests: measured 2026-09-29 on Linux / CPython 3.11 at `3736d8f` with
`PYTHONSAFEPATH=1 python -m pytest tests/` (the CI selection; `addopts` deselects
`requires_network`), extension absent: **2100 passed, 13 skipped, 54 deselected,
4 xfailed** (~5 min). Not measured with the extension built. v0.9.0's release note records
2 local failures with it on Windows -- the documented MSYS2 gfortran drift against the
Linux-built goldens (CLAUDE.md, Validation Status), not defects. `ci / Fortran parity` is the
authoritative parity check.

The strict xfails: only the four 2012 PeakfqSA manual comparisons
(`tests/validation/test_big_sandy.py`, a non-reproducible reference). Cains Coulee's
`skew_weighted` xfail and its `compare_engines` twin were resolved in #64; the
`batch.run_multi_site` xfail was fixed in #46.
Run the number, do not carry it forward -- it has been wrong in a commit message and a PR body
already.

**This file went stale once already.** Between 2026-09-25 and 2026-09-29 about forty PRs
(#42–#86) merged without updating it, so it listed as open nearly everything v0.9.0 had
delivered. Update the Open Items entry in the same PR that closes it.

Every P1, P2 and P3 item is done; the Done sections and "P3 — The `var_mom` port" hold the
history. Open work is the roadmap (below) plus the small items after it.
Fortran reference: **vendored** at `vendor/peakfqr/` (peakfq 8.1.0, CC0); the bridge builds
with `python build_fortran/build.py` (gfortran + meson), and CI builds and checks it with
`make parity`.

### Environment constraints -- read before starting work

These bit repeatedly and are not discoverable from the code:

- **Tag pushes and branch deletes return HTTP 403 from a Claude Code *web* session** --
  a credential boundary, not a transient failure: retrying and re-authenticating do not
  help. Branch pushes to the same repo succeed, and the egress proxy records no denial.
  The qualifier matters and this entry lacked it until 2026-09-07: **a CLI session on a
  developer machine is not affected.** Tags `v0.5.0`, `v0.6.0`, `v0.6.1` and `v0.7.0`
  were all created and pushed straight from one, as were several branch deletes. Read
  unqualified, this entry cost a session's worth of unnecessary hand-offs -- work was
  set aside as owner-action that the session could have done itself. The web-session
  finding is from `docs/PHASE1_RUNBOOK.md`, which had the distinction right.
- **NWIS (`nwis.waterdata.usgs.gov`/`waterservices.usgs.gov`) is blocked by the egress
  proxy, but that does not generalize to every USGS host.** `requires_network` tests
  against NWIS cannot run in a Claude Code session; use the committed fixtures under
  `tests/fixtures/`. **StreamStats (`streamstats.usgs.gov`) is a different host and is
  reachable** -- confirmed directly, 2026-09-11: both the raw endpoint probes and the
  actual `pytest tests/test_streamstats.py -m requires_network` suite ran and passed
  from inside a Claude Code session, no hand-off to the user needed. Read this
  unqualified the way the earlier tag-push/403 entry was misread once already: check
  per host, don't assume "USGS" is one block.
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

---

## Open Items (prioritised)

Audited 2026-09-29 against `main` at `3736d8f` and every remote branch. The Phase A /
Wave 1 deliverables that the 2026-09-25 audit listed as open are mostly merged; that
audit, with its live-verification findings, is kept under "Done — Roadmap Phase A and
Wave 1" below. `docs/MASTER_ROADMAP.md`'s status markers are still the 2026-09-25 ones
(e.g. "OGC backend is a stub", "no equations transcribed") and need the same refresh.

### Open pull requests

Each is a single commit on top of `main`; none is merged.

- [ ] **#81 `fix/ema-p3est-fixed-point`** -- the native EMA fixed point iterates as
      `p3est_ema` does; documents the FIXED-threshold residual as `emafitpr`'s own
      ill-conditioning, with a live `requires_fortran` test
      (`tests/fortran_parity/test_fixed_threshold_live.py`).
- [ ] **#83 `feat/regional-skew-at`** -- `regional_skew_at(lat, lon)`,
      `regional_skew_at_huc`, `regional_skew_for_site`, with HUC tables built by
      `tools/build_regional_skew_hucs.py`. The lookup half of #34.
- [ ] **#85 `fix/diel-variation-dst`** -- `diel_variation`'s `expected_obs` follows each
      local day's real length. Moves numbers `diel_variation` has already reported, so it
      belongs in a release allowed to move them (see Sub-daily below).

### Roadmap Phase A and Wave 1 (issues #27–#40) -- what is still open

- [ ] **#29 Water Data OGC API migration -- residue.** Peaks (#42, #49, #63, #70),
      instantaneous values (#45, #77) and daily values (#79) all read the OGC API by
      default, each after a live legacy-vs-OGC parity gate; `nwis-legacy` stays selectable.
      Left open:
      - Big Sandy has no WY1988–2002 peaks in the API, yet its `revision_note` mentions
        a revised WY2002 peak. Never checked against legacy.
      - The API carries no perception thresholds (`peak_since` is null); they remain
        analyst inputs, as with legacy.
- [ ] **#31 Published perception thresholds.** The `.psf` converter is done (#52,
      `flowfreq/psf_convert.py`); retrieving published thresholds is not.
- [ ] **#32 Regulation / urbanization screen.** Unchanged since 2026-09-25: only
      `peak_codes.classify_from_codes`. Open: GAGES-II/NID storage, NLCD impervious
      fraction, and refusing B17C on a regulated record without an override.
- [ ] **#33 National gage catalog.** `catalog.py` and `tools/build_gage_catalog.py` exist,
      and the #29 dependency for auto-discovery is now met, but no catalog has been built:
      `flowfreq/data/gage_attributes.csv` is still the 3-row seed.
- [ ] **#34 Regional skew -- residue.** WA/OR/ID rows are verified (PNW B-GLS, G = −0.07,
      MSE 0.18; #51, #69), and the silent −0.302 default is gone (#44: callers choose a
      regional skew, `station_skew_only=True`, or `use_default_skew=True`). Left open:
      - The location lookup -- PR #83 above.
      - Wiring the table into `Bulletin17C` / `run_ffa` so a site gets its skew without
        the caller looking it up.
      - **MT** stays `pending` on purpose: no B-WLS/B-GLS study covers it, and USGS MT
        practice uses the B17B Plate I map (the CSV row has the full reasoning). The
        roadmap's Montana Bayesian regional skew study is the real fix.
      - **ID Snake River Plain** has no regional skew by the source report's own statement.
- [ ] **#35 Regression library -- residue.** `tools/snapshot_nss.py` and the exact NSS
      equation parser landed (#47, #58, #60); snapshots are under `data/nss_snapshots/`.
      Left open:
      - `streamstats.estimate_flow_statistics` still keeps only `ASEp`/`SE`
        (`_STANDARD_ERROR_CODES`); NSS's raw `SEp` is kept by the snapshot tool
        (`nss_sep`) but not surfaced on `FlowStatisticEstimate`.
- [ ] **#36 Future-flow change factors.** Unchanged: `future_flow.py` ships no factor
      values and `flowfreq/data/future/` does not exist. Roadmap §6.3.1 (national sets
      first) and the module docstring (per-state, wave by wave) still disagree.
- [ ] **#37 WA -- `partial`.** 32 equations in `data/regression/WA.json` (#55). Region 4
      carries no covariance: Table 7's matrix fits only a P/10 basis the report never
      states, and #65 showed the report itself applies it to raw P. Its intervals fall
      back to `sep_log`. The WA section of `docs/FUTURE_FLOW_GUIDANCE.md` is not written.
- [ ] **#38 OR -- `partial`.** Western OR (SIR 2005-5116, #57), eastern OR (OFR SW 06-001,
      #71) and the 2A/2B transition blend (#67) are in. No covariance is published, so
      intervals use `sep_log`. Future-flow guidance section not written.
- [ ] **#39 ID -- `verified`** (SIR 2016-5083, #48, #59). Open: future-flow guidance.
- [ ] **#40 MT -- `verified`** (SIR 2015-5019-F, #56; channel-width equations and method
      weighting from SIR 2020-5142, #73). Open: the regional skew (above) and future-flow
      guidance.
- [ ] **Wave 1 definition of done** (`.github/ISSUE_TEMPLATE/state-rollout.md`,
      roadmap "Definition of done"): region polygons with point-in-polygon selection, the
      at-site B17C appendix re-run, the cross-border check and transposition LOOCV have
      not been done for any Wave 1 state.

### Small open items

Every small item listed on 2026-09-25 is done -- see "Done — small items closed after
2026-09-25". None new has been recorded since.

### Sub-daily metrics follow-ups

`flowfreq/subdaily.py` landed (see `docs/SUBDAILY_METRICS_DESIGN.md`, and the
Unreleased CHANGELOG entry). Three things it left open, all recorded rather
than fixed:

- [ ] **(PR #85 open)** **`diel_variation`'s `expected_obs` is a fixed 1440 minutes / median step,
      so it mismarks the two daylight-saving transition days of every year** --
      the spring day (23 local hours) reads as incomplete at any
      `min_completeness_frac` above ~0.96, and the autumn day (25 hours) reads
      as more than complete. `daily_extreme_timing` and `ramping_rates` use each
      day's actual local length and do not have this. Left alone deliberately:
      changing it moves numbers `diel_variation` has already reported. Fix it
      together with the next release that is allowed to move those, not as a
      drive-by.
- [ ] **Paired stage-and-discharge retrieval.** `download_instantaneous_flow`
      and `download_instantaneous_stage` are separate calls returning separate
      frames. Anyone comparing a cfs/hr limit against a ft/hr one joins them on
      the index by hand. A single call retrieving both parameters onto one index
      is the natural next step and is not built.
- [ ] **No validation against a published ramping-rate or peak-timing figure for
      a real gage.** Every number in `subdaily.py` is checked against synthetic
      series with known analytic answers and against internal consistency, which
      catches the circular-mean and local-`dt` errors the module exists to avoid
      -- and is *not* the same as reproducing a published hydropeaking statistic.
      That needs live NWIS, blocked from a Claude Code session (see the
      environment constraints above). Treat the metrics as verified in arithmetic
      and unverified against the literature until someone runs a real gage
      through them.

### Next — StreamStats Phase 2 follow-ups

Phase 2 itself (NSS flow-statistics estimation) is done -- see "Done -- StreamStats
module Phase 2" below. Two things it found are still open:

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
      instead of a dict, or should `batch.py` adapt?) -- not done here.
- [x] `cli.py` (57, now larger with the `compare` subcommand) -- `tests/test_cli.py` covers
      `validate`/`benchmark` and `compare`'s extension-agnostic paths (CSV validation, the
      `ImportError`-to-`ClickException` mapping, exit codes) via a monkeypatched
      `compare_engines` so it runs the same with or without the built extension.
      `compare`'s real, Fortran-backed happy path is covered end-to-end in
      `tests/fortran_parity/test_live_cli_compare.py` instead (`requires_fortran`).

### Downstream

- [ ] **Bump the app's pin to v0.9.0.** `pinhead001/flowfreq-app`'s `requirements.txt`
      was `flowfreq @ git+https://github.com/pinhead001/flowfreq@v0.7.0` when last checked
      (2026-09-25); not re-checked since, because that repo was outside this session's
      scope. v0.8.0 was never picked up, so the bump goes straight to v0.9.0, and **FFA
      output will not be identical**. v0.9.0 is breaking in three ways the app reaches
      through `workflow.run_ffa`:
      - no silent −0.302 regional skew: `run_ffa` needs a regional skew,
        `station_skew_only=True`, or `use_default_skew=True` (which reproduces old results);
      - `download_peak_flow` reads the Water Data OGC API by default (`peak_date` is a UTC
        date; water years are unchanged), and partial-date peaks (e.g. Big Sandy's three
        historic ones) are now kept;
      - NWIS peak codes are applied by default (`apply_peak_codes=False` to opt out).
      The native EMA fixes in v0.9.0 (zero-flow rows, perception thresholds, exact LP3
      quantiles, the B17B skew-MSE switch, near-zero-skew bounds) also move numbers at
      sites that hit them. Verify as before: `pip show flowfreq` confirms the tag, and the
      app's suite runs before and after, with any changed expectation explained.

      History: bumped through v0.5.0, v0.6.0, v0.6.1 and v0.7.0, each verified rather than
      assumed. The app's `plot_peak_timeseries` copy was deleted in the v0.6.1 switch to
      `plot_peak_flows_with_thresholds`. See `flowfreq-app/CLAUDE.md`'s "one edit,
      deliberate" note.

### P3 — The `var_mom` port, now complete

**Resolved 2026-09-26: no `xfail` remains from this item.** Cains Coulee's `skew_weighted`
rung passes, at 6e-6 (quantiles to 0.0012%). The cause was not anything in the `var_mom` port:
`emafit.f:707-711` switches the at-site skew MSE to the plain Bulletin 17B `mseg(n, G)` when
MGBT computes the threshold and finds low outliers, and the native engine used ADJE
unconditionally. The "unexplained ~3x `as_G_mse` discrepancy" recorded below, and in the
Cains Coulee rows of the tables that follow, is that switch. Those passages are kept as history;
read them with this in mind.

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
      source label). A fourth, later-found gap: the app also let the user toggle the y-axis
      between linear and log per-plot (a real sidebar control, "Log: Peak Flow Time Series"),
      while the library function was always log. Added `yscale: str = "log"` (same convention
      `plot_frequency_curve` already uses), gating the power-of-10 tick formatting behind
      `yscale == "log"`. Verified live, not assumed: the quantile-line values this function
      computes analytically from `lp3_params` via `_lp3_quantiles` match `run_ffa`'s actual
      `quantile_df["Flow (cfs)"]` to 0.0% on both an uncensored site (Powder River, all nine of
      the app's quantile options including 1.5/200/500-yr) and a censored one (Big Sandy);
      `core.log_pearson3_cdf` (used for the max-peak annotation) matches the app's own
      `scipy.stats.pearson3.cdf`-based formula to ~1e-15. App-side switch (calling this instead
      of `plot_peak_timeseries`, then deleting the app's copy) landed in the v0.6.1 pin bump --
      done, see that repo's `TODO.md` for the session record.

- [x] **`FrequencyComparator` compares every parameter by percent difference.** That was the
      wrong metric for skew, which legitimately crosses zero: Big Sandy's reference at-site
      skew is 0.0066, so an absolute gap of 0.016 read as 249% and dominated `max_diff_pct`. Done
      in commit `dc563c5` (skews now go to their own `skew_diffs` dict, compared in skew units
      against `skew_tolerance_abs` and excluded from `max_diff_pct`) but this checkbox was never
      updated to say so. Verified still correct and covered by
      `tests/validation/test_comparisons.py::TestSkewComparedInSkewUnits`, including the exact
      Big Sandy numbers above (`test_a_near_zero_skew_no_longer_dominates_max_diff`).

- [x] **`origin/dev` deleted.** The read-and-judge pass found only `extra_curves` worth
      keeping (`flowfreq/setup.py` was stale packaging contradicting `pyproject.toml`;
      everything else was superseded or older than main) and it was already ported to main
      -- see the commit "Port extra_curves from dev". Confirmed gone from `git ls-remote
      --heads` on 2026-09-08; tip was `86cb147` if it's ever needed back.

### Blocked

- [x] **Tag pushes and branch deletes: resolved, and the constraint was narrower than
      recorded.** The HTTP 403 is specific to a Claude Code *web* session; a CLI session on
      a developer machine has neither limitation. Stated unqualified, this entry sent work
      to the owner that a session could have done -- see the Environment constraints note
      above for the corrected wording.

      Everything this item tracked is now done: `v0.4.0` was tagged from a local clone,
      `v0.5.0` through `v0.7.0` were tagged and pushed directly from a CLI session,
      `typecheck` and `tests-engine-report` are gone, and `parity-12363000` was deleted on
      2026-09-07 (tip `576d92a`, recoverable by SHA; it was fully contained in `main`, so
      nothing was lost). It was on `flowfreq` itself, not `hydrolib` as this item said for
      most of its life. Historically this class of block also hit `v0.2.0` and an
      `archive/dev-2026-02` tag.

---

## Done

### Done — Roadmap Phase A and Wave 1 (the 2026-09-25 audit, kept as history)

Kept verbatim from the 2026-09-25 audit: the "open" wording below is as of that date, and
the live-verification findings (OGC API site-ID prefix, `parameter_code` filter, paging,
UTC vs local `peak_date`, qualifier-token mapping, IV limits) still hold. What closed each
item: #29 peaks #42/#49/#63/#70, IV #45/#77, daily #79, API key #84; #30 #52/#68/#72;
#31 converter #52; #34 #44/#51/#69; #35 #47/#58/#60; #37 #55/#65; #38 #57/#67/#71;
#39 #48/#59; #40 #56/#73. Remaining work is in Open Items above.


`docs/MASTER_ROADMAP.md` requires a TODO.md entry for every roadmap item that has
started. PR #26 (992a1c2) started all of these. It is **scaffolding only**: the
types, loaders and tests are real (93 tests across the seven modules), but no
endpoint has been live-verified and no published coefficient, skew or change
factor has been transcribed. The roadmap's `[~]` means "scaffolded", not
"partly delivered". For every item below, the substantive deliverable is still
open. Audited 2026-09-25.

**Phase A — data foundation (#27, epic)**

- [ ] **#29 Peak-data backend adapter / Water Data OGC API migration.**
      `flowfreq/peak_sources.py` has the `PeakDataBackend` protocol,
      `validate_peak_frame`, and `LegacyNwisBackend` (which wraps `USGSgage`).
      `WaterDataApiBackend.fetch_peaks` always raises `NotImplementedError`, and
      `DEFAULT_BACKEND` stays `nwis-legacy`. `USGSgage` itself does not route
      through the adapter. `api.waterdata.usgs.gov` **is reachable** from a Claude
      Code session. The egress block on the legacy NWIS hosts does not apply to it.

      **Live-verified 2026-09-25 (read-only):** go on `fetch_peaks`. Every
      fixture peak matches exactly: Big Sandy 47/47 including the three historic
      peaks, and Orestimba 82/82 including the 12 zeros. Use `/ogcapi/v1/`; `v0`
      returns identical bytes, and every link points at v1. No API key is needed.
      No rate-limit headers were returned. CloudFront caches responses for up to
      an hour. Implementation requirements, each of which fails silently if
      missed:
      - **Site ID needs the `USGS-` prefix.** `monitoring_location_id=03606500`
        returns 200 with 0 features. Build `USGS-{site_no}`, and raise on an empty
        result.
      - **Filter `parameter_code=00060`.** Gage-height (00065) peaks are separate
        rows in the same collection, and without the filter every water year
        appears twice, so `validate_peak_frame` raises.
      - **Page explicitly.** `limit` defaults to **10** and caps at 50,000. Follow
        `links[rel=next]` (a cursor; there is no `numberMatched`). Order is
        arbitrary without `sortby=water_year`.
      - **Take `water_year` from its own field**, never from `time`. `time` is
        the *UTC* date when the time of day is known, while legacy `peak_dt` is
        the local date, so evening peaks differ by a day. For example, the WY2023
        2100 cfs peak shows as 2023-01-04, but it was 2023-01-03 20:45 CST. A Sep 30
        evening peak would land in the wrong water year.
      - **Unknown day or month is a placeholder date, not null.** For example,
        Big Sandy 1897 appears as `1897-03-01` with the `DAYUNKNOWN` qualifier.
        Use the null `month`/`day` fields or the `DAYUNKNOWN`/`MONTHUNKNOWN`
        tokens to detect it.
      - **Translate `qualifier` tokens to NWIS codes before `peak_codes`.** The
        tokens are words. As a list, `peak_interval(100, ["LESSTHAN"])` returns
        systematic, so censoring is silently lost. As the CSV string
        `"DAYUNKNOWN,HISTORIC"`, `parse_codes` splits it into characters that
        include O and C, so the historic peak is **removed**. The following
        mapping was verified against vendored WATSTORE files: MAXDAILYMEAN=1,
        ESTIMATED=2, DAMFAILURE=3, LESSTHAN=4, UNKNOWNREGULATION=5, REGULATED=6,
        HISTORIC=7, GREATERTHAN=8, EVENT=9, URBAN=C, OPPORTUNISTIC=O,
        REVISED=R. No tokens were observed for A or F. Keep unknown tokens and
        log them; never drop them.

      Also found:
      - [x] **Legacy partial-date bug: fixed.** `download_peak_flow` ran
            `pd.to_datetime(peak_dt, errors="coerce")`. NWIS RDB encodes an
            unknown day or month as `00` (`1897-03-00`), which coerced to NaT, so
            the row was dropped. That lost Big Sandy's three historic peaks and
            every partial-date year. Now `usgs._parse_peak_dt` follows peakfq
            8.1.0's reader (`vendor/peakfqr/R/DataReaderFunctions_shinyapp.R`):
            month 00 becomes January and day 00 becomes the 1st. Fixed in the
            same pass: `peak_cd` is read as a string, because an all-numeric
            code column was being inferred as float (`"7"` became `"7.0"`).
            Tested against a synthetic RDB fixture (`PEAK_PARTIAL_DATES`). The
            `-00` encoding is confirmed by peakfq's reader, which handles it
            explicitly, but no real NWIS response has been captured, because
            NWIS is blocked from Claude Code sessions.
      - [ ] Decide how `peak_date` is represented (UTC vs. local date, and
            placeholder dates) before `DEFAULT_BACKEND` switches.
      - Big Sandy has no peaks for WY1988–2002 in the API, yet the site's
        `revision_note` mentions a revised WY2002 peak. This has not been checked
        against legacy.
      - The API carries no perception thresholds (`peak_since` is null), so they
        remain analyst inputs, as with legacy.

      **Daily:** no blocker. Filter `statistic_id=00003`, and page or `sortby=time`.
      Omitting `time` returns the full period of record, which makes the
      `DEFAULT_START_DATE` workaround unnecessary. Ice days are numeric with
      `["ESTIMATED","ICE"]`, where legacy coerced "Ice" text to NaN.

      **Continuous (IV):** blocked on three design decisions:
      - [ ] **Multi-sensor.** Series are identified by a 32-character hex
            `time_series_id` and returned *interleaved* at identical timestamps.
            The legacy DD number that `ts_id` means is not in the API. Discovery
            has to go through the `time-series-metadata` collection
            (`sublocation_identifier`, `primary`; at 03612600 both HEADWATER and
            TAILWATER are Primary). `_download_instantaneous`'s
            `index.duplicated(keep="first")` would silently merge two sensors,
            so it must refuse multiple series instead.
      - [ ] **No local time.** Timestamps are RFC 3339 UTC, and the API returns
            no `tz_cd`. `datetime_local`/`tz_cd` would have to be derived from
            the monitoring location's `time_zone_abbreviation` and
            `uses_daylight_savings`, or be dropped from the contract.
      - [ ] **Limits.** The `time` interval is capped at 1100 days (a hard 400)
            and pages at 50,000 rows. The `chunk_years=1` default fits. An empty
            window returns 200 with 0 features, not 400, so
            `_is_no_data_response` becomes dead code. With no `time`, the
            service returns the last year, not the period of record.

      **Monitoring locations:** `drainage_area`, a 12-digit `hydrologic_unit_code`
      (take the first 8 digits for HUC8), `state_code` as a FIPS code (needs a
      lookup to an abbreviation), and lon/lat in WGS84. `items/03606500` without
      the prefix returns 404.

      Next: implement `fetch_peaks` to the requirements above, confirm the
      legacy date bug, then write the parity test before `DEFAULT_BACKEND`
      switches.
- [ ] **#30 Peak qualification codes → B17C treatment.** `flowfreq/peak_codes.py`
      (`parse_codes`, `peak_interval` following peakfq's `siteQT` for codes
      4/8/3/O/6/C) is implemented and tested, but the analysis path does not use it.
      Only `catalog.py` imports it. `PeakRecord` carries no code flags, and
      `usgs.py` just copies `peak_cd` into a `qualification_code` string. Open: parse
      codes at ingestion, and derive EMA intervals from them.
- [ ] **#31 PeakFQ `.psf` reader.** `flowfreq/psf.py` (`parse_psf`, `read_psf`,
      `StationSpec`) is implemented and checked against the vendored WY/MT file.
      Open: the converter from a parsed station to `Bulletin17C`/`EMAParameters`
      arguments, which `psf.py`'s own docstring calls "the open half", and the
      retrieval of published perception thresholds.
- [ ] **#32 Regulation / urbanization screen.** Only step 1 exists:
      `peak_codes.classify_from_codes`, which by design has no `REFERENCE` class.
      Open: GAGES-II/NID storage, NLCD impervious fraction, and refusing B17C on a
      regulated record without an override.
- [ ] **#33 National gage catalog.** `flowfreq/catalog.py` (schema and loader) and
      `tools/build_gage_catalog.py` exist. The tool needs network access and a
      hand-supplied `--sites` list; auto-discovery waits on #29. No catalog has been
      built, so `flowfreq/data/gage_attributes.csv` is still the 3-row seed.
- [ ] **#34 Regional skew table and lookup.** `flowfreq/regional_skew.py`
      (`load_table`, `regional_skew_for`, which never returns a `pending` row) is
      implemented. `flowfreq/data/regional_skew.csv` has four Wave 1 rows, all
      `pending`. Open:
      - Populate the rows. MT can be cross-checked against GenSkew/SkewSE in
        `vendor/peakfqr/inst/testdata/wymt_ffa_2022A.psf`.
      - Add a lat/lon lookup (roadmap §1.3's `regional_skew_at`).
      - Wire it into `Bulletin17C`.
      - **Resolve the conflict with the silent default.** `workflow.B17C_DEFAULT_SKEW =
        -0.302` is the default argument of `run_ffa`/`compare_engines` and is also
        the CLI fallback. Roadmap §1.3 says to raise rather than default, and so
        does `regional_skew.py`'s own docstring. Changing the default is
        user-visible, so it needs a decision, not a drive-by fix.

**Wave 1 prerequisites and states (#28, epic)**

- [ ] **#35 Regression equation schema and offline evaluator.**
      `flowfreq/regression/` (`RegressionEquation`, `Citation`, `Variable`,
      `evaluate` with range checks and prediction intervals, `evaluate_weighted`,
      and the table of 56 jurisdictions and their waves) is implemented and tested,
      but only against synthetic equations.
      - [ ] **`tools/snapshot_nss.py` does not exist** in any commit or branch.
        It is still referenced by `regression/library.py`, the roadmap, and all
        four Wave 1 JSON files. NSS scenario templates carry **no equation
        field**: equation strings appear only in `Scenarios/Estimate` results. So
        the tool must:
        1. POST synthetic in-range inputs for every region.
        2. Parse strings such as `3.846*DRNAREA^0.745*10^(0.032*PRECPRIS10)/...`
           into the intercept and coefficients.

        First check whether NSS's `RegressionRegions` or `apiconfig` resources
        expose equations directly. Reuse `streamstats.py`: `list_regions`,
        `list_statistic_groups`, the template fetch and estimate loop in
        `estimate_flow_statistics`, `_fetch_citations`, and
        `_request_with_backoff`.
      - [ ] NSS cannot supply everything the schema needs. `Citation.table` (the
        report table number), the covariance, `model_error_variance` and `n_sites`
        have to come from the published report. A snapshot alone therefore cannot
        produce a `verified` state file.
      - [ ] `estimate_flow_statistics` discards NSS's `sep` field, which appears to
        be the log10 standard error of prediction and maps to `sep_log`. Keep it.
- [ ] **#36 Future-flow `ChangeFactorSet` framework.** `flowfreq/future_flow.py`
      (`ChangeFactorSet`, `select_factor_set` with state-over-national precedence,
      and `apply_change_factors`, which keeps both current and future columns) is
      implemented. **It ships no factor values**, and `flowfreq/data/future/` does
      not exist. Roadmap §6.3.1 wants the national sets (HEC-17, NCHRP 15-61)
      first, whereas the module docstring says sets arrive per-state wave by wave.
      Reconcile the two.
- [ ] **#37 WA (pilot).** Everything is pending: `data/regression/WA.json`, the
      skew row, and the WA section of `docs/FUTURE_FLOW_GUIDANCE.md`. The source
      report is known and already live in NSS: Mastin, Konrad, Veilleux & Tecca
      (2016), SIR 2016-5118, with 4 peak-flow regions (GC1750–GC1753).
      `tests/fixtures/streamstats_responses.py` has two live WA equation strings,
      which make a ready parser test for `snapshot_nss.py`. Follow the
      definition of done in `.github/ISSUE_TEMPLATE/state-rollout.md`.
- [ ] **#38 OR**, **#39 ID**, **#40 MT.** Everything is pending, as for WA. The
      current peak-flow (PFS) reports were identified 2026-09-25 from live NSS
      (region IDs: ID=16, MT=30, OR=41, WA=51) and USGS pubs pages. Nothing has been
      transcribed yet.
      - **ID (#39): SIR 2016-5083** (Wood, Fosness, Skinner & Veilleux, 2016,
        doi 10.3133/sir20165083; NSS citation 39). This is the best schema fit,
        because every field is published:
        - 6 regions (GC1735–GC1740) and 11 AEPs (80 % to 0.2 %).
        - Table 4 has the equations, MEV, AVP, SEP and per-region n. Table 5 has
          the variable ranges. Table A5 has the covariance matrices.
        - A quirk: NSS writes a `(MINBELEV/1000)^b` term. `Variable.scale`
          (`scale=0.001`) now stores it as published; the limits stay in feet.
        - Appendix B gives a Pacific Northwest regional skew of −0.07, relevant to
          #34. Its MSE is not yet confirmed.
      - **MT (#40): SIR 2015-5019-F** (Sando, Sando, McCarthy & Dutton, 2016,
        ver. 1.1 2018; doi 10.3133/sir20155019F; NSS citation 88).
        - 8 regions and 10 AEPs.
        - Tables: equations in 1–4, ranges in 3, covariance in 1–5, and the
          appendix workbook `sir20155019F_tables.xlsx`.
        - The drainage variable is `CONTDA`, not `DRNAREA`. The W region uses
          `(FOREST+1)^b` (`log10_plus1`). The NW region is WLS, not GLS.
        - NSS also serves **SIR 2020-5142** (Chase et al., 2021; channel-width
          equations; 24 regions = 8 × 3 methods; citation 163). That report
          weights by SEP and cross-correlation, which `evaluate_weighted` does
          not model.
        - Only chapter F of SIR 2015-5019 is the ungaged-site peak-flow chapter.
      - **OR (#38): SIR 2005-5116** (Cooper, 2005; NSS citation 118). It covers
        **western Oregon only**.
        - 3 regions: 1 coastal, 2A (≥ 3,000 ft) and 2B (< 3,000 ft).
        - 7 AEPs, with no 0.5 %.
        - Equations are in Tables 10–12.
        - Whether the covariance or model error is published: UNCONFIRMED.
        - Eastern Oregon has no NSS equations. OWRD Open File Report SW 06-001
          (Cooper, 2006) exists outside NSS.
        - NSS gates the regions on the `ORREG2`/`ELEV` selector parameters, which
          have no limits. Cooper's eq. 7 blends 2A and 2B by interpolating on
          elevation, not by area.
      - **Do not reuse hydrolib's regression branch**
        (`origin/claude/regression-roi-comparison-6BgBi`, unmerged). Its
        WA/OR/ID/MT coefficients are self-described "representative illustrative
        coefficients", and its citations are wrong: the OR one is a low-flow
        report, and the MT one is superseded and misnumbered. Its
        `RegressionTable` API and blank-template writer are usable as design
        references only.


### Done — small items closed after 2026-09-25

#13 in #54 (`flowfreq.donor_similarity`); `run_multi_site`, `make clean` on Windows and
`as_G_PRL_o` (`FrequencyResults.pseudo_record_length`) in #46; CLI `compare --historical`
/ `--threshold` in #50. The entries below are as recorded on 2026-09-25.

Each verified still present on 2026-09-25. Where the finding was first recorded
elsewhere in this file, that entry keeps the history; this list is the open tracker.

- [x] **#13 Screen donors on basin similarity, not drainage area alone.** Roadmap §5.1
      (multiple-donor weighting and similarity ranking). The issue was filed against
      functions that are not in this repo (`screen_donor_ratios`, `assemble_target_fdc`,
      `PREFERRED_RATIO_BAND`); flowfreq's own area screen is `transpose._check_areas`
      (area ratio only), and `qppq.rank_donors` ranks by concurrent-flow correlation, which
      needs a record at the target. No predictor-similarity screen exists in either.
- [x] **`batch.run_multi_site` cannot analyze real `fetch_nwis_batch` output.** Dicts
      in, `PeakRecord` expected, every site swallowed into `{"error": ...}`. Pinned by
      `tests/test_batch.py::TestAnalyzeSites::test_real_fetch_output_shape_is_analyzable`,
      `xfail(strict=True)`. The fix needs a decision on which side adapts; see "Modules with
      no tests".
- [x] **`make clean` misses the Windows extension.** `Makefile` `clean` removes
      `flowfreq/peakfqr/_emafort*.so` only, not `_emafort.cp3xx-win_amd64.pyd` or its four
      MinGW DLLs, so `make clean-verify` on a Windows machine that has built the extension
      still tests with it present. See "Done — the Fortran as a selectable engine".
- [x] **CLI `compare` has no historical-peak or perception-threshold flags**
      (`flowfreq/cli.py`, said so in the command's docstring/`--help`). Such records go
      through `workflow.compare_engines` directly.
- [x] **`as_G_PRL_o` (pseudo effective record length) is not surfaced.** Only
      `validation.reference.ReferenceResult.pseudo_record_length` carries it (from the golden
      file or a live `emafitpr`); neither engine puts it in `FrequencyResults`. See P3's
      `VAR_EMAB` entry.

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

## Resolved Questions

- PeakfqSA: not a standalone binary. The vendored `vendor/peakfqr/src/` Fortran is the
  reference; `flowfreq/peakfqsa/` was deleted (P2 above).
- Reference material: vendored into `vendor/peakfqr/` (CC0). No external workspace needed.
- Big Sandy 2012 manual values: not reproducible by peakfq 8.1.0 — the HWN skew weighting
  postdates the manual and diverges by design on censored records.
- MGBT: verified line-by-line against the Fortran (`GGBCRITP`/`FP_TNC_CDF`), validated on
  Orestimba Creek (USGS 11274500).
- Python: 3.11–3.14, all four in the CI matrix; `requires-python >= 3.11`.
- Regional skew default: `workflow.B17C_DEFAULT_SKEW = -0.302`, SE 0.55 (MSE 0.3025), is
  still the silent default in `run_ffa`/`compare_engines` and the CLI. Whether to keep it
  is an open decision, not a resolved one; see #34 above and roadmap §1.3.
- FrequencyAnalyzer API: added `validate()` and `to_comparison_dict()` to the `Bulletin17C`
  facade.
