# TODO — FlowFreq Hybrid 17C Implementation

> Long-range plan (national 17C + regional regression coverage, transposition,
> nonstationarity, future quantiles): see `docs/MASTER_ROADMAP.md`. Roadmap items
> in progress are tracked under "Roadmap Phase A, Wave 1 and Wave 2" in Open Items below.

## Status
Last updated: 2026-11-21. Version **0.28.0** is released (tag `v0.28.0`). `main` carries
**0.29.0.dev0**; 0.28.x patches come from `release/0.28.x` (CLAUDE.md, Branches and Releases).
flowfreq-app pins v0.28.0 (app #55, merged 2026-11-21); no app PR is open.
51 jurisdictions are `verified`, and 5 have no equations.

**StreamStats blocked this machine's address (HTTP 403) from 2026-11-09 14:40 to 2026-11-16, on purpose: issue
#365, now closed.** A StreamStats developer traced an NSS `bylocation` slowdown to `tools/border_qa.py` with three
workers: NSS keeps a query running after the client disconnects, and the client's 60-second timeout and retries
ran one slow basin's query up to nine times at once. The client side is fixed (no retry on a `bylocation`
timeout, a 200 s timeout, one call at a time with a pause; released in 0.27.0). The owner replied on 2026-11-15,
and StreamStats lifted the block on 2026-11-16, adding that their web services are getting rate limiting. One
test request on 2026-11-17 (a snap at TN 03465500) returned HTTP 200. The owner's reply promised to wait for
StreamStats' go-ahead before any larger run, so `border_qa` and the other live checks below still wait for the
owner's word. Run them with one worker.

### Review of 2026-11-19

Three read-only reviews found about 37 things. #396 fixed most; the follow-up PR did the eight recommendations
left over (transposition notes, four data-retrieval leads, report consistency and the `sections.py` split, the
slow test, the old engine's deprecation, an opt-in coverage requirement, historical notes on stale documents,
CI cache). CHANGELOG.md, Unreleased, has the detail. Still open from it:

1. **Null-value peaks: no change needed** (checked 2026-11-21). peakfq's `siteQT` drops a row with no
   `peak_va` itself ("Remove GH only peaks", `vendor/peakfqr/R/readInputs.R`) before filling in every year of
   the range, so the year is missing data exactly as flowfreq treats it: censored below a perception threshold
   that covers it, otherwise a gap. The warning in `peak_sources.py` said otherwise and now says this.
2. Done 2026-11-21 (#405): the deprecated engine (`engine.py`, `batch.py`, `plots.py`) is removed for 0.29.
   It was not a Bulletin 17B implementation, only B17B's first step (a moments fit with station skew).
   `core.log_pearson3_ppf` stays (it is correct; its docstring says `p` is the non-exceedance probability).
   The B17B procedure that replaces the idea is item 7.
3. Closed 2026-11-21, no rename: the names that differ name different things (`target_name` is the ungaged
   site, `annual_table` holds minima or durations, `detail` defaults follow whether the parts stay in the
   report). The `reporting/sections.py` module docstring says so.
4. Done 2026-11-21 (#403): a report fits each gage once. `gaged` memoises the fit on its options (it was 8
   fits counting nested sections and renders); `gaged.clear_fit_cache()` forgets them. A cache across sections was not
   added: the builders mutate the options they are given and a fit is about a second.
5. Done 2026-11-21 (#401): the chat tools take perception thresholds, peaks typed by hand, a regional skew by
   state, `ungaged` `aeps`, `weighting_space` and `min_coverage_percent`, and document a summary location's
   `basin`.
6. Done 2026-11-21 (#402): the eight class-scoped fixtures are staticmethods, and `PytestRemovedIn10Warning`
   is an error, so the `pytest<10` pin can move when pytest 10 is out.
7. **Bulletin 17B procedure, for comparing B17C with the superseded B17B. Designed and planned, not started**
   (owner, 2026-11-21: do not start building yet).
   - Design (approved): `docs/superpowers/specs/2026-11-21-bulletin17b-design.md`. Plan, 8 tasks:
     `docs/superpowers/plans/2026-11-21-bulletin17b.md` (#408). Execution: subagent-driven (owner's choice),
     one implementer per task and a fresh reviewer before the next.
   - Owner's decisions: a new `Bulletin17B` class (`Bulletin17C(method=MOM)` unchanged); B17B is computed and
     reported **only when requested** (a `Bulletin17B` call or a `b17b_ffa` section); the B17B side uses the
     B17C side's regional skew by default, with `generalized_skew` for the old study's Plate I value (three-
     column comparison).
   - The transcribed bulletin (Appendix 12 examples 1-4, the Appendix 6 Big Sandy example, Appendix 4 K_N,
     Tables 1, 9-1, 11-1, the procedure with page labels) is on branch `wip/b17b-transcription` under
     `tests/fixtures/b17b/`; plan Task 1 starts from it (its source path in the plan was a scratch folder).
   - Before starting: a fresh `spec/` or feature branch from `main`, then Task 1.
8. Done 2026-11-21 (#407): agent ops run on the shared `agentops` package (`D:/a/agentops`, `pip install -e` into
   `.venv`); `tools/agentops/` is gone, the config is `config/agents.toml`, `run.bat` runs it, the push guards
   are tighter. The features it needs came from pinhead001/agentops#6 (merged 2026-11-21). All 12 finished
   cards are archived.

### Next steps (2026-11-15)

**The batches in `docs/VERIFICATION_PLAN.md` are done. Read its "What `verified` rests on here".**

Owner's:

1. Done: every `release/0.N.x` through 0.28.x is pushed, and app #47-#55 are merged (2026-11-17 to 2026-11-21).
2. **Overlap calls:** none open. Urban and NRCS-regulated basins (AL, MO, TN, OK) need `region=` on purpose; say
   if you want a default instead.

   The #340 choices and the method-choice orders are confirmed (2026-11-08). A drop to a drainage-area-only set
   for a missing characteristic now says so in `region_choice_source`, which the app's caption shows.
3. **`allow_extrapolation`** stays allowed everywhere by your decision, including GA, NC and SC above 10 percent
   impervious. The app cautions before proceeding and marks the result (app #37), and the report names the inputs.
4. **Known-miss cautions** cover the three traced cases (OK, MN, LA; `known_misses` in each state's file).
   Add a case by tracing it first, as in `docs/VERIFICATION_PLAN.md`, Checks of 2026-11-05. The 12 of 54
   cross-border gages differing by 2 times or more aren't cautioned, since none is traced to a miss.
5. **MD** stays at the Third Edition, with a caution on every estimate that the Sixth (2023) replaces it for state
   practice (`region_cautions`). A Piedmont boundary basin still needs `region=`.
6. **TX:** the regional sets come before Statewide, and the Colorado sets split at 32 mi2. SIR 2009-5087 (the authors'
   preferred report, not in NSS) is stored as region `Statewide OmegaEM` (2026-11-23), used only with `region=`. Your
   call: whether it should come before Statewide. **It needs a channel slope that StreamStats does not return in
   Texas** (none at the 8 gages saved from the cross-border run; OmegaEM is there), so a StreamStats basin can't run it
   until a slope is supplied. The slope is the main channel's fall over its whole length, not 10-85.
7. **LA's southern gap** uses the 1998 Non-Pine Hills set (TR 60, NSS GC879) when NSS locates no current region,
   with a caution and no interval (`fallback_region`). OH 2006-5312, TN 03-4176 and LA's Delta stay dropped.
8. Issues #293 and #365 are closed. **StreamStats is reachable again** (above): say when the cross-border run
   and the other live checks may go, and whether to tell StreamStats on #365 first.
9. Open the app in a browser: the resizable period grips (app #24), and app #42-#45 (raster hydrograph, basin
   changes, the Bulletin 17C curve, the characteristics form) have never run in one here. The same holds for
   app #47-#52 (merged), which were checked with Streamlit's AppTest only.
10. Have a person spot-check the highest-use states against their reports. Every set was checked only by agents.

Library follow-ups the batches found, none changed:

11. Multi-region weighting: 15 states' reports average in flow and their files set `weighting_space` to
    `linear`: SD, WY, KY, and since 2026-11-15 AZ, CA, CO (also the set awaiting NSS), FL, GA, IA, MT, NC, NJ,
    NV, NY, SC. FL and NJ give the rule in words only. Still open:
    - MT limits the rule to AEPs above 4 percent for a basin running from the Northwest region into the
      Northwest Foothills (SIR 2015-5019-F p. 13). Done 2026-11-16: the library still weights at every AEP, and
      an estimate weighting NW with NWF carries a caution naming the AEPs at 4 percent and rarer
      (`weighting_cautions` in `MT.json`). The channel-width sets (SIR 2020-5142) weren't read for the rule.
    - Read 2026-11-15 and silent on such basins, so left in log: AK, AL (2020 and 2010), AR, CT, DE (2022 and
      2006), HI, ID, KS, LA, MA, ME, MI, MN (2023 and 2019), MO (2014 and 2010), MS, ND, NE, NH, NM (2008), OH,
      OK, OR (SIR 2005-5116), PA, PR, RI, TN, TX, UT (SIR 2007-5158), VA and WV (SIR 2025-5110), VT, WA, WI.
      The reading was a text search of each PDF for weighting language, not a full read.
    - Not read: IL, MD and eastern OR (Cooper, 2006). IN picks a region by point location.
    - SIR 2013-5086 (p. 34) says StreamStats itself weights in flow by drainage area for a site draining more
      than one region, citing Ries and others (2008). If that holds in every state, a silent state's log
      average differs from what StreamStats returns. Whether to make `linear` the default everywhere is the
      maintainer's call.
12. `power` at base 0 (IL stores an offset of 1e-12, which moves a term by at most 1.1e-6 log10 units, so it
    was left). Done: IL's adjusted sets use the report's degrees of freedom (`Variable.fitted`); a rounded AVP now gives way to the SEP in AL
    and TN, IA's CCM above 3.87 goes to the drainage-area-only set, and LA's MAP slope is `CSL10_85` by the
    report's definition.
13. KY regions 1 and 4 and Jefferson have no live slope. TN Memphis has no live CHANCOND. IN prints no
    intervals.
14. **CO and WY have new peak-flow reports coming.** NSS already serves them with placeholder citations: Kohn
    (CO, "SIR 2023-XXXX", GC1949-GC1952) and Taylor and others (WY, "SIR 2027-XXXX", GC2151-GC2154). When
    each is published, verify it and decide whether it replaces the stored set. The live tests ignore those
    regions until then.
15. NH and NY drainage-area-only sets: done (`companion_regions`). The cross-border run estimated 17 New Hampshire
    and 51 New York gages live on 2026-11-17 with the companions in place, all on the full sets; no basin there
    fell outside a full set's ranges, so the drop to drainage-area-only itself has still run only in unit tests.
16. FL, MI and NE can't run through `estimate_at_site` (no StreamStats grid, no NSS region lookup).
    `estimate_from_characteristics` and `flowfreq ungaged-chars` now take the region and characteristics from
    the caller, and the app's ungaged page has a form for it (app #45, on since the v0.26.0 pin).

17. The cross-border check is finished and its offsets traced (`docs/BORDER_QA.md`; `docs/VERIFICATION_PLAN.md`
    has the traces). 1,689 of 1,691 gages attempted; NM 09455800 and 09443000 timed out once at 200 s and were
    left. 15 pairs are systematic. Rhode Island against Connecticut compares unlike equations, and Rhode
    Island's values match its report at 13 of 13 gages. Illinois against Indiana and Colorado against New
    Mexico were the tool's doing (a region taken by list order) and are out. Open:
    - NSS returns every Indiana region, and every New Mexico region, at equal weight for a basin or circle
      near those borders. Why wasn't looked into; it is why `estimate_at_site` needs `region=` there.
    - Rhode Island 01115098 (Peeptoad Brook): StreamStats' stream density is 2.35 against the report's 1.68, so
      today's estimate is 27 percent above the report's printed one. No other gage of the 13 differs by more
      than 3 percent.
    - Run any further StreamStats job on 0.27.1 or later (#389).
18. The two 2026 urban reports are stored as region `Urban 2026` (2026-11-18): SIR 2026-5044 in `AL.json` and
    `MS.json`, SIR 2025-5104 in `TN.json`, from their data releases through `tools/build_urban_2026.py`. Open:
    - NSS served neither on 2026-11-18 (it still serves Alabama's 2010 set and Tennessee's 1984 sets). When it
      does, key the sets to its region codes and decide whether the older sets go.
    - Alabama's and Mississippi's covariance is derived, not printed (the owner accepted this, 2026-11-17).
    - SIR 2025-5104 is stored for Tennessee only. It is the one candidate for urban basins above 1 mi2 in GA,
      NC and SC, where an estimate above 10 percent impervious is refused today; the report doesn't say it
      applies there in place of those states' equations, so using it there is the owner's call.
    - Whether StreamStats computes `L3_PIEDMNT` and `L3_RDGVLY` for Tennessee is unknown, and whether the
      2026-5044 set applies in Alabama's flood region 3 or the Mississippi Alluvial Plain isn't stated.
    - SIR 2025-5034 (Levin, 2025), nonstationary flood frequency by regression in the north-central United
      States, is a reference for `future_flow` and `stationarity`, not an equation set. Not read.

App improvements suggested 2026-11-09 (app #47-#52) were merged 2026-11-17, and the app pins v0.28.0. None
has run in a browser (item 9 above). Not done: the neighbor's estimate one click away from a border caution.

Reports through chat (2026-11-14 and 2026-11-15, unreleased on `main`). `flowfreq.reporting` builds a report from
a JSON specification in Markdown or DOCX (#375), with section kinds `gaged_ffa`, `ungaged`, `weighted` (#380),
`transposed` (#382), `engine_check` (#383) and `text`. `flowfreq-mcp` (#381, the `mcp` extra) gives a chat
client twenty-seven tools over a `ReportSession`; `docs/CHAT_REPORTS.md` has the client configuration. The chat
page in the app was dropped, since only the owner chats. Not yet done:

- Connect `flowfreq-mcp` to Claude Code or Claude Desktop and write a report with it. The tests drive the
  server in process and over stdio from an `mcp.Client`; no chat client has been attached.
- Open a DOCX in Word. The tests read each one back with python-docx, which checks structure, not layout.
- Run the tools that reach a service: `find_gage` and `run_ffa` against the Water Data API with no peaks
  supplied, and `estimate_ungaged` against StreamStats (refused, as above). The tests stub every request.
- The section kinds `low_flow`, `flow_duration` and `future_flow` (2026-11-16): the tests use a synthetic
  daily record and made-up factors. `low_flow` and `flow_duration` were also built once from the Water Data
  API at USGS 03606500 (80 years used, 18 left out as incomplete; 7Q10 37.8 ft3/s), through `build_report`
  and not through a chat client, and nobody has checked those numbers against another program.
- `monthly_statistics`, `summary_hydrograph` and `baseflow` (2026-11-22): built once from the Water Data API at USGS
  03606500 and read in Word. HYSEP is checked against WRIR 96-4040's French Creek example (#412, ported to 0.28.x in
  #413); `lyne_hollick` has not been checked against another program.
  `ih_smoothed_minima` was, on 23 gages against the USGS BFI program's published indexes (OFR 03-146,
  2026-11-23): 0.016 high on average (median 0.012, most 0.061), mostly because the program joins the logarithms of the
  turning points and ours joins the flows. Joining logarithms gives a mean difference of -0.001. Your call: add it as a
  method (or option) named for the program, make it the default, or leave it. The method's documentation no longer
  says it is the Wahl and Wahl program's. Each daily section stores its
  own copy of the record.

Python 3.11–3.14 (`requires-python >= 3.11`; CI matrix in `.github/workflows/ci.yml`).

Tests, measured 2026-10-31 on Windows / CPython 3.12 with the Fortran extension built and
`PYTHONSAFEPATH=1 PYTHONUTF8=1 pytest tests/` (the CI selection; `addopts` deselects
`requires_network`): **4040 passed, 115 deselected, 6 xfailed, 4 failed** (227 s). The four
failures were local-only: two Fortran oracle checks from the documented gfortran toolchain
drift (CLAUDE.md, Validation Status) and two agentops tests from a local profile edit, which #299
removed. `ci / Fortran parity` is the authoritative
parity check.

Open GitHub issues (2026-10-31): #27 (Phase A epic) and #36 (future flow, kept open by
decision). #203-#218 and #246 are closed.

The strict xfails are listed in CLAUDE.md (Validation Status). There are four 2012 PeakfqSA
manual comparisons in `tests/validation/test_big_sandy.py`, which run everywhere and are the
4 xfailed above. There are also two in `tests/fortran_parity/test_fixed_threshold_live.py`,
which need the extension and are where `emafitpr` itself is ill-conditioned (a FIXED threshold
at near-zero skew) or unconverged (`MN2MVARB` on Big Sandy with a 6000 cfs threshold), both
from #81. Run the number, do not carry it forward -- it has been wrong in a commit message
and a PR body already.

Every P1, P2 and P3 item is done; the Done sections and "P3 — The `var_mom` port" hold the
history. Open work is the roadmap: "Carried into 0.11.0" below, then Phase A, Wave 1 and
Wave 2, and the fixed-window ramp metric.
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
- **Legacy NWIS (`nwis.waterdata.usgs.gov`/`waterservices.usgs.gov`) is blocked by the
  egress proxy of a Claude Code *web* session only.** From a CLI session on a developer
  machine it is reachable: the legacy-vs-OGC parity tests for peaks, IV and daily
  (`tests/test_*_backend_parity.py`, `requires_network`) all ran live against it. On
  2026-10-03, from a CLI session, the peak service answered 200, and `waterservices` was
  reached but answered with USGS's own 503 "try again later". In a web session, use the
  committed fixtures under `tests/fixtures/`. The Water Data API (`api.waterdata.usgs.gov`)
  and StreamStats (`streamstats.usgs.gov`) are separate hosts and are reachable from both.
  Check per host, and don't assume "USGS" is one block.
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

Finished work that used to fill the rest of this file is in `docs/HISTORY.md`: the StreamStats, transposition and Fortran-engine "Done" blocks, the P1 and P2 cleanup, and the original Phase 0 to 16 build checklist (which describes a `flowfreq/peakfqsa/` subsystem that was later removed). The `var_mom` port (P3) and the peakfqr reference notes stay here, because code comments and `docs/FORTRAN_UPLOAD.md` point at them.

### Nationwide rollout (2026-10-25 plan)

Plan: `docs/superpowers/plans/2026-10-25-nationwide-rollout.md`; design:
`docs/superpowers/specs/2026-10-25-nationwide-rollout-design.md` ("NSS sweep" and "A. The
`nss_snapshot` tier").

- [x] **A. `nss_snapshot` tier and the 41-jurisdiction sweep.** Three NSS parser extensions
      (mixed-case codes, `logN` in the exponent of 10, negated affine bases), Indiana gating,
      `snapshot_to_library`, `tools/build_nss_snapshot_library.py`, opt-in
      `load_state(..., allow_nss_snapshot=True)`, and 41 jurisdictions with 1817 equations
      (1441 with `sep_log` at harvest, 1473 since #213). A snapshot is not "supported": every jurisdiction still owes the
      roadmap 3.2 definition of done.
- [x] **24 region-indicator equations** (#211) in GA, NC and SC (8 each, SIR 2023-5006) parse
      since `Variable.times` (an interaction term `T(x) * Y`) and match NSS's Estimate values to
      its 3-figure rounding. SC #208 is the same study (GC1935).
- [x] **New York's 2 urban regions** (#212) (GC1077, GC1078, WRI84_4350 Method 1 and 2 Urban
      Peak): out of scope. They are Sauer and others' (1983) nationwide urban adjustments
      (WRI 84-4350 eqs. 2A-2F and 3A-3F), each a power of the *rural* peak `(Q_T)r` from
      another equation set, so NSS's Estimate cannot evaluate them on their own. The report
      evaluates them for New York without recommending them (bias +28 to +53 percent), and
      they cover only the 2- to 100-year floods. They stay in `StateLibrary.unparsed`;
      supporting them would need an equation form that consumes another region's estimate.
- [x] **ASEp below 1** (#213) on 32 LA, AR and MS equations: NSS stores the report's SEP percent
      divided by 100 (SIR 2021-5046 Table 2; SIR 2024-5031 Table 3). `ASEP_AS_FRACTION` in
      `nss_library` names the four regions, and the rebuilt libraries carry `sep_log`
      (1473 of the snapshot equations now have one).
- [x] **The advisory warning** (#214, PR #230): `load_state(..., warn=False)` logs an opt-in
      set's advisory at DEBUG, and the coverage tool uses it.

- [x] **B. `estimate_at_site` and `flowfreq ungaged`** (0.16.0, #197), with
      `OverlappingRegionsError` and `region=`. The flowfreq-app ungaged page is app PR #7
      (pinned to v0.16.0); the app is not merged from here.
- [x] **C. Regional skew rows for all 56 jurisdictions** (0.16.0, #195/#196): 77 rows,
      34 verified, with the state-level fallback.
- [x] **E. Nationwide smoke run and docs** (#200, #201, #202). Docs done (README Coverage
      and Ungaged sites, MASTER_ROADMAP §0/§4, this file). The smoke run
      (`tools/nationwide_smoke.py`, report `docs/NATIONWIDE_SMOKE.md`) first ran on 2026-10-25
      (55 of 56 gaged, 30 of 56 ungaged). The 2026-10-30 rerun, after the follow-ups below,
      gets 56 of 56 gaged and 42 of 56 ungaged, with `tools/nationwide_smoke_points.csv`
      supplying a region (`region` column) or a pour point for 18 states.

Follow-ups from the rollout reviews (issues opened 2026-10-26):

- [x] **Multi-code `region=`** (#203, PR #239). IA and MT use it in the rerun.
- [x] **HUC membership rows for the new skew studies** (#204, PRs #251 and this one). Rows for
      every study the reports define by hydrologic unit: Great Lakes and Ohio, Eastern NY and
      PA, SIR 2025-5110 Skew Regions 1 and 2, New England, Arkansas-Louisiana, and Alaska's
      Skew Areas 1 and 2 (SIR 2016-5024 Table 5; the rest of HU 19 is the exclusion area, and
      Table 6's drainage-area ranges apply through `drainage_area_sqmi=`). A HUC may carry one
      row per overlapping study, chosen by the site's state. Statewide studies need no rows.
- [x] **A nationwide HUC8-states table** (#205): `tools/build_huc8_states.py` pages the WBD HU8
      layer; `data/regional_skew_huc8_states.csv` lists every HUC8 (2,456 rows on 2026-10-30, the
      417 old rows kept verbatim), so a border HU12 resolves its state without `state=`.
- [x] **Catalog defect: gage 15388960 filed under DC** (#206, PRs #234 and #243). The catalog's DC
      rows are rebuilt, and DC's smoke gage is now 01651770.
- [x] **NSS region mismatches** (#207, PR #240; #211). `region=` names a library region where NSS
      locates no peak-flow region (WY, NV), and NC's GC1936 parses since #211.
- [x] **SC GC1587** (#208): not a range error; the gage belongs to GC1935 (#211). SC passes.
- [x] **Region choice where sets overlap** (#209): decided 2026-10-30. Keep the explicit
      `region=`, with no per-state default. The rerun's picks, with the reason for each, are in
      `tools/nationwide_smoke_points.csv`.
- [x] **Unsnappable smoke points** (#210, PR #232): `--points` with hand-checked pour points.
      IL, NH and PR snap a few tens of metres from the gage. At ME, AL, VT and HI the gage
      coordinates snapped to a tributary (0.03 to 1.9 mi2), so their points were checked
      against the catalog drainage area.

Still failing in the 2026-10-30 rerun (14 of 56; WY and MN fixed since, live recheck pending):

- [ ] **No library**: DC, VI, GU, MP and AS. NSS serves no peak-flow equations for these; expected.
- [x] **MN** (#246): the gage is in the Rainy River Basin, which SIR 2023-5079 leaves to SIR
      2019-5012 (GC1775, NSS region `RRB`). GC1775 is now in MN's library and located through
      `NSS_EXTRA_REGIONS`; the live smoke rerun waits on StreamStats (ss-hydro 502).
- [ ] **OR**: GC729 to GC731 need characteristics StreamStats did not return at the gage
      (I24H2Y, JANMAXT2K, WATCAPORC, SOILPERM, BSLOPD, JANMINT2K).
- [x] **WY**: region 1 needs `LONG_OUT` (outlet longitude), which StreamStats does not return;
      `estimate_at_site` now takes it from the snapped pour point (`outlet_coordinates`). The
      live smoke rerun waits on StreamStats (ss-hydro returned 502 on 2026-10-30).
- [x] **MT**: refused as out of range, not extrapolated (SW CONTDA 2616 > max 2472.17).
      Pinned in `tests/test_smoke_failures.py`; the smoke point moved to gage 06035000
      (Willow Creek near Harrison, 85.9 mi2, SW only).
- [x] **FL, MI, NE**: StreamStats' snap service has no stream grid for these regions
      (404 "Stream grid raster not found" / "No such region", 2026-10-31; ss-delineate
      fails the same way). Now `RegionNotDelineatedError`; a known gap in the smoke report,
      since no point there can be delineated through the API.
- [x] **ND**: StreamStats flags the Cannonball basin REGULATED; now
      `RegulatedDelineationError`. The smoke point moved to gage 06352000 (Cedar Creek near
      Haynes, 557 mi2, reference).
- [x] **AK**: gage 15072000 does not snap (nor 10 of the next 11 AK candidates). The
      smoke point moved to gage 15290000 (Little Susitna River near Palmer, 62 mi2) in
      GC1656, which `assemble` now picks by drainage area from AK's two overlapping statewide
      classes (under and over 1,000 mi2).
- [x] **Live recheck** of the moved MT, ND and AK points, with WY and MN: all estimate
      (2026-11-03, `docs/NATIONWIDE_SMOKE.md`).
- Large gaged-versus-ungaged differences with a correct delineation, checked 2026-10-31:
  - **AZ 09484000** (ratio 0.10) was the wrong region. The smoke point forced GC1618 (High
    Elevation), which NSS uses only at a mean basin elevation of 7,500 ft or more; Sabino
    Creek's is about 5,980 ft (3DEP, 59 points), so region 5 (GC1622, Southeastern Basin and
    Range, like neighbouring Rincon Creek) applies: 8,604 cfs at 1% against 14,403 gaged
    (0.60, inside the 90% interval). `assemble` now chooses it from `ELEV` through AZ's
    region gates (`StateLibrary.region_gates`), so AZ has no row in the points file.
  - **NV 10329500** (0.36): Northwest Region 2 is the report's own region for the gage, the
    equation is double-entered from WSP 2433 Table 6, and 0.36 is about 1.65 standard errors
    of prediction. No defect found.
  - **TX 08190000** (0.20): the NSS snapshot's GC2146 (SIR 2008-5084, drainage area and
    precipitation only) is stored exactly as NSS serves it and reproduces the smoke value. The
    Nueces at Laguna is on the Balcones Escarpment (616,000 cfs in 1935), where a statewide
    equation with no regional term underestimates; a limitation of the published set.
- [x] **Stationarity by default**: `run_ffa` attaches the screen dictionary; the figure is
      opt-in (`stationarity_figure=True`).
- [ ] **Regional skew `pending` rows** (41): each waits for a published B-WLS/B-GLS study or
      a reading of the report's printed value; the row's notes say which. Re-checked
      2026-11-15: a Publications Warehouse search for 2025 and 2026 found no regional skew
      study beyond SIR 2025-5110 (VA and WV), which is already used.

The methow-flow items and the `pending` skew rows have no issues; #203-#218 are closed.

The parked minors found by the rollout reviews and fixed after 0.16.0 (see CHANGELOG,
Unreleased): the `aeps=[]` guard, `region=` case-insensitivity and its ambiguity refusal, a
shared `variance_log` refusal, a public `PERCENT_WEIGHT_TOL`, the Tennessee Skew Region 2
notes trimmed, and `regional_skew_at_huc`'s docstring and membership predicate.

### Agent ops (2026-10-26)

- [x] **`tools/agentops/`**, ported from flowfield2d: runner, governor profiles, usage meter,
      hooks (`docs/AGENTOPS.md`). The hooks and status line still need merging into
      `.claude/settings.json` from `tools/agentops/settings.example.json` (owner action).
- [x] **First and second card batches** (2026-10-26 and 2026-10-30), all merged on Sonnet after
      the agy fixes (host allow rules, per-call `--settings` #233, `AGENTS.md` seeded into
      worktrees #229): i214 #230, i210 #232, i206 #234, i205 #235, m12 #236, m4 #237, m7 #238,
      i203 #239, i207 #240. The decisions they waited on are made and done: #211 (#244),
      #209 (explicit `region=`), #213 (#245), #204 (#251, #252), #218 (#242); the DC catalog
      was rebuilt live (#243). The card queue is empty.
- [x] **Methow item 7, part 1 merged** (`m7-monthly-transfer`, PR #238). Owner decisions,
      2026-10-26: transfer by the runoff-volume ratio `(At Rt)/(Ag Rg)` with exponent 1
      (as production, METHOW_METHODS §12); sampling variance and covariance from concurrent
      water years (delta method for the mean, bivariate-normal asymptotics for the median);
      zero-flow gage-months excluded and reported. Statistics: the mean and the median of
      the yearly monthly means (the hydraulic handoff uses May and September medians).
      Part 2: monthly leave-one-out and blend; part 3: duration quantiles (q10, q50, q90).

### Priorities for methow-flow (set 2026-10-21)

Roadmap items ordered by what they do for the methow-flow analysis (gage transfer in Peak
Region 2, `flowfreq.transfer_variance` since v0.13.0). Roadmap section in brackets.

Tier 1 -- unblocks the transfer-variance work:
1. [x] **EMA variance of each quantile estimate, then B17C Appendix 9 gage weighting**
       (§2). `FrequencyResults` does not carry it, so methow-flow's
       `reach_flows.weight_with_regression` back-solves a variance from the 90 % interval
       width, and `transfer_variance` takes `V_g`/`V_i` from WREG's formula, which ignores
       censoring and historic periods. Check against `emafitpr`'s `var_est`, then move the
       weighting upstream.
2. [x] **Nested-gage finder** (§1.2 NLDI, §5.2 same-stream check). The missing input to
       `build_pairs`. Upstream methow-flow's `flow_relation` (polygon containment) and add
       NLDI navigation for gages without a cached polygon.
3. [x] **Regression-weighted drainage-area ratio** (§5.1, SIR 2016-5118 eqs. 11-12, Ries 2007).
       Method 2 in the transfer-variance validation.

Tier 2 -- fixes known Methow biases or moves general code upstream:
4. [x] **Interpolation between two gages** (§5.1; `transpose.interpolate_between_gages`, #237): upstream methow-flow's bracketing
       (`choose_bracketing_gages`, `bracketed_peak_flows`) and reconcile it with
       `transfer_variance.combine`.
5. [x] **GLS validation against SIR 2016-5118** (§3.3, Northwest task 4; done in
       `docs/WA_REGION4_REFIT.md`), and why Peak Region 2 is low at the Methow
       (`docs/WA_REGION2_METHOW.md`, `tools/wa_region2_elevation.py`): mean basin elevation,
       absent from Table 6's design, explains about half the bias at the 50% AEP (t = 4.1,
       Methow ratio 2.08 to 1.46; regional, it holds without the Methow HUC8's gages) and
       fades by the 1% AEP. No Methow-specific effect remains after it. For methow-flow: an
       ELEV refit for the 50%/20% AEPs, the published equation for rarer floods.
6. [x] **Regulation guardrail in `transpose` and `qppq`** (§5.2), with a caller-supplied
       flag for diversions GAGES-II does not see (Beaver Creek):
       `donor_regulation`/`target_regulation` on every transfer, and
       `classify_site(known_diversions=...)`.

Tier 3 -- ordered for the per-reach monthly mean flows (WY1992-2025) that every EDT flow
attribute in methow-flow waits on (`docs/EDT_ATTRIBUTE_CROSSWALK.md` §6 there):
7. [x] **Transfer variance for monthly and duration flows** (§5.1); decisions made
       2026-10-26. Part 1 `flowfreq.monthly_transfer` (#238); part 2 monthly leave-one-out
       and blend (#260); part 3 `flowfreq.duration_transfer`, median annual q10/q50/q90
       (2026-10-31). Open from the original list: the exponent per month stays 1 (no
       duration-indexed exponent set is published for the Methow). Extend
       `transfer_variance` from flood quantiles to a monthly mean per calendar month and
       to duration quantiles: pairs and `s2_d` fitted per month, the GLS blend and
       leave-one-out unchanged. Open questions: the exponent `b` per month (no
       duration-indexed exponent set is published for the Methow), sampling covariance of
       monthly means across nested gages, and months with zero flow (the Mazama gage's
       October-March dry days), where log space fails. Needs item 6 first, since
       diversions move the low-flow months most.
8. [x] **Stationarity screening** (§6.1; `flowfreq.stationarity`, 2026-10-31: Mann-Kendall
       with Hamed-Rao, Sen's slope, Spearman, Pettitt, CUSUM, `field_significance` (FDR, Walker),
       `attribute_change_point` and `stationarity_screen`; then `lombard`,
       `field_significance_bootstrap` and `basin_change_events` (NID/NLCD); then `pot_counts`,
       `pot_count_trend` and `stationarity_report`). The EDT attributes compare two periods (§15
       Q2, §16 low flow, §18 TQmean), and a monthly transfer over WY1992-2025 assumes the
       donor relation holds across it.
9. [x] **Mixed snowmelt/rain-on-snow populations** (§2, warn only):
       `stationarity.mixed_population_screen` (cool-season share of annual peaks and the
       concentration of their dates; warns inside a 20-80 percent split).
10. [x] **Scheduled live-network CI** (§7): the weekly `Live services` workflow
       (`.github/workflows/`), which opens an issue on failure (#218 was one).
11. [x] **Naturalization add-back** (§5.2; `flowfreq.naturalize`, 2026-10-31: `Diversion`,
       `naturalize`; no routing lag or storage). Gaged flow plus the diversions upstream of a
       site, net of return flow (canal seepage), as a daily or monthly series a transfer can
       start from. `qppq` says it does not naturalize, so a regulated donor carries its
       regulation into every transfer. Inputs are caller-supplied (diversion rates,
       return fractions, periods in effect); estimating them is out of scope. In
       methow-flow it is the §16 low-flow template (the 2002 canal inventory).
12. [x] **Uncertainty in choosing the dissimilarity power** (#181 follow-up; #236:
       `bootstrap_band(profile=...)` re-chooses the power in every replicate).
       `bootstrap_band` holds a `powers=` choice fixed, so the band interval leaves out
       the uncertainty in choosing it. Re-run `profile_power` inside each replicate and
       report how often each grid power wins alongside the interval.

Raised by methow-flow's Peak Region 2 run (methow-flow #29, merged 2026-10-22):
- [x] **Non-squared dissimilarity terms** (#181, PR #184): `powers=` on
      `fit_dissimilarity`, `bootstrap_band` and `leave_one_out`, `solve_band` for any power,
      and `profile_power` to choose one by likelihood. Next, in methow-flow: profile the
      `log_ar` power on the Region 2 pairs and, if one power fits the full 1-3,000x range,
      drop the two-model stopgap (within-4x band fit plus full-range fit) and its own fold
      loop. Still open here: the bootstrap holds the power fixed, so the band interval
      leaves out the uncertainty in choosing it (item 12).
- [x] **Indefinite blend covariance from a negative donor rho** (#177). `leave_one_out`
      floors fold estimates at 0 (`min_donor_correlation`), reports the raw estimate and a
      clipped flag, and records a gage it cannot score with a `note` instead of aborting.
      `combine` names the `-1/(k-1)` bound. methow-flow can drop its own floored full-data
      `rho` workaround.

### Carried into 0.11.0

Open at v0.10.4 and carried into the 0.11.0 (countrywide) cycle. Details are in the
Phase A / Wave entries below.

- [x] **Partial states.** None left: NV is `verified` since #157.
      - ~~WA region 4 covariance (#65)~~: stored in the `P/10` basis after a GLS refit
        (`docs/WA_REGION4_REFIT.md`); WA is `verified`.
      - ~~UT region 2 (#96)~~: NSS's 4-digit precipitation bases are stored; each rounds to
        SIR 2007-5158's printed 3-digit base. UT is `verified` (#144, merged 2026-10-13).
      - ~~NV region 6 and the hybrid regions (#104)~~: region 6's 2-year `Q=0` is a
        `form="constant"` equation, and the 11 hybrid errors are `published_errors` that never
        feed an interval (#157, `docs/NEVADA_REGION6_SCHEMA.md`).
        Its equivalent years of record are stored too, with the report's gaged-site
        weighting (eq. 1), and `tools/snapshot_nss.py` parses NSS's `0*(...)^(0)` as the
        constant form.
- [x] **Oregon closed basins (1712)**: provisional PNW member after a B-GLS applicability
      test (`docs/OREGON_UNRESOLVED_SKEW.md`). The Klamath (1801) stays unresolved (2 gages).
- [x] **Montana skew v2** (provisional, opt-in): east-of-Divide G = -0.04, MSE 0.385
      (`docs/MONTANA_REGIONAL_SKEW_V2.md`). PRECIP is the significant covariate but is not
      representable in the table. A published USGS study is still wanted.
- [x] **Snake River Plain skew**: not feasible, 2 of 47 gages (`docs/SNAKE_RIVER_PLAIN_SKEW.md`).
- [ ] **Regional skew** for MT (a published USGS Bayesian study; the provisional #91 study
      stays out of the table) and for the rest of Wave 2 (roadmap §1.3). Arizona is done
      (SIR 2014-5211 app. 1, G = -0.09, MSE 0.08, 84 HUC8s valid in AZ only). CO, UT, WY
      and NV were checked 2026-10-12 and have no B-WLS/B-GLS study; each has a `pending` row
      saying what USGS uses instead (CO and WY: the B17B Plate I map; UT: the Perica and
      Stayner 2004 contour map; NV: WSP 2433's zero, MSE 0.31). NM waits on Tillery 2026.
- [x] **Colorado: SIR 2025-5047** (#151, #156): stored in `awaiting_nss/CO.json`, a new
      opt-in status for USGS-published sets NSS does not serve yet
      (`load_state("CO", allow_awaiting_nss=True)`). It uses NSS staging citation 200's codes,
      region bounds and limits, with Table 1's coefficients. Once NSS publishes citation 200,
      it replaces the Mountain, Northwest, Rio Grande and Southwest regions in `CO.json`.
      SWE is in inches, as stored, although Table 1.2 lists SNODAS SWE in meters
      (`docs/COLORADO_SIR_2025_5047.md`, checked 2026-10-18).
      Original entry: (Kohn, Mast and Gross, 2026) replaces Capesius and
      Stephens (2009) for central and western Colorado: four elevation regions (Plateau,
      Mid-Elevation, Sub-Alpine, Alpine), 418 gages through WY2019, tables as CSV on the pubs
      page. NSS still serves the 2009 equations (citations 11, 20, 21, checked 2026-10-12),
      so `CO.json` matches StreamStats today. Transcribe when NSS adopts it, or earlier as a
      second CO set if wanted.
- [x] **Northwest studies plan** (`docs/superpowers/plans/2026-10-11-northwest-regional-studies.md`), done 2026-10-12
      (WA merged as #132 after the maintainer's go-ahead). Follow-on:
      `docs/superpowers/plans/2026-10-12-oregon-regression-update.md`. Oregon StreamStats
      returned only 9-13 basin characteristics on 2026-10-12 (a missing-raster error), so it is
      not evaluating Cooper's western-Oregon variables (I24H2Y, JANMAXT2K). Re-check it, and
      report it to StreamStats. Oregon plan Tasks 2-5 done 2026-10-13: a provisional
      six-region set in `provisional/OR.json`, less precise than the current equations
      (`docs/OREGON_REGRESSION_UPDATE.md`). Soil Ksat (STATSGO2, #153) and young-volcanic
      geology and lake area (SGMC/NHD, #158) followed. Only region E's stored equations changed
      (SEP at 1 % AEP 94 -> 86 %), and every region is still less precise than the current
      equations. #164 then added one variable set chosen across the 50/10/1/0.2 % AEPs
      (`--select multi`, with a leverage limit of 0.5), and a young-volcanic fraction from
      Oregon's state geologic map (OGDC-8, `VOLCQP_OGDC`). W-171003's SEP is now
      47/42/46 % (it was 47/49/54) and E's 85/85/99 % (it was 90/87/101). Lake area for
      Waldo Lake is still not significant (t = -1.9), and Waldo's own leverage is 0.55-0.96.
      #166 (merged 2026-10-21) takes the three follow-ups:
      - **High Cascades base flow.** Wolock's (2003) BFI grid (`BFI_WOLOCK2003`) was tested
        and not adopted. It is selected in no region, and the High Cascades 1 % residuals
        stay at -0.30 to -0.39.
      - **Lake-dominated basins.** A `LAKEAREA_NHD <= 5 %` applicability limit, through a new
        `RegressionEquation.applicability` field, adopted with #166's merge. The five gages
        above 5 % average -0.38 at 1 %.
      - **W-170800 leverage.** W-170800 is refitted under the leverage limit as DRNAREA,
        BSLOPD. Its 1 % SEP rises from 51.2 to 53.7 %.

      Still open after #166:
      - a High Cascades groundwater variable other than BFI;
      - whether to refit without the five gages above the lake limit (1 % SEP -6 points in
        W-171003, -3 in W-170900).

      #135 (06163400 vs PeakFQ) fixed: lowercase peak codes. A 2026-10-11 search
      found no newer USGS report for WA, MT, ID or OR. The studies follow published USGS
      procedures, and their results enter `regional_skew.csv` only as `provisional`, which
      callers must opt into (`allow_provisional=True`).
- [ ] **New Mexico: newer USGS report** (Wave 2). Tillery and others (2026) supersedes
      Waltemeyer (2008), which `NM.json` uses. It has nine flood regions, data through
      WY2020, and a regional skew of -0.145 (SE 0.454). Data release:
      https://doi.org/10.5066/P9QWEXCY. Blocked as of 2026-10-12: the report is not in the
      USGS Publications Warehouse, NM StreamStats still cites Waltemeyer (2008), and the data
      release was withdrawn for revision on 2026-09-29 (its files are no longer public). The
      NM skew row is `pending` with both values in its notes (Waltemeyer 2008 p. 5: zero, MSE
      0.31, B17B-era). Re-check the Publications Warehouse and the data release; transcribe
      the skew, its region and the nine regions' equations from the published report.
      Re-checked 2026-10-31: still not in the Publications Warehouse.
      Re-checked 2026-11-15: still not in the Warehouse. The data release is titled "(Under
      Revision)", is marked temporarily unavailable and lists no files (ScienceBase, last
      updated 2026-11-09).
- [ ] **Future flow (#36), kept open** (decided 2026-10-11). No state tabulates per-AEP
      factors. The issue stays open for WDFW-style per-site tools rather than closing on
      the survey.

### Roadmap Phase A, Wave 1 and Wave 2 (issues #27–#40, #93–#104)

`docs/MASTER_ROADMAP.md` requires a TODO.md entry for every roadmap item that has
started. PR #26 scaffolded Phase A and Wave 1. Since then almost every item has been
delivered; this section keeps one line of history per item and the open remainder.
Audited 2026-10-03 against `git log v0.9.0..origin/main`, `CHANGELOG.md` and
`gh issue list --state all`.

**Phase A: data foundation (#27, epic, open while #32 and #33 are open)**

- [x] **#29 Water Data OGC API migration** (closed). The OGC API is the default backend
      for peaks (#42, parity #49, default #63, `USGSgage.download_peak_flow` routed
      through `peak_sources` #70), instantaneous values (#45, default #77) and daily
      values (#79), each switched only after a live parity test against legacy NWIS. #84
      added `waterdata.request`, with an `X-Api-Key` from `USGS_API_KEY` and 429/503
      backoff. `backend="nwis-legacy"` keeps the old services. The design questions are
      settled: `peak_date` is the API's UTC date (0.9.0 breaking-change note), `ts_id` is
      the 32-hex `time_series_id` on the OGC backend, and a site with several IV series
      and no `ts_id` raises rather than merging them. The legacy partial-date fix (#41) keeps peaks with
      an unknown day or month. The API carries no perception thresholds, so they remain
      analyst inputs.
- [x] **#30 Peak qualification codes** (closed). `peak_codes` follows peakfq's `siteQT`.
      Applied by default in `run_ffa`, `compare_engines`, `analyze_gage` and
      `flowfreq compare` (#52, #68). Codes 4/8 are fitted as intervals (#72), and code 7
      with 4/8 as historic intervals (#90).
- [x] **#31 PeakFQ `.psf` reader and converter** (closed). `psf.py`, `psf_convert.py`
      (#52) and `watstore.py`. All 24 WY/MT stations in `wymt_ffa_2022A.psf` convert and
      match live `emafitpr` (CLAUDE.md, Validation Status).
- [x] **#32 Regulation / urbanization screen** (closed 2026-10-06; v0.10.2 added the opt-in
      NID/NLCD current-data screen, see CHANGELOG.md). Done in #95: `flowfreq.regulation`
      combines the packaged GAGES-II screen (9,322 gages) with peak code 6, and
      `analyze_gage` (and `run_ffa` given `site_no=`) refuses a regulated gage unless
      `allow_regulated=True`. Open:
      - Current NID storage. The national CSV is live
        (`nid.sec.usace.army.mil/api/nation/gpkg`, #155), but summing it upstream needs basin
        polygons and point-in-polygon, which flowfreq does not depend on.
      - The NLCD impervious **time series** (step 3). GAGES-II carries NLCD 2006 only.
      - Gages outside GAGES-II, which have peak-code evidence only.
- [x] **#33 National gage catalog** (closed 2026-10-24). Built in #105:
      `flowfreq/data/gage_catalog.csv.gz` covers 20,745 USGS peak sites with at least 10
      water years, nationally, rebuilt by `tools/build_gage_catalog.py` with two Water
      Data API bulk queries per state. `regression_region` is filled for all of Waves 1
      and 2. WA/OR/ID/MT and CO/UT/AZ use point-in-polygon. WY, NV and NM have no region
      polygons, so they use the region each report prints for its gages
      (`tools/extract_station_regions.py`, `data/catalog_regions/`): 251 of 504 WY,
      151 of 342 NV and 193 of 405 NM gages (2026-10-18). The rest were not in the
      reports and stay blank. A point-in-polygon region comes from the gage point,
      so a basin that crosses a boundary needs
      `streamstats.locate_regression_regions`. Western Oregon's 2A/2B split is by basin
      mean elevation, so those rows list `GC729;GC730;GC731`.
- [x] **#34 Regional skew table and lookup** (closed). Pacific Northwest rows for
      WA/OR/ID (#51, #69); lookup by HUC or location, `regional_skew_at` /
      `regional_skew_at_huc` / `regional_skew_for_site` (#83); no silent `-0.302`
      default (#44); B-WLS/B-GLS tooling (`skew_study`, #91). Still open from it: the
      Montana skew (below).

**Wave 1: Columbia River basin (#28, epic, closed)**

- [x] **#35 Regression schema and offline evaluator** (closed). `tools/snapshot_nss.py`
      and the exact NSS equation parser (#47), `Variable.scale`/`offset` (#58, #60).
      NSS's `sep` is kept verbatim as `nss_sep` and never mapped to `sep_log`, because it
      does not behave like a log10 SEP.
- [ ] **#36 Future-flow `ChangeFactorSet`** (open). Framework, `data/future/*.factors.json`
      loader and the national source review are done (#92). HEC-17 (2016) and NCHRP 15-61
      (2019), checked against the primary documents, give methods and single-site worked
      examples but **no tabulated national factors**, so no national set ships
      (`docs/FUTURE_FLOW_GUIDANCE.md`). The Wave 1 and Wave 2 state survey found no state
      that tabulates per-AEP factors either (WA requires one, read per site from WDFW's
      tool), so no state set ships. NCHRP 15-61's Ch. 6/7/8 procedures are implemented and
      reproduce its Tables 6.8-6.9, 7.1 and 8.1; `apply_with_precedence` reports the
      national result alongside a state override. Open: NOAA Atlas 15 Volume 2 once
      published (plus a cited precipitation-to-flow elasticity), and the unretrieved NMDOT
      (2018) and NDOT (2006) drainage manuals.
- [x] **#37 WA (pilot)**, `verified` since the Region 4 refit (Northwest studies task 4;
      `docs/WA_REGION4_REFIT.md`). A GLS refit of the report's own data and Table 6's Sp both
      put Table 7's Region 4 matrix in the `P/10` basis, so it is stored with PRECPRIS10
      `scale = 0.1`. The report's workbook and Table 8 apply it to raw P (#65); flowfreq's
      Region 4 intervals are about 3x narrower than the workbook's.
- [x] **#38 OR**, `verified` since #107: western Oregon SIR 2005-5116 (#57), the 2A/2B
      blend (#67), eastern Oregon OWRD OFR SW 06-001 (#71). No covariance is published.
- [x] **#39 ID**, `verified` (#48, #59).
- [x] **#40 MT**, `verified`: SIR 2015-5019-F (#56) and the SIR 2020-5142 channel-width
      equations with eqs. 5-12 weighting (#73).
- [ ] **Montana regional skew.** `regional_skew.csv` keeps MT `pending`. No published
      B-WLS/B-GLS study covers the state, and USGS Montana does not adopt the PNW study in
      HUC 1701. #91 built a **provisional, non-USGS** statewide study
      (`docs/MONTANA_REGIONAL_SKEW_PROVISIONAL.md`), which is deliberately not in the
      table. This stays open for a published USGS study (roadmap §1.3).

**Wave 2: Colorado River basin (#93–#104)**

Status per the definition in `flowfreq/regression/library.py` (#107): `verified` means
everything the cited reports publish for the stored regions is stored, double-entered and
cross-checked, against live NSS where NSS serves the region. An unpublished covariance
does not block it. `partial` means something the reports *do* publish is missing or
unresolved.

- [x] **#93 CO**, `verified` (#94).
- [x] **#96 UT**, `verified` (#97, region 2 bases 2026-10-12). Table 5 prints region 2's
      precipitation bases to 3 significant figures, NSS uses 4 (1.07 against 1.074). No
      published source gives the fourth digit, but every NSS base rounds to the printed
      one and the base the report leaves unrounded (PK200) is identical, so NSS's are
      stored, with both values in each equation's citation. Region 2 now matches live NSS.
- [x] **#98 WY**, `verified` (#99).
- [x] **#100 NM**, `verified` (#101).
- [x] **#102 AZ**, `verified` since v0.10.2: region 5's `10^(a - b DRNAREA^-c)` form is
      stored with `Variable(transform="power")`.
- [x] **#104 NV**, `verified` since #157. Before that, `partial` (#106). v0.10.2 added the transition-zone equations
      (`flowfreq.regression.nevada`) and replaced NSS's limits with the ranges of each
      region's gages. Still open: region 6's 2-year equation is printed `Q=0`, and
      hybrid-method errors are not comparable standard errors; neither is stored.

Arizona has a verified regional skew (SIR 2014-5211 app. 1). CO, UT, WY and NV have no B-GLS study (pending rows); NM waits on its 2026 report. The wave's cross-border consistency check
(roadmap §3.2) ran 2026-10-06 (`docs/WAVE_QA_W1_W2.md`). Their §6.3.2 future-flow guidance is recorded: none tabulates factors
(`docs/FUTURE_FLOW_GUIDANCE.md`).

### Small open items

- [x] **Unmapped Water Data peak qualifier tokens.** Done in #149. `peak_sources.QUALIFIER_TOKEN_CODES`
      warns on `YEARUNKNOWN` (NWIS code A, year not exact) and on gage-height flags
      (`BACKWATER`, `NOTMAXGH`, `EQUIP`, `ICE`, `BLWMIN`), plus `OTHERAGENCY`. None changes
      a B17C treatment (peakfq ignores A), but a bulk download logs one warning per row:
      hundreds in the Oregon dataset build. Map A like `DAYUNKNOWN` and the rest like
      `GHNOTASSCPKQ`.
- [x] **Live-service tests (#140).** Re-run 2026-10-13: 106 passed, and StreamStats computes
      PRECIP in Washington again. NID's national CSV (`api/nation/csv`) returns 403 from NID's
      own proxy, and its S3 object is access-denied, so `load_nid` now reads the national
      GeoPackage (`api/nation/gpkg`) with stdlib `sqlite3` (#155). Its `dams` table holds the
      CSV's values as the same text. All four live regulation tests pass, including the end-to-end
      classification, which had been skipping.
- [x] **v0.10.4 GitHub release.** Deleted 2026-10-13 to re-attach the dists; the workflow
      re-run could not recreate it (PyPI refuses the duplicate upload; `skip-existing` landed
      in 0.11.0, after the tag). It exists again: published 2026-10-14 with both dists
      attached (`gh release view v0.10.4`, checked 2026-11-15).

Earlier items, all closed; each verified in the code on 2026-10-03.

- [x] **#13 Screen donors on basin similarity, not drainage area alone.**
      `flowfreq.donor_similarity` (#54): `screen_donors` (area plus per-predictor bands)
      and `rank_donors_by_similarity` (Burn 1990 region-of-influence distance). Opt-in;
      `transpose`'s area screen is unchanged.
- [x] **`batch.run_multi_site` cannot analyze real `fetch_nwis_batch` output.** Fixed in
      #46; the strict xfail in `tests/test_batch.py` now passes as a plain test.
- [x] **`make clean` misses the Windows extension.** Fixed in #46: `clean` also removes
      `_emafort*.pyd` and the MinGW DLLs.
- [x] **CLI `compare` has no historical-peak or perception-threshold flags.**
      `--historical` and `--threshold` (#50); peak codes are applied by default (#68).
- [x] **`as_G_PRL_o` (pseudo effective record length) is not surfaced.**
      `FrequencyResults.pseudo_record_length` on the native EMA path (#46), matching the
      golden files (54.378 vs 54.373 on Big Sandy, 98.0 exactly on 12363000).

### Sub-daily metrics follow-ups

`flowfreq/subdaily.py` landed in 0.9.0 (see `docs/SUBDAILY_METRICS_DESIGN.md`). Of the
four follow-ups below, three are done (#85, #87, #89); the fixed-window ramp metric
is open.

- [x] **`diel_variation`'s `expected_obs` is a fixed 1440 minutes / median step,
      so it mismarks the two daylight-saving transition days of every year.**
      **Done for v0.10.0:** it now uses each local day's actual length, like
      `daily_extreme_timing` and `ramping_rates`. Only `expected_obs`/`complete`
      on the two transition days moved; CHANGELOG Fixed says exactly which.
- [x] **Paired stage-and-discharge retrieval.** `download_instantaneous_flow`
      and `download_instantaneous_stage` are separate calls returning separate
      frames. Anyone comparing a cfs/hr limit against a ft/hr one joins them on
      the index by hand. A single call retrieving both parameters onto one index
      is the natural next step and is not built.
      **Done:** `USGSgage.download_instantaneous_flow_and_stage` and
      `join_flow_and_stage`: an outer join on UTC with NaN, never interpolation,
      per-parameter qualification codes, and mismatched local time refused.
- [x] **No validation against a published ramping-rate or peak-timing figure for
      a real gage.** **Done (2026-09-30), for daily extremes:** Exelon (2012), *Final Study Report: Downstream Flow Ramping and Stranding Study, RSP 3.8*, Conowingo Hydroelectric Project, FERC No. 405,
      SS4.1.1/4.3.1, states discharge figures read from USGS 01578310
      (Susquehanna at Conowingo). `daily_extreme_timing` on the approved 15-minute
      record, grouped on America/New_York days, reproduces 141,000 cfs (2010-10-03),
      the morning peaks 26,100 / 46,200 / 80,000 cfs (Studies 9-11) and the Study 1
      prior-day peak 36,500 cfs **exactly**. The 80,000 falls at 08:45 EST, inside
      the report's "0600-0900 h". These are pinned in
      `tests/test_subdaily_published.py` (`requires_network`). The Study 3 prior-day
      peak does not agree (reported 80,900; the record's 2010-05-12 maximum is 81,100)
      and is recorded, not pinned. Searched without finding a usable figure:
      Bevelhimer et al. (2015, RRA) gives no per-gage values and uses 2003-2009
      hourly data; Hells Canyon and Glen Canyon give limits, not observed rates;
      HESS 26:3393 (2022) gives no per-gage diel values. Abernethy et al. (2021,
      Ecosphere) and "Vanishing weekly hydropeaking cycles" (Nat. Commun., 2021,
      doi:10.1038/s41467-021-27465-4) were not accessible (403 or captcha).
      Nothing was found for Methow 12449950 diel timing.
- [x] **No fixed-window ramp metric.** Done in #152: `ramping_rates(window_hours=1)` gives
      2.76 / 4.20 / 5.60 / 3.29 ft on the survey mornings, and Studies 9-11 are pinned. Study 12's
      published 3.1 is the window one reading earlier (3.09 ft), a convention that gives 3.97
      against 4.2 for Study 10, so the report is inconsistent there. Original entry: The same report states hydropeaking as a
      stage decline over "the initial one-hour following reduction to minimum
      flow" (2.8 / 4.2 / 5.6 / 3.1 ft, Studies 9-12), as compliance limits
      usually are ("1 ft per hour"). `ramping_rates` reports only per-interval
      (15-minute) rates, which are not that number. A plain rolling one-hour
      maximum decline on 01578310 gage height gives 2.76 / 4.20 / 5.60 / 3.29 ft,
      so three of four agree at 0.1 ft. A `window_hours=` option on
      `ramping_rates` would make this reproducible from the library. Study 12's
      3.29 vs 3.1 would need explaining first.

### Downstream

- [x] **The app's pin.** flowfreq-app pins v0.17.0 (app PR #8, merged 2026-10-31; before it
      #7 pinned v0.16.0). Each bump is "one edit, deliberate" (`flowfreq-app/CLAUDE.md`):
      review what the release changes in the numbers the app shows, confirm the tag installs,
      and run the app's suite before and after. The app is merged only on the owner's say-so.

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
