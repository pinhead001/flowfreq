# Verifying the NSS snapshot states

**Goal, set by the owner on 2026-11-03:** as many states as possible with `verified` equations by
mid-November 2026.

**Where it stands (2026-11-04, batch 1 done, on `main` and not yet released):** 22 jurisdictions are
`verified` (AK, AZ, CA, CO, CT, HI, ID, KS, MA, ME, MT, NM, NV, OR, PA, PR, UT, VA, VT, WA, WV, WY).
29 run on NSS snapshots, and 5 have no equations (AS, DC, GU, MP, VI). So
"verified everywhere" is not reachable by then: the plan is ordered to verify the most states
soonest, not to finish a region at a time. It replaces the wave order in
`docs/MASTER_ROADMAP.md` §3.2 for this goal, and the rule there that a wave's cross-border check
passes before the next wave starts does not apply here. Cross-border checks come after.

**Where it stands (2026-11-04, batches 1-5 done, on `main` and not yet released):** 50 jurisdictions
are `verified`: AK, AL, AR, AZ, CA, CO, CT, FL, GA, HI, IA, ID, IL, IN, KS, KY, LA, MA, MD, ME, MI, MN,
MO, MS, MT, NC, ND, NE, NH, NJ, NM, NV, NY, OH, OK, OR, PA, PR, RI, SC, SD, TN, TX, UT, VA, VT, WA, WI,
WV and WY. DE is `partial`, because of a conflict in its report (below), and 5 have no equations (AS,
DC, GU, MP, VI). The whole order below is finished, ahead of mid-November. This plan replaced the wave
order in `docs/MASTER_ROADMAP.md` §3.2 for this goal; the cross-border checks it deferred are still open.

## What `verified` rests on here

Read this before relying on a batch 2-5 state.

- **Who transcribed:** a separate agent per state did each transcription, and the lead reran a separate
  check per state (range, monotonicity and NSS comparisons) before merging. No person has reviewed them.
- **What "entered twice" meant:** text parsed from the PDF compared with the page images or with NSS,
  not two independent typings.
- **Live checks:** `docs/NATIONWIDE_SMOKE.md` reaches 47 of the 50 at one point each. FL, MI and NE have
  no StreamStats stream grid, so their sets have never run through a live delineation. That one point
  shows the set loads and returns a value, not that every region is right.
- **Calls that were judgment, not proof:** RI's page 25 example does not reproduce (it prints 3,590;
  its equation gives 3,418, and Table 13 prints 3,410 for the same basin); two Table 4 cells at IA's
  05448500 are off by more than rounding; NJ's and MO's weighted columns are not reproduced and not
  stored; 34 NE and 13 NY per-gage flows fall outside tolerance and are pinned as exceptions; MN has six
  A/B differences with no second source; TX region 3's covariance matrix is from another fit and stored
  as printed. Each state's `notes` and tests hold the detail.

## Still open

1. **DE:** Table 6's Coastal Plain model error variances are 12-16 % below Table 8 and the data
   release. While DE is `partial`, its snapshot wins under `allow_nss_snapshot`, the app's path.
2. **Overlapping sets:** tested live on 2026-11-04 at three catalog gages per state. Every point that
   delineated in CT (1), ME (3), VT (2) and TN (3) located overlapping sets and stops with
   `OverlappingRegionsError` unless the caller passes `region=`. NH (2), NY (3) and OH (3) located one set
   at each. Iowa's full-set preference could extend to the first four.
3. **Cross-border checks,** then the rest of the roadmap's definition of done (see Out of scope).
4. **Library follow-ups** the batches found: `TODO.md`, Next steps.

## Order

All done. Kept as the record of what each batch was expected to cost. A region count is what NSS
serves, which can include superseded sets; several sets it serves were superseded and left out.

| Batch | Why it's cheap or not | States (NSS regions; main report) | Count |
|---|---|---|---|
| 1, done | One recent USGS report | MA (1; SIR 2016-5156), CT (2; 2020-5054), ME (2; 2020-5092), KS (2; 2017-5063), PR (2; 2021-5062), VT (2; 2025-5088), PA (5; 2019-5094), VA and WV (11 and 6; both 2025-5110) | 9 |
| 2, done | One older USGS report, few regions | RI (1; SIR 2012-5109), NH (2; 2008-5206), ND (3; 2015-5096), NJ (5; 2009-5167) | 4 |
| 3, done (DE `partial`) | Two reports, or a shared report | IA (6; SIR 2013-5086, 2015-5055), MN (7; 2023-5079, 2019-5012), MO (4; 2014-5165, 2010-5073), DE (3; 2022-5005, 2006-5146), AL (6; 2020-5032, 2010-5012), OK (6; 2019-5143, 2015-5134), WI (6; 2022-5118, WRI 86-4005), AR and MS (6 and 4; 2016-5081 and 2018-5148, with shared 2021-5046), GA, NC and SC (5, 4 and 4; shared 2014-5030, 2023-5006) | 12 |
| 4, done | Many regions, old, or not a USGS report | NY (12; SIR 2006-5112; 2 equations unparsed), IL (13; Illinois report 23-019), MD (4; University of Maryland panel report), IN (8; Purdue report), SD (7; WRI 98-4055), NE (8; WRI 99-4032), MI (5; WRI 84-4207) | 7 |
| 5, done | Several old reports per state | KY (8), OH (9), LA (5; 2 regions cite no report), TN (19), TX (10), FL (7) | 6 |

FL, MI and NE have no StreamStats stream grid, so a verified set there is reachable only by a
caller who supplies the region and characteristics.

## Pace

Batch 1 took one working session on 2026-11-03 and 11-04. Batches 2-5 took one overnight run on
2026-11-04, far faster than this plan expected, because each state went to its own agent in parallel.
That speed is why the caveats under "What `verified` rests on here" matter.

## What the checks have found

Don't skip steps 3 and 4 below: batch 1 found faults in both the reports and NSS.

- **A wrong variable order:** CT's Table 9 states an order for its covariance matrix that gives a
  prediction variance of 1.14. Only one of the 24 orders reproduces the worked example's 0.00286.
- **Labels that disagree:** ME's appendix table 1.8 swaps its 2 and 1 percent labels.
- **NSS errors:** PA region 2 has a transposed coefficient at 4 percent (-0.0831 for -0.0813) and a
  drainage-area limit of 1,160 mi2 for 1,610; PA has two wrong SEPs; PR region 2 has region 1's SEP at
  0.2 percent. KS serves region 2 to 14,600 mi2, where the report's gages stop at 3,555.

- **More NSS errors, batches 2-5:** MN region F serves PMPE limits in inches labeled mm, so it refuses
  every real site; WI region 4 has a limit of 1,304 for 1,034; DE urban has a limit of 117 for 319; AR's
  35 SEPs come from the previous study; NY serves each full region's standard error on its
  drainage-area-only region (60 equations); TN GC2021 has 528 and 0.64 for 540 and 0.643; TX has two wrong citations.

Where a report prints per-gage regression estimates (ME, KS, PA), check them all rather than one
example. That is how PA's coefficient was settled.

## Procedure for each state

"Verified" means what `flowfreq/regression/library.py` says it does: everything the cited reports
publish for the stored regions is stored, double-entered and cross-checked.

1. Get the report, and decide which sets are current. A superseded set NSS still serves stays out
   of scope and is recorded as such, as CA did with WSP 2433.
2. Transcribe the coefficients, ranges, SEP or AVP, and the model error variance, covariance and
   site counts where they are printed. Enter them twice, from the PDF, independently.
3. Compare every equation with the NSS snapshot, and explain each difference.
4. Reproduce the report's worked example to its printed precision, and its per-gage estimates where
   it prints them. Check any covariance against the example, or against the printed average error
   where there is no example, because NSS carries no covariance to compare with.
5. Write `<STATE>.json` with `status: "verified"` and tests that pin the worked example. Keep the
   snapshot file as evidence, as AK, CA and HI did.
6. Run the state's live smoke point. Live data found three real bugs on 2026-11-03 that no unit
   test had.
7. Update `CHANGELOG.md`, the generated `docs/REGRESSION_COVERAGE.md` and this table.

## Releases

Everything since v0.21.0 is unreleased, and a verified set reaches the app only through a release and
a pin bump. v0.22.0 carries all five batches; the owner decides when to cut it.

## Out of scope for this goal

The roadmap's per-jurisdiction definition of done also asks for a regional skew, polygons, an
at-site Bulletin 17C re-run and future-flow guidance. A state reaching `verified` under this plan
has its equations checked and nothing more. Those items stay open in `docs/MASTER_ROADMAP.md`.
