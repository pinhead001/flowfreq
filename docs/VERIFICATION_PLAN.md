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

**Where it stands (2026-11-05, batches 1-5 done, released in v0.22.0):** 50 jurisdictions
are `verified`: AK, AL, AR, AZ, CA, CO, CT, FL, GA, HI, IA, ID, IL, IN, KS, KY, LA, MA, MD, ME, MI, MN,
MO, MS, MT, NC, ND, NE, NH, NJ, NM, NV, NY, OH, OK, OR, PA, PR, RI, SC, SD, TN, TX, UT, VA, VT, WA, WI,
WV and WY. DE is `partial`, because of a conflict in its report (below), and 5 have no equations (AS,
DC, GU, MP, VI). The whole order below is finished, ahead of mid-November. This plan replaced the wave
order in `docs/MASTER_ROADMAP.md` §3.2 for this goal; the cross-border checks it deferred are still open.

**Where it stands (2026-11-06, released in v0.23.0):** DE is `verified`, using Table 8 by the owner's decision, making
51. 5 have no equations. AR's Mississippi Alluvial Plain slope was corrected, and the checks below ran.

## What `verified` rests on here

Read this before relying on a batch 2-5 state.

- **Who transcribed:** a separate agent per state did each transcription, and the lead reran a separate
  check per state (range, monotonicity and NSS comparisons) before merging. No person has reviewed them.
- **What "entered twice" meant:** text parsed from the PDF compared with the page images or with NSS,
  not two independent typings.
- **Live checks:** `docs/NATIONWIDE_SMOKE.md` reaches 48 of the 51 at one point each, and the check below
  adds three gages in each of 34 states. FL, MI and NE have no StreamStats stream grid, so their sets
  have never run through a live delineation. None of this shows that every region is right.
- **Calls that were judgment, not proof:** RI's page 25 example does not reproduce (it prints 3,590;
  its equation gives 3,418, and Table 13 prints 3,410 for the same basin); two Table 4 cells at IA's
  05448500 are off by more than rounding; NJ's and MO's weighted columns are not reproduced and not
  stored; 34 NE and 13 NY per-gage flows fall outside tolerance and are pinned as exceptions; MN has six
  A/B differences with no second source; TX region 3's covariance matrix is from another fit and stored
  as printed; DE's Coastal Plain model error variance is Table 8's, not Table 6's, by the owner's decision.
  Each state's `notes` and tests hold the detail.

## Still open

1. **Refusals still on real sites,** from the check below: 18 of 102 gages, mostly points that don't
   snap. IN still needs `region=` (or a region looked up by point), and an urban
   or NRCS-regulated basin does on purpose.
2. **Cross-border checks:** a sample ran (below). The full check by AEP on every border hasn't.
3. **Library follow-ups** the batches found: `TODO.md`, Next steps.

## Wider live check, 2026-11-04

Three catalog gages in each of 34 newly verified states (102 gages: unregulated or unclassified, 20-500
mi2, at least 30 peaks, spread across the state's range of area). Not FL, MI or NE. The gage's Bulletin 17C
fit with station skew was compared with the regression at the gage's coordinates. It ran before the region
rules, at fc24b45, and after them, on the branch merged as #340. StreamStats' ss-hydro failed during the
second run, and the 21 gages it stopped were retried until all had answered.

**Values, after the rules.** With no `region=`, and where the delineation matched the catalog's area within
25 percent (79 gages; OK's area check could not run because OK uses CONTDA, so its gages are counted), the
median ratio of regression to at-site Q1% was 0.92, with the 10th and 90th percentiles at 0.62 and 1.53. The
at-site Q1% fell inside the regression's prediction interval at 69 of 76 gages with an interval. Those gages
probably include many used to fit the equations, so this is a check for gross errors and not an independent
validation. Large differences, none traced to a transcription fault: 05457778 MN (site 16,190, regression
4,086), 06483430 IA (8,763, 2,863), 07377782 LA (5,290, 13,669) and 07301420 OK (2,479, 24,492; Sweetwater
Creek). No gage that already returned a value moved by more than 1 percent.

**Sweetwater Creek, 07301420 (investigated 2026-11-05).** The tenfold gap is real, and no fault was found
in the stored data. Nothing was changed.

- **The regression is what the report's equation gives.** StreamStats locates Region 1 at 100 percent, and
  CONTDA 410.41 mi2 with PRECIP 23.48 in gives 24,492 ft3/s. USGS lists 410 mi2 contributing. The NRCS-regulated
  set doesn't apply: `DAUNREG` equals `CONTDA`.
- **The at-site value holds up.** 40 peaks, 1986-2025, maximum 1,940; Q1% 2,479, 90 percent interval 1,403 to
  11,490. The upstream gage, 07301410 at Kelton, TX (297 mi2), has 64 peaks, 1962-2025, maximum 2,890 in 1977;
  Q1% 2,716, interval 1,938 to 8,212. Their 40 shared years correlate at 0.86 in log. The downstream gage's
  peak code 5 (39 of 40 peaks) doesn't explain it: the Kelton record has no codes.
- **It is far outside what the report saw.** Neither gage is in SIR 2019-5143 Table 1 or SIR 2015-5134 Table
  1. The gap is 0.995 log units, 5.9 standard errors of prediction, and the regression's interval (12,785 to
  46,918) misses the at-site value and its whole interval. Among the 57 Region 1 rows of Table 1, the same
  equation's largest overestimate is 0.563 log units (07327420). Nearby gages don't share the problem: Elm
  Fork near Carl (07303400, 438 mi2) prints 48,800 against the equation's 28,249, and North Fork Red near
  Shamrock (07301300, 817 mi2) prints 25,300 against 35,331.
- **The cause is not established.** A candidate is the High Plains aquifer: the gage is about 2 miles east of
  the 100th meridian, and SIR 2015-5134 reports peaks declining with irrigation, with low gages of this size
  (07157500, 822 mi2, 3,810; 07142300, 475 mi2, 4,880). But Kelton was already low before that report's 1978
  irrigated period (16 peaks, maximum 2,890), and the share of this basin on the aquifer wasn't checked. The
  Panhandle equation would give about 9,840, still four times the at-site value, and StreamStats doesn't
  locate it here.

Near this stream, use the gage record rather than the regression alone. The report's gaged-site weighting and
transfer (equations 18-20) aren't stored. Whether results should warn near this boundary is an open question
in `TODO.md`.

**Refusals.** Without `region=`, 84 gages returned an estimate, up from 53.

| Stopped by | Before | After | Note |
|---|---:|---:|---|
| Overlap: full or drainage-area-only (CT, ME, VT, TN) | 12 | 0 | `region_preference` |
| Overlap: rural or urban (AL, MO), adjusted or not (IL) | 8 | 1 | The AL gage snapped to a 0.04 mi2 tributary, where the small-stream set does apply |
| Overlap: rural or small basin (GA, NC, SC) | 7 | 1 | The NC gage snapped to a 0.18 mi2 tributary, where both apply |
| Overlap: unregulated or NRCS-regulated (OK) | 3 | 0 | `DAUNREG` equals `CONTDA` at all three |
| Overlap: IN (all eight regions), AR (D against MAP) | 3 | 3 | No rule set; IN 2, AR 1 |
| StreamStats "may contain karst" warning (KY) | 3 | 0 | Now a recorded caveat |
| `UnsnappablePointError`: PR (3), TX (2), GA, IL, ND, NH, NJ | 10 | 10 | Gage coordinates; moving the point may work |
| `OutOfRangeError`: MN (LAKEAREA 14.02 against 13.46), WV (0.11 mi2 tributary) | 2 | 2 | |
| StreamStats refused the point (SD) | 1 | 1 | |

Before the rules, the same refusals occurred at v0.21.0 with snapshots (nine gages checked: MO, AL, GA, NC, SC,
IL, OK, IN, KY), so the verified sets had not caused them. Separately, an estimate from two regions returns no
interval (NJ, VA here).

**Rules that can change a result silently.** Where a characteristic is missing, CT, ME, VT, TN and IA fall back
to the drainage-area-only set, as intended. ss-hydro's answers vary between runs, so the same point can get a
different set on another day: TN's smoke gage went from 57,091 (drainage-area-only) on 2026-11-04 to 51,414
(full). The earlier value implies CLIMFAC2YR was missing from that run, though that run did not record it.
`provenance["missing_characteristics"]` records it whenever it happens.

## Checks of 2026-11-05

None of these found a transcription fault in the equations. One found a wrong input code, now fixed.

**Arkansas's alluvial plain slope.** NSS names `CSL1085LFP`, but SIR 2021-5046's slope is `CSLBlue` in
Arkansas: at the 24 gages whose live area matched Table 1 within 15 percent, the median difference was 2
percent against 19, and 17 of 24 were within 5 percent against 5 of 23. Corrected in #345; at 07047942 the old
code cut Q1% by about half. Mississippi's `CSL10_85fm` matches Table 1 exactly at all 5 gages compared.
Louisiana has only 2 comparable gages and one code, and neither matches, so it stays open.

**GA, NC and SC small basins.** SIR 2023-5006 presents "updated methods" for rural streams and defines rural as
under 10 percent impervious, saying it doesn't apply above that. It cites SIR 2014-5030 only for earlier urban
studies, and doesn't say it replaces it. So at 1 mi2 or less both sets apply at or below 10 percent, and only the
2014 sets apply above it. Above 1 mi2, HR1-5's stored `LC19IMP` range already refuses an urban basin. An earlier
version of this section said nothing did, which was wrong. The owner chose HR1-5 first at or below 10 percent,
the 2014 sets above it, and refusal above 1 mi2.

**The other outliers.**

- **MN 05457778, Little Cedar River (at-site 16,190, regression 4,086): real, on the Iowa line.** The at-site
  interval (9,853 to 48,565) and the regression's (1,799 to 9,282) don't overlap. The same gage in Iowa's report
  prints 10,200 from Iowa's equations, and Iowa and Minnesota gages of 37-160 mi2 nearby fit 9,100 to 15,300.
  Minnesota's region D looks low here.
- **IA 06483430, Otter Creek (8,763, 2,863): not a real conflict.** The report prints 2,880 at this gage. The
  at-site interval is 3,741 to 83,432 and contains the regression. The fit leans on two peaks in 1953 and 1954,
  in a record that ended in 1987.
- **LA 07377782, White Bayou (5,290, 13,669): real, cause unknown.** The report prints 12,500, and live slope 6.87
  against its 5.54 explains the rest. The five largest peaks all lie between 4,140 and 4,730, an unusual ceiling.
  The record ends in 2012, before the 2016 flood that nearby Comite River near Olive Branch recorded at 78,000
  against a previous maximum of 25,300.

**FL, MI and NE, without a StreamStats grid.** 38 gages sampled across their regions: stored equations on the
report's printed characteristics against a live Bulletin 17C fit with station skew. Of 34 compared, the median
regression/at-site Q1% was FL 1.25, MI 0.98 and NE 0.82. The largest differences were NE 06465680 (6,595 against 541,
a 13-year record; the report prints 6,590) and FL 02245050 (7,860 against 2,939). NSS located no region for any
of the 31 NLDI basins it was given in these states, though the same request works in KS and OK, so region lookup
can't be tested here either. This checks the equations against gage records. It doesn't test live characteristics,
which still can't be obtained.

**Cross-border sample.** For each adjacent pair of `verified` states, the longest-record gage on each side within
20 km of the line (25+ peaks, 20-800 mi2, not classed regulated or urban): 191 gages, 99 of 107 pairs. Asking
StreamStats in the neighboring state rarely works, with an estimate from both states at only 4 gages. So, as
`docs/WAVE_QA_W1_W2.md` did, the neighbor's equations were evaluated on the home basin's characteristics, in the
neighbor region nearest the gage, allowing extrapolation. That gave 54 gages in 42 pairs. Of the 127 gages with a home
estimate on this pass, 57 lacked a characteristic the neighbor needs, 14 located no neighbor region, and 2 delineated
a basin far from the catalog's area. So this is a sample and not a pass. Rows are in `docs/cross_border_sample_2026-11-05.csv`.

| Neighbor against home, Q1% | Gages |
|---|---:|
| Within 1.5 times | 34 |
| 1.5 to 2 times | 8 |
| 2 to 3 times | 5 |
| 3 times or more | 7 |

The median factor was 1.24. The home equation was nearer the at-site value at 32 of 54, with a median factor from
it of 1.29 against 1.39. Differences are documented, not smoothed. The clearest pattern is Missouri's region 3 against
Kentucky and Tennessee, 3.6 to 5 times lower at the gages near that line. Missouri's equation fits its own gage
(7,220 against 7,373), and Kentucky's and Tennessee's are nearer theirs, so this looks like real contrast across the
Mississippi River rather than an error. Others of 3 times or more: 09414900 AZ/UT, 05052500 ND/MN, 07154400 NM/OK
and 01064500 NH/ME. 4 of the 54 needed a neighbor input out of range.

## Order

All done. Kept as the record of what each batch was expected to cost. A region count is what NSS
serves, which can include superseded sets; several sets it serves were superseded and left out.

| Batch | Why it's cheap or not | States (NSS regions; main report) | Count |
|---|---|---|---|
| 1, done | One recent USGS report | MA (1; SIR 2016-5156), CT (2; 2020-5054), ME (2; 2020-5092), KS (2; 2017-5063), PR (2; 2021-5062), VT (2; 2025-5088), PA (5; 2019-5094), VA and WV (11 and 6; both 2025-5110) | 9 |
| 2, done | One older USGS report, few regions | RI (1; SIR 2012-5109), NH (2; 2008-5206), ND (3; 2015-5096), NJ (5; 2009-5167) | 4 |
| 3, done | Two reports, or a shared report | IA (6; SIR 2013-5086, 2015-5055), MN (7; 2023-5079, 2019-5012), MO (4; 2014-5165, 2010-5073), DE (3; 2022-5005, 2006-5146), AL (6; 2020-5032, 2010-5012), OK (6; 2019-5143, 2015-5134), WI (6; 2022-5118, WRI 86-4005), AR and MS (6 and 4; 2016-5081 and 2018-5148, with shared 2021-5046), GA, NC and SC (5, 4 and 4; shared 2014-5030, 2023-5006) | 12 |
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

A verified set reaches the app only through a release and
a pin bump. v0.22.0 carries all five batches.

## Out of scope for this goal

The roadmap's per-jurisdiction definition of done also asks for a regional skew, polygons, an
at-site Bulletin 17C re-run and future-flow guidance. A state reaching `verified` under this plan
has its equations checked and nothing more. Those items stay open in `docs/MASTER_ROADMAP.md`.
