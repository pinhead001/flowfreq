# Verifying the NSS snapshot states

**Goal, set by the owner on 2026-11-03:** as many states as possible with `verified` equations by
mid-November 2026.

**Where it stands at v0.21.0:** 13 jurisdictions are `verified` (AK, AZ, CA, CO, HI, ID, MT, NM,
NV, OR, UT, WA, WY). 38 run on NSS snapshots, and 5 have no equations (AS, DC, GU, MP, VI). So
"verified everywhere" is not reachable by then: the plan is ordered to verify the most states
soonest, not to finish a region at a time. It replaces the wave order in
`docs/MASTER_ROADMAP.md` §3.2 for this goal, and the rule there that a wave's cross-border check
passes before the next wave starts does not apply here. Cross-border checks come after.

## Order

Grouped by likely cost per state, read from each snapshot's regions and citations. No report in
this table has been opened yet, so a state may move once its report has. A region count is what
NSS serves, which can include superseded sets.

| Batch | Why it's cheap or not | States (NSS regions; main report) | Count |
|---|---|---|---|
| 1 | One recent USGS report | MA (1; SIR 2016-5156), CT (2; 2020-5054), ME (2; 2020-5092), KS (2; 2017-5063), PR (2; 2021-5062), VT (2; 2025-5088), PA (5; 2019-5094), VA and WV (11 and 6; both 2025-5110) | 9 |
| 2 | One older USGS report, few regions | RI (1; SIR 2012-5109), NH (2; 2008-5206), ND (3; 2015-5096), NJ (5; 2009-5167) | 4 |
| 3 | Two reports, or a shared report | IA (6; SIR 2013-5086, 2015-5055), MN (7; 2023-5079, 2019-5012), MO (4; 2014-5165, 2010-5073), DE (3; 2022-5005, 2006-5146), AL (6; 2020-5032, 2010-5012), OK (6; 2019-5143, 2015-5134), WI (6; 2022-5118, WRI 86-4005), AR and MS (6 and 4; 2016-5081 and 2018-5148, with shared 2021-5046), GA, NC and SC (5, 4 and 4; shared 2014-5030, 2023-5006) | 12 |
| 4 | Many regions, old, or not a USGS report | NY (12; SIR 2006-5112; 2 equations unparsed), IL (13; Illinois report 23-019), MD (4; University of Maryland panel report), IN (8; Purdue report), SD (7; WRI 98-4055), NE (8; WRI 99-4032), MI (5; WRI 84-4207) | 7 |
| 5 | Several old reports per state | KY (8), OH (9), LA (5; 2 regions cite no report), TN (19), TX (10), FL (7) | 6 |

FL, MI and NE have no StreamStats stream grid, so a verified set there is reachable only by a
caller who supplies the region and characteristics. They are last in their batches for that reason.

## Pace

AK, CA and HI took one working session on 2026-11-01, which is the only measure so far. If that
holds, batches 1 and 2 are reachable by mid-November, for 26 verified jurisdictions, and batch 3
would be a stretch that could reach 38. Report quality can change this, and so can the usage
limits on each session.

## Procedure for each state

"Verified" means what `flowfreq/regression/library.py` says it does: everything the cited reports
publish for the stored regions is stored, double-entered and cross-checked.

1. Get the report, and decide which sets are current. A superseded set NSS still serves stays out
   of scope and is recorded as such, as CA did with WSP 2433.
2. Transcribe the coefficients, ranges, SEP or AVP, and the model error variance, covariance and
   site counts where they are printed. Enter them twice, from the PDF, independently.
3. Compare every equation with the NSS snapshot, and explain each difference.
4. Reproduce the report's worked example to its printed precision. Where it prints none, reproduce
   its per-gage estimates.
5. Write `<STATE>.json` with `status: "verified"` and tests that pin the worked example. Keep the
   snapshot file as evidence, as AK, CA and HI did.
6. Run the state's live smoke point. Live data found three real bugs on 2026-11-03 that no unit
   test had.
7. Update `CHANGELOG.md`, the generated `docs/REGRESSION_COVERAGE.md` and this table.

## Releases

Cut a minor release when each batch finishes, or sooner if a week passes, because a verified set
reaches the app only through a release and a pin bump.

## Out of scope for this goal

The roadmap's per-jurisdiction definition of done also asks for a regional skew, polygons, an
at-site Bulletin 17C re-run and future-flow guidance. A state reaching `verified` under this plan
has its equations checked and nothing more. Those items stay open in `docs/MASTER_ROADMAP.md`.
