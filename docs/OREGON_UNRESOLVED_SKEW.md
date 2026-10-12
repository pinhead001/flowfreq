# Oregon's closed basins and the Klamath: does the Pacific Northwest skew apply?

Northwest studies plan, task 9 (`docs/superpowers/plans/2026-10-11-northwest-regional-studies.md`).
This is not a USGS study.

## The question

SIR 2016-5118 (pp. 1, 23) defines the Pacific Northwest (PNW) regional skew study region as:
- the Columbia basin;
- the Oregon-Washington coast;
- the rest of Washington.

That leaves out HUC4 **1712** (Oregon closed basins) and **1801** (Klamath). Yet 5 of the
study's 290 gages are there (SIR 2016-5083 Table B1). `regional_skew_hucs.csv` therefore left
both units `unresolved`, and the lookup refused them.

## Test

`tools/pnw_edge_skew_test.py`; inputs in `data/skew_study/pnw_edge_inputs.csv`, results in
`pnw_edge_results.json`.

- Every Oregon peak-flow gage in 1712 and 1801 (68), from the Water Data OGC API,
  water years up to 2025.
- Screens as in the Montana v2 study: the pre-regulation record, at least 25 peaks, GAGES-II
  regulated/urban gages out unless the record was cut at a coded regulation start, nested
  pairs out. Excluded: 53 short records, 4 redundant, 1 with no drainage area.
- A CONSTANT B-WLS/B-GLS skew is fitted with the PNW correlation-distance model (SIR
  2016-5083 eq. B11). It is compared with G = -0.07 by t = (β + 0.07) / SD(β).
- Decision rule, fixed in the plan beforehand: |t| < 2 means consistent.

| Sample | Gages | β (SD) | t vs -0.07 | Decision |
|---|---|---|---|---|
| 1712 + 1801 | 10 | -0.107 (0.188) | -0.20 | consistent |
| 1712, closed basins | 8 | -0.195 (0.210) | -0.59 | consistent |
| 1801, Klamath | 2 | not fitted | | too few |

The gages are:
- **Closed basins:** Twentymile, Deep and Honey creeks, Chewaucan, Donner und Blitzen, Bridge
  and Silver creeks, Trout Creek.
- **Klamath:** Sprague River near Beatty, Williamson River below Sprague.

## Outcome

- **1712 (Oregon closed basins):** a provisional member of a new region, "Pacific Northwest
  (Oregon closed basins, provisional)". It carries the PNW value G = -0.07, MSE 0.18, is
  valid in Oregon only, and is returned only with `allow_provisional=True`.
  **Caveat:** with SD(β) ≈ 0.2, the test cannot detect a difference smaller than about 0.4.
  "Consistent" means the data give no reason to reject -0.07, not that the value is confirmed.
- **1801 (Klamath):** stays `unresolved`. Two gages cannot test anything. This departs from
  the plan's combined-sample rule, deliberately: the combined t is driven by the closed basins.
- **Defaults are unchanged.** `regional_skew_for("OR")` still returns the verified
  state-keyed PNW row, as the OR row's notes intend. A provisional row is invisible to
  lookups that do not opt in.
