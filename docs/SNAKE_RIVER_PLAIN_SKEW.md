# Snake River Plain regional skew: feasibility (Idaho)

Northwest studies plan, task 8 (`docs/superpowers/plans/2026-10-11-northwest-regional-studies.md`).
This is not a USGS study.

## Decision: not feasible

**A B-WLS/B-GLS regional skew for the Snake River Plain cannot be built from the available
record.** Only **2** gages pass the screens, against the minimum of **20** that was set before
the data were looked at. Idaho's Snake River Plain row in `regional_skew.csv` stays `pending`
with no provisional value. Inside region 0, the guidance remains SIR 2016-5083's (p. 52): the
Pacific Northwest skew does not apply. Consult the USGS Idaho Water Science Center
(Bulletin 17C p. 31).

## Data

`tools/snake_river_plain_skew.py`; inputs in `data/skew_study/srp_inputs.csv`, results in
`srp_results.json`.

- **Idaho peak-flow sites:** 772, from the Water Data OGC API, water years up to 2025.
- **In region 0:** 47.
  - 19 have a 12-digit HUC that `regional_skew_hucs.csv` marks `excluded`.
  - 28 are in a `boundary` HU, placed on the Plain by the gage location against the SIR
    2016-5083 region polygons (StreamStats).
- **Screens:** the same as the Montana v2 study.

| Outcome | Gages |
|---|---|
| Fewer than 25 pre-regulation peaks | 39 (median record length 0, since most never had a peak before a regulation-coded one; only 11 have ≥ 10) |
| No drainage area on record (canals, spillbacks, bypasses, a dry bed) | 6 |
| **Used** | **2**: 13068500 Blackfoot River near Blackfoot (30 peaks, station skew 0.54); 13075983 Spring Creek at Sheepskin Rd (45 peaks, 0.58) |

The Plain's flood records are short, regulated or both. SIR 2016-5083 excluded it for those
reasons, and the data confirm they rule out a study of its own.
