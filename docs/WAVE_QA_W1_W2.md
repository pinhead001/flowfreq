# Wave 1 and 2 QA: cross-border consistency, transposition LOOCV, vignette

Roadmap §3.2 sets a per-wave definition of done: a cross-border consistency check, the
transposition LOOCV of §5.2 on the wave's gage network, and a wave vignette. None of
the three had been run for Wave 1 (WA, OR, ID, MT) or Wave 2 (CO, UT, WY, NM, AZ, NV).
This report covers all three. It is QA only: **no equation was changed.**

- Tool: `tools/wave_qa.py` (stages `states candidates b17c delineate locate evaluate dar
  report vignette`, each cached under `.cache/wave_qa/`, so a finished run re-runs
  offline).
- Library code: `flowfreq/validation/wave_qa.py`. Tests: `tests/test_wave_qa.py`
  (offline, with a trimmed real basin in `tests/fixtures/wave_qa/`, plus one
  `requires_network` NSS test).
- Generated tables: `docs/wave_qa/tables.md`, with every row in the CSVs beside it
  (`flagged_pairs.csv`, `flagged_vs_b17c.csv`, `unevaluable.csv`, `dropped_sites.csv`,
  `border_by_pair.csv`, `dar_loocv_by_*.csv`, ...; per-target DAR rows stay in the cache).
- Run 2026-10-04 to 2026-10-06 against `main` at #122, with the NV (#120), UT (#121)
  and AZ region 5 (#116) libraries and the Wave 2 catalog regions (#122).

## Summary

| Check | Ran | Outcome |
|---|---|---|
| Cross-border consistency (per site) | 233 border and straddling site-pairs; **89** evaluable on both sides at AEP 0.01 | **1** site-pair beyond the threshold, at one AEP: Tongue River near Dayton, WY (06298000), Q2 WY/MT = 10.5 against a threshold of 7.3. |
| Cross-border consistency (systematic) | same sample, by state pair | One-sided offsets larger than 26 % at three or more sites: ID above MT (1.83x at 1 %), ID above UT (1.65x at 10 %, 3.1x at 50 %), ID above OR (1.27x/1.48x at 10 %/50 %), WY above MT (1.32x/1.56x at 10 %/50 %), WA below ID (0.68x at 50 %). None is flagged site by site. |
| Equations vs. gage B17C | 260 basins with a usable B17C fit (of 306 evaluated), both states' equations | Home-state equations: bias -4 % to +35 % at AEP 0.01. A neighbour's equations carried across the line: MT -44 % (n 25), WA -57 % (n 12), OR -33 % (n 10), CO +44 % (n 10). 76 site-equation pairs exceed the threshold. |
| Transposition LOOCV (§5.2, drainage-area ratio) | 980 targets in 7 states | Median absolute error 0.15-0.20 log10 (41-58 %) by AEP; RMSE 0.30-0.43 log10. WY, NM and NV not covered (no region geometry for the DAR exponent). |
| Wave vignette | `docs/vignettes/wave1_wave2_ungaged.py` | 9 of 10 states. **NV has no example**: no NV peak-flow region can be located (below). |

**For Wave 3.** By the per-site test the check passes: one flagged pair in 89, at one
AEP, with an identified cause. It cannot be called complete for every border, though.
**Nevada's equations could not be evaluated anywhere**, and Wyoming's and New Mexico's only
on their own side of the line. Those gaps, and the systematic offsets, are listed under
"What blocks a full pass".

## Method

**Sample.** The national catalog (`flowfreq/data/gage_catalog.csv.gz`) supplied every gage
in the ten states with at least 20 years of peaks and a regulation class other than
`regulated` or `urban`: 2,281 gages. Census TIGERweb state outlines (~200 m
generalised) were used for the distances. A gage is a candidate for a state pair if
either:

- it lies within 20 km of the neighbouring state (`border`); or
- its basin could reach the line, i.e. its distance is under 1.2 sqrt(drainage area)
  (`straddle?`).

The 17 shared lines are WA-OR, WA-ID, OR-ID and ID-MT; CO-UT, CO-WY, CO-NM, UT-WY,
UT-NV, UT-AZ, NM-AZ and AZ-NV; and, between the waves, ID-UT, ID-NV, ID-WY, MT-WY and
OR-NV. The Four Corners point pairs (CO-AZ, UT-NM) are left out. That gave 440
candidate site-pairs at 399 gages.

**Basins.** Each candidate was delineated in its own state's StreamStats region; on
failure the neighbour's region was tried. The basin's share in the neighbouring state
comes from a 400-point grid sample against the TIGERweb outline. A basin with any share
there is `straddles`. A `straddle?` candidate that turned out not to cross is `inland`:
it is kept for the B17C comparison and left out of the border tables. 93 gages were
dropped before evaluation (`dropped_sites.csv`):

- 54: delineation failed in both regions (unsnappable points, "OUTSIDE BOUNDARY OF
  DATASET", "Flow statistics not enabled").
- 33: StreamStats `DRNAREA` differed from the catalog's drainage area by more than a
  factor of 1.2. The catalog point had snapped to the wrong stream, e.g. 0.04 mi² for
  the 856 mi² WY 13027500.
- 6: CO basins returned no `DRNAREA`.

**Regions.** For each basin, each state's peak-flow regions were found, with area shares:

- WA, OR, ID, MT, CO, UT and AZ: the StreamStats `nss/regions` MapServer peak layers
  (40, 31, 11, 22, 5, 37, 3). Those polygons are what StreamStats applies; several
  extend past their state line.
  - Arizona's High Elevation region 1 is an overlay, taken when `ELEV` is at least
    7,500 ft, as NSS does.
  - Oregon's western-interior polygon is Region 2, blended between 2A and 2B on mean
    elevation by `regression.oregon`.
- NM: the `HIGHREG` characteristic from an NM delineation. NSS's nine NM regions share
  one statewide polygon.
- WY: the `WYPK_IND` index from a WY delineation. NSS has no WY peak-region geometry.
  `bylocation` returns an unrelated statewide set, GC2098-GC2115, everywhere.
- NV: none (below).

A basin outside every polygon of the neighbouring state takes the polygon nearest its
outlet: "the equation just across the line". This happened in 45 of the 233 cases
(`neighbor_how`).

**Characteristics.** Both states' equations were evaluated on the characteristics of the
one delineation. Where that region does not compute a code the other state's equation
reads, only a stand-in for the *same quantity* from another data vintage is used
(`wave_qa.PROXIES`), and every use is recorded per row:

- `CONTDA`/`DRNAREA`;
- PRISM mean annual precipitation, `PRECIP` / `PRECPRIS10` / `PRECPRIS20`;
- NLCD forest, `FOREST` / `LC11FOREST` / `LC16FOREST` / `LC21FOREST`.

`LAT_OUT`, `LAT_GAGE` and `LONG_OUT` come from the snapped outlet. Nothing else is
substituted. Inputs outside a region's calibrated range are evaluated anyway and listed
(`*_oor`).

**The threshold.** Two regional estimates of one quantile, each with standard error of
prediction SEP (log10), are flagged when

    |log10 Q_a - log10 Q_b| > 1.645 * sqrt(SEP_a^2 + SEP_b^2)

- 1.645 is the two-sided 90 % level that the library's own prediction intervals use.
- Treating the two errors as independent is conservative. Two regressions fitted partly
  on the same border gages have positively correlated errors, which makes their real
  spread smaller. So a flag is a real inconsistency, not noise.
- At typical SEPs of 0.2-0.3 log10, the threshold is a factor of about 2.5-3.5.
- An equation is compared with the gage's B17C quantile in the same way, using
  SEP^2 + var(log Q_B17C), the latter from the B17C 90 % confidence limits.
- Nevada's hybrid regions 6 and 10 publish no SEP. Where they apply, the other side's
  SEP alone sets the threshold (stricter).

A single site can sit inside that band while the whole border is offset. So a second,
descriptive screen marks a state pair **systematic** when all of these hold:

- n is at least 3;
- the median difference exceeds 0.1 log10 (26 %, about half a typical SEP);
- at least 75 % of sites fall on the same side.

**B17C at the gages.** Each gage was fitted with EMA through `workflow.run_ffa` on the bulk
Water Data API peaks, with peak codes applied.

- Skew: the verified PNW regional skew where its HUC8 is covered (WA, OR, ID, advisory
  in western MT). Station skew elsewhere: no Wave 2 state has a verified regional skew
  yet (§1.3).
- 54 of the 2,281 fits failed (below). 56 more are **degenerate**: Q2 < 1 cfs, or Q1%/Q2 > 1,000.
  - These are short, mostly zero-flow records. One 21-year Columbia Plateau tributary
    gives Q2 = 1e-11 cfs and Q1% = 1e15 cfs.
  - Degenerate fits are not used as a reference. Thirty-six of them are in NV.
- B17C quantiles below 1 cfs are not compared either.

## Cross-border consistency

AEP 0.01; all AEPs are in `border_by_pair.csv`. Each pair is listed home-neighbour: the
gage's state first. `ratio` is home/neighbour at the median.

| pair | n | median ratio | median abs diff (log10) | median threshold | share home higher | flagged | systematic |
|---|---|---|---|---|---|---|---|
| ID-MT | 9 | 1.83 | 0.272 | 0.625 | 0.78 | 0 | yes (also 10 %) |
| MT-ID | 22 | 0.94 | 0.143 | 0.625 | 0.41 | 0 | no |
| WA-ID | 11 | 0.67 | 0.234 | 0.615 | 0.36 | 0 | at 50 % only (0.68x, 10 of 11) |
| ID-OR | 4 | 1.10 | 0.058 | 0.417 | 0.75 | 0 | at 10 % and 50 % (1.27x, 1.48x) |
| ID-UT | 3 | 1.27 | 0.105 | 0.928 | 0.67 | 0 | at 10 % and 50 % (1.65x, 3.09x) |
| NM-AZ | 11 | 0.90 | 0.108 | 0.507 | 0.27 | 0 | no |
| NM-CO | 6 | 1.14 | 0.093 | 0.559 | 0.67 | 0 | no |
| WY-MT | 6 | 0.96 | 0.066 | 0.519 | 0.33 | 0 at 1 %, 1 at 50 % | at 10 % and 50 % (1.32x, 1.56x) |
| WY-UT | 4 | 1.02 | 0.023 | 0.466 | 0.75 | 0 | no |
| WY-CO | 3 | 1.05 | 0.114 | 0.545 | 0.67 | 0 | no |
| UT-AZ | 3 | 1.00 | 0.017 | 0.593 | 0.33 | 0 | no |
| AZ-UT, OR-WA, UT-ID, OR-ID | 1-2 each | | | | | 0 | too few |

### The flagged site-pair

**06298000 Tongue River near Dayton, WY (WY-MT, border, 17 km)**

| AEP | WY (region 1, Rocky Mountains) | MT (Southeast Plains) | B17C (92 peaks) |
|---|---|---|---|
| 0.5 | 1,439 cfs | 137 cfs | 1,600 cfs |

WY/MT is 10.5 at AEP 0.5, against a threshold of 7.3; it is inside the threshold at
every rarer AEP.

Cause: **region boundary mismatch.** Montana's SEP plains polygon extends south over this
Bighorn Mountains basin, so the MT side applies a plains equation to a mountain
snowmelt basin. The MT estimate's `FOREST` (61 %, as NLCD 2016) is also above SEP's
calibrated maximum of 57.6 %. The WY mountain equation matches the gage to 10 %.

### Systematic offsets and likely causes

- **ID above MT in the Idaho panhandle (ID gages; 1.83x at AEP 0.01).**
  - Seven of the nine are ID region 1_2 basins (Clark Fork, Coeur d'Alene, St Joe). There
    the MT side is MT's West region, from the region-layer polygon or the nearest one.
  - Against the gages, ID is close (median home-vs-B17C about +0.03 log10) and MT West
    is low (-0.1 to -0.3, and -1.0 at Trapper Creek).
  - Likely cause: **report vintage and calibration domain**. MT West (SIR 2015-5019F)
    is calibrated on drier western Montana; carried to the wetter panhandle it
    under-predicts. On Montana's own side of the same line (MT-ID, n 22) the two states
    agree (median 0.94x).
  - So the jump lies at the panhandle part of the line, not along all of it.
- **WY above MT (10 % and 50 % AEPs).**
  - Mostly the Yellowstone-area basins, where WY region 1 meets MT UYCM and SW. Also
    the Bighorn front (Tongue River above).
  - Cause: **region boundary mismatch** of the MT plains polygons, as above.
  - At AEP 0.01 the pair agrees (0.96x).
- **WA below ID (AEP 0.5, 10 of 11 basins).**
  - These are WA region 1 (Asotin Creek, Palouse) against ID region 3. ID region 3 is a
    one-variable (drainage-area) equation from 13 gages.
  - The WA equation uses `CANOPY_PCT` and `PRECPRIS10`, so the difference at the
    frequent end is the **variable set and calibration sample**.
  - At AEP 0.01 the median is 0.67x but split by basin (36 % WA-higher).
- **ID above UT and OR (10 % and 50 % AEPs, 3-4 basins each).**
  - The UT side is the nearest UT region for southern Idaho basins (Goose and Trapper
    Creeks, region 3, and the Cub River, region 1).
  - The OR side is eastern Oregon region E3, the nearest OR region to the Weiser and
    Salmon basins.
  - Cause: **region boundary mismatch** (the nearest region is the wrong hydrologic
    analogue), plus few gages.
  - They converge at AEP 0.01: 1.27x (UT) and 1.10x (OR).

No Wave 2-internal pair (NM-AZ, NM-CO, WY-CO, WY-UT, UT-AZ) shows a systematic offset.

## Equations vs. gage B17C

This is not a leave-one-out. The published equations are fixed and many of these gages
were in their calibration. It shows how each state's equations reproduce B17C at
border gages, on both sides of the line. AEP 0.01; all AEPs and regions are in
`equation_vs_b17c_by_state.csv` and `regression_by_region.csv`.

| equations | applied in | n | bias | RMSE (log10) | median abs (log10) | beyond threshold |
|---|---|---|---|---|---|---|
| AZ | home | 9 | +11 % | 0.32 | 0.20 | 1 |
| CO | home | 16 | +23 % | 0.27 | 0.19 | 1 |
| ID | home | 50 | +8 % | 0.31 | 0.12 | 6 |
| MT | home | 66 | -4 % | 0.26 | 0.16 | 7 |
| NM | home | 31 | -3 % | 0.24 | 0.15 | 4 |
| OR | home | 12 | 0 % | 0.21 | 0.10 | 1 |
| UT | home | 17 | +35 % | 0.40 | 0.18 | 2 |
| WA | home | 41 | +11 % | 0.25 | 0.17 | 4 |
| WY | home | 18 | +7 % | 0.26 | 0.12 | 2 |
| AZ | neighbour | 22 | 0 % | 0.30 | 0.19 | 2 |
| CO | neighbour | 10 | +44 % | 0.30 | 0.16 | 2 |
| ID | neighbour | 78 | +5 % | 0.24 | 0.13 | 7 |
| MT | neighbour | 25 | **-44 %** | 0.37 | 0.26 | 6 |
| OR | neighbour | 10 | **-33 %** | 0.24 | 0.14 | 2 |
| UT | neighbour | 17 | -13 % | 0.38 | 0.14 | 0 |
| WA | neighbour | 12 | **-57 %** | 0.41 | 0.32 | 5 |

Home-state equations are within their SEPs of the gages. Several states' equations carry
poorly across the line:

- **MT** in north Idaho and Wyoming: the region-layer issue above.
- **WA region 4** on west-slope Oregon Cascade basins (Bull Run, Little Sandy, Youngs
  River): 0.09-0.13x. The nearest WA region to those basins is region 4, the dry east
  Cascades, not region 3.
- **OR E3** on the Clearwater (0.25x).

The 76 site-equation pairs beyond the threshold are listed, with the AEPs and the worst
ratio, in `flagged_vs_b17c.csv` and the generated tables.

Home-state regions whose median absolute error at AEP 0.01 is large:

| state | region | n | bias at AEP 0.01 | note |
|---|---|---|---|---|
| UT | 3 | 3 | +151 % | |
| UT | 6 | 5 | +104 % | Castle Creek near Moab is 23x at AEP 0.5 |
| ID | 6_8 | 17 | +61 % | Pahsimeroi, Medicine Lodge and Little Lost: spring-fed Snake River Plain margin basins |
| ID | 7 | 5 | +83 % | |
| CO | Southwest | 11 | +27 % | |

The UT and ID cases are known hazards of those regions: groundwater-dominated basins, and
few calibration gages (UT region 3 has n = 14). They are not transcription errors; every
stored equation reproduces its report and NSS (status files).

## Not evaluable, and why

At AEP 0.01, among the 233 border and straddling site-pairs (`unevaluable.csv`):

- **NV: no region can be located, on either side (17).**
  - Nevada's regions (WSP 2433) have no geometry anywhere this check can reach. NSS
    `bylocation` returns nothing in the state; there is no `nss/regions` layer; the
    catalog leaves NV blank (#122).
  - **No NV equation was evaluated against any neighbour.** Locating NV regions needs
    the WSP 2433 region map digitised or obtained from USGS.
- **WY as the neighbour (39).**
  - WY's regions are known only through `WYPK_IND`, which StreamStats computes on a WY
    delineation. A basin delineated in MT, ID, UT or CO therefore has no WY region.
  - NSS carries no WY peak-region polygons.
- **NM as the neighbour (14).** The same limitation through `HIGHREG`.
- **Characteristics a state's own StreamStats does not compute (home side).**
  - **OR**: eastern-Oregon variables `JANAVPRE2K`, `JULAVPRE2K`, `WATCAPORC`, `BSLOPD`,
    `STATSGODEP` and `ASPECT`. Also `ELEV`, which the Region 2A/2B split needs. 13
    OR-home basins.
  - **WY**: `JANAVPRE`, which region 5 (Overthrust Belt) reads. 9 WY-home basins.
  - Region 5 is therefore unusable from StreamStats alone, and these are worth an issue
    against the state rollouts.
- **Characteristics the neighbour's region does not compute (neighbour side).**
  - WA's `CANOPY_PCT` (8).
  - CO's `EL7500`, `STATSCLAY` and others (13).
  - MT's `EL6000` (5).
  - ID's `MINBELEV` (3).
  - OR's variables (32).
  - These are **characteristic definition differences**. A state-specific variable has
    no same-quantity stand-in, and none was invented.

## Transposition LOOCV (§5.2)

**Method.** Drainage-area-ratio transposition, as `transpose_frequency` does it:

    log Q_t = log Q_d + b(p) * log10(A_t / A_d)

- `b(p)` is the drainage-area exponent of the target's own regional equation at that AEP.
  The region is the catalog's gage-point region; for AZ, region 1 is skipped for lack of
  mean elevation, and for OR Region 2 the 2B exponent is used.
- Each gage in turn is the target. Its donor is the nearest gage in the same HUC8 with
  0.5 <= A_t/A_d <= 1.5 and within 100 km (median distance 18 km).
- **980 targets** in WA, OR, ID, MT, CO, UT and AZ. WY, NM and NV have no catalog regions,
  so they have no exponent.

| AEP | n | median abs error (log10) | RMSE (log10) |
|---|---|---|---|
| 0.5 | 980 | 0.157 | 0.341 |
| 0.2 | 980 | 0.146 | 0.305 |
| 0.1 | 980 | 0.148 | 0.302 |
| 0.04 | 980 | 0.149 | 0.315 |
| 0.02 | 980 | 0.156 | 0.333 |
| 0.01 | 980 | 0.168 | 0.356 |
| 0.002 | 980 | 0.196 | 0.425 |

By state at AEP 0.01, median absolute error (log10):

| state | median abs error |
|---|---|
| ID | 0.09 |
| OR | 0.15 |
| WA | 0.16 |
| CO | 0.18 |
| MT | 0.20 |
| AZ | 0.23 |
| UT | 0.23 |

The mean bias is near zero by construction and says little. Mutually nearest donors use
the same exponent, so their errors cancel pairwise.

The useful comparison is with the regression equations at the border sample:

| AEP 0.01 | median abs error (log10) |
|---|---|
| DAR from a same-HUC8 donor in the area band | about 0.17 |
| Home-state equations at the border sample | 0.10-0.20 |

So a DAR donor in the band is about as good as the regression equations, not better. It
gets worse toward the rare end (0.20 at AEP 0.002). A weighting of the two is what §5.1's
regression-weighted DAR would test next.

## Vignette

`docs/vignettes/wave1_wave2_ungaged.py` runs offline from `wave1_wave2_sites.json`.

- For each state it takes one real border-area basin, captured from StreamStats by the
  `vignette` stage. The selection prefers a basin that spans two regions.
- It treats the basin as ungaged, evaluates each region's equation with its 90 %
  prediction interval, and area-weights them (flows, as NSS's `areaave` does).
- It prints the gage's B17C Q1% as a check: e.g. the Little Snake near Slater, CO
  (Mountain 0.83 / Northwest 0.17): 3,833 cfs weighted against B17C 4,674.

Nine states have an example. **Nevada has none**, for the reason above.

## Defects and data issues found on the way

None of these was fixed here.

1. **`run_ffa` fails opaquely when peak codes remove every peak.**
   - 29 gages (16 CO, 5 NM, 3 AZ, 3 WA, 1 NV, 1 OR; e.g. AZ 09483042, 09483045,
     09483250) carry code C (urbanization) on every peak. `peak_code_kwargs` removes
     them all.
   - The fit then dies with `zero-size array to reduction operation minimum` instead of
     saying that no peaks remain.
   - The catalog classed these gages `unknown`, not `urban`.
2. **The bulk Water Data API `peaks` query returns duplicate water years at 25 gages.**
   - Examples: 09447000, 13113000 and 13296000.
   - Sometimes the duplicate is a second peak, sometimes the same value twice.
   - `Bulletin17C` rightly refuses a duplicated year, so this tool records a failure.
   - The per-site `download_peak_flow` path was not checked for the same rows.
3. **33 catalog coordinates snap to the wrong stream in StreamStats.** Any automated
   ungaged-site workflow keyed on catalog coordinates needs the drainage-area check used
   here.

## What blocks a full pass

1. **Nevada.** No NV equation could be evaluated: there is no region geometry. The NV
   borders (UT-NV, AZ-NV, ID-NV, OR-NV) are unchecked.
2. **WY and NM region location away from their own delineations.** Their side of a
   border is checked only where the gage is theirs.
3. The **systematic offsets** above, chiefly ID-MT in the panhandle and the MT plains
   polygons over Wyoming mountain basins. They are documented, not smoothed, as §3.2
   asks.
4. **Regional skew** for MT and all of Wave 2 (§1.3): the B17C side used station skew
   there.

By the per-site threshold the waves pass (1 of 89 evaluable pairs, at one AEP). Items 1
and 2 are coverage gaps, not failures. Whether they hold up Wave 3 is a scope decision:
Wave 3's CA borders NV and OR, so NV region geometry will be needed there in any case.

## Reproduce

    python tools/wave_qa.py all                      # ~2.5 h cold, mostly StreamStats
    python tools/wave_qa.py evaluate dar report vignette   # offline from the cache
    python docs/vignettes/wave1_wave2_ungaged.py

Water Data API traffic is about 10 bulk `peaks` requests, well inside the
1,000-per-hour anonymous limit. StreamStats delineations run serially or 2-3 at a time
(the service allows 4). TIGERweb and the StreamStats region layers are fetched once and
cached.
