# Oregon peak-flow regression update (provisional, not USGS)

Oregon regression update plan, tasks 2-5
(`docs/superpowers/plans/2026-10-12-oregon-regression-update.md`), carried out
2026-10-13. **This is not a USGS study.** The equations are stored as a
`provisional` set in `flowfreq/data/regression/provisional/OR.json`, which
`load_state("OR")` never returns. You get it only with
`load_state("OR", allow_provisional=True)`, and it then carries the
`PROVISIONAL (not USGS)` advisory.

## Bottom line

- **Prefer the current equations wherever you can evaluate them.** The provisional
  equations are less precise in every region. At the 1-percent AEP:
  - West (by sub-region): SEP 36-51 %, against Cooper (2005)'s 29-36 %;
  - East: 67-85 %, against OWRD's 41-75 %.
- **What the set offers.** Bulletin 17C station statistics (EMA/MGBT, PNW regional
  skew, peaks through WY2025), and inputs that can all be computed today.
- **Why that matters.** Oregon StreamStats was re-checked on 2026-10-12 and
  2026-10-13 and still does not compute the current western equations' I24H2Y,
  JANMAXT2K, WATCAPORC or SOILPERM. So `OR.json` regions 1, 2A and 2B cannot
  currently be evaluated through StreamStats.
- **A soil-permeability variable helps a little, not enough.** The largest negative
  residuals are spring-fed High Cascades basins:
  - Waldo Lake outlet;
  - Thielsen Creek;
  - Clearwater River above Trap Creek;
  - North Umpqua at Toketee Falls.

  Their floods are small for their size and precipitation. Cooper captured this
  with SOILPERM and WATCAPORC. STATSGO2 saturated hydraulic conductivity
  (`SOILKSAT_STATSGO2`, added 2026-10-13) enters three regions with the expected
  negative sign and lowers their 1-percent SEP by 2-3 points (below). It barely
  moves the High Cascades residuals (Thielsen Creek -0.58 to -0.56 log units,
  Clearwater -0.45 to -0.41, Waldo Lake unchanged at -0.52). All four basins sit
  in the same STATSGO2 map unit (4.2 in./hr, 2.5 times the median), which also
  covers basins that flood normally. What sets them apart is groundwater storage
  in young High Cascades lava, and Waldo Lake's own storage. Neither is a soil
  property at 1:250,000.
- **Geology and lake storage help more, but still not enough.** Two variables were
  added 2026-10-14 (details below):
  - `VOLCQP_SGMC`: percent of the basin on volcanic rock no older than Pliocene;
  - `LAKEAREA_NHD`: percent of the basin under lakes and reservoirs.

  Basins on young volcanic rock:
  - Thielsen 92 %;
  - Clearwater 69 %;
  - Toketee 52 %;
  - Waldo 51 %;
  - the median basin 0 %.

  Waldo is also 32 % lake.

  The results:
  - **Region E:** both variables enter, and SEP falls at every AEP (1 %: 91 to
    86). The stored E equations are replaced.
  - **Southern Oregon Coastal (W-171003):** they enter and nearly halve two High
    Cascades residuals at the 1 % AEP:
    - Thielsen -0.56 to -0.35 log units;
    - Toketee -0.36 to -0.18;
    - Clearwater -0.41 to -0.35.

    The 1 % SEP falls from 49 to 44. But the 50 % AEP SEP rises from 47 to 57,
    and leverage reaches 0.67, so the stored W-171003 equations were kept
    (`--hold`). That was superseded on 2026-10-18 (next bullet).
  - **Willamette:** neither variable is selected, and Waldo Lake stays at -0.52.
    `LAKEAREA` has t = -1.9 there and would cut Waldo's residual to -0.33.
- **Selecting at four AEPs, with Oregon's own geologic map, fixes W-171003
  (2026-10-18).** Two things changed (details in
  [Multi-AEP selection, OGDC-8 and Waldo Lake](#multi-aep-selection-ogdc-8-and-waldo-lake-2026-10-18)):
  - `--select multi` picks one variable set across the 50, 10, 1 and 0.2 % AEPs,
    with a leverage limit of 0.5;
  - `VOLCQP_OGDC` classifies young volcanic rock in Oregon from the DOGAMI state
    map (OGDC-8) instead of SGMC. It puts Waldo at 95 % (SGMC 51 %), Clearwater
    at 87 % (69 %) and Toketee at 64 % (52 %).

  The results:
  - **W-171003** becomes DRNAREA, PRECIP, VOLCQP_OGDC. Its SEP is 47 / 42 / 46
    (50 / 1 / 0.2 % AEP), against 47 / 49 / 54 stored, with LOOCV 49 / 37 and a
    maximum leverage of 0.15. The stored equations are replaced.
  - **E** becomes DRNAREA, I24H2Y, ELEV, BSLOPD, VOLCQP_OGDC: 85 / 85 / 99,
    against 90 / 86 / 101. Replaced.
  - **The other four regions** select their stored sets again, or keep them
    under the rule that a new set must be no worse at every one of those AEPs
    (W-170800).
  - **The High Cascades residuals barely move.** At the 1 % AEP, Thielsen is
    -0.39 log units, Clearwater -0.30 and Toketee -0.36. The gain is spread
    over the region's other basins.
  - **Waldo Lake is not resolved.** `LAKEAREA` stays below |t| = 2 in the
    Willamette (t = -1.9), and Waldo itself carries leverage 0.55 there.
- **Base flow, lake storage and leverage (2026-10-20).** Details in
  [Base-flow index, lake limit and W-170800 leverage](#base-flow-index-lake-limit-and-w-170800-leverage-2026-10-20):
  - **Wolock's (2003) base-flow index is not adopted.** It is high in the High
    Cascades basins (76-81 %, median 59 %). But in W-171003 it is not
    significant at the 1 % AEP (t = -0.1), and the three residuals do not move.
  - **A lake-area limit is adopted** (approved with #166's merge, 2026-10-21). Every equation now
    carries `LAKEAREA_NHD <= 5` percent as an applicability limit. The five
    calibration basins above 5 % average -0.38 log units at the 1 % AEP. Up to
    5 %, there is no trend.
  - **W-170800** becomes DRNAREA, BSLOPD. Its SEP is 53 / 54 / 54, against
    50 / 51 / 51 stored, and its maximum leverage is 0.48 instead of 0.59. The
    stored set's SOILKSAT term rested on Burnt Bridge Creek alone.

## Data

| Step | Tool | Output |
|---|---|---|
| Station statistics (task 2) | `tools/build_oregon_regression_dataset.py` | `data/regression_study/oregon_stations.csv` |
| Basin characteristics (task 3) | `tools/build_oregon_basin_characteristics.py` | `oregon_basin_characteristics.csv` (+ `.sources.json`) |
| Regions, GLS, diagnostics (tasks 4-5) | `tools/oregon_regression_fit.py --store` | `oregon_regression_results.json`, `provisional/OR.json` |

### Gages (task 2)

The candidates are:
- Oregon stream gages;
- WA, ID, NV and CA gages within 25 miles of the Oregon boundary (U.S. Census
  TIGERweb).

The screens:
- the record before the first peak coded 5, 6 or C;
- the GAGES-II regulated or urban class, unless the record was already cut at a
  coded regulation peak;
- at least 10 nonzero peaks;
- drainage area at most 20,000 mi² (this excludes the Columbia at The Dalles).

Each record is fitted by EMA/MGBT, with the skew weighted by the regional skew for
its HUC12:
- the PNW skew at 600 gages;
- the provisional Oregon closed-basin skew (HUC4 1712) at 17;
- station skew alone at 43, all in the Klamath or NV/CA.

**Result:** 660 gages: OR 520, WA 93, ID 24, CA 18 and NV 5.

**Check against USGS.** On the 29 gages it shares with SIR 2016-5083 Table A2,
refitted on that report's terms (WY2013, skew -0.07), the median flowfreq/USGS
quantile ratio is 0.9994-1.0009 at all eight AEPs. The 4 gages off by more than
5 % all have Table A3 perception thresholds, which the check does not apply.

### Basin characteristics (task 3)

Oregon StreamStats returns usable drainage areas and polygons, so the polygon is
used and every other characteristic is computed over it from one national source:

| Code in the equations | Source | Units |
|---|---|---|
| `DRNAREA` | StreamStats, or the NLDI basin (below) | mi² |
| `PRECIP_PRISM9120` | PRISM 1991-2020 normals, annual precipitation, 30 arc-second | in |
| `JANMAXT_PRISM9120` | PRISM 1991-2020 normals, January mean maximum temperature | °F |
| `I24H2Y_NA2` | NOAA Atlas 2 2-year 24-hour precipitation, HDSC western U.S. grid (NAD27) | in |
| `ELEV_3DEP` | USGS 3DEP 1/3 arc-second DEM, basin mean (3DEP ImageServer) | ft |
| `BSLOPD_3DEP` | 3DEP `Slope Degrees` raster function, basin mean | degrees |
| `SOILKSAT_STATSGO2` | NRCS STATSGO2 via Soil Data Access, rasterized at 15 arc-seconds (below) | in./hr |
| `VOLCQP_SGMC` | USGS State Geologic Map Compilation (Horton and others, 2017), young volcanic units, 15 arc-second mask (below) | percent of basin |
| `VOLCQP_OGDC` | In Oregon, DOGAMI Oregon Geologic Data Compilation release 8 (OGDC-8; Darin and others, 2025), young volcanic units; elsewhere the `VOLCQP_SGMC` mask; 15 arc-second (below) | percent of basin |
| `LAKEAREA_NHD` | NHD high-resolution waterbodies, lakes/ponds and reservoirs ≥ 0.05 km², 5 arc-second mask | percent of basin |
| `BFI_WOLOCK2003` | USGS base-flow index grid (Wolock, 2003), 1 km, resampled to 30 arc-seconds (`bfi` stage, added 2026-10-20) | percent |

**Soil permeability (`soils` stage).** `SOILKSAT` is computed as Wolock's (1997)
STATSGO PERMAVE was:
1. Each component's horizon `ksat_r` is averaged, weighted by thickness, over the
   whole profile. Bedrock horizons coded 0 count; horizons with no Ksat do not.
2. Each map unit's value is the mean of its components' values, weighted by
   `comppct_r`. Components with no horizon data (rock outcrop, water) are left out.
3. The result is rasterized from SDA's `gsmmupolygon` polygons, generalized to
   0.002°. Eight slivers that generalized to lines or geometry collections are
   dropped.

gNATSGO was not used:
- it is a 30 m mukey raster plus tables, a far larger download;
- where SSURGO has gaps, which includes much of the Cascade national-forest land,
  it fills them from STATSGO2 or raster soil surveys. Whether that gives the High
  Cascades basins more detail than STATSGO2 has not been checked.

**Geology and lakes (`geology` stage).** Both come from ArcGIS feature services:
- SGMC: energy.usgs.gov, `State_Geological_Map_Compilation`;
- NHD: hydro.nationalmap.gov, `nhd` MapServer layer 12.

Each is rasterized as a 0/1 mask, so a basin's zonal mean is the fraction of its
area covered.

`VOLCQP` takes the SGMC units with:
- `generalize` = "Igneous, volcanic";
- an oldest age (`age_max`) in the Pliocene or Quaternary.

That is a national rule on SGMC's standardized fields. Two limits:
- **Age coding.** Some High Cascades units carry a Miocene `age_max` despite a
  "QT" map label. Oregon's "Mafic vent complexes" (`QTmv`), just south of Waldo
  Lake, are an example. These units are left out, so `VOLCQP` probably
  undercounts young volcanic rock there. Reading the state maps' labels instead
  was not done, because it would not carry across state lines.
- **State lines.** SGMC units are not reconciled across states, so the WA, ID, CA
  and NV basins are classified from their own state maps. Holocene pumice and
  ash deposits mapped as unconsolidated are not counted.

**The Oregon state map (`ogdc` stage, added 2026-10-18).** `VOLCQP_OGDC` answers
the age-coding limit. Its source is the Oregon Geologic Data Compilation,
release 8 (Darin and others, 2025; DOGAMI Digital Data Series OGDC-8), which
supersedes OGDC-7. It is read from DOGAMI's map service,
`gis.dogami.oregon.gov` Public/OGDC MapServer layer 2, whose service
description is "OGDCv8". The rule uses the service's thematic fields:
- `ThematicRockType` is "volcanic rocks", "vent and pyroclastic rocks" or
  "volcaniclastic rocks". The last includes the Mazama pumice and ash.
  Intrusive rocks are left out.
- `ThematicAge` begins no earlier than the Pliocene: Holocene, Pleistocene,
  Quaternary, Pliocene, Pliocene/Pleistocene, Pliocene/Quaternary, or
  Tertiary/Quaternary. In this layer the volcanic Tertiary/Quaternary units are
  27 km² of late High Cascade lavas.
- Ranges that reach back into the Miocene are left out. That covers
  Miocene/Pliocene (the early High Cascade lavas and the Deschutes Formation),
  Miocene/Quaternary and Miocene/Pleistocene.

The service returned 12,491 polygons.

OGDC maps Oregon only. Outside the state (TIGERweb boundary, generalized to
0.001°), each cell takes the `VOLCQP_SGMC` mask's value, so the WA, ID, CA and
NV basins keep the SGMC classification. The two variables correlate at 0.88
over the 579 basins, and 52 basins differ by more than 10 points:

| Basin | `VOLCQP_SGMC` % | `VOLCQP_OGDC` % |
|---|---|---|
| Thielsen Creek (14312700) | 92 | 99 |
| Clearwater above Trap Creek (14314500) | 69 | 87 |
| North Umpqua at Toketee Falls (14315500) | 52 | 64 |
| Waldo Lake outlet (14147000) | 51 | 95 |

Inside the Waldo basin, OGDC maps 94 % as Quaternary "lavas of Cascade
platform", 8 % as Quaternary vents and 6 % as Holocene sediment. SGMC counts
only 51 %. The SGMC units left out there were not inventoried; the `QTmv`
example above suggests they are Miocene-coded.

The largest drops are small basins on Miocene/Pliocene early High Cascade
lava, which SGMC calls Pliocene. One is 14209100 (4 mi²), 100 to 0 %. Another
is 14158950, which is 87 % Quaternary glacial and alluvial sediment over
platform lava. OGDC counts neither.

The OGDC-7 bundle (`OGDC-7_publication-bundle.zip`, Esri file geodatabase) was
downloaded too but not used: OGDC-8 supersedes it, and the venv has no GDAL
vector bindings to read a file geodatabase.

`LAKEAREA` counts NHD FTYPE 390 (lake/pond) and 436 (reservoir) of at least
0.05 km². It enters the regressions as log10(LAKEAREA + 1), stored as a `log10`
term with `offset` 1.

**Why the codes name their source.** A same-named StreamStats characteristic
(PRECIP, ELEV, I24H2Y) comes from other grids, and the codes keep one from being
substituted silently.

**Zonal means.** Each mean covers the cells whose centres fall in the polygon. Two
limits:
- 3DEP is queried at a 30 m pixel, coarsened above ~7,000 mi² (4 basins, where
  slope is biased low);
- the Salmon River at White Bird (13,420 mi²) gets no 3DEP statistics and drops
  out.

**StreamStats mis-snaps.** StreamStats snaps a gage onto its own stream grid and
lands on another stream at roughly a quarter of these gages. For example, Fort Goff
Creek (13 mi²) came back as the Klamath mainstem (10,638 mi²). So a basin is used
only when its area is within 10 % of the NWIS drainage or contributing area:
- StreamStats' polygon first: 470 gages;
- otherwise USGS NLDI's basin for the gage's indexed NHDPlus flowline: 109 gages;
- 81 gages have neither and are left out.

### Regions (task 4)

NSS cannot place a basin in Oregon's current regions:
- every western basin intersects all of Regions 1, 2A and 2B at 100 %;
- OWRD's six eastern regions have no geometry.

So the regions are sets of hydrologic units, assignable from a basin's outlet HUC.
The starting split is West (HUC4 1708, 1709, 1710 and HUC6 180101) and East.
Then, repeatedly:
1. Fit each region.
2. Find the HUC6 with the largest mean GLS residual among those with |t| ≥ 2 and
   at least 25 gages, leaving at least 25 behind.
3. Split that HUC6 off as its own region.

The splits, in order:

| Split | From | Gages | t |
|---|---|---|---|
| 170701 Middle Columbia | E | 58 | 3.6 |
| 171003 Southern Oregon Coastal (Umpqua, Rogue, south coast) | W | 107 | 3.4 |
| 170800 Lower Columbia | W | 52 | 2.7 |
| 170900 Willamette | W | 130 | 2.3 |

The same four splits come out with and without `SOILKSAT` among the candidates;
only their order and t values changed. (Without it, the order was 171003, 170701,
170800, 170900, with t = 4.0, 3.6, 2.9, 2.8.) `--splits` adopts a given list
without re-deriving it.

The West remainder is the coast:
- 171001 Washington Coastal;
- 171002 Northern Oregon Coastal;
- 180101 Northern California Coastal.

Cooper's Region 1 is likewise a coastal region, which is a reassuring
coincidence. No eastern HUC6 besides 170701 qualifies; 170501 (Owyhee/Malheur,
21 gages, t = 1.98) comes closest.

### GLS (task 5)

The WREG procedure (`flowfreq.regression.gls`), specified as in the Washington
refit:
- station SD smoothed by OLS on the region's design;
- skew MSE 0.18 for regional-skew fits and 0.302 for station-only fits;
- n = peaks;
- the PNW concurrent-peak correlation model (SIR 2016-5083 eq. B11) between basin
  centroids, with concurrent years capped at each record's peak count;
- the sampling covariance clipped to the nearest positive-definite matrix where
  the estimator was not.

**Variables.** log10 DRNAREA is forced. The others enter by forward selection at
AEP 0.01 while every added term has |t| ≥ 2 and every VIF < 5. The selected set
is used at all eight AEPs. `--hold` keeps a region's stored variables instead
of selecting them. Since 2026-10-18 the stored set comes from `--select multi`
instead (see
[Multi-AEP selection, OGDC-8 and Waldo Lake](#multi-aep-selection-ogdc-8-and-waldo-lake-2026-10-18)):

```
python tools/oregon_regression_fit.py --splits 170701,171003,170800,170900 --select multi --store
```

The stored set (2026-10-18; W-170800 replaced 2026-10-20):

| Region | Gages | Variables | SEP % (50 / 1 / 0.2 % AEP) | LOOCV SEP % (50 / 1 %) | Pseudo-R² (1 %) |
|---|---|---|---|---|---|
| W (coastal) | 53 | DRNAREA, ELEV, PRECIP | 28 / 36 / 40 | 28 / 33 | 0.98 |
| W-170900 Willamette | 130 | DRNAREA, BSLOPD, JANMAXT, PRECIP | 38 / 43 / 46 | 39 / 40 | 0.96 |
| W-171003 Southern Oregon Coastal | 107 | DRNAREA, PRECIP, VOLCQP_OGDC | 47 / 42 / 46 | 49 / 37 | 0.96 |
| W-170800 Lower Columbia | 52 | DRNAREA, BSLOPD | 53 / 54 / 54 | 56 / 59 | 0.94 |
| E-170701 Middle Columbia | 58 | DRNAREA, ELEV, PRECIP | 89 / 67 / 82 | 93 / 59 | 0.84 |
| E (rest of the east) | 177 | DRNAREA, I24H2Y, ELEV, BSLOPD, VOLCQP_OGDC | 85 / 85 / 99 | 89 / 87 | 0.80 |

Until 2026-10-18 the W-171003 equations used DRNAREA, BSLOPD, JANMAXT, PRECIP
and SOILKSAT (SEP 47 / 49 / 54, LOOCV 50 / 47). The E equations used DRNAREA,
ELEV, VOLCQP_SGMC, I24H2Y and LAKEAREA (SEP 90 / 86 / 101, LOOCV 93 / 89). The
subsections below describe those sets.

**What `VOLCQP` and `LAKEAREA` changed** (2026-10-14; regions held at the four
splits above with `--splits`):

| Region | SEP before | SEP after | LOOCV before | LOOCV after | Stored |
|---|---|---|---|---|---|
| E | 94 / 91 / 105 | 90 / 86 / 101 | 97 / 94 | 93 / 89 | **replaced** |
| W-171003 | 47 / 49 / 54 | 57 / 44 / 48 | 50 / 47 | 63 / 46 | kept (`--hold`); replaced 2026-10-18 |
| W, W-170800, W-170900, E-170701 | unchanged | | | | neither selected |

SEP is at the 50 / 1 / 0.2 % AEPs and LOOCV at 50 / 1 %. In both columns,
"before" is the set that used `SOILKSAT`.

**Region E (the 2026-10-14 set, since replaced).** Its 1 % equation was:

> log10 Q = 2.00 + 0.737 log10 DRNAREA − 0.084 ELEV/1000 − 0.605 VOLCQP/100
> + 1.26 log10 I24H2Y − 0.388 log10(LAKEAREA+1)

- every t is at least 2.4 in absolute value;
- the largest VIF is 1.6 and the maximum leverage 0.28;
- `SOILKSAT` drops out (t = -0.7 once `VOLCQP` is in).

**W-171003.** The selected set is DRNAREA, VOLCQP, I24H2Y, LAKEAREA and BSLOPD.
At the 1 % AEP:
- `VOLCQP` has t = -6.0;
- the SEP falls from 49 to 44 %.

At the 50 % AEP, though, the SEP rises from 47 to 57 % and the LOOCV from 50 to
63 %, and the maximum leverage is 0.67. Because the set does not improve SEP at
every AEP, the stored equations stay as they were.

Hand-built alternatives with `VOLCQP` gave the same trade-off. For example,
DRNAREA, BSLOPD, PRECIP and VOLCQP give 50 / 46 / 50, with LOOCV 53 / 41. A
selection at more than one AEP may be the way out.

**High Cascades residuals** (1 % AEP, log10 units, observed minus predicted):

| Basin | Region | Before | With `VOLCQP`/`LAKEAREA` selected |
|---|---|---|---|
| Thielsen Creek (14312700) | W-171003 | -0.56 | -0.35 |
| Clearwater River above Trap Creek (14314500) | W-171003 | -0.41 | -0.35 |
| North Umpqua at Toketee Falls (14315500) | W-171003 | -0.36 | -0.18 |
| Waldo Lake outlet (14147000) | W-170900 | -0.52 | -0.52 (not selected; -0.33 with `LAKEAREA` forced, t = -1.9) |

The W-171003 "after" values come from the held-back set, so the stored equations
still give the "before" values.

**What `SOILKSAT` changed.** Before it was a candidate, the three regions above
that now use it had these SEPs (50 / 1 / 0.2 % AEP):
- W-171003: 48 / 51 / 55, with DRNAREA, BSLOPD, JANMAXT, PRECIP;
- W-170800: 53 / 54 / 54, with DRNAREA, BSLOPD;
- E: 99 / 94 / 109, with DRNAREA, ELEV, I24H2Y.

Its coefficient is about -0.3 to -0.4 at the 1 % AEP, with t = -2.4 to -2.7.
Elsewhere it is not selected:
- in the Willamette it has t = -0.2 alone, and VIF > 5 beside JANMAXT and PRECIP,
  so it carries the same High Cascades signal as the climate variables;
- in W and E-170701 it has |t| < 2.

Before `VOLCQP` and `LAKEAREA` were added, the best subset of the seven
candidates (every term |t| ≥ 2, VIF < 5) beat the forward selection in no
region.

**Comparison.** The current equations' published SEPs at the same three AEPs:
- Region 1: 27 / 29 / 33;
- Region 2A: 39 / 34 / 38;
- Region 2B: 33 / 36 / 39;
- OWRD E1-E5: 49-63 / 41-49 / 48-55;
- OWRD E6: 104 / 75 / 83.

The comparison is not like for like:
- the gage sets differ;
- the current western equations are B17B-era;
- OWRD's are not USGS.

The plan's LOOCV "against the current equations at the same gages" was not
possible, because those equations' inputs cannot be computed. The leave-one-out
SEPs above are the substitute.

**Leverage.** Highest in W-170800 (0.57 with SOILKSAT, up from 0.48; Lower
Columbia, where 52 gages span 0.2-1,800 mi²). Everywhere else it is 0.28 or less.

## Multi-AEP selection, OGDC-8 and Waldo Lake (2026-10-18)

This section covers the three next steps that "Limits and next steps" listed on
2026-10-14. Regions are held at the four splits with `--splits`.

### Selection at four AEPs (`--select multi`)

The forward selection above works at the 1 % AEP alone. That is why W-171003's
geology set won at the 1 % AEP and lost at the 50 % AEP.

`--select multi` chooses one set for all eight AEPs, as USGS studies do:
1. Every subset of up to five candidates is added to log10 DRNAREA.
2. Sets with a VIF of 5 or more are dropped.
3. Each remaining set is fitted at the 50, 10, 1 and 0.2 % AEPs.
4. The criterion is the mean, over those four AEPs, of the set's average
   variance of prediction divided by the DRNAREA-only model's. So each AEP
   weighs the same, however large its variances.

A set must pass three rules:
- every added term has |t| ≥ 2 at the 1 % AEP;
- every added term has |t| ≥ 2 at more than half the four AEPs;
- no site's GLS leverage exceeds 0.5 at any of them.

The leverage guard is aimed at the W-171003 set of 2026-10-14. Its 0.67 belongs
to Eel Creek at Lakeside (14323300), a 10.6 mi² basin that is 10 % Eel Lake.
The next-largest lake share in W-171003 is Toketee's 1.7 %. So that one basin
gives `LAKEAREA` its slope in the region.

A region already in the stored file keeps its stored variables unless the
selected set's SEP is no worse at every one of the 50, 1 and 0.2 % AEPs.
`VOLCQP_OGDC` is among the candidates (next subsection).

| Region | Selected | SEP selected (50 / 1 / 0.2 %) | SEP stored before | Outcome |
|---|---|---|---|---|
| W-171003 | DRNAREA, PRECIP, VOLCQP_OGDC | 47.3 / 41.9 / 45.9 | 47.4 / 49.0 / 53.7 | **replaced** |
| E | DRNAREA, I24H2Y, ELEV, BSLOPD, VOLCQP_OGDC | 84.9 / 85.3 / 99.3 | 89.8 / 86.5 / 100.8 | **replaced** |
| W-170800 | DRNAREA, BSLOPD | 53.0 / 53.8 / 53.8 | 50.4 / 51.2 / 50.9 | stored kept (worse everywhere) |
| W, W-170900, E-170701 | their stored sets | | | unchanged |

W-170800's stored set has a maximum leverage of 0.59, above the guard, so the
selection cannot pick it. The high-leverage gage is Burnt Bridge Creek at
Vancouver, WA (14211900): 23 mi², the region's flattest basin (1.7°) and its
highest Ksat (6.7 in./hr). The stored set stays because the guard only blocks
new sets. This is an open question. (Resolved on 2026-10-20: the stored set now
faces the guard too, and W-170800 becomes DRNAREA, BSLOPD; see
[W-170800 under the leverage limit](#w-170800-under-the-leverage-limit).)

**The new W-171003 equations.** The 1 % equation is:

> log10 Q = 0.747 + 0.964 log10 DRNAREA + 0.862 log10 PRECIP − 0.900 VOLCQP_OGDC/100

At the 1 % AEP, the t values are:
- DRNAREA 44.0;
- PRECIP 5.7;
- VOLCQP_OGDC -10.9.

The VIF is 1.04 and the maximum leverage 0.15. The LOOCV SEP is 49.5 % at the
50 % AEP and 37.1 % at the 1 % AEP. Before, it was 50 and 47.

BSLOPD, JANMAXT and SOILKSAT drop out. In the old set with VOLCQP_OGDC added,
at the 1 % AEP:
- BSLOPD has t = 0.2;
- JANMAXT has t = -1.1;
- SOILKSAT has t = -2.1, but -1.5 at the 50 % AEP.

That set also has a VIF of 6.5.

The same set with the SGMC variable gives 49.3 / 46.8 / 51.2. So most of the
gain comes from the state map.

**The new E equations.** The 1 % equation is:

> log10 Q = 2.28 + 0.739 log10 DRNAREA + 1.56 log10 I24H2Y − 0.075 ELEV/1000
> − 0.373 log10 BSLOPD − 0.896 VOLCQP_OGDC/100

VOLCQP_OGDC (t = -5.0 at the 1 % AEP) replaces
VOLCQP_SGMC, and `LAKEAREA` drops out. Added to DRNAREA, I24H2Y, ELEV and
VOLCQP_OGDC, LAKEAREA has t = -1.7 at the 1 % AEP and -1.1 at the 50 % AEP. BSLOPD enters with a negative
coefficient: t = -2.4 at the 1 % AEP and -1.7 at the 50 % AEP. ELEV's
coefficient changes sign between the 50 % AEP (t = +3.4) and the 1 % AEP
(t = -2.3), as it did in the stored E set. The maximum leverage is 0.18 and the
LOOCV 88.5 / 86.5, against 93 / 89 before.

### OGDC-8 against SGMC for the High Cascades

With the same companion variables, the OGDC-8 variable fits better than the SGMC
one:

| Region | Set | SEP with VOLCQP_SGMC | SEP with VOLCQP_OGDC |
|---|---|---|---|
| W-171003 | DRNAREA, PRECIP, + volcanic | 49.3 / 46.8 / 51.2 | 47.3 / 41.9 / 45.9 |
| W-171003 | DRNAREA, BSLOPD, PRECIP, + volcanic | 49.7 / 45.7 / 50.4 | 46.9 / 42.8 / 47.3 |
| W-171003 | stored 2026-10-14 set + volcanic | 45.5 / 45.4 / 50.2 | 44.8 / 42.8 / 47.2 (VIF 6.5) |
| E | DRNAREA, I24H2Y, ELEV, BSLOPD, + volcanic | 91.2 / 87.9 / 102.8 | 84.9 / 85.3 / 99.3 |
| W-170900 | stored set + volcanic | 38.4 / 43.6 / 46.6 (t = -1.2) | 38.3 / 43.6 / 46.6 (t = -0.4) |

The four High Cascades residuals did not follow (1 % AEP, log10 units):

| Basin | Stored before (2026-10-14) | DRNAREA, PRECIP, VOLCQP_SGMC | DRNAREA, PRECIP, VOLCQP_OGDC (stored now) |
|---|---|---|---|
| Thielsen Creek | -0.56 | -0.31 | -0.39 |
| Clearwater above Trap Creek | -0.41 | -0.35 | -0.30 |
| North Umpqua at Toketee Falls | -0.36 | -0.39 | -0.36 |
| Waldo Lake outlet (W-170900) | -0.52 | not in the region | -0.52 (VOLCQP_OGDC not selected) |

Thielsen was already 92 % young volcanic under SGMC, so OGDC hardly changes its
value. Toketee's residual falls only with `LAKEAREA`, to about -0.17: Diamond
and Lemolo Lakes lie upstream. That set fails the leverage guard on Eel Creek.
The better region-wide SEP comes from basins whose young-volcanic share OGDC
corrects, not from these four. At about -0.3 to -0.4 log units, these basins'
1 % floods are still only 40-50 % of what the equations predict. Whatever else
sets them apart, such as groundwater-fed baseflow or snowmelt timing, is not in
the candidate set.

### Waldo Lake and `LAKEAREA` in the Willamette

| Set (W-170900, 130 gages) | SEP (50 / 1 / 0.2 %) | t of the lake term (50 / 1 %) | Max leverage | Waldo residual (1 %) |
|---|---|---|---|---|
| stored: DRNAREA, BSLOPD, JANMAXT, PRECIP | 38.5 / 43.5 / 46.4 | | 0.24 | -0.52 |
| + log10(LAKEAREA+1) | 38.2 / 43.2 / 46.1 | -1.9 / -1.9 | 0.55 (Waldo) | -0.33 |
| + LAKEAREA/100 (lake fraction, untransformed) | 37.0 / 41.4 / 43.8 | -3.4 / -3.6 | 0.96 (Waldo) | +0.02 |
| + VOLCQP_OGDC + log10(LAKEAREA+1) | 38.0 / 43.3 / 46.3 | -1.9 / -1.9 | 0.56 (Waldo) | -0.33 |

**`LAKEAREA` is not significant in the Willamette.** On the log scale it stays at
t = -1.9 and lowers the SEP by only 0.3 points. It also puts leverage 0.55 on
Waldo itself.

**The untransformed lake fraction is significant, but only through Waldo.** It
has t = -3.6, but Waldo's leverage is 0.96: the term exists to fit Waldo.
Waldo is 32 % lake; the next-largest share in the region is 4.2 % (North Fork
of the Middle Fork Willamette, 14147500). Its residual then goes to zero by
construction.

Neither set is adopted. Waldo is better treated as a special case, for example
by flagging basins with large natural lake storage as outside the equations'
range, than by a regional term that one gage determines.

### Sources for this section

- Darin, M.H., McClaughry, J.D., Azzopardi, C.J.M., Franczyk, J.J., and Madin,
  I.P., 2025, Oregon Geologic Data Compilation, release 8 (OGDC-8): Oregon
  Department of Geology and Mineral Industries Digital Data Series OGDC-8
  (https://www.oregon.gov/dogami/pubs/Pages/dds/p-OGDC-8.aspx; NGMDB
  https://ngmdb.usgs.gov/Prodesc/proddesc_118650.htm). It was read through the
  map service
  https://gis.dogami.oregon.gov/arcgis/rest/services/Public/OGDC/MapServer (layer
  2) on 2026-10-18.
- Franczyk, J.J., Madin, I.P., Duda, C.J.M., and McClaughry, J.D., Oregon
  Geologic Data Compilation, release 7 (OGDC-7)
  (https://pubs.oregon.gov/dogami/dds/p-OGDC-7.htm). Superseded by OGDC-8 and
  not used.
- U.S. Census Bureau TIGERweb, State_County MapServer layer 0: the Oregon
  boundary.

## Base-flow index, lake limit and W-170800 leverage (2026-10-20)

This section covers three of the open questions that "Limits and next steps"
listed on 2026-10-18. Regions are held at the four splits, and the stored set
comes from the same command as before (`--select multi --store`).

### A base-flow variable for the High Cascades (`BFI_WOLOCK2003`)

**Source.** The `bfi` stage of `tools/build_oregon_basin_characteristics.py`
reads the USGS base-flow index grid for the conterminous United States
(Wolock, 2003, Open-File Report 03-263):
- the 1-km Albers ArcInfo grid `bfi48grd`, from
  `water.usgs.gov/GIS/dsdl/bfi48grd.zip`;
- warped bilinearly to 30 arc-seconds (`bfi_wolock2003_30s.tif`);
- BFI in percent: the base-flow share of streamflow, interpolated between
  gages.

**Values.** The four High Cascades basins sit at the top of the range:

| Basin | BFI % |
|---|---|
| Thielsen Creek (14312700) | 79 |
| Clearwater above Trap Creek (14314500) | 76 |
| North Umpqua at Toketee Falls (14315500) | 76 |
| Waldo Lake outlet (14147000) | 81 |
| median of the 579 basins | 59 |

BFI correlates with ELEV at 0.71 and with `VOLCQP_OGDC` at 0.50.

**Result: not adopted.** BFI is a candidate in every region, and `--select multi`
does not pick it in any of them:

| Region | Set | SEP % (50 / 10 / 1 / 0.2 % AEP) | t of BFI (50 / 1 %) | 1 % residuals: Thielsen / Clearwater / Toketee |
|---|---|---|---|---|
| W-171003 | stored: DRNAREA, PRECIP, VOLCQP_OGDC | 47.3 / 39.3 / 41.9 / 45.9 | | -0.39 / -0.30 / -0.36 |
| W-171003 | stored + BFI | 43.7 / 39.0 / 42.8 / 47.2 | -4.1 / -0.1 | -0.39 / -0.29 / -0.36 |
| W-171003 | DRNAREA, PRECIP, BFI | 49.3 / 50.4 / 56.5 / 60.5 | -10.5 / -5.9 | -0.81 / -0.65 / -0.51 |
| E | stored: DRNAREA, I24H2Y, ELEV, BSLOPD, VOLCQP_OGDC | 84.9 / 71.9 / 85.3 / 99.3 | | |
| E | DRNAREA, I24H2Y, ELEV, VOLCQP_OGDC, BFI (selected, then refused) | 85.1 / 71.9 / 84.3 / 97.3 | -1.5 / -2.6 | |
| W-170900 | stored + BFI | 38.6 / 39.7 / 43.6 / 46.2 (VIF 8.0) | -1.1 / 1.3 | Waldo -0.52 to -0.47 |

- **W-171003.** BFI is significant at the 50 % AEP but not at the 1 % AEP
  (t = -0.1). It lowers the 50 % SEP by 3.6 points and raises the 1 and 0.2 %
  SEPs by about 1 point each. The three residuals do not move. Used in place of
  `VOLCQP_OGDC` it is much worse, so the young-volcanic map carries the
  signal that BFI carries, and more.
- **E.** The selection's best set swaps BSLOPD for BFI. Its SEP is lower at the
  1 and 0.2 % AEPs (by 1.0 and 2.0 points) and higher at the 50 % AEP (by 0.2). So
  the stored set is kept under the no-worse rule. The swapped set's LOOCV SEP
  is also worse at both AEPs checked: 89.0 / 86.8 against 88.5 / 86.5 (50 / 1 %).
- **W, W-170800, E-170701.** |t| < 2 at the 1 % AEP.

**Why it does not help.** The grid is interpolated between gages at 1 km. The
three basins sit at 76-79 %, and so do other basins in the region whose
floods are normal. Their deficit at the 1 % AEP, about -0.3 to -0.4 log units,
is not explained by a smoothed base-flow share. `BFI` stays in the
characteristics file and among the candidates. No other groundwater proxy was
tried: the BFI grid could be obtained, and it was the first candidate.

### Lake-dominated basins: an applicability limit

**Adopted 2026-10-21** (proposed 2026-10-20, approved with #166's merge). Every stored equation carries
`LAKEAREA_NHD <= 5` percent as an applicability limit, which `evaluate` checks
like a variable's range. For a basin over 5 % lake area, `evaluate` raises
`OutOfRangeError`, unless `allow_extrapolation=True`, which computes the value
and records the violation. The value is
`LAKE_LIMIT_PERCENT` in `tools/oregon_regression_fit.py`. To change it, rerun
the fit with `--lake-limit X`; `--lake-limit 0` removes it.

The equation schema gained an `applicability` field for this
(`RegressionEquation.applicability`): ranges on characteristics an equation
does not use. A limit carries no coefficient and no covariance entry.

**Evidence.** The stored equations' residuals were binned by `LAKEAREA_NHD`
(NHD lakes, ponds and reservoirs of at least 0.05 km²), over the 577
calibration basins:
- the table gives each class's mean residual, in log10 units (observed minus
  predicted);
- t is the mean over its standard error, with sites treated as independent.

| Lake area | Basins | 50 % AEP | 10 % AEP | 1 % AEP | t (1 %) |
|---|---|---|---|---|---|
| ≤ 0.5 % | 515 | +0.02 | +0.01 | +0.04 | 3.1 |
| 0.5-1 % | 20 | -0.05 | -0.03 | -0.00 | -0.0 |
| 1-2 % | 28 | -0.01 | -0.02 | -0.03 | -0.6 |
| 2-5 % | 9 | -0.01 | -0.02 | -0.02 | -0.3 |
| **> 5 %** | **5** | **-0.32** | **-0.36** | **-0.38** | **-2.9** |

Up to 5 % there is no trend. Above it, every basin is overpredicted at the
1 % AEP:

| Gage | Region | Lake % | Residual (50 / 1 %) | Standardized (1 %) |
|---|---|---|---|---|
| Waldo Lake outlet (14147000) | W-170900 | 32.1 | -0.48 / -0.52 | -2.6 |
| Odell Creek near Crescent (14055500) | E | 15.0 | -0.25 / -0.19 | -0.6 |
| Grizzly Creek near Lilyglen (11513000) | E | 10.5 | +0.06 / -0.09 | -0.2 |
| Eel Creek at Lakeside (14323300) | W-171003 | 10.1 | -0.64 / -0.79 | -4.0 |
| Cultus Creek above Crane Prairie Reservoir (14051000) | E | 6.7 | -0.27 / -0.32 | -1.0 |

The standardized residual divides by the site's own prediction error (model
error plus sampling variance). Region E's model error is large. So its three
lake basins are each within the noise there, while Waldo and Eel Creek are not.

**Why 5 %.** No calibration basin lies between 4.2 % and 6.7 %:
- 4.2 % is the North Fork of the Middle Fork Willamette (14147500), residual
  +0.03;
- 6.7 % is Cultus Creek.

Any cutoff in that gap splits the gages the same way, and 5 % is a round number
inside it. Five of the 577 gages exceed it. A lower cutoff gains nothing:
- at 3 % it would add two Willamette basins with residuals near zero;
- at 2 % it would add seven more.

**Dropping these gages from the calibration was tested, not stored.** Each
region was refitted with its stored variables (SEP % at 50 / 10 / 1 / 0.2 %):

| Region | Gages | All gages | Without the > 5 % basins |
|---|---|---|---|
| W-171003 | 107 → 106 (Eel Creek) | 47.3 / 39.3 / 41.9 / 45.9 | 43.7 / 33.9 / 35.8 / 40.0 |
| W-170900 | 130 → 129 (Waldo) | 38.5 / 39.3 / 43.5 / 46.4 | 36.7 / 37.2 / 40.9 / 43.4 |
| E | 177 → 174 | 84.9 / 71.9 / 85.3 / 99.3 | 85.5 / 72.2 / 86.0 / 100.3 |

Leaving out the out-of-range gages is consistent with the limit. It lowers
the 1 % SEP by 6 points in W-171003 and by 3 in W-170900. It is not done here,
because it ties the equations to the cutoff, which is not yet approved.

**USGS practice.** USGS reports handle lake storage in two ways:
- **As a variable.** Minnesota's peak-flow regressions use lake and storage area
  as explanatory variables. Ohio StreamStats limits percent storage (NLCD) to
  0-19 % in its valid ranges, and warns outside them.
- **As a scope limit.** Reports define their equations for unregulated streams,
  without large reservoirs, as Cooper (2005) does for western Oregon.

The limit here follows the first pattern: a percent-storage range, checked
like a variable's. Lake area itself does not enter the equations. These
precedents come from a web search on 2026-10-20 and were not checked
against the reports' own tables.

### W-170800 under the leverage limit

**The problem.** On 2026-10-18 the stored W-170800 set (DRNAREA, BSLOPD,
SOILKSAT) had a maximum leverage of 0.59, above the 0.5 limit. The gage was
Burnt Bridge Creek at Vancouver, WA (14211900):
- 23 mi²;
- slope 1.7°, the region's flattest;
- Ksat 6.7 in./hr, the region's most permeable (the next-largest is 4.6).

Without that one gage, SOILKSAT's t falls from -2.4 to -1.0 at the 1 % AEP. So
the term rests on one gage, as `LAKEAREA` rested on Waldo in the Willamette.

**The change.** `choose` now holds a stored set to the same rules as a new one.
A stored set that fails them is replaced by the selection, even when the
selection's SEP is higher. W-170800 is the only region this changes:

| Set | SEP % (50 / 10 / 1 / 0.2 %) | LOOCV SEP % (50 / 1 %) | Max leverage | Pseudo-R² (1 %) |
|---|---|---|---|---|
| DRNAREA, BSLOPD, SOILKSAT (stored before) | 50.4 / 51.6 / 51.2 / 50.9 | 51.8 / 57.2 | 0.59 | 0.95 |
| **DRNAREA, BSLOPD (stored now)** | **53.0 / 53.6 / 53.7 / 53.8** | **56.1 / 59.0** | **0.48** | **0.94** |
| DRNAREA, PRECIP, SOILKSAT (next best) | 53.2 / 54.3 / 54.3 / 54.1 | 56.0 / 57.7 | 0.27 | |

**The cost.** SEP rises by 2.5-2.9 points and the 1 % LOOCV by 1.8. The new
1 % equation is:

> log10 Q = 1.488 + 0.826 log10 DRNAREA + 0.860 log10 BSLOPD

The t values are 20.4 for DRNAREA and 5.1 for BSLOPD. Burnt Bridge Creek still
has the region's largest leverage, 0.48. Its 52 gages span 0.2-1,800 mi², so
the Lower Columbia has the least leverage to spare of any region.

### Stored changes (2026-10-20)

- **W-170800:** all eight AEPs, now DRNAREA and BSLOPD (above).
- **Every equation:** gains the `LAKEAREA_NHD` applicability limit (adopted 2026-10-21).
- **Unchanged:** every other region's coefficients, covariances and error
  statistics are bit-identical to the 2026-10-18 set.

### Sources for this section

- Wolock, D.M., 2003, Base-flow index grid for the conterminous United States:
  U.S. Geological Survey Open-File Report 03-263, 1-km raster
  (https://water.usgs.gov/GIS/metadata/usgswrd/XML/bfi48grd.xml). Retrieved
  2026-10-20 from https://water.usgs.gov/GIS/dsdl/bfi48grd.zip.
- Cooper, R.M., 2005, Estimation of peak discharges for rural, unregulated
  streams in western Oregon: U.S. Geological Survey Scientific Investigations
  Report 2005-5116 (https://pubs.usgs.gov/sir/2005/5116).

## Using the set

```python
from flowfreq.regression import evaluate, load_state

lib = load_state("OR", allow_provisional=True)    # logs the PROVISIONAL warning
eq = lib.equation("W-170900", 0.01)               # Willamette, 1 % AEP
est = evaluate(eq, {"DRNAREA": 120.0, "BSLOPD_3DEP": 18.0,
                    "JANMAXT_PRISM9120": 45.0, "PRECIP_PRISM9120": 70.0,
                    "LAKEAREA_NHD": 0.4})
```

`LAKEAREA_NHD` is required by the lake-area limit (every equation's
`applicability`), even where it is not a variable. A basin over 5 % lake area
raises `OutOfRangeError`.

**Picking the region.** Use the HUC6 of the basin outlet: 171003, 170701, 170800
and 170900 have their own regions. Otherwise use W for HUC4 1708, 1709 and 1710
and HUC6 180101, and E for everything else. The mapping is also in the stored
file's `notes`.

**Inputs.** Characteristics must come from the sources in the table above.

## Limits and next steps

- The High Cascades residual is only partly resolved. The 2026-10-14 next steps
  were done on 2026-10-18:
  - multi-AEP selection;
  - the DOGAMI OGDC-8 volcanic variable;
  - `LAKEAREA` in the Willamette.

  W-171003's 1 % SEP fell from 49 to 42 %. Thielsen, Clearwater and Toketee
  remain at -0.30 to -0.39 log units, and Waldo at -0.52. Three open questions
  were taken up on 2026-10-20:
  - Wolock's BFI grid does not separate the spring-fed basins (not adopted).
  - Lake-dominated basins are flagged by a 5 % `LAKEAREA_NHD` applicability
    limit, adopted on 2026-10-21. Whether to also drop
    the five gages above it from the calibration is open: that lowers the 1 %
    SEP by 6 points in W-171003 and by 3 in W-170900.
  - W-170800's stored set now faces the leverage guard and was replaced.

  Still open:
  - What else separates Thielsen, Clearwater and Toketee, such as snowmelt
    timing or a finer groundwater map than a 1-km interpolated grid.
  - Whether OGDC-8's Quaternary glacial sediment over platform lava should
    count as young volcanic.
- Regions follow outlet HUC6s. A basin is in one region, and none straddles
  (plan Review Focus 5 is satisfied by construction).
- East region SEPs stay near 90 %. Eastern Oregon's crest-stage records are short
  and zero-heavy (45 of the 660 dataset gages have zero-flow years), and a log-space GLS cannot
  represent a conditional-probability adjustment.
- A published USGS update for Oregon supersedes all of this.
