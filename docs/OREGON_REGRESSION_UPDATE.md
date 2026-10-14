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
  - East: 67-86 %, against OWRD's 41-75 %.
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
    and leverage reaches 0.67, so the stored W-171003 equations are kept
    (`--hold`).
  - **Willamette:** neither variable is selected, and Waldo Lake stays at -0.52.
    `LAKEAREA` has t = -1.9 there and would cut Waldo's residual to -0.33.

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
| `LAKEAREA_NHD` | NHD high-resolution waterbodies, lakes/ponds and reservoirs ≥ 0.05 km², 5 arc-second mask | percent of basin |

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
of selecting them; W-171003 is held (below).

| Region | Gages | Variables | SEP % (50 / 1 / 0.2 % AEP) | LOOCV SEP % (50 / 1 %) | Pseudo-R² (1 %) |
|---|---|---|---|---|---|
| W (coastal) | 53 | DRNAREA, ELEV, PRECIP | 28 / 36 / 40 | 28 / 33 | 0.98 |
| W-170900 Willamette | 130 | DRNAREA, BSLOPD, JANMAXT, PRECIP | 38 / 43 / 46 | 39 / 40 | 0.96 |
| W-171003 Southern Oregon Coastal | 107 | DRNAREA, BSLOPD, JANMAXT, PRECIP, SOILKSAT | 47 / 49 / 54 | 50 / 47 | 0.95 |
| W-170800 Lower Columbia | 52 | DRNAREA, BSLOPD, SOILKSAT | 50 / 51 / 51 | 52 / 57 | 0.95 |
| E-170701 Middle Columbia | 58 | DRNAREA, ELEV, PRECIP | 89 / 67 / 82 | 93 / 59 | 0.84 |
| E (rest of the east) | 177 | DRNAREA, ELEV, VOLCQP, I24H2Y, LAKEAREA | 90 / 86 / 101 | 93 / 89 | 0.79 |

**What `VOLCQP` and `LAKEAREA` changed** (2026-10-14; regions held at the four
splits above with `--splits`):

| Region | SEP before | SEP after | LOOCV before | LOOCV after | Stored |
|---|---|---|---|---|---|
| E | 94 / 91 / 105 | 90 / 86 / 101 | 97 / 94 | 93 / 89 | **replaced** |
| W-171003 | 47 / 49 / 54 | 57 / 44 / 48 | 50 / 47 | 63 / 46 | kept (`--hold`) |
| W, W-170800, W-170900, E-170701 | unchanged | | | | neither selected |

SEP is at the 50 / 1 / 0.2 % AEPs and LOOCV at 50 / 1 %. In both columns,
"before" is the set that used `SOILKSAT`.

**Region E.** E's 1 % equation is:

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

## Using the set

```python
from flowfreq.regression import evaluate, load_state

lib = load_state("OR", allow_provisional=True)    # logs the PROVISIONAL warning
eq = lib.equation("W-170900", 0.01)               # Willamette, 1 % AEP
est = evaluate(eq, {"DRNAREA": 120.0, "BSLOPD_3DEP": 18.0,
                    "JANMAXT_PRISM9120": 45.0, "PRECIP_PRISM9120": 70.0})
```

**Picking the region.** Use the HUC6 of the basin outlet: 171003, 170701, 170800
and 170900 have their own regions. Otherwise use W for HUC4 1708, 1709 and 1710
and HUC6 180101, and E for everything else. The mapping is also in the stored
file's `notes`.

**Inputs.** Characteristics must come from the sources in the table above.

## Limits and next steps

- The High Cascades residual is only partly resolved. `VOLCQP` separates the
  spring-fed basins, but in W-171003 it costs precision at the 50 % AEP, so it
  is not stored there. Next steps:
  - selecting variables at more than one AEP;
  - a High Cascades volcanic flag built from the DOGAMI state map's `Q`/`QT`
    labels, which would avoid SGMC's age coding (see Data);
  - testing `LAKEAREA` in the Willamette with Waldo Lake's outlet storage in mind.
- Regions follow outlet HUC6s. A basin is in one region, and none straddles
  (plan Review Focus 5 is satisfied by construction).
- East region SEPs stay near 90 %. Eastern Oregon's crest-stage records are short
  and zero-heavy (45 of the 660 dataset gages have zero-flow years), and a log-space GLS cannot
  represent a conditional-probability adjustment.
- A published USGS update for Oregon supersedes all of this.
