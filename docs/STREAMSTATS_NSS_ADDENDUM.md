# Addendum: National Streamflow Statistics (NSS) flow-statistics estimation

**Status:** live-verification pass complete, 2026-09-11; implemented (Phase 2). Region
selection by watershed polygon verified and implemented 2026-09-27 (S5).
**Context:** `docs/STREAMSTATS_MODULE_DESIGN.md` scoped NSS flow-statistics regression
estimation out of Phase 1 ("conflating [characteristics and regression evaluation]
invites an uncited exponent") and required a fresh live-verification pass before any
Phase 2 code, for the same reason that document exists: the `ff-idea02` PDF/py
transcript's NSS-adjacent details were never independently verified, and its
`ss-hydro` endpoint guess had already turned out wrong once. This addendum is that
pass. Every endpoint and failure mode below was exercised against the live service on
2026-09-11, using Phase 1's own `flowfreq.streamstats.delineate_and_get_characteristics`
to supply real basin characteristics (Goat Creek, WA: `48.57426, -120.37893`).

---

## 1. The transcript's own endpoint guesses: wrong, in the same pattern as Phase 1's

- `GET /nssservices/regions/{region}` does **not** carry `statisticGroups`, contrary to
  what the PDF/py's `list_available_statistics` assumed. Confirmed live: it returns
  only `{"id": 51, "name": "Washington", "code": "WA"}`.
- `POST /nssservices/estimate` — the exact path the PDF/py posted to — **does not
  exist**. Confirmed 404 on both `GET` and `POST`, and on every plausible nested
  variant tried (`/nssservices/regions/{region}/estimate`,
  `/nssservices/regions/{region}/statisticgroups/{group}/estimate`, etc.). There is no
  bare `estimate` resource in this API at all.
- No OpenAPI/Swagger spec is published for `nssservices`
  (`/nssservices/openapi.json` 404s), unlike `ss-delineate`/`ss-hydro`. The actual
  machine-readable reference is `GET /nssservices/apiconfig` — a custom, non-OpenAPI
  format describing every resource, method, URI, and parameter. Found via the
  Angular docs viewer at `https://streamstats.usgs.gov/docs/nssservices`, whose
  `assets/config.json` names it (`"apiConfig": ".../nssservices/apiconfig"`). The
  human-readable per-endpoint docs it links to (`code.usgs.gov/.../Docs/...`) live on
  an internal, sign-in-gated USGS GitLab and are **not fetchable without
  authentication** — `apiconfig`'s machine-readable parameter lists, plus live
  experimentation, are the only usable reference.

## 2. The real protocol, verified end to end

```
GET  /nssservices/regions
        -> [{id, name, code}, ...]                         (confirmed identical to
                                                              flowfreq.streamstats.list_regions()'s
                                                              existing assumption)

GET  /nssservices/statisticgroups
        -> [{id, name, code, defType}, ...]                 all possible groups, not
                                                              region-filtered (PC, PFS,
                                                              FVS, LFS, FDS, AFS, MFS,
                                                              SFS, monthly *FDS, GFS, BF,
                                                              FPS, BNKF, YIELD, RCHG, DPS,
                                                              LFCS, LFBS, PROB, UPFS,
                                                              RPFS, BURNFS, ...)

GET  /nssservices/regions/{region}/statisticgroups
        -> [{id, name, code, defType}, ...]                 groups actually valid for
                                                              this region -- WA: only
                                                              PFS (Peak-Flow Statistics)
                                                              and LFS (Low-Flow
                                                              Statistics)

GET  /nssservices/regions/{region}/Scenarios
        ?statisticgroups={code}&unitsystem={1|2|3}
        -> [ { statisticGroupID, statisticGroupName,
               regressionRegions: [
                 { id, name, code, citationID,
                   parameters: [
                     { id, name, description, code, unitType, value: -999.99,
                       limits: {min, max} }, ...
                   ]
                 }, ...
               ],
               links: [{"rel":"Citations","href":".../citations?regressionregions=...","method":"GET"}]
             }, ... ]
        -- a "scenario template": one scenario per requested statistic group, each
           containing one or more independently-calibrated regressionRegions (named
           geographic sub-areas of the state -- WA Peak-Flow has four: GC1750-GC1753,
           each its own equation set and parameter range). unitsystem: 1=Metric,
           2=US Customary, 3=Universal.

POST /nssservices/Scenarios/Estimate
     body: a BARE JSON ARRAY of scenario objects from the GET above, with each
           parameter's "value" replaced by the caller's real basin-characteristic
           value (matched by "code" -- exactly the codes
           WatershedCharacteristics.characteristics already carries: DRNAREA,
           PRECPRIS10, CANOPY_PCT, ...)
        -> the same array shape, each regressionRegion now carrying a "results" array:
           [ { id, name, code, description, value, unit,
               equation: "3.846*DRNAREA^0.745*10^(0.032*PRECPRIS10)/10^(0.0078*CANOPY_PCT)",
               errors: [{id, name, code: "ASEp", value}],   (prediction-error
                                                               STATISTICS, not failures --
                                                               confusingly named)
               intervalBounds: {lower, upper},
               equivalentYears
             }, ... ]

GET  /nssservices/citations?regressionregions={comma-separated ids}
        -> [{id, title, author, citationURL, lastYearOfData}, ...]
        -- resolves a regressionRegion's citationID to a real, DOI-linked reference.
           Confirmed for WA's Peak-Flow regions: Mastin, Konrad, Veilleux & Tecca
           (2016), USGS SIR 2016-5118, DOI 10.3133/sir20165118.
```

**The critical, non-obvious detail: the POST body is a bare array, not
`{"scenarioList": [...]}`.** The `apiconfig` metadata names the body parameter
`scenarioList` (type `array(scenario)`), which reads as an envelope key. It is not —
that name is the server-side model-binding parameter name. Wrapping the array in an
object with that key produces `500 Internal Server Error` with no other diagnostic
(`{"code":500,"message":"An error occured while processing your request. See messages
for more information.","content":"Internal Server Error Occured"}`); posting the same
array bare returns `200` with real, correct results. Confirmed by isolating every
other variable (query params, field subset, single vs. multiple regressionRegions) --
only the envelope-vs-bare-array distinction changed the outcome.

## 3. Failure modes: NSS's own "a 200 is not an answer"

Design doc S4 established the governing principle for Phase 1 from `ss-delineate`'s
behavior. NSS has its own version, found live, and in one respect it is **worse**:

- **No automatic geographic region resolution, and no error when you get it wrong.**
  A state can have several independently-calibrated `regressionRegions` per statistic
  group (WA Peak-Flow: 4). `Scenarios`/`Scenarios/Estimate` compute a result for
  *every* region handed to them, with no check that the region geographically
  contains the point. Confirmed live: submitting Goat Creek's real characteristics
  against all 4 WA peak-flow regions returned 4 different "successful" 50-percent-AEP
  estimates (4,370 / 3,290 / 5,280 / 7,970 cfs) with no error on any of them, though
  only one region is actually correct for that location.
- **`Scenarios/ByLocation`, the one geometry-aware call, does not filter on a bare
  Point.** Tested directly with a GeoJSON `Point` (the only geometry Phase 1 can
  produce -- no watershed polygon is available, per
  `docs/STREAMSTATS_MODULE_DESIGN.md`'s own addendum): all 4 regions came back
  unfiltered, identical to the non-geometry call. A real Polygon/MultiPolygon may be
  required for actual spatial filtering; untested here, since Phase 1 produces none.
  **This is a compounding gap, not a new one**: Phase 2 inherits Phase 1's missing
  polygon as a correctness blocker for automatic region selection, not just a nice-to-have.
- **Out-of-range extrapolation returns HTTP 200 with a plausible, silently-wrong
  number.** Confirmed live: submitting `DRNAREA=999999` (the calibrated range for the
  region used is `[0.25, 3310]` square miles, right there in the same request's own
  `limits`) still returned a computed value (1,450,000 cfs) with **no error, no
  warning**. The one soft signal found: the result's own `errors` list (the ASEp
  prediction-error entry, populated for an in-range estimate) came back **empty** for
  the out-of-range one -- but the value itself is returned as if legitimate regardless.
  **A client must validate every submitted parameter against that same request's own
  `limits.min`/`limits.max` before trusting a result; the service will not do it for
  you.** This is a stricter requirement than Phase 1's `WarningMsg` check: there is no
  message to scan for here at all, only the absence of one signal that isn't
  documented as meaningful.

## 4. What Phase 2 needs to resolve before implementing

1. **Region selection without a polygon.** Given Phase 1 provides no watershed
   geometry and `ByLocation` doesn't filter on a bare point, Phase 2 must either (a)
   require the caller to supply the regression region explicitly (mirroring FR-8's
   "region is a required argument, never inferred" precedent from Phase 1, extended
   one level deeper), or (b) find another resolution path not yet tried here (the
   `RegressionRegions` resource returns metadata but no geometry in what was fetched;
   worth one more look before committing to (a)).
2. **Client-side range validation is mandatory, not optional**, using each request's
   own echoed-back `limits.min`/`limits.max` -- this is the FR-3 analog for this
   module and should be implemented with the same "raise rather than return a
   plausible wrong number" posture.
3. **Citation resolution should be automatic, not a caller afterthought** -- every
   returned estimate should carry (or make trivially available) the resolved
   citation from `/nssservices/citations`, not just the raw `citationID`, closing the
   "uncited exponent" concern design doc FR-9 raised.
4. Low-Flow Statistics (`LFS`, the other group WA supports) was not exercised live
   here -- only Peak-Flow (`PFS`). Worth one confirmation pass before assuming the
   shape generalizes, though there is no structural reason to expect it differs.
   *Done 2026-09-28, S6.* The shape does generalize except for the error code.

## 5. Region selection by watershed polygon (verified live 2026-09-27)

S4 item 1 is resolved. Phase 1 now returns the watershed polygon
(`docs/STREAMSTATS_MODULE_DESIGN.md` S10), and USGS's own workflow notebook
(`https://s3.us-east-1.amazonaws.com/streamstats.usgs.gov/StreamStatsFlowStatisticsWorkflow.ipynb`,
linked from the `ss-delineate` OpenAPI description) resolves regions by POSTing that
polygon to `regressionregions/bylocation`. Everything below was exercised live before
any client code was written.

**Protocol.**

```
POST /nssservices/regions/{region}/regressionregions/bylocation
     body: the bare GeoJSON geometry (Polygon or MultiPolygon), WGS84 [lon, lat]
  -> [ {id, name, code, description, citationID, statusID,
        percentWeight,     <- % of the basin in this region, rounded to a whole percent
        area}, ... ]       <- overlap, mi^2, unrounded
```

Also listed in `apiconfig` with an optional `statisticgroups` query parameter.
`statisticgroups=PFS` and `=LFS` each restricted the answer to that group's regions.
The client sends no filter and matches the answer to the scenario template by `code`.

**Evidence.**

| basin | answer (region-scoped) |
|---|---|
| Goat Creek, WA (412 mi²) | `GC1751` Peak_Region_2, 100% -- **not** `GC1750`, the region the Phase 2 fixtures and live test used |
| Methow nr Pateros, WA (1793) | `GC1751`, 100% |
| Wenatchee at Monitor, WA (1302) | `GC1751`, 100% |
| Green R at Auburn, Snoqualmie nr Carnation, Cedar at Renton, Nooksack at Deming | `GC1752` Peak_Region_3 and `GC1557` Low_Flow_Western_2, 100% each |
| Cowlitz at Randle, Soleduck nr Beaver, Skookumchuck nr Bucoda, Klickitat nr Pitt | `GC1753` Peak_Region_4 (+ `GC1557` west of the crest), 100% |
| Nespelem at Nespelem, Little Nespelem, Antoine Cr | `GC1750` Peak_Region_1, 100% |
| **Ogeechee R nr Louisville, GA (805)** | SIR 2014-5030 regions 1/3/4 at **44/32/24%** (urban `GC1539/1541/1542`; rural-under-1 mi² `GC1572/1573` at 44/32), and `GC1934` (2023 SE US method) 100% |

For Goat Creek the correct answer changes the number: the 50%-AEP peak is **3,290 cfs**
from `GC1751`, not the 4,370 cfs `GC1750` gives.

The Ogeechee weights agree with ss-hydro's own `PCTREG1/3/4` for the same basin
(44.05/31.87/24.08%). The WA regions appear to follow drainage divides: none of the 15 WA
gage basins tried spans two WA regions. A few are less than 0.5% short of 100%, an edge
sliver that NSS still rounds to 100.

**Failure modes, each confirmed live:**

- **A bare Point is refused.** The response is `400 {"message": "Geometry is not of
  type: Polygon,MultiPolygon"}`, so this call cannot do point-only selection.
- **The unscoped `/nssservices/regressionregions/bylocation` mixes studies.** It also
  returns national studies (Bieger 2015, Crippen-Bue, Castro-Jackson) and
  neighbouring states' regions: a Klickitat basin picked up Oregon's SIR 2008-5126
  low-flow regions. The client uses only the region-scoped path.
- **Some regions have no geometry and can never be located.**
  `GET /nssservices/regressionregions/{id}?includeGeometry=true` returns a MultiPolygon
  for WA's four peak regions and for `GC1557`. It returns **no `location`** for the
  low-flow regions `GC1434` (Nooksack, `ELEV1000`), `GC1556` and `GC1558`, all
  `statusID` 3; the located ones are `statusID` 4. Estimating by location therefore
  never returns those three, and `include_unlocated=True` is the way to still get them.
- **Point-only route, verified but not implemented.** The same
  `?includeGeometry=true` call serves each region's polygon (3.5-4.3 MB each), and a
  point-in-polygon test on `GC1751` places Goat Creek and the Methow in it, and Seattle
  out of it. It isn't needed while every delineation yields a polygon, and it would
  miss area weighting.

**NSS's own area averaging.** Verified on the GA rural regions `GC1572`/`GC1573`:

- A `Scenarios/Estimate` body whose regions carry `percentWeight` gets back an extra
  region, `{"id": 0, "code": "areaave", "name": "Area-Averaged"}`. Its values are
  `sum(w_i * Q_i) / 100`: with weights 60/40, PK50AEP is 0.6*116 + 0.4*33.3 = 82.92
  exactly, and the equation reads `"Weighted Average"`.
- It is **silent** when the weights do not sum to 100: 60 + 30 gives no `areaave` and no
  error.
- It is also silent when only some regions carry a weight: 60 + none, or 100 + none,
  gives no `areaave`.

**Client** (`locate_regression_regions`, and `estimate_flow_statistics(...,
watershed_polygon=..., include_unlocated=False)`):

1. The `bylocation` answer is validated against the polygon: every overlap is at most
   the basin's own area plus 2%, and every `percentWeight` is within 1.5 points of
   `100 * area / basin area`. That allows 0.5 for NSS's rounding and 1.0 for projected
   vs spherical area; live the gap is under 0.3%. An answer about some other polygon
   fails this.
2. Only located regions are estimated. Each result carries `located=True` and
   `percent_weight`, and the rest go to `skipped` with the reason. With
   `include_unlocated=True`, the other in-range regions come back too, labelled
   `located=False`, and no average is requested.
3. When one statistic group's located regions are all in range and their weights sum to
   100 (within 0.5), they are sent with `percentWeight`. NSS's `areaave` comes back as
   `RegionFlowEstimates(area_averaged=True)`, and every value is checked against the
   weighted mean to 0.5%. If it is missing, or wrong, the client raises. If a located
   region was out of range, or the weights do not sum to 100, no weights are sent and
   `skipped["{group}:areaave"]` says why.
4. With no polygon the old behaviour is unchanged: every in-range region is returned
   with `located=None`. `batch_estimate_flow_statistics` passes each point's own
   polygon unless `select_by_location=False`.

## 6. Low-Flow Statistics, live (2026-09-28)

**WA's four LFS regions** (`GET /nssservices/regions/WA/Scenarios?statisticgroups=LFS`):

| code | name | parameters (limits) | status |
|---|---|---|---|
| `GC1434` | Low_Flow_Nooksack_Basin_2009_5170 | DRNAREA 1.1-786, ELEV1000 0.091-4.27 ("Elevation in Thousands", 1000 ft) | 3, no geometry |
| `GC1556` | Low_Flow_Western_1_var_2012_5078 | DRNAREA 0.1-48.9 | 3, no geometry |
| `GC1557` | Low_Flow_Western_2_var_2012_5078 | DRNAREA 0.1-48.9, PRECIP 25.1-143 | 4, has geometry |
| `GC1558` | Low_Flow_Western_3_var_2012_5078 | DRNAREA 0.1-48.9, PRECIP 25.1-143, TAU_ANN_G 19-129 | 3, no geometry |

Citations resolve to Curran and Olsen (2009), SIR 2009-5170 (`GC1434`), and Curran, Eng
and Konrad (2012), SIR 2012-5078 (the others).

**`ELEV1000` is unobtainable from StreamStats in WA.** WA's ss-hydro computes 11
characteristics (`GET /ss-hydro/v1/basin-characteristics/WA`), and `ELEV1000` is not
among them. Requesting it returns the -999 "not found" sentinel (design doc S11).
`TAU_ANN_G` is not among them either. So `GC1434` and `GC1558` can only be estimated
from caller-supplied values. This module does not derive `ELEV1000` from `ELEV`,
because whether it means mean-basin or outlet elevation was not verified.

**Live estimate.** Skookumchuck River near Vail (USGS 12025700 site; snapped 46.77233,
−122.59394). Characteristics came live through `BCs=DRNAREA;PRECIP`: DRNAREA 39.9 mi²,
PRECIP 71.07 in.

| region | statistic | NSS value | equation | `evaluate_expression` | error |
|---|---|---:|---|---:|---|
| `GC1556` | M7D10Y | 16.2 ft³/s | `0.15*DRNAREA^1.27` | 16.193 | `SE` 133 |
| `GC1557` | M7D10Y | 12.0 ft³/s | `0.000848*DRNAREA^1.17*PRECIP^1.23` | 11.998 | `SE` 114 |

Both reproduce through `flowfreq.regression.nss`'s own parser to NSS's
three-significant-figure rounding. Neither result carries `intervalBounds`.

**The one shape difference, and a bug it exposed.** The error statistic's code is
`"SE"` ("Average standard error (of either estimate or prediction)"), where PFS uses
`"ASEp"`. The client read only `ASEp`, so every low-flow `standard_error_pct` was
`None`, the same value S3 treats as a soft out-of-range signal. The client now reads
`ASEp`, then `SE`, and records which one in `FlowStatisticEstimate.standard_error_code`.
