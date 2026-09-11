# Addendum: National Streamflow Statistics (NSS) flow-statistics estimation

**Status:** live-verification pass complete, 2026-09-11. Implementation not yet started.
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
