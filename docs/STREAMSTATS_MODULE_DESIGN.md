# Requirements: a StreamStats module for flowfreq

**Status:** Phase 1 (SS3-SS7, this document's actual scope) implemented --
`flowfreq/streamstats.py`, PR #23, 2026-09-10/11. See `TODO.md`'s "Done --
StreamStats module Phase 1" entry for what shipped, and one documented
departure from S6's literal cache-key wording (keyed on the requested, not
snapped, coordinate, so a cache hit needs no network call at all -- NFR-5
taken as the binding requirement over NFR-1's literal phrasing). Phase 2
(NSS flow-statistics estimation, out of scope for this document by design --
see S2) is tracked separately in `TODO.md`'s "Next" section; it needs its
own live-verification pass before any code, for the same reason this
document exists.

**Two implementation bugs the live `requires_network` tests caught and fixed
(2026-09-11), neither in this document's own verified protocol**: the
`pourpoint` snap response's `output` turned out to be a GeoJSON Point
(`{"coordinates": [lon, lat]}`), not flat fields, and an extra
`ss-delineate/v1/delineate/features/{region}` call the implementation added
on top of this document's three-call protocol (to fetch a watershed polygon
FR-3 asks for) was found, live, to return an unrelated zero-area point
feature rather than a basin polygon. Removed rather than fixed -- this
document's own protocol never called it, and getting the real polygon is
still an open question. `WatershedCharacteristics.polygon_geojson` is
therefore always `None` in the shipped implementation; FR-3's `WarningMsg`
check is done by scanning the `sshydro` response instead of checking
polygon geometry.
**Date:** 2026-09-10
**Context:** written from the Methow sub-basin work, where 113 EDT reach pour
nodes need basin characteristics that are currently transcribed by hand from
PDF exports. Every endpoint and failure mode below was exercised against the
live service on 2026-09-09/10; the evidence is in the appendix.

---

## 1. Purpose

Provide programmatic access to USGS StreamStats basin delineation and basin
characteristics, so that a caller with a list of pour points can obtain
drainage area, mean annual precipitation, canopy and the other regression
predictors without manual delineation through the web application.

### Why this belongs in flowfreq, not in a basin repository

Nothing about it is basin-specific. StreamStats covers most of the United
States, the three-call protocol is identical in every region, and the failure
modes below are properties of the service rather than of any watershed. It
sits naturally beside `flowfreq.usgs.USGSgage`, which already wraps NWIS: same
role, different USGS service.

What stays with the caller is which points to delineate, which region code
applies, and what to do about a refused point.

Proposed home: `flowfreq/streamstats.py` — flat package, per the repository
layout in `CLAUDE.md`, beside `flowfreq/usgs.py`.

### Read this first: the same defect is live in the neighbouring module

`flowfreq.usgs.USGSgage.download_daily_flow` has the exact bug this document
is written to prevent, in v0.7.0 today. See §4 item 1. It should be fixed
before or alongside this work — the fix is two lines — but the reason it
survived matters more here than the fix.

`TODO.md` records that NWIS is blocked by the egress proxy in Claude sessions,
and `pyproject.toml` deselects `requires_network` by default. So the one call
that would have caught it is the one call the suite never makes. A StreamStats
module built the same way will hide the same class of defect for the same
reason, which is why §8 asks for recorded fixtures of the *failure* responses
rather than only the happy path.

## 2. Scope

**In scope**
- Delineate a watershed from a coordinate.
- Compute basin characteristics for a delineated watershed.
- Snap a coordinate to the stream network, and refuse when it will not snap.
- Cache, so a delineation is paid for once.
- Record provenance sufficient to reproduce or audit a result later.

**Out of scope for the first version**
- Flow statistics and regression equations. StreamStats exposes these through
  separate services (`nssservices`, `gagestatsservices`); worth a later
  module, but the predictors are the blocking need.
- Editing or refining a delineation, which the web application supports
  interactively and an API caller cannot sensibly do.
- The `batchprocessor` service — see §9.

---

## 3. The service, as it actually is

The historical `streamstatsservices` base path is **retired** and returns 404.
Anything written against `watershed.geojson?includeparameters=true` no longer
works. The web application now calls two services:

| Service | Version | OpenAPI |
|---|---|---|
| `ss-delineate` | 1.2.0 | `https://streamstats.usgs.gov/ss-delineate/openapi.json` |
| `ss-hydro` | 1.4.0 | `https://streamstats.usgs.gov/ss-hydro/openapi.json` |
| `pourpoint` (snap) | — | `https://streamstats.usgs.gov/pourpoint/docs` |

The protocol is two calls, plus a third that should precede them:

```
GET  /pourpoint/v1/snap/str900?region={r}&lat={lat}&lon={lon}
        -> {"output": {...}, "couldSnap": true|false}

GET  /ss-delineate/v1/delineate/sshydro/{region}?lat={lat}&lon={lon}
        -> {"stateAbbreviation": ..., "bcrequest": {...}}

POST /ss-hydro/v1/basin-characteristics/calculate-using-ssdelineate/
        ?region={r}&lat={lat}&lon={lon}
     body: the bcrequest object from the delineate response
        -> [ {name, description, code, unit, value, msg}, ... ]
```

No API key. Roughly 6–9 s for a delineation and 9–11 s for characteristics, so
15–20 s per point end to end. The `region` path parameter is required and is a
state or study-area code (`WA` for Washington); the `region` query parameter on
the `ss-hydro` POST is required too, and omitting it returns a 422 naming the
missing fields.

---

## 4. The governing design principle

**A 200 is not an answer.** This is the requirement everything else serves,
and it is not hypothetical — it is the third instance of the same pattern in
this project inside a week:

1. NWIS returns **only the most recent day** when a daily-values request
   carries no date range, and `USGSgage.download_daily_flow` sends none unless
   given one — `start_date` and `end_date` both default to `None` and are
   omitted from `params`. A one-row record builds a degenerate flow-duration
   curve rather than raising. Measured on 12449500: **1 row** with no range,
   **26,958 rows** (1919-06-01 to 2025-09-30) with one. Cost downstream: a
   full analysis run producing numbers that meant nothing, caught three tiers
   later by a monotonicity check that had no idea what it was really catching.
2. A donor tie in `assemble_target_fdc` was broken by argument order, moving
   the peak tier 17% with no signal at all.
3. **StreamStats delineates a hillslope sliver for an unsnappable point and
   returns HTTP 200** with a structurally valid GeoJSON polygon. The only
   indication is a string in
   `features[0].properties.WarningMsg`:

   > `", Point not snappable using ss-pourpoint API service; results may be inaccurate."`

   Fed onward, that yields a drainage area near zero that is indistinguishable
   from data. Across 113 reaches, a handful of nodes landing off the flowline
   is not a possibility, it is a certainty.

So: **every response is validated against what it is supposed to contain, and
a result that fails validation raises rather than returning.** Convenience
that hides a wrong number is worse than no module at all.

---

## 5. Functional requirements

**FR-1 — Snap before delineating.** Snap the coordinate and refuse when
`couldSnap` is false. Do not delineate an unsnappable point and warn; the
caller cannot act on a warning buried in a log across 113 points.

**FR-2 — Return the snapped coordinate.** The delineation is of the snapped
point, not the requested one, and the distance between them is a data-quality
signal the caller must be able to see and threshold. Report both, plus the
separation in metres.

**FR-3 — Validate the delineation response.** At minimum: `WarningMsg` is
absent or empty, the `globalwatershed` feature exists and carries a polygon,
and the polygon is not degenerate. A response failing any of these raises.

**FR-4 — Never send a request with unstated parameters.** Every query
parameter that affects the result is passed explicitly, even where the service
has a default. This is the direct lesson of the NWIS defect: a default that
changes, or that was never what you assumed, is silent.

**FR-5 — Return characteristics as a typed, indexable result**, keyed by
StreamStats code (`DRNAREA`, `PRECPRIS10`, `CANOPY_PCT`, …), each carrying
value, unit, description and any service message. Do not flatten to a bare
dict of floats: the units and messages are how a caller notices that
`DRNAREA` came back in the wrong units or with a caveat attached.

**FR-6 — Surface per-characteristic messages.** The service attaches a `msg`
to individual values (for example a note about local versus total area). These
must reach the caller.

**FR-7 — Fail one point without failing the batch.** A caller delineating 113
points needs the 110 that worked and a clear account of the 3 that did not,
not an exception on point 4. Provide a batch entry point that returns results
and failures side by side, with the reason per failure.

**FR-8 — Region is a required argument, never inferred.** Do not guess a
region from a coordinate. A wrong region silently returns a delineation
computed against the wrong data layers.

**FR-9 — Do not compute regression flows.** The module returns predictors.
Evaluating a regional regression is a separate concern with its own citation
requirements, and conflating them invites an uncited exponent.

## 6. Non-functional requirements

**NFR-1 — Cache by default, keyed on (region, snapped lat, snapped lon,
service versions).** Delineation is expensive on USGS servers and slow for the
caller. A cached result must be reproducible and must not be silently reused
across a service version change.

**NFR-2 — Be a polite client.** Serial by default with a configurable small
concurrency ceiling; exponential backoff on 5xx and timeouts; a descriptive
`User-Agent` identifying the tool and a contact. This is a free public service
with no key, and 113 delineations is a real load.

**NFR-3 — Generous, configurable timeouts.** Observed 6–11 s per call under
good conditions; a large or complex basin will be slower. A default that is
too tight turns a working delineation into a failure.

**NFR-4 — Record provenance with every result.** Service versions, the request
URLs and parameters, the UTC timestamp, and the snapped coordinate. StreamStats
updates its underlying data layers, so a characteristic obtained today is not
guaranteed to be reproducible from the same coordinate next year. Without a
stamped provenance record a published number cannot be defended.

**NFR-5 — Work offline from cache.** A caller with a populated cache must be
able to complete a run with no network, so an analysis can be re-run and
audited where the service is unreachable.

**NFR-6 — No hard dependency on a GIS stack.** Return the watershed polygon as
GeoJSON. Do not require `geopandas` or `shapely` to obtain a drainage area.

## 7. Errors

Distinct, catchable exception types, because the caller's response differs:

| Condition | Caller's response |
|---|---|
| Point will not snap | Move the node; a data problem, not a transient one |
| Point outside the region / region unsupported | Fix the region code |
| Delineation returned a degenerate or warned result | Investigate the node |
| Transport failure, timeout, 5xx | Retry |
| Malformed or unexpected response shape | Bug in this module or a service change |

The last row matters: a service that changes its response shape must produce a
loud failure, not a silently empty result.

## 8. Testing

- **No network in CI**, consistent with `requires_network` being deselected by
  `addopts`. Recorded fixtures under `tests/fixtures/` covering a good
  delineation, an unsnappable point, a 422, and a malformed response. The
  failure fixtures are the ones that matter: a suite that only records the
  happy path is how the `download_daily_flow` defect stayed alive.
- **The unsnappable-point fixture is mandatory** and asserts that the module
  raises. That test is the whole point of the module.
- **One opt-in live test** marked `requires_network`, checking the two Methow
  points below against their known values — so a service change is detected
  deliberately rather than discovered in someone's results. It will not run in
  a Claude session; that is expected, and is why the fixtures above carry the
  real assertions.
- Cache round-trip, provenance completeness, and batch partial-failure.

## 9. Open questions

1. **`batchprocessor`.** `streamstats.usgs.gov/batchprocessor` exists and has a
   docs endpoint; it was not evaluated. If it accepts a list of points and
   returns characteristics, much of §5 could be a thin wrapper instead. Worth
   an hour before building anything.
2. **Region codes.** Where does the authoritative list live, and should the
   module validate against it or pass through?
3. **Snap tolerance.** Should a snap that moves a point more than some distance
   be refused rather than accepted? `couldSnap` is binary and the service will
   happily move a point a long way.
4. **`bcLabels`.** The delineate response requests `"*"` (all characteristics).
   Should a caller be able to request a subset for speed?

## Appendix — verification evidence

Exercised 2026-09-09/10 against the live service, region `WA`.

Both hand-delineated Methow test points reproduce their StreamStats PDF
exports exactly:

| point | lat, lon | DRNAREA | PRECPRIS10 | CANOPY_PCT |
|---|---|---:|---:|---:|
| `tp_goat_creek` | 48.57426, −120.37893 | 412.0 / 412 | 45.62 / 45.62 | 45.242 / 45.242 |
| `tp_lost_river` | 48.65041, −120.51172 | 252.0 / 252 | 47.99 / 47.99 | 43.090 / 43.09 |

An off-network point (48.584, −120.370, about 1 km onto the hillslope):

- `pourpoint` snap returns `"couldSnap": false`.
- `ss-delineate` returns **HTTP 200**, a 1,297-byte body against 184,417 for
  the valid point, a structurally valid polygon of a hillslope sliver, and the
  `WarningMsg` quoted in §4.

Retired endpoints confirmed 404: `/streamstatsservices/watershed.geojson`,
`/streamstatsservices/parameters.json`.
