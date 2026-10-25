# Nationwide rollout by mid-November 2026: design

**Goal.** By 2026-11-15, flowfreq gives a documented flood-frequency estimate at any gaged
or ungaged site in the 56 jurisdictions, and the Streamlit app exposes both paths. "Any
site" means: Bulletin 17C at a gage (done today), a regional-regression estimate at an
ungaged site from an offline library covering every jurisdiction NSS serves, and an honest
"not available, and why" everywhere else.

**Decisions taken 2026-10-25** (user approved all four):

1. A new opt-in library status, `nss_snapshot`, carries equations harvested from the USGS
   NSS web service and parsed by flowfreq, but not checked against the source report. It is
   "available with caveat", never `verified`.
2. The methow-flow tier 2 items #4 (two-gage interpolation) and #5 (Region 2 refit) and all
   of tier 3 are parked until after the deadline.
3. The app gets the v0.15.0 pin bump (reviving the closed PR #5 work) and an ungaged-site
   page. The user merges and deploys app PRs.
4. 0.15.0 is released from main as-is (PR #191).

## Baseline (2026-10-25)

- Library v0.14.1 released; main 0.15.0.dev0, release PR #191 open.
- Regression library: 10 of 56 jurisdictions `verified` (Waves 1 and 2); 46 `not started`.
- Regional skew: verified WA/OR/ID/AZ; provisional MT; pending CO/UT/WY/NV/NM; no rows
  for the other 46.
- No `estimate_at_site` entry point (roadmap §4 all open). The pieces exist:
  `streamstats.delineate_and_get_characteristics`, `streamstats.locate_regression_regions`,
  `regression.evaluate` / `evaluate_weighted`, `gage_weighting.weight_with_regression`.
- flowfreq-app pins v0.7.0; the v0.10.1 bump (branch `bump-flowfreq-v0.10.0`, PR #5,
  +951 lines) was closed unmerged on 2026-10-11. No regression or ungaged UI.

## NSS sweep (2026-10-25, `tools/snapshot_nss.py`, group PFS)

41 of the 46 jurisdictions have NSS peak-flow regions. DC, VI, GU, MP and AS have none.
Indiana's 8 regions resolve to no equations because its `HIGHREG` gating parameter is
bounded (1005-1012) and the tool only guesses gating values for unbounded parameters.

| Equations | Count |
|---|---|
| harvested | 1793 |
| parsed and round-tripped | 1573 |
| refused by the parser | 220 |

Refusals, by cause:

| Cause | States (equations) | Fix |
|---|---|---|
| identifier with lowercase letters (`CSLBlue`, `CSL10_85fm`, `BSLDEM10ff`, `SOILCorD`) | AR 7, TX 34, MS 8, OK 14, MD 30 | relax `_VARIABLE` |
| `logN(X, 10)` / `LogN(...)` in the exponent of 10, with `(X+1)` arguments, identity terms and scaled identity terms | MN 48, WI 32, MA 8 | parse `logN` in `_linear` |
| negated affine base `(-1*(BDF-13))^(-b)` | FL 3, KY 6, OH 6 | affine base with `scale=-1`, `offset=13` |
| region-indicator interactions `DRNAREA^(0.646+0.00224*PCTREG2)` | GA 8, NC 8, SC 8 | **out of scope**; recorded as unparsed |

The snapshot files are in the session scratchpad (`nss_sweep/`) and must be copied to
`data/nss_snapshots/` by the data task; Indiana is re-harvested after the gating fix.

## Design

### A. The `nss_snapshot` tier

Follows the existing opt-in pattern in `flowfreq/regression/library.py` exactly:

- Files live in `flowfreq/data/regression/nss_snapshot/<CODE>.json` with
  `status: "nss_snapshot"`. `STATUSES` and `OPT_IN_SUBDIRS` gain the status.
- `load_state(code, allow_nss_snapshot=True)` returns it, with `advisory` set and a logged
  warning. Precedence when several flags are set and files exist: `awaiting_nss`, then
  `provisional`, then `nss_snapshot`. A generator never writes an `nss_snapshot` file for a
  jurisdiction that has a published file, so in practice there is no overlap.
- `StateLibrary.nss_snapshot` property; `nss_snapshot_states()` listing;
  `StateLibrary.unparsed: List[Dict[str, Any]]` from a new optional JSON key `unparsed`
  (entries: `region_code`, `region_name`, `statistic_code`, `aep`, `equation`, `reason`).
  Unresolved regions are listed there too, with an empty `equation` and the reason.
- Advisory text (exact):
  `NSS SNAPSHOT (not report-verified): the {code} equations were harvested from the USGS
  NSS web service on {date} and parsed by flowfreq. They have not been checked against the
  source report and carry no covariance or model-error variance; intervals use NSS's
  average standard error of prediction where NSS reports one. {n_unparsed} published
  equations or regions could not be parsed and are listed in StateLibrary.unparsed.`
- Each equation is `nss.to_regression_equation(...)` with `region_code` = the NSS code
  (`GC2015`), `region_name` = the NSS region name, `citation` =
  `Citation(publication=f"{author} {title}", table=f"NSS region {region_id} ({region_code}),
  citation {citation_id}", year=<first four digits of title, else None>, url=citation_url)`,
  limits and units from the region's parameters. Then `sep_log` is set from `asep_pct` where
  present and positive: `sqrt(ln(1 + (asep/100)^2)) / ln(10)`, the same conversion WY.json
  documents for a printed SE_P in percent. NSS's raw `nss_sep` is not stored (it is not a
  log10 SEP; `data/nss_snapshots/` keeps it).
- `source_reports`: one entry per NSS citation:
  `{"citation": f"{author} {title}", "url": citation_url, "nss": {"citation_id": id,
  "last_year_of_data": y}}`. `issue`: null. `notes`: a generated sentence naming the
  snapshot file, `retrieved_utc` and the counts.
- Generator: `flowfreq/regression/nss_library.py::snapshot_to_library(snapshot: dict) ->
  Optional[dict]` (None when the snapshot has no regions), driven by
  `tools/build_nss_snapshot_library.py --snapshot <file>... --out <dir>`. The tool skips a
  jurisdiction with a published file and says so.
- Coverage: `tools/gen_regression_coverage.py` reports `nss snapshot (N eqs, M unparsed)` in
  the Equations column for these states, and `tests/test_catalog.py`'s staleness check keeps
  `docs/REGRESSION_COVERAGE.md` current.

### B. `estimate_at_site`

New module `flowfreq/ungaged.py`. One call that decides and documents:

```
estimate_at_site(lat, lon, *, state, aeps=None, at_site=None,
                 allow_nss_snapshot=False, allow_provisional=False,
                 allow_extrapolation=False, weighting_space="log",
                 cache=None, timeout=60.0) -> SiteEstimate
```

- Delineate with `streamstats.delineate_and_get_characteristics(state, lat, lon, cache=,
  timeout=)`; locate regions with `streamstats.locate_regression_regions(state,
  ws.polygon_geojson, timeout=)`; load the library with `load_state(state,
  allow_nss_snapshot=, allow_provisional=)`.
- Map each located NSS region to a library region with a new
  `StateLibrary.region_for_nss_code(code) -> Optional[str]`: equality with `region_code`
  first, then a `region_name` containing the NSS code (verified files write
  `"... (NSS GC1239)"`). A located region with no library match is recorded in
  `SiteEstimate.unmatched_regions` and the estimate proceeds on the rest, renormalised; if
  none match, `EquationsUnavailable`.
- Per AEP (default: the AEPs every matched region publishes, largest first): evaluate each
  region's equation on the characteristic values (`regression.evaluate(eq, inputs,
  allow_extrapolation=)`); one region gives its own interval; several regions give
  `evaluate_weighted(..., space=weighting_space)` with no interval and
  `interval_method="none: area-weighted across N regions"`.
- `at_site`: a `FrequencyResults` the caller fitted (the app already has one). When given
  and the basin lies in one region, the result also carries
  `gage_weighting.weight_with_regression(at_site, estimates)`; `method` becomes
  `"at-site weighted with regression"`, otherwise `"regression"`.
- `SiteEstimate` fields: `lat`, `lon`, `state`, `method`, `quantiles` (DataFrame: `aep`,
  `flow_cfs`, `lower`, `upper`, `interval_level`, `interval_method`, `extrapolated`,
  `regions`), `regions` (the `RegressionRegionLocation`s used), `unmatched_regions`,
  `characteristics` (code to value), `library_status`, `advisory`, `estimates` (AEP to list
  of `(fraction, RegressionEstimate)`), `weighted` (DataFrame or None), `provenance`
  (dict: snap, delineation, library file status, out-of-range codes, timestamps),
  `to_markdown()`.
- No transposition branch (roadmap §5 is out of scope for the deadline). No NSS online
  fallback: a jurisdiction with no library file raises `EquationsUnavailable` whose message
  says whether NSS serves the state at all (DC and four territories do not).
- CLI: `flowfreq ungaged --lat --lon --state [--allow-nss-snapshot] [--aep ...]` prints
  `to_markdown()`.

### C. Regional skew rows for every jurisdiction

`flowfreq/data/regional_skew.csv` gets one `pending` row per jurisdiction that has none (46
rows), each citing the state's current USGS peak-flow report (or, for a jurisdiction with no
B-GLS study, the Bulletin 17B Plate I map or the report that says which map it uses), with
the printed statewide or regional values in `notes`. The row stays `pending` until a second
agent re-reads the cited source and confirms every printed value; rows whose values are
confirmed are flipped to `verified` with `skew`, `skew_mse`, `effective_record_years` and
`table` filled. A row with a citation that cannot be opened is `pending` with that noted.
`load_table()` validation and `tests/test_catalog.py`'s coverage check gate the data.

### D. The app

1. New branch from `main` with the three commits of `bump-flowfreq-v0.10.0` cherry-picked,
   pinned to `v0.15.0`, drift fixed, `make check` green. PR for the user to merge.
2. After 0.16.0 is released with A and B: a sidebar mode `Gaged site` / `Ungaged site`.
   `ungaged_ui.py` at the repo root renders the ungaged page: lat, lon, jurisdiction select
   (56, from `flowfreq.regression.jurisdictions.JURISDICTIONS`), "Allow NSS snapshot
   equations (not report-verified)" checkbox, AEP multiselect; output: the quantile table,
   the advisory as a warning banner, a coverage badge (`verified` / `nss snapshot` / `none`
   from `available_states()` and `nss_snapshot_states()`), the regions and characteristics
   used, and a Markdown download of `to_markdown()`. Errors (`UnsnappablePointError`,
   `UnsupportedRegionError`, `EquationsUnavailable`, transport errors) are shown, never
   raised. Pin bumped to `v0.16.0`.

### E. QA and release

- `tools/nationwide_smoke.py`: for each jurisdiction, the longest-record unregulated catalog
  gage; run `run_ffa(station_skew_only=True, site_no=...)` and
  `estimate_at_site(lat, lon, state=..., allow_nss_snapshot=True)` at its location; write
  `docs/NATIONWIDE_SMOKE.md` with pass/fail, the 1 % AEP at-site and regression values, and
  the failure reason. Network-bound, run by hand; its table is the deadline's evidence.
- Docs: README coverage section, `docs/MASTER_ROADMAP.md` §3.2 and §4 status, `TODO.md`
  entry, CHANGELOG, `docs/vignettes/ungaged_nationwide.py`.
- Releases: 0.16.0 after A, B and C merge; a patch or 0.17.0 after the smoke fixes; app pin
  bumps follow each.

## Out of scope until after 2026-11-15

Region-indicator interaction equations (GA/NC/SC, 24 equations); transposition inside
`estimate_at_site`; an NSS online backend behind the same interface; region polygons for the
41 new states (region selection is by StreamStats polygon, which needs the network);
`verified` transcription of any new state; the methow-flow tier 2 #4/#5 and tier 3 items.

## Risks

- StreamStats does not delineate everywhere (Alaska partially, no territories). The smoke
  run documents which jurisdictions the ungaged path reaches.
- NSS `asep_pct` as `sep_log` is a judgement; the advisory says so and the snapshot keeps the
  raw field.
- Hallucinated skew citations. Mitigated by the two-agent transcribe-then-confirm rule and
  the requirement that every citation URL was opened.
