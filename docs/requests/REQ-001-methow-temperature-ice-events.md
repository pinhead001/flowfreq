# REQ-001 — from methow-temperature: ice events from flagged days, cold air and gauge qualifiers

- From: methow-temperature (its REQ-014, `D:\a\methow-temperature\docs\requests\REQ-014-flowfreq-ice-events.md`) · Date: 2026-11-21 ·
  Status: **answered** 2026-11-22. Placed here uncommitted; answer in this file under a "Reply" heading.
- Decision behind it (Chris, 2026-11-21): the generic fetch, event and QA/QC code for an icing analysis belongs in
  flowfreq, built by a flowfreq session. methow-temperature then takes flowfreq as a git dependency pinned to the
  merge commit.
- **How you build it is yours to decide** (Chris, 2026-11-21): branch name, cards or native, task split, where the
  code sits in the package, and whether this file is committed in flowfreq.

## What to read

Both files are in `D:\a\methow-temperature`; read them in place, nothing needs copying.

- Plan, **Part A** (Tasks A1 to A4), with "Global Constraints" and "Review Focus" items 2 and 3:
  `docs\superpowers\plans\2026-11-21-icing-38-logger-events.md`. It gives signatures, test code and the steps as
  this repo would have built them. Treat the steps as a proposal.
- Design, sections "Icing event" and "Where the code goes / flowfreq":
  `docs\superpowers\specs\2026-11-21-icing-38-logger-events-design.md`.

## Asked for

1. **`flowfreq.ice`**, counting episodes of flagged days and rating nothing:
   `winter_of`, `count_events`, `cold_air_days`, `cold_spells_per_winter`, `events_in_months` (plan Task A1);
   `winter_event_counts`, `decade_scaled` (A2); `ice_days_from_qualifiers`, `gauge_ice_qaqc` (A3).
2. **An opt-in on the daily download** (`keep_null` on `parse_daily_features` and `download_daily`) that keeps a day
   whose value is null, so an `ICE` day with no flow is not lost. Default behaviour unchanged (A3).
3. **`flowfreq.ghcn.fetch_ghcn_daily`**: NCEI daily minimum and maximum air temperature for one station, in °C,
   NaN where missing, through `waterdata.request` (A4).
4. **A `flowfreq ice-events` command** that writes `air_<station>.csv`, `daily_<site>.csv`, `gauge_ice_events.csv`,
   `gauge_ice_qaqc.csv` and `ice_events_run.json` (A4).
5. A PR that Chris merges to `main`.

## What we rely on

Part B here is written against these, so please keep them or tell us what changed:

- the function names, arguments and returned columns in each task's **Interfaces** block;
- the command's options and its five output files, with the column names given in A4;
- every threshold a required argument or option, with no default in flowfreq (our rule is that settings live in
  this repo's config);
- the event rule: two flagged days join one event when at most `merge_gap_days` unflagged days lie between them;
  days outside the given months are dropped before events are counted; an event belongs to the winter (water year)
  and month of its first day; missing air temperature counts as not cold;
- `ESTIMATED` without `ICE` is not an ice day;
- `decade_scaled` returns NaN, not 0 and not an error, for the scaled values when the air index sums to 0 over the
  site's winters;
- no Methow names in flowfreq code or docs (station numbers may appear as test fixture values).

The call we will make after the merge, for reference:

```
flowfreq ice-events --gauge 12447383 --gauge 12447390 --gauge 12448000 --gauge 12448500 --gauge 12448998
  --gauge 12449500 --gauge 12449950 --air-station USC00459376 --start 2006-10-01 --end 2026-09-30
  --months 11,12,1,2,3 --air-c -10 --lookback-days 2 --merge-gap-days 2 --min-cold-spell-days 5 --out <dir>
```

It runs under pandas 3.1.1, numpy 2.5.4 and pyarrow 26.0.1 on Python 3.12 (flowfreq's lock has pandas 3.0.5).

## Reply wanted

1. The **merge commit sha** on flowfreq `main`. We pin to it; no release tag is asked for.
2. Any name, argument, column or file that differs from the plan's Interfaces.
3. Anything in the event rule you implemented differently, and why.
4. Whether the `ice-events` command was run live against the services, or only tested on recorded responses.

Not blocking us today: Task B1 here (daily minimum, maximum and hour counts) proceeds without flowfreq. Tasks B2
onward wait for the reply.

## Reply (flowfreq, 2026-11-22)

Built on branch `feat/ice-events`, PR pinhead001/flowfreq#410 (Chris merges).

1. **Merge commit sha:** `e953e2fd21fb1ce829cfaa0b282e71ed8ba34043` (PR #410, merged 2026-11-22T01:21Z; all CI green).
2. **Differences from the plan's Interfaces:** none removed or renamed. Additions:
   - `waterdata.request(..., send_api_key=False)` (keyword-only, default True). `fetch_ghcn_daily` uses it, so the USGS
     `X-Api-Key` is never sent to NCEI.
   - `ice_events_run.json` also carries `air_tmin_days` (valid TMIN days) and per gauge `daily_rows` and `ice_days`.
   - `fetch_ghcn_daily` returns every day from start to end (index named `date`), NaN where NCEI has no value.
3. **Event rule, implemented as written,** with these calls of ours:
   - `cold_air_days` looks back by calendar day (a gap in the TMIN index does not stretch the window).
   - `gauge_ice_qaqc.csv` covers only days inside `--months` and days the gauge has a record (no "cold without ice"
     before a gauge's record starts or on days it has no row).
   - The command refuses (non-zero exit, nothing written) an air station with no valid TMIN in the window, `--start`
     after `--end`, and a gauge whose download fails (named).
   - Under "missing air counts as not cold", an ICE day during a gap in the air record shows up as `ice_without_cold`
     in the QA/QC file; read those rows with that in mind.
4. **Live runs:** `fetch_ghcn_daily('USC00459376', '2006-10-01', '2026-09-30')` was run live once: 7305 days in one
   request, 24 s, 7263 valid TMIN. The `ice-events` command as a whole was tested only on recorded/synthetic
   responses, not run live against the Water Data API.
