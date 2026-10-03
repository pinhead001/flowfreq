# Requirements: sub-daily timing and ramping-rate metrics for flowfreq

**Status:** implemented -- `flowfreq/subdaily.py`, branch
`claude/flowfreq-hour-max-ramping-rate-z0omag`. This document is the design
that was built, not a forward-looking proposal; where the implementation
departed from a first draft of it, the section says so and why.

**Date:** 2026-09-12
**Context:** written from the Methow sub-basin work. Two things the existing
`diel_variation` cannot answer are being asked of the instantaneous record:
*when* in the day the flow peaks, and *how fast* it changes. `diel_variation`
reports a day's amplitude and CV, which say how much flow moved but nothing
about the shape or timing of the movement -- and those are the two properties
that separate a snowmelt signal from a hydropeaking one, and the two that
salmon-habitat work actually consumes (redd dewatering and fry stranding are
driven by down-ramp rate, not by daily range).

---

## S1. Scope

In scope:

1. **Hour of daily maximum and hour of daily minimum**, per local calendar
   day, plus a period-of-record summary using **circular** statistics.
2. **Ramping rate** -- the rate of change between consecutive instantaneous
   observations -- per local calendar day, reported as up-ramp, down-ramp and
   mean absolute rate, in three unit systems: **cfs/hr**, **%/hr**, and
   **ft/hr** (stage).
3. **Gage-height (NWIS parameter 00065) retrieval**, without which ft/hr has
   no data source.

Out of scope, deliberately:

- **Any natural-versus-regulated classification.** This follows
  `diel_variation`'s existing stance for the same reason: an afternoon
  snowmelt peak and an afternoon load-following peak are not distinguishable
  from the hydrograph alone. These metrics are descriptive inputs to that
  judgement, not the judgement.
- **Frequency analysis of ramping rates.** An annual-maximum-ramp series is a
  plausible future thing to fit, but the ramping literature reports
  exceedance counts against a stated licence limit, not return periods, and a
  down-ramp series is not plausibly LP3. Not attempted.
- **Sub-daily flow transposition.** `flowfreq.transpose` and
  `flowfreq.qppq` work on daily and duration statistics; nothing here feeds
  them.

## S2. Where the code lives, and why it is a new module

`diel_variation` and `diel_variation_summary` are today in
`flowfreq/regime.py`, whose module docstring already stretches to accommodate
them ("plus ... metrics computed from an instantaneous series"). The obvious
move would be to add the new functions beside them.

**They go in a new `flowfreq/subdaily.py` instead, and the two diel functions
move with them.** The deciding reason is `pyproject.toml`'s mypy override
list: `flowfreq.regime` is on it, exempted as pre-existing debt, and that
file's own comment states the rule -- *"Nothing new gets added to this list --
a new module is checked by default."* Six new public functions added to
`regime.py` would be six functions the type checker never looks at, in a
package that ships `py.typed` and therefore promises downstream checkers that
its annotations mean something. A new module is checked from its first commit.
That is the whole argument; it also happens to keep `regime.py` from passing
1,600 lines.

The move is source-compatible and is meant to stay that way:
`flowfreq.regime` re-exports `diel_variation` and `diel_variation_summary`
from `flowfreq.subdaily`, so `from flowfreq.regime import diel_variation`
keeps working, as does `from flowfreq import diel_variation`. The existing
`tests/test_regime.py` imports them from `regime` and was **not** edited, which
is the back-compatibility test: if the shim broke, that file would fail.

Moving them has one consequence worth stating plainly: `diel_variation`'s
body is now mypy-checked where it previously was not. Two annotations were
tightened to satisfy that. No behaviour changed, and
`tests/test_regime.py::TestDielVariation` passes unedited.

## S3. Hour of daily maximum and minimum

### S3.1 What is computed

`daily_extreme_timing(iv_data, tz, ...)` returns one row per local calendar
day with:

| column | meaning |
| --- | --- |
| `date` | local calendar date |
| `hour_of_max`, `hour_of_min` | fractional local hour in `[0, 24)` of the day's extreme |
| `max_flow_cfs`, `min_flow_cfs` | the extreme values themselves |
| `n_tied_max`, `n_tied_min` | how many observations equal that extreme |
| `n_obs`, `expected_obs`, `complete` | coverage, as in `diel_variation` |

Fractional, not integer, hour: `14.25` is 14:15 local. Rounding a 15-minute
record to the hour throws away three quarters of the timing resolution the
record has, and the circular concentration in S3.3 is sensitive to it.

### S3.2 Two cases where a timing number would be a fabrication

**A flat day has no timing.** If the day's every valid observation is the
same value, `max_flow_cfs == min_flow_cfs` and the "hour of the maximum" is
whichever observation happened to be first. `hour_of_max` and `hour_of_min`
are **NaN** for such a day, not `0.0`. This matters more than it sounds: a
regulated reach held at a constant release produces runs of exactly these
days, and a `0.0` would enter the circular mean in S3.3 as a real
midnight peak and drag it.

**Ties are reported, not hidden.** When the extreme is not unique but the day
is not flat (two equal afternoon peaks, or a sensor resolution coarse enough
to report the same value twice at the crest), the *first* occurrence is
reported and `n_tied_max` says how many there were. A reader can then see
that "hour_of_max = 15.5" was one of four equally high readings rather than a
single identified moment. Returning the first occurrence silently, with no
column saying so, is the failure mode this avoids.

### S3.3 The summary must be circular, and this is the load-bearing part

`extreme_timing_summary(daily, which="max")` pools a timing column across
days. **The arithmetic mean of an hour column is wrong**, and wrong in a way
that looks plausible on paper. A site that peaks reliably just before and
just after midnight -- hours 23.5 and 0.5 -- has an arithmetic mean hour of
12.0: noon, the one time of day it never peaks. The error is largest for
exactly the sites whose timing is most consistent.

So the summary uses the standard circular (directional) statistics, mapping
hour to angle by `theta = 2*pi*h/24`:

- `C = mean(cos theta)`, `S = mean(sin theta)`
- **mean hour** `= atan2(S, C) mod 2pi`, back-converted to hours
- **`concentration`** `r = sqrt(C^2 + S^2)`, in `[0, 1]`: the resultant
  length. `r = 1` is identical timing every day; `r = 0` is timing spread
  evenly round the clock. This is the number that says whether the mean hour
  means anything, and it is reported next to it for that reason.
- **`circular_std_hours`** `= sqrt(-2 ln r)` in radians, converted to hours --
  the usual circular standard deviation. Infinite as `r -> 0`, which is
  correct: scattered timing has no spread in hours to quote.
- **`rayleigh_p`** -- p-value of the Rayleigh test for uniformity, using
  Zar's approximation
  `p = exp(sqrt(1 + 4n + 4(n^2 - R^2)) - (1 + 2n))`, `R = n*r`.
  A large `p` says the data are consistent with no preferred hour at all, in
  which case the mean hour should not be reported as a diel peak time.

`circular_hour_statistics(hours)` is public and separately testable, because
the near-midnight case above is precisely the thing a future edit might
"simplify" back into a `.mean()`.

Only rows marked `complete` are pooled, matching `diel_variation_summary`.
`which="min"` pools `hour_of_min` on the same machinery; the min/max pair is
the point -- a morning minimum with a late-afternoon maximum is the snowmelt
signature, and either half alone is ambiguous.

### S3.4 Daylight saving

Grouping is on local calendar days and the hour is read from local wall-clock
time, so on a fall-back day two observations genuinely share a fractional
hour and the day is 25 hours long; on a spring-forward day the 02:00-03:00
hour does not exist and the day is 23 hours. Both are real properties of
local time, not errors, and hour-of-peak in local time is what a diel
question is asking about.

One thing does follow from it. `expected_obs` here is computed from the
**actual local-day length** (`day_hours * 60 / step_minutes`), not from a
fixed 1440 minutes. `diel_variation` originally used the fixed 1440 and so
marked the two transition days of every year as `complete=False` (spring) or
as more complete than a full day (autumn) for no real reason. That was left
alone at first because fixing it moves numbers already reported from that
function; it was fixed in v0.10.0, and all three per-day functions now share
the actual-length convention. Only `expected_obs` and `complete` on the two
transition days changed -- range, CV and `n_obs` never depended on it.

## S4. Ramping rate

### S4.1 Definition

For each consecutive pair of valid observations,

```
dt_hours = (t_i - t_{i-1}) in hours       # measured on the UTC index
rate     = (v_i - v_{i-1}) / dt_hours
```

**`dt` is measured on the UTC index, never on local wall-clock time.** This
is not a stylistic choice. Across a fall-back transition the local clock goes
01:59 -> 01:00, so a local-time difference is *negative* for a real forward
hour, and a 15-minute interval spanning the transition would compute as
`-45 minutes`: a modest real ramp becomes a large ramp of the wrong sign.
`download_instantaneous_flow` already keeps the index on a single UTC axis for
this family of reasons; this uses it.

### S4.2 Three unit systems

- **`cfs/hr`** -- absolute. Physically comparable within a site; scales with
  basin size, so not comparable between a headwater and a mainstem gage.
- **`%/hr`** -- `100 * (v_i - v_{i-1}) / (v_{i-1} * dt_hours)`, normalized by
  the value the ramp started from. Comparable across sites of different size,
  and the form most hydropeaking limits are written in. **NaN when the
  preceding value is not strictly positive** -- a percentage change from zero
  is not a number, and an intermittent reach produces these.
- **`ft/hr`** -- the same calculation on a gage-height series. This is the
  form most FERC and HCP ramping-rate conditions are actually written in, and
  the quantity that drives stranding: fish respond to the water's edge
  retreating, which is a stage rate.

Column names carry the unit (`max_down_ramp_cfs_per_hr`,
`max_down_ramp_ft_per_hr`), inferred from `value_col` via a small known
mapping and overridable with `units=`. A rate column whose unit is only
recorded in a docstring gets misread eventually.

### S4.3 A guard: %/hr on stage is meaningless

Gage height is measured from an **arbitrary local datum**. A gage reading 4.0
ft and rising to 4.2 ft has not "risen 5%" in any physical sense -- the same
water-surface change at a gage whose datum sits two feet lower reads as a
2.9% rise, and the true depth of water is not what either number divides by.
A percent-per-hour figure computed on stage is a number with no referent that
will nonetheless look like a regulatory ramping rate.

So `ramping_rates` **raises `ValueError`** when relative rates are requested
for a non-`cfs` series, naming the datum as the reason. `relative` defaults
to True for a true-zero quantity (discharge) and False otherwise, so the
common calls need no argument and the mistake is not available. This follows
the precedent of `transpose`'s `_require_kind` and `assemble_target_fdc`'s
`index_gage=`: where two quantities share a plausible-looking shape and only
one is correct, the library refuses rather than returning the plausible
wrong number.

### S4.4 Gaps are not ramps

An NWIS unit-value record has outages. A `dQ` computed across a six-hour hole
in a 15-minute record is not a ramping rate -- it is the net change over six
hours, which will read as a gentle ramp however violently flow actually moved
inside the gap, and can equally hide a large real ramp. Averaged in, these
systematically bias the mean absolute rate **downward**.

Intervals longer than `max_gap_hours` are therefore excluded from every rate
and counted in a `n_intervals_gapped` column, so the exclusion is visible
rather than silent. The default is `3 x` the record's own median time step,
inferred the way `diel_variation` infers it -- a multiple rather than a fixed
hour count, so the default is sane for a 15-minute record and for an hourly
one without the caller having to know which they have.

### S4.5 Which day an interval belongs to

An interval spans two instants and can straddle local midnight. It is
assigned to the local day of its **ending** timestamp by default
(`interval_label="end"`), so an overnight recession lands on the day whose
low flow it produced -- which is the day a stranding question is about.
`interval_label="start"` is accepted for the other convention. The parameter
exists because the choice is real and shifts every midnight-straddling
interval by a day; leaving it implicit would make two runs disagree with no
visible cause.

### S4.6 Per-day columns

| column | meaning |
| --- | --- |
| `max_up_ramp_<u>_per_hr` | largest positive rate; `0.0` if intervals exist but none rose |
| `max_down_ramp_<u>_per_hr` | **signed** minimum, so always `<= 0` |
| `mean_abs_ramp_<u>_per_hr` | mean of `abs(rate)` |
| `n_reversals` | sign changes in the rate sequence -- the standard hydropeaking count |
| `n_intervals`, `n_intervals_gapped` | valid intervals used, and intervals dropped per S4.4 |
| `n_obs`, `expected_obs`, `complete` | as in S3.1 |

Down-ramp is kept **signed**. "A maximum down-ramp of 50" is ambiguous about
sign in exactly the context where sign matters, and a magnitude that compares
as `>` against a limit expressed as a negative number is a bug waiting to
happen. A value that is always `<= 0` cannot be misread.

`0.0` versus NaN for `max_up_ramp` is a real distinction and is kept: `0.0`
means "intervals were measured and the largest rise was zero", NaN means "no
valid interval was available at all". Collapsing them would make a
fully-gapped day indistinguishable from a monotonically falling one.

### S4.7 Summary and exceedance

`ramping_rate_summary(daily, limits=None)` pools complete days: mean and
most extreme up- and down-ramp, mean absolute rate, mean reversals per day.
Passing `limits={"max_down_ramp_cfs_per_hr": -100.0}` adds, per named column,
the count and fraction of complete days exceeding that limit -- exceedance
against a stated licence condition being what ramping numbers are actually
used for. Comparison direction is taken from the sign of the limit, so a
negative limit counts days at or below it and a positive limit counts days at
or above it.

## S5. Gage-height retrieval (NWIS 00065)

`USGSgage.download_instantaneous_stage()` mirrors
`download_instantaneous_flow()` -- same chunking, same UTC axis, same
multi-sensor `ts_id` disambiguation, same `NoInstantaneousDataError` -- and
returns a `gage_height_ft` column in place of `flow_cfs`.

Implementation is a generalization rather than a copy: `_parse_iv_rdb`,
`_resolve_flow_column` (now `_resolve_value_column`) and `_empty_iv_frame`
take a parameter code and an output column name. Duplicating the RDB parser
for a second parameter code would mean the tz-mapping and multi-sensor rules
exist twice and drift.

Two limits belong in the docstring where a caller will hit them:

- **Stage is often shorter than discharge at the same gage, or absent.** A
  site with fifteen years of unit-value discharge may have a partial stage
  record or none. `NoInstantaneousDataError` names the parameter so the
  failure is not mistaken for the site having no unit values at all.
- **Stage is not comparable between gages, and not convertible to depth.**
  Arbitrary datum, per S4.3, and a datum can be reset during a gage
  rebuild, which puts a step in the series that is not a hydrologic event.

## S6. Testing

Per `CLAUDE.md`, every public function gets a happy path and at least one
error case. The tests that carry the actual risk:

- **Circular mean near midnight.** Hours `23.5` and `0.5` must summarize to
  `~0.0`, not `12.0`. The single most important assertion in this work.
- **Uniform timing gives low concentration and a non-significant Rayleigh
  p**, so `r` is doing its job of telling the caller the mean is meaningless.
- **Known-phase synthetic diel series.** Reusing `test_regime.py`'s
  `_pacific_diel_series` shape: a sinusoid with its trough at 04:00 and peak
  at 16:00 local, stored in UTC, must recover `hour_of_min ~ 4` and
  `hour_of_max ~ 16` -- which fails if local-day grouping is dropped.
- **DST fall-back day.** A 25-hour local day must produce a positive `dt` for
  every interval, which is the assertion that fails if `dt` is ever computed
  on local wall-clock time (S4.1).
- **Gapped interval excluded, not averaged in.** A long hole must raise
  `n_intervals_gapped` and leave `mean_abs_ramp` equal to what the
  gap-free intervals alone give.
- **Flat day yields NaN timing**, not `0.0` (S3.2).
- **`relative=True` on a stage series raises** (S4.3).
- **Stage RDB parse** against a captured-shape 00065 payload, including the
  `gage_height_ft` column name and a multi-sensor `ts_id` case.

## S7. Open, and not attempted here

- **Stage-discharge pairing** is built:
  `USGSgage.download_instantaneous_flow_and_stage` / `join_flow_and_stage`, an
  outer join on the UTC instant with NaN where one sensor did not report (no
  interpolation), each parameter keeping its own `qualification_code`.
- **Validation against a published figure** (2026-09-30), for daily extremes
  only. Exelon (2012), *Final Study Report: Downstream Flow Ramping and
  Stranding Study, RSP 3.8*, Conowingo Hydroelectric Project, FERC No. 405,
  SS4.1.1/4.3.1, gives discharges at USGS 01578310 that
  `daily_extreme_timing` reproduces exactly: 141,000 cfs on 2010-10-03, morning
  peaks of 26,100 / 46,200 / 80,000 cfs, and 36,500 cfs on 2010-04-28. One
  figure disagrees by 200 cfs (80,900 reported vs 81,100 on 2010-05-12). See
  `tests/test_subdaily_published.py` for the verbatim quotations. The report's
  one-hour stage declines are a fixed-window metric that `ramping_rates` does
  not compute (TODO.md). **Ramping rates themselves therefore remain verified
  in arithmetic only.**
