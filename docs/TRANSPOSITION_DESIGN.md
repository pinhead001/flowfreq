# Design: transposing computed flows to an ungaged site

Moving a computed streamflow statistic from a gaged donor basin to a nearby
ungaged target basin by drainage-area ratio, with the exponent taken from the
applicable published USGS regional regression.

**Status: the flood series is implemented** in `flowfreq/transpose.py`
(`RegressionExponents`, `transpose_frequency`, `TransposedResults`), tested in
`tests/test_transpose.py`. Sections 6 and 7 — flow-duration and low-flow
transposition — remain design only.

---

## 1. Why

The motivating case is the same one behind the rest of this library: a
CLOMR/LOMR or a permit application needs a discharge at a site with no gage.
There is a gage 4 miles downstream on the same stream. The defensible move is
to transpose the gage's fitted statistics to the site of interest, state the
method, and show the arithmetic.

Today a user does that by hand, in a spreadsheet, with an exponent typed in
from a PDF. That is exactly the kind of step that produces a number nobody can
reproduce two years later when a reviewer asks. The library already computes
the donor's statistics with full provenance; transposition is the last
undocumented hop between that and the number that goes in the submittal.

This is deliberately **not** an attempt to reimplement StreamStats or NSS. It
does one thing those do not: it transposes *the user's own at-site fit*, made
from their own record with their own perception thresholds and PILF choices,
rather than a regression estimate computed from basin characteristics alone.

---

## 2. What already exists

| Piece | Where | State |
|---|---|---|
| Donor flood statistics | `FrequencyResults.quantiles` / `.confidence_limits` | Works |
| Donor low-flow statistics | `LowFlowResults`, `flowfreq.lowflow` | Works, incl. `p_zero`/`n_zero_years` |
| Donor drainage area | `usgs.GageAttributes.drainage_area_sqmi` | Works |
| Baseflow index, flashiness | `flowfreq.regime` | Works — makes a donor-similarity screen feasible |
| Flow-duration statistics | `Hydrograph.plot_flow_duration_curve` | **Plot side effect only.** Returns a table at fixed percentiles (1,5,10,20,50,80,90,95,99); there is no standalone FDC function to build on |
| Regional regression coefficients | — | **Nothing.** No lookup, no table, no schema |

Two gaps, then: the exponents have nowhere to come from, and flow-duration
statistics can only be obtained by drawing a figure.

---

## 3. The method, and where the risk is

For a statistic `Q` at exceedance probability `p`:

```
Q_target(p) = Q_donor(p) * (A_target / A_donor) ** b(p)
```

`b(p)` is the drainage-area exponent from the applicable regional regression —
that is, the exponent on drainage area in the published equation
`Q_p = a * A^b * (other basin characteristics)`. `b = 1.0` is the naive
area-ratio and is almost never what the published regression says.

### The failure mode

A wrong exponent does not raise. It produces a plausible discharge — right
order of magnitude, monotone across return periods, passes every internal
consistency check — that is simply wrong by 10–30%. Inside a LOMR that is
worse than refusing to compute it, because it manufactures confidence in a
number nobody will question.

Three ways to get it wrong, all silent:

1. Using an exponent from the wrong region or the wrong (superseded) report.
2. Applying a flood exponent to a low-flow or flow-duration statistic
   (section 7 — this is a category error, not an approximation).
3. Extrapolating `b(p)` past the range the regression was published over,
   which is **guaranteed for every quantile more frequent than Q2** —
   six of the fourteen in the standard set (section 4).

The module cannot detect (1). It must therefore *record* the citation as
mandatory provenance and refuse to run without one, rather than accept a bare
float and lose the only evidence of where the number came from.

### Applicability limits

The drainage-area ratio method is defensible over a limited range. Commonly
cited USGS guidance is an area ratio between roughly 0.5 and 1.5, on the same
stream, within the same hydrologic region. Outside that band the published
regression should be used directly instead.

Recommendation: **raise by default** outside the band, with an explicit
`allow_out_of_range=True` override that records the violation in the result's
provenance. Precedent: `compare_engines` raises rather than falling back
(`FORTRAN_ENGINE_DESIGN.md` §9 q4), for the same reason — a report that
sometimes means one thing and sometimes another is a report a reviewer cannot
take at face value.

---

## 4. The exponent set, and the frequent end

Published flood regressions give `b` at a fixed handful of AEPs — typically
50%, 20%, 10%, 4%, 2%, 1%, 0.5%, 0.2% (Q2 through Q500).
`Bulletin17C.STANDARD_AEP` runs from **0.995 to 0.002** — fourteen points,
Q1.005 through Q500. The whole set transposes: `transpose_frequency` scales
every quantile the donor carries, and the six more frequent than Q2 are the
ones the published regression has nothing to say about.

So `b(p)` has to be interpolated between published points and, below AEP 0.50,
extrapolated beyond them.

### Interpolating

Interpolate `b` against the standard normal deviate `z = Φ⁻¹(1 - AEP)`, not
against AEP directly. That is the abscissa the LP3 curve itself is nearly
linear in, published points are roughly evenly spaced along it, and it treats
both tails symmetrically. `log10(T)` is an acceptable alternative and gives
very similar answers over the published range.

`b(p)` typically varies smoothly and monotonically across the published range,
often decreasing as recurrence interval increases — but **the direction is
region-specific and the module must not assume one**. Use a shape-preserving
interpolant (PCHIP) or plain linear interpolation; do not fit a spline that
can overshoot between published points.

### Extrapolating to the frequent end

This is the load-bearing decision, and it is worth being precise about how far
the extrapolation actually reaches. In normal-deviate space the published set
spans `z = 0` (AEP 0.50) to `z = 2.88` (AEP 0.002). AEP 0.995 sits at
`z = -2.58`. **The frequent end is as far below the published range as Q500 is
above it** — it is a long extrapolation, not a short one, and the six
quantiles more frequent than Q2 all live there.

Two defenses, both implemented:

- Default `extrapolation="clamp"`: hold the endpoint exponent. `"linear"`
  extends the slope through the two nearest published points for a caller who
  has a reason; `"error"` refuses outright. Clamp is the default because
  extending a trend 2.6 normal deviates past the last observation is a curve
  fit, not an estimate.
- A hard `b_bounds` clamp (default `(0.5, 1.0)`) applied after either mode,
  with the affected AEPs flagged separately. This is a backstop against a
  pathological exponent set, not a routine correction: on a normally-sloped
  published set, even full linear extrapolation to AEP 0.995 only reaches
  about 0.89, so the bound does not bite. `test_linear_extrapolation_to_0995
  _stays_under_the_default_bound` measures that rather than assuming it.

Every extrapolated AEP is recorded in the provenance, logged once, and named
in `to_markdown()`'s caveats, so a report says which quantiles rest on a
published exponent and which on an assumption of the analysis.

**The honest caveat, which flagging does not remove:** a flood regression's
exponent describes flood response. Extrapolated to AEP 0.9 or 0.995 it is
being asked about flows that are not floods, where baseflow and channel
geometry control the scaling and the regression has no observations at all.
Clamping keeps the arithmetic sane and the flagging keeps it honest, but the
better answer, where the data exist, is to **supply exponents covering the
frequent end** — from a flow-duration regression (section 6) — in the same
`RegressionExponents`. Doing that removes the extrapolation entirely, and the
API deliberately makes it a matter of passing more points rather than calling
a different function.

---

## 5. Proposed API

Three layers, each usable without the ones above it.

```python
from flowfreq.transpose import RegressionExponents, transpose_frequency

exponents = RegressionExponents(
    aeps=[0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002],
    exponents=[0.83, 0.81, 0.80, 0.78, 0.77, 0.76, 0.75, 0.74],
    citation="USGS SIR 2015-5164 table 6, Tennessee hydrologic region 3",
    region="TN-3",
    valid_area_range_sqmi=(1.0, 3000.0),
)

transposed = transpose_frequency(
    results,                      # the donor's own FrequencyResults
    donor_area_sqmi=178.0,
    target_area_sqmi=131.0,
    exponents=exponents,
)

transposed.quantiles            # same shape as FrequencyResults.quantiles
transposed.provenance.area_ratio        # 0.736
transposed.provenance.extrapolated_aeps # [0.995 ... 0.67] -- flagged
transposed.to_markdown()        # the arithmetic, per quantile, for an appendix
```

`RegressionExponents` requires `citation` — not optional, no default. A
transposition whose exponent has no source is not a result worth having.

`TransposedResults` wraps rather than impersonates `FrequencyResults`. It
carries the transposed `quantiles`/`confidence_limits` plus a `provenance`
record, and deliberately does **not** expose `mean_log`/`std_log`/`skew_*`:
those describe a fit to the donor's record, and there is no fitted LP3 at the
target. Synthesizing them would be inventing a fit that was never performed —
the same rule the Fortran adapter already follows
(`FORTRAN_ENGINE_DESIGN.md` §4: never synthesize a field the source did not
report).

### Confidence limits

Scale the donor's limits by the same factor and mark them transposed. This is
first-order and **understates** the true uncertainty, because it treats the
transposition itself as exact. Say so in the field name and in the markdown
output. If the caller supplies the regression's standard error of prediction,
a combined interval becomes possible — a Phase 2 item, not a first pass, and
never a fabricated default.

---

## 6. Exceedance (flow-duration) transposition

Different problem, and the right answer depends on what the caller wants.

### If they want duration *statistics* (Q10%, Q50%, Q95%)

Same area-ratio machinery, different exponent set: `b` varies substantially
across the duration curve — near-linear at the wet end (high flows are
largely a function of contributing area) and dropping away, sometimes
erratically, at the dry end where geology and storage dominate rather than
area. A single exponent applied across the whole FDC is wrong at one end or
the other.

So: `transpose_duration()`, taking a `RegressionExponents` indexed by
*exceedance* probability from a published duration regression, with the same
interpolation, guardrails and provenance.

Prerequisite: extract a real FDC function out of
`Hydrograph.plot_flow_duration_curve`. Duration statistics should not be a
by-product of drawing a figure, and transposition needs them at arbitrary
percentiles, not the nine hardcoded ones.

### If they want a daily *time series* at the target

Then area-ratio scaling of the whole hydrograph is the wrong tool, and the
recommendation is the **QPPQ / probability-preserving transfer**: for each
day, map the donor's flow to its nonexceedance probability on the donor FDC,
then map that probability back through the *target's* estimated FDC. Timing
comes from the donor; magnitude comes from the target's own duration curve.
It degrades gracefully where a constant ratio does not, because it never
assumes the two basins scale by the same factor at every flow level.

This is well established in the literature; the specific citations should be
pinned down and verified before any of it appears in a report, rather than
carried over from this document on trust.

Scope call: QPPQ is a larger piece of work than the statistic-level
transposition and depends on having a target FDC from somewhere. Recommend
**deferring it** to a second pass, and shipping the duration-statistic
transposition first.

### Zero flows

If the donor's Q95 (or any transposed duration statistic) is zero, the ratio
method is degenerate: `0 * anything` is `0`, which asserts a fact about the
target that the donor's record does not support. Return `None` with a stated
reason, never a zero. `LowFlowResults` already treats zero-flow years as a
first-class concept; this must be consistent with that.

---

## 7. Low-flow transposition

**Do not transpose low flows with a flood exponent.** This is the single most
important sentence in this document. 7Q10 is controlled by baseflow storage —
surficial geology, aquifer transmissivity, soil permeability — not by
drainage area. Two adjacent basins of identical area can differ severalfold
in 7Q10, and one of them can be zero. That is why published low-flow
regressions carry a baseflow or geology term and flood regressions do not.

Recommended approach, in descending order of defensibility:

1. **Use the regional low-flow regression directly** at the target, with the
   target's own basin characteristics. Transposition is not the best tool
   here; it is the fallback when the target's characteristics are unavailable.
2. **Area-ratio with a low-flow-specific exponent** from a published low-flow
   regression, gated behind a donor-similarity screen (below).
3. **Baseflow-normalized ratio** — scale by the ratio of baseflow yields
   rather than raw areas, when both basins have enough daily record to compute
   a BFI. `regime.baseflow_index()` already exists and takes
   `drainage_area_sqmi`. Note `regime`'s own warning: a BFI is only comparable
   to another BFI computed by the same method with the same parameters, so
   the screen must compute both sides itself rather than accept a supplied
   number.

### Guardrails, stricter than for floods

- Tighter area-ratio band (recommend 0.7–1.3 default, versus 0.5–1.5).
- Refuse when the donor statistic is zero or when `p_zero` is nonzero above a
  threshold — a donor that goes dry cannot tell you what a target does.
- Require an explicit caller assertion that donor and target share a
  physiographic/hydrogeologic setting, recorded in provenance. The module
  cannot verify this and should not pretend to.
- Warn when the donor's record does not span the target's drought of record,
  or is short enough that the low-flow tail is poorly determined.

### Separate function, not a flag

`transpose_low_flow()`, distinct from `transpose_frequency()`. The exponent
source differs, the guardrails differ, the failure modes differ, and the zero
handling differs. Folding them into one function behind a `kind="low"` flag
would make the strict path reachable by accident from the lenient one.

---

## 8. Correctness plan

| Risk | Test |
|---|---|
| Exponent interpolation wrong | Reproduce published `b` values exactly at the published AEPs (interpolation must be an identity there) |
| Extrapolation silently wrong | Assert every AEP above 0.50 is flagged in `extrapolated_aeps`, that `clamp` holds the endpoint, and that `error` raises |
| Area-ratio band not enforced | Ratios of 0.4 and 1.6 raise; 0.5/1.5 pass; override records the violation |
| Arithmetic wrong | Hand-computed transposition on a fixture (Big Sandy, halved area) checked against the closed form per quantile |
| Monotonicity lost | Transposed quantiles must remain monotone in AEP for any valid exponent set |
| Ratio of 1.0 is not the identity | `A_target == A_donor` must return the donor's own numbers bitwise, for every exponent |
| Zero low flow mishandled | Donor 7Q10 of 0 returns `None` with a reason, never 0.0 |
| Provenance lost | A `RegressionExponents` without a citation fails to construct |

The identity test (ratio 1.0) is cheap and catches an entire class of
exponent-plumbing bug that no eyeball review would.

---

## 9. Estimate

| Piece | Estimate | State |
|---|---|---|
| `RegressionExponents` + interpolation/extrapolation + provenance | 0.5 day — the risk lives in extrapolation policy | **Done** |
| `transpose_frequency` + guardrails + markdown output | 0.5 day | **Done** |
| Standalone FDC extraction out of `hydrograph.py` | 2 h | Not started |
| `transpose_duration` | 2 h once the FDC function exists | Not started |
| `transpose_low_flow` + similarity screen | 0.5 day | Not started |
| QPPQ daily-series transfer, if pursued later | ~1 week | Deferred |

---

## 10. Open questions

1. **Where do exponents live?** Caller-supplied is the only option with no
   maintenance liability, and is assumed above. A small bundled table for a
   few well-documented regions would be convenient and would immediately go
   stale — every state revises its regressions, and a stale exponent is the
   silent failure in section 3. Recommend caller-supplied only for the first
   pass; revisit if the same three regions keep getting retyped.
2. **Should `transpose_frequency` accept a bare float exponent?** Convenient
   for `b = 1.0` sensitivity checks; also the easiest way to lose the
   citation. Suggest allowing it only via an explicit
   `RegressionExponents.constant(b, citation=...)` constructor, so the
   citation requirement survives.
3. **Weighted donor/regression estimate.** USGS practice for a site on a
   gaged stream weights the transposed gage estimate against the regression
   estimate at the ungaged site, with the weight depending on how far the
   area ratio departs from 1. That needs the *full* regression equation, not
   just its area exponent, and a decision about whose basin characteristics
   to use. Deliberately out of scope for the first pass; note it so the API
   does not foreclose it.
4. **Does the target need its own AEP set?** Transposing preserves AEPs by
   construction, so no — but a caller who wants an AEP at the target and did
   not compute it at the donor cannot get it by transposition. Worth an
   explicit error rather than a quiet empty row.
