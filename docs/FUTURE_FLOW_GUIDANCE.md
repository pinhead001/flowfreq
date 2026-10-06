# Future-flow guidance by jurisdiction

`docs/MASTER_ROADMAP.md` §6.3 is the plan this file serves. The national framework
(`flowfreq.future_flow`, §6.3.1) applies everywhere. This file records **state-specific**
guidance, one section per jurisdiction, added as each state's rollout wave reaches it
(§6.3.2).

Rules:
- Every entry cites its source. If a state has no guidance, record "None published as of
  <date>". The national default then applies and is labelled as such.
- Quantitative guidance is encoded as a `ChangeFactorSet` under `flowfreq/data/future/`,
  with a test that reproduces one published example. Qualitative guidance is recorded here
  only; no numbers are made up for it.

Fields to record for each entry: source and citation · legal status (required /
recommended / informational) · method type (multiplier, projected-precipitation
regression, or simulation) · scenario · horizon · applicable geography · encoded factor
set (if any).

---

## National (applies to all jurisdictions)

Status: **reviewed 2026-09-30; no national factor set exists to encode** (#36).
`flowfreq.future_flow.available_factor_sets(level="national")` is therefore empty, and a
caller who wants future-condition quantiles must supply a site-specific `ChangeFactorSet`
derived by one of the methods below. The review is machine-readable in
`flowfreq/data/future/national_sources.json` (`future_flow.national_source_review()`),
with every table checked and why it is not a factor set.

| Source | What it provides | Tabulated national factors? |
|---|---|---|
| FHWA HEC-17, 2nd ed. (Kilgore, Herrmann, Thomas & Thompson, June 2016, FHWA-HIF-16-018) | Method: a five-level analysis framework (Ch. 7, Tables 7.1-7.4) choosing trend analysis, or rainfall-runoff / regression driven by downscaled climate projections, per project | **No.** Tables 4.1 (p. 4-14), 8.1 (p. 8-3), 8.2-8.3 (pp. 8-8 to 8-9) and the MnDOT ranges (p. 8-5) are single-site case studies |
| NCHRP Project 15-61, *Design Practices* (Kilgore et al., March 15, 2019; TRB online, preliminary unedited) and its Final Report | Method: trend projection (Ch. 5), ensemble GCM precipitation ratios for rainfall-runoff (Ch. 6), projected climate variables in USGS regression equations (Ch. 7), flood-index scaling (Ch. 8), continuous simulation (Ch. 9). Every factor is derived per site | **No.** Tables 6.6-6.9 (pp. 54-56, Denver precipitation ratios), 7.1 (p. 80, hypothetical MAP increase) and 8.1 (p. 87, Pajarito Creek NM) are worked examples |
| NOAA Atlas 15 (NWS Office of Water Prediction) | Future precipitation frequency (Volume 2) | **Not published.** CONUS estimates in peer review September 2026, publication planned 2027; a 2024 Montana pilot only. Needs a documented precipitation-to-flow elasticity (roadmap §6.3.1) before it can become a flow factor set |

Fields: legal status informational (both federal documents are guidance, not regulation) ·
method type: per-site derivation (multiplier from projections, projected-precipitation
regression, or simulation) · scenario / horizon / geography: chosen by the user per site.

Revisit when NOAA Atlas 15 Volume 2 is published, or when FHWA or NCHRP publishes a
tabulated national product.

### NCHRP 15-61 per-site procedures (implemented)

The Guide gives no national factors, but three of its per-site procedures are mechanical once
the projection inputs exist. They are implemented in `flowfreq.future_flow`. Each takes
user-supplied projections (nothing is downloaded) and reproduces the Guide's own worked
example in `tests/test_future_flow.py`:

| Function | Guide | What it does | Worked example reproduced |
|---|---|---|---|
| `ensemble_ratio_summary` | Ch. 6, Steps 8-9, Eqs. 6.1-6.3 (pp. 48-50) | Ensemble mean, SD and confidence limits of per-model future/baseline quantile ratios; projected quantile = historical x ratio, with the 0.1-AEP ratio substituted for rarer AEPs (Eq. 6.3) | Tables 6.8 and 6.9 (p. 56) from Tables 6.6 and 6.7 (p. 54), Denver CO, 24-hour precipitation |
| `regression_change_factors` | Ch. 7, Section 7.1 (pp. 77-79) | Evaluates a regional regression equation with historical and projected inputs; refuses any input outside the calibrated range (no override, per Section 7.3) and any equation without a projected variable | Table 7.1 (p. 80), from Lumia et al. (2006), SIR 2006-5112, Table 1, Region 5 |
| `index_flood_projection` | Ch. 8, Section 8.1 (pp. 85-86) | Scales a historical curve by the projecting method's projected/historical index-flood (0.1 AEP) ratio | Table 8.1 (p. 87), Pajarito Creek NM, with the index flood from Eq. 8.1 through `regression_change_factors` |

Two findings from reproducing them:

- **Table 7.1's 0.2-AEP column is internally inconsistent.** Region 5's printed 5-year
  equation gives 5,270 and 6,580 cfs at the example's inputs; the table prints 5,060 and 6,320.
  Their ratio (1.249) is the equation's, so the flows, not the factor, are off. The other five
  columns agree at the table's three significant figures. The test pins this.
- **The confidence limits in Tables 6.8-6.9 are `mean -/+ z * sd / sqrt(n)`.** The Guide says
  only that they assume a normal distribution; this is the convention that reproduces them.

What stays with the engineer, because the Guide makes it a judgement: choosing GCMs, scenarios,
periods and grid cells (Ch. 3, Ch. 6 Steps 2-4); discarding anomalous models or cells (Ch. 6,
Step 8); and the rainfall-runoff model that turns Ch. 6's projected precipitation into flow.
The trend projection of Ch. 5 and the continuous simulation of Ch. 9 are not implemented.

Each function returns a site-level `ChangeFactorSet` (`level="site"`) or a `FutureQuantiles`
built on one. Site sets never take part in national/state precedence.

### Precedence and side-by-side output

`apply_with_precedence(quantiles, national, state_sets, state=...)` applies the governing set.
When an applicable state set overrides the national default, the result is the state result,
with the national result alongside: `FutureQuantiles.national`, plus the table columns
`national_change_factor`, `national_future_flow_cfs` and `national_extrapolated`.
`provenance["basis"]` names the rule that applied, and `provenance["national_alongside"]` names
the national set. Since no national set ships, `national` may be `None`; a state set must then
apply.

## Per-state status (Waves 1 and 2)

Surveyed 2026-10-06 against the primary documents. **No Wave 1 or Wave 2 state publishes a
tabulated per-AEP flood-flow multiplier**, so no `<STATE>.factors.json` ships and
`available_factor_sets(level="state")` is empty. The national default, which is itself empty,
applies, and a user supplies a site-specific set. Washington is the one state that *requires*
a future-flow factor, but it is read per site from a web tool, not from a table.

| State | Guidance found | Method type | Legal status | Encoded |
|---|---|---|---|---|
| Washington (WA) | WSDOT Hydraulics Manual (2026): 2080s 100-year flow = current 100-year flow x the WDFW tool's per-site percent change | Multiplier, one AEP (0.01), per site | Required for WSDOT water crossings (unless the State Hydraulics Office finds it impracticable) | No tabulated factors |
| Oregon (OR) | ODOT Hydraulics Manual (2026): climate resilience listed as a site-assessment factor | Qualitative | Informational | No tabulated factors |
| Idaho (ID) | ITD Bridge Hydraulics Manual (2021): rainfall-runoff modelling suggested for climate-related rainfall change | Qualitative (simulation suggested) | Informational | No tabulated factors |
| Montana (MT) | None in MDT Hydraulics Manual Ch. 9 (2026); USGS trends study is historical | None | None | No tabulated factors |
| Colorado (CO) | CDOT Resilience Improvement Plan (2025): future IDF framework to be tested; not yet in the Drainage Design Manual | Projected precipitation (planned) | Informational | No tabulated factors |
| Utah (UT) | None in UDOT Drainage Manual (2025) | None | None | No tabulated factors |
| Wyoming (WY) | None in WYDOT Road Design Manual culvert section; USGS trends study is historical | None | None | No tabulated factors |
| New Mexico (NM) | None found; NMDOT Drainage Design Manual (2018) not retrievable for this survey | None found | None found | No tabulated factors |
| Arizona (AZ) | ADOT Hydrology Manual (2014): stationarity screening of gauged records only | None | None | No tabulated factors |
| Nevada (NV) | NDEP climate vulnerability study (2024): annual runoff, not flood quantiles; NDOT Drainage Manual (2006) not retrievable | None | None | No tabulated factors |

Revisit a state when its DOT manual is revised, when CDOT publishes its future-IDF framework,
or when NOAA Atlas 15 Volume 2 is published.

## Wave 1: Columbia River basin

### Washington (WA): #37
Status: **no tabulated factors** (Surveyed 2026-10-06). The guidance is quantitative and
**required**, but per site.

- **WSDOT Hydraulics Manual M 23-03.12 (April 2026), Ch. 7 Water Crossings, §7-3.5.5 Asset
  Management Preservation (pp. 7-38 to 7-39).** The procedure: (1) take the "percentage change
  in 100-year flood event" from WDFW's Climate-Adapted Culvert Design tool; (2) apply "the
  projected increase in 2080" to the current 100-year design flow; (3-5) model the 2080
  projected 100-year flow, water-surface elevation and scour. §7-3.6.1, Table 7-3 note a
  (p. 7-40): "The 2080 100-year projected flood event shall be used for the design, unless the
  State Hydraulics Office has determined that the 2080 100-year projected flood event is not
  practicable."
- **WDFW tool.** Mauger, G.S., Liu, M., Adam, J.C., Won, J., Wilhere, G., Dulan, D., Atha, J.,
  Helbrecht, L., and Quinn, T., 2021 (report dated July 11, 2022), *New Culvert Projections for
  Washington State: Improved Modeling, Probabilistic Projections, and an Updated Web Tool*,
  UW Climate Impacts Group, doi:10.6069/31T3-RE28. The tool reports per-site projected percent
  changes in bankfull width, bankfull flow and 100-year peak flow (Figure 15); the report
  tabulates no statewide or regional flow factors.
- **WDFW 2017 report.** Wilhere, G., Atha, J., Quinn, T., Helbrecht, L., and Tohver, I.,
  November 2017, *Incorporating Climate Change into the Design of Water Crossing Structures*,
  Final Project Report. Table 3 (p. 18) gives mean percent change in bankfull discharge and
  width by ecoregion division for the 2040s and 2080s. Not a factor set: it is the bankfull
  flow (recurrence interval 1.2-1.5 years by ecoregion, Table 2, p. 15), not the 100-year flood
  WSDOT designs to, and it averages over grid cells that the procedure reads individually.

Fields: legal status required (WSDOT) · method type multiplier at one AEP (0.01) · scenario
the tool's GCM ensemble · horizon 2080s · geography per site (tool grid cell). A user who reads
the tool's percent change for a site builds
`ChangeFactorSet(aeps=(0.01,), factors=(1 + pct/100,), horizon="2080s", level="site", ...)`;
`apply_change_factors` then holds that factor at every other AEP and flags them extrapolated.

### Oregon (OR): #38
Status: **no tabulated factors** (Surveyed 2026-10-06). Qualitative only.

- ODOT Hydraulics Manual, Ch. 7 Culvert Design (effective March 2026), Table 7-3 "Site
  Assessment Key Factors" (p. 11) lists "Climate Resilience: increased frequency and intensity
  events" among hydrology risks, and Table 7-12 "Alternative Analysis Factors" (p. 30) lists "resilience to climate
  change" under hydraulic performance. No numbers.
- ODOT Hydraulics Manual (2014) Ch. 7 Hydrology: no future-climate provision.

Fields: legal status informational · method type qualitative · scenario/horizon none ·
geography statewide.

### Idaho (ID): #39
Status: **no tabulated factors** (Surveyed 2026-10-06). Qualitative only.

- ITD Bridge Hydraulics Manual (July 2021), §4.4.3 Rainfall-Runoff Hydrograph Methods
  (p. 28): rainfall-runoff methods are suggested "when it is important to estimate the effects
  of future changes such as land use in the watershed or climate-related changes in rainfall."
  No factors.
- ITD Roadway Design Manual §600 Hydraulics: no future-climate provision.

Fields: legal status informational · method type simulation suggested · scenario/horizon
none · geography statewide.

### Montana (MT): #40
Status: **no tabulated factors** (Surveyed 2026-10-06). None published.

- MDT Hydraulics Manual Ch. 9 Hydrology (July 2026): no future-climate provision (future land
  use only).
- Sando, S.K., Barth, N.A., Sando, R., and Chase, K.J., 2025, *Peak streamflow trends in
  Montana and northern Wyoming and their relation to changes in climate, water years
  1921-2020*, USGS SIR 2023-5064-G: historical trends and change points (mostly downward), no
  projections.
- NOAA Atlas 15's 2024 Montana pilot is precipitation, not flow (see National).

## Wave 2: Colorado River basin

### Colorado (CO): #93
Status: **no tabulated factors** (Surveyed 2026-10-06). Projected precipitation is planned,
not published as design guidance.

- CDOT Resilience Improvement Plan (February 10, 2025), §3.1 (p. 19): "CDOT will test a future
  precipitation framework to develop future conditions IDF curves that modify the hydrologic
  and hydraulic design protocols ... The results will be used to update ... [the] Drainage
  Design Manual"; §3.1.1 (p. 20) and §6.1 (p. 46) name the CWCB/FEMA *Statewide Precipitation
  Scaling for Future Conditions* (Michael Baker International, 2024) as a study CDOT will
  review. That study is precipitation, and it could not be located for this survey.
- CDOT Drainage Design Manual (2019): not retrievable here (served through a document viewer).
  The RIP above says the future framework is not yet in it.
- NCHRP 15-61 Tables 6.6-6.9 (Denver) are a worked example, reproduced by
  `ensemble_ratio_summary`, not Colorado guidance.

### Utah (UT): #96
Status: **no tabulated factors** (Surveyed 2026-10-06). None published.

- UDOT Drainage Manual of Instruction (August 2025): no mention of climate change or future
  flows.

### Wyoming (WY): #98
Status: **no tabulated factors** (Surveyed 2026-10-06). None published.

- WYDOT Road Design Manual §3-04 Culvert Design (August 2011): no future-climate provision.
  The WYDOT Bridge Design Manual has no hydraulics chapter.
- USGS SIR 2023-5064-G (see Montana) covers northern Wyoming: historical trends only.

### New Mexico (NM): #100
Status: **no tabulated factors** (Surveyed 2026-10-06). None found, with one gap.

- NMDOT Drainage Design Manual (2018): the NMDOT site serves it through a script and it could
  not be retrieved for this survey, so its text is unchecked. No other NMDOT future-flow
  guidance was found.
- NCHRP 15-61 Table 8.1 (Pajarito Creek at Newkirk, 07225000) is a worked example, reproduced
  by `index_flood_projection`, not New Mexico guidance.

### Arizona (AZ): #102
Status: **no tabulated factors** (Surveyed 2026-10-06). None published.

- ADOT Highway Drainage Design Manual, Hydrology, 2nd ed. (January 2014), §10.2.6 (p. 10-7):
  screen gauged records for nonstationarity by plotting; formal tests "are beyond the scope of
  this manual." No future factors.

### Nevada (NV): #104
Status: **no tabulated factors** (Surveyed 2026-10-06). None found, with one gap.

- NDEP, *Climate Vulnerability Evaluation of Nevada's Watersheds*, Final Draft (December
  2024), Table 10 (p. 36): projected change in **annual runoff** by HUC12, 2061-2090 versus
  1971-2000. Not a flood-quantile factor.
- NDOT Drainage Manual (December 2006): the NDOT site refused automated access (HTTP 403), so
  its text is unchecked. Its date predates current climate guidance.
