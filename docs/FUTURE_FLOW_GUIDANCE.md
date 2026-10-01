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

NCHRP 15-61 Table 7.1 is used in `tests/test_future_flow.py` to check that
`apply_change_factors` reproduces a published per-AEP worked example. It is not shipped as
data: the source labels it a hypothetical precipitation increase for one watershed.

Revisit when NOAA Atlas 15 Volume 2 is published, or when FHWA or NCHRP publishes a
tabulated national product.

## Wave 1: Columbia River basin

### Washington (WA): #37
Status: **pending**. Survey the WSDOT Hydraulics Manual, Ecology, the state climate office,
and USGS/state cooperative studies.

### Oregon (OR): #38
Status: **pending**. Survey the ODOT Hydraulics Design Manual, OWRD/DEQ, the state climate
office, and USGS/state cooperative studies.

### Idaho (ID): #39
Status: **pending**. Survey ITD design guidance, IDWR, the state climate office, and
USGS/state cooperative studies.

### Montana (MT): #40
Status: **pending**. Survey MDT hydraulics guidance, DNRC, the state climate office, and
USGS/state cooperative studies.
