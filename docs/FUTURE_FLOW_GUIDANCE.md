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

Status: **pending** (#36). Candidate sources to transcribe: FHWA HEC-17; NCHRP 15-61;
NOAA Atlas 15 future precipitation once published.

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
