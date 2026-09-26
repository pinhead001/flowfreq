---
name: State rollout
about: One jurisdiction's regression, skew, validation, and future-flow work (MASTER_ROADMAP §3.2)
title: "Wave <N> state: <Name> (<XX>)"
labels: ["roadmap"]
---

Parent wave epic: #<epic>

Definition of done (`docs/MASTER_ROADMAP.md` §3.2):
- [ ] Every published peak AEP is in `flowfreq/data/regression/<XX>.json`, with citation, table number, and limits.
- [ ] ≥1 published worked example from the source report is reproduced to its printed precision.
- [ ] The online (NSS) and offline evaluators agree.
- [ ] Regional skew row(s) in `flowfreq/data/regional_skew.csv` are filled and second-checked.
- [ ] Region polygons are loaded, with a point-in-polygon test on a known site.
- [ ] The report's at-site B17C appendix is re-run and compared.
- [ ] State future-flow guidance is recorded in `docs/FUTURE_FLOW_GUIDANCE.md`. "None published as of <date>" is acceptable.
- [ ] `python tools/gen_regression_coverage.py` is rerun, and its output committed.

Cross-border neighbours to check:
