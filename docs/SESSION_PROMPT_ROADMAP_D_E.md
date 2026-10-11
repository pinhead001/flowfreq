# Session prompt: Roadmap Phases D–E, F screening, MT skew, release

Paste everything below the line into a **local** Claude Code CLI session, started at the
root of a clean, up-to-date `flowfreq` checkout (`git checkout main && git pull`). Local
matters: tag pushes, legacy NWIS and long live-data jobs work there and not from a web
session (TODO.md, Environment constraints). Set `USGS_API_KEY` first if you have one; the
batch and Montana items make thousands of Water Data API requests.

Decisions already made (2026-10-11): one PR per item, merged in dependency order; Montana
skew improved but kept out of `regional_skew.csv`; flowfreq-app pin bumped in its own PR;
release v0.11.0 at the end. Edit the **Merge policy** line if you want to approve each
merge yourself.

---

You are the **orchestrator** of a multi-agent session in the `flowfreq` repository (this
working directory). You plan, dispatch, review and merge; workers write the code. Read
`CLAUDE.md`, `TODO.md` (Status, Environment constraints, Open Items) and
`docs/MASTER_ROADMAP.md` §2, §4, §5, §6.1 and §1.3 before dispatching anything.

## Ground rules (every agent, every item)

- `CLAUDE.md` is binding. In particular: never edit `vendor/`; `pip install -e ".[dev]"`
  (black is pinned `<25`); run `make check` (black, isort, mypy, pytest, in CI's order)
  before every push, with `PYTHONSAFEPATH=1` if invoking pytest by hand; Python 3.11
  syntax only; type hints and NumPy docstrings on public API; `logging`, never `print`, in
  library code; tests for every public function (happy path and one error case).
- **Validate against published numbers, not just synthetic data.** Every method item
  names a published target below. Reproduce it to the source's printed precision and pin
  it in a test. Where something does not reproduce, record the exact gap and its cause in
  the PR and the relevant doc. Do not widen a tolerance to pass; a known failure is
  `xfail(strict=True)` and the PR says so.
- **Existing numbers must not move.** If a change alters any existing quantile, skew,
  bound or test expectation, stop and report to the orchestrator, who asks the user.
- No new runtime dependency without asking the user. Test-only cross-checks may use what
  is already in the `dev` extra.
- Live-service code: mark tests `requires_network`; commit trimmed fixtures so the
  offline suite covers the logic; long data pulls are cached and resumable under a
  gitignored cache directory, as `tools/wave_qa.py` and `tools/build_gage_catalog.py` do.
- Each PR adds its own bullet(s) to `CHANGELOG.md` → Unreleased and ticks or updates only
  its own lines in `TODO.md` and `docs/MASTER_ROADMAP.md`. New modules go in `CLAUDE.md`'s
  Architecture list. Do not reword other items' entries.
- Commit messages and PR bodies: what changed, why, what was validated against what, and
  the measured `make check` result. No model names anywhere in the repo.

## How to run the agents

- Dispatch each item to a worker with the Agent tool, `isolation: "worktree"`, in the
  background. Give each worker a **self-contained** brief: copy the item's section below,
  the ground rules, and the branch name. Workers start cold; they do not see this prompt.
- At most **4 workers at once**, so CI and the Water Data API rate limit are not swamped.
- Worker procedure: branch from the latest `origin/main`; implement; `make check` green;
  push; open a PR (use the repo's PR template if one exists); report back the PR number,
  the validation results, anything not reproduced and any open question.
- Orchestrator on each returned PR: read the diff adversarially (or run `/code-review`),
  check CI, send fixes back to the same worker (SendMessage) rather than starting a new
  one. Before merging, merge `origin/main` into the PR branch (a merge commit, never a
  rebase or force-push) and resolve conflicts in the shared files (`CHANGELOG.md`,
  `TODO.md`, `docs/MASTER_ROADMAP.md`, `CLAUDE.md`, `flowfreq/__init__.py`,
  `pyproject.toml`'s mypy overrides) by keeping both sides.
- **Merge policy: merge yourself once CI is green on the latest head and your review is
  clean.** Always ask the user before: pushing a tag, opening a PR in another repository,
  merging anything that moves existing numbers, or adding a dependency.
- Start a dependent item only after its prerequisites are **merged**, from fresh `main`.
- Keep a running status table (item, branch, PR, state, blocker) and show it to the user
  whenever a wave finishes.

## Dependency order

```
Wave 0  S0 housekeeping                               (orchestrator, first, alone)
Wave 1  A1 quantile variance -> A2 App. 9 weighting   (A2 starts when A1 merges)
        B  MOVE record extension (App. 8)
        C  state at-site batch (Idaho)
        D  stationarity screening
        E1 ramp window   E2 mixed population   E3 Crippen-Bue envelope
        M  Montana skew study (long-running; start early, finishes whenever)
Wave 2  F  regression-weighted DAR + regulation guardrail   (needs A2)
        C2 weighted columns in the Idaho batch              (needs A2 and C)
Wave 3  G  estimate_at_site, WA pilot                       (needs A2, B, F)
Wave 4  H  docs refresh -> release v0.11.0 -> flowfreq-app pin PR   (needs everything)
```

Wave 1 has nine items; run them four at a time, starting A1, C and M first (longest),
then fill slots with B, D, E1–E3. M does not block anything; if it is still running at
Wave 4, release without it and say so.

## Items

### S0 — Housekeeping (orchestrator, branch `docs/roadmap-refresh-2026-10`)
- If branch `claude/todo-list-8mwbt3` is unmerged, open its PR and merge it (TODO.md
  refresh for #110–#125), or fold its change into this branch.
- `docs/MASTER_ROADMAP.md` stale entries: §7's "CI: `requires_network` live checks on a
  scheduled workflow" is done (#114, `.github/workflows/live.yml`); §3.2's Wave 2 row says
  NV has "no published range table", but #120 stored the data-section ranges. Fix both,
  and refresh the "status as of" date.
- Record the baseline: `make check` on fresh `main`, with the test count, in the PR body.

### A1 — At-site quantile variance (branch `feat/quantile-variance`)
`FrequencyResults` carries quantiles and confidence limits but not the sampling variance
of each log10 quantile, which Appendix 9 needs. Surface it per AEP on both engines:
- Fortran engine: `emafitpr`'s `var_est` (see `vendor/peakfqr/src/emafit.f` and
  `R/fortranWrappers.R` for its meaning and units).
- Native engine: the same quantity from the ported `var_mom`/`var_emab` path
  (`flowfreq/_var_mom.py`, `flowfreq/_var_emab.py`).
- Add it as a new column (or field) without changing any existing value. Validate native
  against `emafitpr` on the golden-file sites and all 24 WY/MT `.psf` stations
  (`requires_fortran`), to the agreement the rest of the parity suite achieves; state the
  tolerance and why.

### A2 — Bulletin 17C Appendix 9 weighting (branch `feat/b17c-app9-weighting`, after A1)
New module (e.g. `flowfreq/weighting.py`):
- `log Qw = (log Qg·Vr + log Qr·Vg) / (Vg + Vr)`, `Vw = Vg·Vr / (Vg + Vr)`, per AEP, with
  Vg from A1 and Vr from the regression equation. Make the regression prediction variance
  public (`regression/equations.py::_prediction_sd` today) and report which basis was used
  (site-specific `MEV + x'Ux`, AVP, or `sep_log` fallback).
- Inputs: an at-site `FrequencyResults` and `RegressionEstimate`s (single region or
  `evaluate_weighted` multi-region). Refuse extrapolated regression estimates unless
  explicitly allowed. Output: weighted flow, interval, weights, equivalent years of
  record, full provenance.
- **Validation:** WA SIR 2016-5118 Table 8 (report eq. 10, "WIE") gives at-site,
  regression and weighted estimates with 95 % intervals per gage. `data/regression/WA.json`'s
  notes show Vr was already inverted from that table to 1.000–1.002 of the stored
  covariance for Regions 1–3. Reproduce the weighted estimates for Regions 1–3 gages from
  the report's at-site values and Vs and the stored equations; then at a few gages with
  flowfreq's own at-site fit, explaining differences. Region 4 has no stored covariance
  (sep_log fallback); show and document the effect. Cross-check against the B17C text of
  Appendix 9.

### B — Record extension, MOVE.1 / MOVE.3 (branch `feat/record-extension`)
B17C Appendix 8. New module (e.g. `flowfreq/record_extension.py`): MOVE.1 and MOVE.3
(Hirsch 1982; Vogel & Stedinger 1985) in log space, concurrent-period diagnostics
(correlation, the minimum correlation for an effective gain), effective record length, and
the extended record handed to EMA as B17C Appendix 8 prescribes. Read Appendix 8 for how
extended values enter EMA and the variance; do not invent a convention. Validate against
B17C's Appendix 8 example and at least one other published MOVE example; pin both.

### C — State at-site batch, Idaho first (branch `feat/state-atsite-batch`)
Roadmap §2 "Batch B17C for a whole state" and §7 per-state validation. ID is `verified`
with the PNW regional skew. Build a resumable tool (e.g. `tools/state_atsite.py`, logic in
the library where reusable) that selects the report's gages from the catalog, pulls peaks
(Water Data API, peak codes applied as peakfq would), applies the report's regional skew
and settings, and produces a table in the form of SIR 2016-5083's at-site appendix. Diff
against the published values and categorise each disagreement: different period of
record, perception thresholds or historic info not in NWIS, low-outlier handling, data
revisions, unexplained. Write `docs/STATE_ATSITE_ID.md` with the summary and commit the
comparison CSV. A useful bar: matching sites agree to the report's printed precision, and
every mismatch has a stated cause or is listed as unexplained.

### C2 — Weighted columns in the Idaho batch (branch `feat/state-atsite-weighted`, after A2 and C)
Add A2's weighted estimates to C's table and compare against SIR 2016-5083's weighted
values.

### D — Stationarity screening (branch `feat/stationarity`)
Roadmap §6.1. New module (e.g. `flowfreq/stationarity.py`): Mann-Kendall with tie
correction (and an autocorrelation-corrected variant), Sen's slope with confidence
interval, Pettitt change point, and `stationarity_report(peaks)` returning a structured
summary plus an optional plot. Screening only: it never changes a fit. Validate against
published worked examples (Helsel et al. 2020, USGS TM 4-A3, has MK, Sen and Pettitt
examples) and one USGS report's per-gage trend table; pin them.

### E1 — Fixed-window ramp metric (branch `feat/ramp-window`)
TODO.md Sub-daily, "No fixed-window ramp metric". Add `window_hours=` to
`subdaily.ramping_rates`. Default behaviour unchanged. Reproduce Exelon (2012) RSP 3.8
Studies 9–11 (2.8 / 4.2 / 5.6 ft stage decline in the first hour after reduction to
minimum flow, USGS 01578310 gage height). Explain Study 12 (3.29 vs 3.1 ft) or record it
as not reproduced. Fixture-backed offline test plus `requires_network` live test.

### E2 — Mixed-population detection (branch `feat/mixed-population`)
Roadmap §2: at minimum detect and warn. Use peak timing (circular statistics on peak
dates: a seasonally bimodal distribution) and, where informative, magnitude by season.
Surface it as a warning and a structured flag in `run_ffa`'s result, never as a change to
the fit. Validate on gages a USGS report describes as mixed population (snowmelt/rain in
the West; cite the report) and on gages it treats as single population. Report false
positive and false negative behaviour honestly.

### E3 — Crippen–Bue envelope check (branch `feat/envelope-check`)
Roadmap §2: a flag, not a correction. Transcribe the Crippen & Bue (1977, WSP 1887)
regional envelope curves (regions, equation form and coefficients; double-enter them, as
the regression files were). The caller supplies the region; no polygon dependency. Flag
any quantile above the envelope for its drainage area. Reproduce the curves' published
values at tabulated drainage areas and pin them.

### M — Montana regional skew, improved provisional study (branch `data/montana-skew-v2`)
Keep it **out of** `flowfreq/data/regional_skew.csv`; MT stays `pending` and the roadmap
item stays open for a USGS study. Improve `docs/MONTANA_REGIONAL_SKEW_PROVISIONAL.md`
(v1: CONSTANT G = -0.03, MSE 0.36, σ²δ = 0.35, ERL 23 yr) to a v2 that fixes its known
gaps, using `flowfreq.skew_study` and `tools/build_montana_skew_study.py`:
- historical periods taken from published sources (the vendored
  `wymt_ffa_2022A.psf`, SIR 2025-5019's data release) instead of approximated from NWIS,
  which is why MBV* is 9.3 rather than 10 in the PNW reproduction;
- basin-level explanatory variables from StreamStats (mean basin elevation, precipitation,
  drainage area, and whatever else MT computes), cached and resumable, to test whether
  any explains σ²δ where gage altitude explained only 5 %;
- regulation screening through `flowfreq.regulation` (current NID option) and peak codes.
Report v1 vs v2 side by side, including if nothing improves. The PNW reproduction must
still reproduce SIR 2016-5083 Tables B2/B3 after any tooling change.

### F — Regression-weighted DAR and regulation guardrail (branch `feat/regression-weighted-dar`, after A2)
Roadmap §5.1 and §5.2. In `transpose.py`: `Q_u,w = Q_u,r·[R − (2|A_g−A_u|/A_g)(R−1)]`,
`R = Q_g,w/Q_g,r`, valid for `0.5 ≤ A_u/A_g ≤ 1.5`, using A2's weighted gage estimate.
Each state's exact variant comes from its report; implement WA's first (SIR 2016-5118)
and structure for others. Reproduce the report's ungaged-near-gage worked example. Add the
§5.2 regulation guardrail: refuse a regulated donor for an unregulated target and the
reverse, unless forced, using `regulation.classify_site`, in `transpose` and `qppq`.

### G — `estimate_at_site`, proven on the WA pilot (branch `feat/estimate-at-site`, after A2, B, F)
Roadmap §4 / Phase E. One entry point that, given a location (or site number):
- on a gage: B17C at-site with its regional skew, weighted with regression (A2);
- near a gage on the same stream within the state's area band: regression-weighted DAR (F);
- otherwise: regression (offline library, with NSS cross-check where available), with
  prediction intervals.
Optional MOVE extension (B) for short records. Every result carries provenance: method and
why, donors, equations, skew source, regulation screen, extrapolation flags. Prove it end
to end on three WA sites (one gaged, one near a gage, one ungaged) against the report's
own numbers, and add a vignette under `docs/vignettes/`.

### H — Docs, release, app pin (orchestrator, after everything else merges)
1. Docs refresh PR: TODO.md (Status with a freshly measured test count, Open Items),
   roadmap checkboxes and §9 phase status, CLAUDE.md (Architecture, test runtime), README.
2. Release PR `release/v0.11.0`: version bump, CHANGELOG `[0.11.0]` section listing the
   new APIs and anything that changed behaviour. After it merges, **ask the user**, then
   tag `v0.11.0` on the merge commit and push the tag.
3. `pinhead001/flowfreq-app`: **ask the user** before opening the PR. Bump
   `requirements.txt` from `v0.7.0` to `v0.11.0`, run the app's suite before and after,
   and in the PR body list every change since v0.7.0 that moves the app's numbers
   (TODO.md Downstream lists them through v0.10.2: required regional skew choice, OGC
   peaks with UTC `peak_date`, peak codes on by default, partial dates kept, the EMA
   fixes, one peak per water year), with before/after values for at least one site.
   Follow `flowfreq-app/CLAUDE.md`'s "one edit, deliberate" note.

## Finish

Report: the final status table, every PR with its validation outcome, what did not
reproduce and why, open questions, and what remains on the roadmap.
