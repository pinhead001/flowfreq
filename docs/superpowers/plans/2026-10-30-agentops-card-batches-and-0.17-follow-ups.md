# flowfreq: card batches for the 0.17 follow-ups — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Clear the follow-ups the 0.16.0 nationwide rollout left (issues #203–#218, TODO.md
"Nationwide rollout (2026-10-25 plan)") and ship 0.17.0, running the mechanical code through
the repo's own agentops delegate runner (`tools/agentops/`, `docs/AGENTOPS.md`, landed in #216
and #217) so the Claude Code lead spends its quota on specs, data reading, review and release.

**Architecture:** The lead writes failing tests on `spec/<id>`, a card in `.agentops/cards/`
with `base = "origin/spec/<id>"`, and later reviews `card/<id>` and opens the PR
(`docs/AGENTOPS.md`, "The loop"). Work splits three ways: **delegate cards** (code with a
named test, no data or network), **lead items** (investigation, USGS report reading, spec
writing, review, merge, release) and **user decisions** (hydrologist picks, agy permissions,
app merges, GitHub actions the classifier blocks).

**Tech Stack:** `tools/agentops/` (in-repo port of flowfield2d's; the shared `D:/a/agentops`
package is the same lineage, and moving flowfreq onto it is a separate decision), `agy`
1.4.10, flowfreq (pip/setuptools, Python 3.11–3.14, `make check`), GitHub CLI.

**Spec:** `TODO.md` ("Open Items" → "Nationwide rollout" and "Agent ops"), GitHub issues
#203–#218, `docs/MASTER_ROADMAP.md` §4, `docs/AGENTOPS.md`. There is no separate design
document. This file replaces the 2026-10-29 plan of the same name, which was written against a
checkout seven commits behind `main` and re-planned tooling that already existed.

## Global Constraints

- Cards follow `docs/AGENTOPS.md`: tests first on `spec/<id>`, prompt names every file to read
  and edit, `acceptance = "python -m pytest <files> -q -p no:cacheprovider"`. The runner puts
  the worktree on `PYTHONPATH` itself (`child_env`); do not use `PYTHONSAFEPATH=1 pytest`.
- Delegates get code and synthetic fixtures only: no card downloads data, reads a USGS
  report, or needs `USGS_API_KEY`. Anything that does is a lead item.
- Model names come from `agy models` on the day (2026-10-30: `gemini-3.8-flash-*`,
  `gemini-4-argon-*`, `claude-sonnet-5-5-*`, `claude-opus-5-5-*`, `gpt-oss-120b-medium`).
  Weekly capacity is about 1.2M Gemini tokens and about 0.33M third-party tokens.
- **The host agy settings must carry allow/deny rules before any `cards run`.** Without them
  every call is auto-denied and still billed (459k tokens on 2026-10-30). The proposed file is
  in the session scratchpad as `agy-settings.proposed.json`; the user applies it.
- Known-failing tests are `xfail(strict=True)`, never skipped. A live-service outage (#218) is
  not a known failure.
- The Fortran parity job is pinned to `ubuntu-24.04` (#223); move the pin only with the goldens.
- Release 0.17.0 comes from `main`.

## Review Focus

1. **A `region=` list whose weights sum below 100 %** is refused with the listing, never
   area-weighted over part of a basin. — `spec/i203-multi-region`.
2. **A `region=` naming a library region when NSS located nothing** must say in provenance
   that the caller chose it. — `spec/i207-region-unlocated`.
3. **A HUC8 whose WBD `states` is `CN,MT`** keeps the Canadian code; the reader returns
   `{"CN", "MT"}`. — `spec/i205-huc8-states`.
4. **A Water Data feature without `country_code`** (older cached payloads) is still US; a
   feature with `country_code = "CA"` and `state_code = "11"` is not DC. —
   `spec/i206-catalog-outside-us`.
5. **Silencing the advisory leaves `StateLibrary.advisory` set.** — `spec/i214-snapshot-warning`.

---

## Inventory

| # | Item | Source | Owner | Where |
|---|---|---|---|---|
| 1 | agy permission rules on the host | AGENTOPS.md step 3 | **user** | scratchpad `agy-settings.proposed.json` |
| 2 | Claude Code hooks and status line | AGENTOPS.md step 2 | **user** | `tools/agentops/settings.example.json` |
| 3 | Cards i214, i210, i203, i206 (country_code), m12, m4, m7 | #217 | runner, then lead review | `.agentops/cards/`, queue pending |
| 4 | Cards i205-huc8-states, i207-region-unlocated | this plan | runner, then lead review | `.agentops/cards/`, queue pending |
| 5 | Rebuild the catalog for DC; run `tools/build_huc8_states.py` live | after i206, i205 merge | lead | — |
| 6 | #211 region-indicator forms: schema decision (interaction terms) | issue comment 2026-10-30 | **user** then lead | #211, #208 |
| 7 | #209 per-state default set; MT and NM picks | decision packet 2026-10-30 | **user** | #209 |
| 8 | #213 `sep_log` LA/AR/MS; #204 HUC rows for 15 skew studies | reports | lead | — |
| 9 | #218 live tests (WA PRECIP raster) | rerun 2026-10-30 failed | lead, after two more weekly runs | #218 |
| 10 | Dependabot #220 (needs `workflow` scope) | PR | **user** | GitHub UI |
| 11 | v0.10.4 GitHub release recreation | TODO small items | **user** | `gh release create v0.10.4 --verify-tag` |
| 12 | Release 0.17.0; app pin bump PR | flowfreq-state | lead; **user** merges the app PR | — |
| 13 | Move flowfreq onto the shared `D:/a/agentops` package | agentops README | **user decision** | — |

---

### Task 1: Run the queue once the agy rules are in place

- [ ] `python tools/agentops_cli.py --config agentops.toml gov profile max` (the committed default
  is `pro`, which `tests/agentops/test_flowfreq_setup.py` pins; the switch is a local, uncommitted
  edit), then `gov plan` shows `dispatch` for i203 (Sonnet) and the first Gemini card.
- [ ] `cards run --card i214-snapshot-warning` first, alone. Expected: `done`, ledger `ok: true`,
  acceptance tail `passed`. If the ledger shows `ok: false` with an empty response, stop: the
  permission rules are not in effect.
- [ ] `cards run` for the rest (concurrency from the governor). i207 is not queued: `depends_on`
  only orders dispatch, and its worktree is cut from `spec/i207-region-unlocated`, which lacks
  i203's code. After `card/i203` merges, recut that spec branch from `main` (rebase its two test
  commits), then `cards add` i207.
- [ ] For each `done` card: `git diff main...card/<id>`, the `reviewer` subagent, push, PR with
  `Closes #<n>`, CI, merge in order: i214, i210, i206, i203, i207, i205, then m12, m4, m7.

### Task 2: Lead follow-through after merges

- [ ] After i206: `USGS_API_KEY=... python tools/build_gage_catalog.py --states DC --refresh`;
  confirm 15388960 is gone and the DC row of the smoke report says "no USGS gage in the
  catalog"; commit the catalog and `gage_catalog.meta.json`.
- [ ] After i205: `python tools/build_huc8_states.py`; expect about 2,300 rows, every current
  row preserved with the same `states`; drop the "lists only some HUC8s" caveats in
  `flowfreq/regional_skew.py` (two comments) and the roadmap; PR.
- [ ] After i210: pick downstream points for AR, GA, ND, PA, PR, TX on the gage pages into a
  `points.csv`; rerun the smoke for those states plus FL, MI, NE (`--resume --states ...`).
- [ ] After i207: rerun the smoke for WY and NV with `--region` (WY: the report's region at the
  gage; NV: region 6 per `docs/NEVADA_REGION6_SCHEMA.md`).

### Task 3: Reports (lead)

- [ ] #213: SIR 2021-5046 (alluvial plain, LA/AR/MS), SIR 2024-5031 (LA coastal plain), TR 60
  (1998), WRIR 85-4150: read each standard-error table; set `sep_log` on the 32 equations with
  `notes` citing the page; a per-state assertion in `tests/test_regression_nss_snapshot.py`.
- [ ] #204: HUC rows for the 15 skew studies into `regional_skew_hucs.csv`, smoke states first.

### Task 4: Release 0.17.0

- [ ] CHANGELOG Unreleased → `0.17.0`; `pyproject.toml` `0.17.0`; PR; CI; merge; annotated tag
  `v0.17.0` pushed (release.yml publishes; never `gh release create`); `release/0.17.x`;
  `main` → `0.18.0.dev0`; CLAUDE.md Branches section.
- [ ] flowfreq-app PR: pin `@v0.17.0`; ungaged page passes a list to `region=` (i203) and lets
  the user name a library region for WY/NV (i207). The user merges and deploys.

---

## What happened on 2026-10-30 (so the next session does not repeat it)

- A plan and PR (#222, closed) were built from a stale local main; the tooling they added
  already existed (#216, #217). The surviving pieces are this plan, two new spec branches and
  cards, the amended #206 spec (country code), and the findings on issues #206–#211, #218.
- Main CI went red when `ubuntu-latest` moved to 26.04; fixed by pinning the parity job (#223).
- The first `cards run` burned 459k Gemini tokens on auto-denied calls (no host agy allow
  rules) and three zero-token Sonnet calls (stale model name). All cards are reset to pending.
- Dependabot #219 and #221 merged; #220 needs the `workflow` OAuth scope (user).
- The local flowfreq-app checkout was on a stale branch; it is now on `main` at the v0.16.0 pin.
