# Agent ops: a Claude lead with agy delegates

flowfreq's development uses the same multi-agent setup as flowfield2d. The Claude Code lead
plans the work, writes the failing tests and a task card, and later reviews and merges. Headless
delegates, Google's Antigravity CLI (`agy`), implement the cards, each in its own git worktree.
A governor reads the lead's Claude usage meter and the delegates' token ledger, then decides
how much work each side takes on.

The code is in `tools/agentops/`, a port of `flowfield2d.agentops`. Its design record is
flowfield2d's `docs/adr/0002-agy-integration.md`: the agy envelope, the success rule,
confinement and the quota pools. The tool is development tooling only. It is not in the
published `flowfreq` package, and nothing under `flowfreq/` imports it.

## Pieces

| Piece | Where | What it does |
|---|---|---|
| Config | `agentops.toml` | Paths, governor profiles, agy models and pools, `child_env`, hook guards |
| Claude meter | `statusline.py` | The status line writes `rate_limits.five_hour` / `seven_day` to `.agentops/state/claude_usage.json` |
| Governor | `governor.py` | Sets the mode (`normal` / `conserve` / `review-only` / `wait`) from the meter and the active profile, and sets agy concurrency from recent errors, calibration and cooldowns |
| Ledger | `ledger.py` | `.agentops/state/ledger.jsonl`, one record per agy call (tokens, duration, ok) plus quota events |
| Usage report | `usage.py` | Claude 5h/7d use, plus agy calls and tokens per quota pool over 24 h and 7 days, failed calls, quota events, top cards |
| Cooldowns | `cooldown.py` | Starts a per-pool cooldown on a quota hit (parsed from "Resets in …"), then falls back to the next free pool |
| Queue | `queue.py` | Card TOMLs in `.agentops/state/queue/{pending,running,done,escalated}/` |
| Runner | `runner.py` | Per card: worktree, agy, then acceptance; commits to `card/<id>` on success, retries with the failure output, otherwise escalates. Never merges |
| Hooks | `hooks.py` | SessionStart prints the governor status. PreToolUse blocks force-pushes, pushes to `main`/`release/*`, tag deletes, and teammate/workflow spawns outside `normal` |

## One-time setup (developer machine)

1. `pip install -e ".[dev]"`, which now includes `filelock`. Then install agy and sign in.
2. **Wire the hooks and status line.** Merge `tools/agentops/settings.example.json` into
   `.claude/settings.json`, or into `.claude/settings.local.json` to keep it to one
   machine. The commands go through `tools/agentops.sh`. It picks `.venv/Scripts/python.exe`
   (Windows), then `.venv/bin/python`, then `python3`/`python`, so the Windows Store
   `python` shim is never used. A missing dependency fails open: the session keeps working
   and the status line reports the tool as unavailable.
3. **agy permissions** (`~/.gemini/antigravity-cli/settings.json`; Deny beats Ask beats Allow),
   as in ADR-0002 decision 3. Allow file read/write inside the worktrees. Deny `git push`,
   `merge`, `rebase`, `checkout`, `switch`, `reset`, `branch`, `worktree`, `remote` and
   `config`; any command using `../`; and `write_file` on the main checkout, `.git/` and
   `vendor/`. `AGENTS.md` tells the delegate not to use the shell at all, but models try
   anyway, and **without an allow rule a headless `command` is auto-denied and ends the turn**:
   the envelope says `SUCCESS` with an empty response and the call still costs its tokens
   (2026-10-30: eight calls, 459k Gemini tokens, nothing written). So allow the read-only
   commands too (`python -m pytest|black|isort|mypy`, `git status|diff|log|show`, directory
   listings), as the sandbox's `D:/a/sbx/state/gemini/antigravity-cli/settings.json` does.
   Model names change: check `agy models` before writing `model =` on a card
   (`claude-sonnet-4-6` disappeared between 2026-10-26 and 2026-10-30).
4. Select the governor profile for your Claude plan:
   `python tools/agentops_cli.py --config agentops.toml gov profile pro` (or `max`; use
   `overnight` when you are away and only review).

## The loop

```bash
A="python tools/agentops_cli.py --config agentops.toml"   # or: bash tools/agentops.sh
$A gov next                 # mode, queue state, and what to do next
$A cards add .agentops/cards/<id>.toml
$A gov plan                 # which ready cards would be dispatched now
$A cards run                # up to the governor's concurrency; or --card <id>
$A cards list
$A usage                    # Claude meter + agy tokens per pool
$A gov calibrate --agy-remaining 63   # after reading agy's /usage
$A cards reset <id>         # discard the worktree and branch, back to pending
```

**Lead (Claude) per card:**
1. On a branch `spec/<id>` from `main`, write the failing tests and push them. Do not open a
   PR for this branch: it exists only to hold the tests.
2. Write `.agentops/cards/<id>.toml` with `base = "origin/spec/<id>"` (run `git fetch
   origin` before `cards run`, so the remote branch is current), an outcome-level `prompt`
   that names the files to read and edit, and an exact `acceptance` command. Use
   `python -m pytest <the card's test files> -q -p no:cacheprovider`. Do not use
   `PYTHONSAFEPATH=1 pytest`: the runner runs commands through `cmd.exe` on Windows,
   which cannot parse that syntax.
3. Once the runner reports `done`, review `card/<id>`. The `reviewer` pass covers scope,
   tests not weakened, the full suite, `make check` and the CLAUDE.md conventions. Then
   push the branch and open the PR to `main` (the cards in `.agentops/cards/` are the
   committed record). `escalated` means the lead takes the card over or rewrites it.

**Why `child_env`:** the editable install points at the main checkout, so a worktree's
`python -m pytest` would import the main checkout's `flowfreq`. To prevent that, the runner
puts the card's worktree on `PYTHONPATH` (`{worktree}` in `agentops.toml`) for both the
delegate and the acceptance run, and also sets `PYTHONUTF8=1`.

## Card format

```toml
id = "i214-snapshot-warning"           # [A-Za-z0-9._-]
repo = "flowfreq"
title = "Silence the nss_snapshot advisory for sweeps (#214)"
base = "origin/spec/i214-snapshot-warning"   # holds the failing tests; git fetch first
model = "claude-sonnet-5-5-medium"   # optional; omit to take the first free pool in [agy] models
depends_on = []                        # card ids that must be done first
max_retries = 2                        # optional override
prompt = """
Outcome, not steps. Files to read: ... Files to edit: ...
Tests (already written, do not edit): tests/test_x.py::test_y
Constraints: ...
"""
acceptance = "python -m pytest tests/test_x.py -q -p no:cacheprovider"
```

## Caveats

- **Quota is per account, but the ledger is per project.** flowfield2d and flowfreq share the
  same agy pools, yet each project's governor sees only its own ledger. Calibrate from agy's
  `/usage` before a large run. A quota hit costs one failed call and starts the cooldown in
  that project only.
- A cloud (web) session has no status-line meter, so the governor uses the profile's
  `unknown_mode`. In `pro` that is `conserve`, which blocks `Workflow` and `TeamCreate` there.
- PreToolUse command patterns are best effort, the same residual risk ADR-0002 accepts.
  The server-side branch protection on `main` is the real guard.
