# Agent ops: a Claude lead with agy delegates

flowfreq's development uses the same multi-agent setup as flowfield2d. The Claude Code lead
plans the work, writes the failing tests and a task card, and later reviews and merges. Headless
delegates, Google's Antigravity CLI (`agy`), implement the cards, each in its own git worktree.
A governor reads the lead's Claude usage meter and the delegates' token ledger, then decides
how much work each side takes on.

The code is the shared `agentops` package in `D:/a/agentops` (github.com/pinhead001/agentops),
which flowfield2d and the other projects use too; flowfreq kept its own copy in `tools/agentops/`
until 2026-11-21, when the features only that copy had (`child_env`, `--settings`, `AGENTS.md`
seeding, `usage`, the machine-local profile) moved into the package. Its design record is
flowfield2d's `docs/adr/0002-agy-integration.md`: the agy envelope, the success rule,
confinement and the quota pools. The tool is development tooling only. It is not a flowfreq
dependency (CI does not install it), and nothing under `flowfreq/` imports it. flowfreq's part
is `config/agents.toml`, `config/agy-settings.json`, `run.bat`, `tools/agentops.sh` (the same
launcher for the hooks) and the cards in `.agentops/cards/`.

## Pieces

| Piece | Where | What it does |
|---|---|---|
| Config | `config/agents.toml` | Paths, governor profiles, agy models and pools, `child_env`, `settings`, hook guards (checked in CI by `tests/test_agentops_config.py`) |
| Claude meter | `statusline.py` | The status line writes `rate_limits.five_hour` / `seven_day` to `.agentops/state/claude_usage.json` and prints `5h <used>% (reset <in>) \| 7d <used>% \| gov <mode>`, plus `agy <n> run, <m> pend \| 7d tok <pool> <k>` while cards are queued or delegates ran in the last 24 h |
| Governor | `governor.py` | Sets the mode (`normal` / `conserve` / `review-only` / `wait`) from the meter and the active profile, and sets agy concurrency from recent errors, calibration and cooldowns |
| Ledger | `ledger.py` | `.agentops/state/ledger.jsonl`, one record per agy call (tokens, duration, ok) plus quota events |
| Usage report | `usage.py` | Claude 5h/7d use, plus agy calls and tokens per quota pool over 24 h and 7 days, failed calls, quota events, top cards |
| Cooldowns | `cooldown.py` | Starts a per-pool cooldown on a quota hit (parsed from "Resets in …"), then falls back to the next free pool |
| Queue | `queue.py` | Card TOMLs in `.agentops/state/queue/{pending,running,done,escalated,archived}/` |
| Runner | `runner.py` | Per card: worktree, agy, then acceptance; commits to `card/<id>` on success, retries with the failure output, otherwise escalates. Never merges |
| Hooks | `hooks.py` | SessionStart prints the governor status. PreToolUse blocks force-pushes, pushes to `main`/`release/*`, tag and branch deletes, and teammate/workflow spawns outside `normal`. Each `deny_bash` regex is searched in the whole command line, so the patterns stop at `;`, `&` and `\|` to stay inside the command that pushes |
| Review | `review.py` | `cards show` / `cards diff` a finished card; `cards accept` archives it, `cards reject --reason` sends it back, each with a ledger record |

## One-time setup (developer machine)

1. `pip install -e ".[dev]"`, then `.venv\Scripts\python -m pip install -e D:gentops` (the
   package needs Python 3.12, which the development venv is; it is not in `pyproject.toml`
   because CI tests 3.11-3.14 and has no `D:/a/agentops`). Install it into flowfreq's own
   `.venv`, not a separate one: the runner puts its interpreter's directory first on `PATH`
   for `pre_acceptance` and acceptance commands, which need black and flowfreq. Then install
   agy and sign in.
2. **Wire the hooks and status line** in `.claude/settings.local.json` (one machine; a cloud
   session has neither). The commands go through `tools/agentops.sh` (`statusline`,
   `hook session-start`, `hook pre-tool-use` on `Bash|TeamCreate|Workflow`). It picks
   `.venv/Scripts/python.exe` (Windows), then `.venv/bin/python`, then `python3`/`python`, so
   the Windows Store `python` shim is never used. Without the package installed the hook
   reports a non-blocking `No module named agentops` error and the session keeps working.
3. **agy permissions**, a different file from step 2. The runner passes
   `config/agy-settings.json` (`[agy] settings` in `config/agents.toml`) to every
   delegate call as `--settings`; agy then uses that file's `permissions` block for the call
   instead of its own `~/.gemini/antigravity-cli/settings.json` (its `cli.log` line
   `CLI settings initialized: permissions=...` shows which rules are live). So the repo governs
   what a delegate may run, and the host file needs no hand edits for runner calls; keep the
   host file in step with the example anyway for interactive `agy` use. Deny beats Ask beats
   Allow, as in ADR-0002 decision 3. The Claude Code hooks of step 2 do not belong in either
   file (agy ignores them and starts with `permissions=<nil>`). Allow file read/write inside the worktrees. Deny `git push`,
   `merge`, `rebase`, `checkout`, `switch`, `reset`, `branch`, `worktree`, `remote` and
   `config`; any command using `../`; and `write_file` on the main checkout, `.git/` and
   `vendor/`. `AGENTS.md` tells the delegate not to use the shell at all, but models try
   anyway, and **without an allow rule a headless `command` is auto-denied and ends the turn**:
   the envelope says `SUCCESS` with an empty response and the call still costs its tokens
   (2026-10-30: eight calls, 459k Gemini tokens, nothing written). So allow the read-only
   commands too (`python -m pytest|black|isort|mypy`, `git status|diff|log|show`, directory
   listings), as the sandbox's `D:/a/sbx/state/gemini/antigravity-cli/settings.json` does.
   Any denial ends the delegate's turn, so the allow list must cover what a model actually types:
   on 2026-10-30 the three attempts on one card were denied for bare `pytest tests/...` and for
   `git grep "..."`, neither of which `python -m pytest` and `git log|show` rules match. The example
   allows bare `pytest`, `python ...`, `git grep|ls-files|blame`, and the usual search and listing
   commands; the deny list (push, merge, data roots, `..`) still wins.
   Sonnet goes further and issues compound PowerShell lines (`git status --short; git diff;
   Select-String -Path ... -Pattern "a|b" | select -First 6`), which no single-command rule can
   match, so the example's first rule admits chains of whitelisted read-only commands joined by
   `;`, `&&` or `|`, with quoted arguments, and no redirection. The rule is case-insensitive
   (`(?i)`): PowerShell cmdlets are, and Sonnet's `| Select -First 6` cost one attempt on
   2026-10-30 against a rule that only knew `select`. Test a new rule against the
   command lines in `~/.gemini/antigravity-cli/conversations/*.db` before spending a card on it.
   Model names change: check `agy models` before writing `model =` on a card
   (`claude-sonnet-4-6` disappeared between 2026-10-26 and 2026-10-30).
4. Select the governor profile for your Claude plan: `run.bat gov profile pro` (or `max`; use
   `overnight` when you are away and only review; `gov profile` alone shows the active one and
   where it came from, `gov profile --clear` drops the choice). The choice is this machine's:
   it is saved as `governor_profile` in the gitignored state directory and overrides `profile`
   in `config/agents.toml`, which stays `pro` as the default for a cloud session.

## The loop

```bash
A="run.bat"                 # or: bash tools/agentops.sh (the same, from bash)
$A gov next                 # mode, queue state, and what to do next
$A cards add .agentops/cards/<id>.toml
$A gov plan                 # which ready cards would be dispatched now
$A cards run                # up to the governor's concurrency; or --card <id>
$A cards list
$A usage                    # Claude meter + agy tokens per pool
$A gov calibrate --agy-remaining 63   # after reading agy's /usage
$A cards reset <id>         # discard the worktree and branch, back to pending
$A cards show <id>          # a finished card: result, ledger, branch
$A cards diff <id> --stat   # its branch against the base
$A cards accept <id>        # after the PR merged: archive it
$A cards reject <id> --reason "..."   # back to pending with the reason
```

**Lead (Claude) per card:**
1. On a branch `spec/<id>` from `main`, write the failing tests and push them. Do not open a
   PR for this branch: it exists only to hold the tests. Cut it from current `main` and, if
   `main` moves before the card runs, merge `main` into it: the worktree is cut from this
   branch, so it must carry `AGENTS.md`. On 2026-10-30 seven spec branches predated that
   file; Sonnet, given no rules, wrote files through a PowerShell script in the shell, the
   host's agy permissions denied it, and three attempts escalated with a clean tree. The
   runner now copies the repo's `AGENTS.md` into a worktree whose branch lacks it (and
   drops the copy before the card commit), but that is a fallback, not the procedure.
2. Write `.agentops/cards/<id>.toml` with `base = "origin/spec/<id>"` (run `git fetch
   origin` before `cards run`, so the remote branch is current), an outcome-level `prompt`
   that names the files to read and edit, and an exact `acceptance` command. Use
   `python -m pytest <the card's test files> -q -p no:cacheprovider`. Do not use
   `PYTHONSAFEPATH=1 pytest`: the runner runs commands through `cmd.exe` on Windows,
   which cannot parse that syntax.
3. Once the runner reports `done`, review `card/<id>`. The `reviewer` pass covers scope,
   tests not weakened, the full suite, `make check` and the CLAUDE.md conventions. Then
   push the branch and open the PR to `main` (the cards in `.agentops/cards/` are the
   committed record). After the PR merges, `cards accept <id>` archives the card. `escalated`
   means the lead takes the card over or rewrites it.

**Why `child_env`:** the editable install points at the main checkout, so a worktree's
`python -m pytest` would import the main checkout's `flowfreq`. To prevent that, the runner
puts the card's worktree on `PYTHONPATH` (`{worktree}` in `config/agents.toml`) for both the
delegate and the acceptance run, and also sets `PYTHONUTF8=1`.

## Card format

```toml
id = "i214-snapshot-warning"           # [A-Za-z0-9._-]
repo = "flowfreq"
title = "Silence the nss_snapshot advisory for sweeps (#214)"
base = "origin/spec/i214-snapshot-warning"   # holds the failing tests; git fetch first
model = "claude-sonnet-5-5-medium"   # optional; omit to take the first free pool in [agy] models
depends_on = []   # orders dispatch only; recut this card's spec branch from main after they merge
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
