# Task cards

Source files for the `agentops` delegate queue (`D:/a/agentops`, `docs/REQUIREMENTS.md` there).
A card is TOML with `id`, `repo = "flowfreq"`, `title`, `prompt`, `acceptance`, and optionally
`base` (default `main`), `model` (default: the first free model in `config/agents.toml`),
`print_timeout`, `max_retries`, `depends_on`. The prompt names every file the delegate may read
and edit (see `AGENTS.md`); the acceptance command is `run.bat accept <pytest paths>`, which
lints, type-checks and runs those tests in the card's worktree.

```
run.bat cards add config/cards/<id>.toml      # queue it
run.bat cards run --workers 1 --card <id>     # run one card (omit --card to run every ready card)
run.bat cards list                            # queue by state
run.bat ledger --tail 5                       # the delegate calls it cost
```

The runner works in `D:/a/sbx/state/agentops-flowfreq/worktrees/flowfreq-<id>` on branch
`card/<id>`, commits there when acceptance passes, and never pushes or merges. Results and the
acceptance output are in `D:/a/sbx/state/agentops-flowfreq/results/<id>/`. The lead reviews
with `git diff main..card/<id>`, pushes the branch and opens the PR. An `escalated` card is
finished by the lead in that worktree, not re-queued.
