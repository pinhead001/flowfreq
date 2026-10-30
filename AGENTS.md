# AGENTS.md — instructions for delegate coding agents (agy)

You implement exactly one task card in this worktree. The card's prompt and acceptance command define done.

- Make the acceptance command pass. Do not edit, skip, xfail or weaken existing tests; the card's tests are
  already written. Add tests only if the card asks.
- Do not run tests, builds or other long commands yourself, and never start background processes: the
  runner runs the acceptance command after you finish and sends you its output if it fails. Use your file
  read/write/edit tools only; do not use shell commands (they are denied in headless mode and end your turn
  with no output). Read files, edit files, then end your turn with a one-line summary.
- Never launch subagents or research agents; read the named files yourself.
- Read only the files the card names (plus this file and the parts of `CLAUDE.md` the card points to). Do not
  browse the repository, list directories recursively, or open other modules for context: every file you
  read costs quota.
- Stay inside this worktree. Never push, merge, rebase or switch branches.
- **Never edit anything under `vendor/`** (a verbatim USGS reference copy) or `tools/agentops/`,
  `agentops.toml`, `AGENTS.md`, `.claude/`, `.github/`.
- No monkeypatching, import-time side effects or test-specific special cases in library code. If a test
  looks wrong (for example it fails for a reason unrelated to the card), stop and say so instead of working
  around it.
- Python 3.11–3.14: no syntax or standard-library features newer than 3.11. No new dependencies.
- Style: black and isort (line length 100) run automatically after you finish; type hints on all
  signatures; NumPy-format docstrings on public API; no bare `except:`; no `print()` in `flowfreq/` (use
  `logging.getLogger(__name__)`).
- Notation: flow-duration quantiles are lowercase (q10, q95); flood quantiles keep the capital (Q100) or are
  given as AEP. Data are in log10 space when interfacing with Fortran conventions.
- Tests are offline: no network access, fixtures under `tests/fixtures/` only.
- Engineering judgment (regions, skews, tolerances, which equation applies) is the project owner's. If the
  card is ambiguous or impossible, stop and say why in your final response instead of guessing.
