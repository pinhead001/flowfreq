# AGENTS.md — instructions for delegate coding agents (agy)

You implement exactly one task card in this worktree. The card's prompt and acceptance command define done.

- Make the acceptance command pass. Do not edit, skip or weaken existing tests.
- Do not run tests, builds or other long commands yourself, and never start background processes: the
  runner runs the acceptance command after you finish and sends you its output if it fails. Use your file
  read/write/edit tools only; do not use shell commands (they are denied in headless mode and end your turn
  with no output). Read files, edit files, then end your turn with a one-line summary.
- Read only the files the card names (plus this file). Do not browse the repository, list directories
  recursively, or open other modules for context: every file you read costs quota.
- Never read `vendor/`, `flowfreq/data/*.csv.gz` or anything under `.cache/`.
- No monkeypatching, import-time side effects or test-specific special cases in source code. If a test looks
  wrong (for example it fails for a reason unrelated to the card), stop and say so instead of working around it.
- Never launch subagents or research agents; read the named files yourself.
- Stay inside this worktree. Never push, merge, rebase or switch branches.
- Do not read or write outside this worktree except the files the card names. Use only the synthetic
  fixtures in `tests/fixtures/`. Never download data.
- Python 3.11-compatible syntax; the runner checks black, isort, mypy and pytest with `run.bat accept`.
  NumPy docstrings on public functions, type hints on every signature, `logging.getLogger(__name__)` not
  `print`, no bare `except`.
- If the card is ambiguous or impossible, stop and say why in your final response instead of guessing.
