"""Fixtures for tools/agentops (ported from flowfield2d tests/unit/agentops)."""

import sys
from pathlib import Path

import pytest

TOOLS = Path(__file__).resolve().parents[2] / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from agentops.config import AgentOpsConfig, load_config  # noqa: E402

CONFIG = """
[paths]
state_dir = "state"
worktree_root = "wt"
worktree_prefix = "wt-"

[repos]
lib = "repos/lib"

[governor]
normal_min_remaining = 60
conserve_min_remaining = 25
weekly_pace_slack = 10
unknown_mode = "conserve"

[agy]
exe = "agy"
models = ["m1", "m2"]
pools = { a = ["m1"], b = ["m2"] }
default_cooldown_min = 60
print_timeout = "1m"
extra_args = []
max_retries = 1
concurrency_start = 2
error_window_min = 60
stop_below_pct = 10
tokens_per_pct = 0

[hooks]
deny_bash = ['\\bgit\\s+push\\b', '\\bnetsh\\b']
deny_substrings = ["_secret"]
teammate_tools = ["TeamCreate"]
task_acceptance = "python -c \\"import sys; sys.exit(0)\\""
"""


@pytest.fixture
def cfg(tmp_path: Path) -> AgentOpsConfig:
    p = tmp_path / "agents.toml"
    p.write_text(CONFIG, encoding="utf-8")
    return load_config(p)
