"""flowfreq's agentops config (config/agents.toml), checked without the agentops package.

The tool itself is the shared package in D:/a/agentops (docs/AGENTOPS.md), which CI does not
install, so this reads the committed config with tomllib and applies the push guards the way
the package's ``hooks.check`` does: ``re.search(pattern, command, re.IGNORECASE)`` over the
whole command line.
"""

from __future__ import annotations

import json
import re
import tomllib
from pathlib import Path
from typing import Any, Dict

import pytest

REPO = Path(__file__).resolve().parent.parent
CONFIG = REPO / "config" / "agents.toml"
PUSH = "git " + "push"


@pytest.fixture(scope="module")
def config() -> Dict[str, Any]:
    return tomllib.loads(CONFIG.read_text(encoding="utf-8"))


def denied(config: Dict[str, Any], command: str) -> bool:
    return any(re.search(p, command, re.IGNORECASE) for p in config["hooks"]["deny_bash"])


def test_committed_profile_is_pro(config: Dict[str, Any]) -> None:
    # A machine's own choice lives in the gitignored state dir, never in this file.
    assert config["governor"]["profile"] == "pro"
    assert set(config["governor"]["profiles"]) == {"pro", "max", "overnight"}


def test_delegate_settings(config: Dict[str, Any]) -> None:
    agy = config["agy"]
    assert "--dangerously-skip-permissions" not in agy["extra_args"]
    assert agy["child_env"]["PYTHONPATH"] == "{worktree}"
    settings = CONFIG.parent / agy["settings"]
    perms = json.loads(settings.read_text(encoding="utf-8"))["permissions"]
    assert perms["allow"]
    assert any(rule.startswith("command(regex:git (push|") for rule in perms["deny"])


def test_paths_point_at_this_repo(config: Dict[str, Any]) -> None:
    assert (CONFIG.parent / config["repos"]["flowfreq"]).resolve() == REPO
    assert (CONFIG.parent / config["paths"]["state_dir"]).resolve() == REPO / ".agentops" / "state"


@pytest.mark.parametrize(
    "command,blocked",
    [
        (f"{PUSH} -u origin feat/agentops", False),
        (f"{PUSH} origin docs/status-2026-10-26", False),
        (f"{PUSH} origin main", True),
        (f"{PUSH} origin HEAD:main", True),
        (f"{PUSH} origin release/0.16.x", True),
        (f"{PUSH} --force origin feat/x", True),
        (f"{PUSH} -f origin feat/x", True),
        (f"{PUSH} --force-with-lease origin feat/x", True),
        (f"{PUSH} origin --delete feat/x", True),
        (f"{PUSH} origin :refs/tags/v0.16.0", True),
        (f"{PUSH} origin --tags", True),
        ("git tag -d v0.16.0", True),
        ("git tag v0.17.0", False),
        ("git status", False),
        # A later command that merely ends in "main" is not a push to main.
        (f"{PUSH} origin feat && git switch main", False),
        (f"{PUSH} -u origin feat; git checkout main", False),
        (f"{PUSH} origin feat && rm -f x", False),
        (f"{PUSH} origin feat/main", False),
        (f"{PUSH} origin main-fix", False),
        (f"{PUSH} origin main:feat", False),
        (f"git tag v1 && {PUSH} origin v1", False),
        # Pushes to main, force pushes and deletes in the forms the old guards missed.
        (f"{PUSH} -u origin main", True),
        (f"{PUSH} origin main --dry-run", True),
        (f"{PUSH} origin feat:main", True),
        (f"{PUSH} origin HEAD:refs/heads/main", True),
        (f"{PUSH} origin +feat", True),
        (f"{PUSH} origin :feat", True),
        (f"{PUSH} -d origin feat", True),
        ("git -C D:/x push origin main", True),
        (f"git switch main && {PUSH} origin main", True),
        (f"{PUSH} origin feat && {PUSH} origin main", True),
    ],
)
def test_lead_push_guards(config: Dict[str, Any], command: str, blocked: bool) -> None:
    assert denied(config, command) is blocked
