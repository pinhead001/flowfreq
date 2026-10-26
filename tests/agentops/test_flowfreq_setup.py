"""flowfreq's own agentops wiring: the committed config, child_env, usage report, launcher."""

import io
import json
import subprocess
import sys
from datetime import UTC, datetime, timedelta
from pathlib import Path

import pytest
from agentops import cooldown, governor, hooks, ledger, runner, usage
from agentops.config import AgentOpsConfig, load_config

from .test_runner import _card, repo_cfg  # noqa: F401  (fixture re-export)

REPO = Path(__file__).resolve().parents[2]
NOW = datetime(2026, 10, 26, 12, 0, tzinfo=UTC)


@pytest.fixture
def repo_toml() -> AgentOpsConfig:
    return load_config(REPO / "agentops.toml")


def test_committed_config_loads(repo_toml: AgentOpsConfig) -> None:
    assert repo_toml.repos == {"flowfreq": REPO}
    assert repo_toml.governor.profiles == ["pro", "max", "overnight"]
    assert repo_toml.agy.child_env["PYTHONPATH"] == "{worktree}"
    assert "--dangerously-skip-permissions" not in repo_toml.agy.extra_args
    assert cooldown.pool_of(repo_toml, "gemini-3.8-flash-low") == "gemini"
    assert cooldown.pool_of(repo_toml, "gpt-oss-120b-medium") == "third_party"


@pytest.mark.parametrize(
    "cmd,denied",
    [
        ("git push -u origin feat/agentops", False),
        ("git push origin docs/status-2026-10-26", False),
        ("git push origin main", True),
        ("git push origin HEAD:main", True),
        ("git push origin release/0.16.x", True),
        ("git push --force origin feat/x", True),
        ("git push -f origin feat/x", True),
        ("git push --force-with-lease origin feat/x", True),
        ("git push origin --delete feat/x", True),
        ("git push origin :refs/tags/v0.16.0", True),
        ("git push origin --tags", True),
        ("git tag -d v0.16.0", True),
        ("git tag v0.17.0", False),
        ("git status", False),
    ],
)
def test_lead_push_guards(repo_toml: AgentOpsConfig, cmd: str, denied: bool) -> None:
    assert (hooks.check(repo_toml, "Bash", {"command": cmd}, now=NOW) is not None) is denied


def test_overnight_profile_is_review_only(tmp_path: Path) -> None:
    text = (REPO / "agentops.toml").read_text(encoding="utf-8")
    text = text.replace('profile = "pro"', 'profile = "overnight"')
    text = text.replace('state_dir = ".agentops/state"', f'state_dir = "{tmp_path.as_posix()}"')
    p = tmp_path / "agentops.toml"
    p.write_text(text, encoding="utf-8")
    cfg = load_config(p)
    t = int(NOW.timestamp())
    cfg.claude_usage_path.write_text(
        json.dumps(
            {
                "five_hour": {"used_percentage": 1, "resets_at": t + 3600},
                "seven_day": {"used_percentage": 1, "resets_at": t + 6 * 86400},
            }
        )
    )
    assert governor.status(cfg, now=NOW).mode == "review-only"


def test_child_env_puts_worktree_on_pythonpath(repo_cfg: AgentOpsConfig) -> None:  # noqa: F811
    repo_cfg.agy.child_env = {"PYTHONPATH": "{worktree}", "AGENTOPS_TEST": "1"}
    check = (
        'python -c "import os,pathlib,sys; '
        "ok = os.path.samefile(os.environ['PYTHONPATH'], os.getcwd()) "
        "and os.environ['AGENTOPS_TEST'] == '1' "
        "and pathlib.Path('hello.txt').read_text() == 'hi'; sys.exit(0 if ok else 1)\""
    )
    res = runner.run_card(repo_cfg, _card(repo_cfg, acceptance=check))
    assert res.state == "done", res.detail


def test_usage_summarizes_pools(cfg: AgentOpsConfig) -> None:
    def rec(hours_ago: float, model: str, tok: int, ok: bool = True, card: str = "c1") -> None:
        ts = (NOW - timedelta(hours=hours_ago)).isoformat()
        ledger.append(
            cfg.ledger_path,
            {"ts": ts, "kind": "agy", "card": card, "model": model, "ok": ok, "total_tokens": tok},
        )

    rec(1, "m1", 100)
    rec(30, "m1", 200, ok=False)
    rec(2, "m2", 50, card="c2")
    rec(24 * 8, "m1", 9999)  # outside the week
    ledger.append(
        cfg.ledger_path,
        {"ts": NOW.isoformat(), "kind": "quota", "card": "c1", "model": "m1", "pool": "a",
         "until": (NOW + timedelta(hours=2)).isoformat()},
    )  # fmt: skip
    cooldown.set_cooldown(cfg, "a", NOW + timedelta(hours=2))
    u = usage.summarize(cfg, now=NOW)
    pools = {p.pool: p for p in u.pools}
    assert (pools["a"].calls_24h, pools["a"].tokens_24h) == (1, 100)
    assert (pools["a"].calls_7d, pools["a"].tokens_7d, pools["a"].failed_7d) == (2, 300, 1)
    assert pools["a"].cooling_until is not None and pools["b"].tokens_7d == 50
    assert u.by_card_7d == {"c1": 300, "c2": 50}
    assert len(u.quota_events_7d) == 1
    text = usage.render(u)
    assert "quota events" in text and "c1 300" in text


def test_usage_cli_json(cfg: AgentOpsConfig, capsys: pytest.CaptureFixture[str]) -> None:
    from agentops.cli import main

    assert main(["--config", str(cfg.source), "usage", "--json"]) == 0
    d = json.loads(capsys.readouterr().out)
    assert [p["pool"] for p in d["pools"]] == ["a", "b"]


def test_launcher_fails_open_for_hooks(tmp_path: Path) -> None:
    launcher = REPO / "tools" / "agentops_cli.py"
    p = subprocess.run(
        [sys.executable, str(launcher), "--config", str(tmp_path / "missing.toml"), "hook",
         "session-start"],
        input="{}", capture_output=True, text=True, cwd=tmp_path,
    )  # fmt: skip
    assert p.returncode == 0 and "hook skipped" in p.stderr


def test_session_start_hook_reports_governor(repo_toml: AgentOpsConfig, tmp_path: Path) -> None:
    repo_toml.state_dir = tmp_path
    out = io.StringIO()
    hooks.session_start(repo_toml, io.StringIO("{}"), out, now=NOW)
    ctx = json.loads(out.getvalue())["hookSpecificOutput"]["additionalContext"]
    assert ctx.startswith("governor: conserve (no Claude meter data)")
