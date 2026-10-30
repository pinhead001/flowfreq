import io
import json
from datetime import UTC, datetime

from agentops import governor, ledger, statusline
from agentops.config import AgentOpsConfig

NOW = datetime(2026, 10, 23, 12, 0, tzinfo=UTC)
T = int(NOW.timestamp())


def test_config_resolves_relative_paths(cfg: AgentOpsConfig) -> None:
    assert cfg.state_dir.is_absolute() and cfg.state_dir.name == "state"
    assert cfg.repos["lib"].is_absolute()
    assert cfg.agy.model == "m1"


def test_statusline_writes_usage(cfg: AgentOpsConfig) -> None:
    payload = {
        "rate_limits": {
            "five_hour": {"used_percentage": 23.5, "resets_at": T + 3600},
            "seven_day": {"used_percentage": 41.2},
        }
    }
    out = io.StringIO()
    statusline.run(cfg, io.StringIO(json.dumps(payload)), out, now=NOW)
    data = json.loads(cfg.claude_usage_path.read_text())
    assert data["five_hour"] == {"used_percentage": 23.5, "resets_at": T + 3600}
    assert data["seven_day"] == {"used_percentage": 41.2, "resets_at": "unknown"}
    assert "5h 24%" in out.getvalue()


def test_statusline_survives_garbage(cfg: AgentOpsConfig) -> None:
    out = io.StringIO()
    statusline.run(cfg, io.StringIO("not json"), out, now=NOW)
    data = json.loads(cfg.claude_usage_path.read_text())
    assert data["five_hour"]["used_percentage"] == "unknown"


def _usage(
    cfg: AgentOpsConfig,
    five: float | str,
    seven: float | str,
    seven_reset: int | str = T + 3 * 86400,
) -> None:
    cfg.claude_usage_path.parent.mkdir(parents=True, exist_ok=True)
    cfg.claude_usage_path.write_text(
        json.dumps(
            {
                "five_hour": {"used_percentage": five, "resets_at": T + 3600},
                "seven_day": {"used_percentage": seven, "resets_at": seven_reset},
            }
        )
    )


def test_mode_thresholds(cfg: AgentOpsConfig) -> None:
    # 3 days left of 7 -> elapsed 4/7 = 57%; pace limit 67%
    for five, mode in [
        (10, "normal"),
        (40.1, "conserve"),
        (74, "conserve"),
        (80, "review-only"),
        (100, "wait"),
    ]:
        _usage(cfg, five, 30)
        assert governor.status(cfg, now=NOW).mode == mode, five


def test_weekly_pace_demotes_normal(cfg: AgentOpsConfig) -> None:
    _usage(cfg, 10, 70)
    s = governor.status(cfg, now=NOW)
    assert s.mode == "conserve" and not s.weekly_pace_ok


def test_unknown_usage_uses_unknown_mode(cfg: AgentOpsConfig) -> None:
    assert governor.status(cfg, now=NOW).mode == "conserve"
    _usage(cfg, "unknown", "unknown")
    assert governor.status(cfg, now=NOW).mode == "conserve"


def test_expired_window_counts_as_reset(cfg: AgentOpsConfig) -> None:
    cfg.claude_usage_path.parent.mkdir(parents=True, exist_ok=True)
    cfg.claude_usage_path.write_text(
        json.dumps(
            {
                "five_hour": {"used_percentage": 100, "resets_at": T - 10},
                "seven_day": {"used_percentage": 30, "resets_at": T + 3 * 86400},
            }
        )
    )
    assert governor.status(cfg, now=NOW).mode == "normal"


def test_agy_concurrency_halves_on_errors_and_stops_low(cfg: AgentOpsConfig) -> None:
    _usage(cfg, 10, 30)
    assert governor.status(cfg, now=NOW).agy_concurrency == 2
    ledger.append(
        cfg.ledger_path, {"ts": NOW.isoformat(), "kind": "agy", "ok": False, "total_tokens": 5}
    )
    assert governor.status(cfg, now=NOW).agy_concurrency == 1
    governor.calibrate(cfg, agy_remaining_pct=5, now=NOW)
    s = governor.status(cfg, now=NOW)
    assert s.agy_concurrency == 0 and s.agy_remaining_pct == 5


def test_agy_estimate_uses_tokens_since_calibration(cfg: AgentOpsConfig) -> None:
    governor.calibrate(cfg, agy_remaining_pct=50, now=NOW)
    later = NOW.replace(hour=13)
    ledger.append(
        cfg.ledger_path, {"ts": later.isoformat(), "kind": "agy", "ok": True, "total_tokens": 3000}
    )
    cfg.agy.tokens_per_pct = 100
    assert governor.status(cfg, now=later).agy_remaining_pct == 20


def test_statusline_shows_time_to_five_hour_reset(cfg: AgentOpsConfig) -> None:
    payload = {
        "rate_limits": {
            "five_hour": {"used_percentage": 23.5, "resets_at": T + 3600 + 600},
            "seven_day": {"used_percentage": 41.2, "resets_at": T + 86400},
        }
    }
    out = io.StringIO()
    statusline.run(cfg, io.StringIO(json.dumps(payload)), out, now=NOW)
    assert "5h 24% (reset 1h10m)" in out.getvalue()


def test_statusline_omits_agy_when_idle(cfg: AgentOpsConfig) -> None:
    out = io.StringIO()
    statusline.run(cfg, io.StringIO("{}"), out, now=NOW)
    assert "agy" not in out.getvalue()


def _queue_card(cfg: AgentOpsConfig, state: str, cid: str) -> None:
    d = cfg.queue_dir / state
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{cid}.toml").write_text(
        f'id = "{cid}"\nrepo = "lib"\ntitle = "t"\nprompt = "p"\nacceptance = "true"\n',
        encoding="utf-8",
    )


def test_statusline_shows_agy_queue_and_pool_tokens(cfg: AgentOpsConfig) -> None:
    _queue_card(cfg, "running", "c1")
    _queue_card(cfg, "running", "c2")
    _queue_card(cfg, "pending", "c3")
    ledger.append(
        cfg.ledger_path,
        {"ts": NOW.isoformat(), "kind": "agy", "card": "c1", "model": "m1", "ok": False,
         "total_tokens": 45507},
    )  # fmt: skip
    out = io.StringIO()
    statusline.run(cfg, io.StringIO("{}"), out, now=NOW)
    text = out.getvalue()
    assert "agy 2 run, 1 pend" in text
    assert "7d tok a 46k" in text  # 45,507 tokens on model m1, pool a


def test_render_shows_time_to_five_hour_reset(cfg: AgentOpsConfig) -> None:
    _usage(cfg, 23.5, 41.2)
    text = governor.render(governor.status(cfg, now=NOW))
    assert "(in 1h00m)" in text
    assert "5h used 24%" in text
