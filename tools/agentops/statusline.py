"""Claude Code status-line command: records the rate-limit meter, prints one line.

Reads the status-line JSON on stdin and writes `claude_usage.json` with
`rate_limits.five_hour` and `seven_day` (`used_percentage`, `resets_at` in
epoch seconds). Fields Claude Code leaves out are recorded as "unknown".
Never raises: a broken status line must not break the session.

The line is `5h <used>% (reset <in>) | 7d <used>% | gov <mode>`, plus, while
agentops is active (cards running or pending, or delegate calls in the last
24 h), `agy <n> run, <m> pend | 7d tok <pool> <k>...`.
"""

from __future__ import annotations

import json
from datetime import UTC, datetime
from typing import Any, TextIO

from agentops.config import AgentOpsConfig
from agentops.ledger import write_json_atomic

UNKNOWN = "unknown"


def _window(rl: Any, name: str) -> dict[str, Any]:
    w = rl.get(name) if isinstance(rl, dict) else None
    w = w if isinstance(w, dict) else {}
    out: dict[str, Any] = {}
    for key in ("used_percentage", "resets_at"):
        v = w.get(key)
        out[key] = v if isinstance(v, int | float) and not isinstance(v, bool) else UNKNOWN
    return out


def run(cfg: AgentOpsConfig, stdin: TextIO, stdout: TextIO, now: datetime | None = None) -> None:
    now = now or datetime.now(UTC)
    try:
        payload = json.loads(stdin.read() or "{}")
    except json.JSONDecodeError:
        payload = {}
    rl = payload.get("rate_limits") if isinstance(payload, dict) else None
    usage: dict[str, Any] = {
        "updated_at": now.isoformat(),
        "five_hour": _window(rl, "five_hour"),
        "seven_day": _window(rl, "seven_day"),
    }
    try:
        write_json_atomic(cfg.claude_usage_path, usage)
    except OSError:
        pass
    parts = []
    for label, key in (("5h", "five_hour"), ("7d", "seven_day")):
        v = usage[key]["used_percentage"]
        text = f"{label} {v:.0f}%" if v != UNKNOWN else f"{label} ?"
        reset = usage[key]["resets_at"]
        if label == "5h" and reset != UNKNOWN:
            from agentops.governor import fmt_duration

            text += f" (reset {fmt_duration(max(0, int(reset - now.timestamp())))})"
        parts.append(text)
    try:
        from agentops import governor

        parts.append(f"gov {governor.status(cfg, now=now).mode}")
    except Exception:  # noqa: BLE001 - status line must not fail
        parts.append("gov ?")
    try:
        parts.extend(_agy_parts(cfg, now))
    except Exception:  # noqa: BLE001 - status line must not fail
        pass
    stdout.write(" | ".join(parts) + "\n")


def _k(tokens: int) -> str:
    return f"{tokens / 1000:.0f}k" if tokens < 1_000_000 else f"{tokens / 1e6:.1f}M"


def _agy_parts(cfg: AgentOpsConfig, now: datetime) -> list[str]:
    """The delegate side, only while agentops is active: queue counts and 7-day pool tokens."""
    from agentops import usage
    from agentops.queue import Queue

    q = Queue(cfg)
    running, pending = len(q.cards("running")), len(q.cards("pending"))
    u = usage.summarize(cfg, now=now)
    if not (running or pending or any(p.calls_24h for p in u.pools)):
        return []
    parts = [f"agy {running} run, {pending} pend"]
    toks = [f"{p.pool} {_k(p.tokens_7d)}" for p in u.pools if p.tokens_7d]
    if toks:
        parts.append("7d tok " + ", ".join(toks))
    return parts
