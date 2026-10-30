"""Usage report: Claude meter plus delegate tokens per quota pool, from the ledger.

Delegate quota has no headless readout (flowfield2d ADR-0002), and the observed limit is a
weekly window per pool, so the useful view is ledger tokens per pool over the last day and
the last seven days, the calls that failed, and the quota events with their reset times.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime, timedelta
from typing import Any

from agentops import cooldown, governor, ledger
from agentops.config import AgentOpsConfig


@dataclass
class PoolUsage:
    pool: str
    calls_24h: int = 0
    tokens_24h: int = 0
    calls_7d: int = 0
    tokens_7d: int = 0
    failed_7d: int = 0
    cooling_until: str | None = None


@dataclass
class Usage:
    governor: dict[str, Any]
    pools: list[PoolUsage] = field(default_factory=list)
    quota_events_7d: list[dict[str, Any]] = field(default_factory=list)
    by_card_7d: dict[str, int] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def summarize(cfg: AgentOpsConfig, now: datetime | None = None) -> Usage:
    now = now or datetime.now(UTC)
    recs = ledger.read(cfg.ledger_path)
    week = ledger.since(recs, now - timedelta(days=7))
    day_start = now - timedelta(days=1)
    cooling = cooldown.active(cfg, now)

    names = list(cfg.agy.pools) + [
        p for p in (cooldown.pool_of(cfg, m) for m in cfg.agy.models) if p not in cfg.agy.pools
    ]
    pools: dict[str, PoolUsage] = {n: PoolUsage(n) for n in names}
    by_card: dict[str, int] = {}
    for r in week:
        if r.get("kind") != "agy":
            continue
        name = cooldown.pool_of(cfg, str(r.get("model") or ""))
        pu = pools.setdefault(name, PoolUsage(name))
        tok = int(r.get("total_tokens") or 0)
        pu.calls_7d += 1
        pu.tokens_7d += tok
        pu.failed_7d += 0 if r.get("ok") else 1
        if datetime.fromisoformat(str(r["ts"])) >= day_start:
            pu.calls_24h += 1
            pu.tokens_24h += tok
        card = str(r.get("card") or "?")
        by_card[card] = by_card.get(card, 0) + tok
    for name, until in cooling.items():
        pools.setdefault(name, PoolUsage(name)).cooling_until = until.isoformat()

    quota = [
        {k: r.get(k) for k in ("ts", "card", "model", "pool", "until")}
        for r in week
        if r.get("kind") == "quota"
    ]
    return Usage(
        governor=governor.status(cfg, now=now).as_dict(),
        pools=list(pools.values()),
        quota_events_7d=quota,
        by_card_7d=dict(sorted(by_card.items(), key=lambda kv: -kv[1])),
    )


def render(u: Usage) -> str:
    g = u.governor

    def pct(v: Any) -> str:
        return "?" if v is None else f"{v:.0f}%"

    lines = [
        f"claude: 5h used {pct(g['five_hour_used'])}, 7d used {pct(g['seven_day_used'])}, "
        f"governor {g['mode']} ({g['reason']})",
        f"{'pool':<16}{'24h calls':>10}{'24h tok':>12}{'7d calls':>10}{'7d tok':>12}"
        f"{'7d failed':>10}  cooling",
    ]
    for p in u.pools:
        lines.append(
            f"{p.pool:<16}{p.calls_24h:>10}{p.tokens_24h:>12,}{p.calls_7d:>10}"
            f"{p.tokens_7d:>12,}{p.failed_7d:>10}  {(p.cooling_until or '-')[:16]}"
        )
    if u.quota_events_7d:
        lines.append("quota events (7d):")
        lines += [
            f"  {e['ts'][:16]} {e['pool']} ({e['model']}) on {e['card']}, until {str(e['until'])[:16]}"
            for e in u.quota_events_7d
        ]
    if u.by_card_7d:
        top = list(u.by_card_7d.items())[:5]
        lines.append("top cards (7d tokens): " + ", ".join(f"{c} {t:,}" for c, t in top))
    return "\n".join(lines)
