"""Agent ops: Claude meter, delegate runner, governor and hooks.

Ported from flowfield2d's ``flowfield2d.agentops`` (the same design, see its ADR-0002) as a
dev-only tool: it is not part of the published ``flowfreq`` package. All paths and thresholds
come from the config file (``agentops.toml`` at the repo root).
"""
