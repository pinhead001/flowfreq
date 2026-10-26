"""Launcher for the `agentops` dev tool (tools/agentops/).

Puts this directory on sys.path explicitly, so it works under PYTHONSAFEPATH too, and
fails open for the hook and status-line entry points: a missing dev dependency must never
break a Claude Code session.

    python tools/agentops_cli.py --config agentops.toml gov status
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

if __name__ == "__main__":
    try:
        from agentops.cli import main
    except ImportError as e:
        if len(sys.argv) > 1 and {"hook", "statusline"} & set(sys.argv[1:]):
            print(f"agentops: unavailable ({e})", file=sys.stderr)
            raise SystemExit(0)
        raise
    raise SystemExit(main())
