"""Build ``nss_snapshot`` equation libraries from NSS snapshots.

Each snapshot written by ``tools/snapshot_nss.py`` is turned into
``<out>/<CODE>.json`` by :func:`flowfreq.regression.nss_library.snapshot_to_library`.
A jurisdiction is skipped when it already has a published set
(``flowfreq/data/regression/<CODE>.json``), or when its snapshot has no regions
(NSS serves no peak-flow equations there).

Usage::

    python tools/build_nss_snapshot_library.py \
        --snapshot data/nss_snapshots/TN_2026-10-25.json \
        --out flowfreq/data/regression/nss_snapshot
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from flowfreq.regression.nss_library import snapshot_to_library  # noqa: E402

PUBLISHED_DIR = REPO_ROOT / "flowfreq" / "data" / "regression"
DEFAULT_OUT = PUBLISHED_DIR / "nss_snapshot"


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--snapshot", action="append", required=True, type=Path, help="snapshot file; repeat"
    )
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT, help=f"default {DEFAULT_OUT}")
    args = ap.parse_args(argv)

    for path in args.snapshot:
        snapshot = json.loads(path.read_text(encoding="utf-8"))
        code = str(snapshot["state"]).upper()
        if (PUBLISHED_DIR / f"{code}.json").exists():
            print(f"skipped {code}: published set exists")
            continue
        library = snapshot_to_library(snapshot)
        if library is None:
            print(f"skipped {code}: no NSS peak-flow regions")
            continue
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / f"{code}.json").write_text(
            json.dumps(library, indent=1) + "\n", encoding="utf-8"
        )
        print(f"{code}: {len(library['equations'])} equations, {len(library['unparsed'])} unparsed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
