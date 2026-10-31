"""Build ``nss_snapshot`` equation libraries from NSS snapshots.

Each snapshot written by ``tools/snapshot_nss.py`` is turned into
``<out>/<CODE>.json`` by :func:`flowfreq.regression.nss_library.snapshot_to_library`.
A jurisdiction is skipped when it already has a published set
(``<published-dir>/<CODE>.json``), when its snapshot has no regions (NSS
serves no peak-flow equations there), or when none of its regions yields an equation
(every region unresolved or refused): a library with no equations is not written.

A snapshot that cannot be read (invalid JSON, no ``state``), or a second snapshot
for a jurisdiction already given, is reported with its path and skipped; the rest
are still built, and the exit code is 1.

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
from typing import Any, Dict, Optional, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from flowfreq.regression.nss_library import merge_snapshots, snapshot_to_library  # noqa: E402
from flowfreq.streamstats import NSS_EXTRA_REGIONS  # noqa: E402

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
    ap.add_argument(
        "--published-dir",
        type=Path,
        default=PUBLISHED_DIR,
        help=f"published equation sets, never overridden (default {PUBLISHED_DIR})",
    )
    args = ap.parse_args(argv)

    status = 0
    seen: Dict[str, Path] = {}
    snapshots: Dict[str, Dict[str, Any]] = {}
    for path in args.snapshot:
        try:
            snapshot = json.loads(path.read_text(encoding="utf-8"))
            code = str(snapshot["state"]).upper()
        except (OSError, ValueError, KeyError, TypeError) as exc:
            print(f"error {path}: cannot read snapshot ({type(exc).__name__}: {exc})")
            status = 1
            continue
        if code in seen:
            print(f"error {code}: two snapshots given, {seen[code]} and {path}; skipped {path}")
            status = 1
            continue
        seen[code] = path
        snapshots[code] = snapshot

    # An NSS region that holds a study inside a jurisdiction (RRB in MN) is merged into
    # that jurisdiction's snapshot, never written on its own.
    owner = {extra: state for state, extras in NSS_EXTRA_REGIONS.items() for extra in extras}
    for extra in [c for c in snapshots if c in owner]:
        state = owner[extra]
        snap = snapshots.pop(extra)
        if state not in snapshots:
            print(
                f"error {extra}: an extra NSS region of {state}, given without {state}'s snapshot"
            )
            status = 1
            continue
        snapshots[state] = merge_snapshots(snapshots[state], [snap])
        print(f"{extra} merged into {state}")

    for code, snapshot in snapshots.items():
        if (args.published_dir / f"{code}.json").exists():
            print(f"skipped {code}: published set exists")
            continue
        library = snapshot_to_library(snapshot)
        if library is None:
            print(f"skipped {code}: no NSS peak-flow regions")
            continue
        if not library["equations"]:
            print(f"skipped {code}: no parseable equations")
            continue
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / f"{code}.json").write_text(
            json.dumps(library, indent=1) + "\n", encoding="utf-8"
        )
        print(f"{code}: {len(library['equations'])} equations, {len(library['unparsed'])} unparsed")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
