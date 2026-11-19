"""Remove build and test artifacts, without a Unix shell.

``make clean`` and ``make clean-verify`` call this. They used ``rm``, ``find`` and
``xargs``, which a Windows ``cmd.exe`` does not have, and make is the supported route on
Windows. Standard library only.

    python tools/clean.py            # the Fortran build, the f2py extension, the caches
    python tools/clean.py --all      # also build/, dist/, *.egg-info and every __pycache__
    python tools/clean.py --dry-run  # say what would go, and remove nothing
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from typing import Iterable, List, Sequence

ROOT = Path(__file__).resolve().parent.parent
BUILD = (
    "build_fortran/mbuild",
    "build_fortran/native.ini",
    "build_fortran/_emafort*.so",
    "build_fortran/_emafort*.pyd",
    "flowfreq/peakfqr/_emafort*.so",
    "flowfreq/peakfqr/_emafort*.pyd",
    "flowfreq/peakfqr/*.dll",
    ".pytest_cache",
    ".mypy_cache",
)
EVERYTHING = ("build", "dist", "*.egg-info")
#: Never walked for ``__pycache__``: the vendored reference, the environments, the caches.
SKIPPED = {"vendor", ".venv", "venv", ".git", ".cache", "node_modules"}


def targets(root: Path, everything: bool) -> List[Path]:
    """The paths a clean removes, in a stable order."""
    found: List[Path] = []
    for pattern in (*BUILD, *(EVERYTHING if everything else ())):
        found.extend(sorted(root.glob(pattern)))
    if everything:
        found.extend(_pycaches(root))
    return found


def _pycaches(root: Path) -> Iterable[Path]:
    for child in sorted(root.iterdir()):
        if not child.is_dir() or child.name in SKIPPED:
            continue
        if child.name == "__pycache__":
            yield child
        else:
            yield from _pycaches(child)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--all", action="store_true", help="also build/, dist/, egg-info, __pycache__"
    )
    parser.add_argument("--dry-run", action="store_true", help="list what would be removed")
    args = parser.parse_args(argv)
    for path in targets(ROOT, args.all):
        print(("would remove " if args.dry_run else "removing ") + str(path.relative_to(ROOT)))
        if args.dry_run:
            continue
        if path.is_dir():
            shutil.rmtree(path, ignore_errors=True)
        else:
            path.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
