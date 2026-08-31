"""Locate the peakfqr reference test data, wherever it happens to live.

The peakfqr R package's own test data (``results_WHIST.csv``, the
``wymt_ffa_2022A_*`` expected outputs) is the authoritative reference for
several fixtures here, but it is not yet part of this repository. Two
locations are searched, in order:

1. ``vendor/peakfqr/inst/testdata`` -- the in-repo home once the reference
   material is vendored (see ``docs/FORTRAN_UPLOAD.md``). Tests that depend
   on it start running automatically as soon as it appears.
2. ``../_shared/peakfqr/inst/testdata`` -- the reference tree beside the
   repository on the original development machine.

``TESTDATA_DIR`` is ``None`` when neither exists, and ``TESTDATA_AVAILABLE``
is then ``False``; dependent tests skip rather than fail, since a missing
reference tree is a missing input, not a defect.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional


def _repo_root() -> Path:
    """Walk up to the directory holding pyproject.toml.

    Counting ``parents[]`` broke silently once: this module moved from
    ``tests/peakfqsa/fixtures/`` to ``tests/fixtures/`` and the hardcoded
    ``parents[3]`` started pointing above the repository, which made
    TESTDATA_AVAILABLE False and skipped three passing tests rather than
    failing anything. Anchor on a file that marks the root instead.
    """
    for parent in Path(__file__).resolve().parents:
        if (parent / "pyproject.toml").is_file():
            return parent
    raise RuntimeError("could not locate the repository root from " + __file__)


REPO_ROOT = _repo_root()

_CANDIDATES = (
    REPO_ROOT / "vendor" / "peakfqr" / "inst" / "testdata",
    REPO_ROOT.parent / "_shared" / "peakfqr" / "inst" / "testdata",
)


def _resolve() -> Optional[Path]:
    for candidate in _CANDIDATES:
        if candidate.is_dir():
            return candidate
    return None


TESTDATA_DIR: Optional[Path] = _resolve()
TESTDATA_AVAILABLE: bool = TESTDATA_DIR is not None

SKIP_REASON = (
    "peakfqr reference test data not present; searched "
    + " and ".join(str(c) for c in _CANDIDATES)
    + " (see docs/FORTRAN_UPLOAD.md)"
)


def testdata_path(name: str) -> Path:
    """Return the full path to *name* inside the reference test data tree."""
    if TESTDATA_DIR is None:
        raise FileNotFoundError(SKIP_REASON)
    return TESTDATA_DIR / name
