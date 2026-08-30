"""Generate golden files recording the vendored Fortran's output.

The f2py extension builds only where gfortran and ``vendor/peakfqr/src`` are
present, but parity tests must run in CI where it is not. This script captures
``emafitpr``'s full output once and commits it as JSON, so
``tests/fortran_parity/test_native_vs_golden.py`` can assert against the
reference implementation anywhere.

Usage::

    python build_fortran/build.py          # produces hydrolib/peakfqr/_emafort*
    python tools/gen_fortran_golden.py     # writes tests/fortran_parity/golden/*.json

Regenerate whenever ``vendor/peakfqr`` changes, and commit the result: a golden
file and the sources that produced it must move together, or the parity tests
silently start comparing against a different reference.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

GOLDEN_DIR = REPO_ROOT / "tests" / "fortran_parity" / "golden"
DESCRIPTION = REPO_ROOT / "vendor" / "peakfqr" / "DESCRIPTION"


def peakfq_version() -> str:
    """Read the vendored package version, so a golden file names its source."""
    if not DESCRIPTION.is_file():
        return "unknown"
    match = re.search(r"^Version:\s*(\S+)", DESCRIPTION.read_text(), re.MULTILINE)
    return match.group(1) if match else "unknown"


def source_fingerprint() -> str:
    """Short git description of the vendored sources, for traceability."""
    try:
        out = subprocess.run(
            ["git", "log", "-1", "--format=%h", "--", "vendor/peakfqr/src"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        return out.stdout.strip() or "unknown"
    except (subprocess.CalledProcessError, OSError):
        return "unknown"


def generate(name: str) -> Path:
    from tests.fortran_parity.cases import CASES, build_emafit_inputs, call_emafitpr

    case = CASES[name]()
    inputs = build_emafit_inputs(case)
    outputs = call_emafitpr(case)

    document = {
        "meta": {
            "site_no": case.site_no,
            "description": case.description,
            "peakfq_version": peakfq_version(),
            "vendor_src_commit": source_fingerprint(),
            "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "generator": "tools/gen_fortran_golden.py",
            "units": "all flow quantities are log10(cfs), exactly as emafitpr returns them",
        },
        "inputs": {
            "n": len(inputs["ql"]),
            "aeps": list(case.aeps),
            "eps": case.eps,
            "weight_opt": case.weight_opt,
            "gbthrsh0": case.gbthrsh0,
            "regional_skew": case.regional_skew,
            "regional_skew_mse": case.regional_skew_mse,
            "reg_m": inputs["reg_m"],
            "reg_m_mse": inputs["reg_m_mse"],
            "reg_sd": inputs["reg_sd"],
            "reg_sd_mse": inputs["reg_sd_mse"],
            "ql": inputs["ql"].tolist(),
            "qu": inputs["qu"].tolist(),
            "tl": inputs["tl"].tolist(),
            "tu": inputs["tu"].tolist(),
            "dtype": inputs["dtype"].tolist(),
        },
        "outputs": outputs,
    }

    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    path = GOLDEN_DIR / f"{name}.json"
    path.write_text(json.dumps(document, indent=2, sort_keys=False) + "\n")
    return path


def main() -> int:
    from tests.fortran_parity.cases import CASES

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "case", nargs="*", default=sorted(CASES), help="case names to regenerate (default: all)"
    )
    args = parser.parse_args()

    try:
        import hydrolib.peakfqr  # noqa: F401
    except ImportError as exc:
        print(
            f"Fortran extension not importable: {exc}\n"
            "Run `python build_fortran/build.py` first "
            "(see docs/FORTRAN_UPLOAD.md).",
            file=sys.stderr,
        )
        return 1

    for name in args.case:
        if name not in CASES:
            print(f"Unknown case {name!r}; known: {', '.join(sorted(CASES))}", file=sys.stderr)
            return 2
        path = generate(name)
        print(f"wrote {path.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
