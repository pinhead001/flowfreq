"""Build the f2py extension wrapping peakfq's Fortran EMA implementation.

Compiles ``vendor/peakfqr/src`` into ``_emafort`` via ``numpy.f2py``. Override
the source location with the ``PEAKFQR_SRC`` environment variable; on Windows,
override the MinGW toolchain location with ``MINGW_BIN``.

The signature file is not optional. ``_emafort.pyf`` restricts wrapping to
``emafitpr``; without it f2py tries to wrap every public symbol in the sources,
including QUADPACK's ``dqag``, which takes a user function as a callback. The
callback wrapper f2py generates for it does not compile against current NumPy
(``unknown type name 'f_t'`` in the generated ``_emafortmodule.c``), so the
build fails with an error that points at generated C rather than at anything
in this repository.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

BUILD_DIR = Path(__file__).resolve().parent
REPO_ROOT = BUILD_DIR.parent
SRC = Path(os.environ.get("PEAKFQR_SRC", REPO_ROOT / "vendor" / "peakfqr" / "src"))
PYF = BUILD_DIR / "_emafort.pyf"

# Order matters: emafit.f references symbols defined in the others.
SOURCE_NAMES = ("emafit.f", "dcdflib1.f90", "imslfake.f", "probfun.f")


def main() -> int:
    if not SRC.is_dir():
        sys.exit(
            f"Fortran sources not found at {SRC}.\n"
            "See docs/FORTRAN_UPLOAD.md -- vendor/peakfqr/src must be populated."
        )
    if not PYF.is_file():
        sys.exit(f"Signature file missing: {PYF}")

    sources = [SRC / name for name in SOURCE_NAMES]
    missing = [str(p) for p in sources if not p.is_file()]
    if missing:
        sys.exit("Missing Fortran sources:\n  " + "\n  ".join(missing))

    if sys.platform == "win32":
        mingw_bin = os.environ.get("MINGW_BIN", r"C:\msys64\mingw64\bin")
        if Path(mingw_bin).is_dir():
            os.environ["PATH"] = mingw_bin + os.pathsep + os.environ["PATH"]

    for tool in ("gfortran", "meson"):
        if shutil.which(tool) is None:
            sys.exit(f"ERROR: {tool} not found on PATH")

    # The .pyf comes first and supplies the module name, so no -m flag.
    cmd = [
        sys.executable,
        "-m",
        "numpy.f2py",
        "-c",
        str(PYF),
        *(str(p) for p in sources),
        "--backend",
        "meson",
        "--build-dir",
        str(BUILD_DIR / "mbuild"),
    ]
    print("Running:", " ".join(cmd))
    rc = subprocess.run(cmd, cwd=BUILD_DIR, env=os.environ).returncode
    if rc != 0:
        return rc

    # f2py drops the extension in the build directory, but flowfreq/peakfqr/__init__.py
    # imports it as flowfreq.peakfqr._emafort. Move it into the package so the bridge
    # actually loads. Both locations are gitignored; this is build output.
    built = sorted(BUILD_DIR.glob("_emafort*.so")) + sorted(BUILD_DIR.glob("_emafort*.pyd"))
    if not built:
        print("WARNING: build reported success but no extension was produced", file=sys.stderr)
        return 1
    dest_dir = REPO_ROOT / "flowfreq" / "peakfqr"
    for artifact in built:
        dest = dest_dir / artifact.name
        shutil.copy2(artifact, dest)
        print(f"Installed {dest.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
