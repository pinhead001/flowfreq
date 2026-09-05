"""Fortran EMA bridge via f2py-compiled peakfqr routines.

This package provides direct Python access to the USGS Fortran EMA
implementation (emafitpr) from the peakfqr R package, compiled via
numpy.f2py.

Usage::

    from flowfreq.peakfqr import emafitpr

The ``emafitpr`` function signature matches the Fortran subroutine
documented in ``vendor/peakfqr/src/emafit.f``.
"""

import os
import sys
from pathlib import Path

# Add DLL directory so Windows can find the bundled mingw runtime DLLs
_pkg_dir = Path(__file__).parent
if sys.platform == "win32" and hasattr(os, "add_dll_directory"):
    os.add_dll_directory(str(_pkg_dir))

try:
    from flowfreq.peakfqr._emafort import emafitpr  # noqa: E402
except ImportError as exc:  # pragma: no cover - depends on the build
    # The extension is gitignored and built on demand, so this is the normal
    # state for anyone who pip-installed flowfreq. The bare ImportError names
    # a private module nobody asked for and gives no hint what to do about it.
    raise ImportError(
        "flowfreq.peakfqr requires the f2py Fortran extension, which is built "
        "on demand rather than shipped. From a source checkout run:\n"
        "    python build_fortran/build.py    (or: make fortran)\n"
        "It needs gfortran and meson, plus the vendored sources under "
        "vendor/peakfqr/. Nothing else in flowfreq depends on this module -- "
        "the native Python EMA in flowfreq.bulletin17c is the default path, "
        "and flowfreq.validation.reference falls back to the committed golden "
        "files when the bridge is absent."
    ) from exc

__all__ = ["emafitpr"]
