"""Fortran EMA bridge via f2py-compiled peakfqr routines.

This package provides direct Python access to the USGS Fortran EMA
implementation (emafitpr) from the peakfqr R package, compiled via
numpy.f2py.

Usage::

    from hydrolib.peakfqr import emafitpr

The ``emafitpr`` function signature matches the Fortran subroutine
documented in ``_shared/peakfqr/src/emafit.f``.
"""

import os
import sys
from pathlib import Path

# Add DLL directory so Windows can find the bundled mingw runtime DLLs
_pkg_dir = Path(__file__).parent
if sys.platform == "win32" and hasattr(os, "add_dll_directory"):
    os.add_dll_directory(str(_pkg_dir))

from hydrolib.peakfqr._emafort import emafitpr  # noqa: E402

__all__ = ["emafitpr"]
