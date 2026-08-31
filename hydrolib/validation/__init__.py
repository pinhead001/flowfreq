"""
Validation module for HydroLib.

Provides comparison engines, benchmarks, and reporting for validating the
native Python EMA against the vendored USGS Fortran.

The reference side used to be :mod:`hydrolib.peakfqsa`, a subprocess wrapper
around a PeakfqSA executable that does not exist. It has been removed;
:class:`~hydrolib.validation.reference.ReferenceResult` takes its place and
reads from references that do exist -- a committed golden file, or a live call
through the f2py bridge in :mod:`hydrolib.peakfqr`.
"""

from hydrolib.validation.comparisons import ComparisonResult, FrequencyComparator
from hydrolib.validation.reference import ReferenceResult, cmoms_to_parameters

__all__ = [
    "ComparisonResult",
    "FrequencyComparator",
    "ReferenceResult",
    "cmoms_to_parameters",
]
