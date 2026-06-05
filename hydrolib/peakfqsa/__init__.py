"""
PeakfqSA integration module for HydroLib.

Provides subprocess wrapping, I/O conversion, and output parsing for the
USGS PeakfqSA Fortran reference implementation of Bulletin 17C flood
frequency analysis.
"""

from hydrolib.peakfqsa.config import PeakfqSAConfig, PeakfqSANotFoundError, find_peakfqsa
from hydrolib.peakfqsa.io_converters import DataFile, SpecificationFile
from hydrolib.peakfqsa.parsers import PeakfqSAResult
from hydrolib.peakfqsa.wrapper import (
    PeakfqSAExecutionError,
    PeakfqSAParseError,
    PeakfqSATimeoutError,
    PeakfqSAWrapper,
)

__all__ = [
    "PeakfqSAConfig",
    "PeakfqSANotFoundError",
    "find_peakfqsa",
    "DataFile",
    "SpecificationFile",
    "PeakfqSAResult",
    "PeakfqSAWrapper",
    "PeakfqSAExecutionError",
    "PeakfqSAParseError",
    "PeakfqSATimeoutError",
]
