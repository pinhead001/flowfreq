"""
Validation module for HydroLib.

Provides comparison engines, benchmarks, and reporting for validating
native Python EMA against reference implementations.
"""

from hydrolib.validation.comparisons import ComparisonResult, FrequencyComparator

__all__ = [
    "ComparisonResult",
    "FrequencyComparator",
]
