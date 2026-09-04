"""Fixtures shared by the Fortran parity tests."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

GOLDEN_DIR = Path(__file__).parent / "golden"


def load_golden(name: str):
    """Load a golden document by case name, or None when it has not been generated.

    Returning None rather than raising lets a caller skip: a missing golden is
    a missing input, not a failure.
    """
    path = GOLDEN_DIR / f"{name}.json"
    if not path.is_file():
        return None
    return json.loads(path.read_text())


@pytest.fixture(scope="session")
def golden_big_sandy() -> dict:
    """The committed Fortran output for Big Sandy (USGS 03606500)."""
    path = GOLDEN_DIR / "big_sandy_03606500.json"
    if not path.is_file():
        pytest.skip(f"golden file missing: {path} (run tools/gen_fortran_golden.py)")
    return json.loads(path.read_text())


@pytest.fixture(scope="session")
def native_big_sandy():
    """flowfreq's own EMA fit of the same record."""
    from flowfreq.bulletin17c import Bulletin17C
    from tests.fixtures.big_sandy import (
        HISTORICAL_PEAKS,
        REGIONAL_SKEW,
        REGIONAL_SKEW_SD,
        SYSTEMATIC_PEAKS,
        THRESHOLDS,
    )

    analysis = Bulletin17C(
        peak_flows=list(SYSTEMATIC_PEAKS.values()),
        water_years=list(SYSTEMATIC_PEAKS.keys()),
        regional_skew=REGIONAL_SKEW,
        regional_skew_mse=REGIONAL_SKEW_SD**2,
        historical_peaks=[(y, f) for y, f in HISTORICAL_PEAKS.items()],
        perception_thresholds={(t["start"], t["end"]): t["lower"] for t in THRESHOLDS},
    )
    analysis.run_analysis()
    return analysis


def aep_index(golden: dict, aep: float) -> int:
    """Index of *aep* within the golden file's quantile vectors."""
    aeps = np.asarray(golden["inputs"]["aeps"], dtype=float)
    return int(np.argmin(np.abs(aeps - aep)))


def asymmetry(lower: float, point: float, upper: float) -> float:
    """Upper/lower half-width ratio of a confidence interval, in log10 space."""
    lo, mid, hi = np.log10([lower, point, upper])
    return (hi - mid) / (mid - lo)
