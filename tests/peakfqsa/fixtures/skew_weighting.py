"""
Test fixtures from peakfqr test-skewweight.R (Halloween skew weighting).

Tests the determinant ratio weighting factor (detrat) against Greg Schwarz's
independent SAS implementation. Data from results_WHIST.csv.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from typing import Optional

# Path to the WHIST CSV (shared test data from peakfqr package)
_WHIST_CSV = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "_shared",
    "peakfqr",
    "inst",
    "testdata",
    "results_WHIST.csv",
)


@dataclass
class SkewWeightCase:
    """A single row from the WHIST skew weighting test data."""

    nu: float
    sigma2: float
    freq: float
    tau: float
    skew: float
    probut: float
    herl_factor: float
    jerl_factor: float


def load_whist_cases() -> list[SkewWeightCase]:
    """Load all WHIST test cases from the CSV file.

    Returns
    -------
    list[SkewWeightCase]
        All test cases with inputs and expected HERL_factor values.
    """
    csv_path = os.path.normpath(_WHIST_CSV)
    cases: list[SkewWeightCase] = []
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cases.append(
                SkewWeightCase(
                    nu=float(row["NU"]),
                    sigma2=float(row["SIGMA2"]),
                    freq=float(row["FREQ"]),
                    tau=float(row["TAU"]),
                    skew=float(row["SKEW"]),
                    probut=float(row["PROBUT"]),
                    herl_factor=float(row["HERL_factor"]),
                    jerl_factor=float(row["JERL_factor"]),
                )
            )
    return cases


# Pre-computed inputs for the detrat function, matching the R test:
# n = 100 (total observations)
# nobs = [n*(1-NU), n*NU]  (censored, uncensored counts)
# tl = [SIGMA*TAU + MU, -20]  (lower perception thresholds)
# tu = [20, 20]  (upper perception thresholds)
# m = MU = 2 (chosen mean)
# sd = sqrt(SIGMA2)
# g = SKEW
# nthresh = 2
DETRAT_N: int = 100
DETRAT_MU: float = 2.0
DETRAT_NTHRESH: int = 2
TOLERANCE: float = 0.001
