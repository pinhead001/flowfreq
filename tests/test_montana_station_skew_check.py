"""flowfreq station skews vs USGS PeakFQ 7.5.1 on the SIR 2025-5019 examples.

Reads data/skew_study/sir2025_5019_station_skews.csv, written by
tools/check_montana_station_skews.py from the report's data release
(https://doi.org/10.5066/P1WHRK8H): the .psf specifications and peaks are fitted
natively and compared with Table 1-5's PeakFQ results.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

CSV = Path(__file__).resolve().parents[1] / "data" / "skew_study" / "sir2025_5019_station_skews.csv"

#: Analyses where flowfreq and PeakFQ 7.5.1 disagree by more than 0.2. Empty since
#: #135: 06163400.00 differed (3.0 against 0.216) because flowfreq upper-cased its
#: lowercase "c" peak codes into code C and dropped 39 peaks; siteQT matches "C"
#: case-sensitively and keeps them.
KNOWN_DISAGREEMENTS: dict = {}


@pytest.fixture(scope="module")
def checks() -> pd.DataFrame:
    return pd.read_csv(CSV, dtype={"analysis_id": str, "site_no": str})


def test_every_analysis_compared(checks):
    assert len(checks) == 63
    assert checks["usgs_station_skew"].notna().all()
    assert checks["flowfreq_station_skew"].notna().all()


def test_station_skews_match_peakfq(checks):
    diff = checks["skew_diff"].abs()
    assert diff.median() < 0.001
    assert set(checks.loc[diff > 0.2, "analysis_id"]) == set(KNOWN_DISAGREEMENTS)


def test_montana_basin_analyses_match_peakfq(checks):
    mt = checks[checks["basin"].isin(["Missouri (10)", "Columbia (17)"])]
    mt = mt[~mt["analysis_id"].isin(KNOWN_DISAGREEMENTS)]
    assert len(mt) >= 40
    assert mt["skew_diff"].abs().max() < 0.1
    assert (mt["skew_diff"].abs() < 0.01).mean() > 0.95


def test_lowercase_c_records_match(checks):
    # The two analyses whose WATSTORE peaks carry lowercase "c" codes (#135).
    rows = checks.set_index("analysis_id").loc[["06151000.00", "06163400.00"]]
    assert (rows["skew_diff"].abs() < 0.001).all()
    assert (rows["flowfreq_n"] > 90).all()
