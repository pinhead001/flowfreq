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

#: Analyses where flowfreq and PeakFQ 7.5.1 disagree by more than 0.2, with why.
KNOWN_DISAGREEMENTS = {
    # PeakFQ 0.216 (57 peaks); flowfreq native 3.0 and the vendored emafitpr
    # (peakfq 8.1.0) 3.26 on the same converted inputs (70 rows). The EMA engines
    # agree; the difference is in the .psf conversion of this record (a 0-1 cfs
    # interval peak in 1968 and stacked 6000 cfs historical thresholds) or between
    # PeakFQ 7.5.1 and 8.1.0. Not yet resolved.
    "06163400.00": "conversion or PeakFQ-version difference; EMA engines agree",
}


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


@pytest.mark.parametrize("analysis_id", sorted(KNOWN_DISAGREEMENTS))
@pytest.mark.xfail(strict=True, reason="documented PeakFQ disagreement; see KNOWN_DISAGREEMENTS")
def test_known_disagreements(checks, analysis_id):
    row = checks.set_index("analysis_id").loc[analysis_id]
    assert abs(row["skew_diff"]) < 0.2
