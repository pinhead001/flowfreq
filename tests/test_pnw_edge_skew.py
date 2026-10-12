"""PNW skew applicability in Oregon's closed basins and the Klamath (task 9)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from flowfreq.regional_skew import RegionalSkewUnavailable, regional_skew_at_huc

DATA = Path(__file__).resolve().parents[1] / "data" / "skew_study"
REGION = "Pacific Northwest (Oregon closed basins, provisional)"


@pytest.fixture(scope="module")
def results():
    return json.loads((DATA / "pnw_edge_results.json").read_text(encoding="utf-8"))


def test_decision_follows_the_rule(results):
    closed = results["by_basin"]["Oregon closed basins"]
    assert closed["n"] == 8 and abs(closed["t_vs_pnw"]) < 2
    assert closed["decision"] == "consistent"
    assert results["by_basin"]["Klamath"] == {"n": 2}


def test_closed_basin_needs_opt_in():
    # Donner und Blitzen River near Frenchglen, OR (10396000), HUC12 171200030108.
    with pytest.raises(RegionalSkewUnavailable, match="provisional"):
        regional_skew_at_huc("171200030108", "OR")
    rs = regional_skew_at_huc("171200030108", "OR", allow_provisional=True)
    assert rs.provisional and rs.skew_region == REGION
    assert (rs.skew, rs.skew_mse) == (-0.07, 0.18)


def test_closed_basin_provisional_is_oregon_only():
    with pytest.raises(RegionalSkewUnavailable, match="NV is outside"):
        regional_skew_at_huc("171200030108", "NV", allow_provisional=True)


def test_klamath_stays_unresolved():
    # Williamson River below Sprague River, OR (11502500).
    with pytest.raises(RegionalSkewUnavailable, match="unresolved"):
        regional_skew_at_huc("180102010605", "OR", allow_provisional=True)
