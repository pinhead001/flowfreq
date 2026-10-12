"""Snake River Plain skew feasibility (tools/snake_river_plain_skew.py, task 8)."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from flowfreq.regional_skew import RegionalSkewUnavailable, regional_skew_for

DATA = Path(__file__).resolve().parents[1] / "data" / "skew_study"


@pytest.fixture(scope="module")
def results():
    return json.loads((DATA / "srp_results.json").read_text(encoding="utf-8"))


def test_decision_follows_the_rule(results):
    assert results["min_sites"] == 20
    assert results["screened_in"] < results["min_sites"]
    assert results["decision"] == "not feasible"
    assert "constant" not in results


def test_inputs_list_every_region0_site(results):
    df = pd.read_csv(DATA / "srp_inputs.csv", dtype={"site_no": str})
    assert len(df) == results["region0_peak_sites"] == 47
    assert set(df["region0"]) == {"excluded", "boundary-in"}
    used = df[df["excluded"].fillna("") == ""]
    assert set(used["site_no"]) == {"13068500", "13075983"}


@pytest.mark.parametrize("allow", [False, True])
def test_no_value_for_the_plain(allow):
    with pytest.raises(RegionalSkewUnavailable):
        regional_skew_for("ID", "Snake River Plain", allow_provisional=allow)
