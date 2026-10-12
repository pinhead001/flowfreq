"""Montana regional skew study v2 (tools/build_montana_skew_study_v2.py; tasks 5 and 7).

Reads the committed inputs and results. NOT a USGS study: the table row it
produces is provisional and returned only on opt-in.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from flowfreq.regional_skew import RegionalSkewUnavailable, regional_skew_at_huc

DATA = Path(__file__).resolve().parents[1] / "data" / "skew_study"
REGION = "Montana east of the Divide (provisional)"


@pytest.fixture(scope="module")
def inputs() -> pd.DataFrame:
    return pd.read_csv(DATA / "montana_v2_inputs.csv", dtype={"site_no": str, "huc": str})


@pytest.fixture(scope="module")
def results():
    return json.loads((DATA / "montana_v2_results.json").read_text(encoding="utf-8"))


def _used(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["excluded"].fillna("") == ""]


def test_neighbour_gages_in_shared_basins(inputs):
    used = _used(inputs)
    assert set(used["state"]) == {"MT", "WY", "ID"}
    nb = used[used["state"] != "MT"]
    assert nb["huc"].str[:4].isin(["1002", "1007", "1008", "1701"]).all()


def test_every_used_gage_has_a_basin(inputs):
    used = _used(inputs)
    assert len(used) == 300
    assert used[["lat_cent", "lon_cent", "ELEV", "PRECIP", "div_dist_mi"]].notna().all().all()


def test_regulated_without_coded_start_excluded(inputs):
    used = _used(inputs)
    flagged = used[used["reg_class"].isin(["regulated", "urban"])]
    assert flagged["regulated_from_wy"].notna().all()


def test_divide_distance_sign(inputs):
    used = _used(inputs)
    assert (used.loc[used["huc"].str[:2] == "17", "div_dist_mi"] < 0).all()
    assert (used.loc[used["huc"].str[:2].isin(["10", "09"]), "div_dist_mi"] > 0).all()


def test_selection_rule(results):
    models = results["models"]
    eligible = {"CONSTANT": models["CONSTANT"]["avp_new"]}
    for name, m in models.items():
        if "t_added" in m and all(t >= 2 for t in m["t_added"]):
            eligible[name] = m["avp_new"]
    assert results["selected"] == min(eligible, key=lambda k: eligible[k]) == "PRECIP"


def test_table_row_is_the_east_only_constant(results):
    from flowfreq.regional_skew import load_table

    row = load_table().set_index(["state", "skew_region"]).loc[("MT", REGION)]
    m = results["models"]["CONSTANT_EAST_ONLY"]
    assert row["status"] == "provisional"
    assert row["skew"] == pytest.approx(m["beta_constant"], abs=0.005)
    assert row["skew_mse"] == pytest.approx(m["avp_new"], abs=0.0005)
    assert row["effective_record_years"] == pytest.approx(m["effective_record_length"], abs=0.5)


def test_fort_benton_needs_opt_in():
    # Missouri River at Fort Benton, MT (06090800), HUC12 100301021602.
    with pytest.raises(RegionalSkewUnavailable, match="provisional"):
        regional_skew_at_huc("100301021602", "MT")
    rs = regional_skew_at_huc("100301021602", "MT", allow_provisional=True)
    assert rs.provisional and rs.skew_region == REGION and rs.skew == -0.04


def test_provisional_row_does_not_apply_in_wyoming():
    # A Bighorn HUC8 shared by MT and WY: WY gages were used, but the row is MT's.
    with pytest.raises(RegionalSkewUnavailable, match="WY is outside"):
        regional_skew_at_huc("10080010", "WY", allow_provisional=True)


def test_west_of_divide_keeps_verified_pnw():
    rs = regional_skew_at_huc("170102040104", "MT", allow_provisional=True)
    assert rs.skew_region == "Pacific Northwest" and not rs.provisional
