"""The Oregon regression update's station dataset (plan task 2).

Reads data/regression_study/oregon_stations.csv and oregon_sir2016_5083_check.csv,
written by tools/build_oregon_regression_dataset.py.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

DATA = Path(__file__).resolve().parents[1] / "data" / "regression_study"
AEPS = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)


@pytest.fixture(scope="module")
def stations() -> pd.DataFrame:
    return pd.read_csv(DATA / "oregon_stations.csv", dtype={"site_no": str, "huc": str})


@pytest.fixture(scope="module")
def used(stations: pd.DataFrame) -> pd.DataFrame:
    return stations[stations["excluded"].isna()]


@pytest.fixture(scope="module")
def check() -> pd.DataFrame:
    return pd.read_csv(DATA / "oregon_sir2016_5083_check.csv", dtype={"site_no": str})


def test_dataset_covers_oregon_and_its_neighbours(used):
    counts = used["state"].value_counts()
    assert counts["OR"] >= 500
    assert set(counts.index) == {"OR", "WA", "ID", "CA", "NV"}
    assert (used.loc[used["state"] != "OR", "border_miles"] <= 25).all()


def test_quantiles_rise_with_rarer_aeps(used):
    q = used[[f"q_{a}" for a in AEPS]].to_numpy(float)
    assert (q > 0).all()
    assert (q[:, 1:] >= q[:, :-1]).all()


def test_every_gage_has_ten_nonzero_peaks(used):
    assert ((used["n_peaks"] - used["n_zeros"]) >= 10).all()


def test_columbia_mainstem_excluded(stations):
    dalles = stations.set_index("site_no").loc["14105700"]
    assert str(dalles["excluded"]).startswith("drainage area >")


def test_regulated_mainstems_keep_only_pre_regulation_years(used):
    # Willamette at Salem and Albany: peaks coded for regulation from 1942.
    rows = used.set_index("site_no").loc[["14191000", "14174000"]]
    assert (rows["regulated_from_wy"] == 1942).all()
    assert (rows["last_wy"] < 1942).all()


def test_closed_basins_use_the_provisional_skew_in_oregon_only(used):
    closed = used[used["huc"].str[:4] == "1712"]
    in_or = closed["state"] == "OR"
    assert in_or.sum() >= 10
    assert closed.loc[in_or, "skew_source"].str.contains("closed basins, provisional").all()
    assert closed.loc[~in_or, "skew_source"].str.startswith("at-site only").all()


def test_klamath_is_at_site_only(used):
    klamath = used[used["huc"].str[:4] == "1801"]
    assert len(klamath) >= 20
    assert klamath["skew_source"].str.startswith("at-site only").all()
    assert klamath["skew_weighted"].isna().all()


def test_matches_sir_2016_5083_station_quantiles(check):
    ratios = check[[f"ratio_{a}" for a in AEPS]]
    assert len(check) >= 25
    assert (ratios.median() - 1).abs().max() < 0.01


def test_sir_2016_5083_outliers_all_have_perception_thresholds(check):
    # The check applies no Table A3 settings; every gage off by more than 5 %
    # at some AEP is one that has them.
    ratios = check[[f"ratio_{a}" for a in AEPS]]
    off = (ratios - 1).abs().max(axis=1) > 0.05
    assert off.sum() <= 5
    assert check.loc[off, "table_a3"].all()
