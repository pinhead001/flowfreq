"""The provisional Oregon regression update (plan tasks 3-5).

Reads data/regression_study/oregon_basin_characteristics.csv and
oregon_regression_results.json (tools/build_oregon_basin_characteristics.py,
tools/oregon_regression_fit.py) and the stored provisional set,
flowfreq/data/regression/provisional/OR.json.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from flowfreq.regression import evaluate, load_state

DATA = Path(__file__).resolve().parents[1] / "data" / "regression_study"
AEPS = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)
#: Stored code -> characteristics CSV column.
COLUMNS = {
    "DRNAREA": "DRNAREA",
    "PRECIP_PRISM9120": "PRECIP",
    "I24H2Y_NA2": "I24H2Y",
    "JANMAXT_PRISM9120": "JANMAXT",
    "ELEV_3DEP": "ELEV",
    "BSLOPD_3DEP": "BSLOPD",
}


@pytest.fixture(scope="module")
def basins() -> pd.DataFrame:
    return pd.read_csv(DATA / "oregon_basin_characteristics.csv", dtype={"site_no": str})


@pytest.fixture(scope="module")
def results() -> dict:
    return json.loads((DATA / "oregon_regression_results.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def provisional():
    return load_state("OR", allow_provisional=True)


def test_accepted_basins_match_the_nwis_area(basins):
    ok = basins[basins["basin_source"].notna()]
    assert len(ok) >= 550
    assert set(ok["basin_source"]) == {"StreamStats", "NLDI"}
    assert ((ok["area_ratio"] - 1).abs() <= 0.10).all()
    rejected = basins[basins["basin_source"].isna()]
    assert rejected["basin_error"].str.startswith("no basin within 10%").all()


def test_characteristics_are_physically_plausible(basins):
    ok = basins[basins["basin_source"].notna()]
    assert ok["PRECIP"].between(5, 200).all()
    assert ok["JANMAXT"].between(20, 60).all()
    assert ok["I24H2Y"].dropna().between(0.5, 8).all()
    assert ok["ELEV"].dropna().between(0, 10000).all()
    assert ok["BSLOPD"].dropna().between(0, 45).all()


def test_every_characteristic_records_its_source():
    sources = json.loads((DATA / "oregon_basin_characteristics.sources.json").read_text())
    assert set(sources) == {"DRNAREA", "PRECIP", "JANMAXT", "I24H2Y", "ELEV", "BSLOPD"}


def test_default_oregon_lookup_is_the_published_set():
    lib = load_state("OR")
    assert lib.status == "verified" and not lib.provisional
    assert {e.region_code for e in lib.equations} >= {"1", "2A", "2B", "E1"}


def test_provisional_set_is_opt_in_and_labelled(provisional):
    assert provisional.provisional
    assert provisional.advisory.startswith("PROVISIONAL (not USGS)")
    assert len(provisional.equations) == 8 * len({e.region_code for e in provisional.equations})


def test_variable_codes_name_their_source(provisional):
    # A StreamStats PRECIP or ELEV is computed from other grids; a same-named code
    # would let one be substituted silently.
    codes = {v.code for e in provisional.equations for v in e.variables}
    assert codes <= set(COLUMNS)
    assert codes & {"PRECIP", "ELEV", "I24H2Y", "JANMAXT", "BSLOPD"} == set()


def test_regions_have_enough_gages(results):
    assert all(r["n"] >= 25 for r in results["regions"].values())
    assert results["n_sites"] >= 550


def test_stored_equations_reproduce_the_fit(provisional, basins, results):
    ok = basins[basins["basin_source"].notna()].dropna(subset=list(COLUMNS.values()))
    site = ok.iloc[len(ok) // 2]
    for e in provisional.equations:
        if e.aep != 0.01:
            continue
        chars = {v.code: float(site[COLUMNS[v.code]]) for v in e.variables}
        est = evaluate(e, chars, allow_extrapolation=True)
        x = [1.0]
        for v in e.variables:
            raw = v.scale * chars[v.code] + v.offset
            x.append(np.log10(raw) if v.transform == "log10" else raw)
        beta = results["regions"][e.region_code]["fits"]["0.01"]["beta"]
        assert est.log_flow == pytest.approx(float(np.dot(x, beta)), abs=1e-9)


def test_covariances_are_positive_definite(provisional):
    for e in provisional.equations:
        assert np.all(np.linalg.eigvalsh(np.asarray(e.covariance)) > 0)
        assert e.sep_log == pytest.approx(np.sqrt(e.avp))


def test_quantiles_rise_with_rarer_aeps_at_each_regions_median_basin(provisional, results):
    df = pd.read_csv(DATA / "oregon_basin_characteristics.csv", dtype={"site_no": str})
    med = {code: float(df[col].median()) for code, col in COLUMNS.items()}
    for region in results["regions"]:
        flows = [
            evaluate(provisional.equation(region, a), med, allow_extrapolation=True).flow_cfs
            for a in AEPS
        ]
        assert all(b > a for a, b in zip(flows, flows[1:])), region
