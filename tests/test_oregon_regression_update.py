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
    "SOILKSAT_STATSGO2": "SOILKSAT",
    "VOLCQP_SGMC": "VOLCQP",
    "LAKEAREA_NHD": "LAKEAREA",
    "VOLCQP_OGDC": "VOLCQP_OG",
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
    # STATSGO2 Ksat, in./hr: every accepted basin has one
    assert ok["SOILKSAT"].between(0.05, 50).all()
    # percentages of the basin
    assert ok["VOLCQP"].between(0, 100).all()
    assert ok["LAKEAREA"].between(0, 100).all()


def test_high_cascades_basins_are_on_young_volcanics(basins):
    # Thielsen Creek, Clearwater above Trap Creek, North Umpqua at Toketee Falls,
    # Waldo Lake outlet: all at least half on Pliocene-Quaternary volcanic rock,
    # against a dataset median of zero; Waldo is a third lake.
    ok = basins[basins["basin_source"].notna()].set_index("site_no")
    assert (ok.loc[["14312700", "14314500", "14315500", "14147000"], "VOLCQP"] > 50).all()
    assert ok["VOLCQP"].median() == 0
    assert ok.loc["14147000", "LAKEAREA"] > 25


def test_high_cascades_soils_are_more_permeable_than_the_median(basins):
    # Waldo Lake outlet, Clearwater River above Trap Creek: pumice and ash soils.
    ok = basins[basins["basin_source"].notna()].set_index("site_no")
    assert (ok.loc[["14147000", "14314500"], "SOILKSAT"] > 1.5 * ok["SOILKSAT"].median()).all()


def test_every_characteristic_records_its_source():
    sources = json.loads((DATA / "oregon_basin_characteristics.sources.json").read_text())
    assert set(sources) == {
        "DRNAREA",
        "PRECIP",
        "JANMAXT",
        "I24H2Y",
        "ELEV",
        "BSLOPD",
        "SOILKSAT",
        "VOLCQP",
        "LAKEAREA",
        "VOLCQP_OG",
    }


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
    assert codes & {"PRECIP", "ELEV", "I24H2Y", "JANMAXT", "BSLOPD", "SOILPERM"} == set()


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


def test_the_state_map_counts_more_young_volcanic_rock_in_the_high_cascades(basins):
    # OGDC-8 dates the High Cascades platform lavas Quaternary where SGMC
    # counts much less: Waldo Lake outlet 51 % (SGMC) -> 95 % (OGDC-8).
    ok = basins[basins["basin_source"].notna()].set_index("site_no")
    hc = ["14147000", "14314500", "14315500"]
    assert (ok.loc[hc, "VOLCQP_OG"] > ok.loc[hc, "VOLCQP"] + 10).all()
    assert ok.loc["14147000", "VOLCQP_OG"] > 90
    assert ok["VOLCQP_OG"].between(0, 100).all()
    assert ok["VOLCQP_OG"].median() == 0


def _sep_percent(e) -> float:
    return 100.0 * float(np.sqrt(np.exp(np.log(10.0) ** 2 * e.sep_log**2) - 1.0))


@pytest.mark.parametrize(
    "region, codes, before",
    [
        # (50 %, 1 %, 0.2 % AEP) SEPs of the equations these replaced on 2026-10-18
        ("W-171003", {"DRNAREA", "PRECIP_PRISM9120", "VOLCQP_OGDC"}, (47.43, 48.98, 53.65)),
        (
            "E",
            {"DRNAREA", "I24H2Y_NA2", "ELEV_3DEP", "BSLOPD_3DEP", "VOLCQP_OGDC"},
            (89.81, 86.48, 100.77),
        ),
    ],
)
def test_multi_aep_selection_replaced_a_region_only_where_no_worse(
    provisional, results, region, codes, before
):
    assert results["select_mode"] == "multi"
    for aep, sep_before in zip((0.5, 0.01, 0.002), before):
        e = provisional.equation(region, aep)
        assert {v.code for v in e.variables} == codes
        assert _sep_percent(e) <= sep_before
    fits = results["regions"][region]["fits"]
    assert max(fits[a]["max_leverage"] for a in ("0.5", "0.1", "0.01", "0.002")) <= 0.5


def test_w171003_volcanic_term_lowers_floods(provisional):
    for aep in AEPS:
        e = provisional.equation("W-171003", aep)
        coef = dict(zip((v.code for v in e.variables), e.coefficients))
        assert coef["VOLCQP_OGDC"] < 0
        assert coef["PRECIP_PRISM9120"] > 0
