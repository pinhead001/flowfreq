"""The 2026 urban equations: SIR 2026-5044 (Alabama, Mississippi) and SIR 2025-5104 (Tennessee;
Georgia, North Carolina and South Carolina by the owner's decision of 2026-11-23).

Neither report prints a coefficient or a worked example, so the check is each data
release's own per-gage table, at every gage and AEP: regression flow to three significant
figures and variance of prediction to four decimals. ``urban_2026_wls.csv`` (42 gages,
doi:10.5066/P13GPMON) and ``urban_2026_gls.csv`` (136, doi:10.5066/P1TE9KTH) are columns of
the releases' ``WLSresults.csv`` and ``GLSresults.csv``, and ``urban_2026_wreg.json`` is
their ``coefs.txt``, ``CoVarMat.txt`` and ``PerformanceMetrics.txt``.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict

import numpy as np
import pandas as pd
import pytest

from flowfreq import estimate_from_characteristics
from flowfreq.regression import evaluate
from flowfreq.regression.equations import OutOfRangeError
from flowfreq.regression.library import load_state
from flowfreq.ungaged import _shares_used, region_for_nss_code

ROOT = Path(__file__).resolve().parent
FIXTURES = ROOT / "fixtures" / "regression"
WLS = pd.read_csv(FIXTURES / "urban_2026_wls.csv", dtype={"site_no": str})
GLS = pd.read_csv(FIXTURES / "urban_2026_gls.csv", dtype={"site_no": str})
WREG = json.loads((FIXTURES / "urban_2026_wreg.json").read_text("utf-8"))
AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
REGION = "Urban 2026"
GLS_STATES = ["TN", "GA", "NC", "SC"]


def _inputs(gage: Any) -> Dict[str, float]:
    return {k: float(v) for k, v in gage.items() if k[:2] in ("DR", "LC", "L3")}


def _gage(table: pd.DataFrame, site_no: str) -> Any:
    return table.set_index("site_no").loc[site_no]


def _three_figures(value: float, printed: float) -> bool:
    step = 10.0 ** (math.floor(math.log10(printed)) - 2)
    return abs(value - printed) <= 0.5 * step * (1 + 1e-9)


@pytest.mark.parametrize("state", ["AL", "MS", "TN", "GA", "NC", "SC"])
def test_the_set_is_stored_at_eight_aeps(state):
    lib = load_state(state)
    assert lib.status == "verified" and REGION in lib.regions
    assert lib.aeps(REGION) == AEPS
    # NSS serves no region for either report, so nothing locates the set.
    assert all("GC" not in e.region_name for e in lib.equations if e.region_code == REGION)
    assert region_for_nss_code(lib, "GC1504") != REGION
    assert lib.equation(REGION, 0.01).n_sites == (42 if state in ("AL", "MS") else 136)


def test_alabama_and_mississippi_hold_the_same_equations():
    al, ms = load_state("AL"), load_state("MS")
    for aep in AEPS:
        a, m = al.equation(REGION, aep), ms.equation(REGION, aep)
        assert (a.intercept, a.coefficients, a.covariance) == (
            m.intercept,
            m.coefficients,
            m.covariance,
        )
        assert [v.code for v in a.variables] == ["DRNAREA", "LC11DEV"]
        assert [v.transform for v in a.variables] == ["log10", "identity"]


@pytest.mark.parametrize("aep", AEPS)
@pytest.mark.parametrize("state", ["AL", "MS"])
def test_alabama_mississippi_release_gages_reproduce(state, aep):
    eq = load_state(state).equation(REGION, aep)
    assert len(WLS) == 42
    for _, gage in WLS.iterrows():
        est = evaluate(eq, _inputs(gage))
        assert _three_figures(est.flow_cfs, gage[f"WLS_{aep:g}"]), gage["site_no"]
        assert est.prediction_variance == pytest.approx(gage[f"VarPred_{aep:g}"], abs=5.1e-5)
        assert est.interval is not None and est.interval[0] < est.flow_cfs < est.interval[1]


@pytest.mark.parametrize("aep", AEPS)
@pytest.mark.parametrize("state", GLS_STATES)
def test_the_southeast_release_gages_reproduce(state, aep):
    eq = load_state(state).equation(REGION, aep)
    assert len(GLS) == 136
    assert any((GLS["L3_RDGVLY"] > 0) & (GLS["L3_RDGVLY"] < 100))  # partly Ridge and Valley
    for _, gage in GLS.iterrows():
        est = evaluate(eq, _inputs(gage))
        assert _three_figures(est.flow_cfs, gage[f"GLS_{aep:g}"]), gage["site_no"]
        assert est.prediction_variance == pytest.approx(gage[f"GLSvar_{aep:g}"], abs=5.1e-5)


@pytest.mark.parametrize("aep", AEPS)
@pytest.mark.parametrize("state", ["AL", "MS", *GLS_STATES])
def test_the_stored_numbers_are_the_releases_wreg_output(state, aep):
    eq = load_state(state).equation(REGION, aep)
    wreg = WREG["wls" if state in ("AL", "MS") else "gls"][f"{aep:g}"]
    assert [eq.intercept, *eq.coefficients] == wreg["coefficients"]
    assert eq.model_error_variance == wreg["ModErrVar"] and eq.avp == wreg["AVP"]
    sep = math.sqrt(math.log(1.0 + (wreg["Sp"] / 100.0) ** 2)) / math.log(10.0)
    assert eq.sep_log == pytest.approx(sep, abs=5.1e-7)
    if state in GLS_STATES:
        assert [list(row) for row in eq.covariance] == wreg["covariance"]
    errors = np.sqrt(np.diag(np.array(eq.covariance)))
    assert errors == pytest.approx(wreg["standard_errors"], rel=1e-5)


def test_the_release_counts_by_state():
    counts = GLS["state"].value_counts()
    assert [int(counts[s]) for s in ("GA", "NC", "SC")] == [32, 38, 16]


@pytest.mark.parametrize("state", ["GA", "NC", "SC"])
def test_the_southeast_states_hold_tennessees_equations(state):
    tn, lib = load_state("TN"), load_state(state)
    for aep in AEPS:
        a, b = tn.equation(REGION, aep), lib.equation(REGION, aep)
        assert (a.intercept, a.coefficients, a.covariance) == (
            b.intercept,
            b.coefficients,
            b.covariance,
        )
        assert (a.variables, a.applicability) == (b.variables, b.applicability)
        assert "Tennessee (" not in b.region_name


@pytest.mark.parametrize("state", GLS_STATES)
def test_the_set_is_used_only_when_named_and_cautions(state):
    lib = load_state(state)
    tiers = lib.region_preference.tiers if lib.region_preference else []
    assert all(REGION not in tier for tier in tiers)
    assert REGION not in lib.region_limits
    caution = lib.region_cautions[REGION]
    for words in (
        "not valid in the Middle Atlantic or Southern Coastal Plain",
        "not known outside the study area",
        "Mississippi Alluvial Plain, Central Appalachians or Blue Ridge",
        "percent of the basin, 0 to 100, not fractions",
        "inappropriate",
    ):
        assert words in caution, words
    assert ("South Carolina gages" in caution) is (state == "SC")
    # What else the USGS publishes for an urban basin there.
    assert ("SIR 2014-5030" in caution) is (state != "TN")
    assert ("the largest urban basin it used was 378, 5.67 and 43.9 mi2" in caution) is (
        state != "TN"
    )
    assert ("does not say it applies in place" in caution) is (state != "TN")


def test_the_south_carolina_caution_has_the_releases_own_residuals():
    """Observed (at-site EMA) over regression flow at the report's 16 South Carolina gages."""
    lib = load_state("SC")
    text = lib.region_cautions[REGION]
    gages = GLS[GLS["state"] == "SC"]
    assert len(gages) == 16
    for aep, words in ((0.1, "10 percent"), (0.01, "1 percent")):
        eq = lib.equation(REGION, aep)
        ratio = np.array(
            [g[f"EMA_{aep:g}"] / evaluate(eq, _inputs(g)).flow_cfs for _, g in gages.iterrows()]
        )
        mean = 10.0 ** float(np.mean(np.log10(ratio)))
        below = int((ratio < 1).sum())
        assert (
            f"{mean:.2f} at {words} AEP ({ratio.min():.2f} to {ratio.max():.2f}, {below} of 16 "
            "below 1)"
        ) in text


def test_a_georgia_urban_estimate_above_one_square_mile():
    """The basin the state's other sets refuse: over 1 mi2 and 10 percent impervious."""
    gage = _gage(GLS, "02203655")
    inputs = _inputs(gage)
    est = estimate_from_characteristics(
        "GA", inputs, region=REGION, lat=33.68, lon=-84.36, aeps=[0.01]
    )
    assert _three_figures(float(est.quantiles.loc[0, "flow_cfs"]), gage["GLS_0.01"])
    assert "SIR 2025-5104" in est.advisory


def test_the_derived_covariance_is_the_coefficients_own():
    """Alabama's and Mississippi's covariance is computed, not printed: its diagonal must be
    the squared standard errors the release prints for AEP 0.5 (coefs.txt)."""
    eq = load_state("MS").equation(REGION, 0.5)
    printed = [0.0705162760581306, 0.0389039993566192, 0.000784597132654258]
    assert np.sqrt(np.diag(np.array(eq.covariance))) == pytest.approx(printed, rel=1e-8)
    assert eq.intercept == pytest.approx(2.31965140648468, abs=1e-12)
    assert eq.coefficients == pytest.approx((0.6272601078372, 0.00539890185097196), abs=1e-12)
    assert eq.model_error_variance == pytest.approx(0.0128111922158582)
    assert eq.avp == pytest.approx(0.0138448440201676)


def test_tennessee_terms_and_limits():
    eq = load_state("TN").equation(REGION, 0.5)
    assert [v.code for v in eq.variables] == ["DRNAREA", "LC11DEV", "L3_PIEDMNT", "L3_RDGVLY"]
    assert eq.intercept == pytest.approx(2.23292096882482, abs=1e-12)
    assert eq.coefficients == pytest.approx(
        (0.651760453342856, 0.00619776358880669, -0.00265954922636201, -0.00211564620997548),
        abs=1e-12,
    )
    assert [(v.code, v.minimum) for v in eq.applicability] == [("LC11IMP", 10.0)]
    assert (eq.variables[0].minimum, eq.variables[0].maximum) == (0.164, 93.4)


@pytest.mark.parametrize(
    "state,inputs",
    [
        ("AL", {"DRNAREA": 5.0, "LC11DEV": 15.0}),
        ("MS", {"DRNAREA": 80.0, "LC11DEV": 50.0}),
        (
            "TN",
            {"DRNAREA": 5.0, "LC11DEV": 20.0, "L3_PIEDMNT": 0.0, "L3_RDGVLY": 0.0, "LC11IMP": 12.0},
        ),
        # Under 10 percent impervious the report sends a Tennessee site to the rural equations.
        (
            "TN",
            {"DRNAREA": 5.0, "LC11DEV": 40.0, "L3_PIEDMNT": 0.0, "L3_RDGVLY": 0.0, "LC11IMP": 8.0},
        ),
    ],
)
def test_outside_the_ranges_is_refused(state, inputs):
    with pytest.raises(OutOfRangeError):
        evaluate(load_state(state).equation(REGION, 0.01), inputs)


def test_the_older_urban_sets_stay_with_a_caution():
    al, tn = load_state("AL"), load_state("TN")
    assert al.equation("Urban", 0.01).n_sites == 23
    assert "supersedes" in al.region_cautions["Urban"]
    assert set(tn.region_cautions) == {"Urban statewide", "Memphis urban", REGION}
    for region in ("Urban statewide", "Memphis urban"):
        text = tn.region_cautions[region]
        assert "'should be used to predict streamflows'" in text and "10 percent" in text
        assert "does not name this set" in text
    assert REGION not in al.region_cautions and load_state("MS").region_cautions == {}


def test_an_estimate_names_the_region():
    gage = _gage(WLS, "02473047")
    est = estimate_from_characteristics(
        "MS",
        {"DRNAREA": gage["DRNAREA"], "LC11DEV": gage["LC11DEV"]},
        region=REGION,
        lat=32.3,
        lon=-90.2,
        aeps=[0.01],
    )
    assert _three_figures(float(est.quantiles.loc[0, "flow_cfs"]), gage["WLS_0.01"])
    assert est.quantiles.loc[0, "lower"] < est.quantiles.loc[0, "flow_cfs"]


def test_alabamas_old_urban_estimate_says_it_is_superseded():
    est = estimate_from_characteristics(
        "AL", {"CONTDA": 5.0, "URBAN": 50.0}, region="Urban", lat=33.5, lon=-86.8, aeps=[0.01]
    )
    assert "SIR 2026-5044" in est.advisory


# -- the two ecoregion percentages ----------------------------------------------------------

BASIN = {"DRNAREA": 20.0, "LC11DEV": 60.0, "LC11IMP": 30.0}


def _shares(piedmont, ridge, state="GA"):
    return estimate_from_characteristics(
        state,
        {**BASIN, "L3_PIEDMNT": piedmont, "L3_RDGVLY": ridge},
        region=REGION,
        lat=33.7,
        lon=-84.4,
        aeps=[0.01],
    )


@pytest.mark.parametrize("state", GLS_STATES)
def test_shares_of_one_basin_cannot_pass_100_percent(state):
    with pytest.raises(ValueError, match="L3_PIEDMNT and L3_RDGVLY.*200 percent"):
        _shares(100.0, 100.0, state)
    assert _shares(50.4, 50.0, state).quantiles.loc[0, "flow_cfs"] > 0  # rounding


@pytest.mark.parametrize("piedmont,ridge", [(1.0, 0.0), (0.6, 0.4), (0.0, 0.98)])
def test_a_share_that_looks_like_a_fraction_is_cautioned(piedmont, ridge):
    est = _shares(piedmont, ridge)
    assert "given as fractions" in est.advisory
    assert "given as fractions" in est.to_markdown()


@pytest.mark.parametrize("piedmont,ridge", [(100.0, 0.0), (0.0, 0.0), (99.0, 1.0), (40.0, 0.5)])
def test_a_share_in_percent_is_not(piedmont, ridge):
    assert "given as fractions" not in _shares(piedmont, ridge).advisory


def test_the_check_is_only_for_the_states_that_take_the_shares():
    assert load_state("GA").basin_shares == [("L3_PIEDMNT", "L3_RDGVLY")]
    est = estimate_from_characteristics(
        "AL",
        {"DRNAREA": 5.0, "LC11DEV": 50.0, "L3_PIEDMNT": 1.0, "L3_RDGVLY": 100.0},
        region=REGION,
        lat=33.5,
        lon=-86.8,
        aeps=[0.01],
    )
    assert "given as fractions" not in est.advisory


def test_another_set_of_the_state_ignores_the_shares():
    """StreamStats may return the codes; a set that does not use them is not refused over them."""
    est = estimate_from_characteristics(
        "AL",
        {"CONTDA": 5.0, "URBAN": 50.0, "L3_PIEDMNT": 100.0, "L3_RDGVLY": 100.0},
        region="Urban",
        lat=33.5,
        lon=-86.8,
        aeps=[0.01],
    )
    assert est.quantiles.loc[0, "flow_cfs"] > 0


def test_a_rural_set_of_the_same_state_ignores_the_shares():
    lib = load_state("GA")
    shares = {"L3_PIEDMNT": 100.0, "L3_RDGVLY": 100.0}
    assert _shares_used(lib, {REGION}, shares) == [("L3_PIEDMNT and L3_RDGVLY", 200.0)]
    assert _shares_used(lib, {"HR1-5", "HR1"}, shares) == []


def test_the_share_limits_are_100_percent_and_1():
    with pytest.raises(ValueError, match="102 percent"):
        _shares(60.0, 42.0)
    assert "given as fractions" not in _shares(1.5, 0.0).advisory
