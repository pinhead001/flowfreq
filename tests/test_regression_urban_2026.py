"""The 2026 urban equations: SIR 2026-5044 (Alabama, Mississippi) and SIR 2025-5104 (Tennessee).

Neither report prints a coefficient or a worked example, so the check is each data
release's own per-gage table: regression flow to three significant figures and variance of
prediction to four decimals. The fixture copies three gages from each
(``WLSresults.csv`` of doi:10.5066/P13GPMON, ``GLSresults.csv`` of doi:10.5066/P1TE9KTH);
``tools/build_urban_2026.py`` checks all 42 and all 136 before it writes the state files.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from flowfreq import estimate_from_characteristics
from flowfreq.regression import evaluate
from flowfreq.regression.equations import OutOfRangeError
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
GAGES = json.loads((ROOT / "fixtures" / "regression" / "urban_2026_gages.json").read_text("utf-8"))
AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
REGION = "Urban 2026"


def _three_figures(value: float, printed: float) -> bool:
    step = 10.0 ** (math.floor(math.log10(printed)) - 2)
    return abs(value - printed) <= 0.5 * step * (1 + 1e-9)


@pytest.mark.parametrize("state", ["AL", "MS", "TN"])
def test_the_set_is_stored_at_eight_aeps(state):
    lib = load_state(state)
    assert lib.status == "verified" and REGION in lib.regions
    assert lib.aeps(REGION) == AEPS
    # NSS serves no region for either report, so nothing locates the set.
    assert all("GC" not in e.region_name for e in lib.equations if e.region_code == REGION)
    assert region_for_nss_code(lib, "GC1504") != REGION
    assert lib.equation(REGION, 0.01).n_sites == (136 if state == "TN" else 42)


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


@pytest.mark.parametrize("gage", GAGES["alms"], ids=lambda g: g["site_no"])
def test_alabama_mississippi_release_gages_reproduce(gage):
    lib = load_state(gage["state"])
    inputs = {"DRNAREA": gage["DRNAREA"], "LC11DEV": gage["LC11DEV"]}
    for aep, printed in gage["flows"].items():
        est = evaluate(lib.equation(REGION, float(aep)), inputs)
        assert _three_figures(est.flow_cfs, printed), (aep, est.flow_cfs, printed)
        assert est.prediction_variance == pytest.approx(gage["var"][aep], abs=5.1e-5)
        assert est.interval is not None and est.interval[0] < est.flow_cfs < est.interval[1]


@pytest.mark.parametrize("gage", GAGES["tn"], ids=lambda g: g["site_no"])
def test_tennessee_release_gages_reproduce(gage):
    lib = load_state("TN")
    inputs = {k: gage[k] for k in ("DRNAREA", "LC11DEV", "L3_PIEDMNT", "L3_RDGVLY", "LC11IMP")}
    for aep, printed in gage["flows"].items():
        est = evaluate(lib.equation(REGION, float(aep)), inputs)
        assert _three_figures(est.flow_cfs, printed), (aep, est.flow_cfs, printed)
        assert est.prediction_variance == pytest.approx(gage["var"][aep], abs=5.1e-5)


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
    assert set(tn.region_cautions) == {"Urban statewide", "Memphis urban"}
    assert all("does not say it supersedes" in t for t in tn.region_cautions.values())
    assert REGION not in al.region_cautions and load_state("MS").region_cautions == {}


def test_an_estimate_names_the_region():
    gage = GAGES["alms"][1]
    est = estimate_from_characteristics(
        "MS",
        {"DRNAREA": gage["DRNAREA"], "LC11DEV": gage["LC11DEV"]},
        region=REGION,
        lat=32.3,
        lon=-90.2,
        aeps=[0.01],
    )
    assert _three_figures(float(est.quantiles.loc[0, "flow_cfs"]), gage["flows"]["0.01"])
    assert est.quantiles.loc[0, "lower"] < est.quantiles.loc[0, "flow_cfs"]


def test_alabamas_old_urban_estimate_says_it_is_superseded():
    est = estimate_from_characteristics(
        "AL", {"CONTDA": 5.0, "URBAN": 50.0}, region="Urban", lat=33.5, lon=-86.8, aeps=[0.01]
    )
    assert "SIR 2026-5044" in est.advisory
