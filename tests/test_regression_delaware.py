"""Delaware's equations: SIR 2022-5005 (equations 3-18, Tables 5, 6, 8), SIR 2006-5146 (urban).

Table 6's Coastal Plain model error variance disagrees with Table 8 and the data release. The
owner chose Table 8 on 2026-11-05; see ``test_coastal_plain_mev_is_table_8_not_table_6``.
"""

from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import numpy as np
import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
NSS = {"Coastal Plain": "GC1881", "Piedmont": "GC1882", "Urban": "GC684"}
#: Characteristics, regression flows and variances at the release's 121 gages (doi 10.5066/P9B7CUVO).
GAGES = json.loads((ROOT / "fixtures" / "regression" / "delaware_gages.json").read_text("utf-8"))


@pytest.fixture(scope="module")
def lib():
    return load_state("DE", warn=False)


def _row(eq, gage):
    logs = [math.log10(gage[v.code] + (v.offset or 0.0)) for v in eq.variables]
    return np.array([1.0, *logs])


def _regressed(region):
    return [g for g in GAGES.values() if g["region"] == region and g["regression"]]


def test_delaware_is_verified_with_three_regions(lib):
    assert lib.status == "verified"
    assert lib.regions == ["Coastal Plain", "Piedmont", "Urban"]
    assert len(lib.equations) == 24
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region


def test_worked_example_horsepen_arm(lib):
    """p. 25: Coastal Plain, 5.21 mi2, 1.64 percent, 0.43 percent soil A, the 1-percent AEP."""
    est = evaluate(
        lib.equation("Coastal Plain", 0.01), {"DRNAREA": 5.21, "BSLDEM3M": 1.64, "SOILA": 0.43}
    )
    assert est.flow_cfs == pytest.approx(1450, abs=0.5)


def test_worked_example_white_clay_creek(lib):
    """p. 26: Piedmont, equation 16, printed to the unit."""
    inputs = {"DRNAREA": 79.2, "SOILD": 14.31, "IMPERV": 6.42, "STORNHD": 0.52}
    est = evaluate(lib.equation("Piedmont", 0.01), inputs)
    assert est.flow_cfs == pytest.approx(16749, abs=0.5)


def test_worked_example_gaged_site_quotes_the_release():
    """p. 25: Qr = 4,097 and Vr = 0.06384 at 01488500, as the fixture carries them."""
    gage = GAGES["01488500"]
    assert gage["flows"][5] == 4097
    assert gage["log_variance"][5] == pytest.approx(0.06384, abs=5e-6)


def test_published_regression_flows_reproduce_at_every_gage(lib):
    """Coefficients print to three decimals and flows to the unit, which bounds the difference."""
    assert len(GAGES) == 121
    assert sum(g["regression"] for g in GAGES.values()) == 94
    largest = 0.0
    for gage in GAGES.values():
        for i, aep in enumerate(lib.aeps(gage["region"])):
            eq = lib.equation(gage["region"], aep)
            inputs = {v.code: gage[v.code] for v in eq.variables}
            flow = evaluate(eq, inputs, allow_extrapolation=True).flow_cfs
            logs = _row(eq, gage)[1:]
            rounding = 0.0005 * (1 + np.abs(logs).sum()) + 0.5 / gage["flows"][i] / math.log(10)
            assert abs(math.log10(flow / gage["flows"][i])) <= rounding + 1e-6
            if gage["flows"][i] >= 200:
                largest = max(largest, abs(flow / gage["flows"][i] - 1))
    assert largest < 0.005


def test_piedmont_covariance_and_mev_give_the_release_variance(lib):
    for i, aep in enumerate(lib.aeps("Piedmont")):
        eq = lib.equation("Piedmont", aep)
        for gage in _regressed("Piedmont"):
            x = _row(eq, gage)
            var = eq.model_error_variance + float(x @ np.array(eq.covariance) @ x)
            assert var == pytest.approx(gage["log_variance"][i], abs=1e-4)


def test_piedmont_average_error_matches_table_8(lib):
    for aep in lib.aeps("Piedmont"):
        eq = lib.equation("Piedmont", aep)
        cov = np.array(eq.covariance)
        var = [
            eq.model_error_variance + float(_row(eq, g) @ cov @ _row(eq, g))
            for g in _regressed("Piedmont")
        ]
        assert math.sqrt(np.mean(var)) == pytest.approx(eq.sep_log, abs=3e-4)


#: Release variance less xUx, the same at all 54 gages. Table 8's SEE and SEP imply these too.
RELEASE_COASTAL_MEV = (
    0.053983,
    0.046343,
    0.045164,
    0.047609,
    0.050724,
    0.056055,
    0.062074,
    0.071523,
)
#: What Table 6 prints, not stored.
TABLE_6_COASTAL_MEV = (
    0.046572,
    0.041077,
    0.03893,
    0.040292,
    0.043251,
    0.047262,
    0.052952,
    0.061914,
)
TABLE_8_COASTAL_SEE = (57.5, 52.8, 52.0, 53.6, 55.6, 58.8, 62.4, 67.9)


def test_coastal_plain_mev_is_table_8_not_table_6(lib):
    """Table 8's SEE, squared, agrees with the release; Table 6 sits 12-16 percent below both."""
    for i, aep in enumerate(lib.aeps("Coastal Plain")):
        eq = lib.equation("Coastal Plain", aep)
        see = math.sqrt(math.log(1 + (TABLE_8_COASTAL_SEE[i] / 100) ** 2)) / math.log(10)
        assert eq.model_error_variance == pytest.approx(see**2, abs=5e-7)
        cov = np.array(eq.covariance)
        rest = [
            g["log_variance"][i] - float(_row(eq, g) @ cov @ _row(eq, g))
            for g in _regressed("Coastal Plain")
        ]
        assert len(rest) == 54
        assert max(rest) - min(rest) < 5e-5
        assert np.mean(rest) == pytest.approx(RELEASE_COASTAL_MEV[i], abs=1e-6)
        assert eq.model_error_variance == pytest.approx(RELEASE_COASTAL_MEV[i], rel=2.5e-3)
        assert eq.sep_log**2 - np.mean(
            [float(_row(eq, g) @ cov @ _row(eq, g)) for g in _regressed("Coastal Plain")]
        ) == pytest.approx(eq.model_error_variance, abs=2.5e-4)
        assert eq.model_error_variance - TABLE_6_COASTAL_MEV[i] > 0.005


def test_worked_example_variance_with_table_8(lib):
    """p. 25: Vr = 0.06384 at 01488500, from the stored MEV plus xUx, to the rounding of Table 8."""
    gage = GAGES["01488500"]
    eq = lib.equation("Coastal Plain", 0.01)
    x = _row(eq, gage)
    var = eq.model_error_variance + float(x @ np.array(eq.covariance) @ x)
    assert var == pytest.approx(0.06384, abs=1e-4)


#: Where NSS differs from Table 14, as printed there and as a percent NSS's value implies.
NSS_URBAN_SEP_ERROR = {0.2: 91.1, 0.04: 93.6}


def test_every_equation_equals_the_nss_snapshot_but_for_its_errors(lib):
    snap = json.loads((DATA / "nss_snapshot" / "DE.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-9)
        assert list(e.coefficients) == s["coefficients"]
        if e.region_code == "Urban" and e.aep in NSS_URBAN_SEP_ERROR:
            assert abs(e.sep_log - s["sep_log"]) > 5e-4
        else:
            assert e.sep_log == pytest.approx(s["sep_log"], abs=1e-5)
        for v, n in zip(e.variables, s["variables"]):
            assert (v.code, v.minimum) == (n["code"], n["minimum"])
            if e.region_code == "Urban" and v.code == "DRNAREA":
                assert (v.maximum, n["maximum"]) == (319.0, 117.0)
            else:
                assert v.maximum == n["maximum"]


def test_urban_set_is_printed_constants_with_the_sep_and_no_covariance(lib):
    """Equations 26-33 and Table 14. The report has no urban example or per-gage estimate."""
    constants = (51.4, 91.8, 126, 179, 225, 277, 334, 420)
    sep = (90.6, 91.1, 91.9, 93.6, 95.4, 97.6, 100, 105)
    for i, aep in enumerate(lib.aeps("Urban")):
        eq = lib.equation("Urban", aep)
        assert 10**eq.intercept == pytest.approx(constants[i], rel=1e-12)
        assert eq.covariance is None
        pct = sep[i] / 100
        assert eq.sep_log == pytest.approx(math.sqrt(math.log(1 + pct**2)) / math.log(10), abs=1e-6)
        assert eq.n_sites == 102
    est = evaluate(lib.equation("Urban", 0.01), {"DRNAREA": 10.0, "HOMEDENS": 1.0})
    assert est.flow_cfs == pytest.approx(277 * 10**0.758 * 2**0.663, rel=1e-12)
