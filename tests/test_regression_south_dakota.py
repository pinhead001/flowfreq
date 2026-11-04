"""South Dakota's verified equations (WRI 98-4055 Tables 4, 5, 6 and 7, and its examples)."""

from __future__ import annotations

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
NSS = dict(zip("ABCDEFG", ["GC667", "GC668", "GC669", "GC670", "GC671", "GC672", "GC673"]))
AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002]
#: Table 6: station and weighted flows, with the characteristics, at the 197 gages.
GAGES = json.loads(
    (ROOT / "fixtures" / "regression" / "south_dakota_gages.json").read_text("utf-8")
)["gages"]

#: Weighted flows that equation 4 misses by more than rounding: gage and recurrence interval
#: to the printed flow. 06478250's 25-year flow is a misprint: equation 4 gives 525.
WEIGHTED_MISPRINT = {("06478250", 0.04): 55.0}
#: Unexplained, and at most 1 percent off.
WEIGHTED_UNEXPLAINED = {
    ("06446250", 0.2),
    ("06446250", 0.1),
    ("06355000", 0.002),
    ("06408500", 0.2),
    ("06354882", 0.02),
    ("06334000", 0.02),
}
#: 06453400's weighted flows imply 13 years of record. Table 6 prints 14 (1961-74).
WEIGHTED_YEARS = {"06453400": 13}


@pytest.fixture(scope="module")
def lib():
    return load_state("SD")


def _inputs(eq, gage):
    return {v.code: gage[v.code] for v in eq.variables}


def _half_unit(value: float, floor: float = 0.0) -> float:
    """Half the last printed unit of a value Table 6 gives to three significant figures."""
    return max(0.5 * 10 ** (math.floor(math.log10(value)) - 2), floor)


def _percent(log_sd: float) -> float:
    return 100.0 * math.sqrt(math.exp((log_sd * math.log(10.0)) ** 2) - 1.0)


def _log_variance(percent: float) -> float:
    return math.log(1.0 + (percent / 100.0) ** 2) / math.log(10.0) ** 2


def test_south_dakota_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == list("ABCDEFG")
    assert len(lib.equations) == 49
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region
    for e in lib.equations:
        assert e.covariance is not None and e.model_error_variance is not None
        assert e.equivalent_years is not None and e.n_sites is not None


def test_worked_example_mosquito_creek(lib):
    """Procedure 1, p. 21: subregion B, 9.2 mi2, PII 0.98."""
    est = evaluate(lib.equation("B", 0.01), {"CONTDA": 9.2, "PII_SD": 0.98})
    assert est.flow_cfs == pytest.approx(1080, abs=5)


def test_worked_example_willow_creek(lib):
    """Procedure 2, p. 21: 86.5 mi2, 66.9 in subregion C and 19.6 in D, averaged by area."""
    c = evaluate(lib.equation("C", 0.01), {"CONTDA": 86.5}).flow_cfs
    d = evaluate(lib.equation("D", 0.01), {"CONTDA": 86.5}).flow_cfs
    assert c == pytest.approx(4730, abs=5)
    assert d == pytest.approx(24400, abs=50)
    assert round(4730 * 66.9 / 86.5 + 24400 * 19.6 / 86.5, -1) == 9190


def test_worked_example_weighted_gage(lib):
    """P. 22: 06477400, subregion B, 12 years, 92 ft3/s at the station."""
    gage = GAGES["06477400"]
    assert (gage["years"], gage["CONTDA"], gage["PII_SD"]) == (12, 0.22, 0.91)
    eq = lib.equation("B", 0.01)
    assert eq.equivalent_years == 12.4
    regression = evaluate(eq, _inputs(eq, gage)).flow_cfs
    assert regression == pytest.approx(100, abs=0.5)
    assert gage["station"][5] == 92
    assert gage["weighted"][5] == 96
    assert (12 * 92 + 12.4 * 100) / 24.4 == pytest.approx(96, abs=0.5)


def test_transfer_examples_read_the_same_table_6_flows():
    """Pp. 22-23: Oak Creek, and the Little White River by both methods."""
    assert GAGES["06354882"]["CONTDA"] == 356
    assert GAGES["06354882"]["weighted"][5] == 10300
    assert 10300 * (269 / 356) ** 0.569 == pytest.approx(8780, abs=5)
    martin, vetal = GAGES["06447500"], GAGES["06449100"]
    assert (martin["CONTDA"], martin["weighted"][5]) == (230, 2110)
    assert (vetal["CONTDA"], vetal["weighted"][5]) == (415, 3500)
    assert 2110 * (325 / 230) ** 0.691 == pytest.approx(2680, abs=5)
    assert 3500 * (325 / 415) ** 0.691 == pytest.approx(2960, abs=5)


def test_gage_counts_equal_table_4(lib):
    """Subregion C's two gages with no flow from 25 years explain its 48 and 46."""
    assert len(GAGES) == 197
    for e in lib.equations:
        i = AEPS.index(e.aep)
        have = [g for g in GAGES.values() if g["region"] == e.region_code and g["weighted"][i]]
        assert len(have) == e.n_sites
    gone = {k for k, g in GAGES.items() if g["weighted"][3] is None}
    assert gone == {"06399300", "06400875"}


def test_ranges_equal_the_extremes_of_table_6(lib):
    """Table 6's largest subregion F slope, 460.1, is what Table 5 rounds to 460."""
    for region in lib.regions:
        gages = [g for g in GAGES.values() if g["region"] == region]
        for v in lib.equation(region, 0.01).variables:
            values = [g[v.code] for g in gages]
            assert v.minimum == min(values)
            assert v.maximum == pytest.approx(max(values), abs=0.5 if v.maximum >= 100 else 1e-9)


def _weighted_misses(lib, years=WEIGHTED_YEARS):
    """Every weighted flow equation 4 misses by more than the rounding of its inputs."""
    missed = {}
    total = 0
    for station, gage in GAGES.items():
        for i, aep in enumerate(AEPS):
            printed = gage["weighted"][i]
            if printed is None:
                continue
            total += 1
            eq = lib.equation(gage["region"], aep)
            est = evaluate(eq, _inputs(eq, gage), allow_extrapolation=True).flow_cfs
            n, en = years.get(station, gage["years"]), eq.equivalent_years
            at_site = gage["station"][i]
            weighted = (n * at_site + en * est) / (n + en)
            # relative error of the regression estimate from the rounding of its inputs
            rel = sum(
                abs(b) * _half_unit(gage[v.code], 0.005) / gage[v.code]
                for v, b in zip(eq.variables, eq.coefficients)
            )
            rel += _half_unit(10**eq.intercept) / 10**eq.intercept
            rel += 0.0005 * abs(math.log(gage["CONTDA"]))
            bound = (n * _half_unit(at_site, 0.5) + en * est * rel) / (n + en)
            bound += _half_unit(printed, 0.5)
            if abs(weighted - printed) > bound:
                missed[(station, aep)] = (printed, weighted)
    return missed, total


def test_weighted_flows_reproduce_at_every_gage_within_rounding(lib):
    """Table 6 prints no regression flows, so equation 4 is the per-gage check."""
    missed, total = _weighted_misses(lib)
    assert total == 1371
    assert set(missed) == set(WEIGHTED_MISPRINT) | WEIGHTED_UNEXPLAINED
    for key, printed in WEIGHTED_MISPRINT.items():
        assert missed[key][0] == printed
    assert missed[("06478250", 0.04)][1] == pytest.approx(525, abs=1)
    for key in WEIGHTED_UNEXPLAINED:
        assert missed[key][1] == pytest.approx(missed[key][0], rel=0.011)


def test_one_gage_is_weighted_on_13_years_not_the_14_printed(lib):
    assert GAGES["06453400"]["years"] == 14
    missed, _ = _weighted_misses(lib, years={})
    assert sum(station == "06453400" for station, _ in missed) >= 5


def test_covariance_reproduces_the_average_standard_error_of_prediction(lib):
    """There is no site-specific example, so Table 7 is checked against Table 4 for all 49."""
    for e in lib.equations:
        u = np.array(e.covariance)
        assert np.array_equal(u, u.T)
        assert np.linalg.eigvalsh(u).min() > 0
        i = AEPS.index(e.aep)
        rows = []
        for g in GAGES.values():
            if g["region"] == e.region_code and g["weighted"][i] is not None:
                rows.append([1.0] + [v.apply(g[v.code]) for v in e.variables])
        x = np.array(rows)
        xux = float(np.einsum("ij,jk,ik->i", x, u, x).mean())
        estimate = round(_percent(math.sqrt(e.model_error_variance)))
        low = 100.0 * math.sqrt(
            math.exp(math.log(10.0) ** 2 * (_log_variance(estimate - 0.5) + xux)) - 1
        )
        high = 100.0 * math.sqrt(
            math.exp(math.log(10.0) ** 2 * (_log_variance(estimate + 0.5) + xux)) - 1
        )
        printed = round(_percent(e.sep_log))
        assert low <= printed + 0.5 and high >= printed - 0.5, (e.region_code, e.aep)


def test_every_equation_equals_the_nss_snapshot(lib):
    snap = json.loads((DATA / "nss_snapshot" / "SD.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-12)
        assert list(e.coefficients) == s["coefficients"]
        assert e.sep_log == pytest.approx(s["sep_log"], abs=1e-5)
        for v, n in zip(e.variables, s["variables"]):
            assert (v.code, v.minimum, v.maximum) == (n["code"], n["minimum"], n["maximum"])


def test_nss_serves_the_report_under_another_title(lib):
    snap = json.loads((DATA / "nss_snapshot" / "SD.json").read_text(encoding="utf-8"))
    assert "A Method for Estimating Magnitude" in snap["source_reports"][0]["citation"]
    assert "Techniques for estimating peak-flow" in lib.equations[0].citation.publication


def test_intervals_use_the_covariance(lib):
    """MEV + xUx at the site, not the average SEP."""
    eq = lib.equation("B", 0.01)
    est = evaluate(eq, {"CONTDA": 9.2, "PII_SD": 0.98})
    assert "site-specific" in est.interval_method
