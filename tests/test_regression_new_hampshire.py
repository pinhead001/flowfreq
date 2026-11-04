"""New Hampshire's verified equations (SIR 2008-5206 equations 3-9 and 20-26, Tables 5 and 7-12)."""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
NSS = {"1": "GC1460", "2": "GC1459"}
AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002]
CODES = ("DRNAREA", "APRAVPRE", "WETLAND", "CSL10_85")
#: Table 7's characteristics and Table 5's at-site flows, at the 117 gages regressed.
GAGES = json.loads(
    (ROOT / "fixtures" / "regression" / "new_hampshire_gages.json").read_text("utf-8")
)
#: Standard error of regression, log10 units: Table 10 (region 1) and Table 12 (region 2).
SE_REGRESSION = {
    "1": [0.132, 0.139, 0.147, 0.159, 0.170, 0.181, 0.210],
    "2": [0.231, 0.238, 0.243, 0.251, 0.258, 0.265, 0.283],
}
#: The matrix as Appendix 3 prints it for the 100-year equation, read from page 50.
APPENDIX_3_MATRIX = [
    [0.060575, -0.0051667, -0.050248, -0.00059817, -0.010587],
    [-0.0051667, 0.0013231, -0.0013914, 0.00010560, 0.0019168],
    [-0.050248, -0.0013914, 0.10452, -0.00036357, -0.0050145],
    [-0.00059817, 0.00010560, -0.00036357, 0.000029723, 0.00028972],
    [-0.010587, 0.0019168, -0.0050145, 0.00028972, 0.0053530],
]


def _sig3(value: float) -> float:
    return float(f"{value:.3g}")


def _design() -> np.ndarray:
    a = np.array([[g[c] for c in CODES] for g in GAGES.values()])
    return np.c_[np.ones(len(a)), np.log10(a[:, 0]), np.log10(a[:, 1]), a[:, 2], np.log10(a[:, 3])]


def _mean_equivalent_years(variance: np.ndarray, aep: float, skew: float) -> float:
    """Equation 13, averaged over the gages."""
    k = float(stats.pearson3.ppf(1 - aep, skew))
    s = np.exp(-(1.31 + 0.134 * _design()[:, 1]))
    return float(np.mean(s**2 * (1 + k * skew + 0.5 * k**2 * (1 + 0.75 * skew**2)) / variance))


@pytest.fixture(scope="module")
def lib():
    return load_state("NH")


def test_new_hampshire_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == ["1", "2"]
    assert len(lib.equations) == 14
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region
        assert lib.aeps(region) == AEPS


def test_worked_example(lib):
    """Appendix 3: 22.0 mi2, 4.05 in., 3.27 percent, 62.6 ft/mi, the 100-year flood."""
    eq = lib.equation("1", 0.01)
    est = evaluate(eq, {"DRNAREA": 22.0, "APRAVPRE": 4.05, "WETLAND": 3.27, "CSL10_85": 62.6})
    assert _sig3(est.flow_cfs) == 3110
    x = np.array([1.0, np.log10(22.0), np.log10(4.05), 3.27, np.log10(62.6)])
    xux = float(x @ np.array(eq.covariance) @ x)
    assert round(xux, 5) == 0.00117
    sep = (eq.model_error_variance + xux) ** 0.5
    assert round(sep, 3) == 0.157
    assert est.prediction_variance == pytest.approx(sep**2)
    assert round(100 * (10**0.157 - 1), 1) == pytest.approx(43.5, abs=0.051)
    assert round(100 * (10**-0.157 - 1), 1) == -30.3
    assert (_sig3(est.interval[0]), _sig3(est.interval[1])) == (1710, 5660)
    s = np.exp(-(1.31 + 0.134 * np.log10(22.0)))
    assert round(float(s), 3) == 0.225
    assert round(float(s**2 * (1 + 0.5 * 2.326**2) / sep**2), 1) == 7.6


def test_worked_example_at_the_ungaged_site(lib):
    """Appendix 3, equations 16-19, with the report's 3,700 ft3/s at the gage."""
    q_g = 3110.0
    est = evaluate(
        lib.equation("1", 0.01),
        {"DRNAREA": 30.0, "APRAVPRE": 4.00, "WETLAND": 3.40, "CSL10_85": 59.0},
    )
    q_u = _sig3(est.flow_cfs)
    assert q_u == 3890
    q_w = 10 ** ((10 * np.log10(3700.0) + 7.6 * np.log10(q_g)) / 17.6)
    assert _sig3(q_w) == 3430
    m = np.log10(q_u / q_g) / np.log10(30.0 / 22.0)
    assert round(float(m), 3) == 0.722
    c = 0.722 + np.log10(q_g / 3430.0) / np.log10(1.5)
    assert round(float(c), 3) == 0.480
    assert _sig3(3430.0 * (30.0 / 22.0) ** 0.480) == 3980


def test_table_11_matches_the_matrix_printed_in_appendix_3(lib):
    """Two printings of the 100-year matrix, including Table 11's odd "0.10452E-00"."""
    assert lib.equation("1", 0.01).covariance == tuple(tuple(r) for r in APPENDIX_3_MATRIX)
    assert lib.equation("1", 0.002).covariance[2][2] == 0.14008


def test_table_11_stated_order_is_the_only_one_that_fits_the_example(lib):
    u = np.array(lib.equation("1", 0.01).covariance)
    logs = [np.log10(22.0), np.log10(4.05), 3.27, np.log10(62.6)]
    near = []
    for order in itertools.permutations(range(4)):
        x = np.array([1.0] + [logs[i] for i in order])
        if round(float(x @ u @ x), 5) == 0.00117:
            near.append(order)
    assert near == [(0, 1, 2, 3)]


@pytest.mark.parametrize("aep", AEPS)
def test_matrix_reproduces_the_printed_average_error(lib, aep):
    """Tables 10 and 11: sqrt(mean(MEV + xUx)) over the 117 gages, printed to 0.001."""
    eq = lib.equation("1", aep)
    x = _design()
    variance = eq.model_error_variance + np.einsum("ij,jk,ik->i", x, np.array(eq.covariance), x)
    assert float(np.mean(variance)) ** 0.5 == pytest.approx(eq.sep_log, abs=5.1e-4)


def test_gage_table_reproduces_table_8():
    assert len(GAGES) == 117
    a = np.array([[g[c] for c in CODES] for g in GAGES.values()])
    assert a.min(axis=0).tolist() == [0.70, 2.79, 0.0, 5.43]
    assert a.max(axis=0).tolist() == [1290.0, 6.23, 21.8, 543.0]
    for mean, printed, half_unit in zip(
        a.mean(axis=0), [83.6, 3.94, 4.86, 114.0], [0.05, 0.005, 0.005, 0.5]
    ):
        assert abs(mean - printed) <= half_unit


@pytest.mark.parametrize("region", ["1", "2"])
@pytest.mark.parametrize("i", range(7))
def test_at_site_residuals_reproduce_the_standard_error_of_regression(lib, region, i):
    """No per-gage regression flow is printed. Table 5 on Table 7 tests every coefficient at once."""
    eq = lib.equation(region, AEPS[i])
    residuals = []
    for gage in GAGES.values():
        inputs = {v.code: gage[v.code] for v in eq.variables}
        est = evaluate(eq, inputs, allow_extrapolation=True)
        residuals.append(np.log10(gage["flows_at_site"][i]) - est.log_flow)
    rms = (np.sum(np.square(residuals)) / (117 - len(eq.variables) - 1)) ** 0.5
    assert rms == pytest.approx(SE_REGRESSION[region][i], abs=1e-3)
    assert abs(np.mean(residuals)) < 0.015


def test_printed_equivalent_years_need_a_skew_the_report_does_not_state(lib):
    """Equation 13 with zero skew, as the report says to use, fits only the 2-year averages.

    A skew of 0.25 fits all 14 printed values to about 0.1 year. That value is inferred. Nothing
    reads ``equivalent_years``, and a site's own value for equation 16 comes from equation 13.
    """
    x = _design()
    for i, aep in enumerate(AEPS):
        full = lib.equation("1", aep)
        da_only = lib.equation("2", aep)
        u = np.array(full.covariance)
        v_full = full.model_error_variance + np.einsum("ij,jk,ik->i", x, u, x)
        v_da = np.full(len(x), da_only.sep_log**2)
        for eq, variance in ((full, v_full), (da_only, v_da)):
            assert _mean_equivalent_years(variance, aep, 0.25) == pytest.approx(
                eq.equivalent_years, abs=0.15
            )
            at_zero = _mean_equivalent_years(variance, aep, 0.0)
            if aep == 0.5:
                assert at_zero == pytest.approx(eq.equivalent_years, abs=0.05)
            elif eq is full:
                assert at_zero < eq.equivalent_years - 0.5


#: Table 9, percent change for a +10 then -10 percent change in the input.
TABLE_9 = {
    "DRNAREA": (
        [9.6, 9.3, 9.1, 8.9, 8.7, 8.6, 8.3],
        [-9.6, -9.3, -9.2, -9.0, -8.8, -8.7, -8.5],
    ),
    "APRAVPRE": (
        [15.4, 18.0, 19.0, 20.0, 20.4, 20.7, 21.4],
        [-14.6, -16.7, -17.5, -18.2, -18.5, -18.8, -19.3],
    ),
    "WETLAND": (
        [-2.7, -2.7, -2.7, -2.8, -2.8, -2.8, -2.9],
        [2.8, 2.8, 2.8, 2.8, 2.9, 2.9, 2.9],
    ),
    "CSL10_85": (
        [1.8, 2.0, 2.0, 2.0, 2.0, 1.9, 1.8],
        [-2.1, -2.2, -2.2, -2.2, -2.1, -2.1, -1.9],
    ),
}
#: Table 9 prints 1.8. The coefficient gives 1.97, and the -10 percent entry agrees with it.
TABLE_9_ERROR = {("CSL10_85", 0.5, 1.1)}


def test_table_9_confirms_the_exponents_but_for_its_error(lib):
    """Each is printed to 0.1, from a coefficient the report rounds, so allow 0.16."""
    for code, rows in TABLE_9.items():
        for factor, row in zip((1.1, 0.9), rows):
            for aep, printed in zip(AEPS, row):
                eq = lib.equation("1", aep)
                b = eq.coefficients[[v.code for v in eq.variables].index(code)]
                if code == "WETLAND":
                    change = 10 ** (b * 4.86 * (factor - 1)) - 1
                else:
                    change = factor**b - 1
                if (code, aep, factor) in TABLE_9_ERROR:
                    assert (printed, round(100 * change, 2)) == (1.8, 1.97)
                else:
                    assert 100 * change == pytest.approx(printed, abs=0.16)


def test_every_equation_equals_the_nss_snapshot(lib):
    """NSS orders slope before wetland; the report's order is the matrix's, so match by code."""
    snap = json.loads((DATA / "nss_snapshot" / "NH.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == s["intercept"]
        assert e.sep_log == pytest.approx(s["sep_log"], abs=2e-4)
        nss = {
            v["code"]: (c, v["transform"], v["minimum"], v["maximum"])
            for v, c in zip(s["variables"], s["coefficients"])
        }
        ours = {
            v.code: (c, v.transform, v.minimum, v.maximum)
            for v, c in zip(e.variables, e.coefficients)
        }
        assert ours == nss


def test_drainage_area_only_uses_the_average_error(lib):
    eq = lib.equation("2", 0.01)
    assert eq.covariance is None
    est = evaluate(eq, {"DRNAREA": 22.0})
    assert est.flow_cfs == pytest.approx(219 * 22.0**0.800)
    assert est.prediction_variance == pytest.approx(0.250**2)


def test_out_of_range_input_raises(lib):
    with pytest.raises(ValueError):
        evaluate(
            lib.equation("1", 0.01),
            {"DRNAREA": 22.0, "APRAVPRE": 7.0, "WETLAND": 3.27, "CSL10_85": 62.6},
        )
