"""Michigan's verified equations (WRI 84-4207, Tables 4-7, 10 and 12-15)."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
NSS = {"1": "GC106", "2": "GC107", "3": "GC108", "4": "GC109", "5": "GC110"}
AEPS = (0.2, 0.1, 0.04, 0.02, 0.01)
#: Characteristics (Tables 5, 6) and regression flows (Table 4, lower lines) at the 185 peak stations.
GAGES = json.loads((ROOT / "fixtures" / "regression" / "michigan_gages.json").read_text("utf-8"))

#: Table 10: constant a, mean-square error, percent error as printed.
TABLE_10 = {
    0.2: (0.6869, 0.01724, 30.1),
    0.1: (0.6688, 0.01908, 32.6),
    0.04: (0.6099, 0.02167, 34.9),
    0.02: (0.5569, 0.02388, 36.7),
    0.01: (0.4936, 0.02620, 38.6),
}
#: Table 7, peak-flow rows by region.
AREAL = {
    "1": (0.9890, 0.9840, 0.9790, 0.9761, 0.9741),
    "2": (1.0209, 1.0238, 1.0301, 1.0354, 1.0426),
    "3": (0.9986, 0.9940, 0.9886, 0.9849, 0.9820),
    "4": (1.0705, 1.0639, 1.0546, 1.0479, 1.0409),
    "5": (0.9171, 0.9303, 0.9447, 0.9539, 0.9605),
}
#: Table 14's peak-flow rows, lower triangle, in its printed order.
TABLE_14_ORDER = (
    "CONTDA",
    "CSL10_85",
    "PLCHSWAMP",
    "OUTWASH",
    "MEDTILL",
    "MORTILL",
    "MUCK",
    "I24H100Y",
    "PCLAY",
    "SLENRAT",
    "TILROCK",
    "CRSTILL",
)
TABLE_14 = (
    (0.05000,),
    (0.04720, 0.1146),
    (-0.005842, -0.001666, 0.03075),
    (-0.009858, -0.0006025, -0.009029, 0.02714),
    (-0.001384, -0.001144, 0.0006744, 0.003859, 0.01284),
    (-0.008871, -0.01554, 0.003113, -0.0008675, 0.001725, 0.03139),
    (-0.007026, 0.004338, -0.002773, 0.008054, 0.005425, 0.002767, 0.06641),
    (0.2367, 0.2304, 0.003678, -0.03797, -0.002807, 0.07515, 0.2009, 10.31),
    (-0.007905, -0.003424, 0.007635, 0.01114, 0.002051, -0.003696, 0.003104, 0.03242, 0.02776),
    (
        -0.007331,
        0.01359,
        -0.01008,
        -0.004133,
        -0.001772,
        0.003870,
        0.008718,
        0.03659,
        -0.01392,
        0.1474,
    ),
    (
        -0.01270,
        -0.02541,
        -0.004047,
        0.01018,
        0.006140,
        0.005099,
        0.0006057,
        -0.08560,
        0.005639,
        -0.002631,
        0.03441,
    ),
    (
        -0.01668,
        -0.01551,
        -0.003168,
        0.002737,
        0.003130,
        0.006707,
        -0.002471,
        -0.01599,
        0.005192,
        0.0007922,
        0.001816,
        0.02693,
    ),
)
#: Table 13's peak-flow rows: mean of logs, sd of logs, minimum and maximum of logs.
TABLE_13 = {
    "CONTDA": (1.934, 0.6143, 0.06070, 2.981),
    "CSL10_85": (0.8634, 0.3290, -0.2757, 1.757),
    "PLCHSWAMP": (1.302, 0.5981, 0.0, 1.987),
    "SLENRAT": (0.8404, 0.2106, 0.3682, 1.362),
    "I24H100Y": (0.6380, 0.02730, 0.5611, 0.7126),
    "OUTWASH": (1.107, 0.6593, 0.0, 2.004),
    "MUCK": (0.1189, 0.3285, 0.0, 1.630),
    "MORTILL": (0.2312, 0.4757, 0.0, 1.641),
    "MEDTILL": (0.9350, 0.7434, 0.0, 2.004),
    "PCLAY": (0.3484, 0.6011, 0.0, 1.868),
    "TILROCK": (0.1717, 0.5015, 0.0, 1.994),
    "CRSTILL": (0.5962, 0.6044, 0.0, 1.986),
}
#: Table 10's Prob>|T| rows, in Table 10's variable order.
PROB_T = {
    0.2: (1e-4, 1e-4, 1e-4, 0.0240, 0.0141, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 0.0648),
    0.1: (1e-4, 1e-4, 1e-4, 0.0210, 0.0078, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 0.0006, 0.0520),
    0.04: (1e-4, 1e-4, 1e-4, 0.0180, 0.0033, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 0.0027, 0.0440),
    0.02: (1e-4, 1e-4, 1e-4, 0.0180, 0.0019, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 0.0069, 0.0420),
    0.01: (1e-4, 1e-4, 1e-4, 0.0180, 0.0011, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 0.0150, 0.0440),
}
#: Table 6 prints no CLAY or OUTWASH here. These restore Table 13's means and the Table 4 flows.
BUTTERNUT = "04147990"
BUTTERNUT_INFERRED = {"PCLAY": 10.8, "OUTWASH": 5.7}
#: Stations beyond 1.2 percent, the largest difference at each.
PINNED_FLOW_DIFFERENCES = {BUTTERNUT: 0.079, "04148800": 0.015, "04148160": 0.0135}


@pytest.fixture(scope="module")
def lib():
    return load_state("MI")


def _table_14() -> np.ndarray:
    c = np.zeros((12, 12))
    for i, row in enumerate(TABLE_14):
        for j, value in enumerate(row):
            c[i, j] = c[j, i] = value
    return c


def _logs(gage: dict, eq) -> np.ndarray:
    return np.array([v.apply(gage[v.code]) for v in eq.variables])


def test_michigan_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == ["1", "2", "3", "4", "5"]
    assert len(lib.equations) == 25
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region
        assert lib.aeps(region) == sorted(AEPS, reverse=True) or set(lib.aeps(region)) == set(AEPS)


def test_constants_are_table_10_times_table_7(lib):
    for region, factors in AREAL.items():
        for aep, factor in zip(AEPS, factors):
            eq = lib.equation(region, aep)
            assert eq.intercept == pytest.approx(TABLE_10[aep][0] + math.log10(factor), abs=1e-12)
            assert eq.model_error_variance == TABLE_10[aep][1]
            assert eq.sep_log == pytest.approx(math.sqrt(TABLE_10[aep][1]), abs=1e-6)
            assert eq.n_sites == 185


def test_published_regression_flows_reproduce_at_every_station(lib):
    """Table 4 prints three significant figures, and Tables 5 and 6 about the same."""
    assert len(GAGES) == 185
    within = 0
    for station, gage in GAGES.items():
        worst = 0.0
        for aep, printed in zip(AEPS, gage["flows"]):
            eq = lib.equation(gage["region"], aep)
            inputs = {v.code: gage[v.code] for v in eq.variables}
            est = evaluate(eq, inputs, allow_extrapolation=True)
            worst = max(worst, abs(printed / est.flow_cfs - 1))
        assert worst < PINNED_FLOW_DIFFERENCES.get(station, 0.012), station
        within += worst < 0.012
    assert within == 182


def test_butternut_creek_is_short_two_cells_in_table_6(lib):
    """With the two cells Table 13 implies, it matches like its neighbours."""
    gage = GAGES[BUTTERNUT]
    assert gage["PCLAY"] == 0.0 and gage["OUTWASH"] == 0.0
    filled = {**gage, **BUTTERNUT_INFERRED}
    for aep, printed in zip(AEPS, gage["flows"]):
        eq = lib.equation(gage["region"], aep)
        est = evaluate(eq, {v.code: filled[v.code] for v in eq.variables}, allow_extrapolation=True)
        assert est.flow_cfs == pytest.approx(printed, rel=0.016)


def test_station_table_reproduces_table_13_and_table_14(lib):
    """The only check on Table 14's labels. It holds element by element in the printed order."""
    filled = {k: ({**g, **BUTTERNUT_INFERRED} if k == BUTTERNUT else g) for k, g in GAGES.items()}
    eq = lib.equation("1", 0.01)
    position = {v.code: i for i, v in enumerate(eq.variables)}
    logs = np.array([_logs(g, eq) for g in filled.values()])
    for code, (mean, sd, _, _) in TABLE_13.items():
        column = logs[:, position[code]]
        assert column.mean() == pytest.approx(mean, abs=1e-3), code
        assert column.std(ddof=1) == pytest.approx(sd, abs=6e-4), code
    order = [position[c] for c in TABLE_14_ORDER]
    centred = logs[:, order] - logs[:, order].mean(axis=0)
    from_stations = np.linalg.inv(centred.T @ centred)
    printed = _table_14()
    scale = np.sqrt(np.outer(np.diag(printed), np.diag(printed)))
    assert np.abs(from_stations - printed).max() < 0.002
    assert (np.abs(from_stations - printed) / scale).max() < 0.02
    leverage = sum(1 / 185 + x @ printed @ x for x in centred)
    assert leverage == pytest.approx(13, abs=0.03)


def test_covariance_is_table_14_uncentred(lib):
    """Model error + x'Ux is the report's S**2 (1 + 1/n + X'CX)."""
    printed = _table_14()
    gage = GAGES["04001000"]
    for region in NSS:
        for aep in AEPS:
            eq = lib.equation(region, aep)
            codes = [v.code for v in eq.variables]
            order = [codes.index(c) for c in TABLE_14_ORDER]
            means = np.array([TABLE_13[c][0] for c in TABLE_14_ORDER])
            x = _logs(gage, eq)[order] - means
            s2 = eq.model_error_variance
            expected = s2 * (1 + 1 / 185 + x @ printed @ x)
            row = np.concatenate([[1.0], _logs(gage, eq)])
            actual = s2 + row @ np.asarray(eq.covariance) @ row
            assert actual == pytest.approx(expected, rel=2e-6)
            est = evaluate(eq, {c: gage[c] for c in codes}, allow_extrapolation=True)
            factor = 10 ** (stats.t.ppf(0.95, 185 - 13) * math.sqrt(expected))
            assert est.interval[1] / est.flow_cfs == pytest.approx(factor, rel=1e-5)
            assert est.interval is not None and "covariance" in est.interval_method


def test_table_14_diagonal_reproduces_prob_t(lib):
    """t = b / sqrt(S**2 Cii) on 172 degrees of freedom. It also settles S**2 at P5."""
    printed = np.diag(_table_14())
    for aep, probs in PROB_T.items():
        eq = lib.equation("1", aep)
        for v, b, prob in zip(eq.variables, eq.coefficients, probs):
            c = printed[TABLE_14_ORDER.index(v.code)]
            p = 2 * stats.t.sf(abs(b) / math.sqrt(eq.model_error_variance * c), 172)
            if prob == 1e-4:
                assert p < 1.5e-4, (aep, v.code)
            else:
                assert p == pytest.approx(prob, abs=8e-4), (aep, v.code)


def test_percent_error_at_p5_is_a_misprint():
    """Table 10's 0.01724 gives 30.9; it prints 30.1. Its other rows all agree."""
    pct = {a: 100 * math.sqrt(10 ** (math.log(10) * s) - 1) for a, (_, s, _) in TABLE_10.items()}
    for aep, (_, _, printed) in TABLE_10.items():
        if aep == 0.2:
            assert round(pct[aep], 1) == 30.9 and printed == 30.1
        else:
            assert round(pct[aep], 1) == printed


def test_worked_example_is_mean_flow_only():
    """Table 15, Washington Creek: not an AEP statistic, so only the variance formula is checked."""
    log_q = (
        -1.4813
        + 1.0217 * math.log10(13.2)
        + 0.5385 * math.log10(90)
        + 1.1703 * math.log10(2.10)
        + math.log10(1.0097)
    )
    assert 10**log_q == pytest.approx(12.5, abs=0.05)
    c = np.array(
        [[0.02757, -0.006255, 0.1748], [-0.006255, 0.2080, 0.3430], [0.1748, 0.3430, 10.62]]
    )
    x = np.array([1.1206 - 2.1073, 1.9542 - 1.7629, 0.3222 - 0.3317])
    xcx = float(x @ c @ x)
    assert xcx == pytest.approx(0.03980, abs=5e-5)
    assert math.sqrt(0.008549 * (1 / 146 + xcx)) == pytest.approx(0.01997, abs=5e-5)
    sd = math.sqrt(0.008549 * (1 + 1 / 146 + xcx))
    assert sd == pytest.approx(0.09459, abs=5e-5)
    assert 10 ** (log_q - 1.96 * 0.01997) == pytest.approx(11.4, abs=0.05)
    assert 10 ** (log_q + 1.96 * 0.01997) == pytest.approx(13.7, abs=0.05)
    assert 10 ** (log_q - 1.96 * sd) == pytest.approx(8.16, abs=0.03)
    assert 10 ** (log_q + 1.96 * sd) == pytest.approx(19.2, abs=0.05)


def test_ranges_are_table_13_as_printed(lib):
    eq = lib.equation("1", 0.01)
    for v in eq.variables:
        low, high = TABLE_13[v.code][2:]
        offset = v.offset or 0.0
        assert math.log10(v.minimum + offset) == pytest.approx(low, abs=1e-3), v.code
        assert math.log10(v.maximum + offset) == pytest.approx(high, abs=1e-3), v.code


def test_table_13_peak_range_is_not_the_peak_stations(lib):
    """Its extremes repeat the flood-volume row, whose sample is smaller and differs."""
    high = {"CSL10_85": 90.5, "MORTILL": 56.1, "PCLAY": 100.0}
    for code, value in high.items():
        assert max(g[code] for g in GAGES.values()) == value
        assert lib.equation("1", 0.01).variables[_index(lib, code)].maximum < value
    assert max(g["I24H100Y"] for g in GAGES.values()) == 5.08
    assert lib.equation("1", 0.01).variables[_index(lib, "I24H100Y")].maximum == 5.16


def _index(lib, code: str) -> int:
    return [v.code for v in lib.equation("1", 0.01).variables].index(code)


def test_every_equation_equals_the_nss_snapshot_but_for_rounding(lib):
    """Every coefficient agrees. NSS rounds 10**a * AREAL to 3 decimals and each range outward."""
    snap = json.loads((DATA / "nss_snapshot" / "MI.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert list(e.coefficients) == s["coefficients"]
        assert 10 ** s["intercept"] == pytest.approx(10**e.intercept, abs=5.1e-4)
        assert s.get("sep_log") is None
        for v, n in zip(e.variables, s["variables"]):
            assert v.code == n["code"]
            assert n["minimum"] <= v.minimum and v.maximum <= n["maximum"]
            assert n["maximum"] != v.maximum or v.maximum == 100.0
