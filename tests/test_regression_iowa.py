"""Iowa's verified equations: SIR 2013-5086 (full sets) and SIR 2015-5055 (drainage area only)."""

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
#: Report tables 3, 4 and 14 at the 394 gages regressed; SIR 2015-5055 tables 1 and 8 at its 80.
GAGES = json.loads((ROOT / "fixtures" / "regression" / "iowa_gages.json").read_text("utf-8"))
FULL = {"1": "GC1561", "2": "GC1562", "3": "GC1563"}
DA_ONLY = {"1": "GC1564", "2": "GC1627", "3": "GC1628"}
AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]

#: Table 3 prints a characteristic that Tables 4 and 14 do not use: (code, printed, the value
#: that reproduces both at all 8 AEPs). The five hold at most three decimals; other gages more.
TABLE_3_VALUE_NOT_USED = {
    "05448600": ("CCM", 1.452, 0.695),
    "0545129280": ("CCM", 0.89, 0.492),
    "0548065350": ("CCM", 0.568, 0.377),
    "05481528": ("CCM", 0.456, 0.282),
    "05460500": ("BSHAPE", 3.875, 3.08),
}
#: Table 4 cells off by more than rounding, with the equation's value. Cause unknown: the
#: gage's other six AEPs reproduce, and its WIE row agrees with the printed cells.
TABLE_4_CELL_DIFFERS = {
    ("05448500", 0.005): (6110.0, 6161.5),
    ("05448500", 0.002): (7260.0, 7206.7),
}
#: Standard error of the model, percent, Tables 9-11. Table 14 was computed from this.
SEM_PERCENT = {
    "1": [39.4, 30.4, 29.4, 30.4, 32.6, 34.7, 37.5, 41.5],
    "2": [45.7, 24.4, 19.3, 17.6, 18.5, 20.3, 22.9, 26.1],
    "3": [41.6, 29.0, 25.6, 24.8, 26.0, 27.2, 28.4, 31.5],
}
#: Table 15 of SIR 2013-5086: exponent, and the constant in log units.
TABLE_15 = {
    "1": [
        (0.641, 1.66),
        (0.592, 2.06),
        (0.569, 2.26),
        (0.548, 2.45),
        (0.535, 2.57),
        (0.524, 2.67),
        (0.516, 2.75),
        (0.506, 2.85),
    ],
    "2": [
        (0.579, 2.15),
        (0.525, 2.57),
        (0.499, 2.78),
        (0.476, 2.97),
        (0.463, 3.08),
        (0.453, 3.18),
        (0.445, 3.26),
        (0.436, 3.35),
    ],
    "3": [
        (0.579, 2.31),
        (0.528, 2.69),
        (0.503, 2.87),
        (0.479, 3.05),
        (0.466, 3.16),
        (0.455, 3.25),
        (0.447, 3.33),
        (0.438, 3.41),
    ],
}
#: Where NSS's limits differ from Table 12: (region, code): (minimum, maximum) as NSS serves them.
NSS_RANGE_DIFFERS = {
    ("GC1561", "DRNAREA"): (0.06, 5464.0),
    ("GC1562", "DRNAREA"): (0.08, 7783.0),
    ("GC1562", "BSHAPE"): (0.806, 13.94),
    ("GC1563", "DRNAREA"): (0.05, 2809.0),
    ("GC1564", "DRNAREA"): (0.06, 5464.0),
    ("GC1627", "DRNAREA"): (0.08, 7783.0),
    ("GC1628", "DRNAREA"): (0.05, 2809.0),
}


@pytest.fixture(scope="module")
def lib():
    return load_state("IA")


def xux(eq, characteristics):
    x = np.array([1.0] + [v.value(characteristics) for v in eq.variables])
    return float(x @ np.array(eq.covariance) @ x)


def sem_variance(percent):
    return (math.sqrt(math.log(1 + (percent / 100) ** 2)) / math.log(10)) ** 2


def table_3_as_used(site):
    gage = dict(GAGES["full_sets"][site])
    if site in TABLE_3_VALUE_NOT_USED:
        code, printed, used = TABLE_3_VALUE_NOT_USED[site]
        assert gage[code] == printed
        gage[code] = used
    return gage


def test_iowa_is_verified(lib):
    """Everything is stored and cross-checked; TABLE_4_CELL_DIFFERS has no known cause."""
    assert lib.status == "verified"
    assert "05448500" in lib.notes
    assert lib.regions == [*FULL.values(), *DA_ONLY.values()]
    assert len(lib.equations) == 48
    for code in lib.regions:
        assert region_for_nss_code(lib, code) == code
        assert f"(NSS {code})" in lib.equation(code, 0.01).region_name
        assert set(lib.aeps(code)) == set(AEPS)


def test_the_region_preference_is_the_snapshots(lib):
    snap = json.loads((DATA / "nss_snapshot" / "IA.json").read_text(encoding="utf-8"))
    assert lib.region_preference.tiers == snap["region_preference"]["tiers"]
    assert lib.region_preference.basis == snap["region_preference"]["basis"]
    assert lib.region_preference.tiers == [list(FULL.values()), list(DA_ONLY.values())]


def test_example_1_north_raccoon_river(lib):
    """05482500, region 1, the 0.2-percent AEP; the interval printed is 17,000 to 68,600 ft3/s."""
    inputs = {"DRNAREA": 1609.35, "I24H10Y": 4.321, "CCM": 1.067}
    eq = lib.equation("GC1561", 0.002)
    est = evaluate(eq, inputs)
    assert round(est.flow_cfs, -2) == 34100
    assert eq.model_error_variance == pytest.approx(0.029998, abs=5e-7)
    # Printed 0.00327303, which takes CCM**0.55 as the 1.036 of the report's arithmetic.
    x = np.array([1.0, math.log10(1609.35), 4.321, 1.036])
    assert float(x @ np.array(eq.covariance) @ x) == pytest.approx(0.00327303, abs=5e-9)
    assert xux(eq, inputs) == pytest.approx(0.00327263, abs=5e-9)
    assert math.sqrt(est.prediction_variance) == pytest.approx(0.182403, abs=2e-6)
    t = 10 ** (stats.t.ppf(0.95, 87) * math.sqrt(est.prediction_variance))
    assert t == pytest.approx(2.0103, abs=1e-4)
    assert round(34100 / t, -2) == 17000 and round(34100 * t, -2) == 68600
    assert est.interval[1] / est.flow_cfs == pytest.approx(t, rel=1e-9)


def test_example_2_english_river(lib):
    """05455500, region 2, the 1-percent AEP; the interval printed is 19,500 to 39,400 ft3/s."""
    inputs = {"DRNAREA": 574.10, "DESMOIN": 0.0, "BSHAPE": 6.155}
    eq = lib.equation("GC1562", 0.01)
    est = evaluate(eq, inputs)
    assert round(est.flow_cfs, -2) == 27700
    assert eq.model_error_variance == pytest.approx(0.007617, abs=5e-7)
    assert xux(eq, inputs) == pytest.approx(0.000902081, abs=5e-10)
    assert math.sqrt(est.prediction_variance) == pytest.approx(0.0922989, abs=5e-6)
    t = 10 ** (stats.t.ppf(0.95, 172) * math.sqrt(est.prediction_variance))
    assert t == pytest.approx(1.4212, abs=1e-4)
    assert round(27700 / t, -2) == 19500 and round(27700 * t, -2) == 39400


def test_the_examples_covariance_tables_match_the_workbook(lib):
    """The report prints two matrices in its text, to nine decimals, and the rest only in xlsx."""
    example_1 = [
        [0.579267620, -0.000733524, -0.141069630, 0.014150952],
        [-0.000733524, 0.000744762, 0.000042403, -0.001270299],
        [-0.141069630, 0.000042403, 0.035474003, -0.007395614],
        [0.014150952, -0.001270299, -0.007395614, 0.020494642],
    ]
    example_2 = [
        [0.041672405, -0.045784820, -0.000039051, -0.000518851],
        [-0.045784820, 0.052558399, 0.000038686, 0.000436101],
        [-0.000039051, 0.000038686, 0.000000287, 0.000000161],
        [-0.000518851, 0.000436101, 0.000000161, 0.000025378],
    ]
    assert np.array(lib.equation("GC1561", 0.002).covariance) == pytest.approx(
        np.array(example_1), abs=6e-10
    )
    assert np.array(lib.equation("GC1562", 0.01).covariance) == pytest.approx(
        np.array(example_2), abs=6e-10
    )
    for code in FULL.values():
        for aep in AEPS:
            eq = lib.equation(code, aep)
            u = np.array(eq.covariance)
            assert u.shape == (len(eq.variables) + 1, len(eq.variables) + 1)
            assert np.array_equal(u, u.T)


def test_t_values_follow_from_the_stored_site_counts(lib):
    """Table 13's t for each region; region 2's 50-percent equation has one variable fewer."""
    printed = {("GC1561", 0.01): 1.6626, ("GC1562", 0.5): 1.6537, ("GC1562", 0.01): 1.6538}
    for (code, aep), t in printed.items():
        eq = lib.equation(code, aep)
        dof = eq.n_sites - len(eq.variables) - 1
        assert stats.t.ppf(0.95, dof) == pytest.approx(t, abs=5e-5)
    assert [lib.equation(c, 0.01).n_sites for c in FULL.values()] == [91, 176, 127]
    assert [lib.equation(c, 0.01).n_sites for c in DA_ONLY.values()] == [91, 176, 127]


def test_regression_flows_reproduce_at_every_gage(lib):
    """Table 4 prints three significant figures, so at most about 0.5 percent."""
    assert len(GAGES["full_sets"]) == 394
    seen = 0
    for site in GAGES["full_sets"]:
        gage = table_3_as_used(site)
        for i, aep in enumerate(AEPS):
            est = evaluate(lib.equation(FULL[gage["region"]], aep), gage, allow_extrapolation=True)
            if (site, aep) in TABLE_4_CELL_DIFFERS:
                printed, computed = TABLE_4_CELL_DIFFERS[(site, aep)]
                assert gage["flows"][i] == printed
                assert est.flow_cfs == pytest.approx(computed, abs=0.1)
                continue
            assert est.flow_cfs == pytest.approx(gage["flows"][i], rel=0.005)
            seen += 1
    assert seen == 394 * 8 - 2


@pytest.mark.parametrize("site", sorted(TABLE_3_VALUE_NOT_USED))
def test_table_3s_printed_value_does_not_give_table_4(lib, site):
    """The pinned values are needed: Table 3's own misses by 2.7 percent or more."""
    gage = GAGES["full_sets"][site]
    for i, aep in enumerate(AEPS[1:], start=1):
        est = evaluate(lib.equation(FULL[gage["region"]], aep), gage, allow_extrapolation=True)
        assert abs(est.flow_cfs / gage["flows"][i] - 1) > 0.027


def test_table_14_variances_confirm_the_covariance_and_its_order(lib):
    """Table 14 prints to 1e-4. It used the MEV that SEM in percent gives, not Table 13's."""
    for site in GAGES["full_sets"]:
        gage = table_3_as_used(site)
        region = gage["region"]
        for i, aep in enumerate(AEPS):
            eq = lib.equation(FULL[region], aep)
            total = sem_variance(SEM_PERCENT[region][i]) + xux(eq, gage)
            assert total == pytest.approx(gage["rre_variance"][i], abs=1e-4)


def test_table_13_mev_is_table_9_to_11s_sem(lib):
    for region, sems in SEM_PERCENT.items():
        for sem, aep in zip(sems, AEPS):
            mev = lib.equation(FULL[region], aep).model_error_variance
            assert mev == pytest.approx(sem_variance(sem), abs=7e-5)


def test_average_prediction_variance_is_the_printed_avp(lib):
    """Mean(MEV + xUx') over a region's gages: the check on the matrix NSS cannot give."""
    for region, code in FULL.items():
        sites = [s for s, g in GAGES["full_sets"].items() if g["region"] == region]
        for aep in AEPS:
            eq = lib.equation(code, aep)
            mean = np.mean([eq.model_error_variance + xux(eq, table_3_as_used(s)) for s in sites])
            assert mean == pytest.approx(eq.avp, abs=5.5e-4)
            assert eq.sep_log**2 == pytest.approx(mean, abs=1.5e-4)


def test_table_12_ranges_are_those_of_table_3():
    """Each limit is its gage extreme rounded as Table 12 prints it, except CCM's maximum.

    Table 12 prints 2.22. SIR 2015-5055 p. 11 says the 2013 equations cannot be used above 3.87,
    which is what NSS serves, so 3.87 is stored.
    """
    gages = GAGES["full_sets"].values()
    ccm = [g["CCM"] for g in gages if g["region"] == "1"]
    assert max(ccm) == pytest.approx(2.2155, abs=1e-4)
    area = {r: max(g["DRNAREA"] for g in gages if g["region"] == r) for r in FULL}
    assert area == pytest.approx({"1": 5463.88, "2": 7782.62, "3": 2809.05}, abs=0.005)


def test_drainage_area_only_flows_reproduce_at_every_gage(lib):
    """SIR 2015-5055 Table 8, from the drainage areas of its Table 1, rounded to StreamStats'."""
    assert len(GAGES["drainage_area_only"]) == 80
    for gage in GAGES["drainage_area_only"].values():
        for i, aep in enumerate(AEPS):
            est = evaluate(lib.equation(DA_ONLY[gage["region"]], aep), gage)
            assert est.flow_cfs == pytest.approx(gage["flows"][i], rel=0.005)


def test_drainage_area_only_equals_table_15_of_the_first_report(lib):
    """The second place these are printed: same exponents, constants in log units to 2 places."""
    for region, rows in TABLE_15.items():
        for (exponent, log_constant), aep in zip(rows, AEPS):
            eq = lib.equation(DA_ONLY[region], aep)
            assert eq.coefficients == (exponent,)
            assert eq.intercept == pytest.approx(log_constant, abs=0.0055)
            assert eq.sep_log**2 == pytest.approx(eq.avp, abs=5.5e-4)


def test_every_equation_equals_the_nss_snapshot_but_for_its_ranges(lib):
    snap = json.loads((DATA / "nss_snapshot" / "IA.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(e.region_code, e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(e.region_code, e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-12)
        assert list(e.coefficients) == s["coefficients"]
        assert e.sep_log == pytest.approx(s["sep_log"], abs=1e-5)
        assert [v.code for v in e.variables] == [n["code"] for n in s["variables"]]
        for v, n in zip(e.variables, s["variables"]):
            assert v.exponent == n.get("exponent")
            if n["transform"] == "log10_plus1":
                # The same term: NSS's alias against offset=1 here.
                assert (v.transform, v.offset) == ("log10", 1.0)
            else:
                assert (v.transform, v.offset) == (n["transform"], 0.0)
            pair = (n["minimum"], n["maximum"])
            if (e.region_code, v.code) in NSS_RANGE_DIFFERS:
                assert pair == NSS_RANGE_DIFFERS[(e.region_code, v.code)]
                assert (v.minimum, v.maximum) != pair
            else:
                assert (v.minimum, v.maximum) == pair
