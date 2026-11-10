"""Kentucky's verified equations (WRI 03-4180, Tables 1-5, and WRI 97-4219, Tables 8, 10, 12-14)."""

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
RURAL_AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
URBAN = "Jefferson County urban"
URBAN_AEPS = RURAL_AEPS[:6]
NSS = {f"Region {i + 1}": f"GC{524 + i}" for i in range(7)}
NSS[URBAN] = "GC531"
SITES = {"Region 1": 28, "Region 2": 68, "Region 3": 24, "Region 4": 17}
SITES.update({"Region 5": 40, "Region 6": 27, "Region 7": 34, URBAN: 17})
FIXTURE = json.loads((ROOT / "fixtures" / "regression" / "kentucky_gages.json").read_text("utf-8"))
RURAL = FIXTURE["rural"]
JEFFERSON = FIXTURE["jefferson"]
#: Table 1's regression row at these basins is an area-weighted mix of two regions. Table 2
#: lists one region and the report prints no fractions, so each one here is inferred from
#: the 8 values: gage -> (region A, region B, fraction in A).
TWO_REGION_GAGES = {
    "03277500": ("Region 2", "Region 3", 0.71),
    "03280000": ("Region 2", "Region 3", 0.68),
    "03310300": ("Region 5", "Region 6", 0.82),
    "03310500": ("Region 5", "Region 6", 0.49),
    "03406500": ("Region 3", "Region 4", 0.141),
}
#: Table 14 prints 1,150 at site RB2, 10-year; Table 10 gives 1,046.
TABLE_14_MISPRINT = ("RB2", 0.1)
#: Table 5's drainage-area-only exponents where Table 3 uses slope. These do not match by design.
TABLE_5_DA_ONLY = {
    ("Region 1", 0.1): 0.642,
    ("Region 1", 0.04): 0.634,
    ("Region 1", 0.02): 0.629,
    ("Region 1", 0.01): 0.625,
    ("Region 1", 0.005): 0.622,
    ("Region 1", 0.002): 0.618,
    ("Region 4", 0.5): 0.824,
    ("Region 4", 0.2): 0.803,
    ("Region 4", 0.1): 0.794,
    ("Region 4", 0.04): 0.786,
    ("Region 4", 0.02): 0.783,
}
#: Every other exponent in Table 5, by region and AEP order.
TABLE_5 = {
    "Region 1": [0.673, 0.651],
    "Region 2": [0.728, 0.721, 0.715, 0.709, 0.704, 0.699, 0.695, 0.690],
    "Region 3": [0.748, 0.712, 0.692, 0.670, 0.656, 0.643, 0.632, 0.620],
    "Region 4": [None, None, None, None, None, 0.780, 0.778, 0.776],
    "Region 5": [0.704, 0.692, 0.686, 0.682, 0.679, 0.677, 0.676, 0.674],
    "Region 6": [0.600, 0.586, 0.578, 0.569, 0.564, 0.559, 0.555, 0.551],
    "Region 7": [0.623, 0.616, 0.613, 0.610, 0.610, 0.609, 0.610, 0.610],
}


@pytest.fixture(scope="module")
def lib():
    return load_state("KY")


def _flow(lib, region, aep, gage):
    eq = lib.equation(region, aep)
    inputs = {v.code: gage[v.code] for v in eq.variables}
    return evaluate(eq, inputs, allow_extrapolation=True).flow_cfs


def test_kentucky_is_verified(lib):
    assert lib.status == "verified"
    assert sorted(lib.regions) == sorted(NSS)
    assert len(lib.equations) == 62
    assert len(lib.source_reports) == 2
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region
    for e in lib.equations:
        assert e.n_sites == SITES[e.region_code]
    for region in NSS:
        assert lib.aeps(region) == (URBAN_AEPS if region == URBAN else RURAL_AEPS)
    assert sum(SITES[r] for r in SITES if r != URBAN) == 238


def test_worked_example_1_two_region_q100(lib):
    """WRI 03-4180, p. 27: 600 mi2, 333 mi2 of it (55.5 percent) in Region 2."""
    q2 = evaluate(lib.equation("Region 2", 0.01), {"DRNAREA": 600.0}).flow_cfs
    q3 = evaluate(lib.equation("Region 3", 0.01), {"DRNAREA": 600.0}).flow_cfs
    assert q2 == pytest.approx(47100, abs=50)
    assert q3 == pytest.approx(67300, abs=50)
    assert 0.555 * q2 + 0.445 * q3 == pytest.approx(56100, abs=50)


def test_worked_example_2_transfer_uses_the_stored_exponents(lib):
    """WRI 03-4180, p. 28: downstream of a 466 mi2 gage whose weighted Q100 is 50,700 ft3/s."""
    wr = 600 / 466 - 1
    assert round(wr, 3) == 0.288
    final = []
    for region, regression, upstream, printed in (
        ("Region 2", 47100, 60500, 56600),
        ("Region 3", 67300, 59600, 61800),
    ):
        b = lib.equation(region, 0.01).coefficients[0]
        qu = 50700 * (600 / 466) ** b
        assert qu == pytest.approx(upstream, abs=50)
        quf = regression * 0.288 + upstream * (1 - 0.288)
        assert quf == pytest.approx(printed, abs=50)
        final.append(printed)
    assert 0.555 * final[0] + 0.445 * final[1] == pytest.approx(58900, abs=50)


def test_worked_example_3_interpolation_between_two_gages():
    """WRI 03-4180, p. 28: 600 mi2 between gages of 466 and 1,101 mi2, Q100 50,700 and 70,000."""
    log_q = math.log10(50700) + (math.log10(70000) - math.log10(50700)) * (
        math.log10(600) - math.log10(466)
    ) / (math.log10(1101) - math.log10(466))
    assert round(log_q, 4) == 4.7462
    assert 10**log_q == pytest.approx(55700, abs=50)


def test_jefferson_county_worked_example(lib):
    """WRI 97-4219, p. 29, whose inputs are site FH1 of Table 8."""
    flow = evaluate(
        lib.equation(URBAN, 0.01), {"CONTDA": 1.66, "CSL10_85": 48.0, "BDF": 7.0}
    ).flow_cfs
    assert flow == pytest.approx(2460, abs=5)


def test_table_1_regression_flows_reproduce_at_single_region_gages(lib):
    """Table 1 prints three significant figures, as Table 2 does its areas and slopes."""
    assert len(RURAL) == 169
    cells = 0
    for number, gage in RURAL.items():
        if number in TWO_REGION_GAGES:
            continue
        for aep, printed in zip(RURAL_AEPS, gage["flows"]):
            assert _flow(lib, gage["region"], aep, gage) == pytest.approx(printed, rel=0.005)
            cells += 1
    assert cells == 164 * 8


@pytest.mark.parametrize("number", sorted(TWO_REGION_GAGES))
def test_table_1_is_an_area_weighted_mix_at_two_region_gages(lib, number):
    a, b, fraction = TWO_REGION_GAGES[number]
    gage = RURAL[number]
    assert gage["region"] in (a, b)
    flows = []
    for aep, printed in zip(RURAL_AEPS, gage["flows"]):
        flows.append(printed)
        mixed = fraction * _flow(lib, a, aep, gage) + (1 - fraction) * _flow(lib, b, aep, gage)
        assert mixed == pytest.approx(printed, rel=0.0035)
    # Computed from the gage's listed region alone, every value is at least 2 percent off.
    alone = [_flow(lib, gage["region"], aep, gage) for aep in RURAL_AEPS]
    assert np.min(np.abs(np.array(alone) / np.array(flows) - 1)) > 0.02


def test_two_region_basins_are_weighted_in_flow_by_default(lib):
    """03277500 at example 1's own split (333 of 466 mi2 in Region 2) gives Table 1's values."""
    from flowfreq import estimate_from_characteristics

    assert lib.weighting_space == "linear"
    gage = RURAL["03277500"]
    chars = {k: v for k, v in gage.items() if isinstance(v, (int, float))}
    share = 333 / 466 * 100
    est = estimate_from_characteristics(
        "KY", chars, region={"Region 2": share, "Region 3": 100 - share}, lat=37.6, lon=-83.3
    )
    q = est.quantiles.set_index("aep")["flow_cfs"]
    for aep, printed in zip(RURAL_AEPS, gage["flows"]):
        assert q.loc[aep] == pytest.approx(printed, rel=0.006)


def test_two_region_fraction_at_hazard_matches_the_worked_example():
    """Example 1 has 333 of 466 mi2 upstream in Region 2, the same gage's area and region mix."""
    assert RURAL["03277500"]["DRNAREA"] == 466.0
    assert TWO_REGION_GAGES["03277500"][2] == pytest.approx(333 / 466, abs=0.01)


def test_table_14_reproduces_from_table_8_but_for_one_misprint(lib):
    assert len(JEFFERSON) == 17
    for site, values in JEFFERSON.items():
        for aep, printed in zip(URBAN_AEPS, values["flows"]):
            flow = _flow(lib, URBAN, aep, values)
            if (site, aep) == TABLE_14_MISPRINT:
                assert printed == 1150
                assert flow == pytest.approx(1046, abs=1)
            else:
                assert flow == pytest.approx(printed, rel=0.004)


#: Table 12, percent change in flow for a percent change in each variable from the mean of the
#: 17 basins (BDF 5), in the printed order of the rows, +50, +25, +10, -10, -25, -50.
TABLE_12_CHANGES = [50, 25, 10, -10, -25, -50]
TABLE_12 = {
    0.5: {
        "A": [29.4, 15.2, 6.2, -6.5, -16.7, -35.6],
        "BDF": [13.4, 5.9, 2.2, -2.1, -4.8, -8.7],
        "SL": [5.3, 2.9, 1.2, -1.4, -3.6, -8.5],
    },
    0.1: {
        "A": [26.2, 13.6, 5.6, -5.9, -15.2, -32.8],
        "BDF": [9.2, 4.1, 1.5, -1.4, -3.4, -6.2],
        "SL": [10.3, 5.7, 2.3, -2.5, -6.7, -15.5],
    },
    0.01: {
        "A": [24.4, 12.8, 5.3, -5.5, -14.3, -31.1],
        "BDF": [7.0, 3.1, 1.2, -1.1, -2.6, -4.8],
        "SL": [13.4, 7.2, 3.0, -3.2, -8.5, -19.3],
    },
}


@pytest.mark.parametrize("aep", sorted(TABLE_12))
def test_table_12_sensitivities_confirm_the_coefficients(lib, aep):
    """A second place for the 2-, 10- and 100-year exponents, so the RB2 misprint is Table 14's."""
    eq = lib.equation(URBAN, aep)
    a, sl, bdf = eq.coefficients
    for change, da, db, ds in zip(
        TABLE_12_CHANGES, TABLE_12[aep]["A"], TABLE_12[aep]["BDF"], TABLE_12[aep]["SL"]
    ):
        factor = 1 + change / 100
        # Printed to 0.1: 46 of the 54 cells round correctly and 8 are off in the last digit,
        # the largest being 10-year SL, +25: 5.57 against 5.7.
        assert 100 * (factor**a - 1) == pytest.approx(da, abs=0.15)
        assert 100 * (factor**sl - 1) == pytest.approx(ds, abs=0.15)
        assert 100 * (((13 - 5 * factor) / 8) ** bdf - 1) == pytest.approx(db, abs=0.15)


def test_sep_and_avp_of_table_3_agree(lib):
    """The AVP is model error + sampling error (p. 26), and the ASEP is its square root."""
    rural = [e for e in lib.equations if e.region_code != URBAN]
    assert len(rural) == 56
    for e in rural:
        assert math.sqrt(e.avp) == pytest.approx(e.sep_log, abs=2e-4)
        assert e.model_error_variance < e.avp
    assert lib.equation("Region 1", 0.02).sep_log == pytest.approx(math.log10(1.529), abs=1e-6)


@pytest.mark.parametrize("log_units,percent", [(0.1688, 40.4), (0.1458, 34.6), (0.1295, 30.4)])
def test_table_11_shows_the_form_of_the_jefferson_percent(log_units, percent):
    """Table 11 prints an RMSE in both units, so Table 10's percents are the exponential form."""
    value = 100 * math.sqrt(math.exp((math.log(10) * log_units) ** 2) - 1)
    assert value == pytest.approx(percent, abs=0.1)


def test_jefferson_stores_the_printed_percent(lib):
    for aep, percent in zip(URBAN_AEPS, [45.4, 40.2, 37.6, 35.4, 34.4, 33.8]):
        eq = lib.equation(URBAN, aep)
        expected = math.sqrt(math.log(1 + (percent / 100) ** 2)) / math.log(10)
        assert eq.sep_log == pytest.approx(expected, abs=1e-5)
        assert eq.avp is None and eq.model_error_variance is None


def test_table_5_confirms_the_drainage_area_exponents(lib):
    for region, row in TABLE_5.items():
        for aep, exponent in zip(RURAL_AEPS, row):
            if exponent is not None:
                eq = lib.equation(region, aep)
                assert [v.code for v in eq.variables] == ["DRNAREA"]
                assert eq.coefficients[0] == exponent
    for (region, aep), exponent in TABLE_5_DA_ONLY.items():
        eq = lib.equation(region, aep)
        assert [v.code for v in eq.variables] == ["DRNAREA", "CSL10_85"]
        assert eq.coefficients[0] > exponent


def test_jefferson_prints_13_minus_bdf(lib):
    for e in (lib.equation(URBAN, aep) for aep in URBAN_AEPS):
        bdf = e.variables[2]
        assert (bdf.code, bdf.transform, bdf.scale, bdf.offset) == ("BDF", "log10", -1.0, 13.0)
        assert (bdf.minimum, bdf.maximum) == (0.0, 7.0)


def test_every_equation_equals_the_nss_snapshot(lib):
    snap = json.loads((DATA / "nss_snapshot" / "KY.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    assert len(by_key) == 62
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-12)
        assert e.coefficients == pytest.approx(s["coefficients"], abs=1e-12)
        for v, n in zip(e.variables, s["variables"]):
            assert v.code == n["code"]
            assert v.minimum == n["minimum"]
            assert v.maximum == n["maximum"]
        if e.region_code == URBAN:
            assert s.get("sep_log") is None  # NSS has no SEP for this set; Table 10 prints it.
        else:
            assert e.sep_log == pytest.approx(s["sep_log"], abs=2.5e-4)


def test_the_abstract_figures_are_not_table_3s(lib):
    """The abstract's extremes (+15.0 in Region 3, +60.3 in Region 5) appear in no Table 3 row."""
    percents = {
        region: [
            round(100 * (10 ** lib.equation(region, aep).sep_log - 1), 1) for aep in RURAL_AEPS
        ]
        for region in ("Region 3", "Region 5")
    }
    assert 15.0 not in percents["Region 3"] and min(percents["Region 3"]) == 20.4
    assert 60.3 not in percents["Region 5"] and max(percents["Region 5"]) == 76.4


def test_notes_record_the_currency_decision(lib):
    assert "SIR 2026-5036" in lib.notes
    assert "RB2" in lib.notes
    assert "figure 6" in lib.notes
