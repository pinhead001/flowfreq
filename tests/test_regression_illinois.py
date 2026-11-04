"""Illinois' verified equations (ICT report 23-019 and USGS data release P9XPWUMI)."""

from __future__ import annotations

import json
import math
import statistics
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

from flowfreq.regression import evaluate
from flowfreq.regression.equations import RegressionEquation
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
AEPS = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)
NSS = {
    "1": "GC2001",
    "2": "GC2002",
    "3": "GC2003",
    "4": "GC2004",
    "5": "GC2005",
    "6": "GC2006",
    "7": "GC2007",
    "1-urban": "GC2008",
    "3-urban": "GC2009",
    "4-urban": "GC2010",
    "5-urban": "GC2011",
    "6-urban": "GC2012",
    "7-urban": "GC2013",
}
#: Characteristics (DR-5) with regression flows (DR-10) and variances (DR-11) at the valid gages.
GAGES = json.loads((ROOT / "fixtures" / "regression" / "illinois_gages.json").read_text("utf-8"))
#: Table 11: (bu)p and its standard error.
BU = dict(zip(AEPS, (0.487, 0.405, 0.391, 0.359, 0.331, 0.312, 0.301, 0.294)))
BU_SE = dict(zip(AEPS, (0.075, 0.080, 0.078, 0.085, 0.086, 0.134, 0.112, 0.115)))
#: Table 12: baseline urban fractions, by region, with the number of gages.
TABLE_12 = {
    "1": (35, 0.005952),
    "2": (119, 0.3029),
    "3": (60, 0.01037),
    "4": (36, 0.01033),
    "5": (36, 0.008841),
    "6": (16, 0.002571),
    "7": (9, 0.0),
}
#: Table 6: intercept and coefficients at the 1-percent AEP, keyed in the report's order.
TABLE_6 = {
    "1": (2.755, (0.5691, -1.320)),
    "2": (2.010, (0.7595, 0.2167, -1.028, 0.4778)),
    "3": (2.919, (0.5813, 0.01216, -0.7065)),
    "4": (2.573, (0.5680, 0.5199)),
    "5": (2.969, (0.5913, 0.004526)),
    "6": (2.898, (0.6239,)),
    "7": (3.138, (0.5090,)),
}
#: Table 8 at the 1-percent AEP: AVP, Sp in percent, MEV.
TABLE_8 = {
    "1": (0.05093, 55.68, 0.04393),
    "2": (0.05468, 57.99, 0.05106),
    "3": (0.03114, 42.37, 0.02745),
    "4": (0.01795, 31.60, 0.01466),
    "5": (0.02382, 36.69, 0.02032),
    "6": (0.06628, 64.89, 0.05626),
    "7": (0.03369, 44.22, 0.02420),
}
#: Rounding of DR-7's coefficients to 4 significant figures, in log10 units of flow.
COEFFICIENT_ROUNDING = 0.0008


@pytest.fixture(scope="module")
def lib():
    return load_state("IL")


def _x(eq, inputs):
    return np.array([1.0, *(v.value(inputs) for v in eq.variables)])


def _xux(eq, inputs):
    x = _x(eq, inputs)
    return float(x @ np.array(eq.covariance) @ x)


def test_illinois_is_verified(lib):
    assert lib.status == "verified"
    assert set(lib.regions) == set(NSS)
    assert len(lib.equations) == 104
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region
        assert set(lib.aeps(region)) == set(AEPS)


@pytest.mark.parametrize("region", sorted(TABLE_6))
def test_table_6_and_table_8_match_the_stored_data_release_values(lib, region):
    """DR-7 and DR-8 are the source; the report prints the 1-percent AEP of each."""
    eq = lib.equation(region, 0.01)
    intercept, coefficients = TABLE_6[region]
    assert eq.intercept == intercept
    assert eq.coefficients == coefficients
    avp, sp, mev = TABLE_8[region]
    assert eq.avp == avp
    assert eq.model_error_variance == mev
    assert eq.n_sites == TABLE_12[region][0]
    assert 100 * math.sqrt(math.exp(math.log(10) ** 2 * avp) - 1) == pytest.approx(sp, abs=0.02)


@pytest.mark.parametrize("region", ["1", "2", "3", "4", "5", "6", "7"])
@pytest.mark.parametrize("aep", AEPS)
def test_every_stored_avp_agrees_with_its_sep_and_covariance_is_symmetric(lib, region, aep):
    eq = lib.equation(region, aep)
    assert eq.sep_log == pytest.approx(math.sqrt(eq.avp), abs=2e-4)
    u = np.array(eq.covariance)
    assert (u == u.T).all()
    assert eq.model_error_variance < eq.avp


def test_table_12_baselines_are_the_medians_of_dr5():
    """Table 12 prints 4 figures; DR-5 prints each gage's median urban fraction."""
    for region, (n, baseline) in TABLE_12.items():
        fractions = [
            g["urban"]
            for g in GAGES.values()
            if g["region"] == region and g["non_redundant"] and g["urban"] is not None
        ]
        assert len(fractions) == n
        assert statistics.median(fractions) == pytest.approx(baseline, abs=1e-5)


def test_table_10_ranges_are_the_extremes_of_the_non_redundant_gages(lib):
    """DR-5 rounds drainage area to 2 decimals, so its small end is checked to that."""
    for region in TABLE_12:
        used = [g for g in GAGES.values() if g["region"] == region and g["non_redundant"]]
        for v in lib.equation(region, 0.01).variables:
            values = [g[v.code] for g in used]
            assert max(values) == v.maximum
            if v.code == "DRNAREA":
                assert min(values) == pytest.approx(v.minimum, abs=0.005)
            else:
                assert min(values) == v.minimum


def test_worked_example_region_2(lib):
    """Figures 14-20: 50 mi2, 20 percent urban, 15 percent poorly drained, 3.67 ft/mi."""
    eq = lib.equation("2", 0.01)
    inputs = {"DRNAREA": 50.0, "FLC16DVLMH": 0.20, "FSSURGDC78": 0.15, "RELRELF": 3.67}
    est = evaluate(eq, inputs)
    assert round(est.flow_cfs) == 1858
    assert _x(eq, inputs) == pytest.approx([1, 1.699, 0.4472, 0.3873, 0.5647], abs=5e-5)
    assert _xux(eq, inputs) == pytest.approx(0.001831, abs=5e-7)
    assert est.prediction_variance == pytest.approx(0.05289, abs=5e-6)
    sp = 100 * math.sqrt(math.exp(math.log(10) ** 2 * est.prediction_variance) - 1)
    assert sp == pytest.approx(56.9, abs=0.05)
    assert stats.t.ppf(0.95, eq.n_sites - 5) == pytest.approx(1.658, abs=5e-4)
    # The report's [773, 4467] rounds S to 0.23; to its 3 figures, [773, 4470].
    assert est.interval[0] == pytest.approx(773, abs=1.5)
    assert est.interval[1] == pytest.approx(4470, abs=5)


def test_figure_16_and_figure_28_are_dr9(lib):
    figure_16 = [
        [34.60, -9.180, -7.487, -6.675, -24.39],
        [-9.180, 3.275, 1.267, -1.495, 6.707],
        [-7.487, 1.267, 9.301, -1.245, 2.409],
        [-6.675, -1.495, -1.245, 37.30, 1.023],
        [-24.39, 6.707, 2.409, 1.023, 21.54],
    ]
    figure_28 = [[19.13, -1.604, -20.03], [-1.604, 0.5936, 0.9224], [-20.03, 0.9224, 24.38]]
    assert np.array(lib.equation("2", 0.01).covariance) * 1000 == pytest.approx(
        np.array(figure_16), abs=1e-9
    )
    assert np.array(lib.equation("4", 0.01).covariance) * 1000 == pytest.approx(
        np.array(figure_28), abs=1e-9
    )


def test_worked_example_region_4_urbanization_adjusted(lib):
    """Figures 24-33: 100 mi2, soil slope 6.0 percent, 20 percent urban."""
    plain = lib.equation("4", 0.01)
    adjusted = lib.equation("4-urban", 0.01)
    inputs = {"DRNAREA": 100.0, "STATSSLPWT": 6.0, "URBTHE2010": 0.20}
    assert round(evaluate(plain, inputs).flow_cfs) == 12989
    est = evaluate(adjusted, inputs)
    assert round(est.flow_cfs) == 14885
    # Figure 26 starts from the rounded 12,989, which moves its 4.17275 by 1.2e-5.
    assert math.log10(12989) + 0.312 * (0.20 - 0.01033) == pytest.approx(4.17275, abs=1e-5)
    assert est.log_flow == pytest.approx(4.17275, abs=2e-5)
    assert _xux(plain, inputs) == pytest.approx(0.001549, abs=5e-7)
    assert evaluate(plain, inputs).prediction_variance == pytest.approx(0.01621, abs=5e-6)
    assert est.prediction_variance == pytest.approx(0.0169, abs=5e-5)
    assert math.sqrt(est.prediction_variance) == pytest.approx(0.130, abs=5e-4)
    sp = 100 * math.sqrt(math.exp(math.log(10) ** 2 * est.prediction_variance) - 1)
    assert sp == pytest.approx(30.6, abs=0.05)
    # To the report's 3 figures, [8,970, 24,700]. The t differs: see the next test.
    assert est.interval[0] == pytest.approx(8970, abs=25)
    assert est.interval[1] == pytest.approx(24700, abs=70)


@pytest.mark.parametrize(
    ("region", "report_df", "report_t"),
    [
        ("4", 33, 1.692),
        ("1", 32, None),
        ("3", 56, None),
        ("5", 33, None),
        ("6", 14, None),
        ("7", 7, None),
    ],
)
def test_adjusted_interval_has_one_degree_of_freedom_fewer_than_the_report(
    lib, region, report_df, report_t
):
    """Pinned: Figure 32 uses n-p, p counting the spatial model only; evaluate counts the urban term."""
    eq = lib.equation(f"{region}-urban", 0.01)
    assert eq.n_sites - len(eq.variables) - 1 == report_df - 1
    if report_t is not None:
        assert stats.t.ppf(0.95, report_df) == pytest.approx(report_t, abs=5e-4)
    if region == "7":
        assert stats.t.ppf(0.95, report_df - 1) == pytest.approx(1.943, abs=5e-4)
        assert stats.t.ppf(0.95, report_df) == pytest.approx(1.895, abs=5e-4)


@pytest.mark.parametrize("region", ["1", "3", "4", "5", "6", "7"])
@pytest.mark.parametrize("aep", AEPS)
def test_adjusted_equation_is_the_unadjusted_one_plus_figures_21_to_23(lib, region, aep):
    plain = lib.equation(region, aep)
    adjusted = lib.equation(f"{region}-urban", aep)
    u0 = TABLE_12[region][1]
    assert adjusted.intercept == plain.intercept
    assert adjusted.coefficients == (*plain.coefficients, BU[aep])
    assert adjusted.variables[:-1] == plain.variables
    urban = adjusted.variables[-1]
    assert (urban.code, urban.transform, urban.offset) == ("URBTHE2010", "identity", -u0)
    assert (adjusted.avp, adjusted.sep_log) == (plain.avp, plain.sep_log)
    assert adjusted.model_error_variance == plain.model_error_variance
    inputs = {v.code: math.sqrt(max(v.minimum, 1e-3) * v.maximum) for v in plain.variables}
    at_baseline = evaluate(adjusted, {**inputs, "URBTHE2010": u0})
    assert at_baseline.log_flow == pytest.approx(evaluate(plain, inputs).log_flow, abs=1e-12)
    assert at_baseline.prediction_variance == pytest.approx(
        evaluate(plain, inputs).prediction_variance, abs=1e-12
    )
    at_half = evaluate(adjusted, {**inputs, "URBTHE2010": 0.5})
    assert at_half.log_flow - evaluate(plain, inputs).log_flow == pytest.approx(
        BU[aep] * (0.5 - u0), abs=1e-12
    )
    assert at_half.prediction_variance - evaluate(plain, inputs).prediction_variance == (
        pytest.approx(BU_SE[aep] ** 2 * (0.5 - u0) ** 2, abs=1e-10)
    )


def _half_unit(flow):
    """Half a unit in the third significant figure, which is how DR-10 prints flows."""
    return 0.5 * 10.0 ** (math.floor(math.log10(flow)) - 2)


def _area_half_unit(area):
    """DR-5 prints drainage area to 2 decimals, or to 4 figures from 100 mi2 up."""
    if area < 100:
        return 0.005
    return 0.05 if area < 1000 else 0.5


def _at_area_ends(eq, gage, get):
    """``get`` of the estimate at both ends of the rounding interval of the printed area."""
    inputs = {v.code: gage[v.code] for v in eq.variables}
    half = _area_half_unit(inputs["DRNAREA"])
    return [
        get(evaluate(eq, {**inputs, "DRNAREA": a}, allow_extrapolation=True))
        for a in (inputs["DRNAREA"] - half, inputs["DRNAREA"] + half)
    ]


def _flow_excess(eq, gage, printed):
    ends = _at_area_ends(eq, gage, lambda e: e.log_flow)
    h = _half_unit(printed)
    return max(min(ends) - math.log10(printed + h), math.log10(printed - h) - max(ends), 0.0)


def _variance_excess(eq, gage, printed):
    ends = _at_area_ends(eq, gage, lambda e: e.prediction_variance)
    return max(min(ends) - printed - 5e-5, printed - 5e-5 - max(ends), 0.0)


def test_there_are_422_gages_with_regression_values():
    assert len(GAGES) == 422
    assert sum(g["non_redundant"] for g in GAGES.values()) == sum(n for n, _ in TABLE_12.values())


def test_regression_flows_reproduce_at_every_gage_outside_regions_1_and_3(lib):
    """DR-10's flows, to their 3 figures, DR-5's rounded areas and DR-7's rounded coefficients."""
    for gage_id, gage in GAGES.items():
        if gage["region"] in ("1", "3"):
            continue
        for i, aep in enumerate(AEPS):
            excess = _flow_excess(lib.equation(gage["region"], aep), gage, gage["flows"][i])
            assert excess <= COEFFICIENT_ROUNDING, (gage_id, aep)


def test_regression_flows_in_regions_1_and_3_show_a_small_tpi_offset(lib):
    """Pinned: the printed TPI forms put these flows above DR-10 by less than 0.5 percent.

    Beyond rounding, the worst value is 0.002 log10 units, and the mean difference is removed
    by a centre of 3.2246 for 3.2259 in region 1 or 1816 for 1823.4 in region 3. So DR-10 and
    DR-11 appear to use slightly different centres than Table 4 prints. Table 4, DR-7 and NSS
    agree, so the printed form is the one stored.
    """
    for gage_id, gage in GAGES.items():
        if gage["region"] in ("1", "3"):
            for i, aep in enumerate(AEPS):
                excess = _flow_excess(lib.equation(gage["region"], aep), gage, gage["flows"][i])
                assert excess <= 0.002, (gage_id, aep)

    def mean_difference(region, transform):
        out = []
        for gage in GAGES.values():
            if gage["region"] != region or gage["DRNAREA"] < 1:
                continue
            for i, aep in enumerate(AEPS):
                eq = lib.equation(region, aep)
                log_q = eq.intercept + sum(
                    b * (transform(gage[v.code]) if v.code == "TPISTATSGO" else v.value(gage))
                    for b, v in zip(eq.coefficients, eq.variables)
                )
                out.append(log_q - math.log10(gage["flows"][i]))
        return float(np.mean(out))

    assert mean_difference("1", lambda t: math.log10(t) - 3.2259) == pytest.approx(
        0.00155, abs=5e-5
    )
    assert abs(mean_difference("1", lambda t: math.log10(t) - 3.2246)) < 2e-4
    assert mean_difference("3", lambda t: (t - 1823.4) / 4758.3) == pytest.approx(0.00087, abs=5e-5)
    assert abs(mean_difference("3", lambda t: (t - 1816.0) / 4758.3)) < 2e-4


def test_covariance_reproduces_the_dr11_variances_at_every_gage(lib):
    """The only check on DR-9, and on its labels: DR-11 prints 4 decimals of MEV + xUx."""
    for gage_id, gage in GAGES.items():
        limit = 7e-5 if gage["region"] == "1" else 3e-5
        for i, aep in enumerate(AEPS):
            excess = _variance_excess(lib.equation(gage["region"], aep), gage, gage["variances"][i])
            assert excess <= limit, (gage_id, aep)


def test_nss_agrees_at_all_104_equations(lib):
    """NSS folds the centres and U0 into the intercept; the flows, ranges and SEP are the same."""
    snapshot = json.loads((DATA / "nss_snapshot" / "IL.json").read_text("utf-8"))
    assert len(snapshot["equations"]) == 104
    for entry in snapshot["equations"]:
        nss = RegressionEquation.from_dict(entry)
        eq = lib.equation(region_for_nss_code(lib, nss.region_code), nss.aep)
        assert [v.code for v in eq.variables] == [v.code for v in nss.variables]
        for ours, theirs in zip(eq.variables, nss.variables):
            assert (ours.minimum, ours.maximum) == (theirs.minimum, theirs.maximum)
        assert eq.sep_log == pytest.approx(nss.sep_log, abs=2e-4)
        for fraction in (0.25, 0.5, 0.75):
            inputs = {
                v.code: (
                    v.minimum ** (1 - fraction) * v.maximum**fraction
                    if v.minimum > 0
                    else v.maximum * fraction
                )
                for v in eq.variables
            }
            assert evaluate(eq, inputs).log_flow == pytest.approx(
                evaluate(nss, inputs).log_flow, abs=1e-11
            )


def test_zero_poorly_drained_fraction_evaluates(lib):
    """Six region 2 gages have 0, Table 10's minimum; NSS's form refuses it."""
    zero = [g for g in GAGES.values() if g["region"] == "2" and g["FSSURGDC78"] == 0.0]
    assert len(zero) == 6
    eq = lib.equation("2", 0.01)
    inputs = {"DRNAREA": 50.0, "FLC16DVLMH": 0.20, "FSSURGDC78": 0.0, "RELRELF": 3.67}
    expected = (
        2.010 + 0.7595 * math.log10(50) + 0.2167 * math.sqrt(0.20) + 0.4778 * math.log10(3.67)
    )
    assert evaluate(eq, inputs).log_flow == pytest.approx(expected, abs=2e-6)
    snapshot = json.loads((DATA / "nss_snapshot" / "IL.json").read_text("utf-8"))
    nss = next(
        RegressionEquation.from_dict(e)
        for e in snapshot["equations"]
        if e["region_code"] == "GC2002" and e["aep"] == 0.01
    )
    with pytest.raises(ValueError, match="not positive"):
        evaluate(nss, inputs)


def test_out_of_range_input_is_refused(lib):
    with pytest.raises(Exception, match="DRNAREA"):
        evaluate(lib.equation("7", 0.01), {"DRNAREA": 1000.0})
