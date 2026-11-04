"""Alabama's verified equations (SIR 2020-5032 Tables 3-5 and appendix 1; SIR 2010-5012 Tables 2-5)."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.equations import OutOfRangeError
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
NSS = {
    "1": "GC1499",
    "2": "GC1500",
    "3": "GC1501",
    "4": "GC1502",
    "Small streams": "GC1503",
    "Urban": "GC1504",
}
AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
#: See the fixture's ``_source``. Neither report prints a worked example or a covariance.
GAGES = json.loads((ROOT / "fixtures" / "regression" / "alabama_gages.json").read_text("utf-8"))


@pytest.fixture(scope="module")
def lib():
    return load_state("AL")


def _flow(lib, region, aep, **inputs):
    return evaluate(lib.equation(region, aep), inputs, allow_extrapolation=True).flow_cfs


def test_alabama_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == ["1", "2", "3", "4", "Small streams", "Urban"]
    assert len(lib.equations) == 48
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region
        assert set(lib.aeps(region)) == set(AEPS)


def test_site_counts_add_up_to_the_reports(lib):
    """217 regional gages (p. 5), 40 small-stream gages, 23 urban gages."""
    regional = sum(lib.equation(r, 0.01).n_sites for r in ("1", "2", "3", "4"))
    assert regional == 217
    assert lib.equation("Small streams", 0.01).n_sites == 40
    assert lib.equation("Urban", 0.01).n_sites == 23


def test_every_equation_equals_the_nss_snapshot(lib):
    snap = json.loads((DATA / "nss_snapshot" / "AL.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-9)
        assert list(e.coefficients) == s["coefficients"]
        assert e.sep_log == pytest.approx(s["sep_log"], abs=1e-6)
        for v, n in zip(e.variables, s["variables"]):
            assert (v.code, v.minimum, v.maximum) == (n["code"], n["minimum"], n["maximum"])


def test_nss_cites_the_urban_report_to_one_of_its_two_authors(lib):
    snap = json.loads((DATA / "nss_snapshot" / "AL.json").read_text(encoding="utf-8"))
    assert snap["source_reports"][1]["citation"].startswith("Hedgecock, T.S., 2010")
    assert lib.equation("Urban", 0.01).citation.publication.startswith(
        "Hedgecock, T.S., and Lee, K.G., 2010"
    )


@pytest.mark.parametrize("region", ["1", "2", "3"])
def test_appendix_1_flows_reproduce_in_regions_1_to_3(lib, region):
    """Three significant figures on the flow, and a drainage area of three more or less."""
    gages = {k: g for k, g in GAGES["rural"].items() if g["region"] == region}
    assert len(gages) == {"1": 68, "2": 20, "3": 18}[region]
    for gage in gages.values():
        for i, aep in enumerate(AEPS):
            flow = _flow(lib, region, aep, CONTDA=gage["CONTDA"])
            assert flow == pytest.approx(gage["flows"][i], rel=0.0056)


#: Region 4 of appendix 1 at most gages, against Table 3: a constant ratio at three AEPs.
#: Table 3 and NSS agree, so Table 3 stands.
REGION_4_APPENDIX_RATIO = {0.5: 0.992, 0.2: 0.989, 0.02: 1.004}
#: The gages whose appendix 1 flows do follow Table 3 at those three AEPs.
REGION_4_FOLLOWING_TABLE_3 = {
    "02427300",
    "02427875",
    "02360000",
    "02377500",
    "02342933",
    "02342150",
    "02342500",
    "02377570",
    "02419000",
    "02373000",
}


def test_region_4_appendix_1_flows_reproduce_at_the_other_aeps(lib):
    gages = {k: g for k, g in GAGES["rural"].items() if g["region"] == "4"}
    assert len(gages) == 46
    for gage in gages.values():
        for i, aep in enumerate(AEPS):
            if aep not in REGION_4_APPENDIX_RATIO:
                flow = _flow(lib, "4", aep, CONTDA=gage["CONTDA"])
                assert flow == pytest.approx(gage["flows"][i], rel=0.0056)


def test_region_4_appendix_1_splits_into_table_3_and_a_near_copy(lib):
    gages = {k: g for k, g in GAGES["rural"].items() if g["region"] == "4"}
    for aep, ratio in REGION_4_APPENDIX_RATIO.items():
        i = AEPS.index(aep)
        for name in REGION_4_FOLLOWING_TABLE_3:
            flow = _flow(lib, "4", aep, CONTDA=gages[name]["CONTDA"])
            assert gages[name]["flows"][i] / flow == pytest.approx(1.0, abs=0.0035)
        others = [
            gages[k]["flows"][i] / _flow(lib, "4", aep, CONTDA=gages[k]["CONTDA"])
            for k in gages
            if k not in REGION_4_FOLLOWING_TABLE_3
        ]
        assert sum(others) / len(others) == pytest.approx(ratio, abs=0.0015)


def test_appendix_1_small_stream_rows_use_region_3s_equations_not_table_4s(lib):
    """Every one of the 15 rows; data release table 3 prints the same flows."""
    assert len(GAGES["small_stream"]) == 15
    for gage in GAGES["small_stream"].values():
        for i, aep in enumerate(AEPS):
            printed = gage["flows"][i]
            assert _flow(lib, "3", aep, CONTDA=gage["CONTDA"]) == pytest.approx(printed, rel=0.0056)
            assert printed / _flow(lib, "Small streams", aep, CONTDA=gage["CONTDA"]) > 1.25


def _small_stream_bias(lib, region, aep, k):
    """Mean of log10(mean + k * sd of the 40 small-stream gages) less the equation's log flow."""
    gages = [g for g in GAGES["moments"].values() if "Small" in g["region"]]
    assert len(gages) == 40
    residuals = [
        g["mean"] + k * g["sd"] - math.log10(_flow(lib, region, aep, CONTDA=g["CONTDA"]))
        for g in gages
    ]
    return sum(residuals) / len(residuals)


@pytest.mark.parametrize("aep, k", [(0.5, 0.0), (0.01, 2.326)])
def test_table_4_fits_the_40_small_stream_gages_and_region_3_does_not(lib, aep, k):
    """The check on Table 4 that appendix 1 cannot give: gage moments, at zero skew."""
    assert abs(_small_stream_bias(lib, "Small streams", aep, k)) < 0.015
    assert _small_stream_bias(lib, "3", aep, k) < -0.2


@pytest.mark.parametrize("region", ["1", "2", "3", "4"])
def test_small_seps_are_as_printed_and_table_3s_mse_is_not_model_error(lib, region):
    """At 10 percent, MSE = AVP + mean at-site sampling variance = mean squared gage residual.

    The sampling variance is the large-sample s**2 * (1 + k**2 / 2) / n, at zero skew.
    """
    k = 1.2816
    eq = lib.equation(region, 0.1)
    mse = GAGES["table3_mse"][region][AEPS.index(0.1)]
    # Data release table 1 labels 02448500 region 2; appendix 1 and Table 5, region 3.
    gages = [
        g
        for name, g in GAGES["moments"].items()
        if (g["region"].split("/")[0] == region) != (name == "02448500" and region in ("2", "3"))
    ]
    assert len(gages) == eq.n_sites
    sampling = sum(g["sd"] ** 2 * (1 + k**2 / 2) / g["n"] for g in gages) / len(gages)
    assert eq.avp + sampling == pytest.approx(mse, abs=0.0015)
    squared = [
        (g["mean"] + k * g["sd"] - math.log10(_flow(lib, region, 0.1, CONTDA=g["CONTDA"]))) ** 2
        for g in gages
    ]
    assert sum(squared) / len(squared) == pytest.approx(mse, abs=0.0015)
    assert eq.model_error_variance is None


def test_table_3_mse_exceeds_avp_almost_everywhere(lib):
    rows = [
        (mse, lib.equation(region, aep).avp)
        for region, values in GAGES["table3_mse"].items()
        for aep, mse in zip(AEPS, values)
    ]
    assert len(rows) == 40
    assert sum(mse > avp for mse, avp in rows) == 39
    assert GAGES["table3_mse"]["2"][0] == lib.equation("2", 0.5).avp == 0.044


def test_the_smallest_sep_is_the_nine_percent_the_report_prints(lib):
    """Text, p. 5: SEP from 9 to 58 percent. The floor in the snapshot's plausibility test."""
    seps = {(e.region_code, e.aep): e.sep_log for e in lib.equations if e.region_code != "Urban"}
    assert min(seps, key=seps.get) == ("4", 0.1)
    assert seps["4", 0.1] == pytest.approx(0.039008, abs=1e-6)
    assert max(seps, key=seps.get) == ("Small streams", 0.5)
    assert seps["Small streams", 0.5] == pytest.approx(0.233866, abs=1e-6)


def test_avp_agrees_with_sep_to_their_rounding(lib):
    """Variance from SEP +/- 0.5 percent overlaps AVP +/- half its last printed digit."""

    def variance(percent):
        return math.log(1 + (percent / 100) ** 2) / math.log(10) ** 2

    for e in lib.equations:
        percent = round(100 * math.sqrt(math.exp((e.sep_log * math.log(10)) ** 2) - 1))
        half = 0.00005 if e.region_code == "Urban" else 0.0005
        assert variance(percent + 0.5) >= e.avp - half, (e.region_code, e.aep)
        assert variance(percent - 0.5) <= e.avp + half, (e.region_code, e.aep)


def _weighted(lib, gage, i, regression_scale):
    """Equation 2 of SIR 2010-5012, weighting Table 2's LP3 flow with the regression's."""
    aep = AEPS[i]
    eq = lib.equation("Urban", aep)
    q_r = _flow(lib, "Urban", aep, CONTDA=gage["CONTDA"], URBAN=gage["URBAN"])
    v_r = eq.avp * regression_scale
    v_s = gage["vp_gage"][i]
    return 10 ** ((v_r * math.log10(gage["lp3"][i]) + v_s * math.log10(q_r)) / (v_s + v_r))


#: The gages whose "weighted-average or best-estimate" flows do follow equation 2.
URBAN_WEIGHTED = {
    "02423400",
    "02471065",
    "02480002",
    "0357526200",
    "0357568980",
    "03575700",
    "0357587090",
    "0357587400",
    "03575890",
    "0357591500",
    "03575980",
}
#: What the report's regression variance is in its weighting, not Table 4's.
REGRESSION_VARIANCE_SCALE = 0.4343


def test_urban_weighted_flows_reproduce_at_11_gages(lib):
    assert len(GAGES["urban"]) == 23
    for name in URBAN_WEIGHTED:
        gage = GAGES["urban"][name]
        for i in range(8):
            flow = _weighted(lib, gage, i, REGRESSION_VARIANCE_SCALE)
            assert flow == pytest.approx(gage["weighted"][i], rel=0.007), (name, AEPS[i])


def test_urban_best_estimates_at_the_other_12_gages_do_not_follow_equation_2(lib):
    others = set(GAGES["urban"]) - URBAN_WEIGHTED
    assert len(others) == 12
    for name in others:
        gage = GAGES["urban"][name]
        worst = max(
            abs(_weighted(lib, gage, i, REGRESSION_VARIANCE_SCALE) / gage["weighted"][i] - 1)
            for i in range(8)
        )
        assert worst > 0.01, name


def test_table_4_avp_is_not_the_regression_variance_table_5_weighted_with(lib):
    """Table 4's matches its own SEP, but Table 5's weighted variances need 0.4343 of it."""
    within = [0, 0]
    total = 0
    for gage in GAGES["urban"].values():
        for i, aep in enumerate(AEPS):
            total += 1
            v_s, v_w = gage["vp_gage"][i], gage["vp_weighted"][i]
            for j, scale in enumerate((1.0, REGRESSION_VARIANCE_SCALE)):
                v_r = lib.equation("Urban", aep).avp * scale
                within[j] += abs(v_s * v_r / (v_s + v_r) - v_w) <= 2.1e-4
    assert total == 184
    assert within[0] <= 6
    assert within[1] >= 178


def test_table_5_gage_variance_out_of_sequence():
    """02473047 at 4 percent prints 0.0046 between 0.0013 and 0.0035; its weighted is 0.0015."""
    gage = GAGES["urban"]["02473047"]
    assert gage["vp_gage"][2:5] == [0.0013, 0.0046, 0.0035]
    assert gage["vp_weighted"][3] == 0.0015


def test_urban_gages_below_the_stated_limit_of_development(lib):
    """The report says the equations do not apply below 20 percent; two of its gages are."""
    assert lib.equation("Urban", 0.01).variables[1].minimum == 20.0
    below = sorted(k for k, g in GAGES["urban"].items() if g["URBAN"] < 20.0)
    assert below == ["02423130", "0357526200"]


def test_table_5_and_appendix_1_put_02448500_in_region_3_against_the_data_release():
    assert GAGES["rural"]["02448500"]["region"] == "3"
    assert GAGES["moments"]["02448500"]["region"] == "2"
    assert GAGES["moments"]["02448500"]["CONTDA"] == 1100.0


def test_region_3_lower_limit_is_table_5s_though_the_release_has_no_such_gage(lib):
    assert lib.equation("3", 0.01).variables[0].minimum == 0.34
    in_3 = [g["CONTDA"] for g in GAGES["moments"].values() if g["region"].split("/")[0] == "3"]
    assert min(in_3) == 0.97


def test_small_streams_and_rural_regions_stay_separate(lib):
    """Both apply from 0.13 to 14 mi2. The report's preference below 5 mi2 is in the notes."""
    assert lib.equation("Small streams", 0.01).variables[0].maximum == 14.0
    flows = (
        evaluate(lib.equation("Small streams", 0.01), {"CONTDA": 2.0}).flow_cfs,
        evaluate(lib.equation("4", 0.01), {"CONTDA": 2.0}).flow_cfs,
    )
    assert flows == pytest.approx((962 * 2.0**0.696, 1380 * 2.0**0.561))


def test_urban_evaluates_as_printed(lib):
    est = evaluate(lib.equation("Urban", 0.01), {"CONTDA": 10.0, "URBAN": 50.0})
    assert est.flow_cfs == pytest.approx(618 * 10.0**0.656 * 50.0**0.223)


def test_out_of_range_urban_is_refused(lib):
    with pytest.raises(OutOfRangeError):
        evaluate(lib.equation("Urban", 0.01), {"CONTDA": 10.0, "URBAN": 10.0})
