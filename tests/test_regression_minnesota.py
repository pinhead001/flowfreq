"""Minnesota's verified equations: SIR 2023-5079 with its data release, and SIR 2019-5012."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

from flowfreq.regression import evaluate
from flowfreq.regression.equations import OutOfRangeError
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
NSS = {
    "A": "GC1928",
    "B": "GC1929",
    "C": "GC1930",
    "D": "GC1931",
    "E": "GC1932",
    "F": "GC1933",
    "B1": "GC1775",
}
AEPS = [0.667, 0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002]
#: Levin (2023) tables 1 and 5 for regions A-F, its table 3, and SIR 2019-5012 table 1.
GAGES = json.loads((ROOT / "fixtures" / "regression" / "minnesota_gages.json").read_text("utf-8"))


@pytest.fixture(scope="module")
def lib():
    return load_state("MN")


def test_minnesota_is_verified(lib):
    assert lib.status == "verified"
    assert sorted(lib.regions) == sorted(NSS)
    assert len(lib.equations) == 56
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region
        assert lib.aeps(region) == sorted(AEPS, reverse=True)


def test_gage_counts_match_the_equations(lib):
    assert {r: len(g) for r, g in GAGES["release_gages"].items()} == {
        "A": 40,
        "B": 66,
        "C": 20,
        "D": 104,
        "E": 23,
        "F": 45,
    }
    assert len(GAGES["rainy_gages"]) == 49
    for e in lib.equations:
        expected = 49 if e.region_code == "B1" else len(GAGES["release_gages"][e.region_code])
        assert e.n_sites == expected


def test_worked_example_flow_little_cannon_tributary(lib):
    """Example 1, region F, the 1-percent AEP. The report rounds the log to 3.025 first."""
    est = evaluate(lib.equation("F", 0.01), {"DRNAREA": 2.216, "PMPE": 133.265})
    assert est.log_flow == pytest.approx(3.025, abs=5e-4)
    assert est.flow_cfs == pytest.approx(1059.254, rel=6e-4)


def test_worked_example_covariance_is_table_4f_but_its_arithmetic_is_not(lib):
    """Steps 4-7 of Example 1 do not follow from Table 4f, or from each other.

    xU's first two elements reproduce Table 4f, and its third (-0.0000655) does not: it is
    -0.0000595. That gives xUx = 0.00729, not the printed 0.0064798, and Table 5f prints
    0.0306 = MEV + 0.00729 for this very gage. The printed 0.1726 follows from neither
    (0.023 + 0.0064798 gives 0.1717). So the library's interval is not the printed one.
    """
    eq = lib.equation("F", 0.01)
    x = np.array([1.0, math.log10(2.216), 133.265])
    u = np.array(eq.covariance)
    xu = x @ u
    assert xu[0] == pytest.approx(0.0157765, abs=3e-6)
    assert xu[1] == pytest.approx(-0.0016085, abs=1e-7)
    assert xu[2] == pytest.approx(-0.0000595, abs=1e-7)
    xux = float(xu @ x)
    assert xux == pytest.approx(0.00729, abs=1e-5)
    assert eq.model_error_variance == 0.023

    printed = GAGES["release_gages"]["F"]["05355100"]["variance_of_prediction"][AEPS.index(0.01)]
    assert printed == pytest.approx(0.0233 + xux, abs=5e-5)
    assert math.sqrt(0.023 + 0.0064798) == pytest.approx(0.1717, abs=5e-5)

    # The printed chain, from its own 0.1726 and t, is what gives the printed limits.
    assert 10 ** (1.6811 * 0.1726) == pytest.approx(1.95055, abs=5e-4)
    assert 1059.254 / 1.95055 == pytest.approx(543.1, abs=0.05)
    assert 1059.254 * 1.95055 == pytest.approx(2066.1, abs=0.05)
    est = evaluate(eq, {"DRNAREA": 2.216, "PMPE": 133.265})
    assert est.interval == pytest.approx((540.0, 2079.2), abs=0.05)


def test_library_t_follows_equation_2_not_table_4(lib):
    """Table 4's critical t has n-p+1 degrees of freedom, where equation 2 says n-p.

    The release's make_tables.R counts one column too few. There is no field for the printed
    t, and this library takes n-p from n_sites, so it uses the equation.
    """
    printed = {"A": 1.687, "B": 1.670, "C": 1.734, "D": 1.660, "E": 1.721, "F": 1.681}
    assert stats.t.ppf(0.95, 43) == pytest.approx(1.6811, abs=5e-5)
    for region, t in printed.items():
        eq = lib.equation(region, 0.01)
        p = len(eq.variables) + 1
        assert stats.t.ppf(0.95, eq.n_sites - p + 1) == pytest.approx(t, abs=5e-4)
    eq = lib.equation("F", 0.01)
    est = evaluate(eq, {"DRNAREA": 2.216, "PMPE": 133.265})
    factor = 10 ** (stats.t.ppf(0.95, 42) * math.sqrt(est.prediction_variance))
    assert est.interval == pytest.approx((est.flow_cfs / factor, est.flow_cfs * factor))


@pytest.mark.parametrize("region", ["A", "B", "C", "D", "E", "F"])
def test_covariance_gives_the_releases_variance_of_prediction_at_every_gage(lib, region):
    """Table 5 prints MEV + xUx per gage, to 4 decimals, with an unrounded MEV.

    So Table 5 minus xUx must be one constant across the gages (a wrong order or basis
    scatters it), and that constant must round to Table 4's MEV. This is what fixes the
    matrices' order, which nothing else in the report does but for region F's example.
    """
    gages = GAGES["release_gages"][region]
    for i, aep in enumerate(AEPS):
        eq = lib.equation(region, aep)
        u = np.array(eq.covariance)
        implied = []
        for gage in gages.values():
            x = np.array([1.0] + [v.apply(gage[v.code]) for v in eq.variables])
            implied.append(gage["variance_of_prediction"][i] - float(x @ u @ x))
        assert max(implied) - min(implied) < 1.5e-4
        assert np.mean(implied) == pytest.approx(eq.model_error_variance, abs=5.5e-4)


@pytest.mark.parametrize("region", ["A", "B", "C", "D", "E", "F"])
def test_printed_sp_and_avp_are_the_mean_variance_of_prediction(lib, region):
    gages = GAGES["release_gages"][region]
    for i, aep in enumerate(AEPS):
        eq = lib.equation(region, aep)
        mean = float(np.mean([g["variance_of_prediction"][i] for g in gages.values()]))
        assert mean == pytest.approx(eq.avp, abs=5.5e-4)
        assert math.sqrt(mean) == pytest.approx(eq.sep_log, abs=5e-5)


def test_covariances_are_symmetric_and_positive_definite(lib):
    for e in lib.equations:
        if e.region_code == "B1":
            assert e.covariance is None and e.model_error_variance is None
            continue
        u = np.array(e.covariance)
        assert np.allclose(u, u.T, rtol=1e-8, atol=1e-14)
        assert np.linalg.eigvalsh((u + u.T) / 2).min() > 0


def test_release_table_3_reproduces_regions_e_and_d_at_its_gages(lib):
    """Twelve region E gages, 7 AEPs, printed to whole cfs from unrounded coefficients.

    Table 2 rounds each coefficient to 3 decimals, so the limit on each row is what half a
    unit there can move the log, plus half a cfs of printing.
    """
    rows = GAGES["release_table3"]
    assert len(rows) == 84
    gages = GAGES["release_gages"]["E"]
    for row in rows:
        for region in ("E", "D"):
            eq = lib.equation(region, row["aep"])
            inputs = {v.code: gages[row["site"]][v.code] for v in eq.variables}
            est = evaluate(eq, inputs, allow_extrapolation=True)
            terms = sum(abs(v.apply(inputs[v.code])) for v in eq.variables)
            limit = 0.0005 * (1 + terms) + math.log10(1 + 0.5 / row[region])
            assert abs(est.log_flow - math.log10(row[region])) <= limit


#: Two values of SIR 2019-5012 Table 1 are printed 1,000 cfs low.
RAINY_TABLE_1_MISPRINTS = {("05137000", 0.01): 2700.0, ("05132000", 0.5): 4060.0}


def test_rainy_table_1_reproduces_at_every_gage_but_for_two_misprints(lib):
    """Table 1 prints three significant figures."""
    worst = 0.0
    for site, gage in GAGES["rainy_gages"].items():
        for i, aep in enumerate(AEPS):
            eq = lib.equation("B1", aep)
            est = evaluate(
                eq,
                {"DRNAREA": gage["DRNAREA"], "LAKEAREA": gage["LAKEAREA"]},
                allow_extrapolation=True,
            )
            printed = gage["flows"][i]
            if (site, aep) in RAINY_TABLE_1_MISPRINTS:
                assert printed == RAINY_TABLE_1_MISPRINTS[(site, aep)]
                assert est.flow_cfs - printed == pytest.approx(1000, abs=5)
                assert gage["flows"][i - 1] < est.flow_cfs < gage["flows"][i + 1]
                continue
            worst = max(worst, abs(est.flow_cfs / printed - 1))
    assert worst < 0.007


#: Where NSS differs from the data release's tables 2a and 2b: (published, NSS).
NSS_INTERCEPT = {
    ("A", 0.5): (29.61, 29.611),
    ("A", 0.04): (50.33, 50.329),
    ("A", 0.02): (54.035, 54.034),
    ("A", 0.01): (57.39, 57.389),
    ("A", 0.002): (64.245, 64.243),
    ("B", 0.5): (65.098, 65.096),
}
NSS_COEFFICIENT = {
    ("A", 0.04, "DRNAREA"): (0.769, 0.768),
    ("A", 0.02, "LAT_OUT"): (-31.051, -31.05),
    ("A", 0.002, "LAT_OUT"): (-36.994, -36.993),
}
#: SEP differences, log10 units (published, NSS). Table 2c's Sp is 43.9 (NSS 43.96), and Table 5's
#: mean variance of prediction gives 43.90. SIR 2019-5012 Table 3's is 42.2 (NSS 44.22), and its
#: AVP of 0.031 agrees with 42.2.
NSS_SEP = {("C", 0.667): (0.18232, 0.182575), ("B1", 0.1): (0.175817, 0.183537)}


def _digits(code: str) -> int:
    return 4 if code in ("LAT_OUT", "LONG_OUT") else 3


def test_every_equation_equals_the_nss_snapshot_but_for_its_differences(lib):
    snap = json.loads((DATA / "nss_snapshot" / "MN.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        key = (e.region_code, e.aep)
        s = by_key[(NSS[e.region_code], e.aep)]
        if key in NSS_INTERCEPT:
            assert (e.intercept, s["intercept"]) == NSS_INTERCEPT[key]
        else:
            assert e.intercept == s["intercept"]
        assert [v.code for v in e.variables] == [n["code"] for n in s["variables"]]
        for v, b, n, nb in zip(e.variables, e.coefficients, s["variables"], s["coefficients"]):
            if key + (v.code,) in NSS_COEFFICIENT:
                assert (b, nb) == NSS_COEFFICIENT[key + (v.code,)]
            else:
                assert b == nb
            assert v.scale == n.get("scale", 1.0)
            if n["transform"] == "log10_plus1":
                assert (v.transform, v.offset) == ("log10", 1.0)
            else:
                assert v.transform == n["transform"]
            if e.region_code == "B1" or v.code == "PMPE":
                continue
            assert v.minimum == round(n["minimum"], _digits(v.code))
            assert v.maximum == round(n["maximum"], _digits(v.code))
        if key in NSS_SEP:
            assert (e.sep_log, round(s["sep_log"], 6)) == NSS_SEP[key]
        else:
            # SIR 2019-5012 prints Sp to 1 decimal; the release prints 3, as NSS has them.
            limit = 2e-4 if e.region_code == "B1" else 5e-5
            assert e.sep_log == pytest.approx(s["sep_log"], abs=limit)


def test_region_b1_avp_agrees_with_its_sp_as_equation_3_says(lib):
    for aep in AEPS:
        eq = lib.equation("B1", aep)
        assert eq.sep_log**2 == pytest.approx(eq.avp, abs=6e-4)
    nss_sp_44_22 = 0.183537
    assert nss_sp_44_22**2 - lib.equation("B1", 0.1).avp > 2.5e-3


def test_nss_pmpe_limits_are_the_reports_millimeters_in_inches(lib):
    """NSS labels 4.587-8.799 as mm. StreamStats serves mm, so those limits exclude everything."""
    snap = json.loads((DATA / "nss_snapshot" / "MN.json").read_text(encoding="utf-8"))
    eq = lib.equation("F", 0.01)
    pmpe = eq.variables[1]
    assert (pmpe.code, pmpe.units, pmpe.minimum, pmpe.maximum) == ("PMPE", "mm", 116.505, 223.496)
    nss = next(
        v
        for s in snap["equations"]
        if s["region_code"] == "GC1933"
        for v in s["variables"]
        if v["code"] == "PMPE"
    )
    assert (nss["units"], nss["minimum"] * 25.4, nss["maximum"] * 25.4) == pytest.approx(
        ("mm", 116.505, 223.496), abs=1e-3
    )
    evaluate(eq, {"DRNAREA": 2.216, "PMPE": 133.265})
    with pytest.raises(OutOfRangeError):
        evaluate(eq, {"DRNAREA": 2.216, "PMPE": 5.25})


def test_region_b1_keeps_table_3s_limits(lib):
    """NSS serves its gages' extremes, 1,841.547 mi2 and 22.268 percent. Table 3 prints these."""
    eq = lib.equation("B1", 0.01)
    assert [(v.code, v.transform, v.minimum, v.maximum) for v in eq.variables] == [
        ("DRNAREA", "log10", 0.037, 1840.0),
        ("LAKEAREA", "identity", 0.0, 22.3),
    ]
    assert evaluate(eq, {"DRNAREA": 100.0, "LAKEAREA": 0.0}).flow_cfs == pytest.approx(
        10 ** (1.975 + 0.754 * 2)
    )


def test_region_b_takes_a_negative_longitude(lib):
    """StreamStats serves LONG_OUT negative. The regression's log.Long is of its magnitude."""
    eq = lib.equation("B", 0.01)
    long_out = eq.variables[3]
    assert (long_out.scale, long_out.minimum, long_out.maximum) == (-1.0, -96.1198, -91.6585)
    assert long_out.apply(-95.0) == pytest.approx(math.log10(95.0))


def test_region_b1_interval_uses_the_avp(lib):
    est = evaluate(lib.equation("B1", 0.01), {"DRNAREA": 100.0, "LAKEAREA": 5.0})
    assert est.prediction_variance == pytest.approx(0.051)
    assert est.interval_method == "average variance of prediction"
