"""Wisconsin's verified equations: SIR 2022-5118 (rural, with its model archive) and WRI 86-4005."""

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
NSS = {
    "1": "GC1913",
    "2": "GC1914",
    "3": "GC1915",
    "4": "GC1916",
    "Urban": "GC1467",
    "Milwaukee County": "GC1466",
}
AEPS = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)
#: Archive: unrounded characteristics, per-gage VP_R and LP3 flows at the 257 regressed gages.
#: WRI 86-4005: Tables 3 and 4 at its 32 stations; ``milwaukee`` is by station name.
FIXTURE = json.loads((ROOT / "fixtures" / "regression" / "wisconsin_gages.json").read_text("utf-8"))
RURAL, URBAN = FIXTURE["rural"], FIXTURE["urban"]

#: SIR 2022-5118 Table 2 as printed: a, then the coefficients in TABLE2_ORDER, then pseudo R2,
#: Sp, SME and AVP. A second entry, from the PDF, against the archive values the file stores.
TABLE2_ORDER = {
    "1": ("DRNAREA", "SSURGOKSAT", "LC01HERB", "LC01WATER", "I24H10Y"),
    "2": ("DRNAREA", "LC01WATER", "I24H10Y", "FOREST", "WETLAND"),
    "3": ("DRNAREA", "LC01WATER", "WETLAND"),
    "4": ("DRNAREA", "SSURGOKSAT"),
}
TABLE2 = {
    "1": """0.936 0.840 -0.668 -0.374 -0.514 0.423 0.92 53.58 50.83 0.05
0.437 0.834 -0.654 -0.440 -0.583 0.610 0.91 55.98 53.05 0.05
0.205 0.830 -0.638 -0.473 -0.621 0.698 0.90 58.47 55.35 0.06
-0.018 0.825 -0.616 -0.507 -0.663 0.784 0.89 61.00 57.64 0.06
-0.148 0.822 -0.599 -0.529 -0.691 0.835 0.89 62.72 59.19 0.06
-0.252 0.819 -0.582 -0.547 -0.717 0.877 0.88 65.24 61.50 0.07
-0.339 0.816 -0.566 -0.564 -0.741 0.913 0.87 67.77 63.80 0.07
-0.433 0.812 -0.544 -0.583 -0.770 0.953 0.86 71.16 66.90 0.08""",
    "2": """-1.980 0.945 -0.720 1.002 -0.007 -0.007 0.95 40.04 38.08 0.03
-2.172 0.939 -0.770 1.121 -0.007 -0.008 0.94 42.88 40.74 0.03
-2.283 0.936 -0.794 1.183 -0.007 -0.009 0.94 44.83 42.54 0.03
-2.405 0.933 -0.819 1.249 -0.007 -0.009 0.93 47.73 45.21 0.04
-2.490 0.931 -0.834 1.293 -0.007 -0.010 0.92 50.55 47.83 0.04
-2.567 0.929 -0.848 1.332 -0.007 -0.010 0.92 53.35 50.42 0.05
-2.641 0.927 -0.860 1.368 -0.007 -0.010 0.91 56.15 53.00 0.05
-2.732 0.925 -0.874 1.412 -0.007 -0.011 0.90 59.88 56.45 0.06""",
    "3": """1.625 0.785 -0.476 -0.013 0.92 44.60 42.75 0.03
1.882 0.772 -0.536 -0.014 0.91 46.05 44.10 0.04
2.012 0.766 -0.565 -0.015 0.90 48.19 46.10 0.04
2.148 0.762 -0.595 -0.016 0.88 51.60 49.28 0.04
2.234 0.759 -0.614 -0.017 0.87 54.49 51.98 0.05
2.310 0.757 -0.630 -0.017 0.86 57.39 54.69 0.05
2.379 0.756 -0.644 -0.018 0.85 60.33 57.43 0.06
2.461 0.754 -0.661 -0.019 0.83 65.00 61.80 0.07""",
    "4": """2.365 0.623 -0.420 0.86 54.80 52.90 0.05
2.892 0.585 -0.560 0.87 49.92 48.08 0.04
3.168 0.565 -0.637 0.87 49.48 47.55 0.04
3.459 0.546 -0.720 0.86 51.56 49.44 0.04
3.645 0.534 -0.773 0.85 53.41 51.13 0.05
3.810 0.525 -0.822 0.83 56.44 53.96 0.05
3.959 0.517 -0.866 0.82 59.15 56.49 0.06
4.138 0.508 -0.920 0.80 63.83 60.90 0.06""",
}
#: Printed constants, exponents and standard errors of estimate (percent), 2- to 100-year.
URBAN_TABLE5 = (
    (4.18, 0.786, 1.02, 37),
    (9.97, 0.739, 0.910, 32),
    (14.7, 0.723, 0.863, 32),
    (21.5, 0.712, 0.818, 33),
    (27.0, 0.707, 0.792, 35),
    (32.8, 0.704, 0.770, 37),
)
MILWAUKEE_TABLE6 = (
    (3.72, 0.743, 1.11, 39),
    (5.73, 0.727, 1.09, 32),
    (7.05, 0.724, 1.09, 32),
    (8.72, 0.725, 1.08, 32),
    (10.0, 0.727, 1.08, 33),
    (11.3, 0.729, 1.07, 34),
)
#: NSS differences, report values winning: (NSS region, AEP) -> {variable: (NSS, report)}.
NSS_COEFFICIENT_ERRORS = {("GC1916", 0.04): {"DRNAREA": (0.54569, 0.5459)}}


@pytest.fixture(scope="module")
def lib():
    return load_state("WI")


@pytest.fixture(scope="module")
def snapshot():
    snap = json.loads((DATA / "nss_snapshot" / "WI.json").read_text(encoding="utf-8"))
    return {(s["region_code"], s["aep"]): s for s in snap["equations"]}


def _row(eq, gage) -> np.ndarray:
    x = [1.0]
    for v in eq.variables:
        value = v.scale * gage[v.code] + v.offset
        x.append(math.log10(value) if v.transform == "log10" else value)
    return np.array(x)


def test_wisconsin_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == ["1", "2", "3", "4", "Milwaukee County", "Urban"]
    assert len(lib.equations) == 44
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region
    assert {lib.equation("1", a).n_sites for a in AEPS} == {66}
    assert lib.equation("2", 0.01).n_sites == 68
    assert lib.equation("3", 0.01).n_sites == 61
    assert lib.equation("4", 0.01).n_sites == 62
    assert lib.equation("Urban", 0.01).n_sites == 32
    assert lib.equation("Milwaukee County", 0.01).n_sites is None


@pytest.mark.parametrize("region", ["1", "2", "3", "4"])
def test_archive_values_round_to_table_2(lib, region):
    order = TABLE2_ORDER[region]
    lines = TABLE2[region].split("\n")
    assert len(lines) == len(AEPS)
    for aep, line in zip(AEPS, lines):
        values = [float(v) for v in line.split()]
        printed, fit = values[: len(order) + 1], values[len(order) + 1 :]
        eq = lib.equation(region, aep)
        by_code = dict(zip((v.code for v in eq.variables), eq.coefficients))
        assert round(eq.intercept, 3) == printed[0]
        assert [round(by_code[c], 3) for c in order] == printed[1:]
        assert eq.sep_log == pytest.approx(
            math.sqrt(math.log(1 + (fit[1] / 100) ** 2)) / math.log(10), abs=2e-5
        )
        assert round(eq.avp, 2) == fit[3]
        # The report's Sp is the AVP's, and its SME the MEV's, the MEV being printed to 4 places.
        sp = 100 * math.sqrt(math.exp(eq.avp * math.log(10) ** 2) - 1)
        assert sp == pytest.approx(fit[1], abs=0.005 + 1e-9)
        sme = 100 * math.sqrt(math.exp(eq.model_error_variance * math.log(10) ** 2) - 1)
        assert sme == pytest.approx(fit[2], abs=0.04)


def test_rural_equations_round_to_the_nss_snapshot_with_its_errors_pinned(lib, snapshot):
    seen = set()
    for eq in lib.equations:
        if eq.region_code not in TABLE2:
            continue
        key = (NSS[eq.region_code], eq.aep)
        seen.add(key)
        s = snapshot[key]
        assert round(eq.intercept, 4) == s["intercept"]
        nss_by_code = dict(zip((v["code"] for v in s["variables"]), s["coefficients"]))
        for v, b in zip(eq.variables, eq.coefficients):
            pinned = NSS_COEFFICIENT_ERRORS.get(key, {}).get(v.code)
            if pinned:
                assert nss_by_code[v.code] == pinned[0]
                assert round(b, 4) == pinned[1]
            else:
                assert round(b, 4) == nss_by_code[v.code]
        assert eq.sep_log == pytest.approx(s["sep_log"], abs=1e-6)
    assert seen == {k for k in snapshot if k[0] in {NSS[r] for r in TABLE2}}
    assert len(seen) == 32


def test_rural_ranges_with_the_nss_error_pinned(lib, snapshot):
    for eq in lib.equations:
        if eq.region_code not in TABLE2:
            continue
        nss = {v["code"]: v for v in snapshot[(NSS[eq.region_code], eq.aep)]["variables"]}
        for v in eq.variables:
            assert v.minimum == nss[v.code]["minimum"]
            if (eq.region_code, v.code) == ("4", "DRNAREA"):
                assert (v.maximum, nss[v.code]["maximum"]) == (1034.01, 1304.01)
            else:
                assert v.maximum == nss[v.code]["maximum"]
            at_gages = [g[v.code] for g in RURAL.values() if g["region"] == eq.region_code]
            lo, hi = round(min(at_gages), 2), round(max(at_gages), 2)
            if (eq.region_code, v.code) == ("4", "DRNAREA"):
                # Table 3 prints 0.30; the archive rounds 05371300's 0.2947 mi2 to 0.29.
                assert (lo, v.minimum) == (0.29, 0.30)
            else:
                assert lo == v.minimum
            assert hi == v.maximum


def test_covariance_reproduces_every_gages_variance_of_prediction(lib):
    """The archive's table 5 prints VP_R at each gage: MEV + x U x', in the variables' order."""
    assert len(RURAL) == 257
    worst = 0.0
    for region in TABLE2:
        gages = [g for g in RURAL.values() if g["region"] == region]
        for i, aep in enumerate(AEPS):
            eq = lib.equation(region, aep)
            u = np.array(eq.covariance)
            computed = [eq.model_error_variance + _row(eq, g) @ u @ _row(eq, g) for g in gages]
            printed = [g["vp_r"][i] for g in gages]
            worst = max(worst, float(np.max(np.abs(np.subtract(computed, printed)))))
            assert float(np.mean(printed)) == pytest.approx(eq.avp, abs=2e-7)
    # 5e-5 is the rounding of the printed MEV; 1e-6 that of the fixture's VP_R.
    assert worst < 5.1e-5


def test_no_per_gage_flows_are_published_so_the_lp3_residuals_stand_in(lib):
    """Regression minus LP3 is centered at every region and AEP, as a GLS fit must be."""
    for region in TABLE2:
        gages = [g for g in RURAL.values() if g["region"] == region]
        for i, aep in enumerate(AEPS):
            eq = lib.equation(region, aep)
            b = np.array([eq.intercept, *eq.coefficients])
            res = [math.log10(g["lp3"][i]) - _row(eq, g) @ b for g in gages]
            assert abs(float(np.mean(res))) < 0.015
            assert 0.15 < float(np.std(res)) < 0.30


def test_example_1_kewaskum(lib):
    """04086200, region 3, the 1-percent AEP: printed 2.87, 740 ft3/s, and 822 ft3/s weighted."""
    printed_terms = 2.310 + 0.757 * 1.73 - 0.630 * 0.536 - 0.017 * 24.1
    assert printed_terms == pytest.approx(2.87, abs=0.005)
    assert 10**2.87 == pytest.approx(740, rel=0.002)
    # The unrounded WETLAND coefficient is -0.01737, for the printed -0.017: 1 percent lower.
    est = evaluate(lib.equation("3", 0.01), {"DRNAREA": 53.9, "LC01WATER": 2.44, "WETLAND": 24.1})
    assert est.flow_cfs == pytest.approx(733, rel=0.001)
    gage = RURAL["04086200"]
    assert gage["lp3"][5] == pytest.approx(862, abs=0.5)
    assert gage["vp_r"][5] == pytest.approx(0.0568, abs=5e-5)
    # Equations 5 and 6, from the printed values.
    logw = (0.0568 * math.log10(862) + 0.0249 * math.log10(740)) / (0.0568 + 0.0249)
    assert logw == pytest.approx(2.915, abs=5e-4)
    assert 0.0249 * 0.0568 / (0.0249 + 0.0568) == pytest.approx(0.0173, abs=5e-5)


def test_example_2_ungaged_upstream(lib):
    """35.1 mi2, 2.49 percent open water, 20.97 percent wetland: printed 588, then 608 ft3/s."""
    assert 2.310 + 0.757 * math.log10(35.1) - 0.630 * math.log10(3.49) - 0.0174 * 20.97 == (
        pytest.approx(2.77, abs=0.005)
    )
    assert 10**2.77 == pytest.approx(588, rel=0.002)
    est = evaluate(lib.equation("3", 0.01), {"DRNAREA": 35.1, "LC01WATER": 2.49, "WETLAND": 20.97})
    assert est.flow_cfs == pytest.approx(595, rel=0.001)
    assert 588 * (0.697 + (1 - 0.697) * 1.111) == pytest.approx(608, rel=0.001)


def test_example_3_prints_a_log_flow_its_own_terms_do_not_give(lib):
    """23.7 mi2, 0 percent open water, 6.44 percent wetland: printed 3.23, 1,698, 703 to 4,100."""
    own_terms = 2.310 + 0.757 * math.log10(23.7) - 0.0174 * 6.44
    assert own_terms == pytest.approx(3.2386, abs=1e-4)
    assert 10**3.23 == pytest.approx(1698, abs=0.5)
    assert 10**own_terms == pytest.approx(1732, abs=0.5)

    eq = lib.equation("3", 0.01)
    est = evaluate(eq, {"DRNAREA": 23.7, "LC01WATER": 0.0, "WETLAND": 6.44})
    assert est.log_flow == pytest.approx(3.2394, abs=1e-4)
    assert est.flow_cfs == pytest.approx(1735.5, abs=0.5)

    # The printed order: Intercept, WETLAND, LC01WATER, DRNAREA, x = {1, 6.44, 0, 1.37}.
    assert [v.code for v in eq.variables] == ["WETLAND", "LC01WATER", "DRNAREA"]
    x = np.array([1.0, 6.44, 0.0, 1.37])
    xux = float(x @ np.array(eq.covariance) @ x)
    assert xux == pytest.approx(0.00294, abs=2e-5)  # 0.00293 here, 0.00294 printed
    assert eq.model_error_variance == 0.0494
    assert est.prediction_variance**0.5 == pytest.approx(0.2288, abs=5e-5)
    assert 10 ** (1.6736 * est.prediction_variance**0.5) == pytest.approx(2.415, abs=5e-4)
    assert 1698 / 2.415 == pytest.approx(703, abs=0.5)
    assert 1698 * 2.415 == pytest.approx(4100, abs=1)
    # flowfreq's critical value uses n-p = 57 where the printed one used 54.
    t = math.log10(est.flow_cfs / est.interval[0]) / est.prediction_variance**0.5
    assert t == pytest.approx(1.6720, abs=1e-4)


@pytest.mark.parametrize(
    "region,table", [("Urban", URBAN_TABLE5), ("Milwaukee County", MILWAUKEE_TABLE6)]
)
def test_urban_equations_equal_the_snapshot_and_the_printed_tables(lib, snapshot, region, table):
    for aep, (a, b, c, se) in zip(AEPS, table):
        eq = lib.equation(region, aep)
        assert 10**eq.intercept == pytest.approx(a, rel=1e-12)
        assert list(eq.coefficients) == [b, c]
        assert eq.sep_log == pytest.approx(
            math.sqrt(math.log(1 + (se / 100) ** 2)) / math.log(10), abs=1e-6
        )
        assert eq.covariance is None and eq.avp is None
        s = snapshot[(NSS[region], aep)]
        assert eq.intercept == pytest.approx(s["intercept"], abs=1e-9)
        assert s["coefficients"] == [b, c]
        assert "sep_log" not in s
    assert not any(e.aep < 0.01 for e in lib.equations if e.region_code == region)


def test_urban_ranges_are_table_4s_with_the_nss_differences_pinned(lib, snapshot):
    assert len(URBAN) == 32
    for region, stations in (
        ("Urban", list(URBAN.values())),
        ("Milwaukee County", [g for g in URBAN.values() if g["milwaukee"]]),
    ):
        eq = lib.equation(region, 0.01)
        for v in eq.variables:
            values = [g[v.code] for g in stations]
            assert (v.minimum, v.maximum) == (min(values), max(values))
    assert len([g for g in URBAN.values() if g["milwaukee"]]) == 22
    # NSS: 0.47 to 124 mi2 and 4.8 to 90.8 percent for both sets. Table 4's smallest is 0.22.
    nss = snapshot[("GC1467", 0.01)]["variables"]
    assert (nss[0]["minimum"], nss[0]["maximum"]) == (0.47, 124.0)
    assert snapshot[("GC1466", 0.01)]["variables"] == nss
    assert URBAN["05429040"]["DRNAREA"] == 0.22


def test_wri_worked_example_fisher_creek(lib):
    """05430403: 1.88 mi2 and the example's 19.0 percent impervious (Table 4 has 18.7): 494 ft3/s."""
    eq = lib.equation("Urban", 0.01)
    assert evaluate(eq, {"DRNAREA": 1.88, "IMPERV": 19.0}).flow_cfs == pytest.approx(494, abs=0.5)
    assert URBAN["05430403"]["IMPERV"] == 18.7


@pytest.mark.parametrize(
    "region,milwaukee_only,computed,printed",
    [
        ("Urban", False, (39.0, 32.4, 31.8, 33.0, 34.6, 36.6), (37, 32, 32, 33, 35, 37)),
        ("Milwaukee County", True, (42.3, 34.6, 32.7, 32.2, 33.0, 34.3), (39, 32, 32, 32, 33, 34)),
    ],
)
def test_printed_standard_errors_from_tables_3_and_4(
    lib, region, milwaukee_only, computed, printed
):
    """No per-gage regression flows are printed. The residuals at Tables 3 and 4's stations are."""
    stations = [g for g in URBAN.values() if g["milwaukee"] or not milwaukee_only]
    for i, aep in enumerate(AEPS[:6]):
        eq = lib.equation(region, aep)
        res = [
            math.log10(g["q"][i])
            - evaluate(
                eq, {"DRNAREA": g["DRNAREA"], "IMPERV": g["IMPERV"]}, allow_extrapolation=True
            ).log_flow
            for g in stations
        ]
        rms = math.sqrt(sum(r**2 for r in res) / (len(stations) - 3))
        pct = 100 * math.sqrt(math.exp((rms * math.log(10)) ** 2) - 1)
        assert pct == pytest.approx(computed[i], abs=0.05)
        assert abs(float(np.mean(res))) < 0.012
        # Pinned: the 2-year equations, and Milwaukee's 5-year, are 2 to 3.3 points above print.
        assert 0 <= abs(pct - printed[i]) <= (3.31 if i < 2 else 0.75)


def test_milwaukee_interval_has_no_printed_site_count(lib):
    est = evaluate(lib.equation("Milwaukee County", 0.01), {"DRNAREA": 20.0, "IMPERV": 30.0})
    z = math.log10(est.flow_cfs / est.interval[0]) / lib.equation("Milwaukee County", 0.01).sep_log
    assert z == pytest.approx(1.6449, abs=1e-3)
