"""Nevada peak-flow regression equations (issue #104).

Source: Thomas, Hjalmarson and Waltemeyer, 1997, USGS WSP 2433, Tables 5, 6, 7, 9 (GLS) and
10, 14 (hybrid method). Every literal number below is copied from the report, or from a
dated live NSS response, never from ``NV.json`` itself.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from flowfreq.regression import EquationsUnavailable, OutOfRangeError, evaluate, load_state

REPO_ROOT = Path(__file__).resolve().parents[1]

AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01]
BASIN = {"DRNAREA": 50.0, "PRECIP": 20.0, "ELEV": 6500.0, "LAT_GAGE": 38.0}


@pytest.fixture(scope="module")
def lib():
    return load_state("NV")


def _sig3(x: float) -> float:
    return float(f"{x:.3g}")


class TestLibraryShape:
    def test_status_and_regions(self, lib):
        assert lib.status == "verified"
        assert lib.issue == 104
        assert lib.regions == ["1", "10", "2", "3", "5", "6"]
        for region in ("1", "2", "3", "5", "6", "10"):
            assert lib.aeps(region) == AEPS
        assert len(lib.equations) == 36

    def test_region_6_two_year_is_the_printed_constant(self, lib):
        """Table 10 (p. 46) prints the 2-year equation as Q = 0."""
        eq = lib.equation("6", 0.5)
        assert (eq.form, eq.constant_cfs, eq.intercept) == ("constant", 0.0, None)
        est = evaluate(eq, BASIN)
        assert est.flow_cfs == 0.0 and est.log_flow == -math.inf
        assert est.interval is None and "printed constant" in est.interval_method
        assert eq.published_errors == ()  # Table 10 prints '--'

    def test_region_6_two_year_is_range_checked(self, lib):
        """The constant carries region 6's ranges, the same as its 5-year equation."""
        eq = lib.equation("6", 0.5)
        assert eq.variables == lib.equation("6", 0.2).variables
        assert set(eq.coefficients) == {0.0}
        with pytest.raises(OutOfRangeError):
            evaluate(eq, {"DRNAREA": 5000.0, "ELEV": 6500.0})
        with pytest.raises(KeyError):
            evaluate(eq, {"DRNAREA": 50.0})

    def test_only_constant_is_region_6_two_year(self, lib):
        from flowfreq.regression.library import DATA_DIR

        constants = [
            (path.stem, e.region_code, e.aep)
            for path in sorted(DATA_DIR.glob("*.json"))
            for e in load_state(path.stem).equations
            if e.form == "constant"
        ]
        assert constants == [("NV", "6", 0.5)]

    def test_constant_entry_has_no_intercept_key(self):
        """A reader that predates ``form`` does ``float(d["intercept"])``: it must fail
        loudly on the constant cell, not misread it."""
        import json

        from flowfreq.regression.library import DATA_DIR

        data = json.loads((DATA_DIR / "NV.json").read_text(encoding="utf-8"))
        (cell,) = [e for e in data["equations"] if e.get("form") == "constant"]
        assert "intercept" not in cell and "coefficients" not in cell
        with pytest.raises(KeyError, match="intercept"):
            float(cell["intercept"])

    @pytest.mark.parametrize(
        "region, n", [("1", 165), ("2", 108), ("3", 35), ("5", 37), ("6", 80), ("10", 104)]
    )
    def test_n_sites(self, lib, region, n):
        assert {e.n_sites for e in lib.equations if e.region_code == region} == {n}

    def test_hybrid_regions_carry_no_prediction_error(self, lib):
        for region in ("6", "10"):
            for eq in (e for e in lib.equations if e.region_code == region):
                assert eq.sep_log is None and eq.model_error_variance is None
                assert eq.avp is None and eq.covariance is None
                est = evaluate(eq, BASIN)
                assert est.interval is None and est.interval_method.startswith("none")

    @pytest.mark.parametrize(
        "region, source, printed",
        [
            # Table 10, p. 46 (2-year: '--'); Table 14, p. 54. "Estimated average standard
            # error of regression, in log units", copied from the rendered pages.
            ("6", "Table 10, p. 46", {5: 1.47, 10: 1.12, 25: 0.796, 50: 1.10, 100: 1.84}),
            (
                "10",
                "Table 14, p. 54",
                {2: 1.14, 5: 0.602, 10: 0.675, 25: 0.949, 50: 0.928, 100: 1.23},
            ),
        ],
    )
    def test_hybrid_published_errors(self, lib, region, source, printed):
        stored = {}
        for eq in (e for e in lib.equations if e.region_code == region and e.form != "constant"):
            (pe,) = eq.published_errors
            assert (pe.kind, pe.units, pe.source) == (
                "hybrid_standard_error_of_regression",
                "log units",
                source,
            )
            stored[round(1 / eq.aep)] = pe.value
        assert stored == printed

    def test_hybrid_published_errors_match_nss_snapshot(self, lib):
        """NSS carries 6 of the 11 hybrid errors, under the GLS regions' 'SE' code."""
        import json

        snap = REPO_ROOT / "data" / "nss_snapshots" / "NV_2026-10-01.json"
        if not snap.exists():
            pytest.skip("NSS snapshot not present")
        cells = []

        def walk(o):
            if isinstance(o, dict):
                if o.get("region_code") in ("GC164", "GC165") and "equations" in o:
                    cells.extend((o["region_code"], e) for e in o["equations"])
                for v in o.values():
                    walk(v)
            elif isinstance(o, list):
                for v in o:
                    walk(v)

        walk(json.loads(snap.read_text(encoding="utf-8")))
        compared = set()
        for code, e in cells:
            se = [x["value"] for x in e["errors"] if x["code"] == "SE"]
            if not se:
                continue
            eq = lib.equation(NSS_CODE_TO_REGION[code], NSS_STAT_TO_AEP[e["statistic_code"]])
            assert eq.published_errors[0].value == se[0], (code, e["statistic_code"])
            compared.add((code, e["statistic_code"]))
        assert len(compared) == 6

    def test_no_other_state_has_published_errors(self):
        from flowfreq.regression.library import DATA_DIR

        for path in sorted(DATA_DIR.glob("*.json")):
            if path.stem != "NV":
                assert not any(e.published_errors for e in load_state(path.stem).equations)

    def test_latitude_term(self, lib):
        v = {v.code: v for v in lib.equation("5", 0.01).variables}["LAT_GAGE"]
        assert (v.transform, v.scale, v.offset) == ("log10", 0.1, -2.8)


class TestPrintedForms:
    @pytest.mark.parametrize(
        "region, aep, printed",
        [
            ("1", 0.5, lambda b: 0.124 * b["DRNAREA"] ** 0.845 * b["PRECIP"] ** 1.44),
            ("2", 0.01, lambda b: 148 * b["DRNAREA"] ** 0.752 * (b["ELEV"] / 1000) ** -0.584),
            ("3", 0.1, lambda b: 1.99 * b["DRNAREA"] ** 0.633 * b["PRECIP"] ** 0.924),
            (
                "5",
                0.01,
                lambda b: 7000
                * b["DRNAREA"] ** 0.782
                * (b["ELEV"] / 1000) ** -2.18
                * ((b["LAT_GAGE"] - 28) / 10) ** 4.6,
            ),
            ("5", 0.1, lambda b: 28.0 * b["DRNAREA"] ** 0.826 * ((b["LAT_GAGE"] - 28) / 10) ** 4.3),
            ("6", 0.2, lambda b: 32 * b["DRNAREA"] ** 0.80 * (b["ELEV"] / 1000) ** -0.66),
            ("10", 0.5, lambda b: 12 * b["DRNAREA"] ** 0.58),
        ],
    )
    def test_power_form(self, lib, region, aep, printed):
        assert evaluate(lib.equation(region, aep), BASIN).flow_cfs == pytest.approx(
            printed(BASIN), rel=1e-12
        )

    @pytest.mark.parametrize("region, aep, sep", [("1", 0.5, 59), ("2", 0.01, 68), ("5", 0.5, 135)])
    def test_sep(self, lib, region, aep, sep):
        eq = lib.equation(region, aep)
        assert 100 * math.sqrt(10 ** (math.log(10) * eq.sep_log**2) - 1) == pytest.approx(
            sep, abs=1e-3
        )


class TestWorkedExample:
    """pp. 39-40: A 57 mi2, ELEV 6,500 ft, in regions 6 and 10."""

    @pytest.mark.parametrize(
        "region, aep, printed",
        [("6", 0.1, 362), ("6", 0.01, 2120), ("10", 0.1, 2450), ("10", 0.01, 13800)],
    )
    def test_regions_6_and_10(self, lib, region, aep, printed):
        q = evaluate(lib.equation(region, aep), {"DRNAREA": 57.0, "ELEV": 6500.0}).flow_cfs
        assert _sig3(q) == printed


# Live NSS responses for BASIN, fetched 2026-10-03 (3 significant figures).
NSS_RECORDED = [
    ("1", 0.5, 253),
    ("1", 0.01, 943),
    ("2", 0.01, 940),
    ("3", 0.01, 704),
    ("5", 0.5, 141),
    ("5", 0.01, 2520),
    ("6", 0.01, 1990),
    ("10", 0.01, 12600),
]


@pytest.mark.parametrize("region, aep, nss", NSS_RECORDED)
def test_worked_example_matches_recorded_nss(lib, region, aep, nss):
    assert _sig3(evaluate(lib.equation(region, aep), BASIN).flow_cfs) == nss


NSS_CODE_TO_REGION = {
    "GC160": "1",
    "GC161": "2",
    "GC162": "3",
    "GC163": "5",
    "GC164": "6",
    "GC165": "10",
}
NSS_STAT_TO_AEP = {
    "PK50AEP": 0.5,
    "PK20AEP": 0.2,
    "PK10AEP": 0.1,
    "PK4AEP": 0.04,
    "PK2AEP": 0.02,
    "PK1AEP": 0.01,
}


@pytest.mark.requires_network
@pytest.mark.parametrize(
    "basin",
    [
        BASIN,
        {"DRNAREA": 5.0, "PRECIP": 12.0, "ELEV": 5800.0, "LAT_GAGE": 37.0},
        {"DRNAREA": 200.0, "PRECIP": 40.0, "ELEV": 7900.0, "LAT_GAGE": 39.4},
    ],
)
def test_matches_live_nss(lib, basin):
    """Every stored region x AEP against NSS region 32 (Nevada), citation 17."""
    from flowfreq.streamstats import Characteristic, estimate_flow_statistics

    chars = {k: Characteristic(k, k, "", v, "") for k, v in basin.items()}
    results, skipped = estimate_flow_statistics("NV", chars, ["PFS"])
    assert not skipped
    compared = 0
    for r in results:
        region = NSS_CODE_TO_REGION[r.region_code]
        for code, est in r.estimates.items():
            if region == "6" and code == "PK50AEP":
                assert est.value == 0  # Table 10: Q = 0
            ours = evaluate(lib.equation(region, NSS_STAT_TO_AEP[code]), basin).flow_cfs
            assert _sig3(ours) == est.value, (region, code, est.value, ours)
            compared += 1
    assert compared == 36


# Data section (pp. 109-195), the gages used in each region: minimum and maximum of each
# explanatory variable, with the station that sets it, read from the rendered rows.
DATA_SECTION_RANGES = [
    ("1", "DRNAREA", 0.30, 1061.0),  # 10205100 Sheep Creek; 09119000 Tomichi Creek
    ("1", "PRECIP", 10.0, 45.0),  # 09204700; 09352500
    ("2", "DRNAREA", 0.16, 1680.0),  # 13172668; 14040500 John Day River at Picture Gorge
    ("2", "ELEV", 2700.0, 9090.0),  # 13228300 Lytle Creek; 10316500 Lamoille Creek
    ("3", "DRNAREA", 2.22, 1450.0),  # 13145700; 13105000 Salmon Falls Creek
    ("3", "PRECIP", 10.0, 46.0),  # 13075600; 13047500 Falls River
    ("5", "DRNAREA", 0.64, 886.0),  # 10336635; 10311000 Carson River near Carson City
    ("5", "ELEV", 5770.0, 10500.0),  # 10311450; 10265700
    ("5", "LAT_GAGE", 36.439, 39.501),  # 10286000; 10342000
    ("6", "DRNAREA", 0.20, 1670.0),  # 10249140; 09418500 Meadow Valley Wash
    ("6", "ELEV", 3340.0, 9960.0),  # 09415800; 10249900 Chiatovich Creek
    ("10", "DRNAREA", 0.01, 1693.0),  # 09429250; 10255885 San Felipe Creek
]


@pytest.mark.parametrize("region, code, lo, hi", DATA_SECTION_RANGES)
def test_limits_from_data_section(lib, region, code, lo, hi):
    for eq in (e for e in lib.equations if e.region_code == region):
        for v in eq.variables:
            if v.code == code:
                assert (v.minimum, v.maximum) == (lo, hi)
                assert "data section" in eq.citation.table


def test_ranges_contain_nss_limits(lib):
    """NSS's limits (snapshot 2026-10-01) lie inside the data-section ranges, to NSS's rounding."""
    nss = {
        ("1", "DRNAREA"): (0.6, 1060.0),
        ("1", "PRECIP"): (11.0, 43.0),
        ("2", "DRNAREA"): (0.8, 1680.0),
        ("2", "ELEV"): (3540.0, 7950.0),
        ("3", "DRNAREA"): (2.2, 1450.0),
        ("3", "PRECIP"): (10.0, 41.0),
        ("5", "DRNAREA"): (4.1, 360.0),
        ("5", "ELEV"): (5770.0, 10500.0),
        ("5", "LAT_GAGE"): (36.44, 39.5),
        ("6", "DRNAREA"): (0.2, 210.0),
        ("6", "ELEV"): (4770.0, 9960.0),
        ("10", "DRNAREA"): (0.1, 1000.0),
    }
    for (region, code), (lo, hi) in nss.items():
        v = {v.code: v for v in lib.equation(region, 0.01).variables}[code]
        assert v.minimum <= lo * 1.01 and hi <= v.maximum * 1.01, (region, code)


class TestTransitionZones:
    """WSP 2433 p. 19, equations 6 and 7; worked examples pp. 39-44."""

    def test_boundary_south_of_41(self):
        from flowfreq.regression.nevada import high_elevation_boundary_ft

        assert high_elevation_boundary_ft(37.0) == 7500.0
        assert high_elevation_boundary_ft(41.0) == 7500.0

    def test_boundary_north_of_41_is_graphical_only(self):
        from flowfreq.regression.nevada import high_elevation_boundary_ft

        with pytest.raises(ValueError, match="figure 5"):
            high_elevation_boundary_ft(41.5)

    @pytest.mark.parametrize(
        "elev, low", [(6700.0, 1.0), (6800.0, 1.0), (7150.0, 0.5), (7500.0, 0.0), (8000.0, 0.0)]
    )
    def test_elevation_weights(self, elev, low):
        from flowfreq.regression.nevada import elevation_transition_weights

        w = elevation_transition_weights(elev, 7500.0)
        assert w["low"] == pytest.approx(low) and w["high"] == pytest.approx(1 - low)

    def test_equation_6_worked_example(self, lib):
        """p. 41: 21 mi2 in Region 6 and 36 mi2 in Region 10, ELEV 6,500 ft:
        Q10(w) = (362x21 + 2,450x36)/57 = 1,680; Q100(w) = (2,120x21 + 13,800x36)/57 = 9,500."""
        from flowfreq.regression.nevada import estimate_area_weighted

        chars = {"DRNAREA": 57.0, "ELEV": 6500.0}
        q10 = estimate_area_weighted(0.1, {"6": 21.0, "10": 36.0}, chars, lib=lib)
        q100 = estimate_area_weighted(0.01, {"6": 21.0, "10": 36.0}, chars, lib=lib)
        assert q10.weights == pytest.approx({"6": 21 / 57, "10": 36 / 57})
        assert q10.blended and q10.method == "area (eq. 6)"
        # The report weights its rounded component flows; the unrounded ones agree to rounding.
        assert (362 * 21 + 2450 * 36) / 57 == pytest.approx(1680, abs=5)
        assert (2120 * 21 + 13800 * 36) / 57 == pytest.approx(9500, abs=50)
        assert round(q10.flow_cfs, -1) == 1680
        assert round(q100.flow_cfs, -2) == 9500

    def test_equation_7_worked_example(self, lib):
        """pp. 42-44: Region 8 site at 7,100 ft, B = 7,500 ft; A 45 mi2, PREC 28 in.
        Region 1 gives Q2 375 and Q50 975; with Region 8's 433 and 2,440 the weighted
        values are 408 and 1,810 ft3/s. Region 8 is not a Nevada region, so its printed
        flows are combined here with the stored Region 1 equations."""
        from flowfreq.regression.nevada import (
            elevation_transition_weights,
            high_elevation_boundary_ft,
        )

        w = elevation_transition_weights(7100.0, high_elevation_boundary_ft(37.0))
        assert w["low"] == pytest.approx(400 / 700)
        r1 = {"DRNAREA": 45.0, "PRECIP": 28.0}
        q2_h = evaluate(lib.equation("1", 0.5), r1).flow_cfs
        q50_h = evaluate(lib.equation("1", 0.02), r1).flow_cfs
        assert (_sig3(q2_h), _sig3(q50_h)) == (375, 975)
        assert _sig3(433 * w["low"] + q2_h * w["high"]) == 408
        assert _sig3(2440 * w["low"] + q50_h * w["high"]) == 1810

    def test_elevation_transition_blends_nevada_regions(self, lib):
        from flowfreq.regression.nevada import estimate_elevation_transition

        chars = {"DRNAREA": 50.0, "PRECIP": 20.0, "ELEV": 6500.0}
        mid = estimate_elevation_transition(
            0.01, "2", chars, site_elevation=7150.0, latitude=38.0, lib=lib
        )
        q2 = evaluate(lib.equation("2", 0.01), chars).flow_cfs
        q1 = evaluate(lib.equation("1", 0.01), chars).flow_cfs
        assert mid.flow_cfs == pytest.approx(0.5 * q2 + 0.5 * q1, rel=1e-12)
        assert mid.weights == {"2": 0.5, "1": 0.5} and mid.method == "site elevation (eq. 7)"
        below = estimate_elevation_transition(
            0.01, "2", chars, site_elevation=6000.0, latitude=38.0, lib=lib
        )
        assert not below.blended and below.flow_cfs == pytest.approx(q2, rel=1e-12)
        above = estimate_elevation_transition(
            0.01, "2", chars, site_elevation=7600.0, latitude=38.0, lib=lib
        )
        assert list(above.estimates) == ["1"] and above.flow_cfs == pytest.approx(q1, rel=1e-12)
        north = estimate_elevation_transition(
            0.01, "2", chars, site_elevation=6600.0, boundary_ft=7000.0, lib=lib
        )
        assert north.weights["2"] == pytest.approx(400 / 700)

    def test_errors(self, lib):
        from flowfreq.regression.nevada import estimate_area_weighted, estimate_elevation_transition

        chars = {"DRNAREA": 50.0, "PRECIP": 20.0, "ELEV": 6500.0}
        with pytest.raises(ValueError, match="figure 5"):
            estimate_elevation_transition(0.01, "2", chars, site_elevation=6600.0, latitude=41.5)
        with pytest.raises(ValueError, match="latitude"):
            estimate_elevation_transition(0.01, "2", chars, site_elevation=6600.0)
        with pytest.raises(ValueError, match="Region 1"):
            estimate_elevation_transition(
                0.01, "1", chars, site_elevation=7000.0, latitude=38.0, lib=lib
            )
        with pytest.raises(ValueError, match="equation 7"):
            estimate_area_weighted(0.01, {"1": 10.0, "2": 40.0}, chars, lib=lib)
        with pytest.raises(ValueError, match="exactly two"):
            estimate_area_weighted(0.01, {"2": 50.0}, chars, lib=lib)
        with pytest.raises(ValueError, match="positive"):
            estimate_area_weighted(0.01, {"2": 0.0, "10": 0.0}, chars, lib=lib)
        with pytest.raises(ValueError, match="Nevada"):
            estimate_area_weighted(0.01, {"2": 1.0, "10": 1.0}, chars, lib=load_state("AZ"))
        # Region 1 has no 200-year equation.
        with pytest.raises(EquationsUnavailable):
            estimate_elevation_transition(
                0.002, "2", chars, site_elevation=7150.0, latitude=38.0, lib=lib
            )

    def test_equation_6_with_region_6_two_year_constant(self, lib):
        """Region 6's Q2 is printed 0, so eq. 6 at AEP 0.5 is Region 10's Q2 x 36/57."""
        from flowfreq.regression.nevada import estimate_area_weighted

        chars = {"DRNAREA": 57.0, "ELEV": 6500.0}
        w = estimate_area_weighted(0.5, {"6": 21.0, "10": 36.0}, chars, lib=lib)
        q2_10 = evaluate(lib.equation("10", 0.5), chars).flow_cfs
        assert q2_10 == pytest.approx(12 * 57**0.58, rel=1e-12)  # Table 14: Q = 12 AREA^0.58
        assert w.flow_cfs == pytest.approx(q2_10 * 36 / 57, rel=1e-12)
        assert w.blended and w.estimates["6"].flow_cfs == 0.0

    def test_equation_7_with_region_6_two_year_constant(self, lib):
        from flowfreq.regression.nevada import estimate_elevation_transition

        chars = {"DRNAREA": 50.0, "PRECIP": 20.0, "ELEV": 6500.0}
        mid = estimate_elevation_transition(
            0.5, "6", chars, site_elevation=7150.0, latitude=38.0, lib=lib
        )
        q1 = evaluate(lib.equation("1", 0.5), chars).flow_cfs
        assert mid.flow_cfs == pytest.approx(0.5 * q1, rel=1e-12)

    def test_constant_cannot_be_weighted_in_log_space(self, lib):
        from flowfreq.regression import evaluate_weighted

        chars = {"DRNAREA": 57.0, "ELEV": 6500.0}
        e6 = evaluate(lib.equation("6", 0.5), chars)
        e10 = evaluate(lib.equation("10", 0.5), chars)
        with pytest.raises(ValueError, match="6 AEP 0.5.*linear space"):
            evaluate_weighted([(0.4, e6), (0.6, e10)], space="log")


# Equivalent years of record, the last column of Tables 5, 6, 7, 9 (pp. 36, 38, 40, 45),
# 10 (p. 46) and 14 (p. 54), 2 to 100 years. Region 6's 2-year cell prints '--'.
EQUIVALENT_YEARS = {
    "1": [0.16, 0.62, 1.34, 2.50, 3.37, 4.19],
    "2": [0.96, 1.80, 3.07, 4.64, 5.47, 6.05],
    "3": [0.29, 0.49, 0.77, 1.23, 1.57, 1.92],
    "5": [0.21, 0.73, 1.69, 2.62, 3.26, 3.80],
    "6": [None, 0.233, 0.748, 2.52, 1.75, 0.794],
    "10": [0.618, 3.13, 3.45, 2.49, 3.22, 2.22],
}


class TestEquivalentYears:
    @pytest.mark.parametrize("region", sorted(EQUIVALENT_YEARS))
    def test_printed(self, lib, region):
        stored = [lib.equation(region, aep).equivalent_years for aep in AEPS]
        assert stored == EQUIVALENT_YEARS[region]

    def test_match_nss_snapshot(self, lib):
        """NSS's equivalentYears (snapshot 2026-10-01) agree in all 35 printed cells; NSS
        gives 0 for region 6's 2-year cell, which the report prints as '--'."""
        import json

        snap = REPO_ROOT / "data" / "nss_snapshots" / "NV_2026-10-01.json"
        if not snap.exists():
            pytest.skip("NSS snapshot not present")
        compared = 0
        for region in json.loads(snap.read_text(encoding="utf-8"))["regions"]:
            for rec in region["equations"]:
                eq = lib.equation(
                    NSS_CODE_TO_REGION[region["region_code"]],
                    NSS_STAT_TO_AEP[rec["statistic_code"]],
                )
                if eq.equivalent_years is None:
                    assert (eq.region_code, eq.aep, rec["equivalent_years"]) == ("6", 0.5, 0.0)
                else:
                    assert eq.equivalent_years == rec["equivalent_years"]
                    compared += 1
        assert compared == 35


# Data section, p. 147: systematic years N (p. 146), drainage area, mean basin elevation,
# and the station (top) and weighted (bottom) 5- to 100-year peaks, cfs.
GAGED = {
    "10243700 Cleve Creek near Ely": (
        29,
        {"DRNAREA": 31.8, "ELEV": 8770.0},
        [88, 130, 201, 270, 356],
        [89, 131, 208, 284, 368],
    ),
    "10244950 Steptoe Creek near Ely": (
        20,
        {"DRNAREA": 11.1, "ELEV": 8940.0},
        [42, 57, 79, 97, 117],
        [42, 58, 86, 109, 129],
    ),
}


class TestGagedSiteWeighting:
    @pytest.mark.parametrize("station", sorted(GAGED))
    def test_equation_1_reproduces_data_section(self, lib, station):
        """Within the rounding of the printed station values (the report weighted the
        unrounded ones): at most 1 cfs."""
        from flowfreq.regression.nevada import weight_gaged_estimate

        n, basin, station_q, weighted = GAGED[station]
        for aep, qs, qw in zip(AEPS[1:], station_q, weighted):
            est = evaluate(lib.equation("6", aep), basin)
            assert weight_gaged_estimate(qs, n, est) == pytest.approx(qw, abs=1.0), aep

    def test_formula(self, lib):
        from flowfreq.regression.nevada import weight_gaged_estimate

        est = evaluate(lib.equation("10", 0.01), {"DRNAREA": 50.0})
        e = 2.22
        assert weight_gaged_estimate(1000.0, 30, est) == pytest.approx(
            (1000.0 * 30 + est.flow_cfs * e) / (30 + e), rel=1e-15
        )
        assert weight_gaged_estimate(1000.0, 30, est, equivalent_years=10.0) == pytest.approx(
            (1000.0 * 30 + est.flow_cfs * 10.0) / 40.0, rel=1e-15
        )

    def test_region_6_two_year_needs_an_explicit_e(self, lib):
        """The data section's weighted 2-year values are below the station's (Cleve Creek
        45 to 40), consistent with Q=0 and an unprinted E of about 3 years."""
        from flowfreq.regression.nevada import weight_gaged_estimate

        est = evaluate(lib.equation("6", 0.5), GAGED["10243700 Cleve Creek near Ely"][1])
        with pytest.raises(ValueError, match="no printed equivalent years"):
            weight_gaged_estimate(45.0, 29, est)
        assert weight_gaged_estimate(45.0, 29, est, equivalent_years=3.6) == pytest.approx(
            45.0 * 29 / 32.6
        )

    @pytest.mark.parametrize(
        "kw, match",
        [
            ({"station_cfs": -1.0}, "station_cfs"),
            ({"station_cfs": math.nan}, "station_cfs"),
            ({"years_of_record": 0}, "years_of_record"),
            ({"equivalent_years": 0.0}, "equivalent_years"),
        ],
    )
    def test_errors(self, lib, kw, match):
        from flowfreq.regression.nevada import weight_gaged_estimate

        args = {"station_cfs": 100.0, "years_of_record": 20, **kw}
        est = evaluate(lib.equation("10", 0.01), {"DRNAREA": 50.0})
        with pytest.raises(ValueError, match=match):
            weight_gaged_estimate(args.pop("station_cfs"), args.pop("years_of_record"), est, **args)
