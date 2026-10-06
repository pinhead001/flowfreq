"""Nevada peak-flow regression equations (issue #104).

Source: Thomas, Hjalmarson and Waltemeyer, 1997, USGS WSP 2433, Tables 5, 6, 7, 9 (GLS) and
10, 14 (hybrid method). Every literal number below is copied from the report, or from a
dated live NSS response, never from ``NV.json`` itself.
"""

from __future__ import annotations

import math

import pytest

from flowfreq.regression import EquationsUnavailable, evaluate, load_state

AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01]
BASIN = {"DRNAREA": 50.0, "PRECIP": 20.0, "ELEV": 6500.0, "LAT_GAGE": 38.0}


@pytest.fixture(scope="module")
def lib():
    return load_state("NV")


def _sig3(x: float) -> float:
    return float(f"{x:.3g}")


class TestLibraryShape:
    def test_status_and_regions(self, lib):
        assert lib.status == "partial"
        assert lib.issue == 104
        assert lib.regions == ["1", "10", "2", "3", "5", "6"]
        for region in ("1", "2", "3", "5", "10"):
            assert lib.aeps(region) == AEPS
        assert lib.aeps("6") == AEPS[1:]
        assert len(lib.equations) == 35

    def test_region_6_two_year_is_not_stored(self, lib):
        """Table 10 prints the 2-year equation as Q = 0."""
        with pytest.raises(EquationsUnavailable):
            lib.equation("6", 0.5)

    @pytest.mark.parametrize(
        "region, n", [("1", 165), ("2", 108), ("3", 35), ("5", 37), ("6", 80), ("10", 104)]
    )
    def test_n_sites(self, lib, region, n):
        assert {e.n_sites for e in lib.equations if e.region_code == region} == {n}

    def test_hybrid_regions_carry_no_error_statistics(self, lib):
        for region in ("6", "10"):
            for eq in (e for e in lib.equations if e.region_code == region):
                assert eq.sep_log is None and eq.model_error_variance is None
        assert evaluate(lib.equation("10", 0.01), BASIN).interval is None

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
                continue
            ours = evaluate(lib.equation(region, NSS_STAT_TO_AEP[code]), basin).flow_cfs
            assert _sig3(ours) == est.value, (region, code, est.value, ours)
            compared += 1
    assert compared == 35


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
        # Region 6's 2-year equation is printed Q = 0 and is not stored.
        with pytest.raises(EquationsUnavailable):
            estimate_area_weighted(0.5, {"6": 21.0, "10": 36.0}, chars, lib=lib)
