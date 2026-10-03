"""Oregon (western) peak-flow regression equations (issue #38).

Source: Cooper, 2005, USGS SIR 2005-5116: Tables 10-12 (equations and error
statistics), Table 15 (ranges, station counts), the report's worked examples (pp.
45-46, Table 14) and Appendices A, D and G (per-gage regression estimates). Every
literal number below is copied from that report or from a dated live NSS response,
never from ``OR.json`` itself: the tests check the file against the source, not
against itself.
"""

from __future__ import annotations

import math

import pytest

from flowfreq.regression import OutOfRangeError, evaluate, load_state

AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002]


@pytest.fixture(scope="module")
def lib():
    return load_state("OR")


def _sig3(x: float) -> float:
    return float(f"{x:.3g}")


def _pct(sep_log: float) -> float:
    return 100 * math.sqrt(10 ** (math.log(10) * sep_log**2) - 1)


class TestLibraryShape:
    def test_status_and_regions(self, lib):
        # verified per library.STATUSES: everything both reports publish is stored
        # and cross-checked against their own per-gage estimates and worked
        # examples (NSS serves only region 1); no covariance is published.
        assert lib.status == "verified"
        assert lib.issue == 38
        assert lib.regions == ["1", "2A", "2B", "E1", "E2", "E3", "E4", "E5", "E6"]
        for region in lib.regions:
            assert lib.aeps(region) == AEPS  # no 200-year (0.5%) equation
        assert len(lib.equations) == 63

    def test_every_equation_is_cited(self, lib):
        for eq in lib.equations:
            if eq.region_code.startswith("E"):
                continue  # eastern Oregon: TestEasternOregon
            assert "2005-5116" in eq.citation.publication
            assert "Table 15" in eq.citation.table
            assert eq.covariance is None and eq.avp is None  # not published
            assert all(v.transform == "log10" for v in eq.variables)
            assert "ORREG2" not in {v.code for v in eq.variables}

    @pytest.mark.parametrize("region, n", [("1", 91), ("2A", 107), ("2B", 178)])
    def test_n_sites(self, lib, region, n):
        assert {e.n_sites for e in lib.equations if e.region_code == region} == {n}

    def test_limits_from_table_15(self, lib):
        def limits(region):
            return {v.code: (v.minimum, v.maximum) for v in lib.equation(region, 0.01).variables}

        assert limits("1") == {
            "DRNAREA": (0.28, 673.0),
            "I24H2Y": (2.52, 5.79),
            "JANMAXT2K": (42.4, 53.9),
            "WATCAPORC": (0.10, 0.23),
            "SOILPERM": (0.72, 4.76),
        }
        assert limits("2A") == {
            "DRNAREA": (0.22, 3940.0),
            "BSLOPD": (6.24, 28.0),
            "I24H2Y": (1.72, 4.34),
            "JANMINT2K": (20.5, 34.0),
            "JANMAXT2K": (33.9, 47.3),
        }
        assert limits("2B") == {
            "DRNAREA": (0.37, 7270.0),
            "BSLOPD": (5.62, 28.3),
            "I24H2Y": (1.53, 4.48),
        }

    def test_out_of_range_raises(self, lib):
        with pytest.raises(OutOfRangeError, match="I24H2Y"):
            evaluate(lib.equation("2B", 0.01), {"DRNAREA": 10.0, "BSLOPD": 10.0, "I24H2Y": 6.0})


class TestAgainstTables10To12:
    @pytest.mark.parametrize(
        "region, aep, sem_pct, pred_pct",
        [
            ("1", 0.5, 25.5, 26.8),
            ("1", 0.002, 30.0, 32.6),
            ("2A", 0.5, 37.1, 38.7),
            ("2A", 0.01, 31.6, 34.4),
            ("2B", 0.04, 33.0, 34.1),
            ("2B", 0.002, 37.7, 39.1),
        ],
    )
    def test_error_statistics(self, lib, region, aep, sem_pct, pred_pct):
        eq = lib.equation(region, aep)
        assert _pct(eq.sep_log) == pytest.approx(pred_pct, abs=1e-3)
        assert _pct(math.sqrt(eq.model_error_variance)) == pytest.approx(sem_pct, abs=1e-3)

    def test_region_2b_power_form(self, lib):
        b = {"DRNAREA": 100.0, "BSLOPD": 15.0, "I24H2Y": 3.0}
        est = evaluate(lib.equation("2B", 0.1), b)
        assert est.flow_cfs == pytest.approx(
            18.49 * 100.0**0.9064 * 15.0**0.4688 * 3.0**0.6937, rel=1e-12
        )


class TestWorkedExamples:
    def test_lobster_creek_region_1_q100(self, lib):
        """Report p. 45, eq. 10: Q100 = 10,200 cfs."""
        site = {
            "DRNAREA": 58.3,
            "I24H2Y": 3.69,
            "JANMAXT2K": 47.6,
            "SOILPERM": 2.51,
            "WATCAPORC": 0.134,
        }
        assert _sig3(evaluate(lib.equation("1", 0.01), site).flow_cfs) == 10200

    def test_quartz_creek_transition_zone(self, lib):
        """Report pp. 45-46, Table 14 (mean elevation 2,970 ft).

        The report gives 7,690 and 8,380 cfs and labels them 2A and 2B; the equations
        give them the other way round (the text swaps the labels). With the labels
        corrected, the report's blend is its equation 8 as printed.
        """
        site = {
            "DRNAREA": 42.1,
            "BSLOPD": 23.8,
            "I24H2Y": 2.84,
            "JANMINT2K": 31.0,
            "JANMAXT2K": 44.3,
        }
        q2a = evaluate(lib.equation("2A", 0.01), site).flow_cfs
        q2b = evaluate(lib.equation("2B", 0.01), site).flow_cfs
        assert q2b == pytest.approx(7690, rel=0.005)
        assert q2a == pytest.approx(8380, rel=0.01)  # Table 14 temperatures are rounded
        e = 2970.0
        qt = q2b * (3125 - e) / 250 + q2a * (e - 2875) / 250
        assert qt == pytest.approx(7950, rel=0.005)


QUARTZ_CREEK = {  # Table 14 (report p. 46)
    "DRNAREA": 42.1,
    "BSLOPD": 23.8,
    "I24H2Y": 2.84,
    "JANMINT2K": 31.0,
    "JANMAXT2K": 44.3,
}


class TestRegion2TransitionZone:
    """Equations 7 and 8 (report pp. 31, 44), via flowfreq.regression.oregon."""

    @pytest.mark.parametrize(
        "elev, w2a",
        [
            (2000.0, 0.0),
            (2875.0, 0.0),  # zone's lower edge: Region 2B alone
            (2970.0, 0.38),  # Quartz Creek: (2,970 - 2,875) / 250
            (3000.0, 0.5),
            (3125.0, 1.0),  # zone's upper edge: Region 2A alone
            (5000.0, 1.0),
        ],
    )
    def test_equation_8_weights(self, elev, w2a):
        from flowfreq.regression.oregon import region2_weights

        w = region2_weights(elev)
        assert w["2A"] == pytest.approx(w2a, abs=1e-15)
        assert w["2B"] == pytest.approx(1 - w2a, abs=1e-15)

    @pytest.mark.parametrize(
        "width, lower, upper",
        [(500.0, 2750.0, 3250.0), (350.0, 2825.0, 3175.0)],  # p. 44 text; Table 13 headnote
    )
    def test_equation_7_zone_edges(self, width, lower, upper):
        from flowfreq.regression.oregon import region2_weights

        assert region2_weights(lower, width)["2A"] == 0.0
        assert region2_weights(upper, width)["2A"] == 1.0
        assert region2_weights((lower + upper) / 2, width)["2A"] == pytest.approx(0.5)

    @pytest.mark.parametrize("width", [0.0, -250.0, math.inf, math.nan])
    def test_invalid_width_raises(self, width):
        from flowfreq.regression.oregon import region2_weights

        with pytest.raises(ValueError, match="width"):
            region2_weights(3000.0, width)

    def test_report_arithmetic_for_quartz_creek(self):
        """Eq. 8 with the report's own component values reproduces its 7,950 cfs.

        The text calls 7,690 the Region 2A value, but eq. 8 puts it in the Q2B slot,
        and the equations agree it is 2B's (see test_quartz_creek_transition_zone).
        """
        from flowfreq.regression.oregon import region2_weights

        w = region2_weights(2970.0)
        assert _sig3(7690 * w["2B"] + 8380 * w["2A"]) == 7950

    def test_quartz_creek_blend(self, lib):
        from flowfreq.regression import evaluate_weighted
        from flowfreq.regression.oregon import estimate_region2

        est = estimate_region2(0.01, QUARTZ_CREEK, mean_elevation=2970.0, lib=lib)
        assert est.blended
        q2a = evaluate(lib.equation("2A", 0.01), QUARTZ_CREEK).flow_cfs
        q2b = evaluate(lib.equation("2B", 0.01), QUARTZ_CREEK).flow_cfs
        assert est.estimates["2A"].flow_cfs == q2a and est.estimates["2B"].flow_cfs == q2b
        # Equation 8 exactly, and evaluate_weighted(space="linear") expresses it exactly.
        assert est.flow_cfs == pytest.approx(
            q2b * (3125 - 2970) / 250 + q2a * (2970 - 2875) / 250, rel=1e-15
        )
        assert est.flow_cfs == evaluate_weighted(
            [(0.62, est.estimates["2B"]), (0.38, est.estimates["2A"])], space="linear"
        )
        # The report's 7,950 used Table 14's rounded temperatures (its 2A value is 8,380
        # against the equation's 8,428), so agreement is to that rounding.
        assert est.flow_cfs == pytest.approx(7950, rel=0.005)

    def test_elevation_from_characteristics(self, lib):
        from flowfreq.regression.oregon import estimate_region2

        a = estimate_region2(0.01, {**QUARTZ_CREEK, "ELEV": 2970.0}, lib=lib)
        b = estimate_region2(0.01, QUARTZ_CREEK, mean_elevation=2970.0, lib=lib)
        assert a.flow_cfs == b.flow_cfs

    def test_missing_elevation_raises(self, lib):
        from flowfreq.regression.oregon import estimate_region2

        with pytest.raises(KeyError, match="ELEV"):
            estimate_region2(0.01, QUARTZ_CREEK, lib=lib)

    def test_rejects_other_state(self):
        from flowfreq.regression.oregon import estimate_region2

        with pytest.raises(ValueError, match="Oregon"):
            estimate_region2(0.01, QUARTZ_CREEK, mean_elevation=2970.0, lib=load_state("WA"))

    def test_mckenzie_river_above_zone_is_region_2a(self, lib):
        """14162500, Appendix G mean elevation 3,850 ft: R = 99,900 cfs (p. 45)."""
        from flowfreq.regression.oregon import estimate_region2

        site = {"DRNAREA": 930.0, "BSLOPD": 16.4, "I24H2Y": 3.20, "JANMINT2K": 25.8}
        site["JANMAXT2K"] = 40.6
        est = estimate_region2(0.01, site, mean_elevation=3850.0, lib=lib)
        assert not est.blended and set(est.estimates) == {"2A"}
        assert est.flow_cfs == pytest.approx(99900, rel=0.01)

    def test_marks_creek_below_zone_is_region_2b(self, lib):
        """14312300, Appendix G mean elevation 752 ft: Appendix D R100 = 176 cfs.

        Appendix A lists this gage as 2A; its elevation puts it in 2B, which is what
        Appendix D's regression estimates use. Only 2B's characteristics are needed.
        """
        from flowfreq.regression.oregon import estimate_region2

        site = {"DRNAREA": 1.31, "BSLOPD": 11.5, "I24H2Y": 1.79}
        est = estimate_region2(0.01, site, mean_elevation=752.0, lib=lib, allow_extrapolation=True)
        assert set(est.estimates) == {"2B"} and est.weights["2B"] == 1.0
        assert est.flow_cfs == pytest.approx(176, rel=0.01)

    def test_continuous_across_zone_edges(self, lib):
        from flowfreq.regression.oregon import estimate_region2

        q = {
            e: estimate_region2(0.01, QUARTZ_CREEK, mean_elevation=e, lib=lib).flow_cfs
            for e in (2874.999, 2875.0, 2875.001, 3124.999, 3125.0, 3125.001)
        }
        assert q[2874.999] == q[2875.0] == pytest.approx(q[2875.001], rel=1e-5)
        assert q[3125.001] == q[3125.0] == pytest.approx(q[3124.999], rel=1e-5)


# Appendix D "R" (regression) estimates at gages, with Appendix G characteristics.
APPENDIX_D = [
    # McKenzie River near Vida (14162500), Region 2A, 100-year: R = 99,900 (p. 45)
    (
        "2A",
        0.01,
        {"DRNAREA": 930.0, "BSLOPD": 16.4, "I24H2Y": 3.20, "JANMINT2K": 25.8, "JANMAXT2K": 40.6},
        99900,
    ),
    # Marks Creek near Roseburg (14312300): Appendix A says 2A, but Appendix D uses 2B
    (
        "2B",
        0.01,
        {"DRNAREA": 1.31, "BSLOPD": 11.5, "I24H2Y": 1.79},
        176,
    ),
]


@pytest.mark.parametrize("region, aep, site, r", APPENDIX_D)
def test_appendix_d_regression_estimates(lib, region, aep, site, r):
    est = evaluate(lib.equation(region, aep), site, allow_extrapolation=True)
    assert est.flow_cfs == pytest.approx(r, rel=0.01)


# Live NSS response (GC729, ORREG2=729) for the Lobster Creek basin, fetched 2026-09-26.
LOBSTER = {"DRNAREA": 58.3, "I24H2Y": 3.69, "JANMAXT2K": 47.6, "WATCAPORC": 0.134, "SOILPERM": 2.51}
NSS_RECORDED = [
    (0.5, 4180),
    (0.2, 5790),
    (0.1, 6850),
    (0.04, 8200),
    (0.02, 9200),
    (0.01, 10200),
    (0.002, 12400),
]


@pytest.mark.parametrize("aep, nss", NSS_RECORDED)
def test_region_1_matches_recorded_nss(lib, aep, nss):
    assert _sig3(evaluate(lib.equation("1", aep), LOBSTER).flow_cfs) == nss


NSS_STAT_TO_AEP = {
    "PK50AEP": 0.5,
    "PK20AEP": 0.2,
    "PK10AEP": 0.1,
    "PK4AEP": 0.04,
    "PK2AEP": 0.02,
    "PK1AEP": 0.01,
    "PK0_2AEP": 0.002,
}


@pytest.mark.requires_network
@pytest.mark.parametrize(
    "basin",
    [
        {"DRNAREA": 5.0, "I24H2Y": 3.0, "JANMAXT2K": 45.0, "WATCAPORC": 0.15, "SOILPERM": 1.5},
        LOBSTER,
        {"DRNAREA": 500.0, "I24H2Y": 5.5, "JANMAXT2K": 53.0, "WATCAPORC": 0.22, "SOILPERM": 4.5},
    ],
)
def test_region_1_matches_live_nss(lib, basin):
    """Region 1 against NSS Oregon (region 41), GC729, citation 118.

    Only GC729 resolves live (with the gating parameter ORREG2 = 729). NSS lists
    Regions 2A and 2B (GC730, GC731) but returns no estimate for them.
    """
    from flowfreq.streamstats import Characteristic, estimate_flow_statistics

    chars = {k: Characteristic(k, k, "", v, "") for k, v in {**basin, "ORREG2": 729}.items()}
    results, _skipped = estimate_flow_statistics("OR", chars, ["PFS"])
    compared = 0
    for r in results:
        if r.region_code != "GC729":
            continue
        for code, est in r.estimates.items():
            ours = evaluate(lib.equation("1", NSS_STAT_TO_AEP[code]), basin).flow_cfs
            assert _sig3(ours) == est.value, (code, est.value, ours)
            compared += 1
    assert compared == 7


# ---------------------------------------------------------------------------
# Eastern Oregon: Cooper, 2006, OWRD Open File Report SW 06-001 (a State of Oregon
# report, not USGS, and not in NSS). Literals below are copied from that report:
# Tables 21-27, the worked examples (pp. 78-81) and Appendices A, D and G.
# ---------------------------------------------------------------------------


class TestEasternOregon:
    def test_citation_is_explicitly_non_usgs(self, lib):
        for eq in lib.equations:
            if not eq.region_code.startswith("E"):
                continue
            pub = eq.citation.publication
            assert "SW 06-001" in pub and "Oregon Water Resources Department" in pub
            assert "not a U.S. Geological Survey publication" in pub and "not in NSS" in pub
            assert eq.citation.year == 2006 and "Table 27" in eq.citation.table
            assert eq.covariance is None and eq.avp is None  # not published
            assert all(v.transform == "log10" for v in eq.variables)
        src = [s for s in lib.source_reports if "SW 06-001" in s["citation"]]
        assert len(src) == 1 and src[0]["nss"] is None

    @pytest.mark.parametrize(
        "region, n", [("E1", 41), ("E2", 58), ("E3", 76), ("E4", 51), ("E5", 28), ("E6", 22)]
    )
    def test_n_sites_from_table_27(self, lib, region, n):
        assert {e.n_sites for e in lib.equations if e.region_code == region} == {n}

    def test_limits_from_table_27(self, lib):
        def limits(region):
            return {v.code: (v.minimum, v.maximum) for v in lib.equation(region, 0.01).variables}

        assert limits("E1") == {
            "DRNAREA": (4.64, 525.0),
            "BSLOPD": (5.03, 19.3),
            "ELEV": (2770.0, 6060.0),
        }
        assert limits("E2") == {
            "DRNAREA": (0.293, 1630.0),
            "JANAVPRE2K": (1.52, 7.97),
            "JULAVPRE2K": (0.234, 1.28),
            "WATCAPORC": (0.102, 0.196),
        }
        assert limits("E3") == {"DRNAREA": (0.240, 3870.0)}
        assert limits("E4") == {
            "DRNAREA": (1.70, 7630.0),
            "BSLOPD": (5.35, 22.6),
            "JANAVPRE2K": (1.37, 5.28),
            "JANMINT2K": (16.4, 23.1),
            "STATSGODEP": (21.3, 50.7),
        }
        assert limits("E5") == {
            "DRNAREA": (0.609, 2992.0),
            "ELEV": (4850.0, 7050.0),
            "JANAVPRE2K": (1.78, 4.46),
            "JULAVPRE2K": (0.197, 0.915),
        }
        assert limits("E6") == {"DRNAREA": (0.152, 453.0), "ASPECT": (133.0, 244.0)}

    def test_scaled_terms_keep_published_units(self, lib):
        elev = {v.code: v for v in lib.equation("E1", 0.01).variables}["ELEV"]
        assert (elev.scale, elev.offset, elev.units) == (0.001, 0.0, "ft")
        aspect = {v.code: v for v in lib.equation("E6", 0.01).variables}["ASPECT"]
        assert (aspect.scale, aspect.offset, aspect.units) == (0.01, 0.0, "degrees")

    def test_out_of_range_raises(self, lib):
        with pytest.raises(OutOfRangeError, match="ASPECT"):
            evaluate(lib.equation("E6", 0.01), {"DRNAREA": 10.0, "ASPECT": 300.0})

    @pytest.mark.parametrize(
        "region, aep, site, printed",
        [
            (
                "E1",
                0.5,
                {"DRNAREA": 50.0, "BSLOPD": 10.0, "ELEV": 4500.0},
                0.7516 * 50.0**0.8787 * 10.0**1.984 * 4.5**-1.069,
            ),
            (
                "E2",
                0.002,
                {"DRNAREA": 50.0, "JANAVPRE2K": 3.0, "JULAVPRE2K": 0.5, "WATCAPORC": 0.15},
                702.7 * 50.0**0.7407 * 3.0**0.5300 * 0.5**-1.348 * 0.15**1.330,
            ),
            ("E3", 0.04, {"DRNAREA": 50.0}, 61.90 * 50.0**0.7415),
            (
                "E4",
                0.2,
                {
                    "DRNAREA": 50.0,
                    "BSLOPD": 10.0,
                    "JANAVPRE2K": 3.0,
                    "JANMINT2K": 20.0,
                    "STATSGODEP": 40.0,
                },
                11.54 * 50.0**0.8576 * 10.0**-0.8858 * 3.0**2.025 * 20.0**1.983 * 40.0**-1.530,
            ),
            (
                "E5",
                0.01,
                {"DRNAREA": 50.0, "ELEV": 6000.0, "JANAVPRE2K": 3.0, "JULAVPRE2K": 0.5},
                0.1174 * 50.0**0.7737 * 6.0**4.505 * 3.0**-1.781 * 0.5**-0.8428,
            ),
            ("E6", 0.1, {"DRNAREA": 50.0, "ASPECT": 180.0}, 112.5 * 50.0**0.7873 * 1.8**-2.664),
        ],
    )
    def test_power_form(self, lib, region, aep, site, printed):
        """Tables 21-26 as printed, (Elev/1,000) and (Aspect/100) included."""
        assert evaluate(lib.equation(region, aep), site).flow_cfs == pytest.approx(
            printed, rel=1e-12
        )

    @pytest.mark.parametrize(
        "region, aep, sem_pct, pred_pct",
        [
            ("E1", 0.5, 46.3, 49.1),
            ("E2", 0.01, 35.6, 41.8),
            ("E3", 0.002, 52.2, 54.3),
            ("E4", 0.1, 40.5, 44.7),
            ("E5", 0.04, 31.7, 37.8),
            ("E6", 0.5, 94.6, 104.0),
        ],
    )
    def test_error_statistics(self, lib, region, aep, sem_pct, pred_pct):
        eq = lib.equation(region, aep)
        assert _pct(eq.sep_log) == pytest.approx(pred_pct, abs=1e-3)
        assert _pct(math.sqrt(eq.model_error_variance)) == pytest.approx(sem_pct, abs=1e-3)

    def test_west_birch_creek_region_2_q100(self, lib):
        """Report pp. 78-79, eq. 9: Q100 = 2,470 cfs."""
        site = {"DRNAREA": 86.25, "JANAVPRE2K": 2.85, "JULAVPRE2K": 0.647, "WATCAPORC": 0.128}
        assert _sig3(evaluate(lib.equation("E2", 0.01), site).flow_cfs) == 2470

    def test_middle_fork_john_day_transfer(self, lib):
        """Report p. 81, eq. 10 with Table 24's Q100 area exponent: 7,350 cfs."""
        eq = lib.equation("E4", 0.01)
        b = dict(zip((v.code for v in eq.variables), eq.coefficients))["DRNAREA"]
        assert b == 0.8375
        assert _sig3(5200 * (791 / 523) ** b) == 7350

    def test_region_5_uses_july_precipitation(self, lib):
        """Table 25's headnote says Mn Jul T, but its equations and Table 27 say Jul P.

        Gage 10370000 (Appendix G: Jul P 0.423 in., Mn Jul T 45.3 deg F): R100 = 2,180.
        """
        site = {"DRNAREA": 66.2, "ELEV": 6240.0, "JANAVPRE2K": 3.82, "JULAVPRE2K": 0.423}
        assert _sig3(evaluate(lib.equation("E5", 0.01), site).flow_cfs) == 2180
        wrong = {**site, "JULAVPRE2K": 45.3}
        est = evaluate(lib.equation("E5", 0.01), wrong, allow_extrapolation=True)
        assert est.flow_cfs < 0.1 * 2180


# Appendix G characteristics and Appendix D 'R' rows (2, 5, 10, 25, 50, 100, 500 years),
# with each gage's Appendix A flood region.
EAST_APPENDIX_D = [
    (
        "E1",
        "11504000",
        {"DRNAREA": 71.8, "BSLOPD": 8.62, "ELEV": 5440.0},
        [377, 551, 670, 822, 938, 1060, 1350],
    ),
    (
        "E2",
        "14010000",
        {"DRNAREA": 61.8, "JANAVPRE2K": 6.88, "JULAVPRE2K": 0.896, "WATCAPORC": 0.145},
        [883, 1310, 1630, 2080, 2430, 2790, 3680],
    ),
    ("E3", "13213800", {"DRNAREA": 53.0}, [437, 711, 912, 1180, 1380, 1590, 2100]),
    (
        # Strawberry Creek above Slide Creek near Prairie City (report p. 78: R100 = 169).
        "E4",
        "14037500",
        {
            "DRNAREA": 6.97,
            "BSLOPD": 22.6,
            "JANAVPRE2K": 5.28,
            "JANMINT2K": 16.8,
            "STATSGODEP": 40.3,
        },
        [79.8, 106, 121, 141, 155, 169, 200],
    ),
    (
        "E5",
        "10370000",
        {"DRNAREA": 66.2, "ELEV": 6240.0, "JANAVPRE2K": 3.82, "JULAVPRE2K": 0.423},
        [360, 680, 960, 1390, 1770, 2180, 3330],
    ),
    (
        "E6",
        "10396000",
        {"DRNAREA": 206.0, "ASPECT": 200.0},
        [391, 813, 1180, 1720, 2170, 2650, 3890],
    ),
]


@pytest.mark.parametrize("region, station, site, r_values", EAST_APPENDIX_D)
def test_eastern_appendix_d_regression_estimates(lib, region, station, site, r_values):
    for aep, r in zip(AEPS, r_values):
        est = evaluate(lib.equation(region, aep), site)
        if region == "E4" and aep == 0.1:
            # Appendix D's Region 4 10-year R values sit a constant ~3.7% above Table 24's
            # printed Q(10) equation at all 51 gages (see OR.json notes); Table 24 is stored.
            assert r / est.flow_cfs == pytest.approx(1.037, abs=0.012), station
            continue
        assert est.flow_cfs == pytest.approx(r, rel=0.01), (station, aep)
