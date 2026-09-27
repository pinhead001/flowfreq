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
        assert lib.status == "partial"
        assert lib.issue == 38
        assert lib.regions == ["1", "2A", "2B"]
        for region in lib.regions:
            assert lib.aeps(region) == AEPS  # no 200-year (0.5%) equation
        assert len(lib.equations) == 21

    def test_every_equation_is_cited(self, lib):
        for eq in lib.equations:
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
