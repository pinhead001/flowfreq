"""Washington peak-flow regression equations (issue #37).

Source: Mastin, Konrad, Veilleux and Tecca, 2016, USGS SIR 2016-5118 ver. 1.2
(November 2017): Table 6 (equations, n, Sp, ranges), Table 7 (model error variance,
t, covariance), and the report's Flood Q Tools workbook (full-precision covariance,
Regions 1 and 2). Every literal number below is copied from that report, from its
worked example (p. 43), or from a dated live NSS response, never from ``WA.json``
itself: the tests check the file against the source, not against itself.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from flowfreq.regression import OutOfRangeError, evaluate, load_state

AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]

# In range for every region (intersection of Table 6's ranges).
BASIN = {"DRNAREA": 100.0, "PRECPRIS10": 45.0, "CANOPY_PCT": 50.0}


@pytest.fixture(scope="module")
def lib():
    return load_state("WA")


def _sig3(x: float) -> float:
    return float(f"{x:.3g}")


class TestLibraryShape:
    def test_status_and_regions(self, lib):
        assert lib.status == "partial"
        assert lib.issue == 37
        assert lib.regions == ["1", "2", "3", "4"]
        for region in lib.regions:
            assert lib.aeps(region) == AEPS
        assert len(lib.equations) == 32

    def test_every_equation_is_cited(self, lib):
        for eq in lib.equations:
            assert "2016-5118" in eq.citation.publication
            assert "ver. 1.2" in eq.citation.publication
            assert "Table 6" in eq.citation.table and "Table 7" in eq.citation.table
            assert eq.model_error_variance is not None and eq.sep_log is not None
            assert eq.avp is None  # the report publishes no AVP

    def test_covariance_present_except_region_4(self, lib):
        for eq in lib.equations:
            if eq.region_code == "4":
                assert eq.covariance is None  # see WA.json notes: basis inconsistent
                continue
            u = np.asarray(eq.covariance)
            assert np.allclose(u, u.T)
            assert np.all(np.linalg.eigvalsh(u) > 0)

    @pytest.mark.parametrize("region, n", [("1", 93), ("2", 89), ("3", 142), ("4", 139)])
    def test_n_sites(self, lib, region, n):
        assert {e.n_sites for e in lib.equations if e.region_code == region} == {n}

    @pytest.mark.parametrize(
        "region, t_printed", [("1", 1.662), ("2", 1.6630), ("3", 1.6559), ("4", 1.656)]
    )
    def test_n_sites_reproduce_table7_t(self, lib, region, t_printed):
        eq = lib.equation(region, 0.01)
        df = eq.n_sites - len(eq.variables) - 1
        assert round(stats.t.ppf(0.95, df), 3) == pytest.approx(round(t_printed, 3))

    def test_limits_from_table_6(self, lib):
        def limits(region):
            return {v.code: (v.minimum, v.maximum) for v in lib.equation(region, 0.01).variables}

        assert limits("1") == {
            "DRNAREA": (0.25, 3304.0),
            "PRECPRIS10": (9.82, 52.45),
            "CANOPY_PCT": (0.0, 77.4),
        }
        assert limits("2") == {
            "DRNAREA": (0.42, 1322.0),
            "PRECPRIS10": (8.86, 84.2),
            "CANOPY_PCT": (0.0, 81.8),
        }
        assert limits("3") == {"DRNAREA": (0.08, 2605.0), "PRECPRIS10": (33.29, 168.0)}
        assert limits("4") == {"DRNAREA": (0.18, 2230.0), "PRECPRIS10": (11.94, 186.6)}

    def test_out_of_range_raises(self, lib):
        with pytest.raises(OutOfRangeError, match="PRECPRIS10"):
            evaluate(lib.equation("3", 0.01), {"DRNAREA": 10.0, "PRECPRIS10": 20.0})


class TestAgainstTable6:
    """Evaluate the report's printed equation forms directly."""

    @pytest.mark.parametrize(
        "region, aep, printed",
        [
            (
                "1",
                0.5,
                lambda b: 3.846
                * b["DRNAREA"] ** 0.745
                * 10 ** (0.032 * b["PRECPRIS10"])
                / 10 ** (0.0078 * b["CANOPY_PCT"]),
            ),
            (
                "1",
                0.002,
                lambda b: 193.20
                * b["DRNAREA"] ** 0.624
                * 10 ** (0.020 * b["PRECPRIS10"])
                / 10 ** (0.0143 * b["CANOPY_PCT"]),
            ),
            (
                "2",
                0.04,
                lambda b: 26.485
                * b["DRNAREA"] ** 0.833
                * 10 ** (0.0174 * b["PRECPRIS10"])
                / 10 ** (0.0105 * b["CANOPY_PCT"]),
            ),
            ("3", 0.1, lambda b: 0.0794 * b["DRNAREA"] ** 0.895 * b["PRECPRIS10"] ** 1.670),
            ("3", 0.002, lambda b: 0.281 * b["DRNAREA"] ** 0.900 * b["PRECPRIS10"] ** 1.542),
            ("4", 0.5, lambda b: 15.346 * b["DRNAREA"] ** 0.911 * 10 ** (0.0073 * b["PRECPRIS10"])),
            (
                "4",
                0.01,
                lambda b: 91.411 * b["DRNAREA"] ** 0.921 * 10 ** (0.0040 * b["PRECPRIS10"]),
            ),
        ],
    )
    def test_power_form(self, lib, region, aep, printed):
        est = evaluate(lib.equation(region, aep), BASIN)
        assert est.flow_cfs == pytest.approx(printed(BASIN), rel=1e-12)

    @pytest.mark.parametrize(
        "region, aep, mev, sp_pct",
        [
            ("1", 0.5, 0.113, 95.04),
            ("1", 0.01, 0.108, 93.55),
            ("2", 0.002, 0.155, 119.59),
            ("3", 0.1, 0.034, 45.6),
            ("4", 0.01, 0.047, 54.21),
        ],
    )
    def test_error_statistics(self, lib, region, aep, mev, sp_pct):
        eq = lib.equation(region, aep)
        assert eq.model_error_variance == mev
        # sep_log is the percent Sp converted to log10 units
        assert 100 * math.sqrt(10 ** (math.log(10) * eq.sep_log**2) - 1) == pytest.approx(
            sp_pct, abs=1e-3
        )


class TestCovariance:
    def test_region_3_from_table_7(self, lib):
        assert lib.equation("3", 0.01).covariance == (
            (6.826e-2, 1.469e-4, -3.471e-2),
            (1.469e-4, 6.309e-4, -5.263e-4),
            (-3.471e-2, -5.263e-4, 1.819e-2),
        )

    def test_regions_1_2_round_to_table_7(self, lib):
        # Full-precision workbook values round to Table 7's printed upper triangle.
        cases = {
            ("1", 0.01): [2.469e-2, -3.136e-3, -7.782e-4, 2.151e-4, 1.849e-3, 1.634e-5]
            + [-1.560e-5, 4.727e-5, -1.970e-5, 1.202e-5],
            ("2", 0.5): [1.635e-2, -2.438e-3, -1.267e-4, -1.153e-4, 1.546e-3, 1.063e-6]
            + [5.239e-7, 5.535e-6, -1.668e-6, 3.360e-6],
        }
        for (region, aep), printed in cases.items():
            u = lib.equation(region, aep).covariance
            upper = [u[i][j] for i in range(4) for j in range(i, 4)]
            assert [float(f"{x:.4g}") for x in upper] == printed

    def test_worked_example_report_p43(self, lib):
        """Region 2, AEP 0.01, USGS 12488500 American River near Nile (report p. 43)."""
        site = {"DRNAREA": 79.22, "PRECPRIS10": 67.6, "CANOPY_PCT": 68.1}
        eq = lib.equation("2", 0.01)
        x = np.array([1.0, math.log10(79.22), 67.6, 68.1])
        assert x @ np.asarray(eq.covariance) @ x == pytest.approx(0.00868359, abs=5e-9)
        est = evaluate(eq, site)
        assert _sig3(est.flow_cfs) == 3190
        assert est.flow_cfs == pytest.approx(3187.8344, rel=1e-5)
        s_i = math.sqrt(0.116 + 0.00868359)
        assert s_i == pytest.approx(0.35311, abs=5e-6)
        assert 10 ** (1.6630 * s_i) == pytest.approx(3.8656, abs=5e-5)
        assert round(est.interval[0]) == 825
        assert _sig3(est.interval[1]) == 12300

    def test_region_4_falls_back_to_average_sep(self, lib):
        est = evaluate(lib.equation("4", 0.01), BASIN)
        assert est.interval_method == "average standard error of prediction"


# Live NSS response for Goat Creek (DRNAREA 412.0, PRECPRIS10 45.62, CANOPY_PCT 45.242),
# recorded in tests/fixtures/streamstats_responses.py (Region 1 = GC1750).
GOAT_CREEK = {"DRNAREA": 412.0, "PRECPRIS10": 45.62, "CANOPY_PCT": 45.242}


@pytest.mark.parametrize("aep, nss", [(0.5, 4370.0), (0.2, 6050.0)])
def test_recorded_nss_goat_creek(lib, aep, nss):
    assert _sig3(evaluate(lib.equation("1", aep), GOAT_CREEK).flow_cfs) == nss


# Live NSS responses for BASIN, fetched 2026-09-26 (NSS rounds to 3 significant figures).
NSS_RECORDED = [
    ("1", 0.5, 1330),
    ("1", 0.01, 4020),
    ("2", 0.5, 896),
    ("2", 0.002, 4090),
    ("3", 0.1, 2820),
    ("3", 0.002, 6280),
    ("4", 0.5, 2170),
    ("4", 0.01, 9620),
]


@pytest.mark.parametrize("region, aep, nss", NSS_RECORDED)
def test_worked_example_matches_recorded_nss(lib, region, aep, nss):
    assert _sig3(evaluate(lib.equation(region, aep), BASIN).flow_cfs) == nss


NSS_CODE_TO_REGION = {"GC1750": "1", "GC1751": "2", "GC1752": "3", "GC1753": "4"}
NSS_STAT_TO_AEP = {
    "PK50AEP": 0.5,
    "PK20AEP": 0.2,
    "PK10AEP": 0.1,
    "PK4AEP": 0.04,
    "PK2AEP": 0.02,
    "PK1AEP": 0.01,
    "PK0_5AEP": 0.005,
    "PK0_2AEP": 0.002,
}


@pytest.mark.requires_network
@pytest.mark.parametrize(
    "basin",
    [
        {"DRNAREA": 5.0, "PRECPRIS10": 35.0, "CANOPY_PCT": 10.0},
        BASIN,
        {"DRNAREA": 1000.0, "PRECPRIS10": 52.0, "CANOPY_PCT": 75.0},
    ],
)
def test_matches_live_nss(lib, basin):
    """Every region x AEP against NSS Washington (region 51), citation 150.

    NSS returns three significant figures, so the check is exact agreement after
    rounding flowfreq's value the same way; the relative difference is bounded by
    that rounding (at most 0.5% for a leading digit of 1).
    """
    from flowfreq.streamstats import Characteristic, estimate_flow_statistics

    chars = {k: Characteristic(k, k, "", v, "") for k, v in basin.items()}
    results, skipped = estimate_flow_statistics("WA", chars, ["PFS"])
    assert not skipped
    compared = 0
    for r in results:
        region = NSS_CODE_TO_REGION[r.region_code]
        for code, est in r.estimates.items():
            ours = evaluate(lib.equation(region, NSS_STAT_TO_AEP[code]), basin).flow_cfs
            assert _sig3(ours) == est.value, (region, code, est.value, ours)
            assert abs(ours - est.value) / est.value < 5e-3
            compared += 1
    assert compared == 32
