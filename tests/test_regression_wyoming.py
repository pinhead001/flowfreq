"""Wyoming peak-flow regression equations (issue #98).

Source: Miller, 2003, USGS WRIR 03-4107, Tables 1-6 (equations, SE_E, SE_P, 95-percent
prediction-interval factors), Table 7 (ranges) and Table 10 (gages by region). Every
literal number below is copied from that report, or from a dated live NSS response,
never from ``WY.json`` itself.
"""

from __future__ import annotations

import math

import pytest
from scipy import stats

from flowfreq.regression import EquationsUnavailable, OutOfRangeError, evaluate, load_state

AEPS = [0.667, 0.5, 0.429, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]

# In range for every region (Table 7 intersection).
BASIN = {
    "DRNAREA": 50.0,
    "ELEV": 8000.0,
    "LONG_OUT": 108.0,
    "SOILINDEX": 3.0,
    "MARAVPRE": 2.0,
    "LAT_OUT": 41.8,
    "JANAVPRE": 3.0,
}


@pytest.fixture(scope="module")
def lib():
    return load_state("WY")


def _sig3(x: float) -> float:
    return float(f"{x:.3g}")


def _pct(sep_log: float) -> float:
    return 100 * math.sqrt(10 ** (math.log(10) * sep_log**2) - 1)


class TestLibraryShape:
    def test_status_and_regions(self, lib):
        assert lib.status == "verified"
        assert lib.issue == 98
        assert lib.regions == ["1", "2", "3", "4", "5", "6"]
        for region in lib.regions:
            assert lib.aeps(region) == AEPS
        assert len(lib.equations) == 60

    def test_unpublished_aep_unavailable(self, lib):
        with pytest.raises(EquationsUnavailable):
            lib.equation("1", 0.25)

    @pytest.mark.parametrize(
        "region, n", [("1", 139), ("2", 74), ("3", 65), ("4", 20), ("5", 37), ("6", 29)]
    )
    def test_n_sites_from_table_10(self, lib, region, n):
        assert {e.n_sites for e in lib.equations if e.region_code == region} == {n}

    @pytest.mark.parametrize(
        "region, limits",
        [
            ("1", {"DRNAREA": (0.52, 2620), "ELEV": (5950, 10700), "LONG_OUT": (105.21, 111.34)}),
            ("2", {"DRNAREA": (0.06, 2070)}),
            ("3", {"DRNAREA": (0.20, 1230), "SOILINDEX": (2.1, 4.0)}),
            ("4", {"DRNAREA": (0.20, 471), "MARAVPRE": (0.71, 4.63), "LAT_OUT": (40.54, 44.57)}),
            ("5", {"DRNAREA": (2.77, 564), "JANAVPRE": (0.94, 8.77)}),
            ("6", {"DRNAREA": (1.26, 1180), "LAT_OUT": (41.02, 42.59)}),
        ],
    )
    def test_limits_from_table_7(self, lib, region, limits):
        eq = lib.equation(region, 0.01)
        assert {v.code: (v.minimum, v.maximum) for v in eq.variables} == limits

    def test_shifted_terms_stored_as_published(self, lib):
        v = {v.code: v for v in lib.equation("1", 0.01).variables}
        assert (v["ELEV"].scale, v["ELEV"].offset) == (0.001, -3.0)
        assert (v["LONG_OUT"].scale, v["LONG_OUT"].offset) == (1.0, -100.0)
        v = {v.code: v for v in lib.equation("6", 0.01).variables}
        assert v["LAT_OUT"].offset == -40.0

    def test_out_of_range_raises(self, lib):
        with pytest.raises(OutOfRangeError, match="LAT_OUT=43.0 > max 42.59"):
            evaluate(lib.equation("6", 0.01), {**BASIN, "LAT_OUT": 43.0})


class TestPrintedForms:
    @pytest.mark.parametrize(
        "region, aep, printed",
        [
            (
                "1",
                0.667,
                lambda b: 0.126
                * b["DRNAREA"] ** 0.885
                * ((b["ELEV"] - 3000) / 1000) ** 2.56
                * (b["LONG_OUT"] - 100) ** 0.032,
            ),
            (
                "1",
                0.01,
                lambda b: 38.6
                * b["DRNAREA"] ** 0.764
                * ((b["ELEV"] - 3000) / 1000) ** 1.00
                * (b["LONG_OUT"] - 100) ** -0.562,
            ),
            ("2", 0.002, lambda b: 728 * b["DRNAREA"] ** 0.425),
            ("3", 0.429, lambda b: 3.10 * b["DRNAREA"] ** 0.403 * b["SOILINDEX"] ** 2.84),
            (
                "4",
                0.002,
                lambda b: 52.5
                * b["DRNAREA"] ** 0.585
                * (b["LAT_OUT"] - 40) ** 0.766
                * b["MARAVPRE"] ** 0.873,
            ),
            ("5", 0.1, lambda b: 8.71 * b["DRNAREA"] ** 0.861 * b["JANAVPRE"] ** 0.529),
            ("6", 0.5, lambda b: 22.2 * b["DRNAREA"] ** 0.608 * (b["LAT_OUT"] - 40) ** -1.24),
        ],
    )
    def test_power_form(self, lib, region, aep, printed):
        assert evaluate(lib.equation(region, aep), BASIN).flow_cfs == pytest.approx(
            printed(BASIN), rel=1e-12
        )

    @pytest.mark.parametrize(
        "region, aep, see, sep",
        [("1", 0.667, 55, 56), ("2", 0.667, 131, 135), ("3", 0.04, 43, 46), ("4", 0.002, 41, 56)],
    )
    def test_error_statistics(self, lib, region, aep, see, sep):
        eq = lib.equation(region, aep)
        assert _pct(eq.sep_log) == pytest.approx(sep, abs=1e-3)
        assert _pct(math.sqrt(eq.model_error_variance)) == pytest.approx(see, abs=1e-3)
        assert eq.covariance is None and eq.avp is None

    @pytest.mark.parametrize(
        "region, aep, upper",
        [("1", 0.667, 2.82), ("2", 0.01, 2.83), ("4", 0.5, 2.87), ("6", 0.002, 3.83)],
    )
    def test_table_10_counts_reproduce_printed_95_percent_factors(self, lib, region, aep, upper):
        """The 95% factor is 10^(t(n - p) SE_P) (p. 19 note 2); SE_P is printed to 2-3
        significant figures, so agreement is to that rounding."""
        eq = lib.equation(region, aep)
        t = stats.t.ppf(0.975, eq.n_sites - len(eq.variables) - 1)
        assert 10 ** (t * eq.sep_log) == pytest.approx(upper, rel=0.01)


class TestWorkedExamples:
    """WRIR 03-4107 pp. 30-33."""

    def test_murphy_creek_region_3_q100(self, lib):
        est = evaluate(lib.equation("3", 0.01), {"DRNAREA": 49.0, "SOILINDEX": 3.5})
        assert round(est.flow_cfs, -2) == 8900
        assert round(est.flow_cfs * 0.382, -2) == 3400  # printed lower 95% limit

    def test_north_fork_crazy_woman_q25_area_weighted(self, lib):
        r1 = evaluate(
            lib.equation("1", 0.04), {"DRNAREA": 170.0, "ELEV": 7300.0, "LONG_OUT": 106.8}
        ).flow_cfs
        r3 = evaluate(lib.equation("3", 0.04), {"DRNAREA": 170.0, "SOILINDEX": 2.9}).flow_cfs
        assert round(r1, -1) == 2080
        assert round(r3, -1) == 4630
        assert round(r1 * 120 / 170 + r3 * 50 / 170, -2) == 2800

    def test_coney_creek_region_1_q50(self, lib):
        q = evaluate(
            lib.equation("1", 0.02), {"DRNAREA": 3.41, "ELEV": 9440.0, "LONG_OUT": 107.32}
        ).flow_cfs
        assert round(q) == 182

    def test_poison_creek_region_2_q100(self, lib):
        q = evaluate(lib.equation("2", 0.01), {"DRNAREA": 500.0}).flow_cfs
        assert round(q, -2) == 6000


# Live NSS responses for BASIN, fetched 2026-10-01 (3 significant figures).
NSS_RECORDED = [
    ("1", 0.667, 264),
    ("1", 0.01, 1190),
    ("1", 0.002, 1580),
    ("2", 0.667, 119),
    ("2", 0.01, 2230),
    ("3", 0.01, 6540),
    ("3", 0.002, 13800),
    ("4", 0.667, 67.1),
    ("4", 0.01, 840),
    ("5", 0.01, 697),
    ("5", 0.002, 857),
    ("6", 0.667, 73.5),
    ("6", 0.002, 1800),
]


@pytest.mark.parametrize("region, aep, nss", NSS_RECORDED)
def test_worked_example_matches_recorded_nss(lib, region, aep, nss):
    assert _sig3(evaluate(lib.equation(region, aep), BASIN).flow_cfs) == nss


NSS_CODE_TO_REGION = {
    "GC1239": "1",
    "GC1240": "2",
    "GC1241": "3",
    "GC1242": "4",
    "GC1243": "5",
    "GC1244": "6",
}
NSS_STAT_TO_AEP = {
    "PK66_7AEP": 0.667,
    "PK50AEP": 0.5,
    "PK42_9AEP": 0.429,
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
        {
            "DRNAREA": 5.0,
            "ELEV": 6500.0,
            "LONG_OUT": 106.0,
            "SOILINDEX": 2.5,
            "MARAVPRE": 1.0,
            "LAT_OUT": 41.2,
            "JANAVPRE": 1.5,
        },
        BASIN,
        {
            "DRNAREA": 400.0,
            "ELEV": 10000.0,
            "LONG_OUT": 110.5,
            "SOILINDEX": 3.8,
            "MARAVPRE": 4.0,
            "LAT_OUT": 42.5,
            "JANAVPRE": 7.0,
        },
    ],
)
def test_matches_live_nss(lib, basin):
    """Every region x AEP against NSS region 54 (Wyoming), citation 157, exact at NSS's
    3 significant figures."""
    from flowfreq.streamstats import Characteristic, estimate_flow_statistics

    chars = {k: Characteristic(k, k, "", v, "") for k, v in basin.items()}
    results, skipped = estimate_flow_statistics("WY", chars, ["PFS"])
    assert not skipped
    compared = 0
    for r in results:
        region = NSS_CODE_TO_REGION[r.region_code]
        for code, est in r.estimates.items():
            ours = evaluate(lib.equation(region, NSS_STAT_TO_AEP[code]), basin).flow_cfs
            assert _sig3(ours) == est.value, (region, code, est.value, ours)
            compared += 1
    assert compared == 60
