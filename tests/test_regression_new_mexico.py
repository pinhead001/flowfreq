"""New Mexico peak-flow regression equations (issue #100).

Sources, as NSS uses them for New Mexico's peak-flow statistic group:

- Waltemeyer, 2008, USGS SIR 2008-5119, Table 1 (ranges), Table 2 (equations and errors)
  and Appendix 1 (gages): flood regions 1-9.
- Waltemeyer, 2006, USGS SIR 2006-5306 (Navajo Nation): regions 8, 11, High Elevation, 6.

Every literal number below is copied from those reports, or from a dated live NSS
response, never from ``NM.json`` itself.
"""

from __future__ import annotations

import math

import pytest

from flowfreq.regression import EquationsUnavailable, OutOfRangeError, evaluate, load_state

AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002]
REGIONS = ["1", "2", "3", "4", "5", "6", "7", "8", "9"]

# In range for regions 1-4 and 6-9 and every Navajo region (region 5 needs ELEV >= 7,570).
BASIN = {
    "DRNAREA": 50.0,
    "BSLDEM30ff": 0.1,
    "PRECIP": 17.0,
    "ELEV": 7000.0,
    "I24H100YA2": 4.0,
}
# In range for region 5 (and 1, 3, 4, 6-9).
HIGH_BASIN = {
    "DRNAREA": 400.0,
    "BSLDEM30ff": 0.15,
    "PRECIP": 20.0,
    "ELEV": 7580.0,
    "I24H100YA2": 5.0,
}


@pytest.fixture(scope="module")
def lib():
    return load_state("NM")


def _sig3(x: float) -> float:
    return float(f"{x:.3g}")


class TestLibraryShape:
    def test_status_and_regions(self, lib):
        assert lib.status == "verified"
        assert lib.issue == 100
        assert lib.regions == REGIONS + ["Navajo11", "Navajo6", "Navajo8", "NavajoHighElev"]
        for region in lib.regions:
            assert lib.aeps(region) == AEPS
        assert len(lib.equations) == 91

    def test_no_200_year_equation(self, lib):
        with pytest.raises(EquationsUnavailable):
            lib.equation("1", 0.005)

    @pytest.mark.parametrize(
        "region, n",
        [
            ("1", 29),
            ("2", 62),
            ("3", 12),
            ("4", 29),
            ("5", 64),
            ("6", 23),
            ("7", 14),
            ("8", 28),
            ("9", 32),
        ],
    )
    def test_n_sites_from_appendix_1(self, lib, region, n):
        assert {e.n_sites for e in lib.equations if e.region_code == region} == {n}

    @pytest.mark.parametrize(
        "region, limits",
        [
            (
                "1",
                {
                    "DRNAREA": (0.536, 10800),
                    "BSLDEM30ff": (0.0205, 0.2023),
                    "PRECIP": (14.50, 20.15),
                },
            ),
            # Table 1 prints "12,7000"; Appendix 1's largest region 2 gage is 12,700 mi2.
            (
                "2",
                {"DRNAREA": (0.059, 12700), "BSLDEM30ff": (0.0231, 0.4447), "ELEV": (4320, 7520)},
            ),
            ("3", {"DRNAREA": (3.35, 716), "ELEV": (4890, 8820)}),
            ("4", {"DRNAREA": (0.140, 4460), "ELEV": (3970, 7620), "BSLDEM30ff": (0.0194, 0.1966)}),
            ("5", {"DRNAREA": (0.606, 3120), "ELEV": (7570, 11500), "I24H100YA2": (3.08, 5.51)}),
            ("6", {"DRNAREA": (0.18, 7220)}),
            ("7", {"DRNAREA": (0.20, 2820)}),
            ("8", {"DRNAREA": (1.70, 1200)}),
            ("9", {"DRNAREA": (0.107, 2150), "BSLDEM30ff": (0.0177, 0.1545)}),
        ],
    )
    def test_limits_from_table_1(self, lib, region, limits):
        eq = lib.equation(region, 0.01)
        assert {v.code: (v.minimum, v.maximum) for v in eq.variables} == limits

    def test_elevation_is_scaled(self, lib):
        v = {v.code: v for v in lib.equation("3", 0.01).variables}["ELEV"]
        assert (v.transform, v.scale) == ("log10", 0.001)

    def test_out_of_range_raises(self, lib):
        with pytest.raises(OutOfRangeError, match="PRECIP=21.0 > max 20.15"):
            evaluate(lib.equation("1", 0.01), {**BASIN, "PRECIP": 21.0})


class TestPrintedForms:
    @pytest.mark.parametrize(
        "region, aep, printed",
        [
            (
                "1",
                0.5,
                lambda b: 8.260e6
                * b["DRNAREA"] ** 0.562
                * b["BSLDEM30ff"] ** 0.351
                * b["PRECIP"] ** -3.81,
            ),
            ("2", 0.5, lambda b: 9.656e1 * b["DRNAREA"] ** 0.428),
            (
                "2",
                0.01,
                lambda b: 1.531e5
                * b["DRNAREA"] ** 0.295
                * b["BSLDEM30ff"] ** 0.457
                * (b["ELEV"] / 1000) ** -1.99,
            ),
            ("3", 0.002, lambda b: 3.430e6 * b["DRNAREA"] ** 0.876 * (b["ELEV"] / 1000) ** -4.71),
            (
                "4",
                0.1,
                lambda b: 9.587e4
                * b["DRNAREA"] ** 0.460
                * (b["ELEV"] / 1000) ** -2.31
                * b["BSLDEM30ff"] ** 0.494,
            ),
            ("6", 0.002, lambda b: 2.354e3 * b["DRNAREA"] ** 0.348),
            ("7", 0.04, lambda b: 5.537e2 * b["DRNAREA"] ** 0.488),
            ("8", 0.5, lambda b: 8.990 * b["DRNAREA"] ** 0.800),
            ("9", 0.002, lambda b: 1.335e4 * b["DRNAREA"] ** 0.446 * b["BSLDEM30ff"] ** 0.904),
        ],
    )
    def test_power_form(self, lib, region, aep, printed):
        assert evaluate(lib.equation(region, aep), BASIN).flow_cfs == pytest.approx(
            printed(BASIN), rel=1e-12
        )

    def test_region_5(self, lib):
        b = HIGH_BASIN
        expect = (
            5.211e4 * b["DRNAREA"] ** 0.642 * (b["ELEV"] / 1000) ** -5.36 * b["I24H100YA2"] ** 4.55
        )
        assert evaluate(lib.equation("5", 0.002), b).flow_cfs == pytest.approx(expect, rel=1e-12)

    @pytest.mark.parametrize(
        "region, aep, reg_log, pred_log",
        [
            ("1", 0.5, 0.298, 0.320),
            ("2", 0.5, 0.419, 0.427),
            ("5", 0.002, 0.304, 0.318),
            ("9", 0.002, 0.345, 0.373),
        ],
    )
    def test_error_statistics_in_log_units(self, lib, region, aep, reg_log, pred_log):
        eq = lib.equation(region, aep)
        assert eq.sep_log == pred_log
        assert eq.model_error_variance == pytest.approx(reg_log**2, rel=1e-12)
        assert eq.covariance is None


class TestWorkedExamples:
    def test_region_2_q100_example(self, lib):
        """p. 11: A 523, S 0.2422, E 6,940 -> printed 10,800 ft3/s (the equation gives
        10,745; the example rounds its intermediate factors)."""
        q = evaluate(
            lib.equation("2", 0.01), {"DRNAREA": 523.0, "BSLDEM30ff": 0.2422, "ELEV": 6940.0}
        ).flow_cfs
        assert q == pytest.approx(10800, rel=0.006)

    @pytest.mark.parametrize(
        "region, site, predicted",
        [
            # Appendix 1 'predicted' lines: map 1 (07153500) and map 227 (08478800, whose
            # predictions use Appendix 2's 0.20 mi2, not Appendix 1's GIS area of 1.80).
            (
                "1",
                {"DRNAREA": 611.0, "BSLDEM30ff": 0.1132, "PRECIP": 16.85},
                {0.5: 3020, 0.2: 7210, 0.01: 38200, 0.002: 73300},
            ),
            ("7", {"DRNAREA": 0.20}, {0.5: 71, 0.2: 131, 0.01: 377, 0.002: 544}),
        ],
    )
    def test_appendix_1_predicted_values(self, lib, region, site, predicted):
        for aep, q in predicted.items():
            assert evaluate(lib.equation(region, aep), site).flow_cfs == pytest.approx(q, rel=0.01)


# Live NSS responses, fetched 2026-10-02 (3 significant figures).
NSS_RECORDED = [
    ("1", 0.5, 680),
    ("1", 0.01, 8680),
    ("2", 0.01, 3530),
    ("3", 0.01, 6020),
    ("4", 0.5, 542),
    ("6", 0.01, 5600),
    ("7", 0.5, 865),
    ("8", 0.01, 1630),
    ("9", 0.01, 5400),
    ("Navajo8", 0.01, 3210),
    ("Navajo11", 0.01, 5770),
    ("NavajoHighElev", 0.5, 219),
]


@pytest.mark.parametrize("region, aep, nss", NSS_RECORDED)
def test_worked_example_matches_recorded_nss(lib, region, aep, nss):
    assert _sig3(evaluate(lib.equation(region, aep), BASIN).flow_cfs) == nss


NSS_CODE_TO_REGION = {f"GC{1091 + i}": str(i) for i in range(1, 10)}
NSS_CODE_TO_REGION.update(
    {"GC1088": "Navajo8", "GC1089": "NavajoHighElev", "GC1090": "Navajo11", "GC1091": "Navajo6"}
)
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
@pytest.mark.parametrize("basin", [BASIN, HIGH_BASIN])
@pytest.mark.parametrize("region", REGIONS)
def test_matches_live_nss(lib, basin, region):
    """SIR 2008-5119 regions are gated in NSS (region 35) by HIGHREG = the region's GC
    number. Every in-range region returned must equal flowfreq at 3 significant figures."""
    from flowfreq.streamstats import Characteristic, estimate_flow_statistics

    gc = f"GC{1091 + int(region)}"
    b = {**basin, "HIGHREG": float(gc[2:])}
    chars = {k: Characteristic(k, k, "", v, "") for k, v in b.items()}
    results, _ = estimate_flow_statistics("NM", chars, ["PFS"])
    for r in results:
        if r.region_code not in NSS_CODE_TO_REGION:
            continue  # WSP 2433 regions, not stored
        mapped = NSS_CODE_TO_REGION[r.region_code]
        assert not (mapped in REGIONS and mapped != region), "HIGHREG gating leaked"
        for code, est in r.estimates.items():
            ours = evaluate(lib.equation(mapped, NSS_STAT_TO_AEP[code]), basin).flow_cfs
            assert _sig3(ours) == est.value, (mapped, code, est.value, ours)
