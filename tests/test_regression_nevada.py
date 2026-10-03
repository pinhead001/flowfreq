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
