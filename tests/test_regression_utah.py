"""Utah peak-flow regression equations (issue #96).

Sources, as NSS combines them for Utah's peak-flow statistic group:

- Kenney, Wilkowske and Wright, 2007, USGS SIR 2007-5158 ver. 4.0 (March 10, 2008),
  Table 5 (equations and errors) and Table 8 (ranges): regions 1-7.
- Waltemeyer, 2006, USGS SIR 2006-5306, Tables 1-2: Navajo Nation region 8.

Every literal number below is copied from those reports, or from a dated live NSS
response, never from ``UT.json`` itself. The report has no worked example, so the
offline worked example is NSS's response for one in-range basin.
"""

from __future__ import annotations

import math

import pytest

from flowfreq.regression import EquationsUnavailable, OutOfRangeError, evaluate, load_state

AEPS_8 = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
AEPS_NAVAJO = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002]

# In range for regions 1-7 (Table 8 intersection); region 1 and 4 need ELEV >= 8,130.
BASIN = {
    "DRNAREA": 30.0,
    "ELEV": 8500.0,
    "PRECIP": 35.0,
    "BSLDEM10M": 25.0,
    "LU92HRBN": 10.0,
    "BSLDEM30ff": 0.2,
}
# In range for Navajo region 8 (ELEV <= 7,750) and regions 2, 3, 5, 6, 7.
LOW_BASIN = {
    "DRNAREA": 50.0,
    "ELEV": 5000.0,
    "PRECIP": 25.0,
    "BSLDEM10M": 20.0,
    "LU92HRBN": 3.0,
    "BSLDEM30ff": 0.25,
}

# Table 5 region 2 bases b in b^PRECIP, as printed, and as NSS evaluates them
# (NSS equation strings, snapshot 2026-10-01).
REGION2_PRINTED = [1.07, 1.07, 1.06, 1.06, 1.06, 1.06, 1.055, 1.05]
REGION2_NSS = [1.074, 1.067, 1.064, 1.059, 1.059, 1.057, 1.055, 1.054]


@pytest.fixture(scope="module")
def lib():
    return load_state("UT")


def _sig3(x: float) -> float:
    return float(f"{x:.3g}")


def _pct(sep_log: float) -> float:
    return 100 * math.sqrt(10 ** (math.log(10) * sep_log**2) - 1)


class TestLibraryShape:
    def test_status_and_regions(self, lib):
        assert lib.status == "partial"
        assert lib.issue == 96
        assert lib.regions == ["1", "2", "3", "4", "5", "6", "7", "Navajo8"]
        for region in "1234567":
            assert lib.aeps(region) == AEPS_8
        assert lib.aeps("Navajo8") == AEPS_NAVAJO
        assert len(lib.equations) == 63

    def test_unpublished_aep_unavailable(self, lib):
        with pytest.raises(EquationsUnavailable):
            lib.equation("Navajo8", 0.005)

    @pytest.mark.parametrize(
        "region, n",
        [
            ("1", 46),
            ("2", 32),
            ("3", 14),
            ("4", 42),
            ("5", 35),
            ("6", 99),
            ("7", 25),
            ("Navajo8", 62),
        ],
    )
    def test_n_sites(self, lib, region, n):
        assert {e.n_sites for e in lib.equations if e.region_code == region} == {n}

    @pytest.mark.parametrize(
        "region, limits",
        [
            ("1", {"DRNAREA": (3.62, 390), "ELEV": (6420, 10500)}),
            ("2", {"DRNAREA": (2.14, 84.1), "PRECIP": (16.5, 53.7)}),
            ("3", {"DRNAREA": (5.72, 66.5)}),
            ("4", {"DRNAREA": (2.95, 667), "ELEV": (8130, 10900), "BSLDEM10M": (9.67, 40.3)}),
            ("5", {"DRNAREA": (0.91, 629), "LU92HRBN": (2.14, 15.6)}),
            ("6", {"DRNAREA": (0.87, 532), "ELEV": (4300, 9380)}),
            ("7", {"DRNAREA": (5.43, 1670)}),
        ],
    )
    def test_limits_from_table_8(self, lib, region, limits):
        eq = lib.equation(region, 0.01)
        assert {v.code: (v.minimum, v.maximum) for v in eq.variables} == limits

    def test_out_of_range_raises(self, lib):
        with pytest.raises(OutOfRangeError, match="ELEV=7000.0 < min 8130"):
            evaluate(lib.equation("4", 0.01), {**BASIN, "ELEV": 7000.0})


class TestPrintedForms:
    """Table 5's printed power and exponential forms, evaluated directly."""

    @pytest.mark.parametrize(
        "region, aep, printed",
        [
            ("1", 0.5, lambda b: 1.52 * b["DRNAREA"] ** 0.677 * 1.39 ** (b["ELEV"] / 1000)),
            ("1", 0.002, lambda b: 85.4 * b["DRNAREA"] ** 0.457 * 1.13 ** (b["ELEV"] / 1000)),
            ("2", 0.5, lambda b: 0.585 * b["DRNAREA"] ** 0.847 * 1.07 ** b["PRECIP"]),
            ("2", 0.005, lambda b: 8.79 * b["DRNAREA"] ** 0.592 * 1.055 ** b["PRECIP"]),
            ("3", 0.01, lambda b: 300 * b["DRNAREA"] ** 0.303),
            (
                "4",
                0.02,
                lambda b: 2.68
                * b["DRNAREA"] ** 0.798
                * 2.72 ** (0.373 * b["ELEV"] / 1000 - 0.028 * b["BSLDEM10M"]),
            ),
            (
                "4",
                0.5,
                lambda b: 0.083
                * b["DRNAREA"] ** 0.822
                * 2.72 ** (0.656 * b["ELEV"] / 1000 - 0.039 * b["BSLDEM10M"]),
            ),
            ("5", 0.1, lambda b: 18.4 * b["DRNAREA"] ** 0.555 * (b["LU92HRBN"] + 1) ** 0.388),
            ("6", 0.01, lambda b: 115000 * b["DRNAREA"] ** 0.391 * (b["ELEV"] / 1000) ** -2.58),
            ("7", 0.002, lambda b: 1620 * b["DRNAREA"] ** 0.280),
        ],
    )
    def test_power_form(self, lib, region, aep, printed):
        est = evaluate(lib.equation(region, aep), BASIN)
        assert est.flow_cfs == pytest.approx(printed(BASIN), rel=1e-12)

    def test_navajo_region_8(self, lib):
        b = LOW_BASIN
        est = evaluate(lib.equation("Navajo8", 0.01), b)
        expect = 1.53e12 * b["DRNAREA"] ** 0.265 * b["BSLDEM30ff"] ** 0.477 * b["ELEV"] ** -2.25
        assert est.flow_cfs == pytest.approx(expect, rel=1e-12)

    def test_exponential_terms_are_identity_transforms(self, lib):
        v = {v.code: v for v in lib.equation("1", 0.01).variables}["ELEV"]
        assert (v.transform, v.scale) == ("identity", 0.001)
        v = {v.code: v for v in lib.equation("4", 0.01).variables}
        assert (v["ELEV"].transform, v["ELEV"].scale) == ("identity", 0.001)
        assert v["BSLDEM10M"].transform == "identity"
        v = {v.code: v for v in lib.equation("5", 0.01).variables}["LU92HRBN"]
        assert (v.transform, v.offset) == ("log10", 1.0)
        v = {v.code: v for v in lib.equation("6", 0.01).variables}["ELEV"]
        assert (v.transform, v.scale) == ("log10", 0.001)

    @pytest.mark.parametrize(
        "region, aep, sep_pct, me_pct",
        [
            ("1", 0.5, 62, 59),
            ("2", 0.005, 51, 47),
            ("3", 0.5, 357, 295),
            ("4", 0.1, 35, 32),
            ("5", 0.002, 86, 76),
            ("6", 0.01, 61, 58),
            ("7", 0.002, 211, 189),
        ],
    )
    def test_error_statistics(self, lib, region, aep, sep_pct, me_pct):
        eq = lib.equation(region, aep)
        assert _pct(eq.sep_log) == pytest.approx(sep_pct, abs=1e-3)
        assert _pct(math.sqrt(eq.model_error_variance)) == pytest.approx(me_pct, abs=1e-3)
        assert eq.avp is None and eq.covariance is None

    def test_interval_uses_sep(self, lib):
        assert evaluate(lib.equation("6", 0.01), BASIN).interval_method == (
            "average standard error of prediction"
        )


class TestRegion2Rounding:
    """Table 5 prints region 2's bases to 3 significant figures; NSS uses 4."""

    def test_stored_bases_are_the_printed_ones(self, lib):
        for aep, b in zip(AEPS_8, REGION2_PRINTED, strict=True):
            eq = lib.equation("2", aep)
            assert eq.coefficients[1] == pytest.approx(math.log10(b), rel=1e-12)

    def test_only_the_unrounded_pk200_base_agrees_with_nss(self):
        assert [p == n for p, n in zip(REGION2_PRINTED, REGION2_NSS)] == [
            False,
            False,
            False,
            False,
            False,
            False,
            True,
            False,
        ]

    def test_size_of_the_discrepancy(self, lib):
        # At PRECIP = 35 in., the printed 2-year equation is 12% below NSS's.
        est = evaluate(lib.equation("2", 0.5), BASIN)
        nss_form = 0.585 * 30**0.847 * 1.074**35
        assert est.flow_cfs / nss_form == pytest.approx((1.07 / 1.074) ** 35, rel=1e-12)
        assert est.flow_cfs / nss_form == pytest.approx(0.878, abs=1e-3)


# Live NSS responses, fetched 2026-10-01 (NSS rounds to 3 significant figures).
NSS_RECORDED = [
    (BASIN, "1", 0.5, 250),
    (BASIN, "1", 0.01, 875),
    (BASIN, "1", 0.002, 1140),
    (BASIN, "2", 0.005, 429),
    (BASIN, "3", 0.01, 841),
    (BASIN, "3", 0.002, 1660),
    (BASIN, "4", 0.5, 136),
    (BASIN, "4", 0.01, 596),
    (BASIN, "5", 0.01, 617),
    (BASIN, "5", 0.002, 923),
    (BASIN, "6", 0.01, 1740),
    (BASIN, "7", 0.01, 2280),
    (LOW_BASIN, "Navajo8", 0.5, 655),
    (LOW_BASIN, "Navajo8", 0.01, 10600),
    (LOW_BASIN, "Navajo8", 0.002, 18400),
]
# NSS's region 2 values at BASIN, which the printed equations do not reproduce.
NSS_REGION2_AT_BASIN = [127, 192, 240, 282, 346, 387, 429, 499]


@pytest.mark.parametrize("basin, region, aep, nss", NSS_RECORDED)
def test_worked_example_matches_recorded_nss(lib, basin, region, aep, nss):
    assert _sig3(evaluate(lib.equation(region, aep), basin).flow_cfs) == nss


def test_recorded_nss_region_2_is_the_unrounded_bases(lib):
    for aep, nss, p, n in zip(AEPS_8, NSS_REGION2_AT_BASIN, REGION2_PRINTED, REGION2_NSS):
        ours = evaluate(lib.equation("2", aep), BASIN).flow_cfs
        assert ours * (n / p) ** BASIN["PRECIP"] == pytest.approx(nss, rel=5e-3), aep


NSS_CODE_TO_REGION = {
    "GC959": "Navajo8",
    "GC960": "1",
    "GC961": "2",
    "GC962": "3",
    "GC963": "4",
    "GC964": "5",
    "GC965": "6",
    "GC966": "7",
}
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
    "basin, n_regions",
    [
        (
            {
                "DRNAREA": 10.0,
                "ELEV": 7000.0,
                "PRECIP": 20.0,
                "BSLDEM10M": 15.0,
                "LU92HRBN": 5.0,
                "BSLDEM30ff": 0.1,
            },
            7,
        ),
        (BASIN, 7),
        (
            {
                "DRNAREA": 60.0,
                "ELEV": 9000.0,
                "PRECIP": 50.0,
                "BSLDEM10M": 35.0,
                "LU92HRBN": 15.0,
                "BSLDEM30ff": 0.3,
            },
            7,
        ),
        (LOW_BASIN, 6),
    ],
)
def test_matches_live_nss(lib, basin, n_regions):
    """Every in-range region x AEP against NSS region 48 (Utah), citations 140 and 11.

    Exact agreement at NSS's 3 significant figures, except region 2 at the seven AEPs
    whose printed base is rounded: there NSS equals the printed equation rescaled by
    (NSS base / printed base)^PRECIP.
    """
    from flowfreq.streamstats import Characteristic, estimate_flow_statistics

    chars = {k: Characteristic(k, k, "", v, "") for k, v in basin.items()}
    results, skipped = estimate_flow_statistics("UT", chars, ["PFS"])
    assert len(results) == n_regions, skipped
    for r in results:
        region = NSS_CODE_TO_REGION[r.region_code]
        for code, est in r.estimates.items():
            aep = NSS_STAT_TO_AEP[code]
            ours = evaluate(lib.equation(region, aep), basin).flow_cfs
            if region == "2" and aep != 0.005:
                i = AEPS_8.index(aep)
                ours *= (REGION2_NSS[i] / REGION2_PRINTED[i]) ** basin["PRECIP"]
            assert _sig3(ours) == est.value, (region, code, est.value, ours)
