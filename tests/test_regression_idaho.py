"""Idaho peak-flow regression equations (issue #39).

Source: Wood and others, 2016, USGS SIR 2016-5083 ver. 1.1 (April 2017), Table 4
(equations, MEV, AVP, SEP, n), Table 5 (variable ranges) and Table A5 (covariance).
Every literal number below is copied from that report, or from a dated live NSS
response, never from ``ID.json`` itself: the tests check the file against the
source, not against itself.

The report gives no worked example, so the offline "worked example" is the value
NSS returned for one in-range basin on 2026-09-25 (region 4: 2026-09-26; NSS rounds
to three significant figures), plus direct evaluation of the report's printed
power-form equations.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from flowfreq.regression import EquationsUnavailable, OutOfRangeError, evaluate, load_state

AEPS = [0.8, 0.667, 0.5, 0.429, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]

# An in-range basin for every region (Table 5 intersection).
BASIN = {"DRNAREA": 100.0, "PRECPRIS10": 35.0, "LC11FOREST": 50.0, "MINBELEV": 4000.0}


@pytest.fixture(scope="module")
def lib():
    return load_state("ID")


def _sig3(x: float) -> float:
    return float(f"{x:.3g}")


class TestLibraryShape:
    def test_status_and_regions(self, lib):
        assert lib.status == "verified"
        assert lib.issue == 39
        assert lib.regions == ["1_2", "3", "4", "5", "6_8", "7"]
        for region in lib.regions:
            assert lib.aeps(region) == AEPS
        assert len(lib.equations) == 66

    def test_unpublished_aep_unavailable(self, lib):
        with pytest.raises(EquationsUnavailable):
            lib.equation("4", 0.3)

    def test_every_equation_is_cited_and_complete(self, lib):
        for eq in lib.equations:
            assert "2016-5083" in eq.citation.publication
            assert "Table 4" in eq.citation.table and "Table A5" in eq.citation.table
            assert eq.model_error_variance is not None and eq.avp is not None
            assert eq.sep_log is not None and eq.covariance is not None
            u = np.asarray(eq.covariance)
            assert np.allclose(u, u.T)
            assert np.all(np.linalg.eigvalsh(u) > 0)

    @pytest.mark.parametrize(
        "region, n", [("1_2", 59), ("3", 13), ("4", 28), ("5", 20), ("6_8", 48), ("7", 24)]
    )
    def test_n_sites(self, lib, region, n):
        assert {e.n_sites for e in lib.equations if e.region_code == region} == {n}

    def test_limits_from_table_5(self, lib):
        limits = {v.code: (v.minimum, v.maximum) for v in lib.equation("7", 0.01).variables}
        assert limits == {"DRNAREA": (0.15, 1400.0), "MINBELEV": (2230.0, 7020.0)}
        limits = {v.code: (v.minimum, v.maximum) for v in lib.equation("1_2", 0.01).variables}
        assert limits == {"DRNAREA": (1.11, 10700.0), "PRECPRIS10": (21.3, 66.8)}
        limits = {v.code: (v.minimum, v.maximum) for v in lib.equation("4", 0.01).variables}
        assert limits == {
            "DRNAREA": (5.36, 13400.0),
            "LC11FOREST": (2.35, 94.4),
            "PRECPRIS10": (18.7, 46.8),
        }

    def test_scaled_terms_stored_as_published(self, lib):
        # Region 4's (F/100 + 1)^b and region 7's (Emin/1,000)^b, with limits in
        # the published units (percent and feet).
        for aep in AEPS:
            v4 = {v.code: v for v in lib.equation("4", aep).variables}
            assert [v.code for v in lib.equation("4", aep).variables] == [
                "DRNAREA",
                "LC11FOREST",
                "PRECPRIS10",
            ]
            f = v4["LC11FOREST"]
            assert (f.transform, f.scale, f.offset) == ("log10", 0.01, 1.0)
            assert all(v.scale == 1.0 and v.offset == 0.0 for c, v in v4.items() if c != f.code)
            e = {v.code: v for v in lib.equation("7", aep).variables}["MINBELEV"]
            assert (e.transform, e.scale, e.offset) == ("log10", 0.001, 0.0)

    def test_out_of_range_raises(self, lib):
        with pytest.raises(OutOfRangeError, match="MINBELEV"):
            evaluate(lib.equation("7", 0.01), {"DRNAREA": 100.0, "MINBELEV": 8000.0})
        # The forest limit is checked in percent, not in the scaled F/100 + 1 space.
        with pytest.raises(OutOfRangeError, match="LC11FOREST=95.0 > max 94.4"):
            evaluate(lib.equation("4", 0.01), {**BASIN, "LC11FOREST": 95.0})


class TestAgainstTable4:
    """Evaluate the report's printed power-form equations directly."""

    @pytest.mark.parametrize(
        "region, aep, printed",
        [
            ("1_2", 0.01, lambda b: 0.0239 * b["DRNAREA"] ** 0.844 * b["PRECPRIS10"] ** 2.11),
            ("1_2", 0.667, lambda b: 0.00141 * b["DRNAREA"] ** 0.945 * b["PRECPRIS10"] ** 2.41),
            ("3", 0.002, lambda b: 171 * b["DRNAREA"] ** 0.764),
            ("5", 0.429, lambda b: 0.000984 * b["DRNAREA"] ** 0.944 * b["PRECPRIS10"] ** 2.63),
            ("6_8", 0.1, lambda b: 0.0534 * b["DRNAREA"] ** 0.713 * b["PRECPRIS10"] ** 1.91),
            ("7", 0.8, lambda b: 1.64 * b["DRNAREA"] ** 0.563 * (b["MINBELEV"] / 1000) ** 1.09),
            ("7", 0.01, lambda b: 4000 * b["DRNAREA"] ** 0.478 * (b["MINBELEV"] / 1000) ** -2.13),
            (
                "4",
                0.8,
                lambda b: 0.000592
                * b["DRNAREA"] ** 0.981
                * (b["LC11FOREST"] / 100 + 1) ** -1.52
                * b["PRECPRIS10"] ** 2.80,
            ),
            (
                "4",
                0.01,
                lambda b: 0.292
                * b["DRNAREA"] ** 0.884
                * (b["LC11FOREST"] / 100 + 1) ** -2.92
                * b["PRECPRIS10"] ** 1.80,
            ),
            (
                "4",
                0.002,
                lambda b: 1.09
                * b["DRNAREA"] ** 0.859
                * (b["LC11FOREST"] / 100 + 1) ** -3.11
                * b["PRECPRIS10"] ** 1.56,
            ),
        ],
    )
    def test_power_form(self, lib, region, aep, printed):
        est = evaluate(lib.equation(region, aep), BASIN)
        assert est.flow_cfs == pytest.approx(printed(BASIN), rel=1e-12)

    @pytest.mark.parametrize(
        "region, aep, mev, avp, sep_pct",
        [
            ("1_2", 0.01, 0.086, 0.093, 79.7),
            ("3", 0.8, 0.135, 0.158, 114),
            ("4", 0.5, 0.044, 0.051, 55.8),
            ("4", 0.002, 0.019, 0.026, 38.2),
            ("5", 0.002, 0.021, 0.026, 38.6),
            ("6_8", 0.5, 0.056, 0.060, 61.3),
            ("7", 0.01, 0.100, 0.121, 94.9),
        ],
    )
    def test_error_statistics(self, lib, region, aep, mev, avp, sep_pct):
        eq = lib.equation(region, aep)
        assert eq.model_error_variance == mev
        assert eq.avp == avp
        # sep_log is the percent SEP converted to log10 units
        assert 100 * math.sqrt(10 ** (math.log(10) * eq.sep_log**2) - 1) == pytest.approx(
            sep_pct, abs=1e-3
        )
        # and it agrees with sqrt(AVP) to the rounding of the printed columns
        assert eq.sep_log == pytest.approx(math.sqrt(avp), rel=0.02)


class TestCovariance:
    def test_region_3_q1_from_table_a5(self, lib):
        eq = lib.equation("3", 0.01)
        assert eq.covariance == (
            (2.0571213000000001e-2, -6.4550781999999996e-3),
            (-6.4550781999999996e-3, 2.6924032e-3),
        )

    def test_region_4_q1_from_table_a5(self, lib):
        """Table A5 region 4 is in the basis (1, log10 A, log10(F/100 + 1), log10 P)."""
        eq = lib.equation("4", 0.01)
        assert eq.covariance == (
            (0.30003542999999999, -1.0941779e-2, 0.14102358000000001, -0.19876046),
            (-1.0941779e-2, 1.7454025e-3, -4.9844045999999998e-3, 5.1330239999999999e-3),
            (
                0.14102358000000001,
                -4.9844045999999998e-3,
                0.36696785999999998,
                -0.13625475000000001,
            ),
            (-0.19876046, 5.1330239999999999e-3, -0.13625475000000001, 0.14174798),
        )

    def test_region_4_interval_uses_scaled_forest_term(self, lib):
        eq = lib.equation("4", 0.01)
        est = evaluate(eq, BASIN)
        x = np.array(
            [
                1.0,
                math.log10(BASIN["DRNAREA"]),
                math.log10(BASIN["LC11FOREST"] / 100 + 1),
                math.log10(BASIN["PRECPRIS10"]),
            ]
        )
        sep_i = math.sqrt(0.019 + x @ np.array(eq.covariance) @ x)  # MEV from Table 4
        ci = stats.t.ppf(0.95, 28 - (3 + 1)) * sep_i
        assert est.interval[1] == pytest.approx(10 ** (est.log_flow + ci), rel=1e-10)

    def test_region_7_covariance_stored_as_published(self, lib):
        """Table A5 region 7 is in the log10(MINBELEV/1000) basis, and ID.json keeps it."""
        report_u = (
            (0.30811845999999998, -2.7867076000000001e-2, -0.37168530999999999),
            (-2.7867076000000001e-2, 6.6210102999999998e-3, 2.3593382e-2),
            (-0.37168530999999999, 2.3593382e-2, 0.48953417999999999),
        )
        eq = lib.equation("7", 0.01)
        assert eq.covariance == report_u
        assert eq.intercept == pytest.approx(math.log10(4000), rel=1e-15)  # Table 4: Q1 = 4,000 ...

    def test_site_specific_interval_follows_report_equations_6_to_8(self, lib):
        eq = lib.equation("3", 0.01)
        est = evaluate(eq, BASIN)
        x = np.array([1.0, math.log10(BASIN["DRNAREA"])])
        u = np.array(eq.covariance)
        sep_i = math.sqrt(0.014 + x @ u @ x)  # MEV from Table 4
        ci = stats.t.ppf(0.95, 13 - (1 + 1)) * sep_i
        assert est.interval[1] == pytest.approx(10 ** (est.log_flow + ci), rel=1e-10)
        assert est.interval[0] == pytest.approx(10 ** (est.log_flow - ci), rel=1e-10)


# Live NSS responses for BASIN, fetched 2026-09-25 and, for region 4, 2026-09-26
# (NSS rounds to 3 significant figures).
NSS_RECORDED = [
    ("1_2", 0.8, 474),
    ("1_2", 0.01, 2110),
    ("1_2", 0.002, 2650),
    ("3", 0.8, 352),
    ("3", 0.01, 3890),
    ("4", 0.8, 617),
    ("4", 0.429, 1050),
    ("4", 0.01, 3150),
    ("4", 0.002, 4140),
    ("5", 0.5, 818),
    ("5", 0.01, 2190),
    ("6_8", 0.01, 1970),
    ("6_8", 0.002, 2470),
    ("7", 0.8, 99.3),
    ("7", 0.01, 1890),
    ("7", 0.002, 2990),
]


@pytest.mark.parametrize("region, aep, nss", NSS_RECORDED)
def test_worked_example_matches_recorded_nss(lib, region, aep, nss):
    assert _sig3(evaluate(lib.equation(region, aep), BASIN).flow_cfs) == nss


NSS_CODE_TO_REGION = {
    "GC1735": "1_2",
    "GC1736": "3",
    "GC1737": "4",
    "GC1738": "5",
    "GC1739": "6_8",
    "GC1740": "7",
}
NSS_STAT_TO_AEP = {
    "PK80AEP": 0.8,
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
        {"DRNAREA": 10.0, "PRECPRIS10": 31.0, "LC11FOREST": 20.0, "MINBELEV": 2500.0},
        BASIN,
        {"DRNAREA": 875.0, "PRECPRIS10": 45.5, "LC11FOREST": 80.0, "MINBELEV": 6800.0},
    ],
)
def test_matches_live_nss(lib, basin):
    """Every region x AEP against NSS region 16 (Idaho), citation 39.

    NSS returns three significant figures, so the check is exact agreement after
    rounding flowfreq's value the same way; the relative difference is bounded by
    that rounding (at most 0.5% for a leading digit of 1).
    """
    from flowfreq.streamstats import Characteristic, estimate_flow_statistics

    chars = {k: Characteristic(k, k, "", v, "") for k, v in basin.items()}
    results, skipped = estimate_flow_statistics("ID", chars, ["PFS"])
    assert not skipped
    compared = 0
    for r in results:
        region = NSS_CODE_TO_REGION[r.region_code]
        for code, est in r.estimates.items():
            ours = evaluate(lib.equation(region, NSS_STAT_TO_AEP[code]), basin).flow_cfs
            assert _sig3(ours) == est.value, (region, code, est.value, ours)
            assert abs(ours - est.value) / est.value < 5e-3
            compared += 1
    assert compared == 66
