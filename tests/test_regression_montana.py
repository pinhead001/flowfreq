"""Montana peak-flow regression equations (issue #40), basin-characteristics method.

Source: Sando, Sando, McCarthy and Dutton, 2016, USGS SIR 2015-5019-F ver. 1.1
(February 2018): Table 1-4 (equations, n, model error variance, MVP, SEP) and
Table 1-5 (covariance), both in sir20155019F_tables.xlsx, and Table 3 (ranges).
Every literal number below is copied from that report, from its worked examples
(pp. 21-23), or from a dated live NSS response, never from ``MT.json`` itself: the
tests check the file against the source, not against itself.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from flowfreq.regression import OutOfRangeError, evaluate, load_state

AEPS = [0.667, 0.5, 0.429, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
REGIONS = ["ECP", "NEP", "NW", "NWF", "SEP", "SW", "UYCM", "W"]

# In range for every region (Table 3).
BASIN = {
    "CONTDA": 100.0,
    "PRECIP": 18.0,
    "FOREST": 45.0,
    "EL5000": 15.0,
    "EL6000": 50.0,
    "SLOP30_30M": 10.0,
    "ET0306MOD": 1.1,
}


@pytest.fixture(scope="module")
def lib():
    return load_state("MT")


def _sig3(x: float) -> float:
    return float(f"{x:.3g}")


class TestLibraryShape:
    def test_status_and_regions(self, lib):
        assert lib.status == "verified"
        assert lib.issue == 40
        assert lib.regions == REGIONS
        for region in lib.regions:
            assert lib.aeps(region) == AEPS
        assert len(lib.equations) == 80

    def test_every_equation_is_cited_and_complete(self, lib):
        for eq in lib.equations:
            assert "2015-5019-F" in eq.citation.publication
            assert "ver. 1.1" in eq.citation.publication
            assert "Table 1-4" in eq.citation.table and "Table 1-5" in eq.citation.table
            assert eq.model_error_variance is not None and eq.avp is not None
            assert eq.sep_log is not None and eq.covariance is not None
            u = np.asarray(eq.covariance)
            assert np.allclose(u, u.T)
            assert np.all(np.linalg.eigvalsh(u) > 0)
            assert eq.variables[0].code == "CONTDA"

    @pytest.mark.parametrize(
        "region, n_q66, n",
        [
            ("W", 113, 113),
            ("NW", 31, 32),
            ("NWF", 28, 31),
            ("NEP", 59, 64),
            ("ECP", 85, 90),
            ("SEP", 67, 68),
            ("UYCM", 87, 91),
            ("SW", 47, 48),
        ],
    )
    def test_n_sites(self, lib, region, n_q66, n):
        assert lib.equation(region, 0.667).n_sites == n_q66
        assert {e.n_sites for e in lib.equations if e.region_code == region and e.aep < 0.6} == {n}

    def test_transforms(self, lib):
        def tr(region):
            return {v.code: v.transform for v in lib.equation(region, 0.01).variables}

        assert tr("W") == {"CONTDA": "log10", "PRECIP": "log10", "FOREST": "log10_plus1"}
        assert tr("ECP") == {
            "CONTDA": "log10",
            "SLOP30_30M": "log10_plus1",
            "ET0306MOD": "log10",
        }
        assert tr("NW") == {"CONTDA": "log10"}

    def test_limits_from_table_3(self, lib):
        def limits(region):
            return {v.code: (v.minimum, v.maximum) for v in lib.equation(region, 0.01).variables}

        assert limits("W") == {
            "CONTDA": (0.60, 2465.66),
            "PRECIP": (14.62, 62.02),
            "FOREST": (20.40, 99.04),
        }
        assert limits("SEP") == {
            "CONTDA": (0.10, 1962.05),
            "FOREST": (0.00, 57.64),
            "ET0306MOD": (0.96, 1.67),
        }
        assert limits("NW") == {"CONTDA": (2.43, 1556.17)}

    def test_out_of_range_raises(self, lib):
        with pytest.raises(OutOfRangeError, match="PRECIP"):
            evaluate(lib.equation("NWF", 0.01), {"CONTDA": 10.0, "PRECIP": 40.0})


class TestAgainstTable14:
    """Evaluate the report's printed power-form equations directly."""

    @pytest.mark.parametrize(
        "region, aep, printed",
        [
            (
                "W",
                0.667,
                lambda b: 0.047
                * b["CONTDA"] ** 0.943
                * b["PRECIP"] ** 2.44
                * (b["FOREST"] + 1) ** -0.840,
            ),
            ("NW", 0.002, lambda b: 1171 * b["CONTDA"] ** 0.614),
            ("NWF", 0.2, lambda b: 0.303 * b["CONTDA"] ** 0.429 * b["PRECIP"] ** 2.03),
            ("NEP", 0.1, lambda b: 62.5 * b["CONTDA"] ** 0.617 * (b["EL5000"] + 1) ** -0.231),
            (
                "ECP",
                0.02,
                lambda b: 497
                * b["CONTDA"] ** 0.454
                * (b["SLOP30_30M"] + 1) ** 0.279
                * b["ET0306MOD"] ** -3.48,
            ),
            (
                "SEP",
                0.002,
                lambda b: 1240
                * b["CONTDA"] ** 0.511
                * (b["FOREST"] + 1) ** 0.012
                * b["ET0306MOD"] ** -4.68,
            ),
            ("UYCM", 0.5, lambda b: 4.73 * b["CONTDA"] ** 0.802 * (b["EL6000"] + 1) ** 0.194),
            # ver. 1.1 corrected this exponent's sign from +0.088 to -0.088
            ("SW", 0.2, lambda b: 13.4 * b["CONTDA"] ** 0.842 * (b["EL6000"] + 1) ** -0.088),
        ],
    )
    def test_power_form(self, lib, region, aep, printed):
        est = evaluate(lib.equation(region, aep), BASIN)
        assert est.flow_cfs == pytest.approx(printed(BASIN), rel=1e-12)

    @pytest.mark.parametrize(
        "region, aep, mev, mvp, sep_pct",
        [
            ("W", 0.01, 0.048, 0.051, 56.03),
            ("NW", 0.667, 0.150, 0.159, 115.0),
            ("NW", 0.01, 0.0, 0.004, 13.6),
            ("SEP", 0.667, 0.294, 0.316, 208.09),
            ("SW", 0.002, 0.083, 0.094, 80.31),
        ],
    )
    def test_error_statistics(self, lib, region, aep, mev, mvp, sep_pct):
        eq = lib.equation(region, aep)
        assert eq.model_error_variance == mev
        assert eq.avp == mvp
        assert 100 * math.sqrt(10 ** (math.log(10) * eq.sep_log**2) - 1) == pytest.approx(
            sep_pct, abs=1e-3
        )


class TestWorkedExamples:
    def test_case_1_southeast_plains_q1(self, lib):
        """Report p. 21-22: A = 27.70, ETSPR = 1.34, F = 24.80 -> Q1 = 970 ft3/s."""
        site = {"CONTDA": 27.70, "FOREST": 24.80, "ET0306MOD": 1.34}
        eq = lib.equation("SEP", 0.01)
        est = evaluate(eq, site)
        assert round(est.flow_cfs) == 970
        u = np.asarray(eq.covariance)
        # The report prints Table 1-5 rounded to 5 decimals ...
        assert np.round(u, 5).tolist() == [
            [0.01061, -0.00045, -0.00306, -0.04999],
            [-0.00045, 0.00173, -0.00074, -0.01283],
            [-0.00306, -0.00074, 0.00724, -0.00778],
            [-0.04999, -0.01283, -0.00778, 0.85177],
        ]
        # ... and with that and its rounded x0 = [1.0, 1.442, 1.412, 0.127] gets 0.00924.
        x_printed = np.array([1.0, 1.442, 1.412, 0.127])
        assert round(float(x_printed @ np.round(u, 5) @ x_printed), 5) == 0.00924
        # Full precision gives 0.00921, and the same SEP0 of 0.276 log units.
        x = np.array([1.0, math.log10(27.70), math.log10(25.80), math.log10(1.34)])
        xux = x @ u @ x
        assert xux == pytest.approx(0.009215, abs=1e-6)
        sep0 = math.sqrt(0.067 + xux)
        assert round(sep0, 3) == 0.276
        # eq. 8 applied to the rounded SEP0, as the report does
        assert round(100 * math.sqrt(math.exp((2.3026 * 0.276) ** 2) - 1), 1) == 70.5
        # t(0.05, 64) = 1.67 and CI = +/-0.461 log units (the report rounds log Q to
        # 2.99 before exponentiating, hence its 338-2,824 ft3/s)
        half = math.log10(est.interval[1] / est.flow_cfs)
        assert round(half, 3) == 0.461

    @pytest.mark.parametrize(
        "region, site, q10",
        [
            ("UYCM", {"CONTDA": 17.27, "EL6000": 0.00}, 339),
            ("ECP", {"CONTDA": 17.27, "SLOP30_30M": 0.05, "ET0306MOD": 1.15}, 420),
        ],
    )
    def test_case_2_q10(self, lib, region, site, q10):
        """Report p. 23: the two Q10 estimates that are area-weighted in case 2."""
        assert round(evaluate(lib.equation(region, 0.1), site).flow_cfs) == q10


# Live NSS responses for BASIN, fetched 2026-09-26 (NSS rounds to 3 significant figures).
NSS_RECORDED = [
    ("W", 0.667, 168),
    ("W", 0.01, 1250),
    ("NW", 0.5, 1240),
    ("NW", 0.002, 19800),
    ("NWF", 0.1, 1330),
    ("NEP", 0.01, 1540),
    ("ECP", 0.04, 3780),
    ("SEP", 0.667, 54.8),
    ("SEP", 0.002, 8740),
    ("UYCM", 0.2, 732),
    ("SW", 0.01, 1260),
    ("SW", 0.002, 1810),
]


@pytest.mark.parametrize("region, aep, nss", NSS_RECORDED)
def test_worked_example_matches_recorded_nss(lib, region, aep, nss):
    assert _sig3(evaluate(lib.equation(region, aep), BASIN).flow_cfs) == nss


NSS_CODE_TO_REGION = {
    "GC1675": "W",
    "GC1768": "NW",
    "GC1677": "NWF",
    "GC1678": "NEP",
    "GC1679": "ECP",
    "GC1680": "SEP",
    "GC1681": "UYCM",
    "GC1682": "SW",
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
            "CONTDA": 5.0,
            "PRECIP": 20.0,
            "FOREST": 30.0,
            "EL5000": 5.0,
            "EL6000": 10.0,
            "SLOP30_30M": 2.0,
            "ET0306MOD": 1.2,
        },
        BASIN,
        {
            "CONTDA": 1000.0,
            "PRECIP": 22.0,
            "FOREST": 55.0,
            "EL5000": 28.0,
            "EL6000": 90.0,
            "SLOP30_30M": 25.0,
            "ET0306MOD": 1.5,
        },
    ],
)
def test_matches_live_nss(lib, basin):
    """Every region x AEP against NSS Montana (region 30), citation 88.

    The channel-width regions (citation 163) need channel-width inputs, so NSS skips
    them; they are out of scope here. NSS returns three significant figures, so the
    check is exact agreement after rounding flowfreq's value the same way.
    """
    from flowfreq.streamstats import Characteristic, estimate_flow_statistics

    chars = {k: Characteristic(k, k, "", v, "") for k, v in basin.items()}
    results, skipped = estimate_flow_statistics("MT", chars, ["PFS"])
    assert not any(k.split(":")[-1] in NSS_CODE_TO_REGION for k in skipped)
    compared = 0
    for r in results:
        region = NSS_CODE_TO_REGION.get(r.region_code)
        if region is None:
            continue
        for code, est in r.estimates.items():
            ours = evaluate(lib.equation(region, NSS_STAT_TO_AEP[code]), basin).flow_cfs
            assert _sig3(ours) == est.value, (region, code, est.value, ours)
            assert abs(ours - est.value) / est.value < 5e-3
            compared += 1
    assert compared == 80
