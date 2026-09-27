"""Montana channel-width peak-flow equations and method weighting (issue #40).

Source: Chase, Sando, Armstrong and McCarthy, 2021, USGS SIR 2020-5142 ver. 1.1
(September 2021): Table 1 (width ranges), Tables 2-4 (equations and error
statistics for active-channel, bankfull and aerial-photo widths), Table 6
(cross-correlations), Tables 7-9 (covariance), equations 5-12 (weighting) and the
worked examples (pp. 39-47). Every literal number below is copied from that report
or from a dated live NSS response, never from ``MT.json`` itself.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from flowfreq.regression import OutOfRangeError, evaluate, load_state
from flowfreq.regression.montana import (
    combine_three,
    combine_two,
    cross_correlation,
    estimate_weighted,
    prediction_sep,
    region_code,
)

AEPS = [0.667, 0.5, 0.429, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
HYDRO = ["ECP", "NEP", "NW", "NWF", "SEP", "SW", "UYCM", "W"]
CW_REGIONS = sorted(f"{r}-{m}" for r in HYDRO for m in ("AC", "BF", "RS"))
WIDTH_CODE = {"AC": "WACTCH", "BF": "WBANKFULL", "RS": "CHANWD_RS"}


@pytest.fixture(scope="module")
def lib():
    return load_state("MT")


def _sig3(x: float) -> float:
    return float(f"{x:.3g}")


def _log_var(pct: float) -> float:
    return math.log(1 + (pct / 100) ** 2) / math.log(10) ** 2


class TestLibraryShape:
    def test_channel_width_regions(self, lib):
        assert lib.status == "verified"
        assert sorted(r for r in lib.regions if "-" in r) == CW_REGIONS
        for region in CW_REGIONS:
            assert lib.aeps(region) == AEPS
        assert len(lib.equations) == 80 + 240

    def test_lookup_is_unambiguous(self, lib):
        """The basin-characteristics and channel-width equations have distinct codes."""
        assert lib.equation("W", 0.01).variables[0].code == "CONTDA"
        assert lib.equation("W-AC", 0.01).variables[0].code == "WACTCH"
        assert region_code("W", "BC") == "W" and region_code("W", "RS") == "W-RS"
        with pytest.raises(ValueError, match="method"):
            region_code("W", "XX")
        with pytest.raises(ValueError, match="region"):
            region_code("ZZ", "AC")

    def test_every_equation_is_cited_and_complete(self, lib):
        for eq in lib.equations:
            if eq.region_code not in CW_REGIONS:
                continue
            assert "2020-5142" in eq.citation.publication
            assert "ver. 1.1" in eq.citation.publication
            method = eq.region_code.split("-")[1]
            table = {"AC": "Table 2", "BF": "Table 3", "RS": "Table 4"}[method]
            cov = {"AC": "Table 7", "BF": "Table 8", "RS": "Table 9"}[method]
            assert table in eq.citation.table and cov in eq.citation.table
            assert [v.code for v in eq.variables] == [WIDTH_CODE[method]]
            assert eq.variables[0].transform == "log10" and eq.variables[0].units == "ft"
            u = np.asarray(eq.covariance)
            assert u.shape == (2, 2) and np.allclose(u, u.T)
            assert np.all(np.linalg.eigvalsh(u) > 0)

    @pytest.mark.parametrize(
        "region, lo_hi",
        [
            ("W", {"AC": (3.0, 213.0), "BF": (5.0, 246.0), "RS": (2.3, 203.8)}),
            ("NW", {"AC": (9.5, 234.5), "BF": (11.5, 273.5), "RS": (8.9, 257.2)}),
            ("ECP", {"AC": (2.0, 91.0), "BF": (3.5, 220.0), "RS": (2.7, 47.4)}),
            ("SW", {"AC": (1.8, 223.0), "BF": (3.5, 260.0), "RS": (2.6, 219.2)}),
        ],
    )
    def test_limits_from_table_1(self, lib, region, lo_hi):
        for method, (lo, hi) in lo_hi.items():
            v = lib.equation(f"{region}-{method}", 0.01).variables[0]
            assert (v.minimum, v.maximum) == (lo, hi)

    def test_out_of_range_raises(self, lib):
        with pytest.raises(OutOfRangeError, match="CHANWD_RS"):
            evaluate(lib.equation("ECP-RS", 0.01), {"CHANWD_RS": 60.0})


class TestAgainstTables2To4:
    @pytest.mark.parametrize(
        "region, aep, width, k, b",
        [
            ("W-AC", 0.667, 12.0, 0.451, 1.90),
            ("NW-AC", 0.002, 50.0, 454, 0.952),
            ("UYCM-AC", 0.5, 20.0, 0.911, 1.76),  # ver. 1.1 revision (1.77 -> 1.76)
            ("NEP-BF", 0.1, 30.0, 6.96, 1.49),
            ("NW-BF", 0.002, 50.0, 316, 0.980),
            ("SW-BF", 0.01, 40.0, 14.0, 1.26),
            ("W-RS", 0.667, 12.0, 1.05, 1.69),  # ver. 1.1 revision (1.70 -> 1.69)
            ("NWF-RS", 0.5, 12.0, 2.35, 1.42),  # ver. 1.1 revision (1.41 -> 1.42)
            ("NWF-RS", 0.002, 12.0, 213, 1.16),  # ver. 1.1 revision (1.17 -> 1.16)
            ("SW-RS", 0.5, 12.0, 1.46, 1.51),  # ver. 1.1 revision (1.56 -> 1.51)
            ("ECP-RS", 0.02, 20.0, 83.1, 1.2),
        ],
    )
    def test_power_form(self, lib, region, aep, width, k, b):
        eq = lib.equation(region, aep)
        est = evaluate(eq, {eq.variables[0].code: width})
        assert est.flow_cfs == pytest.approx(k * width**b, rel=1e-12)

    @pytest.mark.parametrize(
        "region, aep, n, sigma2, mev, mvp, sep_pct",
        [
            ("W-AC", 0.01, 103, 0.048, 0.018, 0.050, 66.1),
            ("NW-AC", 0.01, 28, 0.000, 0.009, 0.003, 25.5),  # WLS, sigma2 printed 0.000
            ("NWF-AC", 0.667, 21, 0.081, 0.018, 0.092, 89.2),  # smaller Q66.7 n
            ("W-BF", 0.01, 100, 0.050, 0.028, 0.053, 72.8),
            ("ECP-BF", 0.1, 80, 0.082, 0.017, 0.088, 86.0),
            ("UYCM-BF", 0.1, 81, 0.058, 0.020, 0.062, 73.8),
            ("SEP-RS", 0.667, 65, 0.411, 0.035, 0.427, 325),
            ("SW-RS", 0.002, 46, 0.172, 0.005, 0.183, 131),
        ],
    )
    def test_error_statistics(self, lib, region, aep, n, sigma2, mev, mvp, sep_pct):
        """model_error_variance is sigma2 + MEV and avp is MVP + MEV (MT.json notes)."""
        eq = lib.equation(region, aep)
        assert eq.n_sites == n
        assert eq.model_error_variance == pytest.approx(sigma2 + mev, abs=1e-12)
        assert eq.avp == pytest.approx(mvp + mev, abs=1e-12)
        assert eq.sep_log == pytest.approx(math.sqrt(_log_var(sep_pct)), abs=1e-6)
        # The report's SEP includes measurement error: SEP^2 = MVP + MEV to the printing.
        assert _log_var(sep_pct) == pytest.approx(mvp + mev, abs=2e-3)


class TestCovariance:
    @pytest.mark.parametrize(
        "region, aep, printed",
        [
            ("W-AC", 0.01, ((0.01258, -0.00687), (-0.00687, 0.00418))),
            ("SW-AC", 0.002, ((0.04725, -0.03021), (-0.03021, 0.02158))),
            ("W-BF", 0.01, ((0.01673, -0.00883), (-0.00883, 0.00507))),
            ("NWF-BF", 0.002, ((0.12534, -0.07883), (-0.07883, 0.05797))),
            ("SEP-RS", 0.667, ((0.07333, -0.06653), (-0.06653, 0.06860))),
            ("UYCM-RS", 0.1, ((0.02715, -0.01628), (-0.01628, 0.01125))),
        ],
    )
    def test_tables_7_to_9(self, lib, region, aep, printed):
        assert lib.equation(region, aep).covariance == printed

    def test_example_1_site_seps(self, lib):
        """Report pp. 40-41: x2'Ux2 = 0.0026, SEP2 = 0.262; x3'Ux3 = 0.00255, SEP3 = 0.283."""
        eq = lib.equation("W-AC", 0.01)
        x = np.array([1.0, math.log10(12)])
        assert x @ np.asarray(eq.covariance) @ x == pytest.approx(0.0026, abs=5e-5)
        assert prediction_sep(eq, {"WACTCH": 12.0}) == pytest.approx(0.262, abs=5e-4)
        assert _sig3(evaluate(eq, {"WACTCH": 12.0}).flow_cfs) == 284
        eq = lib.equation("W-BF", 0.01)
        x = np.array([1.0, math.log10(18)])
        assert x @ np.asarray(eq.covariance) @ x == pytest.approx(0.00255, abs=5e-6)
        # The report prints 0.283 for sqrt(0.081) = 0.2838 (it truncates).
        assert prediction_sep(eq, {"WBANKFULL": 18.0}) == pytest.approx(0.2838, abs=5e-4)
        assert _sig3(evaluate(eq, {"WBANKFULL": 18.0}).flow_cfs) == 317

    def test_site_interval_uses_mev(self, lib):
        """evaluate()'s site interval is sigma2 + MEV + x'Ux, the report's eq. 3 as modified."""
        eq = lib.equation("W-AC", 0.01)
        est = evaluate(eq, {"WACTCH": 12.0})
        assert est.interval_method.startswith("site-specific")
        df = 103 - 2
        half = math.log10(est.interval[1] / est.flow_cfs)
        from scipy import stats

        assert half / stats.t.ppf(0.95, df) == pytest.approx(
            math.sqrt(0.048 + 0.018 + 0.0026), abs=5e-4
        )

    def test_prediction_sep_needs_covariance(self, lib):
        from dataclasses import replace

        eq = replace(lib.equation("W-AC", 0.01), covariance=None)
        with pytest.raises(ValueError, match="covariance"):
            prediction_sep(eq, {"WACTCH": 12.0})


class TestCrossCorrelation:
    @pytest.mark.parametrize(
        "region, aep, m1, m2, r",
        [
            ("W", 0.01, "BC", "AC", 0.553),
            ("W", 0.01, "BC", "BF", 0.518),
            ("W", 0.01, "AC", "BF", 0.907),
            ("ECP", 0.1, "BC", "BF", 0.642),
            ("UYCM", 0.1, "BF", "BC", 0.488),  # either order
            ("SEP", 0.667, "BC", "RS", 0.658),
            ("NW", 0.1, "BC", "RS", 0.009),
            ("SW", 0.002, "BF", "RS", 0.896),
        ],
    )
    def test_table_6(self, region, aep, m1, m2, r):
        assert cross_correlation(region, aep, m1, m2) == r

    def test_unknown_lookups_raise(self):
        with pytest.raises(KeyError):
            cross_correlation("W", 0.3, "BC", "AC")
        with pytest.raises(KeyError):
            cross_correlation("W", 0.01, "AC", "AC")
        with pytest.raises(KeyError):
            cross_correlation("XX", 0.01, "BC", "AC")


class TestWeightingEquations:
    def test_equation_11_12_report_arithmetic(self):
        """Example 2 (p. 44) with the report's rounded inputs: a1 = 0.871, SEPz = 0.237."""
        z, sepz, (a1, a2) = combine_two(math.log10(420), 0.239, math.log10(384), 0.321, 0.642)
        assert round(a1, 2) == 0.87 and a1 + a2 == pytest.approx(1.0)
        assert round(sepz, 3) == 0.237
        assert _sig3(10**z) == 415

    def test_equations_6_to_8_report_arithmetic(self):
        """Example 1 (p. 42): the report rounds A, B, C to 0.067, 0.013, 0.015.

        With those the printed formulas give a1 = 0.623, a2 = 0.332, a3 = 0.045. The
        system is ill-conditioned (A C - B^2 = 0.000836), so the unrounded weights
        differ (0.624, 0.356, 0.019) while Z and SEPz barely move.
        """
        a_, b_, c_ = 0.067, 0.013, 0.015
        s13, s23 = 0.034, 0.067
        det = a_ * c_ - b_**2
        a1 = (c_ * (0.283**2 - s13) - b_ * (0.283**2 - s23)) / det
        a2 = (a_ * (0.283**2 - s23) - b_ * (0.283**2 - s13)) / det
        a1, a2 = round(a1, 3), round(a2, 3)
        assert (a1, a2, round(1 - a1 - a2, 3)) == (0.623, 0.332, 0.045)  # eq. 8 on rounded
        z, sepz, w = combine_three(
            [math.log10(551), math.log10(284), math.log10(317)],
            [0.234, 0.262, 0.283],
            0.553,
            0.518,
            0.907,
        )
        assert _sig3(10**z) == pytest.approx(432, rel=0.005)
        assert round(sepz, 3) == 0.217

    def test_three_way_weights_are_minimum_variance(self):
        """Eqs. 6-8 are the closed form of the GLS weights S^-1 1 / (1' S^-1 1)."""
        sep = np.array([0.234, 0.262, 0.283])
        r = np.array([[1, 0.553, 0.518], [0.553, 1, 0.907], [0.518, 0.907, 1]])
        s = r * np.outer(sep, sep)
        w_gls = np.linalg.solve(s, np.ones(3))
        w_gls /= w_gls.sum()
        _, sepz, w = combine_three([0.0, 0.0, 0.0], list(sep), 0.553, 0.518, 0.907)
        assert np.allclose(w, w_gls, atol=1e-12)
        assert sepz**2 == pytest.approx(1.0 / (np.ones(3) @ np.linalg.solve(s, np.ones(3))))

    def test_invalid_inputs_raise(self):
        with pytest.raises(ValueError):
            combine_two(1.0, 0.0, 1.0, 0.2, 0.5)
        with pytest.raises(ValueError):
            combine_two(1.0, 0.2, 1.0, 0.2, 1.0)
        with pytest.raises(ValueError):
            combine_three([1.0, 1.0], [0.2, 0.2], 0.1, 0.1, 0.1)


class TestWorkedExamples:
    def test_example_1_west_q1_three_methods(self, lib):
        """Report pp. 39-42: 432 cfs, SEPz 0.217.

        Estimate 1 is 558 cfs here against the report's 551: the report multiplies
        rounded factors (21.5 x 7.66 x 93.0 x 0.036); its SEP1 = 0.234 agrees.
        """
        site = {"CONTDA": 12.5, "PRECIP": 31.0, "FOREST": 38.0, "WACTCH": 12.0}
        site["WBANKFULL"] = 18.0
        w = estimate_weighted("W", 0.01, site, ["BC", "AC", "BF"], lib=lib)
        assert w.estimates["BC"].flow_cfs == pytest.approx(551, rel=0.015)
        assert w.seps["BC"] == pytest.approx(0.234, abs=5e-4)
        assert _sig3(w.estimates["AC"].flow_cfs) == 284
        assert _sig3(w.estimates["BF"].flow_cfs) == 317
        assert sum(w.weights.values()) == pytest.approx(1.0)
        assert w.flow_cfs == pytest.approx(432, rel=0.005)
        assert round(w.sep_log, 3) == 0.217
        assert not w.outside_range

    def test_example_2_crossing_regions(self, lib):
        """Report pp. 43-45: ECP 415 (SEP 0.237), UYCM 257 (0.249), area-weighted 358 (0.242).

        UYCM comes out 253 here: the report rounds Z to 2.41 before taking the
        antilogarithm (Z = 2.403 unrounded, even from its own rounded inputs).
        """
        site = {"CONTDA": 17.27, "SLOP30_30M": 0.05, "ET0306MOD": 1.15, "EL6000": 0.0}
        site["WBANKFULL"] = 12.0
        ecp = estimate_weighted("ECP", 0.1, site, ["BC", "BF"], lib=lib)
        assert _sig3(ecp.estimates["BC"].flow_cfs) == 420
        assert _sig3(ecp.estimates["BF"].flow_cfs) == 384
        assert round(ecp.seps["BC"], 3) == 0.239 and round(ecp.seps["BF"], 3) == 0.321
        assert _sig3(ecp.flow_cfs) == 415 and round(ecp.sep_log, 3) == 0.237
        uycm = estimate_weighted("UYCM", 0.1, site, ["BC", "BF"], lib=lib)
        assert _sig3(uycm.estimates["BC"].flow_cfs) == 339
        assert _sig3(uycm.estimates["BF"].flow_cfs) == 189
        assert round(uycm.seps["BC"], 3) == 0.288 and round(uycm.seps["BF"], 3) == 0.288
        assert uycm.weights["BC"] == pytest.approx(0.500, abs=5e-3)
        assert round(uycm.sep_log, 3) == 0.249
        assert uycm.flow_cfs == pytest.approx(257, rel=0.02)
        q = 0.36 * uycm.flow_cfs + 0.64 * ecp.flow_cfs
        assert q == pytest.approx(358, rel=0.005)
        assert 0.36 * uycm.sep_log + 0.64 * ecp.sep_log == pytest.approx(0.242, abs=1.5e-3)

    def test_example_3_southeast_plains_q66_7(self, lib):
        """Report pp. 45-47: a1 = 0.779, 44.7 cfs, SEPz 0.543 (the report says '1.5-percent
        AEP' for Q66.7, the 66.7-percent AEP; its equations are Q66.7's)."""
        site = {"CONTDA": 54.7, "FOREST": 18.2, "ET0306MOD": 1.26, "CHANWD_RS": 15.1}
        w = estimate_weighted("SEP", 0.667, site, ["BC", "RS"], lib=lib)
        assert w.estimates["BC"].flow_cfs == pytest.approx(40.5, rel=0.01)
        assert _sig3(w.estimates["RS"].flow_cfs) == 63.4
        assert round(w.seps["BC"], 3) == 0.553 or round(w.seps["BC"], 3) == 0.554
        assert round(w.seps["RS"], 3) == 0.677
        assert round(w.weights["BC"], 3) == 0.779 and round(w.weights["RS"], 3) == 0.221
        assert w.flow_cfs == pytest.approx(44.7, rel=0.01)
        assert w.sep_log == pytest.approx(0.543, abs=2e-3)

    def test_rejects_bad_method_lists(self, lib):
        site = {"WACTCH": 12.0, "WBANKFULL": 18.0}
        with pytest.raises(ValueError, match="two or three"):
            estimate_weighted("W", 0.01, site, ["AC"], lib=lib)
        with pytest.raises(ValueError, match="two or three"):
            estimate_weighted("W", 0.01, site, ["AC", "AC"], lib=lib)
        with pytest.raises(ValueError, match="Montana"):
            estimate_weighted("W", 0.01, site, ["AC", "BF"], lib=load_state("WA"))


# Live NSS responses (Montana, citation 163) for Wac 25, Wbf 30, Wrs 20 ft, fetched
# 2026-09-27. NSS rounds to 3 significant figures.
MID = {"WACTCH": 25.0, "WBANKFULL": 30.0, "CHANWD_RS": 20.0}
NSS_RECORDED = [
    ("W-AC", 0.667, 204),
    ("W-AC", 0.01, 899),
    ("W-AC", 0.002, 1170),
    ("NW-BF", 0.667, 176),
    ("NW-BF", 0.01, 3560),
    ("NW-BF", 0.002, 8860),
    ("SW-BF", 0.01, 1020),
    ("UYCM-AC", 0.01, 2140),
    ("ECP-RS", 0.01, 4140),
    ("NEP-RS", 0.002, 4890),
]


@pytest.mark.parametrize("region, aep, nss", NSS_RECORDED)
def test_matches_recorded_nss(lib, region, aep, nss):
    eq = lib.equation(region, aep)
    assert _sig3(evaluate(eq, {eq.variables[0].code: MID[eq.variables[0].code]}).flow_cfs) == nss


NSS_STAT_SUFFIX_TO_AEP = {
    "66_7AE": 0.667,
    "50AEP": 0.5,
    "42_9AE": 0.429,
    "20AEP": 0.2,
    "10AEP": 0.1,
    "4AEP": 0.04,
    "2AEP": 0.02,
    "1AEP": 0.01,
    "0_5AEP": 0.005,
    "0_2AEP": 0.002,
}
NSS_METHOD = {"Active": "AC", "Act": "AC", "Bankfull": "BF", "Aerial": "RS"}


@pytest.mark.requires_network
@pytest.mark.parametrize(
    "basin",
    [
        {"WACTCH": 10.0, "WBANKFULL": 12.0, "CHANWD_RS": 9.0},
        MID,
        {"WACTCH": 70.0, "WBANKFULL": 90.0, "CHANWD_RS": 45.0},
    ],
)
def test_matches_live_nss(lib, basin):
    """All 24 channel-width regions x 10 AEPs against NSS Montana, citation 163.

    NSS rounds to 3 significant figures; flowfreq must round to the same value, except at
    an exact tie (e.g. 3.14 x 25^1.5 = 392.5), where either rounding is correct.
    """
    from flowfreq.streamstats import Characteristic, estimate_flow_statistics

    gc_to_region = {
        e.region_name.split("NSS ")[1].split(",")[0]: e.region_code
        for e in lib.equations
        if e.region_code in CW_REGIONS
    }
    chars = {k: Characteristic(k, k, "", v, "") for k, v in basin.items()}
    results, _skipped = estimate_flow_statistics("MT", chars, ["PFS"])
    compared = 0
    for r in results:
        region = gc_to_region.get(r.region_code)
        if region is None:
            continue
        for code, est in r.estimates.items():
            eq = lib.equation(region, NSS_STAT_SUFFIX_TO_AEP[code[4:]])
            ours = evaluate(eq, basin).flow_cfs
            unit = 10 ** (math.floor(math.log10(est.value)) - 2)
            assert abs(ours - est.value) <= 0.5 * unit * (1 + 1e-9), (region, code, ours)
            assert abs(ours - est.value) / est.value < 5e-3
            compared += 1
    assert compared == 240
