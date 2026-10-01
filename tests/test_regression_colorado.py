"""Colorado peak-flow regression equations (issue #93).

Sources, as NSS combines them for Colorado's peak-flow statistic group:

- Capesius and Stephens, 2009, USGS SIR 2009-5136 ver. 1.2 (April 2014), figures 3-6:
  Mountain, Northwest, Rio Grande and Southwest regions.
- Kohn and others, 2016, USGS SIR 2016-5099, figures 9-10 and Appendix 6 (WREG output,
  covariance): Foothills and Plains regions.
- Waltemeyer, 2006, USGS SIR 2006-5306, Tables 1-2: Navajo Nation region 8 and High
  Elevation region.

Every literal number below is copied from those reports, or from a dated live NSS
response, never from ``CO.json`` itself: the tests check the file against the source,
not against itself.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from flowfreq.regression import EquationsUnavailable, OutOfRangeError, evaluate, load_state

AEPS_8 = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
AEPS_NAVAJO = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002]

# One basin inside every region's range (NSS's too).
BASIN = {
    "DRNAREA": 100.0,
    "BSLDEM10M": 15.0,
    "PRECIP": 30.0,
    "EL7500": 50.0,
    "STATSCLAY": 25.0,
    "I6H100Y": 3.5,
    "OUTLETELEV": 6500.0,
    "ELEV": 6500.0,
    "BSLDEM30ff": 0.15,
}


@pytest.fixture(scope="module")
def lib():
    return load_state("CO")


def _sig3(x: float) -> float:
    return float(f"{x:.3g}")


def _pct(sep_log: float) -> float:
    return 100 * math.sqrt(10 ** (math.log(10) * sep_log**2) - 1)


class TestLibraryShape:
    def test_status_and_regions(self, lib):
        assert lib.status == "verified"
        assert lib.issue == 93
        assert lib.regions == [
            "Foothills",
            "Mountain",
            "Navajo8",
            "NavajoHighElev",
            "Northwest",
            "Plains",
            "RioGrande",
            "Southwest",
        ]
        for region in ("Mountain", "Northwest", "RioGrande", "Southwest", "Foothills", "Plains"):
            assert lib.aeps(region) == AEPS_8
        for region in ("Navajo8", "NavajoHighElev"):
            assert lib.aeps(region) == AEPS_NAVAJO
        assert len(lib.equations) == 62
        assert len(lib.source_reports) == 3

    def test_navajo_has_no_200_year_equation(self, lib):
        with pytest.raises(EquationsUnavailable):
            lib.equation("Navajo8", 0.005)

    @pytest.mark.parametrize(
        "region, n",
        [
            ("Mountain", 141),
            ("Northwest", 90),
            ("RioGrande", 44),
            ("Southwest", 78),
            ("Foothills", 89),
            ("Plains", 99),
            ("Navajo8", 62),
            ("NavajoHighElev", 27),
        ],
    )
    def test_n_sites(self, lib, region, n):
        assert {e.n_sites for e in lib.equations if e.region_code == region} == {n}

    @pytest.mark.parametrize(
        "region, limits",
        [
            ("Mountain", {"DRNAREA": (1, 1060), "BSLDEM10M": (7.6, 60.2), "PRECIP": (18, 47)}),
            ("Northwest", {"DRNAREA": (1, 5250), "EL7500": (0, 100), "PRECIP": (8, 49)}),
            ("RioGrande", {"DRNAREA": (2, 517), "PRECIP": (19, 45)}),
            ("Southwest", {"DRNAREA": (1, 4390), "EL7500": (0, 100)}),
            (
                "Foothills",
                {
                    "DRNAREA": (0.60, 2850),
                    "I6H100Y": (2.38, 4.89),
                    "STATSCLAY": (9.87, 37.5),
                    "OUTLETELEV": (4290, 8270),
                },
            ),
            (
                "Plains",
                {"DRNAREA": (0.26, 3560), "BSLDEM10M": (0.41, 21.9), "STATSCLAY": (5.20, 38.5)},
            ),
            (
                "Navajo8",
                {"DRNAREA": (0.06, 4350), "BSLDEM30ff": (0.0231, 0.4447), "ELEV": (4320, 7750)},
            ),
            ("NavajoHighElev", {"DRNAREA": (1.7, 1200)}),
        ],
    )
    def test_limits(self, lib, region, limits):
        eq = lib.equation(region, 0.01)
        assert {v.code: (v.minimum, v.maximum) for v in eq.variables} == limits

    def test_a7500_is_percent_plus_one(self, lib):
        # Figure 2: "Percentage of A above 7,500 feet of elevation plus 1".
        for region in ("Northwest", "Southwest"):
            for aep in AEPS_8:
                v = {v.code: v for v in lib.equation(region, aep).variables}["EL7500"]
                assert (v.transform, v.scale, v.offset) == ("log10", 1.0, 1.0)
        # Zero percent above 7,500 ft is valid and in range.
        est = evaluate(lib.equation("Southwest", 0.01), {"DRNAREA": 100.0, "EL7500": 0.0})
        assert est.flow_cfs == pytest.approx(10**2.91 * 100**0.59, rel=1e-12)

    def test_out_of_range_raises(self, lib):
        with pytest.raises(OutOfRangeError, match="PRECIP"):
            evaluate(lib.equation("RioGrande", 0.01), {"DRNAREA": 100.0, "PRECIP": 50.0})
        with pytest.raises(OutOfRangeError, match="ELEV=7800.0 > max 7750"):
            evaluate(lib.equation("Navajo8", 0.01), {**BASIN, "ELEV": 7800.0})


class TestPrintedForms:
    """Evaluate the reports' printed equations directly."""

    @pytest.mark.parametrize(
        "region, aep, printed",
        [
            (
                "Mountain",
                0.5,
                lambda b: 10**-2.05
                * b["DRNAREA"] ** 0.78
                * b["BSLDEM10M"] ** 0.17
                * b["PRECIP"] ** 2.10,
            ),
            (
                "Mountain",
                0.002,
                lambda b: 10**-0.06
                * b["DRNAREA"] ** 0.74
                * b["BSLDEM10M"] ** 0.15
                * b["PRECIP"] ** 1.14,
            ),
            (
                "Northwest",
                0.01,
                lambda b: 10**0.93
                * b["DRNAREA"] ** 0.74
                * (b["EL7500"] + 1) ** -0.81
                * b["PRECIP"] ** 1.65,
            ),
            ("RioGrande", 0.5, lambda b: 10**-3.00 * b["DRNAREA"] ** 1.00 * b["PRECIP"] ** 2.46),
            (
                "Southwest",
                0.002,
                lambda b: 10**3.21 * b["DRNAREA"] ** 0.58 * (b["EL7500"] + 1) ** -0.39,
            ),
            (
                "Foothills",
                0.01,
                lambda b: 10**13.244
                * b["DRNAREA"] ** 0.572
                * b["I6H100Y"] ** 3.190
                * b["STATSCLAY"] ** 1.013
                * b["OUTLETELEV"] ** -3.631,
            ),
            (
                "Plains",
                0.002,
                lambda b: 10**-0.473
                * b["DRNAREA"] ** 0.412
                * b["BSLDEM10M"] ** 0.567
                * b["STATSCLAY"] ** 2.567,
            ),
            ("Navajo8", 0.5, lambda b: 1.08e7 * b["DRNAREA"] ** 0.457 * b["ELEV"] ** -1.35),
            (
                "Navajo8",
                0.01,
                lambda b: 1.53e12
                * b["DRNAREA"] ** 0.265
                * b["BSLDEM30ff"] ** 0.477
                * b["ELEV"] ** -2.25,
            ),
            ("NavajoHighElev", 0.02, lambda b: 1.15e2 * b["DRNAREA"] ** 0.623),
        ],
    )
    def test_power_form(self, lib, region, aep, printed):
        est = evaluate(lib.equation(region, aep), BASIN)
        assert est.flow_cfs == pytest.approx(printed(BASIN), rel=1e-12)

    @pytest.mark.parametrize(
        "region, aep, sep_pct, sme_pct",
        [
            ("Mountain", 0.5, 49, 48),
            ("Northwest", 0.5, 113, 108),
            ("RioGrande", 0.01, 51, 48),
            ("Southwest", 0.002, 75, 70),
            ("Foothills", 0.5, 117, 111),
            ("Plains", 0.002, 170, 160),
        ],
    )
    def test_percent_error_statistics(self, lib, region, aep, sep_pct, sme_pct):
        eq = lib.equation(region, aep)
        assert _pct(eq.sep_log) == pytest.approx(sep_pct, abs=1e-3)
        if region in ("Foothills", "Plains"):
            # WREG's model error variance (log units, 3 decimals) rounds to the printed SME.
            assert round(_pct(math.sqrt(eq.model_error_variance))) == pytest.approx(sme_pct, abs=1)
        else:
            assert _pct(math.sqrt(eq.model_error_variance)) == pytest.approx(sme_pct, abs=1e-3)
            assert eq.avp is None and eq.covariance is None

    @pytest.mark.parametrize(
        "region, aep, reg_log, pred_log",
        [
            ("Navajo8", 0.5, 0.386, 0.397),
            ("Navajo8", 0.002, 0.237, 0.256),
            ("NavajoHighElev", 0.1, 0.149, 0.160),
            ("NavajoHighElev", 0.002, 0.205, 0.222),
        ],
    )
    def test_navajo_log_unit_errors(self, lib, region, aep, reg_log, pred_log):
        eq = lib.equation(region, aep)
        assert eq.sep_log == pred_log
        assert eq.model_error_variance == pytest.approx(reg_log**2, rel=1e-12)
        assert eq.covariance is None


class TestCovariance2016:
    def test_foothills_q1_from_appendix_6(self, lib):
        eq = lib.equation("Foothills", 0.01)
        assert eq.model_error_variance == 0.099
        assert eq.avp == 0.108
        assert eq.covariance[0] == (
            7.3152209,
            -5.1035335e-2,
            -7.4413241e-1,
            -6.9581139e-2,
            -1.7808901,
        )
        assert eq.covariance[4][4] == 4.4838012e-1

    def test_plains_q10_from_appendix_6(self, lib):
        eq = lib.equation("Plains", 0.1)
        assert eq.covariance == (
            (2.8003530e-01, -8.9194229e-03, -2.1683035e-02, -1.8327544e-01),
            (-8.9194229e-03, 1.8579332e-03, 7.9705313e-05, 3.5475499e-03),
            (-2.1683035e-02, 7.9705313e-05, 1.5275367e-02, 9.4109449e-03),
            (-1.8327544e-01, 3.5475499e-03, 9.4109449e-03, 1.2760268e-01),
        )

    def test_every_2016_matrix_is_symmetric_positive_definite(self, lib):
        for eq in lib.equations:
            if eq.region_code in ("Foothills", "Plains"):
                u = np.asarray(eq.covariance)
                assert np.array_equal(u, u.T)
                assert np.all(np.linalg.eigvalsh(u) > 0)

    def test_variance_of_prediction_at_a_calibration_gage(self, lib):
        """Appendix 5: VPP(g)r at 06655000 (A 196, 6P100 2.59, C 17.322, Eout 4,450 from
        Appendix 6 SiteInfo.txt) is 0.176 at 50% AEP and 0.119 at 1% AEP."""
        site = {"DRNAREA": 196.0, "I6H100Y": 2.59, "STATSCLAY": 17.322, "OUTLETELEV": 4450.0}
        for aep, vp in ((0.5, 0.176), (0.01, 0.119)):
            eq = lib.equation("Foothills", aep)
            x = np.array([1.0] + [math.log10(site[v.code]) for v in eq.variables])
            assert eq.model_error_variance + x @ np.array(eq.covariance) @ x == pytest.approx(
                vp, abs=5e-4
            )

    def test_site_specific_interval(self, lib):
        eq = lib.equation("Plains", 0.1)
        est = evaluate(eq, BASIN)
        x = np.array([1.0] + [math.log10(BASIN[c]) for c in ("DRNAREA", "BSLDEM10M", "STATSCLAY")])
        sd = math.sqrt(0.127 + x @ np.array(eq.covariance) @ x)  # Plains Q10 MEV, Appendix 6
        half = stats.t.ppf(0.95, 99 - 4) * sd
        assert est.interval[1] == pytest.approx(10 ** (est.log_flow + half), rel=1e-10)
        assert est.interval_method.startswith("site-specific")

    def test_2009_and_navajo_intervals_use_sep(self, lib):
        est = evaluate(lib.equation("Mountain", 0.01), BASIN)
        assert est.interval_method == "average standard error of prediction"


class TestWorkedExamples:
    @pytest.mark.parametrize(
        "region, site, values",
        [
            # SIR 2016-5099 Appendix 5, QP(g)r; characteristics from Appendix 6 SiteInfo.txt.
            (
                "Foothills",
                {"DRNAREA": 196.0, "I6H100Y": 2.59, "STATSCLAY": 17.322, "OUTLETELEV": 4450.0},
                [760, 1770, 2770, 4410, 5860, 7610, 9580, 12500],
            ),
            (
                "Plains",
                {"DRNAREA": 897.9, "BSLDEM10M": 5.528, "STATSCLAY": 5.200},
                [68.8, 142, 208, 317, 421, 556, 723, 1010],
            ),
        ],
    )
    def test_appendix_5_regression_estimates(self, lib, region, site, values):
        for aep, q in zip(AEPS_8, values, strict=True):
            assert _sig3(evaluate(lib.equation(region, aep), site).flow_cfs) == q

    def test_navajo_appendix_1_predicted_values(self, lib):
        """SIR 2006-5306 Appendix 1: the 'predicted' line of High Elevation gage 09177500
        (A 15.3 mi2) is Table 2's equations."""
        predicted = [88, 196, 299, 471, 631, 822, 1410]
        for aep, q in zip(AEPS_NAVAJO, predicted, strict=True):
            ours = evaluate(lib.equation("NavajoHighElev", aep), {"DRNAREA": 15.3}).flow_cfs
            assert ours == pytest.approx(q, rel=0.01), (aep, ours, q)

    def test_navajo_appendix_1_region_8_uses_unrounded_coefficients(self, lib):
        """Table 2's rounded region 8 equations sit a constant factor below Appendix 1's
        predicted values (Q10: median ratio 0.959 over the 62 gages, spread < 1%), which a
        refit attributes to coefficient rounding. 09314250: A 946, S 0.2490, E 7,660,
        predicted Q10 3,880."""
        site = {"DRNAREA": 946.0, "BSLDEM30ff": 0.2490, "ELEV": 7660.0}
        ours = evaluate(lib.equation("Navajo8", 0.1), site).flow_cfs
        assert ours / 3880 == pytest.approx(0.956, abs=0.002)

    def test_navajo_worked_example_is_inconsistent(self, lib):
        """p. 9 prints 8,510 ft3/s for region 8 Q100 at A 523, S 0.2422, E 6,940, but the
        equation it restates, 1.53e12 A^0.265 S^0.477 E^-2.25, gives 9,296."""
        site = {"DRNAREA": 523.0, "BSLDEM30ff": 0.2422, "ELEV": 6940.0}
        q = evaluate(lib.equation("Navajo8", 0.01), site).flow_cfs
        assert q == pytest.approx(1.53e12 * 523**0.265 * 0.2422**0.477 * 6940**-2.25, rel=1e-12)
        assert round(q) == 9296
        assert q / 8510 == pytest.approx(1.092, abs=1e-3)


# Live NSS responses for BASIN, fetched 2026-10-01 (NSS rounds to 3 significant figures).
NSS_RECORDED = [
    ("Mountain", 0.5, 648),
    ("Mountain", 0.01, 1580),
    ("Mountain", 0.002, 1910),
    ("Northwest", 0.5, 669),
    ("Northwest", 0.01, 2910),
    ("RioGrande", 0.01, 1900),
    ("Southwest", 0.002, 5060),
    ("Foothills", 0.5, 361),
    ("Foothills", 0.01, 4950),
    ("Plains", 0.01, 20100),
    ("Plains", 0.002, 40400),
    ("Navajo8", 0.5, 631),
    ("Navajo8", 0.01, 5530),
    ("NavajoHighElev", 0.01, 2540),
    ("NavajoHighElev", 0.002, 4000),
]


@pytest.mark.parametrize("region, aep, nss", NSS_RECORDED)
def test_worked_example_matches_recorded_nss(lib, region, aep, nss):
    assert _sig3(evaluate(lib.equation(region, aep), BASIN).flow_cfs) == nss


NSS_CODE_TO_REGION = {
    "GC1204": "Mountain",
    "GC1205": "Northwest",
    "GC1206": "Plains",
    "GC1207": "RioGrande",
    "GC1208": "Southwest",
    "GC1209": "Navajo8",
    "GC1210": "NavajoHighElev",
    "GC1698": "Foothills",
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
    "basin",
    [
        {
            "DRNAREA": 10.0,
            "BSLDEM10M": 10.0,
            "PRECIP": 20.0,
            "EL7500": 10.0,
            "STATSCLAY": 15.0,
            "I6H100Y": 2.6,
            "OUTLETELEV": 5000.0,
            "ELEV": 5000.0,
            "BSLDEM30ff": 0.05,
        },
        BASIN,
        {
            "DRNAREA": 500.0,
            "BSLDEM10M": 20.0,
            "PRECIP": 44.0,
            "EL7500": 95.0,
            "STATSCLAY": 35.0,
            "I6H100Y": 4.5,
            "OUTLETELEV": 8000.0,
            "ELEV": 7500.0,
            "BSLDEM30ff": 0.4,
        },
    ],
)
def test_matches_live_nss(lib, basin):
    """Every region x AEP against NSS region 9 (Colorado), citations 11, 20 and 21.

    NSS returns three significant figures, so the check is exact agreement after
    rounding flowfreq's value the same way.
    """
    from flowfreq.streamstats import Characteristic, estimate_flow_statistics

    chars = {k: Characteristic(k, k, "", v, "") for k, v in basin.items()}
    results, skipped = estimate_flow_statistics("CO", chars, ["PFS"])
    assert not skipped
    compared = 0
    for r in results:
        region = NSS_CODE_TO_REGION[r.region_code]
        for code, est in r.estimates.items():
            ours = evaluate(lib.equation(region, NSS_STAT_TO_AEP[code]), basin).flow_cfs
            assert _sig3(ours) == est.value, (region, code, est.value, ours)
            assert abs(ours - est.value) / est.value < 5e-3
            compared += 1
    assert compared == 62
