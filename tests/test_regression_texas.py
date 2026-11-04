"""Texas's verified equations: seven USGS reports behind NSS's ten regions.

WRI 98-4015 (Colorado small and large), WRI 98-4178 (Brazos 1-3), WRI 86-4069 (Austin urban,
which NSS miscites to 98-4015), WRI 82-18 (Dallas-Fort Worth urban), WRI 96-4072 (Highland
Lakes), WRI 80-17 (Houston urban) and SIR 2008-5084 Table 8 (statewide).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "flowfreq" / "data" / "regression"
GAGES = json.loads((ROOT / "tests/fixtures/regression/texas_gages.json").read_text("utf-8"))
NSS = {
    "Austin urban": "GC1489",
    "Brazos 1": "GC1490",
    "Brazos 2": "GC1491",
    "Brazos 3": "GC1492",
    "Colorado large": "GC1493",
    "Colorado small": "GC1494",
    "Dallas-Fort Worth urban": "GC1495",
    "Highland Lakes": "GC1496",
    "Houston urban": "GC1497",
    "Statewide": "GC2146",
}
AEPS = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01)
#: Sites the 98-4178 text excludes from the regressions (p. 11).
BRAZOS_EXCLUDED = {95, 99, 100, 103, 117, 120, 138, 153, 164}
BEDIAS_CREEK = 24


@pytest.fixture(scope="module")
def lib():
    return load_state("TX")


def _flow(eq, **chars):
    return evaluate(eq, chars, allow_extrapolation=True).flow_cfs


def _wls(rows, slope, index=5):
    """Refit one equation from Table 1, with weights normalized to sum to n as the reports do."""
    x = np.column_stack(
        [np.ones(len(rows)), np.log10([r["A"] for r in rows])]
        + ([np.log10([r["SL"] for r in rows])] if slope else [])
    )
    w = np.array([r["w"] for r in rows], dtype=float)
    w = w * len(w) / w.sum()
    u = np.linalg.inv(x.T @ (w[:, None] * x))
    y = np.log10([r["q"][index] for r in rows])
    b = u @ x.T @ (w * y)
    s = math.sqrt((w * (y - x @ b) ** 2).sum() / (len(rows) - x.shape[1]))
    return b, s, u, x


def _brazos(region, *, bedias_in_3=True):
    out = []
    for r in GAGES["brazos"]:
        if r["site"] in BRAZOS_EXCLUDED:
            continue
        reg = r["region"]
        if bedias_in_3 and r["site"] == BEDIAS_CREEK:
            reg = 3
        if reg == region:
            out.append(r)
    return out


def _printed_matrix(eq):
    """The report's unscaled matrix: the stored one is s^2 times it."""
    return np.array(eq.covariance) / eq.model_error_variance


def test_texas_is_verified(lib):
    assert lib.status == "verified"
    assert sorted(lib.regions) == sorted(NSS)
    assert len(lib.equations) == 65
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region


def test_every_equation_gives_the_nss_flow(lib):
    """No flow, AEP, exponent or limit differs. The (100 + x) forms fold the 100 into NSS's constant."""
    snap = json.loads((DATA / "nss_snapshot" / "TX.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert list(e.coefficients) == s["coefficients"]
        chars = {}
        for v, n in zip(e.variables, s["variables"]):
            assert v.code == n["code"]
            assert (v.minimum, v.maximum) == (n["minimum"], n["maximum"])
            chars[v.code] = math.sqrt(n["minimum"] * n["maximum"])
        log_q = s["intercept"]
        for n, c in zip(s["variables"], s["coefficients"]):
            x = chars[n["code"]]
            if n["transform"] == "power":
                log_q += c * x ** n["exponent"]
            else:
                log_q += c * math.log10(x + n.get("offset", 0.0))
        assert _flow(e, **chars) == pytest.approx(10**log_q, rel=1e-9)


def test_nss_citation_errors_are_corrected(lib):
    assert "86-4069" in lib.equation("Austin urban", 0.01).citation.publication
    assert "98-4015" in lib.equation("Colorado large", 0.01).citation.publication
    assert "Asquith" in lib.equation("Colorado large", 0.01).citation.publication
    assert "SIR 2009-5087" in lib.notes
    assert lib.notes.startswith("SUPERSEDED")


class TestColoradoExample:
    """WRI 98-4015 Appendix: 50 mi2, 10.5 ft/mi, the 100-year peak and its 67-percent interval."""

    CHARS = {"CONTDA": 50.0, "CSLBlue": 10.5}

    def test_small_and_large_estimates(self, lib):
        small = evaluate(
            lib.equation("Colorado small", 0.01),
            self.CHARS,
            allow_extrapolation=True,
            interval_level=0.67,
        )
        large = evaluate(lib.equation("Colorado large", 0.01), self.CHARS, interval_level=0.67)
        assert small.flow_cfs == pytest.approx(16800, rel=0.002)
        assert large.flow_cfs == pytest.approx(26600, rel=0.002)
        assert small.interval == pytest.approx((7710, 36600), rel=0.003)
        assert large.interval == pytest.approx((17100, 41400), rel=0.003)

    def test_leverage_uses_the_printed_matrix(self, lib):
        x = np.array([1.0, 1.6990, 1.0212])
        large = lib.equation("Colorado large", 0.01)
        small = lib.equation("Colorado small", 0.01)
        assert x @ _printed_matrix(large) @ x == pytest.approx(0.27388, abs=3e-5)
        assert x @ _printed_matrix(small) @ x == pytest.approx(0.32152, abs=3e-5)
        assert large.model_error_variance == pytest.approx(0.17**2)
        assert small.model_error_variance == pytest.approx(0.0841)

    def test_equation_3_weighting(self):
        """Stored in no field: equations 3-5 weight by log10 A from 10 to 100 mi2."""
        weight = 2 - math.log10(50)
        assert round(weight, 2) == 0.30  # p. I-1's prose prints the other weight as 0.070
        assert 0.30 * 16800 + 0.70 * 26600 == pytest.approx(23700, rel=0.002)
        assert 0.30 * 36600 + 0.70 * 41400 == pytest.approx(39960)
        assert 0.30 * 7710 + 0.70 * 17100 == pytest.approx(14300, rel=0.002)


class TestColoradoRefit:
    """Table 1's gage data reproduce Tables 3 and 4 and the leverages on p. 11."""

    @pytest.mark.parametrize(
        "region,small,mean_h,max_h",
        [("Colorado small", True, 0.24, 0.54), ("Colorado large", False, 0.14, 0.38)],
    )
    def test_matrix_and_leverage(self, lib, region, small, mean_h, max_h):
        rows = [
            r for r in GAGES["colorado"] if r["site"] not in (13, 16) and (r["A"] < 32) == small
        ]
        assert len(rows) == lib.equation(region, 0.01).n_sites
        _, _, u, x = _wls(rows, True)
        printed = _printed_matrix(lib.equation(region, 0.01))
        # The refit is close, not exact: the largest gap is 0.2 percent, cause not established.
        assert u == pytest.approx(printed, rel=0.0025)
        h = np.einsum("ij,jk,ik->i", x, printed, x)
        assert h.mean() == pytest.approx(mean_h, abs=0.005)
        assert h.max() == pytest.approx(max_h, abs=0.005)

    @pytest.mark.parametrize("region,small", [("Colorado small", True), ("Colorado large", False)])
    @pytest.mark.parametrize("i,aep", list(enumerate(AEPS + (0.004, 0.002))))
    def test_coefficients(self, lib, region, small, i, aep):
        rows = [
            r for r in GAGES["colorado"] if r["site"] not in (13, 16) and (r["A"] < 32) == small
        ]
        b, s, _, _ = _wls(rows, True, i)
        eq = lib.equation(region, aep)
        assert 10 ** b[0] == pytest.approx(10**eq.intercept, rel=0.003)
        assert list(b[1:]) == pytest.approx(list(eq.coefficients), abs=0.0015)
        assert s == pytest.approx(eq.model_error_variance**0.5, abs=0.0055)


class TestBrazos:
    def test_worked_example(self, lib):
        """WRI 98-4178 p. 19: region 3, 10 mi2, the 100-year peak, 90-percent interval."""
        eq = lib.equation("Brazos 3", 0.01)
        est = evaluate(eq, {"CONTDA": 10.0})
        assert est.flow_cfs == pytest.approx(2910, rel=0.002)
        x = np.array([1.0, 1.0])
        assert x @ _printed_matrix(eq) @ x == pytest.approx(0.29634, abs=1e-5)
        # The computation prints 9,670. The sentence after it prints 9,760, a transposition.
        assert est.interval == pytest.approx((876, 9670), rel=0.002)

    def test_table_1_miscounts_bedias_creek(self):
        """Table 1 puts it in region 2, giving 80 and 41 stations against Table 2's 79 and 42."""
        assert len(_brazos(2, bedias_in_3=False)) == 80
        assert len(_brazos(3, bedias_in_3=False)) == 41
        assert len(_brazos(1)) == 56
        assert len(_brazos(2)) == 79
        assert len(_brazos(3)) == 42

    @pytest.mark.parametrize("region", [1, 2])
    def test_matrices_of_regions_1_and_2(self, lib, region):
        """Within 0.3 percent, not exact (region 2 largest), and that fixes the printed order."""
        _, _, u, _ = _wls(_brazos(region), True)
        printed = _printed_matrix(lib.equation(f"Brazos {region}", 0.01))
        assert u == pytest.approx(printed, rel=0.003)

    @pytest.mark.parametrize("region", [1, 2, 3])
    @pytest.mark.parametrize("i,aep", list(enumerate(AEPS)))
    def test_refit_gives_the_printed_equations(self, lib, region, i, aep):
        b, s, _, _ = _wls(_brazos(region), region < 3, i)
        eq = lib.equation(f"Brazos {region}", aep)
        assert 10 ** b[0] == pytest.approx(10**eq.intercept, rel=0.01)
        assert list(b[1:]) == pytest.approx(list(eq.coefficients), abs=0.0025)
        assert s == pytest.approx(eq.model_error_variance**0.5, abs=0.0011)

    def test_region_3_matrix_is_from_a_fit_with_slope(self, lib):
        """Table 3's region 3 matrix does not belong to region 3's printed equation.

        It is the constant and CDA block of a three-term fit with stream slope. Stored as
        printed because the report's own example uses it, and flagged in the notes.
        """
        printed = _printed_matrix(lib.equation("Brazos 3", 0.01))
        _, _, with_slope, _ = _wls(_brazos(3), True)
        assert with_slope[:2, :2] == pytest.approx(printed, abs=6e-6)
        _, _, two_term, _ = _wls(_brazos(3), False)
        assert two_term == pytest.approx(
            np.array([[0.16508, -0.06928], [-0.06928, 0.03398]]), abs=6e-6
        )
        x = np.array([1.0, 1.0])
        assert x @ two_term @ x == pytest.approx(0.0605, abs=1e-4)
        assert "not the matrix of region 3's printed equation" in load_state("TX").notes


def test_austin_tables_reproduce(lib):
    """WRI 86-4069 p. 44's four basins, 45 and 0 percent, and Table 9's Waller Creek at 38th."""
    q5, q100 = lib.equation("Austin urban", 0.2), lib.equation("Austin urban", 0.01)
    for area, row in [
        (12.3, (5380, 2960, 14700, 8520)),
        (4.13, (2650, 1460, 7030, 4070)),
        (3.58, (2420, 1330, 6380, 3690)),
        (6.30, (3490, 1920, 9360, 5410)),
    ]:
        got = (
            _flow(q5, CONTDA=area, IMPERV=45),
            _flow(q5, CONTDA=area, IMPERV=0),
            _flow(q100, CONTDA=area, IMPERV=45),
            _flow(q100, CONTDA=area, IMPERV=0),
        )
        assert got == pytest.approx(row, rel=0.0025)
    for aep, earlier, later, ratio in zip(
        AEPS,
        (738, 1290, 1730, 2360, 2890, 3460),
        (989, 1660, 2200, 2980, 3640, 4360),
        (1.99, 1.82, 1.76, 1.73, 1.73, 1.73),
    ):
        eq = lib.equation("Austin urban", aep)
        assert _flow(eq, CONTDA=2.31, IMPERV=17) == pytest.approx(earlier, rel=0.0025)
        assert _flow(eq, CONTDA=2.31, IMPERV=37) == pytest.approx(later, rel=0.0025)
        assert 1.45 ** eq.coefficients[1] == pytest.approx(ratio, abs=0.0051)


def test_dallas_effects_table_reproduces(lib):
    """WRI 82-18 p. 53, at urbanization index 36 and 9, with one cell pinned."""
    q5 = lib.equation("Dallas-Fort Worth urban", 0.2)
    q100 = lib.equation("Dallas-Fort Worth urban", 0.01)
    for area, row in [
        (17.7, (9790, 3460, 19500, 9640)),
        (7.98, (5500, None, 10700, 5290)),
        (1.25, (1440, 508, 2660, 1310)),
        (52.8, (21600, 7630, 44500, 21900)),
    ]:
        got = (
            _flow(q5, DRNAREA=area, TXURBINDEX=36),
            _flow(q5, DRNAREA=area, TXURBINDEX=9),
            _flow(q100, DRNAREA=area, TXURBINDEX=36),
            _flow(q100, DRNAREA=area, TXURBINDEX=9),
        )
        for g, r in zip(got, row):
            if r is not None:
                assert g == pytest.approx(r, rel=0.0035)
    # Turtle Creek's 5-year minimum prints 1,960 where the equation gives 1,942.
    assert _flow(q5, DRNAREA=7.98, TXURBINDEX=9) == pytest.approx(1942, abs=1)


class TestHouston:
    @pytest.mark.parametrize("gage", GAGES["houston"], ids=lambda g: g["no"])
    def test_table_10_per_gage(self, lib, gage):
        """Table 10's rounding: nearest 5 below 1,000 cfs, else nearest 10. K has 2 figures."""
        for aep, q in zip(AEPS, gage["q"]):
            est = evaluate(
                lib.equation("Houston urban", aep),
                {"DRNAREA": gage["A"], "CONVEY": gage["K"], "URBAN": gage["Ao"]},
            )
            assert abs(est.flow_cfs - q) <= (3.0 if q < 1000 else max(6.0, 0.0012 * q))

    def test_form_of_the_table_9_term(self, lib):
        """Q2 = 2.028 A^0.383 [K(1.0 + 0.01 Ao)]^0.447, the example under Table 9."""
        eq = lib.equation("Houston urban", 0.5)
        got = _flow(eq, DRNAREA=8.81, CONVEY=1.7e5, URBAN=54.0)
        assert got == pytest.approx(2.028 * 8.81**0.383 * (1.7e5 * 1.54) ** 0.447)

    def test_500_year_is_stored_though_table_10_stops_at_100(self, lib):
        assert lib.equation("Houston urban", 0.002).coefficients == (0.451, 0.498, 0.498)

    def test_table_11_factors(self, lib):
        """Rural Ao is 20 percent (p. 26), and rectification multiplies K by 15."""
        for aep, urban, channel in zip(
            AEPS, (1.26, 1.27, 1.28, 1.28, 1.29, 1.29), (3.36, 3.55, 3.65, 3.74, 3.79, 3.82)
        ):
            b = lib.equation("Houston urban", aep).coefficients[1]
            assert (2.0 / 1.2) ** b == pytest.approx(urban, abs=0.0051)
            assert 15.0**b == pytest.approx(channel, abs=0.0051)


class TestStatewide:
    """SIR 2008-5084 Table 8, whose matrices are in Appendix 2's R output."""

    #: Appendix 2's Std. Error columns: constant, A^lambda, log10 P.
    STD = {
        0.5: (0.27228, 0.30614, 0.07183),
        0.2: (0.11743, 0.11665, 0.06857),
        0.1: (0.11028, 0.09729, 0.07091),
        0.04: (0.11304, 0.08916, 0.07593),
        0.02: (0.11793, 0.08744, 0.08037),
        0.01: (0.12386, 0.08748, 0.08509),
    }
    #: Appendix 2's full-precision estimates.
    EST = {
        0.5: (17.3586632, -16.5114820, 0.9428574),
        0.2: (8.0796013, -6.4032928, 0.6065365),
        0.1: (7.199766, -5.106641, 0.439738),
        0.04: (6.8485138, -4.3290741, 0.2729137),
        0.02: (6.7765427, -3.9910263, 0.1706351),
        0.01: (6.77645185, -3.75754823, 0.08211505),
    }

    @pytest.mark.parametrize("aep", AEPS)
    def test_matrix_gives_the_printed_standard_errors(self, lib, aep):
        eq = lib.equation("Statewide", aep)
        assert np.sqrt(np.diag(eq.covariance)) == pytest.approx(self.STD[aep], rel=2.5e-4)

    @pytest.mark.parametrize("aep", AEPS)
    def test_table_8_rounds_appendix_2(self, lib, aep):
        eq = lib.equation("Statewide", aep)
        # Table 8 prints four significant figures, so 17.36 and -16.51 have two decimals.
        tol = 5.1e-3 if aep == 0.5 else 6e-4
        assert eq.intercept == pytest.approx(self.EST[aep][0], abs=tol)
        assert eq.coefficients[0] == pytest.approx(self.EST[aep][1], abs=tol)
        assert eq.coefficients[1] == pytest.approx(self.EST[aep][2], abs=6e-5)

    def test_q100_error_settled_from_appendix_2(self, lib):
        """Table 8 prints 0.3548 for Q100. The R output prints 0.3545, and its own numbers agree."""
        eq = lib.equation("Statewide", 0.01)
        assert eq.model_error_variance == pytest.approx(0.3545**2)
        assert 0.12386 / 0.12205140**0.5 == pytest.approx(0.3545, abs=6e-5)

    def test_power_form(self, lib):
        eq = lib.equation("Statewide", 0.5)
        assert eq.variables[0].exponent == -0.0157
        got = _flow(eq, DRNAREA=100.0, PRECIPCENT=30.0)
        assert got == pytest.approx(10 ** (17.36 - 16.51 * 100.0**-0.0157) * 30.0**0.9429)


def test_highland_lakes_error_is_the_printed_percent(lib):
    """The sheet's 0.036 log10^2 mean square error at 100 years agrees with its 46 percent."""
    eq = lib.equation("Highland Lakes", 0.01)
    assert eq.covariance is None
    assert eq.sep_log**2 == pytest.approx(0.036, abs=5e-4)


def test_out_of_range_is_refused(lib):
    with pytest.raises(ValueError):
        evaluate(lib.equation("Houston urban", 0.01), {"DRNAREA": 500, "CONVEY": 1e5, "URBAN": 50})
