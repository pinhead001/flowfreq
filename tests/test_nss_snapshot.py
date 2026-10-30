"""Tests for flowfreq.regression.nss (NSS equation parser) and tools/snapshot_nss.py.

The offline tests use the real WA equation strings from
``tests/fixtures/streamstats_responses.py`` and synthetic strings for each form
the parser accepts or refuses. ``TestLiveSnapshot`` contacts NSS and is marked
``requires_network``.
"""

from __future__ import annotations

import copy
import importlib.util
import math
from pathlib import Path

import pytest

from flowfreq.regression.equations import Citation, OutOfRangeError, evaluate
from flowfreq.regression.nss import (
    NSSEquationError,
    ParsedConstant,
    evaluate_expression,
    parse_constant_equation,
    parse_equation,
    parse_expression,
    statistic_code_to_aep,
    to_regression_equation,
    variables_in,
)
from tests.fixtures.streamstats_responses import NSS_ESTIMATE_RESPONSE_PFS

REPO_ROOT = Path(__file__).resolve().parent.parent
CIT = Citation(publication="Synthetic test report", table="T1")


def _load_tool():
    spec = importlib.util.spec_from_file_location(
        "snapshot_nss", REPO_ROOT / "tools" / "snapshot_nss.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


tool = _load_tool()


def _roundtrip(equation: str, values: dict) -> float:
    """Relative difference between the parsed form and direct evaluation."""
    parsed = parse_equation(equation)
    direct = evaluate_expression(equation, values)
    return abs(10.0 ** parsed.log10_value(values) - direct) / direct


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


class TestParseWashingtonFixture:
    REGION = NSS_ESTIMATE_RESPONSE_PFS[0]["regressionRegions"][0]
    INPUTS = {p["code"]: p["value"] for p in REGION["parameters"]}

    @pytest.mark.parametrize("result", REGION["results"], ids=lambda r: r["code"])
    def test_reproduces_nss_value(self, result):
        parsed = parse_equation(result["equation"])
        eq = to_regression_equation(
            parsed,
            region_code=self.REGION["code"],
            aep=statistic_code_to_aep(result["code"]),
            citation=CIT,
        )
        flow = evaluate(eq, self.INPUTS).flow_cfs
        # Exact against the string itself...
        direct = evaluate_expression(result["equation"], self.INPUTS)
        assert flow == pytest.approx(direct, rel=1e-12)
        # ...and within NSS's three-significant-figure rounding of its own value.
        assert abs(flow - result["value"]) <= tool._rounding_bound(result["value"])

    def test_structure_of_pk50(self):
        parsed = parse_equation(self.REGION["results"][0]["equation"])
        assert parsed.codes == ("DRNAREA", "PRECPRIS10", "CANOPY_PCT")
        assert parsed.transforms == ("log10", "identity", "identity")
        assert parsed.coefficients == pytest.approx((0.745, 0.032, -0.0078))
        assert parsed.intercept == pytest.approx(math.log10(3.846))


class TestParseForms:
    VALUES = {"A": 37.0, "B": 2350.0, "C": 12.5, "D": 0.4}

    def test_scaled_power_term_folds_into_intercept(self):
        parsed = parse_equation("27.9*A^0.523*(B/1000)^(-0.098)")
        assert parsed.transforms == ("log10", "log10")
        assert parsed.intercept == pytest.approx(math.log10(27.9) + 0.098 * 3.0)
        assert any("scale" in n for n in parsed.notes)
        assert _roundtrip("27.9*A^0.523*(B/1000)^(-0.098)", self.VALUES) < 1e-12

    def test_multiplicative_scale_parenthesized(self):
        s = "0.00953*(A)^(0.392)*(0.001*(B))^(3.36)"
        parsed = parse_equation(s)
        assert parsed.intercept == pytest.approx(math.log10(0.00953) + 3.36 * math.log10(0.001))
        assert _roundtrip(s, self.VALUES) < 1e-12

    def test_plus_one(self):
        s = "5.0*A^0.674*(C+1)^(-0.026)"
        parsed = parse_equation(s)
        assert parsed.transforms == ("log10", "log10_plus1")
        assert parsed.coefficients == pytest.approx((0.674, -0.026))
        assert _roundtrip(s, self.VALUES) < 1e-12

    def test_scaled_plus_one(self):
        s = "0.0000173*(A)^(0.957)*(0.001*(B+1))^(-2.19)"
        assert parse_equation(s).transforms == ("log10", "log10_plus1")
        assert _roundtrip(s, self.VALUES) < 1e-12

    def test_division_by_power_of_ten(self):
        s = "3.846*A^0.745*10^(0.032*C)/10^(0.0078*D)"
        parsed = parse_equation(s)
        assert parsed.transforms == ("log10", "identity", "identity")
        assert parsed.coefficients == pytest.approx((0.745, 0.032, -0.0078))
        assert _roundtrip(s, self.VALUES) < 1e-12

    def test_constant_power_of_ten_and_max_zero(self):
        s = "max(0,1.0176*(10^(-3.0008)*A^1.0485*C^1.7013))"
        parsed = parse_equation(s)
        assert parsed.intercept == pytest.approx(math.log10(1.0176) - 3.0008)
        assert any("max(0" in n for n in parsed.notes)
        assert _roundtrip(s, self.VALUES) < 1e-12

    def test_bare_variable_and_repeated_variable(self):
        parsed = parse_equation("0.145*(A)*(A)^(0.5)")
        assert parsed.codes == ("A",)
        assert parsed.coefficients == pytest.approx((1.5,))

    def test_scientific_notation(self):
        assert parse_equation("1.5E-3*A^2").intercept == pytest.approx(math.log10(1.5e-3))

    @pytest.mark.parametrize(
        "equation",
        [
            "2*(A+C)^0.5",  # affine base in two variables
            "2*(A*C+1)^0.5",  # base not affine
            "2*A^0.5*(A+1)^0.2",  # one variable, two rescalings
            "max(0,1.53*(10^(-4.9)*A^1.09*C^2.6-1))",  # subtracted constant
            "10^(0.1*A)*A^0.5",  # one variable, two transforms
            "A^(C^2)",  # exponent not affine in the variables
            "A^(log(C))",  # exponent not affine in the variables
            "(A^C)^(0.5+0.1*D)",  # an interaction raised to an affine exponent
            "e#^(8.92+(-0.0141*(A)))/(1+e#^(8.92+(-0.0141*(A))))",  # logistic, malformed
            "exp(A)",  # unsupported function
            "-2*A^0.5",  # negative
            "3*",  # malformed
            "",
        ],
    )
    def test_refuses_what_it_cannot_represent(self, equation):
        with pytest.raises(NSSEquationError):
            parse_equation(equation)

    def test_affine_exponent_becomes_interaction_terms(self):
        # The SIR 2023-5006 Southeast US form (NSS GC1934/GC1935/GC1936, #211).
        s = (
            "(LC19IMP^0)*10^(1.82+0.00354*PCTREG1-0.00202*PCTREG3+0.00187*PCTREG5)"
            "*DRNAREA^(0.646+0.00224*PCTREG2)*PCTREG4^(0)"
        )
        parsed = parse_equation(s)
        terms = dict(zip(zip(parsed.codes, parsed.times), parsed.coefficients))
        assert terms == pytest.approx(
            {
                ("LC19IMP", None): 0.0,
                ("PCTREG1", None): 0.00354,
                ("PCTREG3", None): -0.00202,
                ("PCTREG5", None): 0.00187,
                ("DRNAREA", None): 0.646,
                ("DRNAREA", "PCTREG2"): 0.00224,
                ("PCTREG4", None): 0.0,
            }
        )
        assert parsed.intercept == pytest.approx(1.82)
        values = {"DRNAREA": 26.686326086593485, "LC19IMP": 5.0}
        values.update({f"PCTREG{i}": 50.0 for i in range(1, 6)})
        assert _roundtrip(s, values) < 1e-12
        # NSS's own value at those inputs, from the 2026-10-25 GA snapshot (3 s.f.).
        assert 10.0 ** parsed.log10_value(values) == pytest.approx(1180.0, rel=5e-3)

    def test_affine_exponent_on_a_scaled_base(self):
        # A constant inside the base: (3*A)^(0.5 + 0.1*C) adds 0.1*log10(3)*C.
        s = "2*(3*A)^(0.5+0.1*C)"
        parsed = parse_equation(s)
        assert parsed.times == (None, None, "C")
        assert _roundtrip(s, self.VALUES) < 1e-12
        assert parse_equation("A^C").times == (None, "C")

    @pytest.mark.parametrize(
        "equation, code, scale, offset",
        [
            ("79*(A)^(0.93)*(C+0.1)^(0.05)", "C", 1.0, 0.1),
            ("0.000592*A^0.981*(C/100+1)^(-1.52)", "C", 0.01, 1.0),
            ("0.000592*A^0.981*(1+C/100)^(-1.52)", "C", 0.01, 1.0),
            ("max(0,1.04*(10^(-2.1)*A^1.1*(C-10)^1.4))", "C", 1.0, -10.0),
            ("3*(0.5*B+100)^0.25", "B", 0.5, 100.0),
        ],
    )
    def test_affine_base_emits_scale_and_offset(self, equation, code, scale, offset):
        parsed = parse_equation(equation)
        i = parsed.codes.index(code)
        assert parsed.transforms[i] == "log10"
        assert parsed.scales[i] == pytest.approx(scale, rel=1e-15)
        assert parsed.offsets[i] == pytest.approx(offset, rel=1e-15)
        assert _roundtrip(equation, self.VALUES) < 1e-12
        var = to_regression_equation(parsed, region_code="R", aep=0.01, citation=CIT).variables[i]
        assert (var.scale, var.offset) == (parsed.scales[i], parsed.offsets[i])
        d = parsed.to_dict()["variables"][i]
        assert d["scale" if scale != 1.0 else "offset"] == pytest.approx(
            scale if scale != 1.0 else offset
        )

    def test_plus_one_keeps_log10_plus1(self):
        parsed = parse_equation("5.0*A^0.674*(C+1)^(-0.026)")
        assert parsed.scales == (1.0, 1.0) and parsed.offsets == (0.0, 0.0)
        assert parsed.to_dict()["variables"] == [
            {"code": "A", "transform": "log10"},
            {"code": "C", "transform": "log10_plus1"},
        ]

    def test_power_inside_exponent_of_10(self):
        # Arizona region 5's form, 10^(a - b*X^(-c)): power with exponent -c and
        # coefficient -b, a folded into b0.
        s = "10^(6.363-4.386*CONTDA^(-0.060))"
        parsed = parse_equation(s)
        assert parsed.codes == ("CONTDA",) and parsed.transforms == ("power",)
        assert parsed.exponents == (-0.06,)
        assert parsed.coefficients == pytest.approx((-4.386,), rel=1e-15)
        assert parsed.intercept == pytest.approx(6.363, rel=1e-15)
        assert _roundtrip(s, {"CONTDA": 37.0}) < 1e-12
        var = to_regression_equation(parsed, region_code="5", aep=0.5, citation=CIT).variables[0]
        assert (var.transform, var.exponent) == ("power", -0.06)
        assert parsed.to_dict()["variables"] == [
            {"code": "CONTDA", "transform": "power", "exponent": -0.06}
        ]

    def test_power_one_inside_exponent_stays_identity(self):
        parsed = parse_equation("10^(0.5*A^1)")
        assert parsed.transforms == ("identity",) and parsed.exponents == (None,)

    @pytest.mark.parametrize(
        "equation",
        [
            "10^(A^0.5+A)",  # one variable, two powers
            "10^(1-A^0.5*C)",  # product of variables
            "10^(A^C)",  # non-constant power
            "2*(A^0.5+1)^0.3",  # power inside an affine base
        ],
    )
    def test_power_refusals(self, equation):
        with pytest.raises(NSSEquationError):
            parse_equation(equation)

    def test_committed_snapshot_power_refusals_now_parse(self):
        """The committed AZ snapshot refused region 5's eight strings ("exponent of
        10 is not linear"); each now parses to ``power`` and reproduces direct
        evaluation at the snapshot's submitted inputs."""
        import json

        checked = 0
        for path in sorted((REPO_ROOT / "data" / "nss_snapshots").glob("*.json")):
            snap = json.loads(path.read_text(encoding="utf-8"))
            for region in snap["regions"]:
                inputs = {p["code"]: p["submitted"] for p in region.get("parameters", [])}
                for rec in region.get("equations", []):
                    if not (rec.get("parse_error") or "").startswith("exponent of 10 is not"):
                        continue
                    parsed = parse_equation(rec["equation"])
                    assert parsed.transforms == ("power",)
                    direct = evaluate_expression(rec["equation"], inputs)
                    assert 10.0 ** parsed.log10_value(inputs) == pytest.approx(direct, rel=1e-12)
                    checked += 1
        assert checked == 8

    def test_idaho_region4_matches_offline_library_form(self):
        # NSS's string for Idaho region 4 Q80 (snapshot 2026-09-26) parses to the
        # published form: log10 with scale 0.01 and offset 1 on LC11FOREST.
        s = "0.000592*DRNAREA^0.981*(LC11FOREST/100+1)^(-1.52)*PRECPRIS10^2.80"
        parsed = parse_equation(s)
        assert parsed.codes == ("DRNAREA", "LC11FOREST", "PRECPRIS10")
        assert parsed.scales == (1.0, 0.01, 1.0) and parsed.offsets == (0.0, 1.0, 0.0)
        assert parsed.intercept == pytest.approx(math.log10(0.000592), rel=1e-15)

    def test_evaluate_expression_handles_refused_forms(self):
        s = "max(0,1.53*(10^(-4.9)*A^1.09*C^2.6-1))"
        with pytest.raises(NSSEquationError):
            parse_equation(s)
        expected = max(0.0, 1.53 * (10**-4.9 * 37.0**1.09 * 12.5**2.6 - 1))
        assert evaluate_expression(s, self.VALUES) == pytest.approx(expected)

    def test_evaluate_expression_missing_variable(self):
        with pytest.raises(KeyError):
            evaluate_expression("2*A^0.5", {})

    def test_power_is_right_associative_and_binds_tighter_than_minus(self):
        assert evaluate_expression("2^3^2", {}) == 512.0
        assert evaluate_expression("-2^2", {}) == -4.0
        assert parse_expression("10^-1") == ("bin", "^", ("num", 10.0), ("neg", ("num", 1.0)))

    def test_offset_outside_domain_raises(self):
        parsed = parse_equation("2*(C-20)^1.4")
        with pytest.raises(ValueError, match="needs C - 20 > 0"):
            parsed.log10_value({"C": 12.5})

    def test_committed_snapshot_refusals_now_parse(self):
        """Every equation the committed snapshots refused for an additive offset
        now parses and reproduces direct evaluation of the NSS string at the
        snapshot's submitted inputs."""
        import json

        checked = 0
        for path in sorted((REPO_ROOT / "data" / "nss_snapshots").glob("*.json")):
            snap = json.loads(path.read_text(encoding="utf-8"))
            for region in snap["regions"]:
                inputs = {p["code"]: p["submitted"] for p in region.get("parameters", [])}
                for rec in region.get("equations", []):
                    if not (rec.get("parse_error") or "").startswith("additive offset"):
                        continue
                    direct = evaluate_expression(rec["equation"], inputs)
                    value = 10.0 ** parse_equation(rec["equation"]).log10_value(inputs)
                    assert value == pytest.approx(direct, rel=1e-12), rec["equation"]
                    checked += 1
        assert checked == 17  # 11 Idaho region 4, 6 Portland urban (OR)

    def test_variables_in(self):
        assert variables_in("3*A^2*(B+1)^0.5*A") == ["A", "B"]


# NSS codes with lower-case letters (AR, TX, MS, OK, MD).
_MIXED_CASE_INPUTS = {
    "DRNAREA": 50.0,
    "CSLBlue": 20.0,
    "BSHAPELFP": 3.0,
    "CSL10_85fm": 20.0,
    "IMPERV": 5.0,
    "SOILCorD": 30.0,
}


@pytest.mark.parametrize(
    "equation, codes",
    [
        (
            "10^2.183*DRNAREA^0.678*CSLBlue^0.239*BSHAPELFP^(-0.324)",
            ("DRNAREA", "CSLBlue", "BSHAPELFP"),
        ),
        (
            "207*DRNAREA^0.676*CSL10_85fm^0.159*BSHAPELFP^-0.390",
            ("DRNAREA", "CSL10_85fm", "BSHAPELFP"),
        ),
        (
            "5.18*DRNAREA^0.694*(IMPERV+1)^(0.382)*(SOILCorD+1)^(0.414)",
            ("DRNAREA", "IMPERV", "SOILCorD"),
        ),
    ],
)
def test_mixed_case_identifiers_parse(equation, codes):
    parsed = parse_equation(equation)
    assert parsed.codes == codes
    assert parsed.transforms[1] in ("log10", "log10_plus1")
    direct = evaluate_expression(equation, _MIXED_CASE_INPUTS)
    assert 10.0 ** parsed.log10_value(_MIXED_CASE_INPUTS) == pytest.approx(direct, rel=1e-12)


@pytest.mark.parametrize("template", ["2*DRNAREA^0.5*{}^0.3", "2*DRNAREA^0.5*{}"])
@pytest.mark.parametrize("name", ["max", "logN", "LogN", "e", "exp"])
def test_mixed_case_function_name_as_variable_refused(name, template):
    with pytest.raises(NSSEquationError, match="unsupported identifier"):
        parse_equation(template.format(name))


# NSS's logN(X, b) inside 10^(...) (88 equations in MN, WI and MA).
_LOGN_MN = "10^(24.074+0.744*logN(DRNAREA,10)-0.297*logN((LAKEAREA+1),10)-13.720*logN(LAT_OUT,10))"
_LOGN_WI = (
    "10^(0.9357+0.8398*logN(DRNAREA,10)-0.6685*logN(SSURGOKSAT,10)"
    "-0.3740*logN((LC01HERB+1),10)-0.5143*logN((LC01WATER+1),10)+0.4227*I24H10Y)"
)
_LOGN_MA = "10^(1.631+0.801*LogN(DRNAREA,10)+0.589*0.001*(ELEV/3.28084)-1.137*0.01*(LC06STOR))"


def test_logN_terms_parse_to_log10():
    p = parse_equation(_LOGN_MN)
    assert p.intercept == pytest.approx(24.074)
    assert p.codes == ("DRNAREA", "LAKEAREA", "LAT_OUT")
    assert p.transforms == ("log10", "log10_plus1", "log10")
    assert p.coefficients == pytest.approx((0.744, -0.297, -13.720))


def test_identity_terms_with_folded_constants():
    p = parse_equation(_LOGN_MA)
    i = p.codes.index("ELEV")
    assert p.transforms[i] == "identity"
    assert p.coefficients[i] == pytest.approx(0.589 * 0.001 / 3.28084)
    assert any("folded" in n for n in p.notes)


def test_logN_wisconsin_round_trip():
    values = {
        "DRNAREA": 120.0,
        "SSURGOKSAT": 15.0,
        "LC01HERB": 30.0,
        "LC01WATER": 2.5,
        "I24H10Y": 4.1,
    }
    parsed = parse_equation(_LOGN_WI)
    assert parsed.transforms == ("log10", "log10", "log10_plus1", "log10_plus1", "identity")
    direct = evaluate_expression(_LOGN_WI, values)
    assert 10.0 ** parsed.log10_value(values) == pytest.approx(direct, rel=1e-12)


def test_logN_other_base_scales_coefficient():
    s = "10^(1.5+0.7*logN(DRNAREA,2))"
    parsed = parse_equation(s)
    assert parsed.transforms == ("log10",)
    assert parsed.coefficients == pytest.approx((0.7 / math.log10(2.0),))
    assert _roundtrip(s, {"DRNAREA": 37.0}) < 1e-12


def test_logN_affine_base_keeps_scale_and_offset():
    s = "10^(0.5+0.3*logN((2*DRNAREA+3),10))"
    parsed = parse_equation(s)
    assert (parsed.transforms, parsed.scales, parsed.offsets) == (("log10",), (2.0,), (3.0,))
    assert _roundtrip(s, {"DRNAREA": 37.0}) < 1e-12


@pytest.mark.parametrize(
    "equation, match",
    [
        ("10^(1+logN(DRNAREA))", "positive constant base"),  # one argument
        ("10^(1+logN(DRNAREA,ELEV))", "positive constant base"),  # variable base
        ("10^(1+logN(DRNAREA,-10))", "positive constant base"),
        ("10^(1+logN(DRNAREA,1))", "positive constant base"),
        ("10^(1+0.5*logN(DRNAREA,10)+0.1*DRNAREA)", "two different forms"),
        ("10^(1+0.5*logN((A*C+1),10))", "not affine"),
    ],
)
def test_logN_unrepresentable_forms_refused(equation, match):
    with pytest.raises(NSSEquationError, match=match):
        parse_equation(equation)


# A negated affine base, (-1*(X-c))^p (15 NSS equations in FL, KY and OH).
_NEG_FL = "214*(DRNAREA)^(1.13)*(-1*(BDF-13))^(-0.59)*(CSL10_85)^(0.73)"


def test_negated_affine_base():
    p = parse_equation(_NEG_FL)
    i = p.codes.index("BDF")
    assert p.transforms[i] == "log10"
    assert (p.scales[i], p.offsets[i]) == (-1.0, 13.0)
    assert p.coefficients[i] == pytest.approx(-0.59)


def test_negated_base_round_trip():
    inputs = {"DRNAREA": 10.0, "BDF": 4.0, "CSL10_85": 20.0}
    assert 10 ** parse_equation(_NEG_FL).log10_value(inputs) == pytest.approx(
        evaluate_expression(_NEG_FL, inputs), rel=1e-12
    )


@pytest.mark.parametrize(
    "base", ["(-1*(BDF-13))", "(13-BDF)", "((BDF-13)*-1)", "(-(BDF-13))", "((BDF-13)/-1)"]
)
def test_negated_base_spellings_agree(base):
    p = parse_equation(f"2*{base}^(-0.59)")
    assert (p.transforms, p.scales, p.offsets) == (("log10",), (-1.0,), (13.0,))
    assert p.coefficients == pytest.approx((-0.59,))


def test_negated_base_out_of_domain_raises_naming_the_code():
    with pytest.raises(ValueError, match="BDF: log10 argument -1 is not positive"):
        parse_equation(_NEG_FL).log10_value({"DRNAREA": 10.0, "BDF": 14.0, "CSL10_85": 20.0})


@pytest.mark.parametrize(
    "equation, match",
    [
        ("2*(-1*(A*C))^0.5", "negated base that is not affine"),
        ("2*(A-A+13)^0.5", "zero scale"),
    ],
)
def test_unrepresentable_affine_bases_refused(equation, match):
    with pytest.raises(NSSEquationError, match=match):
        parse_equation(equation)


class TestToRegressionEquation:
    def test_limits_and_units_carried(self):
        parsed = parse_equation("2*A^0.5*10^(0.01*C)")
        eq = to_regression_equation(
            parsed,
            region_code="GC1",
            aep=0.01,
            citation=CIT,
            limits={"A": (1.0, 100.0)},
            units={"A": "mi^2"},
        )
        assert eq.variables[0].minimum == 1.0 and eq.variables[0].units == "mi^2"
        assert eq.variables[1].minimum is None
        assert eq.sep_log is None

    def test_bad_aep_raises(self):
        with pytest.raises(ValueError):
            to_regression_equation(parse_equation("2*A"), region_code="R", aep=1.5, citation=CIT)


NV_Q0 = "0*(0*(DRNAREA))^(0)*(0*(ELEV))^(0)"  # NSS GC164 PK50AEP (snapshot 2026-10-01)


class TestParseConstant:
    """NSS's form of a printed ``Q=0``: a product with a literal zero factor."""

    def test_nevada_region_6_two_year(self):
        with pytest.raises(NSSEquationError, match="non-positive factor"):
            parse_equation(NV_Q0)
        parsed = parse_constant_equation(NV_Q0)
        assert parsed == ParsedConstant(equation=NV_Q0, value=0.0, codes=("DRNAREA", "ELEV"))
        assert evaluate_expression(NV_Q0, {"DRNAREA": 50.0, "ELEV": 6500.0}) == 0.0
        assert parsed.to_dict()["form"] == "constant"

    def test_to_regression_equation_is_range_checked(self):
        eq = to_regression_equation(
            parse_constant_equation(NV_Q0),
            region_code="GC164",
            aep=0.5,
            citation=CIT,
            limits={"DRNAREA": (0.2, 210.0), "ELEV": (4770.0, 9960.0)},
            units={"DRNAREA": "mi2", "ELEV": "ft"},
        )
        assert (eq.form, eq.constant_cfs, eq.intercept, eq.coefficients) == (
            "constant",
            0.0,
            None,
            (0.0, 0.0),
        )
        assert evaluate(eq, {"DRNAREA": 50.0, "ELEV": 6500.0}).flow_cfs == 0.0
        with pytest.raises(OutOfRangeError):
            evaluate(eq, {"DRNAREA": 500.0, "ELEV": 6500.0})

    @pytest.mark.parametrize(
        "equation",
        ["0", "2*DRNAREA^0.5", "0/DRNAREA", "0+0*DRNAREA", "0.5*0.1*DRNAREA", "0*bad$"],
    )
    def test_refuses_anything_else(self, equation):
        with pytest.raises(NSSEquationError):
            parse_constant_equation(equation)

    def test_tool_records_the_constant(self):
        stat = {"code": "PK50AEP", "name": "50 percent AEP flood", "equation": NV_Q0, "value": 0}
        region = {
            "code": "GC164",
            "name": "Region 6",
            "parameters": [
                {"code": "DRNAREA", "limits": {"min": 0.2, "max": 210.0}},
                {"code": "ELEV", "limits": {"min": 4770.0, "max": 9960.0}},
            ],
        }
        rec = tool._equation_record(stat, region, {"DRNAREA": 50.0, "ELEV": 6500.0}, CIT)
        assert rec["parse_error"] is None
        assert rec["parsed"]["form"] == "constant" and rec["parsed"]["constant_cfs"] == 0.0
        assert rec["round_trip"]["ok"] and rec["round_trip"]["evaluated"] == 0.0

    def test_committed_snapshot_constant_refusal_now_parses(self):
        import json

        snap = json.loads(
            (REPO_ROOT / "data" / "nss_snapshots" / "NV_2026-10-01.json").read_text("utf-8")
        )
        refused = [
            rec
            for region in snap["regions"]
            for rec in region["equations"]
            if (rec.get("parse_error") or "").startswith("non-positive factor")
        ]
        assert [r["equation"] for r in refused] == [NV_Q0]
        assert parse_constant_equation(NV_Q0).value == refused[0]["nss_value"]


class TestStatisticCodeToAep:
    @pytest.mark.parametrize(
        "code, aep",
        [
            ("PK50AEP", 0.5),
            ("PK66_7AEP", 0.667),
            ("PK0_2AEP", 0.002),
            ("PK42_9AEP", 0.429),
            ("PK1AEP", 0.01),
            ("ACPK66_7AE", 0.667),  # MT prefix, truncated to 10 characters
            ("BWPK0_5AEP", 0.005),
            ("pk4aep", 0.04),
            ("PK100", 0.01),
            ("PK1_5", 1 / 1.5),
        ],
    )
    def test_codes(self, code, aep):
        assert statistic_code_to_aep(code) == pytest.approx(aep)

    @pytest.mark.parametrize("code", ["M7D10Y", "QA", "PK1", "PK100AEP", "D50", ""])
    def test_rejects(self, code):
        with pytest.raises(ValueError):
            statistic_code_to_aep(code)


# ---------------------------------------------------------------------------
# Tool (offline parts)
# ---------------------------------------------------------------------------


class TestToolHelpers:
    def test_synthetic_value(self):
        assert tool.synthetic_value({"min": 1.0, "max": 100.0}) == pytest.approx(10.0)
        assert tool.synthetic_value({"min": 0.0, "max": 77.4}) == pytest.approx(38.7)
        assert tool.synthetic_value({"min": -10.0, "max": 10.0}) == 0.0
        assert tool.synthetic_value({"min": 2.0}) == 4.0
        assert tool.synthetic_value(None) is None
        assert tool.synthetic_value({}) is None

    def test_gating_candidates_use_gc_number(self):
        region = {
            "code": "GC729",
            "parameters": [
                {"code": "DRNAREA", "limits": {"min": 1, "max": 10}},
                {"code": "ORREG2", "limits": None},
            ],
        }
        cands = tool._gating_candidates(region)
        assert cands[0] == {"ORREG2": 729.0}
        assert {"ORREG2": tool.UNBOUNDED_DEFAULT} not in cands
        assert tool._gating_candidates({"code": "GC1", "parameters": []}) == []

    def test_gating_candidates_include_a_bounded_region_id_parameter(self):
        region = {
            "code": "GC1005",
            "parameters": [
                {"code": "HIGHREG", "limits": {"min": 1005, "max": 1012}},
                {"code": "CONTDA", "limits": {"min": 0.27, "max": 13706}},
            ],
        }
        cands = tool._gating_candidates(region)
        assert cands[0] == {"HIGHREG": 1005.0}
        assert tool._fill(region, cands[0])["parameters"][0]["value"] == 1005.0

    def test_gating_candidates_try_unbounded_before_region_id(self):
        region = {
            "code": "GC729",
            "parameters": [
                {"code": "REGID", "limits": {"min": 700, "max": 780}},
                {"code": "ORREG2", "limits": None},
            ],
        }
        cands = tool._gating_candidates(region)
        assert cands[0] == {"ORREG2": 729.0}
        with_id = [i for i, c in enumerate(cands) if "REGID" in c]
        without = [i for i, c in enumerate(cands) if "REGID" not in c]
        assert with_id and max(without) < min(with_id)
        assert cands[with_id[0]] == {"REGID": 729.0, "ORREG2": 729.0}
        assert {"ORREG2": tool.UNBOUNDED_DEFAULT} not in cands

    def test_resolved_unbounded_gate_is_tried_before_bounded_parameter(self):
        # New Mexico GC1092: HIGHREG is unbounded, and a narrow integral range
        # containing 1092 must not pre-empt the gate the region resolves with.
        region = {
            "code": "GC1092",
            "parameters": [
                {"code": "SUBREG", "limits": {"min": 1090, "max": 1100}},
                {"code": "HIGHREG", "limits": None},
            ],
        }
        cands = tool._gating_candidates(region)
        assert cands[0] == {"HIGHREG": 1092.0}
        first_bounded = next(i for i, c in enumerate(cands) if "SUBREG" in c)
        assert cands.index({"HIGHREG": 1092.0}) < first_bounded

    @pytest.mark.parametrize(
        "limits",
        [{"min": 0.536, "max": 10800}, {"min": 1, "max": 13706}, {"min": 1000.5, "max": 1100}],
    )
    def test_basin_characteristic_range_is_not_a_region_id(self, limits):
        region = {
            "code": "GC1092",
            "parameters": [
                {"code": "DRNAREA", "limits": limits},
                {"code": "HIGHREG", "limits": None},
            ],
        }
        cands = tool._gating_candidates(region)
        assert all("DRNAREA" not in c for c in cands)
        all_bounded = {"code": "GC1092", "parameters": [{"code": "DRNAREA", "limits": limits}]}
        assert tool._gating_candidates(all_bounded) == []

    def test_gating_candidates_ignore_bounds_not_containing_gc_number(self):
        region = {"code": "GC5", "parameters": [{"code": "A", "limits": {"min": 10, "max": 20}}]}
        assert tool._gating_candidates(region) == []

    def test_rounding_bound(self):
        assert tool._rounding_bound(4370.0) == pytest.approx(5.0)
        assert tool._rounding_bound(0.0123) == pytest.approx(0.00005)
        assert tool._rounding_bound(0.0) == 0.0

    def test_unresolved_region_is_recorded_not_dropped(self):
        group = {"id": 2, "code": "PFS"}
        tmpl = {
            "id": 479,
            "code": "GC730",
            "name": "Reg_2A",
            "citationID": 118,
            "parameters": [{"code": "ORREG2", "limits": None}],
        }
        rec = tool._region_record(group, tmpl, None, tool._fill(tmpl, {}), {}, [{"X": 1}], {})
        assert rec["status"] == "unresolved"
        assert "omitted" in rec["unresolved_reason"]
        assert rec["equations"] == []

    def test_region_record_round_trips_fixture(self):
        group = {"id": 2, "code": "PFS"}
        answered = NSS_ESTIMATE_RESPONSE_PFS[0]["regressionRegions"][0]
        tmpl = copy.deepcopy(answered)
        del tmpl["results"]
        rec = tool._region_record(group, tmpl, answered, answered, {}, [], {})
        assert rec["status"] == "resolved"
        for e in rec["equations"]:
            assert e["round_trip"]["ok"], e
            assert e["aep"] is not None
            assert e["nss_sep"] is not None and "sep_log" not in e

    def test_diff_reports_changes(self):
        def snap(eq: str, lim: float, extra: bool = False):
            regions = [
                {
                    "statistic_group_code": "PFS",
                    "region_code": "GC1",
                    "region_name": "R1",
                    "status": "resolved",
                    "parameters": [{"code": "A", "min": 1.0, "max": lim, "units": "mi^2"}],
                    "equations": [{"statistic_code": "PK1AEP", "equation": eq, "asep_pct": 50}],
                }
            ]
            if extra:
                regions.append({**regions[0], "region_code": "GC2", "equations": []})
            return {"regions": regions}

        old, new = snap("2*A^0.5", 100.0), snap("2*A^0.6", 200.0, extra=True)
        changes = tool.diff_snapshots(old, new)
        assert any("PK1AEP: 2*A^0.5 -> 2*A^0.6" in c for c in changes)
        assert any("parameter A" in c for c in changes)
        assert any(c.startswith("+ region PFS:GC2") for c in changes)
        assert tool.diff_snapshots(old, old) == []

    def test_previous_snapshot(self, tmp_path):
        assert tool.previous_snapshot(tmp_path, "WA") is None
        for name in ("WA_2026-01-01.json", "WA_2026-03-01.json", "ID_2026-05-01.json"):
            (tmp_path / name).write_text("{}")
        assert tool.previous_snapshot(tmp_path, "wa").name == "WA_2026-03-01.json"


# ---------------------------------------------------------------------------
# Live
# ---------------------------------------------------------------------------


@pytest.mark.requires_network
class TestLiveSnapshot:
    """Snapshot WA from live NSS. Deselect with -m 'not requires_network'."""

    def test_wa_snapshot_round_trips(self):
        snap = tool.snapshot_state("WA", ["PFS"])
        assert snap["summary"]["regions_unresolved"] == 0
        region = next(r for r in snap["regions"] if r["region_code"] == "GC1750")
        assert region["status"] == "resolved"
        inputs = {p["code"]: p["submitted"] for p in region["parameters"]}
        pk1 = next(e for e in region["equations"] if e["statistic_code"] == "PK1AEP")
        parsed = parse_equation(pk1["equation"])
        eq = to_regression_equation(parsed, region_code="GC1750", aep=0.01, citation=CIT)
        flow = evaluate(eq, inputs).flow_cfs
        assert flow == pytest.approx(evaluate_expression(pk1["equation"], inputs), rel=1e-12)
        assert abs(flow - pk1["nss_value"]) <= tool._rounding_bound(pk1["nss_value"])
        assert all(e["round_trip"]["ok"] for r in snap["regions"] for e in r["equations"])
