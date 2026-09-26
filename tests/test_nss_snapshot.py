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

from flowfreq.regression.equations import Citation, evaluate
from flowfreq.regression.nss import (
    NSSEquationError,
    evaluate_expression,
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
            "79*(A)^(0.93)*(C+0.1)^(0.05)",  # offset other than +1
            "0.000592*A^0.981*(C/100+1)^(-1.52)",  # scaled then +1
            "max(0,1.04*(10^(-2.1)*A^1.1*(C-20)^1.4))",  # subtraction in base
            "max(0,1.53*(10^(-4.9)*A^1.09*C^2.6-1))",  # subtracted constant
            "10^(0.1*A)*A^0.5",  # one variable, two transforms
            "A^C",  # non-constant exponent
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

    def test_variables_in(self):
        assert variables_in("3*A^2*(B+1)^0.5*A") == ["A", "B"]


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
