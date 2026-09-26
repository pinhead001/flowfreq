"""Tests for flowfreq.regression: schema, evaluation, library, jurisdictions.

Every equation here is synthetic. No published coefficients are used or implied.
"""

from __future__ import annotations

import json
import math

import pytest

from flowfreq.regression import (
    JURISDICTIONS,
    WAVE_NAMES,
    Citation,
    EquationsUnavailable,
    OutOfRangeError,
    RegressionEquation,
    Variable,
    available_states,
    evaluate,
    evaluate_weighted,
    jurisdiction,
    load_state,
    wave,
)
from flowfreq.regression.library import parse_state

CIT = Citation(publication="Synthetic test report", table="T1")


def _eq(**kw):
    base = dict(
        region_code="R1",
        aep=0.01,
        intercept=2.0,
        variables=(
            Variable("DRNAREA", "log10", "mi2", 1.0, 1000.0),
            Variable("PRECIP", "identity", "in"),
        ),
        coefficients=(0.7, 0.01),
        citation=CIT,
    )
    base.update(kw)
    return RegressionEquation(**base)


class TestSchema:
    def test_citation_required(self):
        with pytest.raises(ValueError):
            Citation(publication="", table="T1")

    @pytest.mark.parametrize(
        "kw",
        [
            {"aep": 1.0},
            {"coefficients": (0.7,)},
            {"covariance": ((1.0,),)},
            {"bias_correction": 0.0},
            {"variables": (Variable("A"), Variable("A")), "coefficients": (1.0, 1.0)},
        ],
    )
    def test_invalid(self, kw):
        with pytest.raises(ValueError):
            _eq(**kw)

    def test_unknown_transform(self):
        with pytest.raises(ValueError):
            Variable("A", transform="sqrt")

    def test_round_trip_from_dict(self):
        d = {
            "region_code": "R1",
            "aep": 0.01,
            "intercept": 2.0,
            "variables": [
                {"code": "DRNAREA", "transform": "log10", "minimum": 1.0, "maximum": 1000.0}
            ],
            "coefficients": [0.7],
            "citation": {"publication": "Synthetic", "table": "T1"},
            "covariance": [[0.01, 0.0], [0.0, 0.001]],
            "model_error_variance": 0.04,
            "n_sites": 50,
        }
        eq = RegressionEquation.from_dict(d)
        assert eq.covariance == ((0.01, 0.0), (0.0, 0.001))


class TestEvaluate:
    def test_value(self):
        est = evaluate(_eq(), {"DRNAREA": 100.0, "PRECIP": 30.0})
        assert est.log_flow == pytest.approx(2.0 + 0.7 * 2.0 + 0.01 * 30.0)
        assert est.flow_cfs == pytest.approx(10**est.log_flow)
        assert not est.extrapolated and est.interval is None

    def test_log10_plus1_transform(self):
        eq = _eq(variables=(Variable("STORAGE", "log10_plus1"),), coefficients=(-1.0,))
        assert evaluate(eq, {"STORAGE": 9.0}).log_flow == pytest.approx(1.0)

    def test_missing_characteristic(self):
        with pytest.raises(KeyError, match="PRECIP"):
            evaluate(_eq(), {"DRNAREA": 100.0})

    def test_out_of_range_raises_then_records(self):
        with pytest.raises(OutOfRangeError, match="DRNAREA"):
            evaluate(_eq(), {"DRNAREA": 5000.0, "PRECIP": 30.0})
        est = evaluate(_eq(), {"DRNAREA": 5000.0, "PRECIP": 30.0}, allow_extrapolation=True)
        assert est.extrapolated and "max" in est.out_of_range[0]

    def test_interval_from_sep(self):
        est = evaluate(_eq(sep_log=0.2), {"DRNAREA": 100.0, "PRECIP": 30.0})
        lo, hi = est.interval
        assert lo * hi == pytest.approx(est.flow_cfs**2)
        assert hi / est.flow_cfs == pytest.approx(10 ** (1.6448536 * 0.2), rel=1e-6)
        assert est.interval_method == "average standard error of prediction"

    def test_site_specific_interval_uses_t(self):
        cov = ((0.01, 0.0, 0.0), (0.0, 0.001, 0.0), (0.0, 0.0, 0.0))
        est = evaluate(
            _eq(covariance=cov, model_error_variance=0.04, n_sites=23),
            {"DRNAREA": 100.0, "PRECIP": 30.0},
        )
        sd = math.sqrt(0.04 + 0.01 + 0.001 * 4.0)
        from scipy import stats

        assert est.interval[1] / est.flow_cfs == pytest.approx(10 ** (stats.t.ppf(0.95, 20) * sd))
        assert est.interval_method.startswith("site-specific")

    def test_weighted(self):
        a = evaluate(_eq(), {"DRNAREA": 10.0, "PRECIP": 0.0})
        b = evaluate(_eq(), {"DRNAREA": 1000.0, "PRECIP": 0.0})
        assert evaluate_weighted([(0.5, a), (0.5, b)]) == pytest.approx(
            math.sqrt(a.flow_cfs * b.flow_cfs)
        )
        assert evaluate_weighted([(0.5, a), (0.5, b)], space="linear") == pytest.approx(
            (a.flow_cfs + b.flow_cfs) / 2
        )
        with pytest.raises(ValueError):
            evaluate_weighted([(0.4, a), (0.4, b)])


class TestJurisdictions:
    def test_56_with_one_wave_each(self):
        assert len(JURISDICTIONS) == 56
        assert len({j.code for j in JURISDICTIONS}) == 56
        assert {j.wave for j in JURISDICTIONS} == set(WAVE_NAMES)
        assert sum(j.is_territory for j in JURISDICTIONS) == 5

    def test_wave1_order(self):
        assert [j.code for j in wave("1")] == ["WA", "OR", "ID", "MT"]

    def test_unknown(self):
        with pytest.raises(KeyError):
            jurisdiction("ZZ")
        with pytest.raises(KeyError):
            wave("9")


class TestLibrary:
    def test_wave1_files_pending_and_empty(self):
        states = available_states()
        # Populated, each with its own tests/test_regression_<state>.py: WA, OR, ID, MT
        for code in ():
            assert states[code] == "pending"
            lib = load_state(code)
            assert lib.equations == []
            with pytest.raises(EquationsUnavailable, match="no equations transcribed"):
                lib.equation("anything", 0.01)

    def test_unstarted_state(self):
        with pytest.raises(EquationsUnavailable, match="wave has not started"):
            load_state("FL")

    def test_populated_library(self, tmp_path):
        eq = {
            "region_code": "R1",
            "aep": 0.01,
            "intercept": 2.0,
            "variables": [{"code": "DRNAREA"}],
            "coefficients": [0.7],
            "citation": {"publication": "Synthetic", "table": "T1"},
        }
        data = {
            "schema_version": 1,
            "state": "wa",
            "status": "partial",
            "equations": [eq, {**eq, "aep": 0.5}],
        }
        (tmp_path / "WA.json").write_text(json.dumps(data))
        lib = load_state("WA", data_dir=tmp_path)
        assert lib.regions == ["R1"] and lib.aeps("R1") == [0.5, 0.01]
        assert lib.equation("R1", 0.01).intercept == 2.0
        with pytest.raises(EquationsUnavailable):
            lib.equation("R1", 0.002)

    @pytest.mark.parametrize(
        "data, match",
        [
            ({"schema_version": 2, "state": "WA", "status": "pending"}, "schema_version"),
            ({"schema_version": 1, "state": "WA", "status": "done"}, "unknown status"),
            ({"schema_version": 1, "state": "WA", "status": "verified"}, "no equations"),
        ],
    )
    def test_parse_rejects(self, data, match):
        with pytest.raises(ValueError, match=match):
            parse_state(data)


class TestScaleOffset:
    """``Variable.scale`` and ``Variable.offset``: the value is ``T(scale*x + offset)``."""

    def _one(self, var, coef=1.0):
        return _eq(variables=(var,), coefficients=(coef,), intercept=0.0)

    def test_defaults(self):
        v = Variable("A")
        assert (v.scale, v.offset) == (1.0, 0.0)

    def test_scale_only(self):
        # (ELEV/1000)**0.5 at ELEV=4000 ft: 0.5 * log10(4)
        eq = self._one(Variable("ELEV", "log10", "ft", scale=0.001), 0.5)
        assert evaluate(eq, {"ELEV": 4000.0}).log_flow == pytest.approx(0.5 * math.log10(4.0))

    def test_offset_only(self):
        # (GUTTER + 0.1)**2
        eq = self._one(Variable("GUTTER", "log10", offset=0.1), 2.0)
        assert evaluate(eq, {"GUTTER": 0.9}).log_flow == pytest.approx(0.0)
        assert evaluate(eq, {"GUTTER": 9.9}).log_flow == pytest.approx(2.0)

    def test_scale_and_offset(self):
        # (FOREST/100 + 1)**-0.8 at FOREST=60 %: -0.8 * log10(1.6)
        eq = self._one(Variable("FOREST", "log10", "%", 0.0, 100.0, scale=0.01, offset=1.0), -0.8)
        est = evaluate(eq, {"FOREST": 60.0})
        assert est.log_flow == pytest.approx(-0.8 * math.log10(1.6), abs=1e-15)
        assert est.inputs == {"FOREST": 60.0}

    def test_identity_with_offset(self):
        eq = self._one(Variable("X", "identity", scale=2.0, offset=-20.0), 0.1)
        assert evaluate(eq, {"X": 15.0}).log_flow == pytest.approx(0.1 * (30.0 - 20.0))

    @pytest.mark.parametrize("x", [0.0, 0.5, 9.0, 123.456])
    @pytest.mark.parametrize("scale", [1.0, 0.01])
    def test_log10_plus1_is_log10_with_offset_1(self, x, scale):
        alias = Variable("S", "log10_plus1", scale=scale)
        explicit = Variable("S", "log10", scale=scale, offset=1.0)
        assert alias.apply(x) == explicit.apply(x)
        a = evaluate(self._one(alias, -1.3), {"S": x})
        b = evaluate(self._one(explicit, -1.3), {"S": x})
        assert a.log_flow == b.log_flow

    def test_log10_plus1_with_extra_offset_adds(self):
        assert Variable("S", "log10_plus1", offset=1.0).apply(8.0) == pytest.approx(1.0)

    @pytest.mark.parametrize(
        "var, x",
        [
            (Variable("A", "log10"), 0.0),
            (Variable("A", "log10"), -5.0),
            (Variable("ELEV", "log10", scale=0.001), 0.0),
            (Variable("X", "log10", offset=-20.0), 20.0),
            (Variable("X", "log10", offset=-20.0), 10.0),
            (Variable("F", "log10", scale=0.01, offset=1.0), -100.0),
            (Variable("S", "log10_plus1"), -1.0),
        ],
    )
    def test_nonpositive_argument_raises(self, var, x):
        with pytest.raises(ValueError, match=f"{var.code}.*needs"):
            var.apply(x)
        # Extrapolation cannot rescue an argument that has no logarithm.
        with pytest.raises(ValueError, match="R1 AEP 0.01"):
            evaluate(self._one(var), {var.code: x}, allow_extrapolation=True)

    def test_error_names_the_argument(self):
        with pytest.raises(ValueError, match=r"0\.01\*F \+ 1 > 0, got -0\.5"):
            Variable("F", scale=0.01, offset=1.0).apply(-150.0)

    @pytest.mark.parametrize(
        "kw", [{"scale": 0.0}, {"scale": -1.0}, {"scale": math.inf}, {"offset": math.nan}]
    )
    def test_invalid_scale_offset(self, kw):
        with pytest.raises(ValueError):
            Variable("A", **kw)

    def test_range_checked_on_raw_value(self):
        # Limits are in published units (percent), not in the scaled (fraction + 1) space.
        eq = self._one(Variable("FOREST", "log10", "%", 0.0, 100.0, scale=0.01, offset=1.0))
        assert not evaluate(eq, {"FOREST": 99.0}).extrapolated
        with pytest.raises(OutOfRangeError, match="FOREST=150.0 > max 100.0"):
            evaluate(eq, {"FOREST": 150.0})

    def test_covariance_in_transformed_basis(self):
        from scipy import stats

        eq = _eq(
            variables=(Variable("F", "log10", scale=0.01, offset=1.0),),
            coefficients=(1.0,),
            covariance=((0.01, 0.0), (0.0, 0.25)),
            model_error_variance=0.04,
            n_sites=30,
        )
        est = evaluate(eq, {"F": 60.0})
        t = math.log10(1.6)
        sd = math.sqrt(0.04 + 0.01 + 0.25 * t * t)
        assert est.interval[1] / est.flow_cfs == pytest.approx(10 ** (stats.t.ppf(0.95, 28) * sd))

    def test_json_round_trip(self):
        v = Variable("FOREST", "log10", "%", 0.0, 100.0, scale=0.01, offset=1.0)
        d = json.loads(json.dumps(v.to_dict()))
        assert d["scale"] == 0.01 and d["offset"] == 1.0
        assert Variable(**d) == v
        eq = RegressionEquation.from_dict(
            {
                "region_code": "R1",
                "aep": 0.01,
                "intercept": 0.0,
                "variables": [d],
                "coefficients": [-0.8],
                "citation": {"publication": "Synthetic", "table": "T1"},
            }
        )
        assert eq.variables == (v,)
        assert evaluate(eq, {"FOREST": 60.0}).log_flow == pytest.approx(-0.8 * math.log10(1.6))

    def test_to_dict_omits_defaults(self):
        d = {"code": "DRNAREA", "transform": "log10", "units": "mi2", "minimum": 1.0}
        assert Variable(**d).to_dict() == d
        assert "scale" not in Variable("A", offset=1.0).to_dict()
        assert "offset" not in Variable("A", scale=2.0).to_dict()

    def test_legacy_variable_dict_unchanged(self):
        # A variable written before scale/offset existed loads with the defaults and
        # evaluates exactly as the old ``TRANSFORMS[transform](x)`` did.
        from flowfreq.regression.equations import TRANSFORMS

        for t, x in [("log10", 250.0), ("log10_plus1", 3.7), ("identity", 42.0)]:
            assert Variable(**{"code": "A", "transform": t}).apply(x) == TRANSFORMS[t](x)

    def test_packaged_files_evaluate_as_before(self):
        # Every packaged equation that does not rescale evaluates exactly as the
        # pre-scale/offset evaluator did.
        from flowfreq.regression.equations import TRANSFORMS
        from flowfreq.regression.library import DATA_DIR

        for path in sorted(DATA_DIR.glob("*.json")):
            for eq in load_state(path.stem).equations:
                if any(v.scale != 1.0 or v.offset != 0.0 for v in eq.variables):
                    continue
                vals = {}
                for v in eq.variables:
                    lo = v.minimum if v.minimum is not None else 1.0
                    hi = v.maximum if v.maximum is not None else lo
                    vals[v.code] = (lo + hi) / 2 if lo <= 0 else math.sqrt(lo * hi)
                old = eq.intercept + sum(
                    b * TRANSFORMS[v.transform](vals[v.code])
                    for b, v in zip(eq.coefficients, eq.variables)
                )
                assert evaluate(eq, vals).log_flow == pytest.approx(old, rel=0, abs=1e-12)
