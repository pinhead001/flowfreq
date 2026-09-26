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
        # ID is populated: tests/test_regression_idaho.py
        for code in ("WA", "OR", "MT"):
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
