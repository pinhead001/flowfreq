"""Tests for flowfreq.peak_codes, pinned to peakfq 8.1.0's siteQT rules."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from flowfreq.peak_codes import (
    Q_MAX,
    Q_MIN,
    PeakTreatment,
    RegulationClass,
    classify_from_codes,
    count_acted_on_codes,
    parse_codes,
    peak_frame_intervals,
    peak_interval,
)
from flowfreq.workflow import peak_code_kwargs, run_ffa


class TestParseCodes:
    @pytest.mark.parametrize(
        "raw, expected",
        [
            ("2,6", {"2", "6"}),
            ("26", {"2", "6"}),
            ("4, c", {"4", "C"}),
            ("", set()),
            (None, set()),
            (math.nan, set()),
            (["7", "A"], {"7", "A"}),
        ],
    )
    def test_forms(self, raw, expected):
        assert parse_codes(raw) == frozenset(expected)

    def test_unknown_code_kept_and_logged(self, caplog):
        assert parse_codes("Z") == frozenset({"Z"})
        assert "Unrecognised" in caplog.text


class TestPeakInterval:
    def test_systematic(self):
        iv = peak_interval(1000.0)
        assert (iv.lower, iv.upper, iv.treatment) == (1000.0, 1000.0, PeakTreatment.SYSTEMATIC)

    def test_code4_less_than(self):
        iv = peak_interval(12.0, "4")
        assert (iv.lower, iv.upper) == (Q_MIN, 12.0)

    def test_code8_greater_than(self):
        iv = peak_interval(500.0, "8")
        assert (iv.lower, iv.upper) == (500.0, Q_MAX)

    @pytest.mark.parametrize("code", ["3", "O"])
    def test_always_removed(self, code):
        iv = peak_interval(500.0, code, include_urban_regulated=True)
        assert (iv.lower, iv.upper, iv.treatment) == (Q_MIN, Q_MAX, PeakTreatment.REMOVED)

    @pytest.mark.parametrize("code", ["6", "C"])
    def test_urban_regulated_removed_by_default(self, code):
        assert peak_interval(500.0, code).treatment is PeakTreatment.REMOVED
        kept = peak_interval(500.0, code, include_urban_regulated=True)
        assert (kept.lower, kept.upper) == (500.0, 500.0)

    def test_removal_wins_over_censoring(self):
        # siteQT applies 4/8 first, then 3/O and 6/C, so removal wins.
        assert peak_interval(10.0, "4,C").treatment is PeakTreatment.REMOVED

    def test_historic_flag(self):
        assert peak_interval(9000.0, "7").is_historic

    def test_negative_raises(self):
        with pytest.raises(ValueError):
            peak_interval(-1.0)


class TestClassify:
    def test_regulated(self):
        cls, frac = classify_from_codes(["", "6", "6,2", ""])
        assert cls is RegulationClass.REGULATED and frac == 0.5

    def test_altered(self):
        assert classify_from_codes(["C", ""])[0] is RegulationClass.ALTERED

    def test_no_evidence_is_not_reference(self):
        cls, frac = classify_from_codes(["", "2"])
        assert cls is RegulationClass.NO_CODE_EVIDENCE and frac == 0.0
        assert "reference" not in {c.value for c in RegulationClass}

    def test_empty(self):
        assert classify_from_codes([])[0] is RegulationClass.NO_CODE_EVIDENCE


class TestParseCodesEncodings:
    """Encodings that used to split into codes nobody wrote."""

    @pytest.mark.parametrize(
        "raw, expected",
        [
            ("7.0", {"7"}),
            (7.0, {"7"}),
            ("2.0,7.0", {"2", "7"}),
            # Live legacy-NWIS values: Big Sandy 03606500 WY1897, Orestimba
            # 11274500 WY1947. These used to parse as {7, B, D} and {B, M}.
            ("7,Bd", {"7", "Bd"}),
            ("Bm", {"Bm"}),
            # WATSTORE run-together forms.
            ("7Bm", {"7", "Bm"}),
            ("4Bd", {"4", "Bd"}),
            ("4,BD", {"4", "Bd"}),
            (["7", "Bd"], {"7", "Bd"}),
        ],
    )
    def test_forms(self, raw, expected, caplog):
        assert parse_codes(raw) == frozenset(expected)
        assert "Unrecognised" not in caplog.text

    @pytest.mark.parametrize("date_code", ["Bd", "Bm"])
    def test_date_precision_codes_do_not_change_treatment(self, date_code):
        for base in ["", "4", "8", "6", "7"]:
            with_code = peak_interval(100.0, f"{base},{date_code}" if base else date_code)
            without = peak_interval(100.0, base)
            assert (with_code.lower, with_code.upper, with_code.treatment) == (
                without.lower,
                without.upper,
                without.treatment,
            )
            assert with_code.is_historic == without.is_historic

    def test_date_precision_codes_are_not_regulation_evidence(self):
        cls, _ = classify_from_codes(["Bd", "Bm", "7,Bd"])
        assert cls is RegulationClass.NO_CODE_EVIDENCE

    def test_float_code_drives_treatment(self):
        assert peak_interval(9000.0, "7.0").is_historic
        assert peak_interval(12.0, 4.0).treatment is PeakTreatment.LESS_THAN


class TestCountActedOnCodes:
    def test_counts(self):
        assert count_acted_on_codes(["", "2", "6,Bd", "7", "4", "6", None]) == {
            "4": 1,
            "6": 2,
            "7": 1,
        }

    def test_none(self):
        assert count_acted_on_codes(["", "2", "Bm", "5"]) == {}


class TestPeakFrameIntervals:
    def _frame(self):
        return pd.DataFrame(
            {
                "water_year": [2001, 2002, 2003, 2004, 2005],
                "peak_flow_cfs": [100.0, 12.0, 500.0, 900.0, float("nan")],
                "qualification_code": ["", "4", "6", "7", ""],
            }
        )

    def test_rows(self, caplog):
        out = peak_frame_intervals(self._frame())
        assert list(out["water_year"]) == [2001, 2002, 2003, 2004]  # NaN row dropped
        assert "no discharge" in caplog.text
        assert list(out["lower"]) == [100.0, Q_MIN, Q_MIN, 900.0]
        assert list(out["upper"]) == [100.0, 12.0, Q_MAX, 900.0]
        assert list(out["treatment"]) == [
            PeakTreatment.SYSTEMATIC,
            PeakTreatment.LESS_THAN,
            PeakTreatment.REMOVED,
            PeakTreatment.SYSTEMATIC,
        ]
        assert list(out["is_historic"]) == [False, False, False, True]

    def test_urban_regulated_kept(self):
        out = peak_frame_intervals(self._frame(), include_urban_regulated=True)
        assert out.loc[out["water_year"] == 2003, "lower"].iloc[0] == 500.0

    def test_no_code_column_is_systematic(self):
        out = peak_frame_intervals(self._frame().drop(columns="qualification_code"))
        assert (out["treatment"] == PeakTreatment.SYSTEMATIC).all()
        assert (out["lower"] == out["upper"]).all()

    def test_input_not_modified(self):
        frame = self._frame()
        peak_frame_intervals(frame)
        assert "lower" not in frame.columns

    def test_missing_column_raises(self):
        with pytest.raises(ValueError, match="missing columns"):
            peak_frame_intervals(self._frame().drop(columns="peak_flow_cfs"))

    def test_negative_raises(self):
        frame = self._frame()
        frame.loc[0, "peak_flow_cfs"] = -1.0
        with pytest.raises(ValueError, match="Negative"):
            peak_frame_intervals(frame)


# --------------------------------------------------------------------------- #
# Missing codes in every spelling, pd.NA above all
# --------------------------------------------------------------------------- #

#: Every way an uncoded peak reaches this module. ``pd.NA`` is what a
#: ``qualification_code`` column read with a nullable ``string`` dtype holds;
#: it used to raise ``TypeError`` (``iter(pd.NA)``).
_MISSING = [pd.NA, None, math.nan, np.nan, np.float64("nan"), ""]


def _nullable_codes(coded: dict, n: int) -> pd.Series:
    """A ``string``-dtype code column of length ``n``: ``pd.NA`` except at ``coded``."""
    out = pd.Series([pd.NA] * n, dtype="string")
    for index, code in coded.items():
        out[index] = code
    return out


class TestMissingCodes:
    @pytest.mark.parametrize("missing", _MISSING, ids=repr)
    def test_parse_codes_reads_each_as_no_code(self, missing):
        assert parse_codes(missing) == frozenset()

    def test_parse_codes_on_a_nullable_string_column(self):
        codes = _nullable_codes({0: "7", 2: "2,6"}, 4)
        assert [parse_codes(c) for c in codes] == [
            frozenset({"7"}),
            frozenset(),
            frozenset({"2", "6"}),
            frozenset(),
        ]

    def test_parse_codes_skips_missing_items_in_an_iterable(self):
        """A missing item is no code, not the characters of ``"<NA>"`` or ``"nan"``."""
        assert parse_codes(["7", pd.NA, None, math.nan, ""]) == frozenset({"7"})
        assert parse_codes(_nullable_codes({1: "4"}, 3)) == frozenset({"4"})

    def test_parse_codes_splits_comma_lists_inside_an_iterable(self):
        assert parse_codes(["2,6", "C"]) == frozenset({"2", "6", "C"})

    @pytest.mark.parametrize("missing", _MISSING, ids=repr)
    def test_peak_interval_is_systematic(self, missing):
        iv = peak_interval(500.0, missing)
        assert iv.treatment is PeakTreatment.SYSTEMATIC
        assert iv.lower == iv.upper == 500.0
        assert iv.codes == frozenset()

    def test_peak_interval_still_rejects_a_negative_value(self):
        with pytest.raises(ValueError, match="Negative"):
            peak_interval(-1.0, pd.NA)

    def test_classify_from_codes(self):
        assert classify_from_codes(_nullable_codes({}, 5)) == (
            RegulationClass.NO_CODE_EVIDENCE,
            0.0,
        )
        cls, frac = classify_from_codes(_nullable_codes({1: "6"}, 4))
        assert cls is RegulationClass.REGULATED
        assert frac == pytest.approx(0.25)

    def test_count_acted_on_codes(self):
        assert count_acted_on_codes(_nullable_codes({}, 5)) == {}
        assert count_acted_on_codes(_nullable_codes({0: "7", 3: "4,7"}, 5)) == {
            "4": 1,
            "7": 2,
        }

    def test_peak_frame_intervals(self):
        frame = pd.DataFrame(
            {
                "water_year": [2001, 2002, 2003, 2004],
                "peak_flow_cfs": [100.0, 200.0, 300.0, 400.0],
                "qualification_code": _nullable_codes({1: "4", 3: "7"}, 4),
            }
        )
        out = peak_frame_intervals(frame)
        assert list(out["treatment"]) == [
            PeakTreatment.SYSTEMATIC,
            PeakTreatment.LESS_THAN,
            PeakTreatment.SYSTEMATIC,
            PeakTreatment.SYSTEMATIC,
        ]
        assert list(out["is_historic"]) == [False, False, False, True]
        assert out["codes"].iloc[0] == frozenset()

    def test_peak_frame_intervals_error_case_unchanged(self):
        frame = pd.DataFrame({"water_year": [2001], "qualification_code": _nullable_codes({}, 1)})
        with pytest.raises(ValueError, match="missing columns"):
            peak_frame_intervals(frame)


_RUN_YEARS = np.arange(2001, 2016)
_RUN_FLOWS = np.array([float(100 + 37 * i % 211) for i in range(len(_RUN_YEARS))])


class TestWorkflowWithNullableCodes:
    """``workflow.peak_code_kwargs``, directly and through ``run_ffa(peak_codes=...)``."""

    def test_all_missing_matches_no_codes(self):
        plain = run_ffa(_RUN_FLOWS, _RUN_YEARS, station_skew_only=True)
        coded = run_ffa(
            _RUN_FLOWS,
            _RUN_YEARS,
            station_skew_only=True,
            peak_codes=_nullable_codes({}, len(_RUN_YEARS)),
        )
        assert coded["error"] is None
        assert coded["parameters"]["mean_log"] == plain["parameters"]["mean_log"]

    def test_coded_peaks_are_applied_alongside_missing_ones(self):
        coded = run_ffa(
            _RUN_FLOWS,
            _RUN_YEARS,
            station_skew_only=True,
            peak_codes=_nullable_codes({3: "6"}, len(_RUN_YEARS)),
        )
        assert coded["error"] is None
        assert coded["parameters"]["peak_codes_applied"] == {"6": 1}
        assert coded["b17c"].results.n_peaks == len(_RUN_YEARS) - 1

    def test_peak_code_kwargs_directly(self):
        n = len(_RUN_YEARS)
        assert peak_code_kwargs(_RUN_FLOWS, _RUN_YEARS, _nullable_codes({}, n)) is None
        kw = peak_code_kwargs(_RUN_FLOWS, _RUN_YEARS, _nullable_codes({0: "7"}, n))
        assert kw is not None
        assert kw["historical_peaks"] == [(2001, _RUN_FLOWS[0])]

    def test_misaligned_nullable_codes_still_raise(self):
        with pytest.raises(ValueError, match="aligned"):
            run_ffa(
                _RUN_FLOWS,
                _RUN_YEARS,
                station_skew_only=True,
                peak_codes=_nullable_codes({}, 3),
            )
