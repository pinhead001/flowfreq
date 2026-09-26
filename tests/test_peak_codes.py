"""Tests for flowfreq.peak_codes, pinned to peakfq 8.1.0's siteQT rules."""

from __future__ import annotations

import math

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
