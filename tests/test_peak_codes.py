"""Tests for flowfreq.peak_codes, pinned to peakfq 8.1.0's siteQT rules."""

from __future__ import annotations

import math

import pytest

from flowfreq.peak_codes import (
    Q_MAX,
    Q_MIN,
    PeakTreatment,
    RegulationClass,
    classify_from_codes,
    parse_codes,
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
