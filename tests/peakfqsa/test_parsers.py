"""Tests for PeakfqSA output parser."""

from __future__ import annotations

import pytest

from hydrolib.peakfqsa.parsers import PeakfqSAResult, parse_peakfqsa_output
from tests.peakfqsa.fixtures.big_sandy import (
    EXPECTED_CONFIDENCE_INTERVALS,
    EXPECTED_PARAMETERS,
    EXPECTED_QUANTILES,
    STATION_NAME,
)
from tests.peakfqsa.fixtures.sample_output import (
    BIG_SANDY_OUTPUT,
    MINIMAL_OUTPUT,
    STATION_SKEW_OUTPUT,
)


class TestPeakfqSAResult:
    """Tests for the PeakfqSAResult dataclass."""

    def test_default_construction(self) -> None:
        """Default result has empty/zero fields."""
        result = PeakfqSAResult()
        assert result.station_name == ""
        assert result.n_peaks == 0
        assert result.quantiles == {}
        assert result.raw_output == ""

    def test_raw_output_stored(self) -> None:
        """raw_output field preserves full text."""
        text = "some output text"
        result = PeakfqSAResult(raw_output=text)
        assert result.raw_output == text


class TestParseOutput:
    """Tests for parse_peakfqsa_output function."""

    def test_parse_station_name(self) -> None:
        """Parser extracts station name."""
        result = parse_peakfqsa_output(BIG_SANDY_OUTPUT)
        assert STATION_NAME in result.station_name

    def test_parse_analysis_period(self) -> None:
        """Parser extracts begin/end years."""
        result = parse_peakfqsa_output(BIG_SANDY_OUTPUT)
        assert result.begyear == 1890
        assert result.endyear == 1973

    def test_parse_peak_counts(self) -> None:
        """Parser extracts peak counts."""
        result = parse_peakfqsa_output(BIG_SANDY_OUTPUT)
        assert result.n_peaks == 47
        assert result.n_systematic == 44
        assert result.n_historical == 3

    def test_parse_parameters(self) -> None:
        """Parser extracts LP3 parameters."""
        result = parse_peakfqsa_output(BIG_SANDY_OUTPUT)
        assert "mean_log" in result.parameters
        assert "std_log" in result.parameters
        assert "skew_weighted" in result.parameters

        # Check values match expected
        assert abs(result.parameters["mean_log"] - EXPECTED_PARAMETERS["mean_log"]) < 0.001
        assert abs(result.parameters["std_log"] - EXPECTED_PARAMETERS["std_log"]) < 0.001
        assert (
            abs(result.parameters["skew_weighted"] - EXPECTED_PARAMETERS["skew_weighted"]) < 0.001
        )

    def test_parse_quantiles(self) -> None:
        """Parser extracts all quantiles from frequency curve."""
        result = parse_peakfqsa_output(BIG_SANDY_OUTPUT)
        assert len(result.quantiles) == 14

        # Spot-check some key quantiles
        for aep, expected_q in EXPECTED_QUANTILES.items():
            if aep in result.quantiles:
                pct_diff = abs(result.quantiles[aep] - expected_q) / expected_q * 100
                assert pct_diff < 1.0, f"AEP {aep}: {result.quantiles[aep]} vs {expected_q}"

    def test_parse_confidence_intervals(self) -> None:
        """Parser extracts confidence intervals."""
        result = parse_peakfqsa_output(BIG_SANDY_OUTPUT)
        assert len(result.confidence_intervals) == 14

        # Check specific CIs from expected
        for aep, (exp_lo, exp_hi) in EXPECTED_CONFIDENCE_INTERVALS.items():
            if aep in result.confidence_intervals:
                lo, hi = result.confidence_intervals[aep]
                assert abs(lo - exp_lo) / exp_lo * 100 < 1.0
                assert abs(hi - exp_hi) / exp_hi * 100 < 1.0

    def test_parse_low_outlier(self) -> None:
        """Parser extracts low outlier information."""
        result = parse_peakfqsa_output(BIG_SANDY_OUTPUT)
        assert result.low_outlier_count == 0
        assert result.low_outlier_threshold == 1234.56

    def test_parse_raw_output_stored(self) -> None:
        """Parser stores raw output text."""
        result = parse_peakfqsa_output(BIG_SANDY_OUTPUT)
        assert result.raw_output == BIG_SANDY_OUTPUT

    def test_parse_minimal_output(self) -> None:
        """Parser handles minimal output."""
        result = parse_peakfqsa_output(MINIMAL_OUTPUT)
        assert result.station_name == "TEST_STATION"
        assert result.n_peaks == 21
        assert len(result.quantiles) == 2
        assert 0.01 in result.quantiles
        assert result.quantiles[0.01] == 5000.0

    def test_parse_station_skew(self) -> None:
        """Parser handles station skew output."""
        result = parse_peakfqsa_output(STATION_SKEW_OUTPUT)
        assert result.station_name == "STATION_SKEW_TEST"
        assert result.begyear == 1950
        assert result.endyear == 2000
        assert result.n_peaks == 51

    def test_parse_empty_output(self) -> None:
        """Parser handles empty/minimal text gracefully."""
        result = parse_peakfqsa_output("")
        assert result.station_name == ""
        assert result.quantiles == {}
        assert result.parameters == {}
