"""Tests for PeakfqSA I/O converters."""

from __future__ import annotations

from pathlib import Path

import pytest

from hydrolib.peakfqsa.io_converters import DataFile, SpecificationFile
from tests.peakfqsa.fixtures.big_sandy import (
    BEGYEAR,
    ENDYEAR,
    HISTORICAL_PEAKS,
    REGIONAL_SKEW,
    REGIONAL_SKEW_SD,
    STATION_NAME,
    SYSTEMATIC_PEAKS,
    THRESHOLDS,
)


class TestSpecificationFile:
    """Tests for .psf file generation."""

    def test_from_analysis_params_creates_spec(self) -> None:
        """from_analysis_params produces a valid SpecificationFile."""
        spec = SpecificationFile.from_analysis_params(
            peaks=SYSTEMATIC_PEAKS,
            historical=HISTORICAL_PEAKS,
            thresholds=THRESHOLDS,
            begyear=BEGYEAR,
            endyear=ENDYEAR,
            regional_skew=REGIONAL_SKEW,
            regional_skew_sd=REGIONAL_SKEW_SD,
            station_name=STATION_NAME,
        )
        assert spec.station_name == STATION_NAME
        assert spec.begyear == BEGYEAR
        assert spec.endyear == ENDYEAR
        assert spec.regional_skew == REGIONAL_SKEW
        assert spec.skew_option == "Weighted"

    def test_to_string_contains_required_fields(self) -> None:
        """Generated PSF contains all required fields."""
        spec = SpecificationFile.from_analysis_params(
            peaks=SYSTEMATIC_PEAKS,
            historical=HISTORICAL_PEAKS,
            thresholds=THRESHOLDS,
            begyear=BEGYEAR,
            endyear=ENDYEAR,
            regional_skew=REGIONAL_SKEW,
            regional_skew_sd=REGIONAL_SKEW_SD,
            station_name=STATION_NAME,
        )
        content = spec.to_string()
        assert f"Station {STATION_NAME}" in content
        assert "PCPT_Thresh" in content
        assert "SkewOpt Weighted" in content
        assert f"GenSkew {REGIONAL_SKEW}" in content
        assert f"SkewSE {REGIONAL_SKEW_SD}" in content
        assert "LOType MGBT" in content
        assert "WeightOpt HWN" in content

    def test_to_string_has_thresholds(self) -> None:
        """Generated PSF has one PCPT_Thresh line per threshold."""
        spec = SpecificationFile.from_analysis_params(
            peaks=SYSTEMATIC_PEAKS,
            historical=HISTORICAL_PEAKS,
            thresholds=THRESHOLDS,
            begyear=BEGYEAR,
            endyear=ENDYEAR,
            regional_skew=REGIONAL_SKEW,
            regional_skew_sd=REGIONAL_SKEW_SD,
            station_name=STATION_NAME,
        )
        content = spec.to_string()
        thresh_lines = [l for l in content.splitlines() if "PCPT_Thresh" in l]
        assert len(thresh_lines) == len(THRESHOLDS)

    def test_infer_years_from_data(self) -> None:
        """Years are inferred from peaks when begyear/endyear are 0."""
        spec = SpecificationFile.from_analysis_params(
            peaks=SYSTEMATIC_PEAKS,
            historical=HISTORICAL_PEAKS,
            station_name="test",
        )
        assert spec.begyear == min(min(SYSTEMATIC_PEAKS), min(HISTORICAL_PEAKS))
        assert spec.endyear == max(max(SYSTEMATIC_PEAKS), max(HISTORICAL_PEAKS))

    def test_station_skew_option(self) -> None:
        """Station skew option when skew_sd <= 0."""
        spec = SpecificationFile.from_analysis_params(
            peaks=SYSTEMATIC_PEAKS,
            regional_skew_sd=0.0,
            station_name="test",
        )
        assert spec.skew_option == "Station"
        content = spec.to_string()
        assert "GenSkew" not in content
        assert "SkewSE" not in content

    def test_validate_success(self) -> None:
        """Valid spec passes validation."""
        spec = SpecificationFile.from_analysis_params(
            peaks=SYSTEMATIC_PEAKS,
            historical=HISTORICAL_PEAKS,
            thresholds=THRESHOLDS,
            begyear=BEGYEAR,
            endyear=ENDYEAR,
            regional_skew=REGIONAL_SKEW,
            regional_skew_sd=REGIONAL_SKEW_SD,
            station_name=STATION_NAME,
        )
        spec.validate()  # Should not raise

    def test_validate_missing_station(self) -> None:
        """Validation fails for missing station name."""
        spec = SpecificationFile(station_name="")
        with pytest.raises(ValueError, match="station_name"):
            spec.validate()

    def test_validate_bad_skew_option(self) -> None:
        """Validation fails for invalid skew option."""
        spec = SpecificationFile(station_name="test", skew_option="INVALID")
        with pytest.raises(ValueError, match="skew_option"):
            spec.validate()

    def test_validate_fixed_lo_requires_threshold(self) -> None:
        """Validation fails when FIXED lo_method has no threshold."""
        spec = SpecificationFile(station_name="test", lo_method="FIXED", lo_threshold=None)
        with pytest.raises(ValueError, match="lo_threshold"):
            spec.validate()

    def test_write_creates_file(self, tmp_path: Path) -> None:
        """write() creates a file on disk."""
        spec = SpecificationFile.from_analysis_params(
            peaks=SYSTEMATIC_PEAKS,
            station_name="test",
        )
        outfile = tmp_path / "test.psf"
        spec.write(outfile)
        assert outfile.exists()
        content = outfile.read_text()
        assert "Station test" in content


class TestDataFile:
    """Tests for data file generation."""

    def test_from_analysis_params(self) -> None:
        """from_analysis_params creates a DataFile."""
        df = DataFile.from_analysis_params(
            peaks=SYSTEMATIC_PEAKS,
            historical=HISTORICAL_PEAKS,
            station_name=STATION_NAME,
        )
        assert df.station_name == STATION_NAME
        assert len(df.peaks) == len(SYSTEMATIC_PEAKS)
        assert len(df.historical) == len(HISTORICAL_PEAKS)

    def test_to_string_has_header(self) -> None:
        """Generated file has RDB comment header."""
        df = DataFile.from_analysis_params(
            peaks=SYSTEMATIC_PEAKS,
            station_name=STATION_NAME,
        )
        content = df.to_string()
        assert content.startswith("#")
        assert STATION_NAME in content

    def test_to_string_has_all_peaks(self) -> None:
        """Every systematic peak appears in output."""
        df = DataFile.from_analysis_params(
            peaks=SYSTEMATIC_PEAKS,
            station_name="test",
        )
        content = df.to_string()
        for year, discharge in SYSTEMATIC_PEAKS.items():
            assert str(year) in content
            assert str(discharge) in content

    def test_to_string_has_historical_with_code7(self) -> None:
        """Historical peaks have qualification code 7."""
        df = DataFile.from_analysis_params(
            peaks=SYSTEMATIC_PEAKS,
            historical=HISTORICAL_PEAKS,
            station_name="test",
        )
        content = df.to_string()
        for year in HISTORICAL_PEAKS:
            # Find the line for this year and check it has code 7
            lines = [l for l in content.splitlines() if str(year) in l]
            assert any("7" in l for l in lines), f"Year {year} missing code 7"

    def test_to_string_column_headers(self) -> None:
        """Output has proper RDB column header."""
        df = DataFile.from_analysis_params(
            peaks={2000: 100.0},
            station_name="test",
        )
        content = df.to_string()
        assert "agency_cd" in content
        assert "site_no" in content
        assert "peak_dt" in content
        assert "peak_va" in content

    def test_validate_success(self) -> None:
        """Valid data passes validation."""
        df = DataFile.from_analysis_params(
            peaks=SYSTEMATIC_PEAKS,
            station_name="test",
        )
        df.validate()  # Should not raise

    def test_validate_missing_station(self) -> None:
        """Validation fails for missing station."""
        df = DataFile(station_name="", peaks={2000: 100.0})
        with pytest.raises(ValueError, match="station_name"):
            df.validate()

    def test_validate_no_data(self) -> None:
        """Validation fails when no data is provided."""
        df = DataFile(station_name="test")
        with pytest.raises(ValueError, match="observation"):
            df.validate()

    def test_validate_negative_discharge(self) -> None:
        """Validation fails for negative discharge."""
        df = DataFile(station_name="test", peaks={2000: -100.0})
        with pytest.raises(ValueError, match="Negative"):
            df.validate()

    def test_write_creates_file(self, tmp_path: Path) -> None:
        """write() creates a file on disk."""
        df = DataFile.from_analysis_params(
            peaks={2000: 100.0},
            station_name="test",
        )
        outfile = tmp_path / "test.dat"
        df.write(outfile)
        assert outfile.exists()
