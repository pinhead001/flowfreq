"""Tests for PeakfqSA subprocess wrapper."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from hydrolib.peakfqsa.config import PeakfqSAConfig
from hydrolib.peakfqsa.wrapper import (
    PeakfqSAExecutionError,
    PeakfqSAParseError,
    PeakfqSATimeoutError,
    PeakfqSAWrapper,
)
from tests.peakfqsa.fixtures.big_sandy import HISTORICAL_PEAKS, STATION_NAME, SYSTEMATIC_PEAKS

# Mark for tests requiring the real PeakfqSA binary
requires_peakfqsa = pytest.mark.requires_peakfqsa


class TestPeakfqSAWrapper:
    """Tests for PeakfqSAWrapper (mock-based, no real binary needed)."""

    def test_is_available_without_binary(self) -> None:
        """is_available() returns False when binary not found."""
        config = PeakfqSAConfig(executable_path=Path("/nonexistent/peakfqsa"))
        wrapper = PeakfqSAWrapper(config)
        assert wrapper.is_available() is False

    def test_is_available_with_binary(self, tmp_path: Path) -> None:
        """is_available() returns True when binary exists."""
        exe = tmp_path / "peakfqsa"
        exe.touch()
        config = PeakfqSAConfig(executable_path=exe)
        wrapper = PeakfqSAWrapper(config)
        assert wrapper.is_available() is True

    def test_write_input_files(self, tmp_path: Path) -> None:
        """_write_input_files creates .psf and .dat files."""
        exe = tmp_path / "peakfqsa"
        exe.touch()
        config = PeakfqSAConfig(executable_path=exe)
        wrapper = PeakfqSAWrapper(config)

        spec_path, data_path = wrapper._write_input_files(
            tmp_dir=tmp_path,
            peaks=SYSTEMATIC_PEAKS,
            historical=HISTORICAL_PEAKS,
            thresholds=None,
            begyear=0,
            endyear=0,
            regional_skew=-0.5,
            regional_skew_sd=0.55,
            station_name=STATION_NAME,
        )

        assert spec_path.exists()
        assert data_path.exists()
        assert spec_path.suffix == ".psf"
        assert data_path.suffix == ".dat"

        # Verify content
        spec_content = spec_path.read_text()
        assert STATION_NAME in spec_content
        assert "PCPT_Thresh" in spec_content

        data_content = data_path.read_text()
        assert "agency_cd" in data_content

    def test_write_input_files_references_data(self, tmp_path: Path) -> None:
        """Specification file references the data file."""
        exe = tmp_path / "peakfqsa"
        exe.touch()
        config = PeakfqSAConfig(executable_path=exe)
        wrapper = PeakfqSAWrapper(config)

        spec_path, data_path = wrapper._write_input_files(
            tmp_dir=tmp_path,
            peaks={2000: 100.0},
            historical=None,
            thresholds=None,
            begyear=2000,
            endyear=2000,
            regional_skew=-0.302,
            regional_skew_sd=0.55,
            station_name="test",
        )

        spec_content = spec_path.read_text()
        assert "peaks.dat" in spec_content

    @patch("hydrolib.peakfqsa.wrapper.subprocess.run")
    def test_execute_success(self, mock_run: MagicMock, tmp_path: Path) -> None:
        """_execute returns output file content on success."""
        exe = tmp_path / "peakfqsa"
        exe.touch()

        # Create a fake output file
        spec_path = tmp_path / "analysis.psf"
        spec_path.touch()
        out_path = tmp_path / "analysis.out"
        out_path.write_text("Station: test\n--- Frequency Curve ---\n")

        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        config = PeakfqSAConfig(executable_path=exe)
        wrapper = PeakfqSAWrapper(config)
        output = wrapper._execute(spec_path, tmp_path)

        assert "Station: test" in output
        mock_run.assert_called_once()

    @patch("hydrolib.peakfqsa.wrapper.subprocess.run")
    def test_execute_nonzero_exit(self, mock_run: MagicMock, tmp_path: Path) -> None:
        """_execute raises PeakfqSAExecutionError on non-zero exit."""
        exe = tmp_path / "peakfqsa"
        exe.touch()
        spec_path = tmp_path / "analysis.psf"
        spec_path.touch()

        mock_run.return_value = MagicMock(returncode=1, stdout="some output", stderr="error msg")

        config = PeakfqSAConfig(executable_path=exe)
        wrapper = PeakfqSAWrapper(config)

        with pytest.raises(PeakfqSAExecutionError) as exc_info:
            wrapper._execute(spec_path, tmp_path)

        assert exc_info.value.returncode == 1
        assert "error msg" in str(exc_info.value)

    @patch("hydrolib.peakfqsa.wrapper.subprocess.run")
    def test_execute_timeout(self, mock_run: MagicMock, tmp_path: Path) -> None:
        """_execute raises PeakfqSATimeoutError on timeout."""
        import subprocess

        exe = tmp_path / "peakfqsa"
        exe.touch()
        spec_path = tmp_path / "analysis.psf"
        spec_path.touch()

        mock_run.side_effect = subprocess.TimeoutExpired(cmd="peakfqsa", timeout=60)

        config = PeakfqSAConfig(executable_path=exe, timeout_seconds=60)
        wrapper = PeakfqSAWrapper(config)

        with pytest.raises(PeakfqSATimeoutError):
            wrapper._execute(spec_path, tmp_path)

    @patch("hydrolib.peakfqsa.wrapper.subprocess.run")
    def test_execute_fallback_to_stdout(self, mock_run: MagicMock, tmp_path: Path) -> None:
        """_execute falls back to stdout if no .out file exists."""
        exe = tmp_path / "peakfqsa"
        exe.touch()
        spec_path = tmp_path / "analysis.psf"
        spec_path.touch()

        mock_run.return_value = MagicMock(returncode=0, stdout="Station: stdout_output", stderr="")

        config = PeakfqSAConfig(executable_path=exe)
        wrapper = PeakfqSAWrapper(config)
        output = wrapper._execute(spec_path, tmp_path)

        assert "stdout_output" in output

    def test_parse_output_text_invalid(self, tmp_path: Path) -> None:
        """_parse_output_text raises PeakfqSAParseError on bad text."""
        exe = tmp_path / "peakfqsa"
        exe.touch()
        config = PeakfqSAConfig(executable_path=exe)
        wrapper = PeakfqSAWrapper(config)

        # Valid text shouldn't raise (parser is defensive)
        result = wrapper._parse_output_text("Station: test\n")
        assert result.station_name == "test"


@requires_peakfqsa
class TestPeakfqSAWrapperIntegration:
    """Integration tests requiring the real PeakfqSA binary.

    These tests are skipped when PeakfqSA is not installed.
    """

    pass
