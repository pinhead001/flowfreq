"""Tests for PeakfqSA configuration and executable detection."""

from __future__ import annotations

from pathlib import Path

import pytest

from hydrolib.peakfqsa.config import (
    PeakfqSAConfig,
    PeakfqSANotFoundError,
    find_peakfqsa,
    validate_peakfqsa,
)


class TestPeakfqSAConfig:
    """Tests for PeakfqSAConfig dataclass."""

    def test_default_config(self) -> None:
        """Default config has reasonable values."""
        config = PeakfqSAConfig()
        assert config.executable_path is None
        assert config.timeout_seconds == 60
        assert config.temp_dir is None
        assert config.keep_temp_files is False

    def test_config_with_path(self, tmp_path: Path) -> None:
        """Config accepts and converts path."""
        exe = tmp_path / "peakfqsa"
        config = PeakfqSAConfig(executable_path=exe)
        assert config.executable_path == exe
        assert isinstance(config.executable_path, Path)

    def test_config_string_path_converted(self, tmp_path: Path) -> None:
        """String paths are converted to Path objects."""
        config = PeakfqSAConfig(executable_path=str(tmp_path / "peakfqsa"))
        assert isinstance(config.executable_path, Path)


class TestFindPeakfqsa:
    """Tests for executable auto-detection."""

    def test_user_path_found(self, tmp_path: Path) -> None:
        """User-supplied path is checked first."""
        exe = tmp_path / "peakfqsa"
        exe.touch()
        result = find_peakfqsa(user_path=exe)
        assert result == exe

    def test_user_path_not_found_raises(self, tmp_path: Path) -> None:
        """Missing user path raises PeakfqSANotFoundError."""
        with pytest.raises(PeakfqSANotFoundError):
            find_peakfqsa(user_path=tmp_path / "nonexistent")

    def test_env_var_override(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """PEAKFQSA_PATH environment variable is respected."""
        exe = tmp_path / "peakfqsa"
        exe.touch()
        monkeypatch.setenv("PEAKFQSA_PATH", str(exe))
        result = find_peakfqsa()
        assert result == exe

    def test_not_found_error_message_quality(self) -> None:
        """Error message includes fix instructions."""
        with pytest.raises(PeakfqSANotFoundError) as exc_info:
            find_peakfqsa(user_path=Path("/nonexistent/peakfqsa"))
        msg = str(exc_info.value)
        assert "PEAKFQSA_PATH" in msg
        assert "executable_path" in msg


class TestValidatePeakfqsa:
    """Tests for executable validation."""

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        """Validation of missing file raises PeakfqSANotFoundError."""
        with pytest.raises(PeakfqSANotFoundError):
            validate_peakfqsa(tmp_path / "nonexistent")
