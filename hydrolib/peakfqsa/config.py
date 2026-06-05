"""
PeakfqSA configuration and executable detection.

Handles locating the PeakfqSA Fortran executable, validating it runs
correctly, and storing configuration for subprocess execution.
"""

from __future__ import annotations

import logging
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class PeakfqSANotFoundError(FileNotFoundError):
    """Raised when the PeakfqSA executable cannot be found.

    Parameters
    ----------
    searched_paths : list of Path
        Paths that were searched for the executable.
    """

    def __init__(self, searched_paths: list[Path] | None = None) -> None:
        paths_str = ""
        if searched_paths:
            paths_str = "\n  ".join(str(p) for p in searched_paths)
        msg = (
            "PeakfqSA executable not found.\n"
            "To fix this, do one of the following:\n"
            "  1. Set the PEAKFQSA_PATH environment variable to the executable path\n"
            "  2. Pass executable_path to PeakfqSAConfig\n"
            "  3. Place the PeakfqSA binary in one of these locations:\n"
            f"  {paths_str}"
        )
        super().__init__(msg)
        self.searched_paths = searched_paths or []


@dataclass
class PeakfqSAConfig:
    """Configuration for PeakfqSA subprocess execution.

    Parameters
    ----------
    executable_path : Path or None
        Path to the PeakfqSA executable. If None, auto-detection is attempted.
    timeout_seconds : int
        Maximum time in seconds to wait for PeakfqSA to complete.
    temp_dir : Path or None
        Directory for temporary input/output files. Uses system temp if None.
    keep_temp_files : bool
        If True, temporary files are not deleted after execution.
    """

    executable_path: Optional[Path] = None
    timeout_seconds: int = 60
    temp_dir: Optional[Path] = None
    keep_temp_files: bool = False

    def __post_init__(self) -> None:
        if self.executable_path is not None:
            self.executable_path = Path(self.executable_path)
        if self.temp_dir is not None:
            self.temp_dir = Path(self.temp_dir)


# Common locations to search for PeakfqSA
_SEARCH_PATHS: list[Path] = [
    Path.home() / "tools" / "peakfqsa",
    Path.home() / ".local" / "bin",
    Path("/usr/local/bin"),
    Path("/usr/bin"),
]

_EXECUTABLE_NAMES: list[str] = ["peakfqsa", "PeakfqSA", "peakfqsa.exe", "PeakfqSA.exe"]


def find_peakfqsa(user_path: Optional[Path] = None) -> Path:
    """Locate the PeakfqSA executable.

    Searches in order:
    1. ``user_path`` argument (if provided)
    2. ``PEAKFQSA_PATH`` environment variable
    3. System PATH via ``shutil.which``
    4. Common installation directories

    Parameters
    ----------
    user_path : Path or None
        Explicit path to check first.

    Returns
    -------
    Path
        Path to the PeakfqSA executable.

    Raises
    ------
    PeakfqSANotFoundError
        If the executable cannot be found in any searched location.
    """
    searched: list[Path] = []

    # 1. User-supplied path
    if user_path is not None:
        user_path = Path(user_path)
        searched.append(user_path)
        if user_path.is_file():
            logger.info("PeakfqSA found at user-supplied path: %s", user_path)
            return user_path

    # 2. Environment variable
    env_path = os.environ.get("PEAKFQSA_PATH")
    if env_path:
        p = Path(env_path)
        searched.append(p)
        if p.is_file():
            logger.info("PeakfqSA found via PEAKFQSA_PATH: %s", p)
            return p

    # 3. System PATH
    for name in _EXECUTABLE_NAMES:
        which_result = shutil.which(name)
        if which_result:
            p = Path(which_result)
            logger.info("PeakfqSA found on system PATH: %s", p)
            return p

    # 4. Common directories
    for directory in _SEARCH_PATHS:
        for name in _EXECUTABLE_NAMES:
            p = directory / name
            searched.append(p)
            if p.is_file():
                logger.info("PeakfqSA found at: %s", p)
                return p

    raise PeakfqSANotFoundError(searched)


def validate_peakfqsa(path: Path) -> str:
    """Validate that a PeakfqSA executable runs without crashing.

    Parameters
    ----------
    path : Path
        Path to the PeakfqSA executable.

    Returns
    -------
    str
        Version string extracted from PeakfqSA output, or "unknown" if
        version cannot be determined.

    Raises
    ------
    PeakfqSANotFoundError
        If the path does not exist.
    RuntimeError
        If the executable crashes on startup.
    """
    import subprocess

    path = Path(path)
    if not path.is_file():
        raise PeakfqSANotFoundError([path])

    try:
        result = subprocess.run(
            [str(path)],
            capture_output=True,
            text=True,
            timeout=10,
        )
        # PeakfqSA with no args typically prints usage/version and exits
        output = result.stdout + result.stderr
        logger.debug("PeakfqSA validation output: %s", output)

        # Try to extract version from output
        for line in output.splitlines():
            line_lower = line.lower()
            if "version" in line_lower or "peakfqsa" in line_lower:
                return line.strip()

        return "unknown"

    except subprocess.TimeoutExpired:
        raise RuntimeError(f"PeakfqSA at {path} timed out during validation")
    except OSError as e:
        raise RuntimeError(f"PeakfqSA at {path} failed to execute: {e}")
