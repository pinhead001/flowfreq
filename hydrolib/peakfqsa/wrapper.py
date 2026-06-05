"""
Subprocess wrapper for PeakfqSA Fortran executable.

Handles file I/O, process execution, and result parsing for running
Bulletin 17C flood frequency analysis via the PeakfqSA binary.

Reference: peakfqr R/fortranWrappers.R for call conventions.
"""

from __future__ import annotations

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from hydrolib.peakfqsa.config import PeakfqSAConfig, find_peakfqsa
from hydrolib.peakfqsa.io_converters import DataFile, SpecificationFile
from hydrolib.peakfqsa.parsers import PeakfqSAResult, parse_peakfqsa_output

logger = logging.getLogger(__name__)


class PeakfqSAExecutionError(RuntimeError):
    """Raised when PeakfqSA exits with a non-zero return code.

    Parameters
    ----------
    returncode : int
        Process exit code.
    stdout : str
        Standard output from the process.
    stderr : str
        Standard error from the process.
    """

    def __init__(self, returncode: int, stdout: str, stderr: str) -> None:
        msg = (
            f"PeakfqSA exited with code {returncode}.\n" f"stdout:\n{stdout}\n" f"stderr:\n{stderr}"
        )
        super().__init__(msg)
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


class PeakfqSATimeoutError(TimeoutError):
    """Raised when PeakfqSA exceeds the configured timeout."""

    pass


class PeakfqSAParseError(ValueError):
    """Raised when PeakfqSA output cannot be parsed."""

    pass


class PeakfqSAWrapper:
    """Subprocess wrapper for the PeakfqSA Fortran executable.

    Handles file I/O, execution, and result parsing for running
    Bulletin 17C flood frequency analysis.

    Parameters
    ----------
    config : PeakfqSAConfig
        Configuration for executable path, timeout, and temp file handling.
    """

    def __init__(self, config: PeakfqSAConfig) -> None:
        self.config = config
        self._executable: Path | None = config.executable_path

    def run(
        self,
        peaks: dict[int, float],
        historical: dict[int, float] | None = None,
        thresholds: list[dict[str, Any]] | None = None,
        begyear: int = 0,
        endyear: int = 0,
        regional_skew: float = -0.302,
        regional_skew_sd: float = 0.55,
        station_name: str = "",
        **kwargs: Any,
    ) -> PeakfqSAResult:
        """Run PeakfqSA analysis on the given peak flow data.

        Parameters
        ----------
        peaks : dict[int, float]
            Systematic annual peaks as {water_year: discharge_cfs}.
        historical : dict[int, float] or None
            Historical flood peaks as {water_year: discharge_cfs}.
        thresholds : list[dict] or None
            Perception thresholds, each with keys: start, end, lower, upper.
        begyear : int
            Beginning year of analysis period.
        endyear : int
            Ending year of analysis period.
        regional_skew : float
            Regional skew coefficient.
        regional_skew_sd : float
            Standard deviation of regional skew.
        station_name : str
            Station name for labeling output.
        **kwargs
            Additional keyword arguments passed to specification file.

        Returns
        -------
        PeakfqSAResult
            Parsed analysis results.

        Raises
        ------
        PeakfqSAExecutionError
            If PeakfqSA exits with non-zero code.
        PeakfqSATimeoutError
            If PeakfqSA exceeds timeout.
        PeakfqSAParseError
            If output cannot be parsed.
        """
        if self._executable is None:
            self._executable = find_peakfqsa()

        # Create temp directory (or use configured one)
        tmp_dir_ctx = None
        if self.config.temp_dir:
            tmp_dir = self.config.temp_dir
            tmp_dir.mkdir(parents=True, exist_ok=True)
        else:
            tmp_dir_ctx = tempfile.TemporaryDirectory(prefix="peakfqsa_")
            tmp_dir = Path(tmp_dir_ctx.name)

        try:
            # Write input files
            spec_path, data_path = self._write_input_files(
                tmp_dir=tmp_dir,
                peaks=peaks,
                historical=historical,
                thresholds=thresholds,
                begyear=begyear,
                endyear=endyear,
                regional_skew=regional_skew,
                regional_skew_sd=regional_skew_sd,
                station_name=station_name,
                **kwargs,
            )

            # Execute PeakfqSA
            output_text = self._execute(spec_path, tmp_dir)

            # Parse output
            result = self._parse_output_text(output_text)
            return result

        finally:
            if tmp_dir_ctx and not self.config.keep_temp_files:
                tmp_dir_ctx.cleanup()

    def _write_input_files(
        self,
        tmp_dir: Path,
        peaks: dict[int, float],
        historical: dict[int, float] | None,
        thresholds: list[dict[str, Any]] | None,
        begyear: int,
        endyear: int,
        regional_skew: float,
        regional_skew_sd: float,
        station_name: str,
        **kwargs: Any,
    ) -> tuple[Path, Path]:
        """Write .psf and .dat input files to the temp directory.

        Returns
        -------
        tuple[Path, Path]
            Paths to the specification file and data file.
        """
        data_filename = "peaks.dat"
        data_path = tmp_dir / data_filename

        # Build data file
        data_file = DataFile.from_analysis_params(
            peaks=peaks,
            historical=historical,
            station_name=station_name,
        )
        data_file.write(data_path)

        # Build specification file
        spec = SpecificationFile.from_analysis_params(
            peaks=peaks,
            historical=historical,
            thresholds=thresholds,
            begyear=begyear,
            endyear=endyear,
            regional_skew=regional_skew,
            regional_skew_sd=regional_skew_sd,
            station_name=station_name,
            data_filename=data_filename,
            **kwargs,
        )
        spec_path = tmp_dir / "analysis.psf"
        spec.write(spec_path)

        logger.info("Wrote input files to %s", tmp_dir)
        return spec_path, data_path

    def _execute(self, spec_file: Path, tmp_dir: Path) -> str:
        """Execute PeakfqSA and return raw output text.

        Parameters
        ----------
        spec_file : Path
            Path to the .psf specification file.
        tmp_dir : Path
            Working directory for execution.

        Returns
        -------
        str
            Raw content of the output file.

        Raises
        ------
        PeakfqSAExecutionError
            If the process exits with non-zero return code.
        PeakfqSATimeoutError
            If execution exceeds timeout.
        """
        cmd = [str(self._executable), str(spec_file)]
        logger.info("Running PeakfqSA: %s", " ".join(cmd))

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=str(tmp_dir),
                timeout=self.config.timeout_seconds,
            )
        except subprocess.TimeoutExpired as e:
            raise PeakfqSATimeoutError(
                f"PeakfqSA timed out after {self.config.timeout_seconds} seconds"
            ) from e

        if result.returncode != 0:
            raise PeakfqSAExecutionError(
                returncode=result.returncode,
                stdout=result.stdout,
                stderr=result.stderr,
            )

        # Look for output file (PeakfqSA typically writes .out with same stem)
        out_file = spec_file.with_suffix(".out")
        if not out_file.exists():
            # Try common alternatives
            for candidate in tmp_dir.glob("*.out"):
                out_file = candidate
                break
            else:
                # Fall back to stdout if no .out file found
                if result.stdout.strip():
                    logger.warning("No .out file found; using stdout")
                    return result.stdout
                raise PeakfqSAParseError(f"No output file found in {tmp_dir} and stdout is empty")

        output_text = out_file.read_text()
        logger.info("Read output file: %s (%d chars)", out_file, len(output_text))
        return output_text

    def _parse_output_text(self, text: str) -> PeakfqSAResult:
        """Parse PeakfqSA output text into a result object.

        Parameters
        ----------
        text : str
            Raw output text from PeakfqSA.

        Returns
        -------
        PeakfqSAResult
            Parsed results.

        Raises
        ------
        PeakfqSAParseError
            If the output cannot be parsed.
        """
        try:
            return parse_peakfqsa_output(text)
        except Exception as e:
            raise PeakfqSAParseError(f"Failed to parse PeakfqSA output: {e}") from e

    def is_available(self) -> bool:
        """Check whether PeakfqSA executable is available.

        Returns
        -------
        bool
            True if the executable can be found and validated.
        """
        try:
            if self._executable is None:
                self._executable = find_peakfqsa()
            return self._executable.is_file()
        except Exception:
            return False
