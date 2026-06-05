"""
I/O converters for PeakfqSA input file formats.

Generates .psf (specification) and .dat (data) files matching the formats
consumed by the PeakfqSA Fortran executable. Format conventions derived
from peakfqr R/readInputs.R and R/fortranWrappers.R.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


class SpecificationFile:
    """PeakfqSA specification (.psf) file generator.

    Produces output matching the peakfqr .psf format. The format is
    line-oriented with keyword-value pairs::

        I RDB {data_filename}
        Station {station_name}
            PCPT_Thresh {start} {end} {lower} {upper} {comment}
            SkewOpt {Weighted|Station|Generalized}
            GenSkew {value}
            SkewSE {value}
            LOType {MGBT|FIXED|NONE}
            WeightOpt {HWN|INV|ERL}

    Parameters
    ----------
    data_filename : str
        Path to the associated data file.
    station_name : str
        Station identifier/name.
    begyear : int
        Beginning year of analysis.
    endyear : int
        Ending year of analysis.
    regional_skew : float
        Regional skew coefficient.
    skew_sd : float
        Standard deviation of regional skew.
    thresholds : list[dict]
        Perception thresholds with keys: start, end, lower, upper.
    skew_option : str
        Skew option: 'Weighted', 'Station', or 'Generalized'.
    lo_method : str
        Low outlier method: 'MGBT', 'FIXED', or 'NONE'.
    lo_threshold : float or None
        Fixed low outlier threshold (required if lo_method is 'FIXED').
    weight_opt : str
        Skew weighting algorithm: 'HWN', 'INV', or 'ERL'.
    intervals : list[dict] or None
        Flow intervals with keys: year, lower, upper, comment.
    """

    def __init__(
        self,
        data_filename: str = "",
        station_name: str = "",
        begyear: int = 0,
        endyear: int = 0,
        regional_skew: float = -0.302,
        skew_sd: float = 0.55,
        thresholds: list[dict[str, Any]] | None = None,
        skew_option: str = "Weighted",
        lo_method: str = "MGBT",
        lo_threshold: float | None = None,
        weight_opt: str = "HWN",
        intervals: list[dict[str, Any]] | None = None,
    ) -> None:
        self.data_filename = data_filename
        self.station_name = station_name
        self.begyear = begyear
        self.endyear = endyear
        self.regional_skew = regional_skew
        self.skew_sd = skew_sd
        self.thresholds = thresholds or []
        self.skew_option = skew_option
        self.lo_method = lo_method
        self.lo_threshold = lo_threshold
        self.weight_opt = weight_opt
        self.intervals = intervals or []

    @classmethod
    def from_analysis_params(
        cls,
        peaks: dict[int, float],
        historical: dict[int, float] | None = None,
        thresholds: list[dict[str, Any]] | None = None,
        begyear: int = 0,
        endyear: int = 0,
        regional_skew: float = -0.302,
        regional_skew_sd: float = 0.55,
        station_name: str = "",
        data_filename: str = "peaks.dat",
        lo_method: str = "MGBT",
        weight_opt: str = "HWN",
        **kwargs: Any,
    ) -> SpecificationFile:
        """Create a SpecificationFile from analysis parameters.

        Parameters
        ----------
        peaks : dict[int, float]
            Systematic peaks {year: discharge}.
        historical : dict[int, float] or None
            Historical peaks {year: discharge}.
        thresholds : list[dict] or None
            Perception thresholds.
        begyear : int
            Start year. Inferred from data if 0.
        endyear : int
            End year. Inferred from data if 0.
        regional_skew : float
            Regional skew coefficient.
        regional_skew_sd : float
            Standard deviation of regional skew.
        station_name : str
            Station name.
        data_filename : str
            Filename for the associated data file.
        lo_method : str
            Low outlier method.
        weight_opt : str
            Skew weighting option.
        **kwargs
            Additional parameters.

        Returns
        -------
        SpecificationFile
            Configured specification file.
        """
        historical = historical or {}
        all_years = list(peaks.keys()) + list(historical.keys())

        if begyear == 0:
            begyear = min(all_years)
        if endyear == 0:
            endyear = max(all_years)

        # Default thresholds if none provided: systematic record covers full period
        if not thresholds:
            thresholds = [
                {"start": begyear, "end": endyear, "lower": 0.0, "upper": 1e20},
            ]

        # Determine skew option based on whether regional skew is meaningful
        skew_option = "Weighted" if regional_skew_sd > 0 else "Station"

        return cls(
            data_filename=data_filename,
            station_name=station_name,
            begyear=begyear,
            endyear=endyear,
            regional_skew=regional_skew,
            skew_sd=regional_skew_sd,
            thresholds=thresholds,
            skew_option=skew_option,
            lo_method=lo_method,
            weight_opt=weight_opt,
        )

    def to_string(self) -> str:
        """Generate the .psf file content as a string.

        Returns
        -------
        str
            Complete .psf file content.
        """
        lines: list[str] = []
        lines.append("' Generated by HydroLib")

        # Input file reference
        if self.data_filename:
            lines.append(f"I RDB {self.data_filename}")

        # Station block
        lines.append(f"Station {self.station_name}")

        # Perception thresholds
        for t in self.thresholds:
            start = int(t["start"])
            end = int(t["end"])
            lower = t["lower"]
            upper = t["upper"]
            comment = t.get("comment", "")
            # Format lower/upper: use scientific notation for large values
            lower_str = _format_threshold(lower)
            upper_str = _format_threshold(upper)
            line = f"    PCPT_Thresh {start} {end} {lower_str} {upper_str}"
            if comment:
                line += f" {comment}"
            lines.append(line)

        # Flow intervals
        for iv in self.intervals:
            year = int(iv["year"])
            lower = iv["lower"]
            upper = iv["upper"]
            comment = iv.get("comment", "")
            line = f"    Interval {year} {lower} {upper}"
            if comment:
                line += f" {comment}"
            lines.append(line)

        # Skew options
        lines.append(f"    SkewOpt {self.skew_option}")
        if self.skew_option in ("Weighted", "Generalized"):
            lines.append(f"    GenSkew {self.regional_skew}")
            lines.append(f"    SkewSE {self.skew_sd}")

        # Low outlier options
        lines.append(f"    LOType {self.lo_method}")
        if self.lo_method == "FIXED" and self.lo_threshold is not None:
            lines.append(f"    LoThresh {self.lo_threshold}")

        # Weighting option
        lines.append(f"    WeightOpt {self.weight_opt}")

        return "\n".join(lines) + "\n"

    def write(self, path: Path) -> None:
        """Write the specification file to disk.

        Parameters
        ----------
        path : Path
            Output file path.
        """
        path = Path(path)
        path.write_text(self.to_string())
        logger.info("Wrote specification file: %s", path)

    def validate(self) -> None:
        """Check internal consistency of the specification.

        Raises
        ------
        ValueError
            If the specification is invalid.
        """
        if not self.station_name:
            raise ValueError("station_name is required")

        if self.begyear > self.endyear:
            raise ValueError(f"begyear ({self.begyear}) must be <= endyear ({self.endyear})")

        valid_skew_opts = ("Station", "Weighted", "Generalized")
        if self.skew_option not in valid_skew_opts:
            raise ValueError(
                f"skew_option must be one of {valid_skew_opts}, got '{self.skew_option}'"
            )

        if self.skew_option in ("Weighted", "Generalized"):
            if self.skew_sd <= 0:
                raise ValueError("skew_sd must be positive for Weighted/Generalized skew")

        valid_lo = ("MGBT", "FIXED", "NONE")
        if self.lo_method not in valid_lo:
            raise ValueError(f"lo_method must be one of {valid_lo}, got '{self.lo_method}'")

        if self.lo_method == "FIXED" and self.lo_threshold is None:
            raise ValueError("lo_threshold is required when lo_method is 'FIXED'")

        valid_weight = ("HWN", "INV", "ERL")
        if self.weight_opt not in valid_weight:
            raise ValueError(f"weight_opt must be one of {valid_weight}, got '{self.weight_opt}'")

        for i, t in enumerate(self.thresholds):
            for key in ("start", "end", "lower", "upper"):
                if key not in t:
                    raise ValueError(f"Threshold {i} missing required key '{key}'")
            if t["start"] > t["end"]:
                raise ValueError(f"Threshold {i}: start ({t['start']}) > end ({t['end']})")


class DataFile:
    """PeakfqSA data (.dat) file generator.

    Produces peak flow data in a tab-delimited RDB-like format compatible
    with peakfqr's ``readRDB()`` function. Columns: agency_cd, site_no,
    peak_dt, peak_va, peak_cd.

    Parameters
    ----------
    station_name : str
        Station identifier for the site_no column.
    peaks : dict[int, float]
        Point observations {water_year: discharge_cfs}.
    historical : dict[int, float] or None
        Historical peaks {water_year: discharge_cfs}.
    intervals : dict[int, tuple[float, float]] or None
        Interval observations {water_year: (lower, upper)}.
    """

    def __init__(
        self,
        station_name: str = "",
        peaks: dict[int, float] | None = None,
        historical: dict[int, float] | None = None,
        intervals: dict[int, tuple[float, float]] | None = None,
    ) -> None:
        self.station_name = station_name
        self.peaks = peaks or {}
        self.historical = historical or {}
        self.intervals = intervals or {}

    @classmethod
    def from_analysis_params(
        cls,
        peaks: dict[int, float],
        historical: dict[int, float] | None = None,
        station_name: str = "",
        **kwargs: Any,
    ) -> DataFile:
        """Create a DataFile from analysis parameters.

        Parameters
        ----------
        peaks : dict[int, float]
            Systematic peaks {year: discharge}.
        historical : dict[int, float] or None
            Historical peaks {year: discharge}.
        station_name : str
            Station name/ID.
        **kwargs
            Additional parameters.

        Returns
        -------
        DataFile
            Configured data file.
        """
        return cls(
            station_name=station_name,
            peaks=peaks,
            historical=historical or {},
        )

    def to_string(self) -> str:
        """Generate the data file content as a tab-delimited RDB string.

        The format matches NWIS RDB peak flow files that peakfqr's
        ``readRDB()`` can parse: comment header, column names, format line,
        then tab-separated data rows.

        Returns
        -------
        str
            Complete data file content.
        """
        lines: list[str] = []

        # RDB comment header
        lines.append(f"# Generated by HydroLib")
        lines.append(f"# Station: {self.station_name}")
        lines.append("#")

        # Column header
        cols = ["agency_cd", "site_no", "peak_dt", "peak_va", "peak_cd"]
        lines.append("\t".join(cols))
        # Format line (RDB convention: field widths)
        lines.append("\t".join(["5s", "15s", "10d", "8s", "27s"]))

        # Systematic peaks
        for year, discharge in sorted(self.peaks.items()):
            date_str = f"{year}-01-01"  # Use Jan 1 as placeholder date
            lines.append(f"USGS\t{self.station_name}\t{date_str}\t{discharge}\t")

        # Historical peaks (code 7 = historic peak)
        for year, discharge in sorted(self.historical.items()):
            date_str = f"{year}-01-01"
            lines.append(f"USGS\t{self.station_name}\t{date_str}\t{discharge}\t7")

        # Interval observations (code 4 = less than, code 8 = greater than)
        for year, (lower, upper) in sorted(self.intervals.items()):
            date_str = f"{year}-01-01"
            if lower == 0:
                # Less-than: discharge < upper
                lines.append(f"USGS\t{self.station_name}\t{date_str}\t{upper}\t4")
            else:
                # Greater-than: discharge > lower
                lines.append(f"USGS\t{self.station_name}\t{date_str}\t{lower}\t8")

        return "\n".join(lines) + "\n"

    def write(self, path: Path) -> None:
        """Write the data file to disk.

        Parameters
        ----------
        path : Path
            Output file path.
        """
        path = Path(path)
        path.write_text(self.to_string())
        logger.info("Wrote data file: %s", path)

    def validate(self) -> None:
        """Check internal consistency of the data file.

        Raises
        ------
        ValueError
            If the data is invalid.
        """
        if not self.station_name:
            raise ValueError("station_name is required")

        all_peaks = {**self.peaks, **self.historical}
        if not all_peaks and not self.intervals:
            raise ValueError("At least one peak or interval observation is required")

        for year, q in all_peaks.items():
            if q < 0:
                raise ValueError(f"Negative discharge {q} in year {year}")

        for year, (lower, upper) in self.intervals.items():
            if lower > upper:
                raise ValueError(f"Interval in year {year}: lower ({lower}) > upper ({upper})")


def _format_threshold(value: float) -> str:
    """Format a threshold value for PSF output.

    Parameters
    ----------
    value : float
        Threshold value.

    Returns
    -------
    str
        Formatted string: '0' for zero, scientific notation for large
        values, plain number otherwise.
    """
    if value == 0:
        return "0"
    if abs(value) >= 1e10:
        return f"{value:.0E}"
    return str(value)
