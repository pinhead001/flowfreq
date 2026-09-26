"""Reader for PeakFQ ``.psf`` specification files.

Follows ``vendor/peakfqr/R/readInputs.R`` (the PSF reader documented there as
``readPSF``/``siteQT`` inputs): a file is a header of ``I``/``O`` option
lines, then one block per ``Station``. Within a block, ``PCPT_Thresh``,
``Interval`` and the deprecated ``Peak`` may repeat; every other keyword is a
single value. Lines whose first token starts with ``'`` are comments.

Where the R reader warns and continues (unknown or deprecated keywords), so
does this one, through :mod:`logging`. It never guesses a value for a field
that is missing.

Roadmap: ``docs/MASTER_ROADMAP.md`` §1.1, issue #31. Converting a parsed
station into ``Bulletin17C`` arguments is the open half of that issue.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Union

logger = logging.getLogger(__name__)

#: Single-value station fields read by peakfq 8.1.0 (``PSFcolumns``).
STATION_FIELDS = ("SkewOpt", "GenSkew", "SkewSE", "LOType", "LoThresh", "Urb/Reg", "WeightOpt")
#: Fields peakfq 8.1.0 recognises as deprecated and ignores (``PSFcolumns_old``).
DEPRECATED_FIELDS = ("Analyze", "BegYear", "EndYear", "HistPeriod")


@dataclass(frozen=True)
class PerceptionThreshold:
    """``PCPT_Thresh StartYear EndYear Lower Upper [comment]``."""

    start: int
    end: int
    lower: float
    upper: float
    comment: str = ""


@dataclass(frozen=True)
class SpecInterval:
    """``Interval Year Lower Upper [comment]``."""

    year: int
    lower: float
    upper: float
    comment: str = ""


@dataclass(frozen=True)
class SpecPeak:
    """Deprecated ``Peak Year Value [comment]``."""

    year: int
    value: float
    comment: str = ""


@dataclass
class StationSpec:
    """Specifications for one station."""

    station_id: str
    thresholds: List[PerceptionThreshold] = field(default_factory=list)
    intervals: List[SpecInterval] = field(default_factory=list)
    peaks: List[SpecPeak] = field(default_factory=list)
    fields: Dict[str, str] = field(default_factory=dict)

    @property
    def site_no(self) -> str:
        """The station ID without PeakFQ's ``.NN`` analysis suffix."""
        return self.station_id.split(".")[0]

    @property
    def skew_option(self) -> Optional[str]:
        """``Station``, ``Weighted`` or ``Generalized``, if given."""
        return self.fields.get("SkewOpt")

    @property
    def generalized_skew(self) -> Optional[float]:
        """``GenSkew`` as a float, if given."""
        v = self.fields.get("GenSkew")
        return float(v) if v is not None else None

    @property
    def generalized_skew_se(self) -> Optional[float]:
        """``SkewSE`` as a float, if given."""
        v = self.fields.get("SkewSE")
        return float(v) if v is not None else None

    @property
    def generalized_skew_mse(self) -> Optional[float]:
        """``SkewSE`` squared: the MSE Bulletin 17C weighting uses."""
        se = self.generalized_skew_se
        return se * se if se is not None else None

    @property
    def include_urban_regulated(self) -> bool:
        """peakfq's ``Urb/Reg = Yes``."""
        return self.fields.get("Urb/Reg", "").upper() == "YES"

    @property
    def analysis_years(self) -> Optional[tuple[int, int]]:
        """First and last water years covered by a perception threshold.

        ``siteQT`` defines the analysis period this way.
        """
        if not self.thresholds:
            return None
        return min(t.start for t in self.thresholds), max(t.end for t in self.thresholds)


@dataclass
class PsfFile:
    """A parsed ``.psf`` file."""

    input_lines: List[List[str]] = field(default_factory=list)
    output_options: Dict[str, str] = field(default_factory=dict)
    stations: Dict[str, StationSpec] = field(default_factory=dict)

    @property
    def confidence_interval(self) -> Optional[float]:
        """``O ConfInterval``, if given."""
        v = self.output_options.get("ConfInterval")
        return float(v) if v is not None else None


def _comment(tokens: List[str], n_fixed: int) -> str:
    return " ".join(tokens[n_fixed:]).strip()


def parse_psf(text: str) -> PsfFile:
    """Parse the text of a ``.psf`` file.

    Parameters
    ----------
    text : str
        File contents.

    Returns
    -------
    PsfFile

    Raises
    ------
    ValueError
        On a station keyword before any ``Station`` line, a duplicated station
        ID, or a malformed repeating field.
    """
    out = PsfFile()
    current: Optional[StationSpec] = None
    for lineno, raw in enumerate(text.splitlines(), start=1):
        tokens = raw.split()
        if not tokens:
            continue
        key = tokens[0]
        if key.startswith("'"):
            continue
        if key == "Station":
            if len(tokens) < 2:
                raise ValueError(f"line {lineno}: Station with no identifier")
            sid = tokens[1]
            if sid in out.stations:
                raise ValueError(f"line {lineno}: duplicate station {sid}")
            current = StationSpec(station_id=sid)
            out.stations[sid] = current
            continue
        if key == "I":
            out.input_lines.append(tokens[1:])
            continue
        if key == "O":
            if len(tokens) >= 3:
                out.output_options[tokens[1]] = " ".join(tokens[2:])
            continue
        if current is None:
            raise ValueError(f"line {lineno}: {key!r} before any Station line")
        try:
            if key == "PCPT_Thresh":
                current.thresholds.append(
                    PerceptionThreshold(
                        int(tokens[1]),
                        int(tokens[2]),
                        float(tokens[3]),
                        float(tokens[4]),
                        _comment(tokens, 5),
                    )
                )
            elif key == "Interval":
                current.intervals.append(
                    SpecInterval(
                        int(tokens[1]), float(tokens[2]), float(tokens[3]), _comment(tokens, 4)
                    )
                )
            elif key == "Peak":
                logger.warning("Deprecated PSF keyword 'Peak' for station %s", current.station_id)
                current.peaks.append(
                    SpecPeak(int(tokens[1]), float(tokens[2]), _comment(tokens, 3))
                )
            elif key in STATION_FIELDS:
                if len(tokens) < 2:
                    raise ValueError(f"{key} with no value")
                current.fields[key] = tokens[1]
            elif key in DEPRECATED_FIELDS:
                logger.warning("Deprecated PSF field name %s will not be read", key)
            else:
                logger.warning("Unknown PSF field name: %s", key)
        except (IndexError, ValueError) as exc:
            raise ValueError(f"line {lineno}: malformed {key!r}: {raw.strip()!r}") from exc
    return out


def read_psf(path: Union[str, Path]) -> PsfFile:
    """Read and parse a ``.psf`` file from disk.

    Parameters
    ----------
    path : str or Path
        Path to the specification file.

    Returns
    -------
    PsfFile
    """
    return parse_psf(Path(path).read_text(encoding="utf-8", errors="replace"))
