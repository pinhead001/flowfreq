"""Reader for legacy WATSTORE (``I ASCI``) peak-flow files.

PeakFQ specification files name their peak data with ``I ASCI <file>`` (this
format) or ``I RDB <file>`` (the NWIS tab-delimited format that
:meth:`flowfreq.usgs.USGSgage.download_peak_flow` parses). This module reads
the first, following ``readWATSTORE`` in ``vendor/peakfqr/R/readInputs.R``:

* fixed-width fields of 16, 8, 7, 12, 8, 4 and 4 characters -- record type
  and station ID, date, discharge, peak codes, gage height, gage-height codes,
  year of last peak;
* ``3`` records are peaks; ``Z``/``H``/``N``/``Y`` are header records
  (``N`` carries the station name);
* the water year is the calendar year for a month before October or with no
  month given, and the calendar year plus one from October on;
* run-together codes are comma-separated (``separateCodes``), with ``Bd`` and
  ``Bm`` kept whole.

The returned frame has the peak-frame columns of
:data:`flowfreq.peak_sources.PEAK_COLUMNS` plus ``site_no`` and the gage-height
fields, so it feeds :func:`flowfreq.peak_codes.peak_frame_intervals` and
:func:`flowfreq.psf_convert.convert_station` directly.

peakfq 8.1.0 itself warns that WATSTORE input is deprecated; it is supported
here because the vendored WY/MT validation data is in this format.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

#: Output columns, in order.
WATSTORE_COLUMNS = (
    "site_no",
    "water_year",
    "peak_date",
    "peak_flow_cfs",
    "qualification_code",
    "gage_height_ft",
    "gage_height_code",
)

# Field boundaries (0-based, end-exclusive), from readWATSTORE's widths
# c(16, 8, 7, 12, 8, 4, 4).
_SITE = slice(1, 16)
_DATE = slice(16, 24)
_VALUE = slice(24, 31)
_PEAK_CD = slice(31, 43)
_GAGE_HT = slice(43, 51)
_GAGE_HT_CD = slice(51, 55)
_NAME = slice(16, 64)

_HEADER_RECORDS = frozenset("ZHNY")
# separateCodes' token set: single characters, except Bd/Bm which stay whole.
_CODE_TOKEN = re.compile(r"B[dm]|\S")


def separate_codes(raw: str) -> str:
    """Comma-separate a run-together WATSTORE code field (``"7Bm"`` -> ``"7,Bm"``).

    Parameters
    ----------
    raw : str
        The code field as it appears in the file.

    Returns
    -------
    str
        Codes joined by commas; empty when the field is blank.
    """
    return ",".join(_CODE_TOKEN.findall(raw.replace(",", " ")))


def _float_or_nan(text: str) -> float:
    text = text.strip()
    if not text:
        return float("nan")
    try:
        return float(text)
    except ValueError:
        return float("nan")


def _parse_date(text: str) -> tuple[Optional[int], Optional[int], Optional[int]]:
    """``YYYYMMDD`` with the month and/or day possibly blank or ``00``."""
    text = text.strip()
    year = int(text[0:4]) if len(text) >= 4 and text[0:4].isdigit() else None
    month = int(text[4:6]) if len(text) >= 6 and text[4:6].isdigit() else None
    day = int(text[6:8]) if len(text) >= 8 and text[6:8].isdigit() else None
    return year, (month or None), (day or None)


def parse_watstore(text: str) -> pd.DataFrame:
    """Parse the text of a WATSTORE peak-flow file.

    Parameters
    ----------
    text : str
        File contents. May hold several stations.

    Returns
    -------
    pandas.DataFrame
        One row per peak record, columns :data:`WATSTORE_COLUMNS`.
        ``peak_date`` follows :func:`flowfreq.usgs._parse_peak_dt`'s
        convention for a partial date: an unknown month is January and an
        unknown day the 1st, placeholders that never change the water year.
        ``peak_flow_cfs`` is NaN for a gage-height-only record. Station names
        from ``N`` records are in ``frame.attrs["station_names"]``.

    Raises
    ------
    ValueError
        On a peak record with no readable year.
    """
    rows: List[dict] = []
    names: Dict[str, str] = {}
    skipped = 0
    for lineno, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        kind = line[0]
        if kind == "N":
            names[line[_SITE].strip()] = line[_NAME].strip()
            continue
        if kind in _HEADER_RECORDS:
            continue
        if kind != "3":
            skipped += 1
            continue
        line = line.ljust(59)
        year, month, day = _parse_date(line[_DATE])
        if year is None:
            raise ValueError(f"line {lineno}: peak record with no year: {line.rstrip()!r}")
        if month is None:
            water_year = year
        else:
            water_year = year + 1 if month >= 10 else year
        peak_date = pd.Timestamp(year=year, month=month or 1, day=day or 1)
        rows.append(
            {
                "site_no": line[_SITE].strip(),
                "water_year": water_year,
                "peak_date": peak_date,
                "peak_flow_cfs": _float_or_nan(line[_VALUE]),
                "qualification_code": separate_codes(line[_PEAK_CD]),
                "gage_height_ft": _float_or_nan(line[_GAGE_HT]),
                "gage_height_code": separate_codes(line[_GAGE_HT_CD]),
            }
        )
    if skipped:
        logger.warning("Skipped %d WATSTORE record(s) of an unrecognised type", skipped)
    frame = pd.DataFrame(rows, columns=list(WATSTORE_COLUMNS))
    frame["water_year"] = frame["water_year"].astype(int)
    frame["peak_flow_cfs"] = frame["peak_flow_cfs"].astype(np.float64)
    before = len(frame)
    frame = frame.drop_duplicates().reset_index(drop=True)
    if len(frame) != before:
        logger.warning("WATSTORE file has %d duplicate row(s); removed", before - len(frame))
    frame.attrs["station_names"] = names
    return frame


def read_watstore(path: Union[str, Path]) -> pd.DataFrame:
    """Read and parse a WATSTORE peak-flow file from disk.

    Parameters
    ----------
    path : str or Path
        Path to the file.

    Returns
    -------
    pandas.DataFrame
        See :func:`parse_watstore`.
    """
    return parse_watstore(Path(path).read_text(encoding="utf-8", errors="replace"))
