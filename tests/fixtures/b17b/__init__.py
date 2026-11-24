"""Bulletin 17B worked examples and tables, as printed, for the B17B tests.

``load`` returns a fixture's JSON untouched. ``record`` builds the input record
the bulletin's procedure starts from. Values are as printed in the bulletin;
nothing here corrects them.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, NamedTuple, Optional

FIXTURES = Path(__file__).parent


class B17BRecord(NamedTuple):
    """The input of one printed example."""

    water_years: list[int]
    peaks_cfs: list[float]
    historic_peaks: list[tuple[int, float]]
    historic_period: Optional[tuple[int, int]]
    generalized_skew: float
    generalized_skew_mse: float


def load(name: str) -> dict[str, Any]:
    """Return the JSON fixture ``name`` (no extension) as printed."""
    with open(FIXTURES / f"{name}.json", encoding="utf-8") as fh:
        data: dict[str, Any] = json.load(fh)
    return data


def _number(text: Any) -> float:
    """Parse a printed number such as ``"71,500"`` or ``".302"``."""
    return float(str(text).replace(",", ""))


def _peak_rows(inp: dict[str, Any]) -> list[dict[str, Any]]:
    for key, rows in inp.items():
        if key.startswith("annual_peaks_table"):
            return list(rows)
    if "systematic_peaks" in inp:
        return list(inp["systematic_peaks"])
    raise KeyError("no annual peaks table in the fixture")


def record(name: str) -> B17BRecord:
    """Build the input record of example ``name``.

    A year printed as ``-`` or with an asterisk (Example 3's 1936 estimate)
    is left out. Example 2's 1953 peak stays in ``peaks_cfs``: the procedure
    finds it as a high outlier. Big Sandy's three historic peaks and their
    period are returned as historic information.
    """
    data = load(name)
    inp = data["input"]
    years: list[int] = []
    peaks: list[float] = []
    for row in _peak_rows(inp):
        raw = row.get("peak_cfs_as_printed", row.get("peak_cfs", row.get("Q_cfs")))
        text = str(raw).strip()
        if text == "-" or "*" in text:
            continue
        years.append(int(row["year"]))
        peaks.append(_number(text))

    historic: list[tuple[int, float]] = []
    period: Optional[tuple[int, int]] = None
    if name == "appendix6_big_sandy":
        historic = [(int(h["year"]), _number(h["Q_cfs"])) for h in inp["historic_peaks"]]
        first, last = inp["historical_period_as_printed"].split(" ")[0].split("-")
        period = (int(first), int(last))

    return B17BRecord(
        water_years=years,
        peaks_cfs=peaks,
        historic_peaks=historic,
        historic_period=period,
        generalized_skew=_number(inp["generalized_skew"]),
        generalized_skew_mse=_number(inp["generalized_skew_mse"]),
    )
