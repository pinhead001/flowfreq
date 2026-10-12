"""Read one worksheet of an .xlsx file with the standard library only.

flowfreq does not depend on openpyxl, and the study builders in ``tools/`` only
need cell values: shared strings, inline strings and numbers, laid out by cell
reference. Formulas are read as their cached values. Dates are not converted.
"""

from __future__ import annotations

import re
import zipfile
from pathlib import Path
from typing import Dict, List
from xml.etree import ElementTree as ET

_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_REF = re.compile(r"([A-Z]+)(\d+)")


def _col_index(letters: str) -> int:
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def _shared_strings(z: zipfile.ZipFile) -> List[str]:
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    root = ET.fromstring(z.read("xl/sharedStrings.xml"))
    return ["".join(t.text or "" for t in si.iter(f"{{{_NS['m']}}}t")) for si in root]


def _sheet_paths(z: zipfile.ZipFile) -> List[str]:
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    target: Dict[str, str] = {r.get("Id", ""): r.get("Target", "") for r in rels}
    out = []
    for s in wb.find("m:sheets", _NS) or []:
        t = target[s.get(f"{{{_REL_NS}}}id", "")]
        out.append(t.lstrip("/") if t.startswith("/xl/") else f"xl/{t}")
    return out


def read_sheet(path: Path, sheet_index: int = 0) -> List[List[str]]:
    """All rows of one worksheet as strings ("" for an empty cell), rectangular.

    Parameters
    ----------
    path : Path
        The .xlsx file.
    sheet_index : int
        Zero-based position of the sheet in the workbook.

    Returns
    -------
    list of list of str
        Row-major, padded to the widest row; row and column positions match
        the sheet's cell references (row 1 is index 0).
    """
    with zipfile.ZipFile(path) as z:
        strings = _shared_strings(z)
        root = ET.fromstring(z.read(_sheet_paths(z)[sheet_index]))
    cells: Dict[int, Dict[int, str]] = {}
    for c in root.iter(f"{{{_NS['m']}}}c"):
        m = _REF.fullmatch(c.get("r", ""))
        if m is None:
            continue
        row, col = int(m.group(2)) - 1, _col_index(m.group(1))
        kind = c.get("t")
        if kind == "inlineStr":
            value = "".join(t.text or "" for t in c.iter(f"{{{_NS['m']}}}t"))
        else:
            v = c.find("m:v", _NS)
            raw = v.text if v is not None and v.text is not None else ""
            value = strings[int(raw)] if kind == "s" and raw else raw
        cells.setdefault(row, {})[col] = value
    if not cells:
        return []
    width = max(max(r) for r in cells.values()) + 1
    return [[cells.get(r, {}).get(k, "") for k in range(width)] for r in range(max(cells) + 1)]
