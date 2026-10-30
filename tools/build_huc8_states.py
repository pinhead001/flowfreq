"""Build the nationwide HUC8-to-states table from the Watershed Boundary Dataset.

Pages the WBD HU8 layer (``huc8``, ``name``, ``states``; no geometry) served by the USGS
National Map and writes the whole table in the shape of
``flowfreq/data/regional_skew_huc8_states.csv`` (``huc8,name,states``). ``states`` is a
comma-joined list of two-letter codes kept verbatim, so ``"CN,MT"`` stays ``"CN,MT"``.

Service: ``https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer/4``.

Usage::

    python tools/build_huc8_states.py [--out PATH] [--page-size N]
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, List, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
WBD_HU8_QUERY = "https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer/4/query"
DEFAULT_OUT = ROOT / "flowfreq" / "data" / "regional_skew_huc8_states.csv"

logger = logging.getLogger(__name__)

Row = Tuple[str, str, str]


def fetch_huc8_pages(
    get: Callable[[str, Dict[str, Any]], Dict[str, Any]], *, page_size: int = 1000
) -> Iterator[Dict[str, Any]]:
    """Yield every page payload of the WBD HU8 layer.

    Parameters
    ----------
    get : callable
        ``get(url, params)`` returning the parsed JSON of one page.
    page_size : int, optional
        Records requested per page. The offset advances by the features actually
        returned, since the service may cap a page below this.

    Yields
    ------
    dict
        One page payload.

    Raises
    ------
    RuntimeError
        If the service answers with an error payload.
    """
    offset = 0
    while True:
        params: Dict[str, Any] = {
            "where": "1=1",
            "outFields": "huc8,name,states",
            "returnGeometry": "false",
            "f": "json",
            "orderByFields": "huc8",
            "resultRecordCount": page_size,
            "resultOffset": offset,
        }
        page = get(WBD_HU8_QUERY, params)
        if "error" in page:
            err = page["error"]
            message = err.get("message", err) if isinstance(err, dict) else err
            raise RuntimeError(f"WBD service error: {message}")
        yield page
        features = page.get("features") or []
        if not features or page.get("exceededTransferLimit") is not True:
            return
        offset += len(features)


def rows_from_features(features: Sequence[Dict[str, Any]]) -> List[Row]:
    """Return ``(huc8, name, states)`` rows from ArcGIS features.

    Parameters
    ----------
    features : sequence of dict
        Features carrying an ``attributes`` mapping.

    Returns
    -------
    list of tuple
        One row per feature with a non-empty ``huc8``.
    """
    rows: List[Row] = []
    for feature in features:
        attrs = feature.get("attributes") or {}
        huc8 = str(attrs.get("huc8") or "").strip()
        if not huc8:
            continue
        rows.append((huc8, str(attrs.get("name") or ""), str(attrs.get("states") or "")))
    return rows


def write_table(rows: Sequence[Row], path: Path) -> None:
    """Write the table sorted by ``huc8``, keeping the first of duplicate HUC8s.

    Parameters
    ----------
    rows : sequence of tuple
        ``(huc8, name, states)`` rows.
    path : pathlib.Path
        Destination CSV (UTF-8, LF line endings).
    """
    unique: Dict[str, Row] = {}
    for row in rows:
        unique.setdefault(row[0], row)
    with Path(path).open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["huc8", "name", "states"])
        for key in sorted(unique):
            writer.writerow(unique[key])


def ensure_not_shrunk(rows: Sequence[Row], path: Path) -> None:
    """Refuse to replace an existing table with a smaller one.

    Parameters
    ----------
    rows : sequence of tuple
        Rows about to be written.
    path : pathlib.Path
        Existing table, if any.

    Raises
    ------
    RuntimeError
        If the unique HUC8 count of ``rows`` is below that of the existing file.
    """
    path = Path(path)
    if not path.exists():
        return
    with path.open(encoding="utf-8", newline="") as fh:
        existing = len({r[0] for r in list(csv.reader(fh))[1:] if r})
    new = len({r[0] for r in rows})
    if new < existing:
        raise RuntimeError(
            f"refusing to shrink {path}: new table has {new} HUC8s, existing has {existing}"
        )


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Fetch the WBD HU8 layer and write the table.

    Parameters
    ----------
    argv : sequence of str, optional
        Command-line arguments (default ``sys.argv[1:]``).

    Returns
    -------
    int
        Exit status, 0 on success.
    """
    import requests

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--page-size", type=int, default=1000)
    args = parser.parse_args(argv)

    def get(url: str, params: Dict[str, Any]) -> Dict[str, Any]:
        resp = requests.get(url, params=params, timeout=300)
        resp.raise_for_status()
        return resp.json()

    rows: List[Row] = []
    for page in fetch_huc8_pages(get, page_size=args.page_size):
        rows.extend(rows_from_features(page.get("features") or []))
    ensure_not_shrunk(rows, args.out)
    write_table(rows, args.out)
    logger.info("wrote %d HUC8 rows to %s", len({r[0] for r in rows}), args.out)
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
