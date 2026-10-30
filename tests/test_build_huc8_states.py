"""tools/build_huc8_states.py: the nationwide HUC8-to-states table from the WBD service (#205).

``flowfreq/data/regional_skew_huc8_states.csv`` (``huc8,name,states``) lists only the HUC8s of a
few skew studies, so a border HU12 outside them cannot resolve its state. The tool pages the
Watershed Boundary Dataset HU8 layer (``huc8``, ``name``, ``states``; no geometry) and writes the
whole table in the same shape ``flowfreq.regional_skew._states_of`` already reads.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import pytest

from flowfreq.regional_skew import _states_of

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "tools"))

import build_huc8_states as bhs  # noqa: E402

FIXTURE = REPO_ROOT / "tests" / "fixtures" / "wbd_huc8_page.json"


def _page() -> Dict[str, Any]:
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


def test_rows_from_features_keep_non_us_codes() -> None:
    rows = bhs.rows_from_features(_page()["features"])
    assert rows[0] == ("09040001", "St. Marys", "CN,MT")
    assert len(rows) == 3


def test_rows_skip_a_feature_without_a_huc8() -> None:
    feats = _page()["features"] + [{"attributes": {"huc8": "", "name": "x", "states": "WA"}}]
    assert len(bhs.rows_from_features(feats)) == 3


def test_write_table_sorted_unique_and_quoted(tmp_path: Path) -> None:
    rows = bhs.rows_from_features(_page()["features"])
    rows.append(("09040001", "St. Marys (duplicate)", "MT"))
    out = tmp_path / "huc8.csv"
    bhs.write_table(rows, out)
    text = out.read_text(encoding="utf-8")
    assert text.splitlines()[0] == "huc8,name,states"
    assert '"CN,MT"' in text
    with out.open(encoding="utf-8", newline="") as fh:
        body = list(csv.reader(fh))[1:]
    assert [r[0] for r in body] == ["09040001", "10020002", "17110019"]
    assert body[0] == ["09040001", "St. Marys", "CN,MT"]  # first occurrence kept


def test_pages_follow_exceeded_transfer_limit() -> None:
    calls: List[Dict[str, Any]] = []

    def fake_get(url: str, params: Dict[str, Any]) -> Dict[str, Any]:
        calls.append(dict(params))
        page = _page()  # three features: fewer than the page size asked for
        page["exceededTransferLimit"] = params["resultOffset"] == 0
        return page

    pages = list(bhs.fetch_huc8_pages(fake_get, page_size=1000))
    assert len(pages) == 2
    # The offset advances by the features actually returned, not by page_size: an ArcGIS
    # service caps a page at its maxRecordCount, and stepping by page_size would skip rows.
    assert [c["resultOffset"] for c in calls] == [0, 3]
    assert calls[0]["outFields"] == "huc8,name,states"
    assert calls[0]["returnGeometry"] in ("false", False)
    assert calls[0]["resultRecordCount"] == 1000
    assert calls[0]["orderByFields"] == "huc8"  # resultOffset paging is stable only with an order


def test_error_payload_raises() -> None:
    def fake_get(url: str, params: Dict[str, Any]) -> Dict[str, Any]:
        return {"error": {"code": 400, "message": "Unable to complete operation."}}

    with pytest.raises(RuntimeError, match="Unable to complete"):
        list(bhs.fetch_huc8_pages(fake_get, page_size=10))


def test_refuses_to_shrink_the_packaged_table(tmp_path: Path) -> None:
    existing = tmp_path / "huc8.csv"
    bhs.write_table(bhs.rows_from_features(_page()["features"]), existing)
    with pytest.raises(RuntimeError, match="3"):
        bhs.ensure_not_shrunk([("09040001", "St. Marys", "CN,MT")], existing)
    with pytest.raises(RuntimeError):
        bhs.ensure_not_shrunk([], existing)
    bhs.ensure_not_shrunk(bhs.rows_from_features(_page()["features"]), existing)  # same size: ok
    bhs.ensure_not_shrunk([("1", "a", "WA")], tmp_path / "absent.csv")  # no file yet: ok


def test_pages_stop_on_an_empty_page() -> None:
    def fake_get(url: str, params: Dict[str, Any]) -> Dict[str, Any]:
        return {"features": [], "exceededTransferLimit": True}

    assert len(list(bhs.fetch_huc8_pages(fake_get, page_size=10))) == 1


def test_written_table_is_read_by_states_of(tmp_path: Path) -> None:
    out = tmp_path / "huc8.csv"
    bhs.write_table(bhs.rows_from_features(_page()["features"]), out)
    assert _states_of("09040001", out) == {"CN", "MT"}
    assert _states_of("171100", out) == {"WA"}  # the reader prefix-matches a shorter HUC


def test_main_default_out_is_the_packaged_table() -> None:
    assert bhs.DEFAULT_OUT == REPO_ROOT / "flowfreq" / "data" / "regional_skew_huc8_states.csv"
    with pytest.raises(SystemExit):
        bhs.main(["--help"])
