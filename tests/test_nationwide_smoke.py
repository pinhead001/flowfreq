"""Offline tests for tools/nationwide_smoke.py: site selection and the report writer."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "tools"))

import nationwide_smoke as ns  # noqa: E402


def _catalog() -> pd.DataFrame:
    rows = [
        ("01000001", "TN", 40, "reference"),
        ("01000002", "TN", 90, "regulated"),
        ("01000003", "TN", 60, "reference"),
        ("01000004", "TN", 80, "unknown"),
        ("02000001", "OH", 30, "regulated"),
        ("02000002", "OH", 20, "urban"),
        ("02000003", "OH", 25, "unknown"),
        ("03000001", "DE", 50, "regulated"),
    ]
    return pd.DataFrame(rows, columns=["site_no", "state", "n_peaks", "regulation_class"])


def test_select_prefers_longest_reference():
    assert ns.select_site(_catalog(), "TN")["site_no"] == "01000003"


def test_select_falls_back_to_unclassified_before_urban():
    assert ns.select_site(_catalog(), "OH")["site_no"] == "02000003"


def test_select_none_when_only_regulated_or_absent():
    assert ns.select_site(_catalog(), "DE") is None
    assert ns.select_site(_catalog(), "ZZ") is None


def _row(state, ffa=True, ung=True, err=None):
    return {
        "state": state,
        "site_no": "01000003",
        "n_peaks": 60,
        "ffa_ok": ffa,
        "ffa_q1pct": 1234.5 if ffa else None,
        "ungaged_ok": ung,
        "ungaged_q1pct": 999.0 if ung else None,
        "library_status": "verified" if ung else None,
        "error": err,
    }


def test_write_report(tmp_path):
    rows = [
        _row("TN"),
        _row("DC", ffa=False, ung=False, err="ungaged: EquationsUnavailable: no library"),
        _row("AK", ung=False, err="ungaged: UnsupportedRegionError: not delineated"),
    ]
    out = tmp_path / "r.md"
    ns.write_report(rows, out, run_date="2026-10-25")
    text = out.read_text(encoding="utf-8")
    assert "2026-10-25" in text
    assert "Known gaps" in text
    assert text.count("\n| TN |") == 1 and "1,234" in text
    assert "Gaged fits succeeded: 2 of 3" in text
    assert "Ungaged estimates succeeded: 1 of 3" in text
    assert "ungaged EquationsUnavailable" in text
    assert "ungaged UnsupportedRegionError" in text
