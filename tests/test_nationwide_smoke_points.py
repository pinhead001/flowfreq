"""tools/nationwide_smoke.py ``--points``: per-state coordinate overrides for the ungaged
estimate (#210).

StreamStats exposes no snap distance, and the smoke script cannot move a point downstream,
so the 2026-10-25 run lost AR, GA, ND, PA, PR and TX to unsnappable gage coordinates. A
points CSV (``state,lat,lon,note``) supplies a hand-picked pour point per state; the gaged
fit still uses the gage, and the report says which estimates used an override.
"""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "tools"))

import nationwide_smoke as ns  # noqa: E402


def _csv(tmp_path: Path, text: str) -> Path:
    p = tmp_path / "points.csv"
    p.write_text(text, encoding="utf-8")
    return p


def test_load_points(tmp_path: Path) -> None:
    p = _csv(
        tmp_path, "state,lat,lon,note\ntx,30.25,-97.75,moved 200 m downstream\nPA,40.5,-77.1,\n"
    )
    pts = ns.load_points(p)
    assert set(pts) == {"TX", "PA"}
    assert pts["TX"] == {"lat": 30.25, "lon": -97.75, "note": "moved 200 m downstream"}
    assert pts["PA"]["note"] == ""


@pytest.mark.parametrize(
    "text,match",
    [
        ("state,lat\nTX,30.0\n", "lon"),
        ("state,lat,lon,note\nTX,95.0,-97.0,x\n", "TX"),
        ("state,lat,lon,note\nTX,30.0,-197.0,x\n", "TX"),
        ("state,lat,lon,note\nTX,30.0,-97.0,a\ntx,31.0,-97.0,b\n", "TX"),
        ("state,lat,lon,note\nTX,abc,-97.0,x\n", "TX"),
    ],
)
def test_load_points_refuses_bad_rows(tmp_path: Path, text: str, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        ns.load_points(_csv(tmp_path, text))


@pytest.fixture
def fakes(monkeypatch: pytest.MonkeyPatch) -> dict:
    calls: dict = {}
    import flowfreq
    import flowfreq.usgs
    import flowfreq.workflow

    class FakeGage:
        def __init__(self, site: str) -> None:
            calls["gage"] = site

        def download_peak_flow(self) -> pd.DataFrame:
            return pd.DataFrame(
                {"peak_flow_cfs": [100.0, 200.0], "water_year": [2000, 2001],
                 "qualification_code": ["", ""]}
            )  # fmt: skip

    def fake_ffa(*a: object, **k: object) -> dict:
        return {"quantile_df": pd.DataFrame({"AEP (%)": [0.01], "Flow (cfs)": [5000.0]})}

    def fake_estimate(lat: float, lon: float, **k: object) -> SimpleNamespace:
        calls["point"] = (lat, lon)
        q = pd.DataFrame({"aep": [0.01], "flow_cfs": [4000.0]})
        return SimpleNamespace(quantiles=q, library_status="verified")

    monkeypatch.setattr(flowfreq.usgs, "USGSgage", FakeGage)
    monkeypatch.setattr(flowfreq.workflow, "run_ffa", fake_ffa)
    monkeypatch.setattr(flowfreq, "estimate_at_site", fake_estimate)
    return calls


ROW = pd.Series(
    {"state": "TX", "site_no": "08158000", "n_peaks": 80, "latitude": 30.24, "longitude": -97.69}
)


def test_smoke_one_uses_gage_coordinates_by_default(fakes: dict) -> None:
    out = ns.smoke_one(ROW)
    assert fakes["point"] == (30.24, -97.69) and fakes["gage"] == "08158000"
    assert out["ungaged_point"] == "gage" and out["ungaged_ok"] and out["ffa_ok"]


def test_smoke_one_with_override_point(fakes: dict) -> None:
    out = ns.smoke_one(ROW, point={"lat": 30.25, "lon": -97.75, "note": "moved downstream"})
    assert fakes["point"] == (30.25, -97.75)
    assert fakes["gage"] == "08158000"  # the gaged fit still uses the gage
    assert out["ungaged_point"] == "override: moved downstream"
    assert out["ungaged_ok"]


def test_report_marks_override_points(tmp_path: Path) -> None:
    base = {
        "site_no": "01000003", "n_peaks": 60, "ffa_ok": True, "ffa_q1pct": 1.0,
        "ungaged_ok": True, "ungaged_q1pct": 2.0, "library_status": "verified", "error": None,
    }  # fmt: skip
    rows = [
        {**base, "state": "TN"},  # an older cached row without ungaged_point
        {**base, "state": "TX", "ungaged_point": "override: moved downstream"},
    ]
    out = tmp_path / "r.md"
    ns.write_report(rows, out, run_date="2026-10-26")
    text = out.read_text(encoding="utf-8")
    assert "| Point |" in text
    tx = next(line for line in text.splitlines() if line.startswith("| TX |"))
    tn = next(line for line in text.splitlines() if line.startswith("| TN |"))
    assert "override: moved downstream" in tx and "| gage |" in tn


def test_main_accepts_points_option() -> None:
    with pytest.raises(SystemExit) as e:
        ns.main(["--help"])
    assert e.value.code == 0


def test_points_option_is_parsed(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    seen: dict = {}
    monkeypatch.setattr(ns, "load_points", lambda p: seen.setdefault("path", p) and {})
    import flowfreq.catalog

    monkeypatch.setattr(flowfreq.catalog, "load_catalog", lambda: pd.DataFrame(
        columns=["state", "site_no", "n_peaks", "regulation_class"]))  # fmt: skip
    p = _csv(tmp_path, "state,lat,lon,note\n")
    ns.main(["--points", str(p), "--states", "ZZ", "--out", str(tmp_path / "o.md"),
             "--cache", str(tmp_path / "c")])  # fmt: skip
    assert seen["path"] == p
