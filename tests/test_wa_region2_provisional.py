"""The provisional Washington set: Peak Region 2 refit with mean basin elevation
(tools/wa_region2_provisional.py; docs/WA_REGION2_METHOW.md)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

from flowfreq.regression.equations import evaluate
from flowfreq.regression.library import load_state, provisional_states

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "tools"))

import wa_region2_elevation as r2  # noqa: E402
import wa_region2_provisional as prov  # noqa: E402

PUB = load_state("WA")
LIB = load_state("WA", allow_provisional=True)


def test_opt_in_only() -> None:
    assert PUB.status == "verified"
    assert LIB.status == "provisional" and LIB.advisory.startswith("PROVISIONAL")
    assert "WA" in provisional_states()


def test_only_region_2_changes_and_it_adds_elev_at_every_aep() -> None:
    assert sorted(LIB.regions) == sorted(PUB.regions)
    for code in PUB.regions:
        assert LIB.aeps(code) == PUB.aeps(code)
        for aep in PUB.aeps(code):
            new, old = LIB.equation(code, aep), PUB.equation(code, aep)
            if code == "2":
                assert [v.code for v in new.variables] == [
                    "DRNAREA",
                    "PRECPRIS10",
                    "CANOPY_PCT",
                    "ELEV",
                ]
                assert new.covariance is not None and np.asarray(new.covariance).shape == (5, 5)
                assert "NOT a USGS study" in new.citation.publication
            else:
                assert new == old


def test_elevation_effect_tapers_to_zero() -> None:
    coef = [LIB.equation("2", a).coefficients[3] for a in LIB.aeps("2")]
    assert coef[0] == pytest.approx(1.6e-4, rel=0.05)  # 50% AEP
    assert all(a >= b for a, b in zip(coef, coef[1:]))
    assert abs(coef[-1]) < 1e-5


def test_frequency_curve_increases_at_every_region_2_gage() -> None:
    d = r2.load_frame()
    d = d[d["ELEV"].notna()]
    for _, row in d.iterrows():
        ch = {
            "DRNAREA": row.DRNAREA,
            "PRECPRIS10": row.PRECIP,
            "CANOPY_PCT": row.CANOPY_PCT,
            "ELEV": row.ELEV,
        }
        q = [
            evaluate(LIB.equation("2", a), ch, allow_extrapolation=True).flow_cfs
            for a in LIB.aeps("2")
        ]
        assert all(b > a for a, b in zip(q, q[1:])), row.site_no


def test_it_raises_the_methow_at_common_floods() -> None:
    d = r2.load_frame().set_index("site_no")
    row = d.loc["12448500"]  # Methow at Winthrop
    ch = {"DRNAREA": row.DRNAREA, "PRECPRIS10": row.PRECIP, "CANOPY_PCT": row.CANOPY_PCT}
    new = evaluate(LIB.equation("2", 0.5), {**ch, "ELEV": row.ELEV}, allow_extrapolation=True)
    old = evaluate(PUB.equation("2", 0.5), ch, allow_extrapolation=True)
    assert new.flow_cfs / old.flow_cfs > 1.2
    assert new.interval is not None


def test_committed_file_is_what_the_tool_builds() -> None:
    # Numerically, not byte for byte: the GLS fit may differ in the last bits by platform.
    built = json.loads(json.dumps(prov.build()))
    committed = json.loads(prov.OUT.read_text(encoding="utf-8"))
    assert built["notes"] == committed["notes"] and built["status"] == committed["status"]
    assert len(built["equations"]) == len(committed["equations"])
    for b, c in zip(built["equations"], committed["equations"]):
        assert (b["region_code"], b["aep"]) == (c["region_code"], c["aep"])
        assert b["intercept"] == pytest.approx(c["intercept"], rel=1e-8, abs=1e-12)
        assert b["coefficients"] == pytest.approx(c["coefficients"], rel=1e-8, abs=1e-12)
        assert b["variables"] == c["variables"]
