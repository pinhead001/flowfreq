"""The committed Washington Peak Region 2 analysis (docs/WA_REGION2_METHOW.md).

Reads the results ``tools/wa_region2_elevation.py`` writes; rerunning the tool is a study
step, not a test. These pin the conclusions the document draws from them.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

DATA = Path(__file__).resolve().parent.parent / "data" / "regression_study"
RESULTS = json.loads((DATA / "wa_2016_5118_r2_elevation.json").read_text(encoding="utf-8"))
BASE, ELEV = "base [1, log A, P, CAN]", "+ELEV"


def test_inputs_cover_region_2_and_the_methow():
    d = pd.read_csv(DATA / "wa_2016_5118_r2_inputs.csv", dtype={"site_no": str, "huc8": str})
    assert int(d["used"].sum()) == 88  # 89 calibration gages, 12480700 has no peaks
    assert set(RESULTS["methow"]) <= set(d["site_no"])
    assert (d.loc[d["used"], "huc8"] == "17020008").sum() == 10
    assert d["peak_years"].str.count(";").ge(9).all()


def test_elevation_is_significant_at_common_floods_and_fades():
    t = {aep: RESULTS["models"][ELEV][aep]["t"]["ELEV"] for aep in ("0.5", "0.1", "0.01")}
    assert t["0.5"] > 4.0 and t["0.1"] < 2.0 and abs(t["0.01"]) < 1.0
    base, elev = RESULTS["models"][BASE]["0.5"], RESULTS["models"][ELEV]["0.5"]
    assert elev["mev"] == pytest.approx(0.0716, abs=5e-4) and elev["mev"] < 0.85 * base["mev"]
    assert base["methow_geomean"] == pytest.approx(2.08, abs=0.01)
    assert elev["methow_geomean"] == pytest.approx(1.46, abs=0.01)


def test_the_elevation_effect_holds_without_the_methow_basin():
    outside = RESULTS["outside_methow_huc8"]
    assert outside[ELEV]["0.5"]["t"]["ELEV"] > 3.0
    assert outside[ELEV]["0.5"]["methow_geomean"] < 0.7 * outside[BASE]["0.5"]["methow_geomean"]


def test_no_methow_specific_effect_remains_after_elevation():
    both = RESULTS["models"]["+ELEV+METHOW_HUC8"]
    assert all(abs(both[a]["t"]["METHOW_HUC8"]) < 1.0 for a in ("0.5", "0.1", "0.01"))
