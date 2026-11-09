"""Offline tests for tools/border_qa.py: adjacency, proxies, row reasons and the summary."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "tools"))

import border_qa as bq  # noqa: E402


def test_each_border_appears_once_and_is_symmetric():
    pairs = [tuple(sorted(p)) for p in bq.ADJACENT]
    assert len(pairs) == len(set(pairs))
    assert "WV" in bq.neighbors("PA") and "PA" in bq.neighbors("WV")
    assert bq.neighbors("HI") == []


def test_proxies_fill_only_what_is_missing_and_are_listed():
    values, subs = bq.fill_proxies({"CSL10_85": 12.0, "BSHAPELFP": 4.0}, 40.0, -90.0)
    assert values["CSL1085LFP"] == 12.0 and values["BSHAPE"] == 4.0
    assert "CSL1085LFP=CSL10_85" in subs and "BSHAPE=BSHAPELFP" in subs
    assert values["CSL10_85"] == 12.0  # an existing value is never replaced
    assert values["LAT_OUT"] == 40.0  # from wave_qa's own proxies


def _row(site, home, nb, aep, hq, nq, at_site=None, why=""):
    return {
        "site_no": site,
        "home": home,
        "neighbor": nb,
        "aep": aep,
        "home_q": hq,
        "neighbor_q": nq,
        "at_site_q": at_site,
        "why_not": why,
    }


def test_summary_counts_factors_and_flags_systematic_pairs():
    rows = pd.DataFrame(
        [
            _row("1", "MO", "KY", 0.01, 1000.0, 4000.0, 1100.0),
            _row("2", "MO", "KY", 0.01, 1000.0, 3500.0, 900.0),
            _row("3", "MO", "KY", 0.01, 1000.0, 3000.0, 2000.0),
            _row("4", "PA", "WV", 0.01, 1000.0, 1100.0),
            _row("5", "PA", "WV", 0.01, None, None, why="neighbour: no region located"),
            _row("1", "MO", "KY", 0.5, 100.0, 120.0),
        ]
    )
    rows["clean"] = True
    s = bq.summarize(rows)
    one = next(b for b in s["by_aep"] if b["aep"] == 0.01)
    assert (one["n"], one["within_1.5x"], one["2_to_3x"], one["3x_or_more"]) == (4, 1, 0, 3)
    mo_ky = next(p for p in s["pairs"] if (p["home"], p["neighbor"]) == ("MO", "KY"))
    assert mo_ky["systematic"] and mo_ky["median_neighbor_over_home"] == 3.5
    pa_wv = next(p for p in s["pairs"] if p["home"] == "PA")
    assert not pa_wv["systematic"]  # one gage is not a pattern
    assert (s["home_nearer"], s["with_fit"]) == (2, 3)
    assert s["why_not"] == {"neighbour: no region located": 1}
    assert math.isclose(10 ** abs(math.log10(3.5)), 3.5)


def test_default_pick_avoids_urban_and_drainage_area_only_sets():
    cands = [
        ("GC1486", "Urban (SIR 2010-5073)", 2010),
        ("GC1564", "Region 1, drainage area only", 2015),
        ("GC1485", "Region 3, rural", 2014),
    ]
    assert bq.default_pick(cands) == "GC1485"
    assert bq.default_pick(cands[:2]) == "GC1564"


def test_clean_excludes_extrapolation_and_urban_default_picks():
    home = {"home_est": {"0.01": [1000.0, False]}}
    assert bq._clean(home, {"0.01": [1200.0, False], "regions": ["GC1:1"]}, "0.01")
    assert not bq._clean(home, {"0.01": [1200.0, True], "regions": ["GC1:1"]}, "0.01")
    urban = {"0.01": [1200.0, False], "regions": ["GC1486:Urban"], "picked": "GC1486"}
    assert not bq._clean(home, urban, "0.01")
    assert not bq._clean({"home_est": {"0.01": [1000.0, True]}}, {"0.01": [1, False]}, "0.01")


def test_reasons_a_row_is_not_compared():
    assert bq._why_not({"error": "UnsnappablePointError: x"}, {}, "0.01").startswith("delineation")
    assert bq._why_not({"mis_snapped": True}, {}, "0.01").startswith("delineated area")
    missing = {"error": "EquationsUnavailable: every matched region needs a basin characteristic"}
    assert bq._why_not({"home_est": {"0.01": [1, False]}}, missing, "0.01") == (
        "neighbour: characteristic not computed in the home region"
    )
    assert bq._why_not({"home_est": {"0.01": [1, False]}}, {"0.01": [2, False]}, "0.01") == ""


def test_a_state_without_the_50_percent_aep_still_compares_at_10_and_1():
    from flowfreq.regression.library import load_state
    from flowfreq.streamstats import RegressionRegionLocation

    loc = [
        RegressionRegionLocation(
            region_id=1, code="GC1009", name="5", percent_weight=100.0, area_sq_mi=1.0
        )
    ]
    values = {"CONTDA": 200.0, "CSL10_85": 4.0, "HIGHREG": 1009.0}
    out = bq._assemble(load_state("IN"), loc, values)
    assert "0.5" not in out and out["0.1"][0] > 0 and out["0.01"][0] > out["0.1"][0]
