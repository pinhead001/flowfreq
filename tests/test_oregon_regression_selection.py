"""Multi-AEP variable selection in tools/oregon_regression_fit.py (``--select multi``).

Synthetic regions, no network and no study data: the quantiles are built from
DRNAREA and PRECIP, and JANMAXT is noise.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import numpy as np
import pandas as pd
import pytest

TOOL = Path(__file__).resolve().parents[1] / "tools" / "oregon_regression_fit.py"


@pytest.fixture()
def tool(monkeypatch) -> ModuleType:
    spec = importlib.util.spec_from_file_location("oregon_regression_fit", TOOL)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    keep = ("log10 DRNAREA", "log10 PRECIP", "JANMAXT")
    monkeypatch.setattr(mod, "CANDIDATES", {k: mod.CANDIDATES[k] for k in keep})
    return mod


def region(seed: int = 1, m: int = 40, outlier: bool = False) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    area = 10 ** rng.uniform(0.5, 3.0, m)
    precip = rng.uniform(20.0, 120.0, m)
    if outlier:
        # one basin far outside the others' PRECIP: it carries the PRECIP slope
        precip[0] = 600.0
    df = pd.DataFrame(
        {
            "site_no": [f"{i:08d}" for i in range(m)],
            "DRNAREA": area,
            "PRECIP": precip,
            "JANMAXT": rng.uniform(30.0, 50.0, m),
            "sd_log": rng.uniform(0.2, 0.3, m),
            "n_peaks": rng.integers(20, 80, m).astype(float),
            "skew_weighted": np.full(m, -0.1),
            "skew_station": np.full(m, -0.1),
            "first_wy": np.full(m, 1950.0),
            "last_wy": np.full(m, 2020.0),
            "lat_cent": rng.uniform(42.0, 46.0, m),
            "lon_cent": rng.uniform(-124.0, -121.0, m),
        }
    )
    noise = rng.normal(0.0, 0.08, m)
    for aep, k in zip((0.5, 0.1, 0.01, 0.002), (0.0, 1.3, 2.3, 2.9)):
        logq = 1.0 + 0.8 * np.log10(area) + 0.9 * np.log10(precip) + 0.25 * k + noise
        df[f"q_{aep}"] = 10**logq
    return df


def ev(t: list, lev: float = 0.2, aeps=(0.5, 0.1, 0.01, 0.002)) -> dict:
    return {a: {"avp": 0.01, "sep_percent": 23.0, "t": np.array(t[i]), "max_leverage": lev}
            for i, a in enumerate(aeps)}  # fmt: skip


def test_criterion_weighs_each_aep_by_its_own_base(tool):
    base = {0.5: {"avp": 0.02}, 0.01: {"avp": 0.08}}
    ev_ = {0.5: {"avp": 0.01}, 0.01: {"avp": 0.06}}
    assert tool.criterion(ev_, base) == pytest.approx((0.5 + 0.75) / 2)


def test_rejection_needs_t_at_one_percent_and_at_most_aeps(tool):
    assert tool.rejection(ev([[3.0], [2.5], [2.2], [2.1]])) is None
    # significant at 3 of 4 AEPs, including 0.01: passes
    assert tool.rejection(ev([[1.0], [2.5], [2.2], [2.1]])) is None
    # not significant at 0.01
    assert "0.01" in tool.rejection(ev([[3.0], [2.5], [1.5], [2.1]]))
    # significant at 0.01 but only at half the AEPs
    assert "half" in tool.rejection(ev([[1.0], [1.5], [2.2], [2.1]]))


def test_rejection_guards_leverage(tool):
    assert tool.rejection(ev([[3.0]] * 4, lev=0.67)).startswith("leverage 0.67")
    assert tool.rejection(ev([[3.0]] * 4, lev=0.67), max_leverage=0.7) is None


def test_select_multi_finds_the_signal_and_not_the_noise(tool):
    names, log = tool.select_multi(region())
    assert names == ["log10 DRNAREA", "log10 PRECIP"]
    assert log[0]["set"] == names
    assert log[0]["criterion"] < 1.0


def test_select_multi_drops_a_set_carried_by_one_high_leverage_site(tool):
    df = region(outlier=True)
    aeps = tool.SELECT_AEPS

    def lev(names: list) -> float:
        return max(e["max_leverage"] for e in tool.evaluate_set(df, names, aeps).values())

    lev_base, lev_full = lev(["log10 DRNAREA"]), lev(["log10 DRNAREA", "log10 PRECIP"])
    assert lev_full > lev_base + 0.1  # the outlier basin takes the PRECIP leverage
    names, log = tool.select_multi(df, max_leverage=(lev_base + lev_full) / 2)
    assert names == ["log10 DRNAREA"]
    assert any(str(e.get("rejected", "")).startswith("leverage") for e in log)
    names, _ = tool.select_multi(df, max_leverage=1.0)
    assert "log10 PRECIP" in names


def test_no_worse_compares_sep_at_every_guard_aep(tool):
    df = region()
    ok, info = tool.no_worse(df, ["log10 DRNAREA", "log10 PRECIP"], ["log10 DRNAREA"])
    assert ok
    assert set(info["sep_percent_new"]) == {"0.5", "0.01", "0.002"}
    ok, _ = tool.no_worse(df, ["log10 DRNAREA"], ["log10 DRNAREA", "log10 PRECIP"])
    assert not ok


def test_choose_keeps_the_stored_set_unless_the_new_one_is_no_worse(tool, monkeypatch):
    df = region()
    monkeypatch.setattr(tool, "stored_variables", lambda r: ["log10 DRNAREA"])
    names, log = tool.choose(df, "X", "multi")
    assert names == ["log10 DRNAREA", "log10 PRECIP"]
    assert log[-1]["replaces"] == ["log10 DRNAREA"]
    # the stored set wins when the selection can only find something worse
    monkeypatch.setattr(tool, "select_multi", lambda d: (["log10 DRNAREA"], []))
    monkeypatch.setattr(tool, "stored_variables", lambda r: ["log10 DRNAREA", "log10 PRECIP"])
    names, log = tool.choose(df, "X", "multi")
    assert names == ["log10 DRNAREA", "log10 PRECIP"]
    assert log[-1]["kept_stored"] == names


def test_choose_replaces_a_stored_set_that_fails_the_leverage_guard(tool, monkeypatch):
    # PRECIP's slope comes from one basin far outside the others: the stored set
    # fits better but breaks the leverage limit, so the selection replaces it.
    df = region(outlier=True)
    full = tool.evaluate_set(df, ["log10 DRNAREA", "log10 PRECIP"], tool.SELECT_AEPS)
    lev_full = max(e["max_leverage"] for e in full.values())
    monkeypatch.setattr(tool, "MAX_LEVERAGE", lev_full - 0.05)
    monkeypatch.setattr(tool, "select_multi", lambda d: (["log10 DRNAREA"], []))
    monkeypatch.setattr(tool, "stored_variables", lambda r: ["log10 DRNAREA", "log10 PRECIP"])
    names, log = tool.choose(df, "X", "multi")
    assert names == ["log10 DRNAREA"]
    assert log[-1]["replaces_failing"] == ["log10 DRNAREA", "log10 PRECIP"]
    assert log[-1]["why"].startswith("leverage")


def test_choose_adopts_the_selection_for_a_region_not_yet_stored(tool, monkeypatch):
    def missing(r: str) -> list:
        raise KeyError(r)

    monkeypatch.setattr(tool, "stored_variables", missing)
    names, log = tool.choose(region(), "NEW", "multi")
    assert names == ["log10 DRNAREA", "log10 PRECIP"]
    assert log[-1]["stored"] is None


def test_run_rejects_an_unknown_select_mode(tool):
    with pytest.raises(ValueError, match="select_mode"):
        tool.run(False, select_mode="best")


def _with_lakes(lakes: list) -> pd.DataFrame:
    g = region(m=len(lakes))
    g["LAKEAREA"] = lakes
    g["region"] = "W"
    return g


def test_lake_screen_lists_gages_above_the_limit_largest_first(tool):
    g = _with_lakes([0.0, 32.1, 4.2, 6.7] + [0.0] * 36)
    got = tool.lake_screen(g, 5.0)
    assert [r["lake_percent"] for r in got] == [32.1, 6.7]
    assert got[0]["site_no"] == g["site_no"][1] and got[0]["region"] == "W"
    assert tool.lake_screen(g, None) == [] and tool.lake_screen(g, 0.0) == []


def test_stored_equations_carry_the_lake_limit(tool):
    from flowfreq.regression import OutOfRangeError, RegressionEquation, evaluate

    g = _with_lakes([0.0, 32.1] + [1.0] * 38)
    names = ["log10 DRNAREA", "log10 PRECIP"]
    res = tool.fit(g, names, 0.01)
    d = tool.equation_json("W", g, names, 0.01, res, lake_limit=5.0)
    assert d["applicability"] == [
        {"code": "LAKEAREA_NHD", "transform": "identity", "units": "percent",
         "minimum": 0.0, "maximum": 5.0}
    ]  # fmt: skip
    eq = RegressionEquation.from_dict(d)
    chars = {"DRNAREA": 50.0, "PRECIP_PRISM9120": 60.0}
    assert evaluate(eq, {**chars, "LAKEAREA_NHD": 4.0}).flow_cfs > 0
    with pytest.raises(OutOfRangeError, match="LAKEAREA_NHD"):
        evaluate(eq, {**chars, "LAKEAREA_NHD": 32.1})
    assert "applicability" not in tool.equation_json("W", g, names, 0.01, res, lake_limit=None)


def test_a_lake_area_variable_takes_the_limit_as_its_maximum(tool, monkeypatch):
    lake = ("LAKEAREA", "LAKEAREA_NHD", "log10", 1.0, "percent")
    monkeypatch.setitem(tool.CANDIDATES, "log10 (LAKEAREA+1)", lake)
    monkeypatch.setitem(tool.OFFSETS, "log10 (LAKEAREA+1)", 1.0)
    g = _with_lakes(list(np.linspace(0.0, 8.0, 40)))
    names = ["log10 DRNAREA", "log10 (LAKEAREA+1)"]
    d = tool.equation_json("W", g, names, 0.01, tool.fit(g, names, 0.01), lake_limit=5.0)
    assert "applicability" not in d
    assert d["variables"][1]["maximum"] == 5.0 and d["variables"][1]["offset"] == 1.0
