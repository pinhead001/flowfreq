"""Lombard change-point tests, bootstrap field significance and basin change events
(roadmap section 6.1; methow item 8)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from flowfreq.regulation import NIDInventory, nid_dams_in_basin
from flowfreq.stationarity import (
    LOMBARD_ABRUPT_CRITICAL,
    _lombard_scores,
    _lombard_stat,
    _mk_p,
    attribute_change_point,
    basin_change_events,
    field_significance_bootstrap,
    lombard,
    mann_kendall,
    stationarity_screen,
)

# ------------------------------------------------------------------ Lombard --


@pytest.mark.parametrize("score", ["wilcoxon", "mood"])
def test_lombard_null_quantiles_match_the_asymptotic_values(score: str) -> None:
    rng = np.random.default_rng(11)
    z = np.vstack([_lombard_scores(x, score) for x in rng.normal(size=(3000, 150))])
    abrupt = np.quantile(_lombard_stat(z, "abrupt"), [0.90, 0.95, 0.99])
    for got, level in zip(abrupt, (0.10, 0.05, 0.01)):
        assert got == pytest.approx(LOMBARD_ABRUPT_CRITICAL[level], rel=0.08)
    assert np.quantile(_lombard_stat(z, "smooth"), 0.95) == pytest.approx(0.040, rel=0.08)


def test_lombard_smooth_statistic_by_brute_force() -> None:
    z = _lombard_scores(np.random.default_rng(2).normal(size=12), "wilcoxon")
    n = z.size
    partial = np.cumsum(z)
    total = sum(np.sum(partial[t1:t2]) ** 2 for t1 in range(1, n) for t2 in range(t1 + 1, n + 1))
    assert float(_lombard_stat(z, "smooth")) == pytest.approx(total / n**5)
    assert float(_lombard_stat(z, "abrupt")) == pytest.approx(np.sum(partial**2) / n**2)


def test_lombard_finds_a_step_and_its_position() -> None:
    rng = np.random.default_rng(3)
    x = np.concatenate([rng.normal(0, 1, 40), rng.normal(2, 1, 40)])
    years = list(range(1940, 2020))
    r = lombard(x, years, seed=1)
    assert r.p_value < 0.01
    assert abs(r.start_index - 40) <= 3
    assert r.start_year == years[r.start_index]
    assert r.end_index == r.start_index


def test_lombard_smooth_brackets_a_ramp() -> None:
    rng = np.random.default_rng(5)
    level = np.concatenate([np.zeros(30), np.linspace(0, 3, 30), np.full(30, 3.0)])
    r = lombard(level + rng.normal(0, 1, 90), kind="smooth", seed=2)
    assert r.p_value < 0.01
    # The level leaves 0 at index 30 and reaches 3 at index 59.
    assert 22 <= r.start_index <= 38
    assert 51 <= r.end_index <= 67


def test_lombard_mood_sees_a_change_in_spread_that_wilcoxon_does_not() -> None:
    rng = np.random.default_rng(8)
    x = np.concatenate([rng.normal(0, 0.3, 50), rng.normal(0, 2.0, 50)])
    assert lombard(x, score="mood", seed=0).p_value < 0.01
    assert lombard(x, score="wilcoxon", seed=0).p_value > 0.05


def test_lombard_is_rank_based_and_reproducible() -> None:
    x = np.random.default_rng(9).normal(size=40)
    a = lombard(x, seed=4)
    b = lombard(np.exp(x), seed=4)
    assert a.statistic == pytest.approx(b.statistic)
    assert a.p_value == b.p_value
    assert 0 < a.p_value <= 1


@pytest.mark.parametrize(
    "kwargs",
    [{"kind": "linear"}, {"score": "ansari"}, {"n_permutations": 0}],
)
def test_lombard_refuses_bad_options(kwargs) -> None:
    with pytest.raises(ValueError):
        lombard(np.arange(10.0), **kwargs)


def test_lombard_refuses_constant_short_or_misaligned_series() -> None:
    with pytest.raises(ValueError, match="constant"):
        lombard(np.ones(10))
    with pytest.raises(ValueError, match="at least 5"):
        lombard([1.0, 2.0, 3.0])
    with pytest.raises(ValueError, match="differ in length"):
        lombard(np.arange(10.0), years=range(9))


def test_screen_reports_lombard_and_flags_a_change_in_variability() -> None:
    rng = np.random.default_rng(12)
    logq = np.concatenate([rng.normal(3, 0.05, 40), rng.normal(3, 0.5, 40)])
    values = pd.Series(10**logq, index=range(1940, 2020))
    out = stationarity_screen(values)
    assert set(out["lombard"]) == {
        "wilcoxon_abrupt",
        "wilcoxon_smooth",
        "mood_abrupt",
        "mood_smooth",
    }
    assert out["lombard"]["mood_abrupt"]["p_value"] < 0.01
    assert any("change in variability" in f for f in out["flags"])
    # Reproducible by default.
    assert stationarity_screen(values)["lombard"] == out["lombard"]


# ------------------------------------------- bootstrap field significance --


def test_mk_p_matches_mann_kendall_with_ties() -> None:
    rng = np.random.default_rng(1)
    x = np.round(rng.normal(size=40), 1)
    assert _mk_p(x) == pytest.approx(mann_kendall(x).p_value)
    assert _mk_p(np.ones(10)) == 1.0


def _frame(x: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame(x, index=range(1950, 1950 + x.shape[0]))


def test_bootstrap_finds_a_regional_trend() -> None:
    rng = np.random.default_rng(2)
    t = np.arange(50)[:, None]
    out = field_significance_bootstrap(
        _frame(0.06 * t + rng.normal(size=(50, 15))), n_bootstrap=200, seed=0
    )
    assert out["n_sites"] == 15
    assert out["n_significant"] >= 10
    assert out["field_significant"]
    assert out["p_value"] < 0.01


def test_bootstrap_widens_the_null_for_correlated_sites() -> None:
    rng = np.random.default_rng(3)
    independent = rng.normal(size=(60, 20))
    common = rng.normal(size=(60, 1))
    correlated = 0.95 * common + np.sqrt(1 - 0.95**2) * rng.normal(size=(60, 20))
    a = field_significance_bootstrap(_frame(independent), n_bootstrap=300, seed=1)
    b = field_significance_bootstrap(_frame(correlated), n_bootstrap=300, seed=1)
    # Independent sites: about the binomial 95% point (3 of 20). Correlated: far more.
    assert a["critical_count"] <= 4
    assert b["critical_count"] >= 2 * a["critical_count"]


def test_bootstrap_blocks_missing_values_and_checks() -> None:
    rng = np.random.default_rng(4)
    df = _frame(rng.normal(size=(30, 3)))
    df.iloc[:25, 2] = np.nan  # five values: left out
    out = field_significance_bootstrap(df, n_bootstrap=50, block_length=5, seed=0)
    assert out["n_sites"] == 2 and out["block_length"] == 5
    assert set(out["site_p_values"]) == {"0", "1"}
    with pytest.raises(ValueError, match="block_length"):
        field_significance_bootstrap(df, block_length=31)
    with pytest.raises(ValueError, match="n_bootstrap"):
        field_significance_bootstrap(df, n_bootstrap=0)
    with pytest.raises(ValueError, match="no site"):
        field_significance_bootstrap(df.iloc[:8])


# ------------------------------------------------- basin change events --

SQUARE = {
    "type": "Polygon",
    "coordinates": [
        [[-120.0, 47.0], [-119.0, 47.0], [-119.0, 48.0], [-120.0, 48.0], [-120.0, 47.0]]
    ],
}


def _nid() -> NIDInventory:
    dams = pd.DataFrame(
        {
            "nid_id": ["WA001", "WA002", "WA003", "WA004", "OUT01"],
            "latitude": [47.5, 47.6, 47.2, 47.9, 49.5],
            "longitude": [-119.5, -119.4, -119.8, -119.1, -119.5],
            "nid_storage_acft": [60000.0, 30000.0, 9000.0, 1000.0, 1e6],
            "year_completed": [1961.0, 1938.0, np.nan, 1990.0, 1950.0],
        }
    )
    return NIDInventory(dams, "2026-10-01", "2026-10-15")


def test_nid_dams_in_basin_largest_first() -> None:
    hit = nid_dams_in_basin(SQUARE, _nid())
    assert list(hit["nid_id"]) == ["WA001", "WA002", "WA003", "WA004"]


def test_basin_change_events_dams_and_half_storage() -> None:
    out = basin_change_events(geometry=SQUARE, nid=_nid())
    assert out["n_dams"] == 4
    labels = [label for label, _ in out["events"]]
    assert ("dam WA001 completed (60,000 acre-ft, 60% of basin NID storage)", 1961) in out["events"]
    assert any("WA002" in label for label in labels)
    # 1938: 30,000 of 100,000; 1961: 90,000, past half.
    assert ("half of today's basin NID storage in place", 1961) in out["events"]
    assert not any("WA004" in label for label, _ in out["events"])  # 1% share
    assert out["notes"] == ["1 of 4 NID dams have no completion year"]


def test_basin_change_events_impervious_crossing_and_notes() -> None:
    out = basin_change_events(impervious={2001: 2.0, 2011: 4.1, 2019: 6.3})
    assert out["events"] == [("basin impervious cover reached 5% (NLCD 2011 to 2019)", 2019)]
    assert out["n_dams"] is None
    early = basin_change_events(impervious={2001: 12.0, 2021: 15.0})
    assert early["events"] == [] and "already 12.0%" in early["notes"][0]
    co = basin_change_events(impervious={2001: 1.0, 2011: 9.0}, region="co")
    assert co["events"] == [] and "CO" in co["notes"][0]
    with pytest.raises(ValueError, match="both geometry and nid"):
        basin_change_events(geometry=SQUARE)


def test_basin_events_feed_attribution() -> None:
    events = basin_change_events(geometry=SQUARE, nid=_nid())["events"]
    out = attribute_change_point(1960, [1959, 1960], ["", ""], events=events)
    assert out["explained"]
    assert out["candidates"][0]["year"] == 1961
