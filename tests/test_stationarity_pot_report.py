"""Peaks-over-threshold count trends and the one-page stationarity report (roadmap
section 6.1; methow item 8)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from flowfreq.stationarity import PotCounts, pot_count_trend, pot_counts, stationarity_report


def _daily(years: range, spikes: dict) -> pd.Series:
    """Flow of 10 every day of the water years, with ``spikes`` (date -> value)."""
    idx = pd.date_range(f"{years[0] - 1}-10-01", f"{years[-1]}-09-30", freq="D")
    s = pd.Series(10.0, index=idx)
    for day, value in spikes.items():
        s[pd.Timestamp(day)] = value
    return s


def test_pot_counts_declusters_and_counts_by_water_year() -> None:
    spikes = {
        "2001-03-01": 100.0,
        "2001-03-04": 150.0,  # 3 days after: same event, peak 150
        "2001-03-20": 120.0,  # 16 days after: a new event
        "2002-10-05": 90.0,  # water year 2003
    }
    out = pot_counts(_daily(range(2001, 2006), spikes), threshold=50.0)
    assert isinstance(out, PotCounts)
    assert out.counts.to_dict() == {2001: 2, 2002: 0, 2003: 1, 2004: 0, 2005: 0}
    assert list(out.events["peak"]) == [150.0, 120.0, 90.0]
    assert list(out.events["water_year"]) == [2001, 2001, 2003]
    # A longer separation merges the March events.
    merged = pot_counts(_daily(range(2001, 2006), spikes), threshold=50.0, min_separation_days=20)
    assert merged.counts[2001] == 1


def test_pot_counts_drops_incomplete_years_and_picks_a_threshold() -> None:
    rng = np.random.default_rng(0)
    idx = pd.date_range("1979-10-01", "2009-09-30", freq="D")
    s = pd.Series(np.exp(rng.normal(3, 1, idx.size)), index=idx)
    s = s.drop(s.index[(s.index >= "1995-01-01") & (s.index < "1995-07-01")])  # WY1995 short
    out = pot_counts(s, events_per_year=2.0)
    assert 1995 not in out.counts.index
    assert len(out.counts) == 29
    assert out.counts.sum() == pytest.approx(2.0 * 29, abs=3)


def test_pot_counts_refuses_bad_input() -> None:
    s = _daily(range(2001, 2006), {})
    with pytest.raises(ValueError, match="indexed by date"):
        pot_counts(s.reset_index(drop=True))
    with pytest.raises(ValueError, match="5 complete"):
        pot_counts(s["2002-10-01":])
    with pytest.raises(ValueError, match="min_separation_days"):
        pot_counts(s, min_separation_days=0)
    with pytest.raises(ValueError, match="more than one value"):
        pot_counts(pd.concat([s, s.iloc[:1]]))


def test_pot_count_trend_recovers_a_poisson_trend() -> None:
    rng = np.random.default_rng(1)
    years = np.arange(1900, 2020)
    b = 0.01
    counts = pd.Series(rng.poisson(2.0 * np.exp(b * (years - years.mean()))), index=years)
    out = pot_count_trend(counts)
    assert out["trend"] == "increasing"
    # One record is one draw (standard error about 0.002 here); the average over many
    # is the estimator's bias.
    slopes = [
        pot_count_trend(
            pd.Series(
                np.random.default_rng(seed).poisson(2.0 * np.exp(b * (years - years.mean()))),
                index=years,
            )
        )["slope_log_per_year"]
        for seed in range(60)
    ]
    assert np.mean(slopes) == pytest.approx(b, abs=0.0008)
    assert out["change_per_decade_pct"] == pytest.approx(
        100 * (np.exp(10 * out["slope_log_per_year"]) - 1)
    )
    assert out["total_events"] == int(counts.sum())
    flat = pd.Series(rng.poisson(2.0, 120), index=years)
    assert pot_count_trend(flat)["p_value"] > 0.01


def test_overdispersion_weakens_the_test() -> None:
    rng = np.random.default_rng(2)
    years = np.arange(1950, 2020)
    # Negative binomial counts with a mild trend: dispersion well above 1.
    mean = 2.0 * np.exp(0.008 * (years - years.mean()))
    counts = pd.Series(rng.negative_binomial(1, 1 / (1 + mean)), index=years)
    out = pot_count_trend(counts)
    assert out["dispersion"] > 1.5


def test_pot_count_trend_refuses_bad_counts() -> None:
    years = range(2000, 2010)
    with pytest.raises(ValueError, match="no events"):
        pot_count_trend(pd.Series(0, index=years))
    with pytest.raises(ValueError, match="non-negative integers"):
        pot_count_trend(pd.Series([1.5] * 10, index=years))
    with pytest.raises(ValueError, match="5 years"):
        pot_count_trend(pd.Series([1, 2, 3], index=range(3)))


def test_report_page_and_figure(tmp_path: Path) -> None:
    rng = np.random.default_rng(3)
    years = range(1940, 2020)
    logq = np.concatenate([rng.normal(3.0, 0.2, 40), rng.normal(3.4, 0.2, 40)])
    peaks = pd.Series(10**logq, index=years)
    idx = pd.date_range("1939-10-01", "2019-09-30", freq="D")
    daily = pd.Series(np.exp(rng.normal(3, 1, idx.size)), index=idx)
    pot = pot_counts(daily)
    rep = stationarity_report(
        peaks,
        title="99999999 Test Creek",
        pot=pot,
        events=[("dam completed", 1979)],
    )
    assert len(rep.figure.axes) == 3
    assert rep.markdown.startswith("# 99999999 Test Creek")
    assert "| Pettitt |" in rep.markdown and "Lombard Mood, smooth" in rep.markdown
    assert "## Peaks over a threshold" in rep.markdown
    assert "## Change-point attribution" in rep.markdown and "dam completed" in rep.markdown
    assert rep.screen["pettitt"]["p_value"] < 0.01
    assert rep.pot is not None and rep.pot["n_years"] == len(pot.counts)
    out = tmp_path / "report.png"
    rep.save(out)
    assert out.stat().st_size > 10_000


def test_report_without_counts_has_two_panels() -> None:
    values = pd.Series(np.random.default_rng(4).lognormal(8, 0.5, 30), index=range(1990, 2020))
    rep = stationarity_report(values)
    assert len(rep.figure.axes) == 2
    assert rep.pot is None
    assert "Peaks over a threshold" not in rep.markdown


def test_snowmelt_record_keeps_a_high_threshold_and_warns_when_it_cannot(caplog) -> None:
    # One smooth freshet a year (April-July, 1,000 to 1,500) over a base flow of 10-60: over
    # any threshold between them it is one event, so two a year are found only in the
    # base flow. One a year keeps the threshold in the freshet.
    rng = np.random.default_rng(5)
    idx = pd.date_range("1979-10-01", "2009-09-30", freq="D")
    s = pd.Series(rng.uniform(10, 60, idx.size), index=idx)
    for wy in range(1980, 2010):
        season = pd.date_range(f"{wy}-04-01", f"{wy}-07-31", freq="D")
        s[season] = 1000 + 500 * np.sin(np.pi * np.arange(1, season.size + 1) / (season.size + 1))
    one = pot_counts(s, events_per_year=1.0)
    assert one.threshold > 900 and one.quantile > 0.75
    assert (one.counts == 1).all()
    with caplog.at_level("WARNING", logger="flowfreq.stationarity"):
        two = pot_counts(s, events_per_year=2.0)
    assert two.quantile < 0.75
    assert any("snowmelt" in r.getMessage() for r in caplog.records)
    rep = stationarity_report(
        pd.Series(rng.lognormal(7, 0.3, 30), index=range(1980, 2010)), pot=two
    )
    assert "a low threshold for floods" in rep.markdown
