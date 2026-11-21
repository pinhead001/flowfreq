"""Tests for flowfreq.ice: event counting, cold air and the air index."""

import numpy as np
import pandas as pd
import pytest

from flowfreq.ice import (
    cold_air_days,
    cold_spells_per_winter,
    count_events,
    decade_scaled,
    events_in_months,
    gauge_ice_qaqc,
    ice_days_from_qualifiers,
    winter_event_counts,
    winter_of,
)


def test_winter_of_is_water_year():
    d = pd.Series(pd.to_datetime(["2014-09-30", "2014-10-01", "2014-12-31", "2015-03-31"]))
    out = winter_of(d)
    assert out.tolist() == [2014, 2015, 2015, 2015]
    assert str(out.dtype) == "int64"


def test_count_events_merges_gaps_up_to_merge_gap_days():
    d = pd.Series(pd.to_datetime(["2015-01-01", "2015-01-02", "2015-01-05", "2015-01-09"]))
    ev = count_events(d, merge_gap_days=2)
    # 01-02 -> 01-05 is a 2-day gap: merged; 01-05 -> 01-09 is 3: split
    assert ev["n_days"].tolist() == [3, 1]
    assert ev["start"].dt.strftime("%m-%d").tolist() == ["01-01", "01-09"]


def test_count_events_empty():
    ev = count_events(pd.Series([], dtype="datetime64[ns]"), 2)
    assert ev.empty
    assert list(ev.columns) == ["start", "end", "n_days"]


def test_count_events_unsorted_and_duplicate_dates():
    d = pd.Series(pd.to_datetime(["2015-01-09", "2015-01-01", "2015-01-01", "2015-01-02"]))
    ev = count_events(d, merge_gap_days=0)
    assert ev["n_days"].tolist() == [2, 1]
    assert ev["end"].dt.strftime("%m-%d").tolist() == ["01-02", "01-09"]


def test_cold_air_days_lookback_and_missing():
    idx = pd.date_range("2015-01-01", periods=5)
    tmin = pd.Series([-12.0, -3.0, -2.0, -1.0, float("nan")], index=idx)
    cold = cold_air_days(tmin, threshold_c=-10.0, lookback_days=2)
    assert cold.tolist() == [True, True, True, False, False]
    assert cold.dtype == bool


def test_cold_spells_per_winter_counts_spells_and_air_days():
    idx = pd.date_range("2014-11-01", "2015-03-31")
    tmin = pd.Series(0.0, index=idx)
    tmin.loc["2014-12-10":"2014-12-12"] = -15.0
    tmin.loc["2015-01-20"] = -11.0
    tmin.loc["2015-02-01":"2015-02-05"] = float("nan")
    out = cold_spells_per_winter(tmin, -10.0, 2, [11, 12, 1, 2, 3])
    assert out.loc[2015, "n_spells"] == 2 and out.loc[2015, "n_air_days"] == len(idx) - 5


def test_events_in_months_drops_days_outside_the_winter():  # Review Focus 2
    d = pd.Series(
        pd.to_datetime(["2014-10-30", "2014-10-31", "2014-11-01", "2015-03-31", "2015-04-01"])
    )
    ev = events_in_months(d, [11, 12, 1, 2, 3], merge_gap_days=2)
    assert ev["n_days"].tolist() == [1, 1] and ev["winter"].tolist() == [2015, 2015]
    assert ev["month"].tolist() == [11, 3]


def test_decade_scaled_rate_and_poisson_share():
    n = pd.Series({2013: 2, 2014: 4})  # 6 events
    air = pd.Series({2012: 6, 2013: 3, 2014: 3, 2015: 6})  # 6 spells in the site's winters
    out = decade_scaled(n, air, reference_winters=[2012, 2013, 2014, 2015])
    assert out["rate_per_spell"] == 1.0
    assert out["events_per_winter_ref"] == 4.5  # 1.0 * mean(6, 3, 3, 6)
    assert out["events_ref_decade"] == 45.0
    assert out["share_winters_ref"] == pytest.approx(np.mean([1 - np.exp(-6), 1 - np.exp(-3)] * 2))
    assert (out["events_per_winter_obs"], out["share_winters_obs"], out["n_winters"]) == (
        3.0,
        1.0,
        2,
    )


def test_decade_scaled_no_event_is_zero_not_nan():
    out = decade_scaled(pd.Series({2013: 0}), pd.Series({2013: 5, 2014: 7}), [2013, 2014])
    assert out["rate_per_spell"] == 0.0 and out["share_winters_ref"] == 0.0


def test_decade_scaled_no_spell_in_site_winters_is_nan():  # Review Focus 1
    out = decade_scaled(pd.Series({2013: 2}), pd.Series({2013: 0, 2014: 7}), [2013, 2014])
    for key in (
        "rate_per_spell",
        "events_per_winter_ref",
        "events_ref_decade",
        "share_winters_ref",
    ):
        assert np.isnan(out[key]), key


def test_winter_event_counts_fills_zero_winters():
    ev = pd.DataFrame({"winter": [2014, 2014, 2016]})
    out = winter_event_counts(ev, [2014, 2015, 2016])
    assert out["n_events"].tolist() == [2, 0, 1]


def test_ice_days_need_the_ice_token():
    daily = pd.DataFrame(
        {"qualification_code": ["A", "A:e", "A:e:ICE", "P:ICE"]},
        index=pd.date_range("2015-01-01", periods=4),
    )
    out = ice_days_from_qualifiers(daily)
    assert list(out.columns) == ["date", "ice", "estimated"]
    assert out["ice"].tolist() == [False, False, True, True]
    assert out["estimated"].tolist() == [False, True, True, False]


def test_ice_days_estimated_word_missing_code_and_tz():
    idx = pd.date_range("2015-01-01", periods=3, tz="UTC")
    daily = pd.DataFrame({"qualification_code": ["A:ESTIMATED", None, np.nan]}, index=idx)
    out = ice_days_from_qualifiers(daily)
    assert out["estimated"].tolist() == [True, False, False]
    assert out["ice"].tolist() == [False, False, False]
    assert out["date"].dt.tz is None
    assert out["date"].dt.strftime("%Y-%m-%d").tolist() == [
        "2015-01-01",
        "2015-01-02",
        "2015-01-03",
    ]


def test_qaqc_flags_both_kinds():
    idx = pd.date_range("2015-01-01", periods=12)
    ice = pd.Series([1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], index=idx).astype(bool)
    cold = pd.Series([0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0], index=idx).astype(bool)
    out = gauge_ice_qaqc(ice, cold, min_cold_spell_days=5, merge_gap_days=0)
    assert list(out.columns) == ["kind", "start", "end", "n_days"]
    assert out["kind"].tolist() == ["ice_without_cold", "cold_without_ice"]
    assert out["n_days"].tolist() == [2, 5]


def test_qaqc_minimum_length_and_empty():
    idx = pd.date_range("2015-01-01", periods=6)
    ice = pd.Series(False, index=idx)
    cold = pd.Series([True, True, False, False, False, False], index=idx)
    out = gauge_ice_qaqc(ice, cold, min_cold_spell_days=3, merge_gap_days=0)
    assert out.empty
    assert list(out.columns) == ["kind", "start", "end", "n_days"]
    assert str(out["n_days"].dtype) == "int64"


def test_qaqc_aligns_on_union_index():
    ice = pd.Series([True], index=pd.to_datetime(["2015-01-01"]))
    cold = pd.Series([True, True], index=pd.to_datetime(["2015-01-02", "2015-01-03"]))
    out = gauge_ice_qaqc(ice, cold, min_cold_spell_days=1, merge_gap_days=0)
    assert out["kind"].tolist() == ["ice_without_cold", "cold_without_ice"]
    assert out["n_days"].tolist() == [1, 2]


def test_cold_air_days_looks_back_by_calendar_day_across_a_gap():
    tmin = pd.Series(
        [-12.0, 0.0], index=pd.DatetimeIndex(["2020-01-01", "2020-01-10"], name="date")
    )
    assert cold_air_days(tmin, -10.0, 2).tolist() == [True, False]


def test_cold_air_days_keeps_callers_order():
    idx = pd.DatetimeIndex(["2020-01-03", "2020-01-01", "2020-01-02"])
    tmin = pd.Series([0.0, -12.0, 0.0], index=idx)
    out = cold_air_days(tmin, -10.0, 1)
    assert out.index.equals(idx)
    assert out.tolist() == [False, True, True]
