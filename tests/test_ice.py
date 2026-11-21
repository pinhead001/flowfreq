"""Tests for flowfreq.ice: event counting, cold air and the air index."""

import pandas as pd

from flowfreq.ice import (
    cold_air_days,
    cold_spells_per_winter,
    count_events,
    events_in_months,
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
