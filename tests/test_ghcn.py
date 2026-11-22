"""NCEI daily air temperature (flowfreq.ghcn), on recorded responses."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from flowfreq import ghcn

CSV = (
    '"STATION","DATE","TMAX","TMIN"\n'
    '"USC00000001","2020-01-01","1.1","-12.2"\n'
    '"USC00000001","2020-01-02","0.0",""\n'
    '"USC00000001","2020-01-04","2.5","-3.0"\n'
    '"USC00000001","2020-01-05","","-8.5"\n'
)


def _serve(monkeypatch, text, seen=None):
    def fake(url, params, timeout, **kw):
        if seen is not None:
            seen.update(url=url, params=params, kw=kw)
        return SimpleNamespace(text=text)

    monkeypatch.setattr("flowfreq.waterdata.request", fake)


def test_fetch_ghcn_daily_parses_recorded_csv(monkeypatch):
    seen: dict = {}
    _serve(monkeypatch, CSV, seen)
    df = ghcn.fetch_ghcn_daily("USC00000001", "2020-01-01", "2020-01-06")
    assert seen["url"] == ghcn.NCEI_URL
    assert seen["params"]["dataTypes"] == "TMIN,TMAX"
    assert seen["params"]["units"] == "metric"
    assert seen["kw"] == {"send_api_key": False}
    assert df.index.name == "date"
    assert len(df) == 6
    assert list(df.columns) == ["tmin_c", "tmax_c"]
    assert df.loc["2020-01-01", "tmin_c"] == -12.2
    assert df.loc["2020-01-01", "tmax_c"] == 1.1
    assert np.isnan(df.loc["2020-01-02", "tmin_c"])  # empty string
    assert df.loc["2020-01-02", "tmax_c"] == 0.0
    assert np.isnan(df.loc["2020-01-03", "tmin_c"])  # day not reported
    assert np.isnan(df.loc["2020-01-05", "tmax_c"])
    assert df.loc["2020-01-05", "tmin_c"] == -8.5
    assert df.dtypes.map(str).tolist() == ["float64", "float64"]


def test_missing_tmax_column_is_nan(monkeypatch):
    _serve(monkeypatch, '"STATION","DATE","TMIN"\n"X","2020-01-01","-1.0"\n')
    df = ghcn.fetch_ghcn_daily("X", "2020-01-01", "2020-01-02")
    assert df["tmin_c"].iloc[0] == -1.0
    assert df["tmax_c"].isna().all()


@pytest.mark.parametrize("text", ["", '"STATION","DATE","TMAX","TMIN"\n'])
def test_empty_response_is_all_nan_with_warning(monkeypatch, caplog, text):
    _serve(monkeypatch, text)
    with caplog.at_level("WARNING"):
        df = ghcn.fetch_ghcn_daily("X", "2020-01-01", "2020-01-03")
    assert len(df) == 3 and df.isna().all().all()
    assert "no daily rows" in caplog.text


def test_start_after_end_raises():
    with pytest.raises(ValueError):
        ghcn.fetch_ghcn_daily("X", "2020-02-01", "2020-01-01")
