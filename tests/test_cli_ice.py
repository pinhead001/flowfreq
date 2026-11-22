"""flowfreq ice-events, offline, with both fetchers replaced by synthetic frames."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
import requests
from click.testing import CliRunner

from flowfreq.cli import cli

START, END = "2020-01-01", "2020-01-20"
ARGS = [
    "ice-events", "--gauge", "01234567", "--air-station", "USC00000001",
    "--start", START, "--end", END, "--months", "11,12,1,2,3", "--air-c", "-10",
    "--lookback-days", "0", "--merge-gap-days", "1", "--min-cold-spell-days", "3",
]  # fmt: skip


def _air() -> pd.DataFrame:
    idx = pd.date_range(START, END, freq="D", name="date")
    tmin = pd.Series(0.0, index=idx)
    tmin["2020-01-05":"2020-01-08"] = -12.0  # cold 5-8
    tmin["2020-01-15":"2020-01-18"] = -12.0  # cold 15-18, gauge has no ice
    return pd.DataFrame({"tmin_c": tmin, "tmax_c": tmin + 5.0})


def _daily(site: str, **kw) -> pd.DataFrame:
    idx = pd.date_range(START, END, freq="D", name="date")
    code = pd.Series("A", index=idx, dtype="object")
    code["2020-01-06":"2020-01-10"] = "A:ICE"  # ice 6-10; cold 5-8 -> event 6-8
    flow = pd.Series(np.where(code == "A", 10.0, np.nan), index=idx)
    return pd.DataFrame({"flow_cfs": flow, "qualification_code": code})


@pytest.fixture
def patched(monkeypatch):
    monkeypatch.setattr("flowfreq.ghcn.fetch_ghcn_daily", lambda *a, **k: _air())
    monkeypatch.setattr("flowfreq.waterdata.download_daily", lambda site, **k: _daily(site))


def test_ice_events_cli_writes_outputs(patched, tmp_path):
    out = tmp_path / "new" / "dir"
    res = CliRunner().invoke(cli, ARGS + ["--out", str(out)])
    assert res.exit_code == 0, res.output
    for name in (
        "air_USC00000001.csv",
        "daily_01234567.csv",
        "gauge_ice_events.csv",
        "gauge_ice_qaqc.csv",
        "ice_events_run.json",
    ):
        assert (out / name).exists(), name
    ev = pd.read_csv(out / "gauge_ice_events.csv", dtype={"site": str})
    assert list(ev.columns) == ["site", "winter", "month", "start", "end", "n_days"]
    assert len(ev) == 1
    row = ev.iloc[0]
    assert (row["site"], row["winter"], row["month"]) == ("01234567", 2020, 1)
    assert (row["start"], row["end"], row["n_days"]) == ("2020-01-06", "2020-01-08", 3)
    qa = pd.read_csv(out / "gauge_ice_qaqc.csv", dtype={"site": str})
    assert list(qa.columns) == ["site", "kind", "start", "end", "n_days"]
    assert set(qa["kind"]) == {"ice_without_cold", "cold_without_ice"}
    assert (qa["site"] == "01234567").all()
    daily = pd.read_csv(out / "daily_01234567.csv")
    assert list(daily.columns) == ["date", "flow_cfs", "qualification_code"]
    assert daily["date"].iloc[0] == "2020-01-01"
    run = json.loads((out / "ice_events_run.json").read_text())
    assert run["arguments"]["air_c"] == -10.0
    assert run["gauges"]["01234567"] == {"daily_rows": 20, "ice_days": 5}
    assert "flowfreq_version" in run and run["utc_time"].endswith("+00:00")


def test_empty_results_still_have_headers(monkeypatch, tmp_path):
    monkeypatch.setattr("flowfreq.ghcn.fetch_ghcn_daily", lambda *a, **k: _air())
    monkeypatch.setattr(
        "flowfreq.waterdata.download_daily",
        lambda site, **k: _daily(site).assign(qualification_code="A"),
    )
    res = CliRunner().invoke(cli, ARGS + ["--out", str(tmp_path)])
    assert res.exit_code == 0, res.output
    first = (tmp_path / "gauge_ice_events.csv").read_text().splitlines()
    assert first == ["site,winter,month,start,end,n_days"]


def test_all_nan_air_fails_and_writes_nothing(monkeypatch, tmp_path):
    monkeypatch.setattr("flowfreq.ghcn.fetch_ghcn_daily", lambda *a, **k: _air() * np.nan)
    monkeypatch.setattr("flowfreq.waterdata.download_daily", lambda site, **k: _daily(site))
    out = tmp_path / "o"
    res = CliRunner().invoke(cli, ARGS + ["--out", str(out)])
    assert res.exit_code != 0
    assert "USC00000001" in res.output and "no daily minimum temperature" in res.output
    assert not out.exists()


def test_air_tmin_days_recorded(patched, tmp_path):
    res = CliRunner().invoke(cli, ARGS + ["--out", str(tmp_path)])
    assert res.exit_code == 0, res.output
    run = json.loads((tmp_path / "ice_events_run.json").read_text())
    assert run["air_tmin_days"] == 20


def test_start_after_end_is_refused_naming_the_dates(patched, tmp_path):
    args = list(ARGS)
    args[args.index("--start") + 1] = "2020-02-01"
    out = tmp_path / "o"
    res = CliRunner().invoke(cli, args + ["--out", str(out)])
    assert res.exit_code != 0
    assert "2020-02-01" in res.output and END in res.output
    assert not out.exists()


def test_no_qaqc_before_the_gauge_record_starts(monkeypatch, tmp_path):
    monkeypatch.setattr("flowfreq.ghcn.fetch_ghcn_daily", lambda *a, **k: _air())
    monkeypatch.setattr(
        "flowfreq.waterdata.download_daily",
        lambda site, **k: _daily(site).loc["2020-01-10":],
    )
    res = CliRunner().invoke(cli, ARGS + ["--out", str(tmp_path)])
    assert res.exit_code == 0, res.output
    qa = pd.read_csv(tmp_path / "gauge_ice_qaqc.csv")
    assert len(qa) > 0  # the cold spell 15-18 is still reported
    assert (pd.to_datetime(qa["start"]) >= "2020-01-10").all()


def _feature(day: str, value, approval: str, qualifier) -> dict:
    return {
        "type": "Feature",
        "properties": {
            "time_series_id": "a" * 32,
            "monitoring_location_id": "USGS-01234567",
            "parameter_code": "00060",
            "statistic_id": "00003",
            "time": day,
            "value": value,
            "unit_of_measure": "ft^3/s",
            "approval_status": approval,
            "qualifier": qualifier,
        },
        "geometry": None,
    }


def test_consumer_settings_merge_lookback_and_off_season(monkeypatch, tmp_path):
    from flowfreq.waterdata import parse_daily_features

    days = pd.date_range("2020-10-01", "2020-12-31", freq="D", name="date")
    tmin = pd.Series(0.0, index=days)
    tmin["2020-10-10":"2020-10-16"] = -12.0  # October spell, no ICE, outside --months
    tmin["2020-12-10":"2020-12-11"] = -12.0
    air = pd.DataFrame({"tmin_c": tmin, "tmax_c": tmin + 5.0})
    ice_days = {"2020-12-10": "10.0", "2020-12-11": None, "2020-12-13": None}
    feats = []
    for d in days.strftime("%Y-%m-%d"):
        if d in ice_days:
            feats.append(_feature(d, ice_days[d], "Provisional", ["ICE"]))
        else:
            feats.append(_feature(d, "50.0", "Provisional", None))
    daily = parse_daily_features(feats, keep_null=True)
    monkeypatch.setattr("flowfreq.ghcn.fetch_ghcn_daily", lambda *a, **k: air)
    monkeypatch.setattr("flowfreq.waterdata.download_daily", lambda site, **k: daily)
    args = [
        "ice-events", "--gauge", "01234567", "--air-station", "USC00000001",
        "--start", "2020-10-01", "--end", "2020-12-31", "--months", "11,12,1,2,3",
        "--air-c", "-10", "--lookback-days", "2", "--merge-gap-days", "2",
        "--min-cold-spell-days", "5", "--out", str(tmp_path),
    ]  # fmt: skip
    res = CliRunner().invoke(cli, args)
    assert res.exit_code == 0, res.output
    ev = pd.read_csv(tmp_path / "gauge_ice_events.csv", dtype={"site": str})
    assert len(ev) == 1  # Dec 12 has no ICE flag, merged across the gap
    row = ev.iloc[0]
    assert (row["start"], row["end"], row["n_days"]) == ("2020-12-10", "2020-12-13", 3)
    qa = pd.read_csv(tmp_path / "gauge_ice_qaqc.csv")
    assert not qa["start"].astype(str).str.startswith("2020-10").any()
    assert (qa["kind"] != "ice_without_cold").all()  # Dec 13 counted via the lookback


def test_months_validation_rejects_13(patched, tmp_path):
    args = [a if a != "11,12,1,2,3" else "13" for a in ARGS]
    res = CliRunner().invoke(cli, args + ["--out", str(tmp_path)])
    assert res.exit_code != 0
    assert "months" in res.output


def test_missing_required_option_fails(patched, tmp_path):
    i = ARGS.index("--air-c")
    res = CliRunner().invoke(cli, ARGS[:i] + ARGS[i + 2 :] + ["--out", str(tmp_path)])
    assert res.exit_code != 0
    assert "air-c" in res.output


def test_gauge_failure_names_the_gauge_and_writes_nothing(monkeypatch, tmp_path):
    monkeypatch.setattr("flowfreq.ghcn.fetch_ghcn_daily", lambda *a, **k: _air())

    def boom(site, **k):
        raise requests.RequestException("down")

    monkeypatch.setattr("flowfreq.waterdata.download_daily", boom)
    out = tmp_path / "o"
    res = CliRunner().invoke(cli, ARGS + ["--out", str(out)])
    assert res.exit_code != 0
    assert "01234567" in res.output
    assert not out.exists()
