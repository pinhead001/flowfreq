"""Tests for flowfreq.peak_sources."""

from __future__ import annotations

import pandas as pd
import pytest

from flowfreq.peak_sources import (
    DEFAULT_BACKEND,
    PEAK_COLUMNS,
    LegacyNwisBackend,
    PeakDataBackend,
    WaterDataApiBackend,
    get_backend,
    validate_peak_frame,
)


def _frame(**overrides):
    base = {
        "water_year": [2000, 2001],
        "peak_date": pd.to_datetime(["2000-03-01", "2001-04-02"]),
        "peak_flow_cfs": [100.0, 200.0],
        "qualification_code": ["", "2"],
        "extra": [1, 2],
    }
    base.update(overrides)
    return pd.DataFrame(base)


def test_validate_restricts_and_orders_columns():
    assert tuple(validate_peak_frame(_frame()).columns) == PEAK_COLUMNS


@pytest.mark.parametrize(
    "overrides, match",
    [
        ({"water_year": [2000, 2000]}, "Duplicate"),
        ({"peak_flow_cfs": [-1.0, 2.0]}, "Negative"),
    ],
)
def test_validate_rejects(overrides, match):
    with pytest.raises(ValueError, match=match):
        validate_peak_frame(_frame(**overrides))


def test_validate_missing_column():
    with pytest.raises(ValueError, match="missing"):
        validate_peak_frame(_frame().drop(columns="qualification_code"))


def test_default_is_legacy_and_protocol_holds():
    b = get_backend()
    assert b.name == DEFAULT_BACKEND == LegacyNwisBackend.name
    assert isinstance(b, PeakDataBackend)


def test_legacy_backend_wraps_usgsgage(monkeypatch):
    monkeypatch.setattr("flowfreq.usgs.USGSgage.download_peak_flow", lambda self: _frame())
    out = LegacyNwisBackend().fetch_peaks("12345678")
    assert list(out["peak_flow_cfs"]) == [100.0, 200.0]


def test_waterdata_backend_not_implemented_until_verified():
    with pytest.raises(NotImplementedError, match="#29"):
        WaterDataApiBackend().fetch_peaks("12345678")


def test_unknown_backend():
    with pytest.raises(KeyError):
        get_backend("nope")
