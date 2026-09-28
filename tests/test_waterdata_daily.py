"""Tests for the Water Data OGC API daily-value backend.

``flowfreq.waterdata.download_daily`` and ``USGSgage.download_daily_flow``'s
``backend`` switch. Offline: HTTP is served from trimmed live captures in
``tests/fixtures/waterdata_ogc/daily_*.json`` (2026-09-27) or synthetic rows in
their shape. The legacy-vs-OGC comparisons -- live, and offline on the paired
``dv_*.rdb`` captures -- are in ``tests/test_daily_backend_parity.py``.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Dict, List, Mapping, Optional, Tuple
from unittest.mock import Mock, patch

import pandas as pd
import pytest
import requests

from flowfreq import waterdata
from flowfreq.usgs import USGSgage
from flowfreq.waterdata import (
    AmbiguousTimeSeriesError,
    daily_time_param,
    download_daily,
    parse_daily_features,
)
from tests.fixtures.nwis_rdb import DV_BASIC
from tests.fixtures.waterdata_ogc_responses import BIG_SANDY, load

DAILY_URL = f"{waterdata.WATERDATA_BASE_URL}/daily/items"
BIG_SANDY_DV = "a61c1b7c739f4d89b8b4c0426586a45a"


def _row(
    time: str,
    value: Optional[str],
    *,
    ts: str = BIG_SANDY_DV,
    approval: str = "Approved",
    qualifier: Optional[List[str]] = None,
    statistic: str = "00003",
    param: str = "00060",
) -> Dict[str, Any]:
    return {
        "type": "Feature",
        "properties": {
            "time_series_id": ts,
            "monitoring_location_id": f"USGS-{BIG_SANDY}",
            "parameter_code": param,
            "statistic_id": statistic,
            "time": time,
            "value": value,
            "unit_of_measure": "ft^3/s",
            "approval_status": approval,
            "qualifier": qualifier,
        },
        "geometry": None,
    }


def _response(payload: Dict[str, Any]) -> Mock:
    response = Mock()
    response.status_code = 200
    response.json.return_value = payload
    response.raise_for_status.return_value = None
    return response


class FakeDailyApi:
    """Serves ``daily/items`` requests; pages by exact next-link URL."""

    def __init__(
        self, first: Dict[str, Any], pages: Optional[Dict[str, Dict[str, Any]]] = None
    ) -> None:
        self.first = first
        self.pages = pages or {}
        self.calls: List[Tuple[str, Optional[Dict[str, Any]]]] = []

    def __call__(
        self, url: str, params: Optional[Mapping[str, Any]] = None, timeout: int = 60
    ) -> Mock:
        self.calls.append((url, dict(params) if params is not None else None))
        if params is None:
            return _response(self.pages[url])
        assert url == DAILY_URL, url
        return _response(self.first)

    @property
    def params(self) -> Dict[str, Any]:
        first = self.calls[0][1]
        assert first is not None
        return first


def _collection(features: List[Dict[str, Any]]) -> Dict[str, Any]:
    return {"type": "FeatureCollection", "features": features, "links": []}


# ----------------------------------------------------------------------------
# parse_daily_features
# ----------------------------------------------------------------------------


class TestParseDailyFeatures:
    def test_frame_shape_and_sorting(self) -> None:
        df = parse_daily_features(
            [
                _row("2024-01-03", "77.9"),
                _row("2024-01-01", "82.5"),
                _row("2024-01-02", "79.6"),
            ]
        )
        assert list(df.columns) == ["flow_cfs", "qualification_code"]
        assert df.index.name == "date"
        assert df.index.tz is None
        assert list(df.index) == list(pd.to_datetime(["2024-01-01", "2024-01-02", "2024-01-03"]))
        assert df["flow_cfs"].tolist() == [82.5, 79.6, 77.9]
        assert df["flow_cfs"].dtype == float

    def test_qualifiers_map_as_on_the_instantaneous_path(self) -> None:
        df = parse_daily_features(
            [
                _row("2022-12-20", "223", qualifier=["ESTIMATED", "ICE"]),
                _row("2024-03-25", "171", qualifier=["EQUIP", "ESTIMATED"]),
                _row("2026-09-17", "2190", approval="Provisional", qualifier=["ESTIMATED"]),
                _row("2026-09-18", "2270", approval="Provisional"),
            ]
        )
        assert df["qualification_code"].tolist() == ["A:e:ICE", "A:EQUIP:e", "P:e", "P"]

    def test_ice_estimate_is_kept_and_ice_null_is_dropped(self) -> None:
        """An ice day with a published estimate is data; one without is not a zero."""
        df = parse_daily_features(
            [
                _row("2026-01-17", "3120", approval="Provisional"),
                _row("2026-01-18", None, approval="Provisional", qualifier=["ICE"]),
                _row("2022-12-20", "223", qualifier=["ESTIMATED", "ICE"]),
            ]
        )
        assert df["flow_cfs"].tolist() == [223.0, 3120.0]
        assert pd.Timestamp("2026-01-18") not in df.index

    def test_empty_gives_typed_empty_frame(self) -> None:
        df = parse_daily_features([])
        assert df.empty and list(df.columns) == ["flow_cfs", "qualification_code"]
        assert df.index.name == "date"

    def test_other_statistic_is_refused(self) -> None:
        with pytest.raises(ValueError, match="statistic '00001'"):
            parse_daily_features([_row("2024-01-01", "1", statistic="00001")])

    def test_other_parameter_is_refused(self) -> None:
        with pytest.raises(ValueError, match="parameter '00065'"):
            parse_daily_features([_row("2024-01-01", "1", param="00065")])

    def test_unsupported_parameter_code_raises(self) -> None:
        with pytest.raises(ValueError, match="Unsupported daily parameter"):
            parse_daily_features([], param_cd="00065")

    def test_two_series_are_refused_not_merged(self) -> None:
        with pytest.raises(AmbiguousTimeSeriesError, match="2 separate") as info:
            parse_daily_features([_row("2024-01-01", "1"), _row("2024-01-01", "2", ts="b" * 32)])
        assert "ts_id" in str(info.value)

    def test_ts_id_accepts_its_series_and_rejects_others(self) -> None:
        dashed = "a61c1b7c-739f-4d89-b8b4-c0426586a45a"
        df = parse_daily_features([_row("2024-01-01", "1")], ts_id=dashed)
        assert len(df) == 1
        with pytest.raises(ValueError, match="refusing to merge"):
            parse_daily_features([_row("2024-01-01", "1", ts="b" * 32)], ts_id=BIG_SANDY_DV)

    def test_identical_duplicate_dates_collapse(self) -> None:
        df = parse_daily_features([_row("2024-01-01", "1"), _row("2024-01-01", "1")])
        assert len(df) == 1

    def test_conflicting_duplicate_dates_raise(self) -> None:
        with pytest.raises(ValueError, match="conflicting"):
            parse_daily_features([_row("2024-01-01", "1"), _row("2024-01-01", "2")])


class TestDailyTimeParam:
    @pytest.mark.parametrize(
        ("start", "end", "expected"),
        [
            ("2024-01-01", "2024-12-31", "2024-01-01/2024-12-31"),
            ("2024-01-01", None, "2024-01-01/.."),
            (None, "2024-12-31", "../2024-12-31"),
            (None, None, None),
            ("2024-1-5", "2024-01-05", "2024-01-05/2024-01-05"),
        ],
    )
    def test_forms(self, start: Optional[str], end: Optional[str], expected: Optional[str]) -> None:
        assert daily_time_param(start, end) == expected

    def test_reversed_raises(self) -> None:
        with pytest.raises(ValueError, match="after end_date"):
            daily_time_param("2025-01-01", "2024-01-01")

    def test_malformed_raises(self) -> None:
        with pytest.raises(ValueError, match="Could not parse"):
            daily_time_param("not-a-date", None)


# ----------------------------------------------------------------------------
# download_daily
# ----------------------------------------------------------------------------


class TestDownloadDaily:
    def test_request_parameters(self) -> None:
        api = FakeDailyApi(load("daily_03606500_00060_equip"))
        with patch("flowfreq.waterdata.requests.get", side_effect=api):
            download_daily(BIG_SANDY, start_date="2024-03-20", end_date="2024-03-31")
        params = api.params
        assert params["monitoring_location_id"] == "USGS-03606500"
        assert params["parameter_code"] == "00060"
        assert params["statistic_id"] == "00003"
        assert params["time"] == "2024-03-20/2024-03-31"
        assert params["limit"] == waterdata.PAGE_LIMIT
        assert params["skipGeometry"] == "true"
        assert "time_series_id" not in params

    def test_no_range_omits_time_for_the_full_record(self) -> None:
        api = FakeDailyApi(load("daily_03606500_00060_equip"))
        with patch("flowfreq.waterdata.requests.get", side_effect=api):
            download_daily(BIG_SANDY)
        assert "time" not in api.params

    def test_ts_id_is_sent(self) -> None:
        api = FakeDailyApi(load("daily_03606500_00060_equip"))
        with patch("flowfreq.waterdata.requests.get", side_effect=api):
            download_daily(BIG_SANDY, ts_id=BIG_SANDY_DV.upper())
        assert api.params["time_series_id"] == BIG_SANDY_DV

    def test_captured_window(self) -> None:
        api = FakeDailyApi(load("daily_03606500_00060_equip"))
        with patch("flowfreq.waterdata.requests.get", side_effect=api):
            df = download_daily(BIG_SANDY, start_date="2024-03-20", end_date="2024-03-31")
        assert len(df) == 12
        assert df.index.is_monotonic_increasing and df.index.is_unique
        assert df.loc["2024-03-25", "qualification_code"] == "A:EQUIP:e"
        assert df.loc["2024-03-20", "qualification_code"] == "A"

    def test_follows_next_links(self) -> None:
        pages = [load(f"daily_03606500_00060_page{n}") for n in (1, 2, 3)]
        links = {}
        for this, nxt in zip(pages, pages[1:]):
            href = next(l["href"] for l in this["links"] if l["rel"] == "next")
            links[href] = nxt
        api = FakeDailyApi(pages[0], links)
        with patch("flowfreq.waterdata.requests.get", side_effect=api):
            paged = download_daily(BIG_SANDY, start_date="2024-03-20", end_date="2024-03-31")
        whole = parse_daily_features(load("daily_03606500_00060_equip")["features"])
        assert len(api.calls) == 3
        pd.testing.assert_frame_equal(paged, whole)

    def test_no_rows_raises_like_legacy(self) -> None:
        api = FakeDailyApi(load("daily_03606500_00060_empty"))
        with patch("flowfreq.waterdata.requests.get", side_effect=api):
            with pytest.raises(ValueError, match="No daily data found for site 03606500"):
                download_daily(BIG_SANDY, start_date="1900-01-01", end_date="1900-01-05")

    def test_request_failure_propagates(self) -> None:
        with patch(
            "flowfreq.waterdata.requests.get", side_effect=requests.ConnectionError("reset")
        ):
            with pytest.raises(requests.RequestException, match="daily values for site"):
                download_daily(BIG_SANDY, start_date="2024-01-01", end_date="2024-01-02")

    def test_reversed_range_sends_nothing(self) -> None:
        with patch("flowfreq.waterdata.requests.get") as get:
            with pytest.raises(ValueError, match="after end_date"):
                download_daily(BIG_SANDY, start_date="2025-01-01", end_date="2024-01-01")
        get.assert_not_called()


# ----------------------------------------------------------------------------
# USGSgage.download_daily_flow backend switch
# ----------------------------------------------------------------------------


class TestDownloadDailyFlowBackend:
    def test_default_is_waterdata(self) -> None:
        api = FakeDailyApi(load("daily_03606500_00060_equip"))
        gage = USGSgage(BIG_SANDY)
        with patch("flowfreq.waterdata.requests.get", side_effect=api):
            with patch.object(
                USGSgage,
                "_download_daily_flow_legacy",
                side_effect=AssertionError("legacy NWIS must not run by default"),
            ):
                df = gage.download_daily_flow("2024-03-20", "2024-03-31")
        assert list(df.columns) == ["flow_cfs"]
        assert df.index.name == "date"
        assert len(df) == 12
        assert gage.daily_data is df

    def test_default_range_matches_legacy(self) -> None:
        """Same start (DEFAULT_START_DATE) and UTC end date as the legacy path."""
        api = FakeDailyApi(load("daily_03606500_00060_equip"))
        fixed_now = datetime(2026, 1, 1, 0, 30, tzinfo=timezone.utc)
        with (
            patch("flowfreq.waterdata.requests.get", side_effect=api),
            patch("flowfreq.usgs.datetime") as mock_datetime,
        ):
            mock_datetime.now.return_value = fixed_now
            USGSgage(BIG_SANDY).download_daily_flow()
        assert api.params["time"] == f"{USGSgage.DEFAULT_START_DATE}/2026-01-01"
        mock_datetime.now.assert_called_with(timezone.utc)

    def test_timeout_is_passed(self) -> None:
        api = FakeDailyApi(load("daily_03606500_00060_equip"))
        with patch("flowfreq.waterdata.requests.get", side_effect=api) as get:
            USGSgage(BIG_SANDY).download_daily_flow(timeout=120)
        assert get.call_args.kwargs["timeout"] == 120

    def test_reversed_range_raises_before_any_request(self) -> None:
        with patch("flowfreq.waterdata.requests.get") as get:
            with pytest.raises(ValueError, match="after end_date"):
                USGSgage(BIG_SANDY).download_daily_flow("2025-01-01", "2020-01-01")
        get.assert_not_called()

    def test_legacy_backend_still_reads_the_rdb_service(self) -> None:
        response = Mock(text=DV_BASIC, status_code=200)
        response.raise_for_status.return_value = None
        with patch("flowfreq.usgs.requests.get", return_value=response) as get:
            df = USGSgage("12449500").download_daily_flow(backend="nwis-legacy")
        assert get.call_args.args[0] == USGSgage.BASE_URL_DAILY
        assert list(df.columns) == ["flow_cfs"] and len(df) == 5

    def test_unknown_backend_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown daily-value backend 'nwis'"):
            USGSgage(BIG_SANDY).download_daily_flow(backend="nwis")
