"""Tests for flowfreq.waterdata, the Water Data OGC API instantaneous-value backend.

Offline: every HTTP call is served from trimmed live captures in
``tests/fixtures/waterdata_ogc/`` (captured 2026-09-25) or small synthetic
payloads built in the shape of those captures. Live checks are marked
``requires_network``.

There is no live legacy-vs-OGC comparison: the legacy NWIS service is not
reachable from where these were written. Parity is checked against the frame
*contract* instead -- the columns, index and dtypes the legacy parser produces
from the committed RDB fixtures.
"""

from __future__ import annotations

import logging
from datetime import timedelta, timezone
from typing import Any, Callable, Dict, List, Mapping, Optional, Tuple
from unittest.mock import Mock, patch

import pandas as pd
import pytest
import requests

from flowfreq import waterdata
from flowfreq.usgs import NoInstantaneousDataError, USGSgage, _parse_iv_rdb
from flowfreq.waterdata import (
    MAX_TIME_ENVELOPE_DAYS,
    AmbiguousTimeSeriesError,
    download_instantaneous,
    fetch_continuous_chunk,
    local_day_bounds,
    map_qualification_code,
    parse_continuous_features,
    plan_chunks,
    resolve_local_zone,
)
from tests.fixtures.nwis_rdb import IV_BASIC, IV_STAGE_BASIC
from tests.fixtures.waterdata_ogc_responses import (
    BIG_SANDY,
    BIG_SANDY_IV_00060,
    OLMSTED,
    OLMSTED_HEADWATER,
    OLMSTED_TAILWATER,
    load,
)

BASE = waterdata.WATERDATA_BASE_URL
ContinuousHandler = Callable[[Mapping[str, Any]], Dict[str, Any]]


def _response(payload: Dict[str, Any]) -> Mock:
    response = Mock()
    response.status_code = 200
    response.json.return_value = payload
    response.raise_for_status.return_value = None
    return response


class FakeApi:
    """Routes ``requests.get`` calls to fixture payloads, recording each call."""

    def __init__(
        self,
        site: str,
        param_cd: str,
        continuous: ContinuousHandler,
        *,
        location: Optional[Dict[str, Any]] = None,
        metadata: Optional[Dict[str, Any]] = None,
        pages: Optional[Dict[str, Dict[str, Any]]] = None,
    ) -> None:
        self.location = location if location is not None else load(f"loc_{site}")
        self.metadata = metadata if metadata is not None else load(f"meta_{site}_{param_cd}")
        self.continuous = continuous
        self.pages = pages or {}
        self.calls: List[Tuple[str, Optional[Dict[str, Any]]]] = []

    def continuous_calls(self) -> List[Dict[str, Any]]:
        return [p for u, p in self.calls if u.endswith("/continuous/items") and p is not None]

    def __call__(
        self, url: str, params: Optional[Mapping[str, Any]] = None, timeout: int = 60
    ) -> Mock:
        self.calls.append((url, dict(params) if params is not None else None))
        if params is None:
            return _response(self.pages[url])
        if "/monitoring-locations/items/" in url:
            return _response(self.location)
        if url.endswith("/time-series-metadata/items"):
            return _response(self.metadata)
        if url.endswith("/continuous/items"):
            return _response(self.continuous(params))
        raise AssertionError(f"unexpected URL {url}")


def _fixed(name: str) -> ContinuousHandler:
    return lambda params: load(name)


def _feature(ts_id: str, time: str, value: Optional[str], **extra: Any) -> Dict[str, Any]:
    props: Dict[str, Any] = {
        "time_series_id": ts_id,
        "parameter_code": "00060",
        "statistic_id": "00011",
        "time": time,
        "value": value,
        "approval_status": extra.pop("approval_status", "Approved"),
        "qualifier": extra.pop("qualifier", None),
    }
    return {"type": "Feature", "properties": props, "geometry": None}


def _collection(features: List[Dict[str, Any]]) -> Dict[str, Any]:
    return {
        "type": "FeatureCollection",
        "features": features,
        "numberReturned": len(features),
        "links": [],
    }


def _inclusive_server(features: List[Dict[str, Any]]) -> ContinuousHandler:
    """A fake ``continuous`` endpoint honouring the real inclusive-both-ends interval."""

    def handler(params: Mapping[str, Any]) -> Dict[str, Any]:
        lo_s, hi_s = str(params["time"]).split("/")
        lo, hi = pd.Timestamp(lo_s), pd.Timestamp(hi_s)
        keep = [f for f in features if lo <= pd.Timestamp(f["properties"]["time"]) <= hi]
        return _collection(keep)

    return handler


def _run(api: FakeApi, param_cd: str = "00060", **kwargs: Any) -> pd.DataFrame:
    site = kwargs.pop("site", BIG_SANDY)
    with patch("flowfreq.waterdata.requests.get", side_effect=api):
        return download_instantaneous(site, param_cd, **kwargs)


# ----------------------------------------------------------------------------
# Single series, the frame contract
# ----------------------------------------------------------------------------


class TestSingleSeries:
    def _frame(self) -> Tuple[pd.DataFrame, FakeApi]:
        api = FakeApi(BIG_SANDY, "00060", _fixed("iv_03606500_00060_dst"))
        return _run(api, start_date="2024-11-02", end_date="2024-11-03"), api

    def test_contract_matches_legacy_parser(self) -> None:
        df, _ = self._frame()
        legacy = _parse_iv_rdb(IV_BASIC)
        assert list(df.columns) == list(legacy.columns)
        assert str(df.index.tz) == "UTC" and df.index.name == legacy.index.name
        assert df.index.is_monotonic_increasing and df.index.is_unique
        for col in legacy.columns:
            assert df[col].dtype.kind == legacy[col].dtype.kind, col

    def test_values_from_capture(self) -> None:
        df, _ = self._frame()
        assert len(df) == 12
        assert df["flow_cfs"].iloc[0] == pytest.approx(75.8)
        assert df.index[0] == pd.Timestamp("2024-11-03T04:15:00Z")
        assert set(df["qualification_code"]) == {"A"}

    def test_request_shape(self) -> None:
        _, api = self._frame()
        (params,) = api.continuous_calls()
        assert params["monitoring_location_id"] == "USGS-03606500"
        assert params["time_series_id"] == BIG_SANDY_IV_00060
        assert params["parameter_code"] == "00060"
        assert params["limit"] == waterdata.PAGE_LIMIT
        # Local days 2024-11-02..03 in Central time: CDT midnight to CST midnight.
        assert params["time"] == "2024-11-02T05:00:00Z/2024-11-04T06:00:00Z"

    def test_stage_column(self) -> None:
        payload = load("iv_03606500_00060_dst")
        for f in payload["features"]:
            f["properties"]["parameter_code"] = "00065"
            f["properties"]["time_series_id"] = "3ea67b8255054104a0150ca1fd852e99"
        api = FakeApi(BIG_SANDY, "00065", lambda p: payload)
        df = _run(api, "00065", start_date="2024-11-02", end_date="2024-11-03")
        assert list(df.columns) == list(_parse_iv_rdb(IV_STAGE_BASIC, param_cd="00065").columns)

    def test_statistic_filter_ignores_daily_and_annual_series(self) -> None:
        series = self._series()
        assert [s.time_series_id for s in series] == [BIG_SANDY_IV_00060]

    def _series(self) -> List[waterdata.TimeSeriesInfo]:
        api = FakeApi(BIG_SANDY, "00060", _fixed("iv_03606500_00060_empty"))
        with patch("flowfreq.waterdata.requests.get", side_effect=api):
            return waterdata.list_instantaneous_series(BIG_SANDY, "00060")


# ----------------------------------------------------------------------------
# Local time
# ----------------------------------------------------------------------------


class TestDstFallBack:
    """Big Sandy 2024-11-03: 02:00 CDT (07:00Z) fell back to 01:00 CST."""

    @pytest.fixture
    def df(self) -> pd.DataFrame:
        api = FakeApi(BIG_SANDY, "00060", _fixed("iv_03606500_00060_dst"))
        return _run(api, start_date="2024-11-02", end_date="2024-11-03")

    def test_repeated_local_hour(self, df: pd.DataFrame) -> None:
        by_utc = df.set_index(df.index.strftime("%H:%M"))
        assert by_utc.loc["06:00", "datetime_local"] == pd.Timestamp("2024-11-03 01:00")
        assert by_utc.loc["06:00", "tz_cd"] == "CDT"
        assert by_utc.loc["07:00", "datetime_local"] == pd.Timestamp("2024-11-03 01:00")
        assert by_utc.loc["07:00", "tz_cd"] == "CST"
        assert by_utc.loc["04:15", "datetime_local"] == pd.Timestamp("2024-11-02 23:15")

    def test_utc_axis_unique_local_repeats(self, df: pd.DataFrame) -> None:
        assert df.index.is_monotonic_increasing
        # 06:00Z and 07:00Z are both 01:00 local: wall-clock time repeats.
        assert not df["datetime_local"].is_unique

    def test_tz_codes(self, df: pd.DataFrame) -> None:
        before = df.index < pd.Timestamp("2024-11-03T07:00Z")
        assert set(df.loc[before, "tz_cd"]) == {"CDT"}
        assert set(df.loc[~before, "tz_cd"]) == {"CST"}

    def test_local_columns_agree_with_subdaily_conversion(self, df: pd.DataFrame) -> None:
        # flowfreq.subdaily groups on index.tz_convert(tz); datetime_local must agree.
        local = df.index.tz_convert("America/Chicago").tz_localize(None)
        assert (pd.DatetimeIndex(df["datetime_local"]) == local).all()


class TestFixedOffsetZone:
    def test_no_dst_site_is_fixed_offset_year_round(self) -> None:
        location = {"properties": {"time_zone_abbreviation": "MST", "uses_daylight_savings": "N"}}
        features = [
            _feature(BIG_SANDY_IV_00060, "2024-01-15T12:00:00+00:00", "10"),
            _feature(BIG_SANDY_IV_00060, "2024-07-15T12:00:00+00:00", "20"),
        ]
        api = FakeApi(BIG_SANDY, "00060", _inclusive_server(features), location=location)
        df = _run(api, start_date="2024-01-01", end_date="2024-12-31")
        assert list(df["tz_cd"]) == ["MST", "MST"]
        assert list(df["datetime_local"]) == [
            pd.Timestamp("2024-01-15 05:00"),
            pd.Timestamp("2024-07-15 05:00"),
        ]
        # Window is Arizona local midnight, UTC-7 in any season.
        assert api.continuous_calls()[0]["time"].startswith("2024-01-01T07:00:00Z")


class TestResolveLocalZone:
    @pytest.mark.parametrize(
        ("abbr", "zone"),
        [
            ("EST", "America/New_York"),
            ("CST", "America/Chicago"),
            ("MST", "America/Denver"),
            ("PST", "America/Los_Angeles"),
            ("AKST", "America/Anchorage"),
            ("HST", "Pacific/Honolulu"),
            ("AST", "America/Puerto_Rico"),
        ],
    )
    def test_dst_zones(self, abbr: str, zone: str) -> None:
        assert resolve_local_zone(abbr, "Y") == zone

    def test_fixed_offset(self) -> None:
        assert resolve_local_zone("MST", "N") == (timezone(timedelta(hours=-7)), "MST")

    @pytest.mark.parametrize(("abbr", "dst"), [("XST", "Y"), ("XST", "N"), ("SST", "Y")])
    def test_unknown_abbreviation_raises(self, abbr: str, dst: str) -> None:
        with pytest.raises(ValueError, match="Unrecognized time_zone_abbreviation"):
            resolve_local_zone(abbr, dst)

    @pytest.mark.parametrize("dst", [None, "", "maybe"])
    def test_bad_dst_flag_raises(self, dst: Optional[str]) -> None:
        with pytest.raises(ValueError, match="uses_daylight_savings"):
            resolve_local_zone("CST", dst)

    def test_unknown_zone_raises_before_any_data_request(self) -> None:
        location = {"properties": {"time_zone_abbreviation": "XST", "uses_daylight_savings": "Y"}}
        api = FakeApi(BIG_SANDY, "00060", _fixed("iv_03606500_00060_dst"), location=location)
        with pytest.raises(ValueError, match="XST"):
            _run(api, start_date="2024-11-02", end_date="2024-11-03")
        assert api.continuous_calls() == []


# ----------------------------------------------------------------------------
# Multi-sensor sites
# ----------------------------------------------------------------------------


class TestMultiSensor:
    """Olmsted (03612600) stage: HEADWATER and TAILWATER, both flagged primary."""

    def test_refuses_without_ts_id(self) -> None:
        api = FakeApi(OLMSTED, "00065", _fixed("iv_03612600_00065_unfiltered"))
        with pytest.raises(AmbiguousTimeSeriesError) as info:
            _run(api, "00065", site=OLMSTED, start_date="2024-06-01", end_date="2024-06-01")
        message = str(info.value)
        for token in (OLMSTED_HEADWATER, OLMSTED_TAILWATER, "HEADWATER", "TAILWATER"):
            assert token in message
        assert "Headwater" in message and "Primary" in message
        assert "2023-04-13T04:00:00Z" in message  # tailwater begin
        assert api.continuous_calls() == []

    def test_is_a_value_error_like_legacy(self) -> None:
        assert issubclass(AmbiguousTimeSeriesError, ValueError)

    def test_refuses_without_ts_id_or_window(self) -> None:
        api = FakeApi(OLMSTED, "00065", _fixed("iv_03612600_00065_unfiltered"))
        with pytest.raises(AmbiguousTimeSeriesError):
            _run(api, "00065", site=OLMSTED)

    def test_selects_by_uuid(self) -> None:
        api = FakeApi(OLMSTED, "00065", _fixed("iv_03612600_00065_tailwater"))
        df = _run(
            api,
            "00065",
            site=OLMSTED,
            start_date="2024-05-31",
            end_date="2024-06-01",
            ts_id=OLMSTED_TAILWATER,
        )
        assert list(df["gage_height_ft"]) == [40.48, 40.48, 40.46]
        assert api.continuous_calls()[0]["time_series_id"] == OLMSTED_TAILWATER

    def test_uuid_match_ignores_case_and_hyphens(self) -> None:
        dashed = "7CA46507-7C26-415C-93FD-0549C697A270"
        api = FakeApi(OLMSTED, "00065", _fixed("iv_03612600_00065_tailwater"))
        df = _run(
            api, "00065", site=OLMSTED, start_date="2024-05-31", end_date="2024-05-31", ts_id=dashed
        )
        assert len(df) == 3

    def test_rows_from_another_series_are_never_merged(self) -> None:
        # Server ignoring the time_series_id filter: both sensors come back.
        api = FakeApi(OLMSTED, "00065", _fixed("iv_03612600_00065_unfiltered"))
        with pytest.raises(ValueError, match="refusing to merge"):
            _run(
                api,
                "00065",
                site=OLMSTED,
                start_date="2024-06-01",
                end_date="2024-06-01",
                ts_id=OLMSTED_TAILWATER,
            )

    @pytest.mark.parametrize("bad", ["45", "0123456789abcdef0123456789abcdef"])
    def test_unknown_ts_id_raises(self, bad: str) -> None:
        api = FakeApi(OLMSTED, "00065", _fixed("iv_03612600_00065_tailwater"))
        with pytest.raises(ValueError, match="does not match any") as info:
            _run(
                api,
                "00065",
                site=OLMSTED,
                start_date="2024-06-01",
                end_date="2024-06-01",
                ts_id=bad,
            )
        assert not isinstance(info.value, AmbiguousTimeSeriesError)

    def test_window_before_second_sensor_needs_no_ts_id(self) -> None:
        # TAILWATER begins 2023-04-13; a 2020 window has only HEADWATER data.
        features = [_feature(OLMSTED_HEADWATER, "2020-06-01T12:00:00+00:00", "30.0")]
        for f in features:
            f["properties"]["parameter_code"] = "00065"
        api = FakeApi(OLMSTED, "00065", _inclusive_server(features))
        df = _run(api, "00065", site=OLMSTED, start_date="2020-06-01", end_date="2020-06-01")
        assert list(df["gage_height_ft"]) == [30.0]
        assert api.continuous_calls()[0]["time_series_id"] == OLMSTED_HEADWATER


# ----------------------------------------------------------------------------
# Duplicate timestamps within one series
# ----------------------------------------------------------------------------


class TestDuplicateTimestamps:
    def test_identical_duplicates_dropped_with_warning(self, caplog: Any) -> None:
        features = [
            _feature("s", "2024-01-01T00:00:00+00:00", "5"),
            _feature("s", "2024-01-01T00:00:00+00:00", "5"),
            _feature("s", "2024-01-01T00:15:00+00:00", "6"),
        ]
        with caplog.at_level(logging.WARNING, logger="flowfreq.waterdata"):
            df = parse_continuous_features(features, "00060", "s")
        assert list(df["flow_cfs"]) == [5.0, 6.0]
        assert "duplicate" in caplog.text

    def test_conflicting_duplicates_raise(self) -> None:
        features = [
            _feature("s", "2024-01-01T00:00:00+00:00", "5"),
            _feature("s", "2024-01-01T00:00:00+00:00", "7"),
        ]
        with pytest.raises(ValueError, match="conflicting values"):
            parse_continuous_features(features, "00060", "s")

    def test_same_value_different_qualifier_is_a_conflict(self) -> None:
        features = [
            _feature("s", "2024-01-01T00:00:00+00:00", "5"),
            _feature("s", "2024-01-01T00:00:00+00:00", "5", approval_status="Provisional"),
        ]
        with pytest.raises(ValueError, match="conflicting"):
            parse_continuous_features(features, "00060", "s")


# ----------------------------------------------------------------------------
# Windows, chunks, limits
# ----------------------------------------------------------------------------


class TestChunks:
    def test_adjacent_chunks_share_boundary(self) -> None:
        chunks = plan_chunks("2020-01-01", "2022-12-31", 1, "America/Chicago")
        assert len(chunks) == 3
        assert chunks[0][0] == pd.Timestamp("2020-01-01T06:00Z")
        for (_, hi), (lo, _) in zip(chunks, chunks[1:]):
            assert hi == lo
        assert chunks[-1][1] == pd.Timestamp("2023-01-01T06:00Z")

    def test_boundary_instant_counted_once(self) -> None:
        boundary = "2021-01-01T06:00:00+00:00"  # local midnight, the chunk seam
        times = ["2020-12-31T23:00:00+00:00", boundary, "2021-01-01T07:00:00+00:00"]
        features = [_feature(BIG_SANDY_IV_00060, t, str(i)) for i, t in enumerate(times)]
        api = FakeApi(BIG_SANDY, "00060", _inclusive_server(features))
        df = _run(api, start_date="2020-01-01", end_date="2021-12-31")
        assert len(api.continuous_calls()) == 2
        assert len(df) == 3 and df.index.is_unique
        assert (df.index == pd.Timestamp(boundary)).sum() == 1

    def test_three_years_fits_the_envelope(self) -> None:
        chunks = plan_chunks("2020-01-01", "2025-12-31", 3, "America/Chicago")
        assert all(hi - lo <= pd.Timedelta(days=MAX_TIME_ENVELOPE_DAYS) for lo, hi in chunks)

    def test_chunk_over_1100_days_raises_before_any_request(self) -> None:
        api = FakeApi(BIG_SANDY, "00060", _fixed("iv_03606500_00060_dst"))
        with pytest.raises(ValueError, match="1100-day"):
            _run(api, start_date="2015-01-01", end_date="2020-12-31", chunk_years=4)
        assert api.calls == []

    def test_short_window_with_large_chunk_years_is_fine(self) -> None:
        assert len(plan_chunks("2020-01-01", "2020-03-01", 10, "America/Chicago")) == 1

    @pytest.mark.parametrize("years", [0, -1])
    def test_chunk_years_below_one_raises(self, years: int) -> None:
        with pytest.raises(ValueError, match="chunk_years"):
            plan_chunks("2020-01-01", "2020-12-31", years, "America/Chicago")

    def test_reversed_window_raises(self) -> None:
        with pytest.raises(ValueError, match="after"):
            local_day_bounds("2021-01-01", "2020-01-01", "America/Chicago")


class TestPaging:
    def test_follows_next_links(self) -> None:
        page1 = load("iv_03606500_00060_page1")
        page2 = load("iv_03606500_00060_page2")
        page3 = load("iv_03606500_00060_page3")
        pages = {
            waterdata._next_link(page1) or "": page2,
            waterdata._next_link(page2) or "": page3,
        }
        assert waterdata._next_link(page3) is None
        api = FakeApi(BIG_SANDY, "00060", lambda p: page1, pages=pages)
        df = _run(api, start_date="2024-11-02", end_date="2024-11-03")

        whole = load("iv_03606500_00060_dst")
        expected = [f["properties"]["time"] for f in whole["features"]]
        assert [t.isoformat() for t in df.index] == expected
        assert len([c for c in api.calls if c[1] is None]) == 2

    def test_repeated_next_link_raises(self) -> None:
        page = load("iv_03606500_00060_page1")
        api = FakeApi(
            BIG_SANDY, "00060", lambda p: page, pages={waterdata._next_link(page) or "": page}
        )
        with pytest.raises(requests.RequestException, match="paging loop"):
            _run(api, start_date="2024-11-02", end_date="2024-11-03")


class TestEmptyWindow:
    def test_chunk_is_an_ordinary_empty_frame(self) -> None:
        api = FakeApi(BIG_SANDY, "00060", _fixed("iv_03606500_00060_empty"))
        with patch("flowfreq.waterdata.requests.get", side_effect=api):
            part = fetch_continuous_chunk(
                BIG_SANDY,
                "00060",
                BIG_SANDY_IV_00060,
                pd.Timestamp("1990-01-01T00:00Z"),
                pd.Timestamp("1990-01-02T00:00Z"),
            )
        assert part.empty and str(part.index.tz) == "UTC"

    def test_whole_request_empty_raises_like_legacy(self) -> None:
        api = FakeApi(BIG_SANDY, "00060", _fixed("iv_03606500_00060_empty"))
        with pytest.raises(NoInstantaneousDataError, match="1990-01-01"):
            _run(api, start_date="1990-01-01", end_date="1990-01-01")

    def test_http_error_is_not_empty(self) -> None:
        bad = Mock()
        bad.raise_for_status.side_effect = requests.HTTPError("400 time envelope too large")
        with patch("flowfreq.waterdata.requests.get", return_value=bad):
            with pytest.raises(requests.RequestException, match="envelope"):
                fetch_continuous_chunk(
                    BIG_SANDY,
                    "00060",
                    BIG_SANDY_IV_00060,
                    pd.Timestamp("1990-01-01T00:00Z"),
                    pd.Timestamp("1990-01-02T00:00Z"),
                )

    def test_no_series_raises(self) -> None:
        api = FakeApi(BIG_SANDY, "00060", _fixed("iv_03606500_00060_empty"))
        api.metadata = _collection([])
        with pytest.raises(NoInstantaneousDataError, match="lists no instantaneous"):
            _run(api, start_date="2024-01-01", end_date="2024-01-02")


class TestDefaultWindow:
    def test_default_start_is_series_begin(self) -> None:
        api = FakeApi(BIG_SANDY, "00060", _fixed("iv_03606500_00060_empty"))
        with pytest.raises(NoInstantaneousDataError):
            _run(api, end_date="2002-06-30")
        (params,) = api.continuous_calls()
        # begin 2002-01-01T12:00Z -> local day 2002-01-01 -> CST midnight 06:00Z.
        assert params["time"] == "2002-01-01T06:00:00Z/2002-07-01T05:00:00Z"

    def test_default_end_is_series_end(self) -> None:
        api = FakeApi(BIG_SANDY, "00060", _fixed("iv_03606500_00060_empty"))
        with pytest.raises(NoInstantaneousDataError):
            _run(api, start_date="2026-01-01")
        # end 2026-09-25T11:00Z -> local day 2026-09-25 -> CDT midnight 2026-09-26T05:00Z.
        assert api.continuous_calls()[-1]["time"].endswith("/2026-09-26T05:00:00Z")

    def test_always_sends_a_time_window(self) -> None:
        api = FakeApi(BIG_SANDY, "00060", _fixed("iv_03606500_00060_empty"))
        with pytest.raises(NoInstantaneousDataError):
            _run(api, start_date="2026-01-01")
        assert all("time" in p for p in api.continuous_calls())


# ----------------------------------------------------------------------------
# Qualification codes
# ----------------------------------------------------------------------------


class TestQualificationCode:
    @pytest.mark.parametrize(
        ("status", "qualifiers", "expected"),
        [
            ("Approved", None, "A"),
            ("Provisional", None, "P"),
            ("Approved", ["ESTIMATED"], "A:e"),
            ("Provisional", ["ICE"], "P:ICE"),
            ("Provisional", ["ICE", "ESTIMATED"], "P:ICE:e"),
            ("Provisional", "ESTIMATED", "P:e"),
            (None, None, ""),
            ("Revised", [], "Revised"),
        ],
    )
    def test_mapping(self, status: Optional[str], qualifiers: Any, expected: str) -> None:
        assert map_qualification_code(status, qualifiers) == expected

    def test_unmapped_token_logged_at_debug(self, caplog: Any) -> None:
        with caplog.at_level(logging.DEBUG, logger="flowfreq.waterdata"):
            map_qualification_code("Provisional", ["ICE"])
        assert "ICE" in caplog.text

    def test_estimated_capture(self) -> None:
        payload = load("iv_06191500_00060_estimated")
        ts_id = payload["features"][0]["properties"]["time_series_id"]
        df = parse_continuous_features(payload["features"], "00060", ts_id)
        assert list(df["qualification_code"]) == ["A:e"]
        assert list(df["flow_cfs"]) == [858.0]

    def test_ice_rows_without_value_dropped_like_legacy(self, caplog: Any) -> None:
        payload = load("iv_06214500_00060_ice")
        ts_id = payload["features"][0]["properties"]["time_series_id"]
        with caplog.at_level(logging.DEBUG, logger="flowfreq.waterdata"):
            df = parse_continuous_features(payload["features"], "00060", ts_id)
        assert df.empty
        assert "P:ICE" in caplog.text

    def test_ice_with_value_keeps_token(self) -> None:
        features = [
            _feature("s", "2024-01-01T00:00:00+00:00", "12", qualifier=["ICE", "ESTIMATED"])
        ]
        df = parse_continuous_features(features, "00060", "s")
        assert list(df["qualification_code"]) == ["A:ICE:e"]


# ----------------------------------------------------------------------------
# Backend switch on USGSgage
# ----------------------------------------------------------------------------


def _legacy_response(text: str) -> Mock:
    response = Mock()
    response.text = text
    response.status_code = 200
    response.raise_for_status.return_value = None
    return response


class TestBackendSwitch:
    def test_default_is_legacy_and_unchanged(self) -> None:
        with patch("flowfreq.usgs.requests.get", return_value=_legacy_response(IV_BASIC)):
            with patch(
                "flowfreq.waterdata.download_instantaneous",
                side_effect=AssertionError("OGC backend must not run by default"),
            ):
                default = USGSgage("12449950").download_instantaneous_flow(
                    "2022-06-15", "2022-06-15"
                )
                explicit = USGSgage("12449950").download_instantaneous_flow(
                    "2022-06-15", "2022-06-15", backend="nwis-legacy"
                )
        pd.testing.assert_frame_equal(default, explicit)
        pd.testing.assert_frame_equal(default, _parse_iv_rdb(IV_BASIC))

    def test_stage_default_is_legacy(self) -> None:
        with patch("flowfreq.usgs.requests.get", return_value=_legacy_response(IV_STAGE_BASIC)):
            df = USGSgage("12449950").download_instantaneous_stage("2022-06-15", "2022-06-15")
        pd.testing.assert_frame_equal(df, _parse_iv_rdb(IV_STAGE_BASIC, param_cd="00065"))

    def test_waterdata_backend_routes_and_applies_tz(self) -> None:
        api = FakeApi(BIG_SANDY, "00060", _fixed("iv_03606500_00060_dst"))
        gage = USGSgage(BIG_SANDY)
        with patch("flowfreq.waterdata.requests.get", side_effect=api):
            df = gage.download_instantaneous_flow(
                "2024-11-02", "2024-11-03", backend="waterdata-ogc", tz="America/Chicago"
            )
        assert str(df.index.tz) == "America/Chicago"
        assert len(df) == 12
        assert gage._instantaneous_data is df

    def test_waterdata_backend_stage(self) -> None:
        api = FakeApi(OLMSTED, "00065", _fixed("iv_03612600_00065_tailwater"))
        with patch("flowfreq.waterdata.requests.get", side_effect=api):
            df = USGSgage(OLMSTED).download_instantaneous_stage(
                "2024-05-31", "2024-05-31", backend="waterdata-ogc", ts_id=OLMSTED_TAILWATER
            )
        assert "gage_height_ft" in df.columns and len(df) == 3

    def test_unknown_backend_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown instantaneous-value backend"):
            USGSgage(BIG_SANDY).download_instantaneous_flow(
                "2024-01-01", "2024-01-02", backend="nwis"
            )


# ----------------------------------------------------------------------------
# Live
# ----------------------------------------------------------------------------


@pytest.mark.requires_network
class TestLiveWaterData:
    """Live Water Data OGC API checks. Deselected by default."""

    @pytest.mark.parametrize(("method", "col"), [("flow", "flow_cfs"), ("stage", "gage_height_ft")])
    def test_big_sandy_few_days(self, method: str, col: str) -> None:
        gage = USGSgage(BIG_SANDY)
        fetch = getattr(gage, f"download_instantaneous_{method}")
        df = fetch("2024-11-02", "2024-11-04", backend="waterdata-ogc")
        assert list(df.columns) == [col, "datetime_local", "tz_cd", "qualification_code"]
        assert str(df.index.tz) == "UTC" and df.index.is_unique
        assert df.index.is_monotonic_increasing
        assert set(df["tz_cd"]) == {"CDT", "CST"}
        assert (df[col] > 0).all()
