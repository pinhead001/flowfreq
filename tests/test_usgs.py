"""Tests for USGS NWIS data retrieval.

Network access is mocked throughout; no test in this module contacts NWIS.
Tests that would require a live service are marked ``requires_network``.
"""

from __future__ import annotations

from unittest.mock import Mock, patch

import pandas as pd
import pytest
import requests

from hydrolib.flowio import load_flow_frame, save_flow_frame
from hydrolib.usgs import (
    NWIS_TZ_OFFSETS,
    NoInstantaneousDataError,
    USGSgage,
    _chunk_date_range,
    _is_no_data_response,
    _parse_iv_rdb,
)
from tests.fixtures.nwis_rdb import (
    IV_BASIC,
    IV_DST_FALL_BACK,
    IV_EMPTY,
    IV_MULTI_SENSOR,
    IV_NO_DATA_400_BODY,
    IV_UNKNOWN_TZ,
    IV_WITH_GAPS,
    IV_WRONG_PARAMETER,
    SITE_SERIES_CATALOG,
    SITE_SERIES_CATALOG_NO_UV,
)


def _mock_response(text: str, status_code: int = 200) -> Mock:
    """Build a stand-in for a requests Response carrying RDB text."""
    response = Mock()
    response.text = text
    response.status_code = status_code
    if status_code >= 400:
        response.raise_for_status.side_effect = requests.HTTPError(f"{status_code} error")
    else:
        response.raise_for_status.return_value = None
    return response


class TestChunkDateRange:
    """Tests for splitting a long request into year-sized chunks."""

    def test_single_chunk_when_range_is_short(self) -> None:
        """A range shorter than the chunk size stays as one request."""
        assert _chunk_date_range("2022-06-01", "2022-09-30", 1) == [("2022-06-01", "2022-09-30")]

    def test_splits_multi_year_range(self) -> None:
        """A three-year range becomes three one-year chunks."""
        chunks = _chunk_date_range("2020-01-01", "2022-12-31", 1)
        assert chunks == [
            ("2020-01-01", "2020-12-31"),
            ("2021-01-01", "2021-12-31"),
            ("2022-01-01", "2022-12-31"),
        ]

    def test_chunks_are_gapless_and_non_overlapping(self) -> None:
        """Consecutive chunks abut exactly: no day is requested twice or missed."""
        chunks = _chunk_date_range("2015-03-17", "2021-08-02", 1)
        for (_, end), (next_start, _) in zip(chunks, chunks[1:]):
            assert pd.Timestamp(next_start) - pd.Timestamp(end) == pd.Timedelta(days=1)

    def test_chunks_cover_exactly_the_requested_range(self) -> None:
        """The union of the chunks is the requested window, not more or less."""
        chunks = _chunk_date_range("2015-03-17", "2021-08-02", 2)
        assert chunks[0][0] == "2015-03-17"
        assert chunks[-1][1] == "2021-08-02"

    def test_multi_year_chunk_size(self) -> None:
        """chunk_years > 1 produces correspondingly longer chunks."""
        chunks = _chunk_date_range("2020-01-01", "2023-12-31", 2)
        assert chunks == [("2020-01-01", "2021-12-31"), ("2022-01-01", "2023-12-31")]

    def test_reversed_range_raises(self) -> None:
        """An end date before the start date is an error, not an empty list."""
        with pytest.raises(ValueError, match="after end_date"):
            _chunk_date_range("2022-12-31", "2022-01-01", 1)


class TestParseInstantaneousRDB:
    """Tests for the instantaneous-value RDB parser."""

    def test_parses_records(self) -> None:
        """A basic payload yields one row per record with the documented columns."""
        df = _parse_iv_rdb(IV_BASIC)
        assert len(df) == 6
        assert list(df.columns) == [
            "flow_cfs",
            "datetime_local",
            "tz_cd",
            "qualification_code",
        ]

    def test_index_is_timezone_aware_utc(self) -> None:
        """The index carries a timezone rather than being naive."""
        df = _parse_iv_rdb(IV_BASIC)
        assert df.index.tz is not None
        assert str(df.index.tz) == "UTC"

    def test_local_time_converted_by_reported_offset(self) -> None:
        """12:00 PDT is 19:00 UTC — the tz_cd offset is actually applied."""
        df = _parse_iv_rdb(IV_BASIC)
        assert df.index[0] == pd.Timestamp("2022-06-15 19:00", tz="UTC")

    def test_local_time_and_zone_preserved_verbatim(self) -> None:
        """Nothing is discarded: the reported wall-clock time survives as a column."""
        df = _parse_iv_rdb(IV_BASIC)
        assert df["datetime_local"].iloc[0] == pd.Timestamp("2022-06-15 12:00")
        assert df["tz_cd"].iloc[0] == "PDT"

    def test_flow_is_float(self) -> None:
        """Discharge is float even when every value in the chunk is a whole number."""
        df = _parse_iv_rdb(IV_BASIC)
        assert df["flow_cfs"].dtype == float
        assert df["flow_cfs"].iloc[0] == pytest.approx(1520.0)

    def test_qualification_code_captured(self) -> None:
        """The NWIS data qualifier is retained per record."""
        df = _parse_iv_rdb(IV_BASIC)
        assert df["qualification_code"].iloc[0] == "A"
        assert df["qualification_code"].iloc[-1] == "P"

    def test_dst_transition_gives_monotonic_utc_index(self) -> None:
        """Across the autumn DST change local time repeats but UTC keeps increasing.

        This is the reason the index is not left in local time: on 2022-11-06
        the wall clock runs 01:45 PDT then 01:00 PST, so a naive index moves
        backwards and any per-day statistic computed from it is quietly wrong.
        """
        df = _parse_iv_rdb(IV_DST_FALL_BACK)
        assert not df["datetime_local"].is_monotonic_increasing
        assert df.index.is_monotonic_increasing

    def test_dst_transition_preserves_distinct_instants(self) -> None:
        """The repeated local hour maps to two different UTC instants, not one."""
        df = _parse_iv_rdb(IV_DST_FALL_BACK)
        repeated = df[df["datetime_local"] == pd.Timestamp("2022-11-06 01:15")]
        assert len(repeated) == 2
        assert repeated.index.nunique() == 2
        assert sorted(repeated["tz_cd"]) == ["PDT", "PST"]

    def test_non_numeric_values_dropped(self) -> None:
        """'Ice' and blank readings are dropped rather than becoming NaN rows."""
        df = _parse_iv_rdb(IV_WITH_GAPS)
        assert len(df) == 2
        assert df["flow_cfs"].tolist() == [210.0, 205.0]

    def test_empty_payload_returns_empty_frame(self) -> None:
        """A header with no records is an empty frame, not an exception."""
        df = _parse_iv_rdb(IV_EMPTY)
        assert df.empty
        assert list(df.columns) == [
            "flow_cfs",
            "datetime_local",
            "tz_cd",
            "qualification_code",
        ]

    def test_multiple_sensors_raises(self) -> None:
        """Two 00060 series at one site must not be resolved by silently picking one."""
        with pytest.raises(ValueError, match="separate 00060 time series"):
            _parse_iv_rdb(IV_MULTI_SENSOR)

    def test_ts_id_selects_one_sensor(self) -> None:
        """ts_id disambiguates a multi-sensor site."""
        df = _parse_iv_rdb(IV_MULTI_SENSOR, ts_id="63680")
        assert df["flow_cfs"].tolist() == [1490.0, 1495.0]

    def test_unmatched_ts_id_raises(self) -> None:
        """A ts_id that matches nothing is an error, not an empty result."""
        with pytest.raises(ValueError, match="does not match any discharge series"):
            _parse_iv_rdb(IV_MULTI_SENSOR, ts_id="99999")

    def test_missing_discharge_column_raises(self) -> None:
        """Records present but no 00060 column is a schema surprise worth raising on."""
        with pytest.raises(ValueError, match="no discharge \\(00060\\) column"):
            _parse_iv_rdb(IV_WRONG_PARAMETER)

    def test_unknown_timezone_code_raises(self) -> None:
        """An unmappable tz_cd raises rather than dropping those records."""
        with pytest.raises(ValueError, match="Unrecognized NWIS time-zone code"):
            _parse_iv_rdb(IV_UNKNOWN_TZ)


class TestTimezoneOffsets:
    """Tests for the NWIS time-zone abbreviation table."""

    def test_covers_conterminous_us_zones(self) -> None:
        """All four conterminous US zones are present in both standard and daylight."""
        for code in ("PST", "PDT", "MST", "MDT", "CST", "CDT", "EST", "EDT"):
            assert code in NWIS_TZ_OFFSETS

    def test_daylight_offset_is_one_hour_ahead_of_standard(self) -> None:
        """Daylight codes sit exactly one hour closer to UTC than their standard twin."""
        for standard, daylight in (("PST", "PDT"), ("MST", "MDT"), ("EST", "EDT")):
            assert NWIS_TZ_OFFSETS[daylight] == NWIS_TZ_OFFSETS[standard] + 1


class TestNoDataResponseDetection:
    """Tests for distinguishing an empty window from a broken request."""

    def test_recognizes_nwis_no_data_body(self) -> None:
        """The service's 400-with-explanation body is recognized as 'no records'."""
        assert _is_no_data_response(IV_NO_DATA_400_BODY)

    def test_does_not_match_arbitrary_error(self) -> None:
        """An unrelated error body is not mistaken for an empty window."""
        assert not _is_no_data_response("Internal Server Error")


class TestDownloadInstantaneousFlow:
    """Tests for USGSgage.download_instantaneous_flow."""

    def test_returns_frame_and_caches_it(self) -> None:
        """The frame is returned and also stored on the gage."""
        gage = USGSgage("12449950")
        with patch("hydrolib.usgs.requests.get", return_value=_mock_response(IV_BASIC)):
            df = gage.download_instantaneous_flow("2022-06-15", "2022-06-15")

        assert len(df) == 6
        assert gage.instantaneous_data is not None
        assert len(gage.instantaneous_data) == 6

    def test_requests_the_instantaneous_endpoint(self) -> None:
        """The unit-value service is used, with no statistic code."""
        gage = USGSgage("12449950")
        with patch("hydrolib.usgs.requests.get", return_value=_mock_response(IV_BASIC)) as mock_get:
            gage.download_instantaneous_flow("2022-06-15", "2022-06-15")

        url = mock_get.call_args.args[0]
        params = mock_get.call_args.kwargs["params"]
        assert url == USGSgage.BASE_URL_IV
        assert params["parameterCd"] == "00060"
        assert "statCd" not in params

    def test_chunks_long_request_by_year(self) -> None:
        """A three-year window is issued as three requests with abutting windows."""
        gage = USGSgage("12449950")
        with patch("hydrolib.usgs.requests.get", return_value=_mock_response(IV_BASIC)) as mock_get:
            gage.download_instantaneous_flow("2020-01-01", "2022-12-31")

        assert mock_get.call_count == 3
        windows = [
            (c.kwargs["params"]["startDT"], c.kwargs["params"]["endDT"])
            for c in mock_get.call_args_list
        ]
        assert windows == [
            ("2020-01-01", "2020-12-31"),
            ("2021-01-01", "2021-12-31"),
            ("2022-01-01", "2022-12-31"),
        ]

    def test_duplicate_timestamps_collapsed(self) -> None:
        """Chunks returning the same instant twice yield one row, not two."""
        gage = USGSgage("12449950")
        with patch("hydrolib.usgs.requests.get", return_value=_mock_response(IV_BASIC)):
            df = gage.download_instantaneous_flow("2020-01-01", "2022-12-31")

        assert len(df) == 6
        assert df.index.is_unique
        assert df.index.is_monotonic_increasing

    def test_empty_chunk_skipped_but_others_kept(self) -> None:
        """A window NWIS has no records for is a gap, not a failure of the whole call."""
        gage = USGSgage("12449950")
        responses = [
            _mock_response(IV_NO_DATA_400_BODY, status_code=400),
            _mock_response(IV_BASIC),
        ]
        with patch("hydrolib.usgs.requests.get", side_effect=responses):
            df = gage.download_instantaneous_flow("2021-01-01", "2022-12-31")

        assert len(df) == 6

    def test_no_data_anywhere_raises_rather_than_returning_empty(self) -> None:
        """A site with no instantaneous record fails loudly, per the brief."""
        gage = USGSgage("12449950")
        with patch(
            "hydrolib.usgs.requests.get",
            return_value=_mock_response(IV_NO_DATA_400_BODY, status_code=400),
        ):
            with pytest.raises(NoInstantaneousDataError, match="No instantaneous discharge"):
                gage.download_instantaneous_flow("2022-01-01", "2022-12-31")

    def test_failed_chunk_raises_naming_the_window(self) -> None:
        """A transport failure never comes back as a silently truncated record."""
        gage = USGSgage("12449950")
        responses = [
            _mock_response(IV_BASIC),
            requests.ConnectionError("connection reset"),
        ]
        with patch("hydrolib.usgs.requests.get", side_effect=responses):
            with pytest.raises(requests.RequestException, match="2022-01-01 to 2022-12-31"):
                gage.download_instantaneous_flow("2021-01-01", "2022-12-31")

    def test_server_error_raises(self) -> None:
        """A 500 is a failure, not an empty window."""
        gage = USGSgage("12449950")
        with patch(
            "hydrolib.usgs.requests.get",
            return_value=_mock_response("Internal Server Error", status_code=500),
        ):
            with pytest.raises(requests.RequestException):
                gage.download_instantaneous_flow("2022-01-01", "2022-06-30")

    def test_tz_argument_converts_index(self) -> None:
        """Passing tz returns the index in that zone, same instants."""
        gage = USGSgage("12449950")
        with patch("hydrolib.usgs.requests.get", return_value=_mock_response(IV_BASIC)):
            df = gage.download_instantaneous_flow(
                "2022-06-15", "2022-06-15", tz="America/Los_Angeles"
            )

        assert str(df.index.tz) == "America/Los_Angeles"
        assert df.index[0] == pd.Timestamp("2022-06-15 12:00", tz="America/Los_Angeles")

    def test_invalid_chunk_years_raises(self) -> None:
        """A non-positive chunk size is rejected before any request is made."""
        gage = USGSgage("12449950")
        with pytest.raises(ValueError, match="chunk_years must be >= 1"):
            gage.download_instantaneous_flow("2022-01-01", "2022-12-31", chunk_years=0)

    def test_defaults_to_instantaneous_period_of_record(self) -> None:
        """With no dates given, the site's unit-value POR bounds the request."""
        gage = USGSgage("12449950")
        gage.site_name = "METHOW RIVER AT PATEROS, WA"
        gage.drainage_area = 1772.0

        def _dispatch(url, **kwargs):
            if url == USGSgage.BASE_URL_SITE:
                return _mock_response(SITE_SERIES_CATALOG)
            return _mock_response(IV_BASIC)

        with patch("hydrolib.usgs.requests.get", side_effect=_dispatch):
            gage.download_instantaneous_flow()

        assert gage.iv_por_start == "2007-10-01"
        assert gage.iv_por_end == "2024-09-30"

    def test_site_without_uv_series_fails_before_requesting_data(self) -> None:
        """A site NWIS lists no unit-value series for never gets a doomed data request."""
        gage = USGSgage("12449950")
        gage.site_name = "METHOW RIVER AT PATEROS, WA"
        gage.drainage_area = 1772.0

        with patch(
            "hydrolib.usgs.requests.get", return_value=_mock_response(SITE_SERIES_CATALOG_NO_UV)
        ) as mock_get:
            with pytest.raises(NoInstantaneousDataError, match="no instantaneous"):
                gage.download_instantaneous_flow()

        assert all(c.args[0] == USGSgage.BASE_URL_SITE for c in mock_get.call_args_list)


class TestFlowFrameIO:
    """Tests for saving and reloading a retrieved flow series."""

    def test_csv_round_trip(self, tmp_path) -> None:
        """CSV preserves the values and the instants."""
        df = _parse_iv_rdb(IV_BASIC)
        path = save_flow_frame(df, tmp_path / "iv.csv")
        loaded = load_flow_frame(path)

        assert len(loaded) == len(df)
        assert loaded["flow_cfs"].tolist() == df["flow_cfs"].tolist()
        assert loaded.index[0] == df.index[0]

    def test_parquet_round_trip(self, tmp_path) -> None:
        """Parquet restores dtypes and the tz-aware index exactly."""
        pytest.importorskip("pyarrow")
        df = _parse_iv_rdb(IV_BASIC)
        path = save_flow_frame(df, tmp_path / "iv.parquet")
        loaded = load_flow_frame(path)

        pd.testing.assert_frame_equal(loaded, df)

    def test_unsupported_extension_raises(self, tmp_path) -> None:
        """An unknown extension is refused rather than guessed at."""
        df = _parse_iv_rdb(IV_BASIC)
        with pytest.raises(ValueError, match="Unsupported flow-frame format"):
            save_flow_frame(df, tmp_path / "iv.xlsx")


@pytest.mark.requires_network
class TestLiveNWIS:
    """Live NWIS checks. Deselect with -m 'not requires_network'."""

    def test_download_instantaneous_flow_live(self) -> None:
        """Retrieve a short real window from NWIS."""
        gage = USGSgage("12449950")
        df = gage.download_instantaneous_flow("2022-06-01", "2022-06-07")

        assert not df.empty
        assert df.index.tz is not None
        assert (df["flow_cfs"] > 0).all()
