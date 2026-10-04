"""Tests for USGS NWIS data retrieval.

Network access is mocked throughout; no test in this module contacts NWIS.
Tests that would require a live service are marked ``requires_network``.
"""

from __future__ import annotations

import json
from datetime import date, datetime, timezone
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd
import pytest
import requests

from flowfreq.flowio import load_flow_frame, save_flow_frame
from flowfreq.peak_sources import DEFAULT_BACKEND
from flowfreq.usgs import (
    NWIS_TZ_OFFSETS,
    NoInstantaneousDataError,
    USGSgage,
    _chunk_date_range,
    _is_no_data_response,
    _parse_iv_rdb,
    _parse_peak_dt,
    fetch_nwis_batch,
    fetch_nwis_peaks,
)
from tests.fixtures.nwis_rdb import (
    DV_BASIC,
    DV_SINGLE_DAY,
    IV_BASIC,
    IV_DST_FALL_BACK,
    IV_EMPTY,
    IV_MULTI_SENSOR,
    IV_NO_DATA_400_BODY,
    IV_STAGE_BASIC,
    IV_STAGE_MULTI_SENSOR,
    IV_UNKNOWN_TZ,
    IV_WITH_GAPS,
    IV_WRONG_PARAMETER,
    PEAK_PARTIAL_DATES,
    SITE_EXPANDED,
    SITE_EXPANDED_NO_COORDS,
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


class TestParseIvRdbStage:
    """Tests for _parse_iv_rdb on gage height (parameter 00065).

    The RDB parser is shared with discharge, so these check the parts that are
    actually parameterized -- column discovery and output naming -- rather than
    re-testing the tz mapping that TestParseIvRdb already covers.
    """

    def test_stage_column_named_gage_height_ft(self) -> None:
        df = _parse_iv_rdb(IV_STAGE_BASIC, param_cd="00065")
        assert list(df.columns) == [
            "gage_height_ft",
            "datetime_local",
            "tz_cd",
            "qualification_code",
        ]
        assert df["gage_height_ft"].dtype == float
        assert df["gage_height_ft"].iloc[0] == pytest.approx(4.52)

    def test_stage_shares_the_utc_axis_rule(self) -> None:
        """12:00 PDT is 19:00 UTC for stage exactly as for discharge."""
        df = _parse_iv_rdb(IV_STAGE_BASIC, param_cd="00065")
        assert df.index[0] == pd.Timestamp("2022-06-15 19:00", tz="UTC")
        assert df["datetime_local"].iloc[0] == pd.Timestamp("2022-06-15 12:00")

    def test_stage_non_numeric_dropped(self) -> None:
        """An 'Ice' stage reading is dropped, not carried as NaN."""
        df = _parse_iv_rdb(IV_STAGE_BASIC, param_cd="00065")
        assert len(df) == 3

    def test_stage_multiple_sensors_raises(self) -> None:
        """Separate primary and backup stage sensors must not be resolved by
        silently picking one."""
        with pytest.raises(ValueError, match="separate 00065 time series"):
            _parse_iv_rdb(IV_STAGE_MULTI_SENSOR, param_cd="00065")

    def test_stage_ts_id_selects_one_sensor(self) -> None:
        df = _parse_iv_rdb(IV_STAGE_MULTI_SENSOR, ts_id="63681", param_cd="00065")
        assert df["gage_height_ft"].tolist() == [4.51, 4.54]

    def test_stage_unmatched_ts_id_raises_naming_the_parameter(self) -> None:
        with pytest.raises(ValueError, match="does not match any gage height series"):
            _parse_iv_rdb(IV_STAGE_MULTI_SENSOR, ts_id="99999", param_cd="00065")

    def test_discharge_payload_has_no_stage_column(self) -> None:
        """Asking a discharge-only payload for 00065 raises rather than
        returning an empty frame that reads as 'the site has no stage'."""
        with pytest.raises(ValueError, match=r"no gage height \(00065\) column"):
            _parse_iv_rdb(IV_BASIC, param_cd="00065")

    def test_unsupported_parameter_code_raises(self) -> None:
        with pytest.raises(ValueError, match="Unsupported instantaneous parameter code"):
            _parse_iv_rdb(IV_BASIC, param_cd="00010")

    def test_discharge_default_unchanged(self) -> None:
        """The default parameter is still discharge, so existing callers that
        pass no param_cd are unaffected by the generalization."""
        assert "flow_cfs" in _parse_iv_rdb(IV_BASIC).columns


class TestDownloadInstantaneousStage:
    """Tests for USGSgage.download_instantaneous_stage."""

    def test_requests_parameter_00065(self) -> None:
        gage = USGSgage("12449950")
        with patch(
            "flowfreq.usgs.requests.get", return_value=_mock_response(IV_STAGE_BASIC)
        ) as mock_get:
            gage.download_instantaneous_stage("2022-06-15", "2022-06-15", backend="nwis-legacy")

        params = mock_get.call_args.kwargs["params"]
        assert mock_get.call_args.args[0] == USGSgage.BASE_URL_IV
        assert params["parameterCd"] == "00065"
        assert "statCd" not in params

    def test_cached_separately_from_discharge(self) -> None:
        """A stage download must not land in instantaneous_data, or a caller
        reading instantaneous_data["flow_cfs"] afterwards gets a KeyError from
        a call it never made.
        """
        gage = USGSgage("12449950")
        with patch("flowfreq.usgs.requests.get", return_value=_mock_response(IV_BASIC)):
            gage.download_instantaneous_flow("2022-06-15", "2022-06-15", backend="nwis-legacy")
        with patch("flowfreq.usgs.requests.get", return_value=_mock_response(IV_STAGE_BASIC)):
            stage = gage.download_instantaneous_stage(
                "2022-06-15", "2022-06-15", backend="nwis-legacy"
            )

        assert "flow_cfs" in gage.instantaneous_data.columns
        assert "gage_height_ft" in gage.instantaneous_stage.columns
        assert gage.instantaneous_stage is stage

    def test_chunks_and_converts_tz_like_discharge(self) -> None:
        gage = USGSgage("12449950")
        with patch(
            "flowfreq.usgs.requests.get", return_value=_mock_response(IV_STAGE_BASIC)
        ) as mock_get:
            df = gage.download_instantaneous_stage(
                "2020-01-01", "2021-12-31", tz="America/Los_Angeles", backend="nwis-legacy"
            )

        assert mock_get.call_count == 2
        assert str(df.index.tz) == "America/Los_Angeles"

    def test_no_stage_record_raises_naming_the_parameter(self) -> None:
        """Stage is absent far more often than discharge, so the error must not
        read as 'this site has no unit values at all'."""
        gage = USGSgage("12449950")
        with patch(
            "flowfreq.usgs.requests.get",
            return_value=_mock_response(IV_NO_DATA_400_BODY, status_code=400),
        ):
            with pytest.raises(NoInstantaneousDataError, match="00065"):
                gage.download_instantaneous_stage("2022-06-15", "2022-06-15", backend="nwis-legacy")

    def test_bad_chunk_years_raises(self) -> None:
        gage = USGSgage("12449950")
        with pytest.raises(ValueError, match="chunk_years must be >= 1"):
            gage.download_instantaneous_stage(
                "2022-06-15", "2022-06-15", chunk_years=0, backend="nwis-legacy"
            )


class TestDownloadInstantaneousFlow:
    """Tests for USGSgage.download_instantaneous_flow."""

    def test_returns_frame_and_caches_it(self) -> None:
        """The frame is returned and also stored on the gage."""
        gage = USGSgage("12449950")
        with patch("flowfreq.usgs.requests.get", return_value=_mock_response(IV_BASIC)):
            df = gage.download_instantaneous_flow("2022-06-15", "2022-06-15", backend="nwis-legacy")

        assert len(df) == 6
        assert gage.instantaneous_data is not None
        assert len(gage.instantaneous_data) == 6

    def test_requests_the_instantaneous_endpoint(self) -> None:
        """The unit-value service is used, with no statistic code."""
        gage = USGSgage("12449950")
        with patch("flowfreq.usgs.requests.get", return_value=_mock_response(IV_BASIC)) as mock_get:
            gage.download_instantaneous_flow("2022-06-15", "2022-06-15", backend="nwis-legacy")

        url = mock_get.call_args.args[0]
        params = mock_get.call_args.kwargs["params"]
        assert url == USGSgage.BASE_URL_IV
        assert params["parameterCd"] == "00060"
        assert "statCd" not in params

    def test_chunks_long_request_by_year(self) -> None:
        """A three-year window is issued as three requests with abutting windows."""
        gage = USGSgage("12449950")
        with patch("flowfreq.usgs.requests.get", return_value=_mock_response(IV_BASIC)) as mock_get:
            gage.download_instantaneous_flow("2020-01-01", "2022-12-31", backend="nwis-legacy")

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
        with patch("flowfreq.usgs.requests.get", return_value=_mock_response(IV_BASIC)):
            df = gage.download_instantaneous_flow("2020-01-01", "2022-12-31", backend="nwis-legacy")

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
        with patch("flowfreq.usgs.requests.get", side_effect=responses):
            df = gage.download_instantaneous_flow("2021-01-01", "2022-12-31", backend="nwis-legacy")

        assert len(df) == 6

    def test_no_data_anywhere_raises_rather_than_returning_empty(self) -> None:
        """A site with no instantaneous record fails loudly, per the brief."""
        gage = USGSgage("12449950")
        with patch(
            "flowfreq.usgs.requests.get",
            return_value=_mock_response(IV_NO_DATA_400_BODY, status_code=400),
        ):
            with pytest.raises(NoInstantaneousDataError, match="No instantaneous discharge"):
                gage.download_instantaneous_flow("2022-01-01", "2022-12-31", backend="nwis-legacy")

    def test_failed_chunk_raises_naming_the_window(self) -> None:
        """A transport failure never comes back as a silently truncated record."""
        gage = USGSgage("12449950")
        responses = [
            _mock_response(IV_BASIC),
            requests.ConnectionError("connection reset"),
        ]
        with patch("flowfreq.usgs.requests.get", side_effect=responses):
            with pytest.raises(requests.RequestException, match="2022-01-01 to 2022-12-31"):
                gage.download_instantaneous_flow("2021-01-01", "2022-12-31", backend="nwis-legacy")

    def test_server_error_raises(self) -> None:
        """A 500 is a failure, not an empty window."""
        gage = USGSgage("12449950")
        with patch(
            "flowfreq.usgs.requests.get",
            return_value=_mock_response("Internal Server Error", status_code=500),
        ):
            with pytest.raises(requests.RequestException):
                gage.download_instantaneous_flow("2022-01-01", "2022-06-30", backend="nwis-legacy")

    def test_tz_argument_converts_index(self) -> None:
        """Passing tz returns the index in that zone, same instants."""
        gage = USGSgage("12449950")
        with patch("flowfreq.usgs.requests.get", return_value=_mock_response(IV_BASIC)):
            df = gage.download_instantaneous_flow(
                "2022-06-15", "2022-06-15", tz="America/Los_Angeles", backend="nwis-legacy"
            )

        assert str(df.index.tz) == "America/Los_Angeles"
        assert df.index[0] == pd.Timestamp("2022-06-15 12:00", tz="America/Los_Angeles")

    def test_invalid_chunk_years_raises(self) -> None:
        """A non-positive chunk size is rejected before any request is made."""
        gage = USGSgage("12449950")
        with pytest.raises(ValueError, match="chunk_years must be >= 1"):
            gage.download_instantaneous_flow(
                "2022-01-01", "2022-12-31", chunk_years=0, backend="nwis-legacy"
            )

    def test_defaults_to_instantaneous_period_of_record(self) -> None:
        """With no dates given, the site's unit-value POR bounds the request."""
        gage = USGSgage("12449950")
        gage.site_name = "METHOW RIVER AT PATEROS, WA"
        gage.drainage_area = 1772.0

        def _dispatch(url, **kwargs):
            if url == USGSgage.BASE_URL_SITE:
                return _mock_response(SITE_SERIES_CATALOG)
            return _mock_response(IV_BASIC)

        with patch("flowfreq.usgs.requests.get", side_effect=_dispatch):
            gage.download_instantaneous_flow(backend="nwis-legacy")

        assert gage.iv_por_start == "2007-10-01"
        assert gage.iv_por_end == "2024-09-30"

    def test_site_without_uv_series_fails_before_requesting_data(self) -> None:
        """A site NWIS lists no unit-value series for never gets a doomed data request."""
        gage = USGSgage("12449950")
        gage.site_name = "METHOW RIVER AT PATEROS, WA"
        gage.drainage_area = 1772.0

        with patch(
            "flowfreq.usgs.requests.get", return_value=_mock_response(SITE_SERIES_CATALOG_NO_UV)
        ) as mock_get:
            with pytest.raises(NoInstantaneousDataError, match="no instantaneous"):
                gage.download_instantaneous_flow(backend="nwis-legacy")

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
        df = _parse_iv_rdb(IV_BASIC)
        path = save_flow_frame(df, tmp_path / "iv.parquet")
        loaded = load_flow_frame(path)

        pd.testing.assert_frame_equal(loaded, df)

    def test_parquet_preserves_timezone_where_csv_does_not(self, tmp_path) -> None:
        """The reason Parquet is preferred: CSV loses the index dtype, Parquet keeps it."""
        df = _parse_iv_rdb(IV_BASIC)

        from_parquet = load_flow_frame(save_flow_frame(df, tmp_path / "iv.parquet"))
        from_csv = load_flow_frame(save_flow_frame(df, tmp_path / "iv.csv"))

        assert from_parquet.index.dtype == df.index.dtype
        assert from_parquet["flow_cfs"].dtype == df["flow_cfs"].dtype
        assert from_parquet["datetime_local"].dtype == df["datetime_local"].dtype
        # CSV agrees on the instants but not on how they are typed.
        assert list(from_csv.index) == list(df.index)
        assert from_csv["datetime_local"].dtype != df["datetime_local"].dtype

    def test_parquet_compression_is_configurable(self, tmp_path) -> None:
        """A different codec still round-trips; the default is not load-bearing."""
        df = _parse_iv_rdb(IV_BASIC)
        path = save_flow_frame(df, tmp_path / "iv.parquet", compression="snappy")
        pd.testing.assert_frame_equal(load_flow_frame(path), df)

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
        df = gage.download_instantaneous_flow("2022-06-01", "2022-06-07", backend="nwis-legacy")

        assert not df.empty
        assert df.index.tz is not None
        assert (df["flow_cfs"] > 0).all()


class TestSiteCoordinates:
    """Decimal-degree coordinates from the siteOutput=expanded call."""

    def test_parses_latitude_and_longitude(self) -> None:
        gage = USGSgage("03606500")
        with patch("flowfreq.usgs.requests.get", return_value=_mock_response(SITE_EXPANDED)):
            gage.fetch_site_info(backend="nwis-legacy")

        assert gage.latitude == pytest.approx(36.0389722)
        assert gage.longitude == pytest.approx(-88.2450000)

    def test_longitude_stays_negative_in_western_hemisphere(self) -> None:
        """Sign is carried through as NWIS reports it, not normalised to positive."""
        gage = USGSgage("03606500")
        with patch("flowfreq.usgs.requests.get", return_value=_mock_response(SITE_EXPANDED)):
            gage.fetch_site_info(backend="nwis-legacy")

        assert gage.longitude < 0

    def test_missing_coordinates_leave_none(self) -> None:
        """Empty NWIS values must not become NaN, 0.0, or an exception."""
        gage = USGSgage("03606500")
        with patch(
            "flowfreq.usgs.requests.get", return_value=_mock_response(SITE_EXPANDED_NO_COORDS)
        ):
            gage.fetch_site_info(backend="nwis-legacy")

        assert gage.latitude is None
        assert gage.longitude is None
        # The rest of the response is still parsed.
        assert gage.site_name == "BIG SANDY RIVER AT BRUCETON, TN"

    def test_default_to_none_before_any_fetch(self) -> None:
        gage = USGSgage("03606500")
        assert gage.latitude is None
        assert gage.longitude is None

    def test_blank_drainage_area_is_none_not_nan(self) -> None:
        """The same NaN trap the coordinates fall into, on the adjacent field.

        NWIS reports a missing numeric as an empty column, pandas reads it as
        NaN, and float(nan) does not raise -- so a naive conversion stores NaN.
        """
        gage = USGSgage("03606500")
        rdb = SITE_EXPANDED_NO_COORDS.replace("\t205\n", "\t\n")
        with patch("flowfreq.usgs.requests.get", return_value=_mock_response(rdb)):
            gage.fetch_site_info(backend="nwis-legacy")

        assert gage.drainage_area is None

    def test_drainage_area_still_parses_when_present(self) -> None:
        gage = USGSgage("03606500")
        with patch("flowfreq.usgs.requests.get", return_value=_mock_response(SITE_EXPANDED)):
            gage.fetch_site_info(backend="nwis-legacy")

        assert gage.drainage_area == pytest.approx(205.0)


class TestDownloadDailyFlow:
    """A range-less NWIS request returns one day, not the record.

    These exist because the omission was invisible: no exception, no warning,
    a well-formed one-row frame that builds a degenerate curve wherever it is
    used. The bug lived in a method with no test, in a module whose live tests
    are deselected by default.
    """

    def test_a_date_range_is_sent_even_when_the_caller_supplies_none(self) -> None:
        with patch("flowfreq.usgs.requests.get") as get:
            get.return_value = _mock_response(DV_BASIC)
            USGSgage("12449500").download_daily_flow(backend="nwis-legacy")

        params = get.call_args.kwargs["params"]
        assert "startDT" in params and "endDT" in params
        assert params["startDT"] and params["endDT"]

    def test_the_default_range_spans_the_period_of_record(self) -> None:
        with patch("flowfreq.usgs.requests.get") as get:
            get.return_value = _mock_response(DV_BASIC)
            USGSgage("12449500").download_daily_flow(backend="nwis-legacy")

        params = get.call_args.kwargs["params"]
        assert params["startDT"] == USGSgage.DEFAULT_START_DATE
        assert params["endDT"] == datetime.now(timezone.utc).date().isoformat()

    def test_default_end_date_uses_utc_not_local_clock(self) -> None:
        """A host clock behind UTC must not silently narrow the requested range.

        Regression test for exactly that: the default used to be computed from
        ``date.today()`` (local wall-clock), so a host behind UTC could omit a
        day NWIS had already published with no error to signal it.
        """
        fixed_now = datetime(2026, 1, 1, 0, 30, tzinfo=timezone.utc)
        with (
            patch("flowfreq.usgs.requests.get") as get,
            patch("flowfreq.usgs.datetime") as mock_datetime,
        ):
            mock_datetime.now.return_value = fixed_now
            get.return_value = _mock_response(DV_BASIC)
            USGSgage("12449500").download_daily_flow(backend="nwis-legacy")

        assert get.call_args.kwargs["params"]["endDT"] == "2026-01-01"
        mock_datetime.now.assert_called_with(timezone.utc)

    def test_the_default_start_precedes_every_usgs_daily_record(self) -> None:
        assert date.fromisoformat(USGSgage.DEFAULT_START_DATE).year <= 1850

    def test_caller_supplied_dates_are_passed_through(self) -> None:
        with patch("flowfreq.usgs.requests.get") as get:
            get.return_value = _mock_response(DV_BASIC)
            USGSgage("12449500").download_daily_flow(
                start_date="1919-06-01", end_date="2025-09-30", backend="nwis-legacy"
            )

        params = get.call_args.kwargs["params"]
        assert params["startDT"] == "1919-06-01"
        assert params["endDT"] == "2025-09-30"

    def test_default_timeout_is_generous(self) -> None:
        """60s, not the 30s a small siteOutput=expanded request is sized for --

        a full period-of-record request (this method's own default) can be tens
        of thousands of RDB rows for a long-running, high-frequency site.
        """
        with patch("flowfreq.usgs.requests.get") as get:
            get.return_value = _mock_response(DV_BASIC)
            USGSgage("12449500").download_daily_flow(backend="nwis-legacy")

        assert get.call_args.kwargs["timeout"] == 60

    def test_timeout_is_configurable(self) -> None:
        with patch("flowfreq.usgs.requests.get") as get:
            get.return_value = _mock_response(DV_BASIC)
            USGSgage("12449500").download_daily_flow(timeout=120, backend="nwis-legacy")

        assert get.call_args.kwargs["timeout"] == 120

    def test_reversed_range_raises(self) -> None:
        """A caller-supplied reversed range must fail before ever reaching NWIS."""
        with patch("flowfreq.usgs.requests.get") as get:
            with pytest.raises(ValueError, match="after end_date"):
                USGSgage("12449500").download_daily_flow(
                    start_date="2025-01-01", end_date="2020-01-01", backend="nwis-legacy"
                )

        get.assert_not_called()

    def test_malformed_date_raises(self) -> None:
        with patch("flowfreq.usgs.requests.get") as get:
            with pytest.raises(ValueError):
                USGSgage("12449500").download_daily_flow(
                    start_date="not-a-date", backend="nwis-legacy"
                )

        get.assert_not_called()

    def test_parses_a_daily_record(self) -> None:
        with patch("flowfreq.usgs.requests.get") as get:
            get.return_value = _mock_response(DV_BASIC)
            frame = USGSgage("12449500").download_daily_flow(backend="nwis-legacy")

        assert len(frame) == 5
        assert frame["flow_cfs"].iloc[0] == 221
        assert frame["flow_cfs"].iloc[-1] == 216

    def test_the_single_day_response_is_what_the_old_default_produced(self) -> None:
        """Documents the defect rather than the fix: this payload is a valid
        200 carrying one row. Nothing downstream can tell it from a record,
        which is why the range is now always sent."""
        with patch("flowfreq.usgs.requests.get") as get:
            get.return_value = _mock_response(DV_SINGLE_DAY)
            frame = USGSgage("12449500").download_daily_flow(backend="nwis-legacy")

        assert len(frame) == 1


class TestDownloadPeakFlowPartialDates:
    """NWIS writes an unknown day or month as ``00``; those peaks must survive.

    ``pd.to_datetime(..., errors="coerce")`` turned ``1897-03-00`` into NaT,
    which made ``water_year`` NaN and dropped the row -- silently losing the
    historic peaks, found while verifying the Water Data OGC API (#29).
    """

    def _download(self) -> pd.DataFrame:
        with patch("flowfreq.usgs.requests.get") as get:
            get.return_value = _mock_response(PEAK_PARTIAL_DATES)
            return USGSgage("03606500").download_peak_flow(backend="nwis-legacy")

    def test_no_peak_is_dropped(self) -> None:
        assert len(self._download()) == 6

    def test_historic_peaks_with_unknown_day_are_kept(self) -> None:
        frame = self._download().set_index("water_year")
        assert frame.loc[1897, "peak_flow_cfs"] == 25000
        assert frame.loc[1919, "peak_flow_cfs"] == 21000
        assert frame.loc[1927, "peak_flow_cfs"] == 18500
        assert list(frame.loc[[1897, 1919, 1927], "qualification_code"]) == ["7", "7", "7"]

    def test_water_years_follow_peakfq(self) -> None:
        """Unknown month -> calendar year; unknown day in October -> next year."""
        assert list(self._download()["water_year"]) == [1897, 1919, 1927, 1930, 1931, 1932]

    def test_placeholder_dates_use_the_first_of_the_month(self) -> None:
        dates = self._download()["peak_date"]
        assert dates.iloc[0] == pd.Timestamp("1897-03-01")
        assert dates.iloc[4] == pd.Timestamp("1931-01-01")

    def test_numeric_only_codes_stay_strings(self) -> None:
        """An all-numeric ``peak_cd`` column must not be read as float ("7.0")."""
        codes = self._download()["qualification_code"]
        assert list(codes) == ["7", "7", "7", "", "", ""]

    def test_a_full_date_is_unchanged(self) -> None:
        assert self._download()["peak_date"].iloc[3] == pd.Timestamp("1930-01-09")

    def test_header_metadata_is_set(self) -> None:
        with patch("flowfreq.usgs.requests.get") as get:
            get.return_value = _mock_response(PEAK_PARTIAL_DATES)
            gage = USGSgage("03606500")
            gage.download_peak_flow(backend="nwis-legacy")
        assert gage.site_name == "Big Sandy River at Bruceton, TN"
        assert gage.drainage_area == 205.0


_FIXTURES = Path(__file__).parent / "fixtures"


class _RoutedGet:
    """``requests.get`` stand-in that answers by URL prefix, in order per prefix.

    ``flowfreq.usgs``, ``flowfreq.peak_sources`` and ``flowfreq.waterdata``
    all call the one ``requests.get``, so a single patch has to serve every
    service a download touches -- and record which it touched.
    """

    def __init__(self, routes: dict) -> None:
        self.routes = {prefix: list(answers) for prefix, answers in routes.items()}
        self.urls: list = []

    def __call__(self, url, params=None, timeout=None):
        assert timeout is not None, "requests must carry a timeout"
        self.urls.append(url)
        for prefix, answers in self.routes.items():
            if url.startswith(prefix):
                answer = answers.pop(0)
                if isinstance(answer, Exception):
                    raise answer
                response = Mock()
                response.raise_for_status.return_value = None
                if isinstance(answer, str):
                    response.text = answer
                else:
                    response.json.return_value = answer
                return response
        raise AssertionError(f"unexpected request to {url}")


def _ogc_pages() -> list:
    data = json.loads((_FIXTURES / "waterdata_peaks_03606500.json").read_text(encoding="utf-8"))
    return [data["page1"], data["page2"]]


def _ogc_location() -> dict:
    return json.loads(
        (_FIXTURES / "waterdata_ogc" / "loc_03606500.json").read_text(encoding="utf-8")
    )


_PEAKS_OGC = "https://api.waterdata.usgs.gov/ogcapi/v1/collections/peaks/items"
_LOCATIONS_OGC = "https://api.waterdata.usgs.gov/ogcapi/v1/collections/monitoring-locations"


class TestDownloadPeakFlowBackends:
    """``download_peak_flow`` routes through :mod:`flowfreq.peak_sources` (#29)."""

    def test_default_is_the_water_data_api(self) -> None:
        get = _RoutedGet({_PEAKS_OGC: _ogc_pages(), _LOCATIONS_OGC: [_ogc_location()]})
        gage = USGSgage("03606500")
        with patch("flowfreq.usgs.requests.get", get):
            frame = gage.download_peak_flow()
        assert not any(u.startswith(USGSgage.BASE_URL_PEAKS) for u in get.urls)
        assert list(frame.columns) == [
            "water_year",
            "peak_date",
            "peak_flow_cfs",
            "qualification_code",
        ]
        assert list(frame["water_year"]) == [1897, 1919, 1927, 1930, 1934, 2010, 2020]
        assert frame is gage.peak_data
        assert gage.period_of_record == (1897, 2020)

    def test_default_sets_site_metadata_from_monitoring_locations(self) -> None:
        get = _RoutedGet({_PEAKS_OGC: _ogc_pages(), _LOCATIONS_OGC: [_ogc_location()]})
        gage = USGSgage("03606500")
        with patch("flowfreq.usgs.requests.get", get):
            gage.download_peak_flow()
        assert gage.site_name == "BIG SANDY RIVER AT BRUCETON, TN"
        assert gage.drainage_area == 205.0

    def test_peak_path_sets_every_site_attribute(self) -> None:
        """Not just name and area: the attributes fetch_site_info sets from the
        same record (HUC, FIPS state, location). loc_12449500.json is the one
        untrimmed monitoring-locations capture."""
        feature = json.loads(
            (_FIXTURES / "waterdata_ogc" / "loc_12449500.json").read_text(encoding="utf-8")
        )
        gage = USGSgage("12449500")
        with patch("flowfreq.waterdata.fetch_monitoring_location_feature", return_value=feature):
            gage._site_metadata_from_waterdata()
        assert gage.drainage_area == 1301.0
        assert gage.huc == "170200080610"
        assert gage.state_code == "53"
        assert gage.latitude == pytest.approx(48.3651, abs=1e-4)
        assert gage.longitude == pytest.approx(-120.1162, abs=1e-4)
        assert gage.site_name

    def test_metadata_failure_still_returns_the_peaks(self, caplog) -> None:
        get = _RoutedGet(
            {_PEAKS_OGC: _ogc_pages(), _LOCATIONS_OGC: [requests.ConnectionError("down")]}
        )
        gage = USGSgage("03606500")
        with patch("flowfreq.usgs.requests.get", get), caplog.at_level("WARNING"):
            frame = gage.download_peak_flow()
        assert len(frame) == 7
        assert gage.site_name is None and gage.drainage_area is None
        assert gage.huc is None and gage.latitude is None
        assert "monitoring-location request failed" in caplog.text

    def test_same_frame_shape_on_both_backends(self) -> None:
        ogc = _RoutedGet({_PEAKS_OGC: _ogc_pages(), _LOCATIONS_OGC: [_ogc_location()]})
        with patch("flowfreq.usgs.requests.get", ogc):
            new = USGSgage("03606500").download_peak_flow()
        with patch("flowfreq.usgs.requests.get", return_value=_mock_response(PEAK_PARTIAL_DATES)):
            old = USGSgage("03606500").download_peak_flow(backend="nwis-legacy")
        assert list(new.columns) == list(old.columns)
        for col in new.columns:
            assert new[col].dtype.kind == old[col].dtype.kind, col
        assert isinstance(new.index, pd.RangeIndex) and isinstance(old.index, pd.RangeIndex)

    def test_legacy_backend_reads_the_rdb_service(self) -> None:
        with patch(
            "flowfreq.usgs.requests.get", return_value=_mock_response(PEAK_PARTIAL_DATES)
        ) as get:
            USGSgage("03606500").download_peak_flow(backend="nwis-legacy")
        assert get.call_args.args[0] == USGSgage.BASE_URL_PEAKS

    def test_unknown_backend_raises(self) -> None:
        with patch("flowfreq.usgs.requests.get") as get:
            with pytest.raises(ValueError, match="Unknown peak backend 'nope'"):
                USGSgage("03606500").download_peak_flow(backend="nope")
        get.assert_not_called()

    def test_fetch_nwis_peaks_passes_the_backend(self) -> None:
        seen = []

        def _fake(self, backend="unset"):
            seen.append(backend)
            self._peak_data = pd.DataFrame(
                {"water_year": [2000], "peak_flow_cfs": [10.0], "qualification_code": [""]}
            )
            return self._peak_data

        with patch.object(USGSgage, "download_peak_flow", _fake):
            assert fetch_nwis_peaks("03606500") == [{"year": 2000, "flow": 10.0, "source": "USGS"}]
            fetch_nwis_peaks("03606500", backend="nwis-legacy")
            results, errors = fetch_nwis_batch(["03606500"], workers=1, backend="nwis-legacy")
        assert seen == [DEFAULT_BACKEND, "nwis-legacy", "nwis-legacy"]
        assert set(results) == {"03606500"} and errors == {}

    def test_fetch_nwis_batch_reports_a_bad_backend_per_site(self) -> None:
        results, errors = fetch_nwis_batch(["03606500"], workers=1, backend="nope")
        assert results == {}
        assert "Unknown peak backend" in errors["03606500"]


@pytest.mark.requires_network
def test_live_default_download_is_the_water_data_api() -> None:
    """The new default, live: peaks plus site metadata from the OGC API."""
    gage = USGSgage("03606500")
    frame = gage.download_peak_flow()
    assert list(frame.columns) == ["water_year", "peak_date", "peak_flow_cfs", "qualification_code"]
    by_year = frame.set_index("water_year")
    assert by_year.loc[1897, "peak_flow_cfs"] == 25000
    assert "7" in by_year.loc[1897, "qualification_code"]
    assert gage.site_name == "BIG SANDY RIVER AT BRUCETON, TN"
    assert gage.drainage_area == pytest.approx(205.0)


class TestParsePeakDt:
    def test_garbage_is_nat_not_an_error(self) -> None:
        parsed = _parse_peak_dt(pd.Series(["1950-05-12", "not a date", "", "1950-13-01"]))
        assert parsed.iloc[0] == pd.Timestamp("1950-05-12")
        assert parsed.iloc[1:].isna().all()

    def test_partial_dates_are_logged(self, caplog: pytest.LogCaptureFixture) -> None:
        with caplog.at_level("INFO", logger="flowfreq.usgs"):
            _parse_peak_dt(pd.Series(["1897-03-00", "1968-00-00"]), "03606500")
        assert "2 peak date(s) with unknown day or month" in caplog.text
