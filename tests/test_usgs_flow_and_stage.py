"""Tests for paired discharge-and-stage retrieval.

:func:`flowfreq.usgs.join_flow_and_stage` and
:meth:`flowfreq.usgs.USGSgage.download_instantaneous_flow_and_stage`. The
offline tests reuse the captured-shape IV fixtures in
``tests/fixtures/nwis_rdb.py``. ``IV_BASIC`` (discharge, six 15-minute readings)
and ``IV_STAGE_BASIC`` (stage, the first four marks, one of them ``Ice``) are
misaligned in the way real pairs are, so the outer-join rule is exercised on
fixture data rather than only on hand-built frames.
"""

from __future__ import annotations

from typing import Any
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from flowfreq.usgs import (
    PAIRED_IV_COLUMNS,
    NoInstantaneousDataError,
    USGSgage,
    _parse_iv_rdb,
    join_flow_and_stage,
)
from tests.fixtures.nwis_rdb import IV_BASIC, IV_DST_FALL_BACK, IV_NO_DATA_400_BODY, IV_STAGE_BASIC
from tests.test_usgs import _mock_response


def _flow() -> pd.DataFrame:
    return _parse_iv_rdb(IV_BASIC)


def _stage() -> pd.DataFrame:
    return _parse_iv_rdb(IV_STAGE_BASIC, param_cd="00065")


def _by_parameter(flow_text: str, stage_text: str, stage_status: int = 200) -> Any:
    """A requests.get stand-in answering 00060 and 00065 with different payloads."""

    def get(url: str, params: Any = None, **kwargs: Any) -> Any:
        if params["parameterCd"] == "00065":
            return _mock_response(stage_text, status_code=stage_status)
        return _mock_response(flow_text)

    return get


class TestJoinFlowAndStage:
    """The join rule itself, on fixture-parsed frames."""

    def test_outer_join_keeps_every_instant_without_interpolating(self) -> None:
        joined = join_flow_and_stage(_flow(), _stage())
        assert list(joined.columns) == list(PAIRED_IV_COLUMNS)
        assert len(joined) == 6  # the union; discharge has every mark here
        assert joined.index.is_monotonic_increasing
        # 12:30 PDT stage was "Ice", so it is not in the stage frame; 13:00 and
        # 13:15 are past the end of the stage record. All three stay NaN.
        stage_missing = joined.index[joined["gage_height_ft"].isna()]
        assert [t.strftime("%H:%M") for t in stage_missing] == ["19:30", "20:00", "20:15"]
        np.testing.assert_allclose(joined["gage_height_ft"].dropna().to_numpy(), [4.52, 4.55, 4.48])
        assert joined["flow_cfs"].notna().all()

    def test_readings_at_different_instants_are_separate_rows(self) -> None:
        """A 12:07 stage and a 12:00 discharge are two rows, never merged."""
        flow = _flow().iloc[:1]
        stage = _stage().iloc[:1].copy()
        stage.index = stage.index + pd.Timedelta(minutes=7)
        stage["datetime_local"] = stage["datetime_local"] + pd.Timedelta(minutes=7)
        joined = join_flow_and_stage(flow, stage)
        assert len(joined) == 2
        assert joined["flow_cfs"].notna().sum() == 1
        assert joined["gage_height_ft"].notna().sum() == 1
        assert joined.dropna(subset=["flow_cfs", "gage_height_ft"]).empty

    def test_qualification_codes_kept_per_parameter(self) -> None:
        joined = join_flow_and_stage(_flow(), _stage())
        at_1245 = joined.loc[pd.Timestamp("2022-06-15 19:45", tz="UTC")]
        assert at_1245["qualification_code_flow"] == "A"
        assert at_1245["qualification_code_stage"] == "P"
        assert "qualification_code" not in joined.columns

    def test_local_time_filled_from_whichever_frame_has_the_instant(self) -> None:
        joined = join_flow_and_stage(_flow().iloc[:2], _stage())
        assert joined["datetime_local"].notna().all()
        assert (joined["tz_cd"] == "PDT").all()
        assert joined["datetime_local"].iloc[-1] == pd.Timestamp("2022-06-15 12:45")

    def test_mismatched_tz_cd_is_refused(self) -> None:
        """Across the fall-back hour, PDT vs PST is exactly the disagreement
        that must not be resolved by picking one."""
        flow = _parse_iv_rdb(IV_DST_FALL_BACK)
        stage = flow.rename(columns={"flow_cfs": "gage_height_ft"}).copy()
        stage["gage_height_ft"] = 3.0
        stage.loc[stage.index[3], "tz_cd"] = "PDT"
        with pytest.raises(ValueError, match="disagree on tz_cd at 1 shared"):
            join_flow_and_stage(flow, stage)

    def test_mismatched_datetime_local_is_refused(self) -> None:
        stage = _stage()
        stage.loc[stage.index[0], "datetime_local"] = pd.Timestamp("2022-06-15 11:00")
        with pytest.raises(ValueError, match="datetime_local"):
            join_flow_and_stage(_flow(), stage)

    def test_matching_fall_back_hour_joins_cleanly(self) -> None:
        flow = _parse_iv_rdb(IV_DST_FALL_BACK)
        stage = flow.rename(columns={"flow_cfs": "gage_height_ft"})
        joined = join_flow_and_stage(flow, stage)
        assert len(joined) == 6
        assert list(joined["tz_cd"]) == ["PDT"] * 3 + ["PST"] * 3

    def test_index_takes_the_flow_frames_zone_and_joins_across_zones(self) -> None:
        flow = _flow()
        flow.index = flow.index.tz_convert("America/Los_Angeles")
        joined = join_flow_and_stage(flow, _stage())  # stage left in UTC
        assert str(joined.index.tz) == "America/Los_Angeles"
        assert len(joined) == 6
        assert joined["gage_height_ft"].notna().sum() == 3

    def test_naive_index_raises(self) -> None:
        stage = _stage()
        stage.index = stage.index.tz_localize(None)
        with pytest.raises(TypeError, match="stage frame has a naive index"):
            join_flow_and_stage(_flow(), stage)

    def test_missing_column_raises(self) -> None:
        with pytest.raises(KeyError, match="gage_height_ft"):
            join_flow_and_stage(_flow(), _flow())

    def test_duplicate_timestamps_raise(self) -> None:
        flow = _flow()
        with pytest.raises(ValueError, match="duplicate timestamps"):
            join_flow_and_stage(pd.concat([flow, flow.iloc[:1]]), _stage())


class TestDownloadInstantaneousFlowAndStage:
    """The download method, with the network mocked."""

    def test_legacy_requests_both_parameters_and_joins(self) -> None:
        gage = USGSgage("12449950")
        with patch(
            "flowfreq.usgs.requests.get", side_effect=_by_parameter(IV_BASIC, IV_STAGE_BASIC)
        ) as mock_get:
            both = gage.download_instantaneous_flow_and_stage(
                "2022-06-15", "2022-06-15", backend="nwis-legacy"
            )
        codes = [c.kwargs["params"]["parameterCd"] for c in mock_get.call_args_list]
        assert codes == ["00060", "00065"]
        pd.testing.assert_frame_equal(both, join_flow_and_stage(_flow(), _stage()))
        # The unjoined frames are cached as the single-parameter calls cache them.
        assert "flow_cfs" in gage.instantaneous_data.columns
        assert "gage_height_ft" in gage.instantaneous_stage.columns

    def test_ts_ids_go_to_their_own_parameter(self) -> None:
        gage = USGSgage("12449950")
        with patch(
            "flowfreq.usgs.requests.get", side_effect=_by_parameter(IV_BASIC, IV_STAGE_BASIC)
        ):
            with patch.object(
                USGSgage, "_download_instantaneous_backend", autospec=True
            ) as backend:
                backend.side_effect = lambda self, p, **kw: (_flow() if p == "00060" else _stage())
                gage.download_instantaneous_flow_and_stage(
                    "2022-06-15",
                    "2022-06-15",
                    ts_id_flow="45",
                    ts_id_stage="46",
                    backend="nwis-legacy",
                )
        by_param = {c.args[1]: c.kwargs["ts_id"] for c in backend.call_args_list}
        assert by_param == {"00060": "45", "00065": "46"}

    def test_waterdata_backend_is_used_for_both(self) -> None:
        gage = USGSgage("12449950")
        with patch(
            "flowfreq.waterdata.download_instantaneous",
            side_effect=lambda site, p, **kw: _flow() if p == "00060" else _stage(),
        ) as ogc:
            both = gage.download_instantaneous_flow_and_stage(
                "2022-06-15", "2022-06-15", backend="waterdata-ogc", tz="America/Los_Angeles"
            )
        assert [c.args[1] for c in ogc.call_args_list] == ["00060", "00065"]
        assert str(both.index.tz) == "America/Los_Angeles"
        assert len(both) == 6

    def test_wrong_form_stage_ts_id_refused_before_any_request(self) -> None:
        gage = USGSgage("12449950")
        with patch("flowfreq.usgs.requests.get") as mock_get:
            with pytest.raises(ValueError, match="NWIS DD number"):
                gage.download_instantaneous_flow_and_stage(
                    "2022-06-15", "2022-06-15", ts_id_stage="60629", backend="waterdata-ogc"
                )
        mock_get.assert_not_called()

    def test_no_stage_record_raises_rather_than_returning_flow_only(self) -> None:
        gage = USGSgage("12449950")
        with patch(
            "flowfreq.usgs.requests.get",
            side_effect=_by_parameter(IV_BASIC, IV_NO_DATA_400_BODY, stage_status=400),
        ):
            with pytest.raises(NoInstantaneousDataError, match="00065"):
                gage.download_instantaneous_flow_and_stage(
                    "2022-06-15", "2022-06-15", backend="nwis-legacy"
                )


@pytest.mark.requires_network
class TestLiveFlowAndStage:
    """One live pair on the default (Water Data OGC) backend."""

    def test_methow_at_pateros_week(self) -> None:
        """12449950 publishes both 00060 and 00065 at 15-minute marks."""
        gage = USGSgage("12449950")
        both = gage.download_instantaneous_flow_and_stage(
            "2024-06-01", "2024-06-03", tz="America/Los_Angeles"
        )
        assert list(both.columns) == list(PAIRED_IV_COLUMNS)
        assert both.index.is_unique and both.index.is_monotonic_increasing
        paired = both.dropna(subset=["flow_cfs", "gage_height_ft"])
        assert len(paired) > 0.9 * len(both)
        assert (paired["flow_cfs"] > 0).all()
        # More water, higher stage: the two series must move together.
        assert paired["flow_cfs"].corr(paired["gage_height_ft"]) > 0.9
