"""Tests for USGS StreamStats delineation and basin-characteristics retrieval.

Network access is mocked throughout; no test in this module contacts StreamStats.
Tests that would require the live service are marked ``requires_network``.
"""

from __future__ import annotations

import copy
import json as json_module
import math
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import Mock, patch

import pytest
import requests

from flowfreq.streamstats import (
    MAX_CONCURRENCY,
    POLYGON_AREA_RTOL,
    SQ_M_PER_SQ_MI,
    Characteristic,
    DegenerateDelineationError,
    FlowStatisticEstimate,
    Provenance,
    RegionFlowEstimates,
    RegressionCitation,
    SnapResult,
    StreamStatsCache,
    StreamStatsResponseError,
    StreamStatsTransportError,
    UnsnappablePointError,
    UnsupportedRegionError,
    WatershedCharacteristics,
    _fill_and_validate_regions,
    _point_in_geometry,
    batch_estimate_flow_statistics,
    batch_get_characteristics,
    delineate_and_get_characteristics,
    estimate_flow_statistics,
    geojson_area_sq_mi,
    list_regions,
    list_statistic_groups,
    locate_regression_regions,
    snap_point,
)
from tests.fixtures.streamstats_responses import (
    DELINEATE_SSHYDRO_GA_OGEECHEE,
    DELINEATE_SSHYDRO_GOOD,
    DELINEATE_SSHYDRO_MALFORMED,
    DELINEATE_SSHYDRO_METHOW_PATEROS,
    DELINEATE_SSHYDRO_WITH_WARNING,
    HYDRO_CHARACTERISTICS_GOOD,
    HYDRO_CHARACTERISTICS_MALFORMED,
    HYDRO_CHARACTERISTICS_METHOW_PATEROS,
    HYDRO_CHARACTERISTICS_MISSING_VALUE,
    HYDRO_CHARACTERISTICS_WRAPPED,
    NSS_BYLOCATION_GA_OGEECHEE,
    NSS_BYLOCATION_WA_GOAT_CREEK,
    NSS_CITATIONS_RESPONSE,
    NSS_ESTIMATE_500_BODY,
    NSS_ESTIMATE_GA_AREA_AVERAGED,
    NSS_ESTIMATE_RESPONSE_PFS,
    NSS_ESTIMATE_WA_GOAT_CREEK_ALL_PFS,
    NSS_ESTIMATE_WA_GOAT_CREEK_GC1751,
    NSS_SCENARIO_TEMPLATE_PFS,
    NSS_SCENARIOS_GA_RURAL_UNDER_1,
    NSS_SCENARIOS_WA_PFS_LFS,
    NSS_STATISTIC_GROUPS_ALL,
    NSS_STATISTIC_GROUPS_WA,
    RESPONSE_422_MISSING_REGION,
    SNAP_GOOD,
    SNAP_METHOW_PATEROS,
    SNAP_UNSNAPPABLE,
    load_capture,
)


def _mock_response(json_data=None, status_code: int = 200, headers=None, text=None) -> Mock:
    """Build a stand-in for a requests Response carrying a JSON (or raw text) body."""
    response = Mock()
    response.status_code = status_code
    response.headers = headers or {}
    response.url = "https://mocked.example/"
    if text is not None:
        response.text = text
    elif json_data is not None:
        response.text = json_module.dumps(json_data)
    else:
        response.text = ""

    if json_data is not None:
        response.json.return_value = json_data
    else:
        response.json.side_effect = ValueError("No JSON object could be decoded")

    if status_code >= 400:
        response.raise_for_status.side_effect = requests.HTTPError(f"{status_code} error")
    else:
        response.raise_for_status.return_value = None
    return response


def _good_get_responses(server: str = "PRODWEBB"):
    """The two GET calls (snap, sshydro) for a fully successful point."""
    return [
        _mock_response(SNAP_GOOD),
        _mock_response(DELINEATE_SSHYDRO_GOOD, headers={"usgswim-hostname": server}),
    ]


class TestSnapPoint:
    """Tests for FR-1: snap before delineating, and never silently past a refusal."""

    def test_good_snap_returns_result(self) -> None:
        with patch("flowfreq.streamstats.requests.get", return_value=_mock_response(SNAP_GOOD)):
            result = snap_point("WA", 48.57430, -120.37890)

        assert result.could_snap is True
        assert result.snapped_lat == pytest.approx(48.57426)
        assert result.snapped_lon == pytest.approx(-120.37893)

    def test_distance_m_is_nonnegative(self) -> None:
        with patch("flowfreq.streamstats.requests.get", return_value=_mock_response(SNAP_GOOD)):
            result = snap_point("WA", 48.57430, -120.37890)

        assert result.distance_m >= 0

    def test_unsnappable_point_raises(self) -> None:
        with patch(
            "flowfreq.streamstats.requests.get", return_value=_mock_response(SNAP_UNSNAPPABLE)
        ):
            with pytest.raises(UnsnappablePointError):
                snap_point("WA", 48.584, -120.370)

    def test_422_raises_unsupported_region(self) -> None:
        with patch(
            "flowfreq.streamstats.requests.get",
            return_value=_mock_response(status_code=422, text=RESPONSE_422_MISSING_REGION),
        ):
            with pytest.raises(UnsupportedRegionError):
                snap_point("ZZ", 48.5, -120.3)

    def test_malformed_json_raises_response_error(self) -> None:
        with patch(
            "flowfreq.streamstats.requests.get",
            return_value=_mock_response(status_code=200, text="not json"),
        ):
            with pytest.raises(StreamStatsResponseError):
                snap_point("WA", 48.5, -120.3)

    def test_missing_output_coordinates_raises(self) -> None:
        payload = {"couldSnap": True, "output": {}}
        with patch("flowfreq.streamstats.requests.get", return_value=_mock_response(payload)):
            with pytest.raises(StreamStatsResponseError):
                snap_point("WA", 48.5, -120.3)


class TestRequestWithBackoff:
    """Tests for NFR-2: polite retries on transient failure, none on a data problem."""

    def test_retries_on_5xx_then_succeeds(self) -> None:
        responses = [
            _mock_response(status_code=500, text="Internal Server Error"),
            _mock_response(SNAP_GOOD),
        ]
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=responses),
            patch("flowfreq.streamstats.time.sleep", return_value=None),
        ):
            result = snap_point("WA", 48.57430, -120.37890)

        assert result.could_snap is True

    def test_exhausts_retries_and_raises_transport_error(self) -> None:
        with (
            patch(
                "flowfreq.streamstats.requests.get",
                return_value=_mock_response(status_code=503, text="Service Unavailable"),
            ),
            patch("flowfreq.streamstats.time.sleep", return_value=None),
        ):
            with pytest.raises(StreamStatsTransportError):
                snap_point("WA", 48.5, -120.3)

    def test_connection_error_retried_then_succeeds(self) -> None:
        with (
            patch(
                "flowfreq.streamstats.requests.get",
                side_effect=[requests.ConnectionError("reset"), _mock_response(SNAP_GOOD)],
            ),
            patch("flowfreq.streamstats.time.sleep", return_value=None),
        ):
            result = snap_point("WA", 48.57430, -120.37890)

        assert result.could_snap is True

    def test_4xx_is_not_retried(self) -> None:
        with patch("flowfreq.streamstats.requests.get") as mock_get:
            mock_get.return_value = _mock_response(status_code=422, text="bad request")
            with pytest.raises(UnsupportedRegionError):
                snap_point("ZZ", 48.5, -120.3)

        assert mock_get.call_count == 1

    def test_user_agent_header_identifies_the_client(self) -> None:
        with patch(
            "flowfreq.streamstats.requests.get", return_value=_mock_response(SNAP_GOOD)
        ) as mock_get:
            snap_point("WA", 48.57430, -120.37890)

        assert "flowfreq" in mock_get.call_args.kwargs["headers"]["User-Agent"]

    def test_timeout_is_always_passed(self) -> None:
        with patch(
            "flowfreq.streamstats.requests.get", return_value=_mock_response(SNAP_GOOD)
        ) as mock_get:
            snap_point("WA", 48.57430, -120.37890, timeout=12.5)

        assert mock_get.call_args.kwargs["timeout"] == 12.5


class TestDelineateAndGetCharacteristics:
    """Tests for the full snap -> delineate -> characteristics pipeline."""

    def test_happy_path_returns_characteristics(self) -> None:
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_good_get_responses()),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(HYDRO_CHARACTERISTICS_GOOD),
            ),
        ):
            result = delineate_and_get_characteristics("WA", 48.57430, -120.37890)

        assert result["DRNAREA"].value == pytest.approx(412.0)
        assert result["PRECPRIS10"].value == pytest.approx(45.62)
        assert "CANOPY_PCT" in result

    def test_server_stickiness_from_delineate_header(self) -> None:
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_good_get_responses("PRODWEBB")),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(HYDRO_CHARACTERISTICS_GOOD),
            ) as mock_post,
        ):
            result = delineate_and_get_characteristics("WA", 48.57430, -120.37890)

        assert result.provenance is not None
        assert result.provenance.server_used == "prodwebb"
        assert "prodwebb.streamstats.usgs.gov" in mock_post.call_args.args[0]

    def test_unsnappable_point_never_reaches_delineate(self) -> None:
        with (
            patch(
                "flowfreq.streamstats.requests.get", return_value=_mock_response(SNAP_UNSNAPPABLE)
            ) as mock_get,
            patch("flowfreq.streamstats.requests.post") as mock_post,
        ):
            with pytest.raises(UnsnappablePointError):
                delineate_and_get_characteristics("WA", 48.584, -120.370)

        assert mock_get.call_count == 1
        mock_post.assert_not_called()

    def test_warning_msg_on_delineation_raises_and_skips_hydro(self) -> None:
        """The mandatory test (design doc S8): this is the whole point of the module.

        A 200 carrying a WarningMsg for an unsnappable point must never reach
        ss-hydro, let alone be returned as if it were a real result.
        """
        responses = [
            _mock_response(SNAP_GOOD),
            _mock_response(DELINEATE_SSHYDRO_WITH_WARNING),
        ]
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=responses),
            patch("flowfreq.streamstats.requests.post") as mock_post,
        ):
            with pytest.raises(DegenerateDelineationError, match="not snappable"):
                delineate_and_get_characteristics("WA", 48.584, -120.370)

        mock_post.assert_not_called()

    def test_malformed_sshydro_response_raises(self) -> None:
        responses = [
            _mock_response(SNAP_GOOD),
            _mock_response(DELINEATE_SSHYDRO_MALFORMED),
        ]
        with patch("flowfreq.streamstats.requests.get", side_effect=responses):
            with pytest.raises(StreamStatsResponseError, match="bcrequest"):
                delineate_and_get_characteristics("WA", 48.57430, -120.37890)

    def test_hydro_422_raises_response_error(self) -> None:
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_good_get_responses()),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(status_code=422, text=RESPONSE_422_MISSING_REGION),
            ),
        ):
            with pytest.raises(StreamStatsResponseError):
                delineate_and_get_characteristics("WA", 48.57430, -120.37890)

    def test_hydro_wrapped_parameters_shape_supported(self) -> None:
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_good_get_responses()),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(HYDRO_CHARACTERISTICS_WRAPPED),
            ),
        ):
            result = delineate_and_get_characteristics("WA", 48.57430, -120.37890)

        assert result["DRNAREA"].value == pytest.approx(412.0)

    def test_malformed_hydro_response_raises(self) -> None:
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_good_get_responses()),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(HYDRO_CHARACTERISTICS_MALFORMED),
            ),
        ):
            with pytest.raises(StreamStatsResponseError):
                delineate_and_get_characteristics("WA", 48.57430, -120.37890)

    def test_characteristic_missing_value_raises(self) -> None:
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_good_get_responses()),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(HYDRO_CHARACTERISTICS_MISSING_VALUE),
            ),
        ):
            with pytest.raises(StreamStatsResponseError):
                delineate_and_get_characteristics("WA", 48.57430, -120.37890)

    def test_bc_labels_defaults_to_wildcard(self) -> None:
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_good_get_responses()),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(HYDRO_CHARACTERISTICS_GOOD),
            ) as mock_post,
        ):
            delineate_and_get_characteristics("WA", 48.57430, -120.37890)

        assert mock_post.call_args.kwargs["params"]["bcLabels"] == "*"

    def test_bc_labels_uses_requested_codes(self) -> None:
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_good_get_responses()),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(HYDRO_CHARACTERISTICS_GOOD),
            ) as mock_post,
        ):
            delineate_and_get_characteristics(
                "WA", 48.57430, -120.37890, characteristic_codes=["DRNAREA", "PRECPRIS10"]
            )

        assert mock_post.call_args.kwargs["params"]["bcLabels"] == "DRNAREA,PRECPRIS10"

    def test_region_is_never_inferred(self) -> None:
        """FR-8: the region passed through is exactly the one the caller supplied."""
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_good_get_responses()),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(HYDRO_CHARACTERISTICS_GOOD),
            ),
        ):
            result = delineate_and_get_characteristics("WA", 48.57430, -120.37890)

        assert result.region == "WA"

    def test_sshydro_call_uses_lat_lon_params(self) -> None:
        """The sshydro delineate call takes lat/lon, matching the rest of the API and
        the design doc's own literally-verified protocol.
        """
        with (
            patch(
                "flowfreq.streamstats.requests.get", side_effect=_good_get_responses()
            ) as mock_get,
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(HYDRO_CHARACTERISTICS_GOOD),
            ),
        ):
            delineate_and_get_characteristics("WA", 48.57430, -120.37890)

        sshydro_call = mock_get.call_args_list[1]
        params = sshydro_call.kwargs["params"]
        assert set(params) == {"lat", "lon"}

    def test_provenance_completeness(self) -> None:
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_good_get_responses()),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(HYDRO_CHARACTERISTICS_GOOD),
            ),
        ):
            result = delineate_and_get_characteristics("WA", 48.57430, -120.37890)

        assert result.provenance is not None
        assert result.provenance.service_versions["ss-delineate"]
        assert result.provenance.service_versions["ss-hydro"]
        assert len(result.provenance.request_urls) == 2
        assert result.provenance.requested_at_utc
        assert result.provenance.server_used


def _run_delineation(snap, sshydro, hydro, **kwargs):
    """Run the pipeline against one snap/sshydro/ss-hydro triple; returns (result, post)."""
    with (
        patch(
            "flowfreq.streamstats.requests.get",
            side_effect=[_mock_response(snap), _mock_response(sshydro)],
        ),
        patch("flowfreq.streamstats.requests.post", return_value=_mock_response(hydro)) as post,
    ):
        lon, lat = snap["input"]["coordinates"]
        return delineate_and_get_characteristics("WA", lat, lon, **kwargs), post


def _outer_ring(sshydro):
    for entry in sshydro["bcrequest"]["wsresp"]["featurecollection"][0]:
        if entry["name"] == "globalwatershed":
            return entry["feature"]["features"][0]["geometry"]["coordinates"][0]
    raise AssertionError("no globalwatershed in capture")


class TestWatershedPolygon:
    """Design doc S10: the globalwatershed polygon, returned only once validated."""

    def test_goat_creek_polygon_returned_and_validated(self) -> None:
        result, _ = _run_delineation(SNAP_GOOD, DELINEATE_SSHYDRO_GOOD, HYDRO_CHARACTERISTICS_GOOD)

        feature = result.polygon_geojson
        assert feature is not None
        assert feature["type"] == "Feature"
        assert feature["geometry"]["type"] == "Polygon"
        ring = feature["geometry"]["coordinates"][0]
        assert ring[0] == ring[-1]
        assert feature["properties"]["GlobalWshd"] == 1
        # Trimmed capture: 411.72 mi^2 against ss-hydro's DRNAREA of 412.0.
        assert result.polygon_area_sq_mi == pytest.approx(411.72, abs=0.01)
        assert result.polygon_area_sq_mi == pytest.approx(412.0, rel=POLYGON_AREA_RTOL)

    def test_hole_is_subtracted_from_area(self) -> None:
        result, _ = _run_delineation(
            SNAP_METHOW_PATEROS,
            DELINEATE_SSHYDRO_METHOW_PATEROS,
            HYDRO_CHARACTERISTICS_METHOW_PATEROS,
        )

        rings = result.polygon_geojson["geometry"]["coordinates"]
        assert len(rings) == 2  # the live capture's exterior ring and one hole
        exterior_only = geojson_area_sq_mi({"type": "Polygon", "coordinates": [rings[0]]})
        assert result.polygon_area_sq_mi < exterior_only
        assert result.polygon_area_sq_mi == pytest.approx(1792.80, abs=0.01)

    def test_include_polygon_false_validates_but_drops_geometry(self) -> None:
        result, _ = _run_delineation(
            SNAP_GOOD,
            DELINEATE_SSHYDRO_GOOD,
            HYDRO_CHARACTERISTICS_GOOD,
            include_polygon=False,
        )

        assert result.polygon_geojson is None
        assert result.polygon_area_sq_mi == pytest.approx(411.72, abs=0.01)

    def test_include_polygon_false_still_rejects_a_bad_polygon(self) -> None:
        with pytest.raises(DegenerateDelineationError, match="DRNAREA"):
            _run_delineation(
                SNAP_GOOD,
                DELINEATE_SSHYDRO_GOOD,
                HYDRO_CHARACTERISTICS_METHOW_PATEROS,
                include_polygon=False,
            )

    def test_area_disagreeing_with_drnarea_raises(self) -> None:
        """Goat Creek's 412 mi^2 polygon against the Methow's 1793 mi^2 DRNAREA."""
        with pytest.raises(DegenerateDelineationError, match="disagrees with DRNAREA"):
            _run_delineation(
                SNAP_GOOD, DELINEATE_SSHYDRO_GOOD, HYDRO_CHARACTERISTICS_METHOW_PATEROS
            )

    def test_drnarea_in_square_kilometres_is_converted(self) -> None:
        hydro = [dict(HYDRO_CHARACTERISTICS_GOOD[0], value=1067.08, unit="square kilometers")]
        result, _ = _run_delineation(SNAP_GOOD, DELINEATE_SSHYDRO_GOOD, hydro)

        assert result.polygon_area_sq_mi == pytest.approx(412.0, rel=POLYGON_AREA_RTOL)

    def test_unrecognised_drnarea_unit_raises(self) -> None:
        hydro = [dict(HYDRO_CHARACTERISTICS_GOOD[0], unit="acres")]
        with pytest.raises(StreamStatsResponseError, match="unrecognised unit"):
            _run_delineation(SNAP_GOOD, DELINEATE_SSHYDRO_GOOD, hydro)

    def test_missing_drnarea_skips_only_the_area_check(self) -> None:
        hydro = [c for c in HYDRO_CHARACTERISTICS_GOOD if c["code"] != "DRNAREA"]
        result, _ = _run_delineation(SNAP_GOOD, DELINEATE_SSHYDRO_GOOD, hydro)

        assert result.polygon_geojson is not None
        assert "DRNAREA" not in result

    def test_snapped_point_outside_polygon_raises_before_ss_hydro(self) -> None:
        """The Methow's pour point is ~60 km from Goat Creek's basin."""
        with pytest.raises(DegenerateDelineationError, match="outside the delineated"):
            _run_delineation(
                SNAP_METHOW_PATEROS, DELINEATE_SSHYDRO_GOOD, HYDRO_CHARACTERISTICS_GOOD
            )

    @staticmethod
    def _goat_shifted_west(metres: float):
        sshydro = load_capture("sshydro_WA_goat_creek_trimmed.json")
        lon, lat = SNAP_GOOD["output"]["coordinates"]
        ring = _outer_ring(sshydro)
        shift = metres / (111320.0 * math.cos(math.radians(lat)))
        ring[:] = [[x - shift, y] for x, y in ring]
        geometry = {"type": "Polygon", "coordinates": [ring]}
        return sshydro, _point_in_geometry(lon, lat, geometry)

    def test_snapped_point_just_outside_boundary_is_accepted(self) -> None:
        """The outlet sits on the boundary by construction; a few metres out is fine."""
        sshydro, inside = self._goat_shifted_west(40.0)  # pour point ends ~24 m outside
        assert not inside

        result, _ = _run_delineation(SNAP_GOOD, sshydro, HYDRO_CHARACTERISTICS_GOOD)

        assert result.polygon_geojson is not None

    def test_snapped_point_beyond_boundary_tolerance_raises(self) -> None:
        sshydro, inside = self._goat_shifted_west(250.0)  # pour point ends ~230 m outside
        assert not inside

        with pytest.raises(DegenerateDelineationError, match="outside the delineated"):
            _run_delineation(SNAP_GOOD, sshydro, HYDRO_CHARACTERISTICS_GOOD)

    def test_unclosed_ring_raises_before_ss_hydro(self) -> None:
        sshydro = load_capture("sshydro_WA_goat_creek_trimmed.json")
        _outer_ring(sshydro).pop()
        with (
            patch(
                "flowfreq.streamstats.requests.get",
                side_effect=[_mock_response(SNAP_GOOD), _mock_response(sshydro)],
            ),
            patch("flowfreq.streamstats.requests.post") as post,
        ):
            with pytest.raises(DegenerateDelineationError, match="not closed"):
                delineate_and_get_characteristics("WA", 48.57426, -120.37893)
        post.assert_not_called()

    def test_ring_with_too_few_positions_raises(self) -> None:
        sshydro = load_capture("sshydro_WA_goat_creek_trimmed.json")
        ring = _outer_ring(sshydro)
        ring[:] = [ring[0], ring[1], ring[0]]
        with pytest.raises(DegenerateDelineationError, match="fewer than 4"):
            _run_delineation(SNAP_GOOD, sshydro, HYDRO_CHARACTERISTICS_GOOD)

    def test_missing_globalwatershed_raises(self) -> None:
        sshydro = load_capture("sshydro_WA_goat_creek_trimmed.json")
        entries = sshydro["bcrequest"]["wsresp"]["featurecollection"][0]
        entries[:] = [e for e in entries if e["name"] != "globalwatershed"]
        with pytest.raises(DegenerateDelineationError, match="no globalwatershed"):
            _run_delineation(SNAP_GOOD, sshydro, HYDRO_CHARACTERISTICS_GOOD)

    def test_point_geometry_in_place_of_polygon_raises(self) -> None:
        """What `delineate/features` returned live: a zero-area Point, not a basin."""
        sshydro = load_capture("sshydro_WA_goat_creek_trimmed.json")
        for entry in sshydro["bcrequest"]["wsresp"]["featurecollection"][0]:
            if entry["name"] == "globalwatershed":
                entry["feature"]["features"][0]["geometry"] = {
                    "type": "Point",
                    "coordinates": SNAP_GOOD["output"]["coordinates"],
                }
        with pytest.raises(DegenerateDelineationError, match="not a Polygon"):
            _run_delineation(SNAP_GOOD, sshydro, HYDRO_CHARACTERISTICS_GOOD)

    def test_missing_wsresp_is_a_response_shape_error(self) -> None:
        sshydro = load_capture("sshydro_WA_goat_creek_trimmed.json")
        del sshydro["bcrequest"]["wsresp"]
        with pytest.raises(StreamStatsResponseError, match="featurecollection"):
            _run_delineation(SNAP_GOOD, sshydro, HYDRO_CHARACTERISTICS_GOOD)

    def test_cache_round_trip_keeps_polygon(self, tmp_path) -> None:
        cache = StreamStatsCache(tmp_path / "cache.json")
        first, _ = _run_delineation(
            SNAP_GOOD, DELINEATE_SSHYDRO_GOOD, HYDRO_CHARACTERISTICS_GOOD, cache=cache
        )
        with patch(
            "flowfreq.streamstats.requests.get",
            side_effect=AssertionError("cache hit must not touch the network"),
        ):
            second = delineate_and_get_characteristics(
                "WA", 48.57426, -120.37893, cache=StreamStatsCache(tmp_path / "cache.json")
            )

        assert second.polygon_geojson == first.polygon_geojson
        assert second.polygon_area_sq_mi == first.polygon_area_sq_mi


class TestGeojsonArea:
    def test_lat_lon_rectangle_matches_spherical_zone_formula(self) -> None:
        """A lon/lat rectangle's exact spherical area is R^2 * dlon * (sin phi2 - sin phi1)."""
        lon1, lon2, lat1, lat2 = -120.0, -119.5, 48.0, 48.5
        ring = [[lon1, lat1], [lon2, lat1], [lon2, lat2], [lon1, lat2], [lon1, lat1]]
        r = 6378137.0
        exact = (
            r**2
            * math.radians(lon2 - lon1)
            * (math.sin(math.radians(lat2)) - math.sin(math.radians(lat1)))
            / SQ_M_PER_SQ_MI
        )

        area = geojson_area_sq_mi({"type": "Polygon", "coordinates": [ring]})

        assert area == pytest.approx(exact, rel=1e-9)

    def test_ring_orientation_does_not_matter(self) -> None:
        ring = [[0.0, 0.0], [0.1, 0.0], [0.1, 0.1], [0.0, 0.1], [0.0, 0.0]]
        cw = geojson_area_sq_mi({"type": "Polygon", "coordinates": [ring]})
        ccw = geojson_area_sq_mi({"type": "Polygon", "coordinates": [ring[::-1]]})
        assert cw == pytest.approx(ccw)

    def test_multipolygon_sums_its_parts(self) -> None:
        a = [[0.0, 0.0], [0.1, 0.0], [0.1, 0.1], [0.0, 0.1], [0.0, 0.0]]
        b = [[1.0, 0.0], [1.1, 0.0], [1.1, 0.1], [1.0, 0.1], [1.0, 0.0]]
        one = geojson_area_sq_mi({"type": "Polygon", "coordinates": [a]})
        both = geojson_area_sq_mi({"type": "MultiPolygon", "coordinates": [[a], [b]]})
        assert both == pytest.approx(2 * one, rel=1e-6)

    def test_non_polygon_raises(self) -> None:
        with pytest.raises(DegenerateDelineationError):
            geojson_area_sq_mi({"type": "Point", "coordinates": [0.0, 0.0]})


class TestStreamStatsCache:
    """Tests for NFR-1/NFR-5: cache round-trip, and offline reuse."""

    def test_round_trip_and_offline_reuse(self, tmp_path) -> None:
        cache_path = tmp_path / "cache.json"
        cache = StreamStatsCache(cache_path)

        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_good_get_responses()),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(HYDRO_CHARACTERISTICS_GOOD),
            ),
        ):
            first = delineate_and_get_characteristics("WA", 48.57430, -120.37890, cache=cache)

        assert cache_path.exists()

        reloaded_cache = StreamStatsCache(cache_path)
        with patch(
            "flowfreq.streamstats.requests.get",
            side_effect=AssertionError("cache hit must not touch the network"),
        ):
            second = delineate_and_get_characteristics(
                "WA", 48.57430, -120.37890, cache=reloaded_cache
            )

        assert second["DRNAREA"].value == first["DRNAREA"].value
        assert second.provenance.server_used == first.provenance.server_used

    def test_missing_cache_file_starts_empty(self, tmp_path) -> None:
        cache = StreamStatsCache(tmp_path / "does_not_exist.json")
        assert cache.get("anything") is None

    def test_corrupt_cache_file_starts_empty_rather_than_raising(self, tmp_path) -> None:
        path = tmp_path / "corrupt.json"
        path.write_text("not json", encoding="utf-8")
        cache = StreamStatsCache(path)
        assert cache.get("anything") is None


class TestWatershedCharacteristicsSerialization:
    def test_to_dict_from_dict_round_trip(self) -> None:
        snap = SnapResult(48.5, -120.3, 48.50001, -120.30001, True, 1.2)
        char = Characteristic("DRNAREA", "Drainage Area", "desc", 412.0, "square miles", "")
        prov = Provenance(
            {"ss-delineate": "1.2.0"}, ["url1"], "2026-09-10T00:00:00+00:00", "prodweba"
        )
        wc = WatershedCharacteristics("WA", snap, None, {"DRNAREA": char}, prov)

        restored = WatershedCharacteristics.from_dict(wc.to_dict())

        assert restored.region == wc.region
        assert restored["DRNAREA"].value == pytest.approx(412.0)
        assert restored.provenance.server_used == "prodweba"


class TestRegionFlowEstimatesSerialization:
    def test_to_dict_from_dict_round_trip(self) -> None:
        estimate = FlowStatisticEstimate(
            code="PK50AEP",
            name="50-percent AEP flood",
            description="",
            value=4370.0,
            unit="ft^3/s",
            equation="3.846*DRNAREA^0.745",
            standard_error_pct=95.0,
            interval_lower=1140.0,
            interval_upper=16700.0,
        )
        citation = RegressionCitation(
            citation_id=150,
            title="2016, Magnitude, frequency, and trends of floods...",
            author="Mastin, M.C., et al.",
            citation_url="http://dx.doi.org/10.3133/sir20165118",
            last_year_of_data=2014,
        )
        region = RegionFlowEstimates(
            region_code="GC1750",
            region_name="Peak_Region_1_2016_5118",
            statistic_group_code="PFS",
            statistic_group_name="Peak-Flow Statistics",
            estimates={"PK50AEP": estimate},
            citation=citation,
        )

        restored = RegionFlowEstimates.from_dict(region.to_dict())

        assert restored.region_code == "GC1750"
        assert restored["PK50AEP"].value == pytest.approx(4370.0)
        assert restored.citation is not None
        assert restored.citation.citation_id == 150
        assert restored.located is None and restored.percent_weight is None

    def test_location_labels_round_trip(self) -> None:
        region = RegionFlowEstimates(
            "areaave", "Area-Averaged", "PFS", "Peak-Flow Statistics", located=True
        )
        region.area_averaged = True
        labelled = RegionFlowEstimates("GC1751", "", "PFS", "", percent_weight=100.0)

        assert RegionFlowEstimates.from_dict(region.to_dict()).area_averaged is True
        assert RegionFlowEstimates.from_dict(labelled.to_dict()).percent_weight == 100.0


class TestListRegions:
    def test_parses_region_list(self) -> None:
        payload = [
            {"id": 53, "name": "Washington", "code": "WA"},
            {"id": 7, "name": "California", "code": "CA"},
        ]
        with patch("flowfreq.streamstats.requests.get", return_value=_mock_response(payload)):
            regions = list_regions()

        assert {"WA", "CA"} <= {r["code"] for r in regions}

    def test_non_list_response_raises(self) -> None:
        with patch(
            "flowfreq.streamstats.requests.get",
            return_value=_mock_response({"not": "a list"}),
        ):
            with pytest.raises(StreamStatsResponseError):
                list_regions()


class TestBatchGetCharacteristics:
    """Tests for FR-7: one bad point never aborts the batch."""

    def test_partial_failure_does_not_abort_batch(self) -> None:
        get_responses = _good_get_responses("PRODWEBA") + [_mock_response(SNAP_UNSNAPPABLE)]
        points = [
            ("node_good", "WA", 48.57430, -120.37890),
            ("node_bad", "WA", 48.584, -120.370),
        ]
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=get_responses),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(HYDRO_CHARACTERISTICS_GOOD),
            ),
        ):
            results, errors = batch_get_characteristics(
                points, concurrency=1, validate_regions=False
            )

        assert "node_good" in results
        assert results["node_good"]["DRNAREA"].value == pytest.approx(412.0)
        assert "node_bad" in errors
        assert "snap" in errors["node_bad"].lower()

    def test_empty_batch_returns_empty_results(self) -> None:
        results, errors = batch_get_characteristics([], validate_regions=False)
        assert results == {}
        assert errors == {}

    def test_concurrency_is_capped_at_documented_rate_limit(self) -> None:
        with patch(
            "flowfreq.streamstats.ThreadPoolExecutor", wraps=ThreadPoolExecutor
        ) as mock_pool_cls:
            batch_get_characteristics([], concurrency=99, validate_regions=False)

        assert mock_pool_cls.call_args.kwargs["max_workers"] == MAX_CONCURRENCY

    def test_concurrency_defaults_to_serial(self) -> None:
        with patch(
            "flowfreq.streamstats.ThreadPoolExecutor", wraps=ThreadPoolExecutor
        ) as mock_pool_cls:
            batch_get_characteristics([], validate_regions=False)

        assert mock_pool_cls.call_args.kwargs["max_workers"] == 1

    def test_bad_region_raises_before_any_network_call(self) -> None:
        with (
            patch(
                "flowfreq.streamstats.list_regions",
                return_value=[{"id": 53, "name": "Washington", "code": "WA"}],
            ),
            patch("flowfreq.streamstats.requests.get") as mock_get,
        ):
            with pytest.raises(UnsupportedRegionError):
                batch_get_characteristics([("n1", "ZZ", 48.5, -120.3)], validate_regions=True)

        mock_get.assert_not_called()


_GOOD_CHARACTERISTICS = {
    "DRNAREA": Characteristic("DRNAREA", "Drainage Area", "", 412.0, "square miles"),
    "PRECPRIS10": Characteristic("PRECPRIS10", "Mean Annual Precip", "", 45.62, "inches"),
    "CANOPY_PCT": Characteristic("CANOPY_PCT", "Percent Canopy", "", 45.242, "percent"),
}


def _nss_get_dispatch(url, **kwargs):
    if url.endswith("/nssservices/statisticgroups"):
        return _mock_response(NSS_STATISTIC_GROUPS_ALL)
    if url.endswith("/nssservices/regions/WA/statisticgroups"):
        return _mock_response(NSS_STATISTIC_GROUPS_WA)
    if url.endswith("/nssservices/regions/WA/Scenarios"):
        return _mock_response(NSS_SCENARIO_TEMPLATE_PFS)
    if url.endswith("/nssservices/citations"):
        return _mock_response(NSS_CITATIONS_RESPONSE)
    raise AssertionError(f"unexpected GET {url}")


class TestFillAndValidateRegions:
    """Tests for the client-side range check NSS itself will not perform."""

    def test_valid_values_are_filled_and_region_kept(self) -> None:
        import copy

        scenario = copy.deepcopy(NSS_SCENARIO_TEMPLATE_PFS[0])
        filled, skipped = _fill_and_validate_regions(scenario, _GOOD_CHARACTERISTICS)

        assert skipped == {}
        assert len(filled["regressionRegions"]) == 1
        values = {p["code"]: p["value"] for p in filled["regressionRegions"][0]["parameters"]}
        assert values == {"DRNAREA": 412.0, "PRECPRIS10": 45.62, "CANOPY_PCT": 45.242}

    def test_missing_characteristic_skips_region(self) -> None:
        import copy

        scenario = copy.deepcopy(NSS_SCENARIO_TEMPLATE_PFS[0])
        incomplete = {k: v for k, v in _GOOD_CHARACTERISTICS.items() if k != "CANOPY_PCT"}

        filled, skipped = _fill_and_validate_regions(scenario, incomplete)

        assert filled["regressionRegions"] == []
        assert "GC1750" in skipped
        assert "CANOPY_PCT" in skipped["GC1750"]

    def test_out_of_range_value_skips_region(self) -> None:
        """The mandatory test for this module's own governing principle, restated
        for NSS: an out-of-range input must never be submitted as if valid.
        """
        import copy

        scenario = copy.deepcopy(NSS_SCENARIO_TEMPLATE_PFS[0])
        out_of_range = dict(_GOOD_CHARACTERISTICS)
        out_of_range["DRNAREA"] = Characteristic("DRNAREA", "", "", 999999.0, "square miles")

        filled, skipped = _fill_and_validate_regions(scenario, out_of_range)

        assert filled["regressionRegions"] == []
        assert "outside this region's valid range" in skipped["GC1750"]


class TestListStatisticGroupsNSS:
    def test_no_region_lists_all_groups(self) -> None:
        with patch(
            "flowfreq.streamstats.requests.get",
            return_value=_mock_response(NSS_STATISTIC_GROUPS_ALL),
        ) as mock_get:
            groups = list_statistic_groups()

        assert mock_get.call_args.args[0].endswith("/nssservices/statisticgroups")
        assert {g["code"] for g in groups} == {"PC", "PFS", "LFS"}

    def test_region_lists_only_valid_groups(self) -> None:
        with patch(
            "flowfreq.streamstats.requests.get",
            return_value=_mock_response(NSS_STATISTIC_GROUPS_WA),
        ) as mock_get:
            groups = list_statistic_groups("WA")

        assert mock_get.call_args.args[0].endswith("/nssservices/regions/WA/statisticgroups")
        assert {g["code"] for g in groups} == {"PFS", "LFS"}

    def test_non_list_response_raises(self) -> None:
        with patch(
            "flowfreq.streamstats.requests.get", return_value=_mock_response({"not": "a list"})
        ):
            with pytest.raises(StreamStatsResponseError):
                list_statistic_groups("WA")


class TestEstimateFlowStatistics:
    """Tests for the Phase 2 NSS pipeline (docs/STREAMSTATS_NSS_ADDENDUM.md)."""

    def test_happy_path_returns_estimates_with_citation(self) -> None:
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_nss_get_dispatch),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(NSS_ESTIMATE_RESPONSE_PFS),
            ),
        ):
            region_estimates, skipped = estimate_flow_statistics(
                "WA", _GOOD_CHARACTERISTICS, statistic_group_codes=["PFS"]
            )

        assert skipped == {}
        assert len(region_estimates) == 1
        result = region_estimates[0]
        assert result.region_code == "GC1750"
        assert result.statistic_group_code == "PFS"
        assert result["PK50AEP"].value == pytest.approx(4370.0)
        assert result["PK50AEP"].equation.startswith("3.846*DRNAREA")
        assert result["PK50AEP"].standard_error_pct == pytest.approx(95.0)
        assert result.citation is not None
        assert "Mastin" in result.citation.author
        assert result.citation.citation_url == "http://dx.doi.org/10.3133/sir20165118"

    def test_estimate_body_is_a_bare_array_not_wrapped(self) -> None:
        """The mandatory regression test: wrapping the body in {"scenarioList": [...]}
        (the api-config's own literal parameter name) is a confirmed-live 500 with no
        other diagnostic. Posting the array bare is what actually works.
        """
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_nss_get_dispatch),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(NSS_ESTIMATE_RESPONSE_PFS),
            ) as mock_post,
        ):
            estimate_flow_statistics("WA", _GOOD_CHARACTERISTICS, statistic_group_codes=["PFS"])

        body = mock_post.call_args.kwargs["json"]
        assert isinstance(body, list)
        assert not isinstance(body, dict)

    def test_default_statistic_groups_discovered_live(self) -> None:
        """With no explicit groups, every group valid for the region is requested."""
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_nss_get_dispatch) as mock_get,
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(NSS_ESTIMATE_RESPONSE_PFS),
            ),
        ):
            estimate_flow_statistics("WA", _GOOD_CHARACTERISTICS)

        template_call = next(c for c in mock_get.call_args_list if c.args[0].endswith("/Scenarios"))
        assert template_call.kwargs["params"]["statisticgroups"] == "PFS,LFS"

    def test_all_regions_invalid_skips_network_call_to_estimate(self) -> None:
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_nss_get_dispatch),
            patch("flowfreq.streamstats.requests.post") as mock_post,
        ):
            region_estimates, skipped = estimate_flow_statistics(
                "WA", {}, statistic_group_codes=["PFS"]
            )

        assert region_estimates == []
        assert "PFS:GC1750" in skipped
        mock_post.assert_not_called()

    def test_estimate_4xx_raises_response_error(self) -> None:
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_nss_get_dispatch),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(status_code=422, text="bad request"),
            ),
        ):
            with pytest.raises(StreamStatsResponseError):
                estimate_flow_statistics("WA", _GOOD_CHARACTERISTICS, statistic_group_codes=["PFS"])

    def test_estimate_persistent_500_raises_transport_error(self) -> None:
        """The confirmed-live real bug (wrapping the body in {"scenarioList": [...]})
        surfaces as a 500 that persists across retries -- exercised here via the
        actual captured error body, expecting the transport-failure path rather than
        the 4xx-data-problem path.
        """
        with (
            patch("flowfreq.streamstats.requests.get", side_effect=_nss_get_dispatch),
            patch("flowfreq.streamstats.time.sleep", return_value=None),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(status_code=500, text=NSS_ESTIMATE_500_BODY),
            ),
        ):
            with pytest.raises(StreamStatsTransportError):
                estimate_flow_statistics("WA", _GOOD_CHARACTERISTICS, statistic_group_codes=["PFS"])

    def test_malformed_template_response_raises(self) -> None:
        def dispatch(url, **kwargs):
            if url.endswith("/Scenarios"):
                return _mock_response({"not": "a list"})
            return _nss_get_dispatch(url, **kwargs)

        with patch("flowfreq.streamstats.requests.get", side_effect=dispatch):
            with pytest.raises(StreamStatsResponseError):
                estimate_flow_statistics("WA", _GOOD_CHARACTERISTICS, statistic_group_codes=["PFS"])

    def test_citation_resolution_failure_does_not_fail_the_whole_call(self) -> None:
        """A citations-lookup failure is not fatal -- the estimates themselves are
        still good; only the citation attachment is best-effort.
        """

        def dispatch(url, **kwargs):
            if url.endswith("/nssservices/citations"):
                return _mock_response(status_code=500, text="error")
            return _nss_get_dispatch(url, **kwargs)

        with (
            patch("flowfreq.streamstats.requests.get", side_effect=dispatch),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(NSS_ESTIMATE_RESPONSE_PFS),
            ),
        ):
            region_estimates, _ = estimate_flow_statistics(
                "WA", _GOOD_CHARACTERISTICS, statistic_group_codes=["PFS"]
            )

        assert region_estimates[0]["PK50AEP"].value == pytest.approx(4370.0)
        assert region_estimates[0].citation is None


def _global_watershed(sshydro):
    for entry in sshydro["bcrequest"]["wsresp"]["featurecollection"][0]:
        if entry["name"] == "globalwatershed":
            return entry["feature"]["features"][0]
    raise AssertionError("no globalwatershed in capture")


GOAT_POLYGON = _global_watershed(DELINEATE_SSHYDRO_GOOD)
OGEECHEE_POLYGON = _global_watershed(DELINEATE_SSHYDRO_GA_OGEECHEE)

GOAT_CHARACTERISTICS = {
    "DRNAREA": Characteristic("DRNAREA", "", "", 412.0, "square miles"),
    "PRECPRIS10": Characteristic("PRECPRIS10", "", "", 45.62, "inches"),
    "CANOPY_PCT": Characteristic("CANOPY_PCT", "", "", 45.242, "percent"),
}


def _wa_pfs_template():
    return [s for s in copy.deepcopy(NSS_SCENARIOS_WA_PFS_LFS) if s["statisticGroupID"] == 2]


class TestLocateRegressionRegions:
    """Addendum S5: which NSS regions a watershed polygon falls in."""

    def test_goat_creek_is_in_peak_region_2(self) -> None:
        with patch(
            "flowfreq.streamstats.requests.post",
            return_value=_mock_response(NSS_BYLOCATION_WA_GOAT_CREEK),
        ) as post:
            located = locate_regression_regions("WA", GOAT_POLYGON)

        assert [(r.code, r.percent_weight) for r in located] == [("GC1751", 100.0)]
        assert located[0].region_id == 718
        assert located[0].citation_id == 150
        url = post.call_args.args[0]
        assert url.endswith("/nssservices/regions/WA/regressionregions/bylocation")
        # The body is the bare geometry, not the Feature (as the USGS notebook posts it).
        assert post.call_args.kwargs["json"] == GOAT_POLYGON["geometry"]

    def test_basin_spanning_regions_reports_area_weights(self) -> None:
        with patch(
            "flowfreq.streamstats.requests.post",
            return_value=_mock_response(NSS_BYLOCATION_GA_OGEECHEE),
        ):
            located = locate_regression_regions("GA", OGEECHEE_POLYGON["geometry"])

        weights = {r.code: r.percent_weight for r in located}
        assert weights["GC1541"] == 32.0  # SIR 2014-5030 region 3
        assert weights["GC1542"] == 24.0  # region 4
        assert weights["GC1572"] == 44.0  # region 1
        assert weights["GC1934"] == 100.0  # the 2023 Southeast US method, whole basin
        assert sum(weights[c] for c in ("GC1539", "GC1541", "GC1542")) == 100.0

    def test_answer_about_another_polygon_is_refused(self) -> None:
        """Goat Creek's 412 mi^2 at 100% cannot describe the 807 mi^2 Ogeechee basin."""
        with patch(
            "flowfreq.streamstats.requests.post",
            return_value=_mock_response(NSS_BYLOCATION_WA_GOAT_CREEK),
        ):
            with pytest.raises(StreamStatsResponseError, match="percentWeight"):
                locate_regression_regions("GA", OGEECHEE_POLYGON)

    def test_overlap_larger_than_basin_is_refused(self) -> None:
        answer = copy.deepcopy(NSS_BYLOCATION_WA_GOAT_CREEK)
        answer[0]["area"] = 900.0
        with patch("flowfreq.streamstats.requests.post", return_value=_mock_response(answer)):
            with pytest.raises(StreamStatsResponseError, match="not about this polygon"):
                locate_regression_regions("WA", GOAT_POLYGON)

    def test_entry_missing_percent_weight_is_refused(self) -> None:
        answer = copy.deepcopy(NSS_BYLOCATION_WA_GOAT_CREEK)
        del answer[0]["percentWeight"]
        with patch("flowfreq.streamstats.requests.post", return_value=_mock_response(answer)):
            with pytest.raises(StreamStatsResponseError, match="missing"):
                locate_regression_regions("WA", GOAT_POLYGON)

    def test_rejected_request_raises(self) -> None:
        body = (
            '{"code":400,"message":"Geometry is not of type: Polygon,MultiPolygon",'
            '"content":"Bad Request Received"}'
        )
        with patch(
            "flowfreq.streamstats.requests.post",
            return_value=_mock_response(status_code=400, text=body),
        ):
            with pytest.raises(StreamStatsResponseError, match="Polygon,MultiPolygon"):
                locate_regression_regions("WA", GOAT_POLYGON)

    def test_point_geometry_is_refused_before_any_request(self) -> None:
        point = {"type": "Point", "coordinates": [-120.37893, 48.57426]}
        with patch("flowfreq.streamstats.requests.post") as post:
            with pytest.raises(DegenerateDelineationError):
                locate_regression_regions("WA", point)
        post.assert_not_called()

    def test_empty_answer_means_no_region(self) -> None:
        with patch("flowfreq.streamstats.requests.post", return_value=_mock_response([])):
            assert locate_regression_regions("WA", GOAT_POLYGON) == []


def _nss_dispatch_for(template):
    def dispatch(url, **kwargs):
        if url.endswith("/Scenarios"):
            return _mock_response(copy.deepcopy(template))
        return _nss_get_dispatch(url, **kwargs)

    return dispatch


def _post_router(bylocation, estimate):
    """Route bylocation and Scenarios/Estimate POSTs; records the estimate bodies."""
    bodies = []

    def post(url, **kwargs):
        if url.endswith("/bylocation"):
            return _mock_response(bylocation)
        if url.endswith("/Scenarios/Estimate"):
            bodies.append(copy.deepcopy(kwargs["json"]))
            return _mock_response(estimate)
        raise AssertionError(f"unexpected POST {url}")

    return post, bodies


class TestEstimateWithLocation:
    """estimate_flow_statistics with a watershed polygon (addendum S5)."""

    def test_only_the_located_region_is_estimated(self) -> None:
        post, bodies = _post_router(NSS_BYLOCATION_WA_GOAT_CREEK, NSS_ESTIMATE_WA_GOAT_CREEK_GC1751)
        with (
            patch(
                "flowfreq.streamstats.requests.get",
                side_effect=_nss_dispatch_for(_wa_pfs_template()),
            ),
            patch("flowfreq.streamstats.requests.post", side_effect=post),
        ):
            results, skipped = estimate_flow_statistics(
                "WA", GOAT_CHARACTERISTICS, ["PFS"], watershed_polygon=GOAT_POLYGON
            )

        submitted = [r["code"] for s in bodies[0] for r in s["regressionRegions"]]
        assert submitted == ["GC1751"]
        assert [r.region_code for r in results] == ["GC1751"]
        assert results[0].located is True
        assert results[0].percent_weight == 100.0
        assert results[0]["PK50AEP"].value == pytest.approx(3290.0)
        for code in ("GC1750", "GC1752", "GC1753"):
            assert "not in this regression region" in skipped[f"PFS:{code}"]
        # A single region: nothing to average, and no weights sent.
        assert "percentWeight" not in bodies[0][0]["regressionRegions"][0]
        assert "PFS:areaave" not in skipped

    def test_include_unlocated_keeps_every_in_range_region_labelled(self) -> None:
        post, bodies = _post_router(
            NSS_BYLOCATION_WA_GOAT_CREEK, NSS_ESTIMATE_WA_GOAT_CREEK_ALL_PFS
        )
        with (
            patch(
                "flowfreq.streamstats.requests.get",
                side_effect=_nss_dispatch_for(_wa_pfs_template()),
            ),
            patch("flowfreq.streamstats.requests.post", side_effect=post),
        ):
            results, _ = estimate_flow_statistics(
                "WA",
                GOAT_CHARACTERISTICS,
                ["PFS"],
                watershed_polygon=GOAT_POLYGON,
                include_unlocated=True,
            )

        labels = {r.region_code: (r.located, r.percent_weight) for r in results}
        assert labels == {
            "GC1750": (False, None),
            "GC1751": (True, 100.0),
            "GC1752": (False, None),
            "GC1753": (False, None),
        }
        assert all("percentWeight" not in r for s in bodies[0] for r in s["regressionRegions"])

    def test_without_polygon_behaviour_is_unchanged(self) -> None:
        with (
            patch(
                "flowfreq.streamstats.requests.get",
                side_effect=_nss_dispatch_for(_wa_pfs_template()),
            ),
            patch(
                "flowfreq.streamstats.requests.post",
                return_value=_mock_response(NSS_ESTIMATE_WA_GOAT_CREEK_ALL_PFS),
            ) as post,
        ):
            results, _ = estimate_flow_statistics("WA", GOAT_CHARACTERISTICS, ["PFS"])

        assert post.call_count == 1  # no bylocation call
        assert {r.region_code for r in results} == {"GC1750", "GC1751", "GC1752", "GC1753"}
        assert all(r.located is None and r.percent_weight is None for r in results)

    def test_no_located_region_estimates_nothing(self) -> None:
        post, bodies = _post_router([], NSS_ESTIMATE_WA_GOAT_CREEK_ALL_PFS)
        with (
            patch(
                "flowfreq.streamstats.requests.get",
                side_effect=_nss_dispatch_for(_wa_pfs_template()),
            ),
            patch("flowfreq.streamstats.requests.post", side_effect=post),
        ):
            results, skipped = estimate_flow_statistics(
                "WA", GOAT_CHARACTERISTICS, ["PFS"], watershed_polygon=GOAT_POLYGON
            )

        assert results == []
        assert bodies == []
        assert len(skipped) == 4


GA_CHARACTERISTICS = {
    "DRNAREA": Characteristic("DRNAREA", "", "", 0.5, "square miles"),
    "LC06IMP": Characteristic("LC06IMP", "", "", 5.0, "percent"),
    "LC06DEV": Characteristic("LC06DEV", "", "", 20.0, "percent"),
}


def _ga_located(weights):
    """A bylocation answer for the Ogeechee polygon with the given percent weights."""
    area = geojson_area_sq_mi(OGEECHEE_POLYGON["geometry"])
    ids = {"GC1572": 107, "GC1573": 108}
    return [
        {
            "id": ids[code],
            "name": code,
            "code": code,
            "citationID": 33,
            "statusID": 4,
            "percentWeight": w,
            "area": area * w / 100.0,
        }
        for code, w in weights.items()
    ]


class TestAreaAveraging:
    """NSS's own area-weighted average for a basin spanning regions (addendum S5)."""

    def _run(self, weights, estimate=NSS_ESTIMATE_GA_AREA_AVERAGED, characteristics=None):
        post, bodies = _post_router(_ga_located(weights), estimate)
        with (
            patch(
                "flowfreq.streamstats.requests.get",
                side_effect=_nss_dispatch_for(NSS_SCENARIOS_GA_RURAL_UNDER_1),
            ),
            patch("flowfreq.streamstats.requests.post", side_effect=post),
        ):
            results, skipped = estimate_flow_statistics(
                "GA",
                characteristics or GA_CHARACTERISTICS,
                ["PFS"],
                watershed_polygon=OGEECHEE_POLYGON,
            )
        return results, skipped, bodies

    def test_weights_sent_and_average_returned(self) -> None:
        results, skipped, bodies = self._run({"GC1572": 60.0, "GC1573": 40.0})

        sent = {r["code"]: r.get("percentWeight") for r in bodies[0][0]["regressionRegions"]}
        assert sent == {"GC1572": 60.0, "GC1573": 40.0}
        average = [r for r in results if r.area_averaged]
        assert len(average) == 1
        assert average[0].region_code == "areaave"
        assert average[0]["PK50AEP"].value == pytest.approx(0.6 * 116.0 + 0.4 * 33.3)
        weights = {r.region_code: r.percent_weight for r in results if not r.area_averaged}
        assert weights == {"GC1572": 60.0, "GC1573": 40.0}
        assert skipped == {}

    def test_tampered_average_is_refused(self) -> None:
        estimate = copy.deepcopy(NSS_ESTIMATE_GA_AREA_AVERAGED)
        estimate[0]["regressionRegions"][2]["results"][0]["value"] = 116.0
        with pytest.raises(StreamStatsResponseError, match="weighted mean"):
            self._run({"GC1572": 60.0, "GC1573": 40.0}, estimate=estimate)

    def test_missing_average_is_refused(self) -> None:
        estimate = copy.deepcopy(NSS_ESTIMATE_GA_AREA_AVERAGED)
        del estimate[0]["regressionRegions"][2]
        with pytest.raises(StreamStatsResponseError, match="no area average"):
            self._run({"GC1572": 60.0, "GC1573": 40.0}, estimate=estimate)

    def test_weights_not_summing_to_100_are_not_sent(self) -> None:
        """Live: NSS silently computes no average for 60 + 30."""
        estimate = copy.deepcopy(NSS_ESTIMATE_GA_AREA_AVERAGED)
        del estimate[0]["regressionRegions"][2]
        results, skipped, bodies = self._run({"GC1572": 60.0, "GC1573": 30.0}, estimate)

        assert all("percentWeight" not in r for r in bodies[0][0]["regressionRegions"])
        assert "sum to 90" in skipped["PFS:areaave"]
        assert not any(r.area_averaged for r in results)

    def test_out_of_range_located_region_blocks_the_average(self) -> None:
        characteristics = dict(GA_CHARACTERISTICS)
        characteristics["LC06DEV"] = Characteristic("LC06DEV", "", "", 1.0, "percent")
        estimate = copy.deepcopy(NSS_ESTIMATE_GA_AREA_AVERAGED)
        estimate[0]["regressionRegions"] = estimate[0]["regressionRegions"][:1]
        results, skipped, bodies = self._run(
            {"GC1572": 60.0, "GC1573": 40.0}, estimate, characteristics
        )

        assert "outside this region's valid range" in skipped["PFS:GC1573"]
        assert "GC1573" in skipped["PFS:areaave"]
        assert [r.region_code for r in results] == ["GC1572"]
        assert "percentWeight" not in bodies[0][0]["regressionRegions"][0]


class TestBatchEstimateFlowStatistics:
    def test_regions_selected_by_each_points_own_polygon(self) -> None:
        def get(url, **kwargs):
            if url.endswith("/pourpoint/v1/snap/str900"):
                return _mock_response(SNAP_GOOD)
            if "/ss-delineate/" in url:
                return _mock_response(DELINEATE_SSHYDRO_GOOD)
            return _nss_dispatch_for(_wa_pfs_template())(url, **kwargs)

        nss_post, bodies = _post_router(
            NSS_BYLOCATION_WA_GOAT_CREEK, NSS_ESTIMATE_WA_GOAT_CREEK_GC1751
        )

        def post(url, **kwargs):
            if "/ss-hydro/" in url:
                return _mock_response(HYDRO_CHARACTERISTICS_GOOD)
            return nss_post(url, **kwargs)

        with (
            patch("flowfreq.streamstats.requests.get", side_effect=get),
            patch("flowfreq.streamstats.requests.post", side_effect=post),
        ):
            estimates, skipped, errors = batch_estimate_flow_statistics(
                [("goat_creek", "WA", 48.57426, -120.37893)],
                statistic_group_codes=["PFS"],
                validate_regions=False,
            )

        assert errors == {}
        assert [(r.region_code, r.located) for r in estimates["goat_creek"]] == [("GC1751", True)]
        assert "PFS:GC1750" in skipped["goat_creek"]

    def test_delineation_failure_reported_as_error(self) -> None:
        with patch(
            "flowfreq.streamstats.requests.get", return_value=_mock_response(SNAP_UNSNAPPABLE)
        ):
            estimates, skipped, errors = batch_estimate_flow_statistics(
                [("bad", "WA", 48.584, -120.370)], validate_regions=False
            )

        assert estimates == {}
        assert "bad" in errors

    def test_successful_point_returns_estimates(self) -> None:
        def dispatch(url, **kwargs):
            if url.endswith("/pourpoint/v1/snap/str900"):
                return _mock_response(SNAP_GOOD)
            if "/ss-delineate/" in url:
                return _mock_response(DELINEATE_SSHYDRO_GOOD)
            return _nss_get_dispatch(url, **kwargs)

        with (
            patch("flowfreq.streamstats.requests.get", side_effect=dispatch),
            patch(
                "flowfreq.streamstats.requests.post",
                side_effect=[
                    _mock_response(HYDRO_CHARACTERISTICS_GOOD),
                    _mock_response(NSS_ESTIMATE_RESPONSE_PFS),
                ],
            ),
        ):
            estimates, skipped, errors = batch_estimate_flow_statistics(
                [("goat_creek", "WA", 48.57430, -120.37890)],
                statistic_group_codes=["PFS"],
                validate_regions=False,
                select_by_location=False,
            )

        assert errors == {}
        assert "goat_creek" in estimates
        assert estimates["goat_creek"][0]["PK50AEP"].value == pytest.approx(4370.0)


@pytest.mark.requires_network
class TestLiveStreamStats:
    """Live StreamStats checks against the design doc's verification points.

    Deselect with -m 'not requires_network'. If ``test_goat_creek_matches_design_doc_appendix``
    fails on the snapped-coordinate assertions, the field-name guesses in
    ``snap_point``'s output parsing (see the module docstring) are the first thing to
    check and correct against the real response.
    """

    def test_goat_creek_matches_design_doc_appendix(self) -> None:
        result = delineate_and_get_characteristics("WA", 48.57426, -120.37893)

        assert result["DRNAREA"].value == pytest.approx(412.0, rel=0.01)
        assert result["PRECPRIS10"].value == pytest.approx(45.62, rel=0.01)
        assert result["CANOPY_PCT"].value == pytest.approx(45.242, rel=0.01)

    def test_lost_river_matches_design_doc_appendix(self) -> None:
        result = delineate_and_get_characteristics("WA", 48.65041, -120.51172)

        assert result["DRNAREA"].value == pytest.approx(252.0, rel=0.01)
        assert result["PRECPRIS10"].value == pytest.approx(47.99, rel=0.01)
        assert result["CANOPY_PCT"].value == pytest.approx(43.090, rel=0.01)

    def test_offnetwork_point_raises(self) -> None:
        with pytest.raises((UnsnappablePointError, DegenerateDelineationError)):
            delineate_and_get_characteristics("WA", 48.584, -120.370)

    def test_goat_creek_polygon_matches_drnarea(self) -> None:
        """Design doc S10: the live polygon, validated, agrees with DRNAREA."""
        result = delineate_and_get_characteristics(
            "WA", 48.57426, -120.37893, characteristic_codes=["DRNAREA"]
        )

        assert result.polygon_geojson is not None
        assert result.polygon_geojson["geometry"]["type"] in ("Polygon", "MultiPolygon")
        # 411.75 mi^2 live on 2026-09-27 against DRNAREA 412.
        assert result.polygon_area_sq_mi == pytest.approx(result["DRNAREA"].value, rel=0.005)

    def test_methow_at_pateros_polygon_with_hole(self) -> None:
        result = delineate_and_get_characteristics(
            "WA", 48.0776, -119.9837, characteristic_codes=["DRNAREA"]
        )

        assert result.polygon_geojson is not None
        assert len(result.polygon_geojson["geometry"]["coordinates"]) >= 2
        # 1791.81 mi^2 live on 2026-09-27 against DRNAREA 1793.
        assert result.polygon_area_sq_mi == pytest.approx(result["DRNAREA"].value, rel=0.005)


@pytest.mark.requires_network
class TestLiveNSS:
    """Live NSS checks (docs/STREAMSTATS_NSS_ADDENDUM.md). Deselect with
    -m 'not requires_network'.
    """

    def test_goat_creek_peak_flow_estimate_is_plausible(self) -> None:
        """Without a polygon every in-range region is returned; GC1750 is one of them
        but is not where Goat Creek is (see the located test below)."""
        watershed = delineate_and_get_characteristics("WA", 48.57426, -120.37893)

        region_estimates, skipped = estimate_flow_statistics(
            "WA", watershed.characteristics, statistic_group_codes=["PFS"]
        )

        assert region_estimates
        by_code = {r.region_code: r for r in region_estimates}
        assert "GC1750" in by_code
        result = by_code["GC1750"]
        assert result["PK50AEP"].value == pytest.approx(4370.0, rel=0.02)
        assert "DRNAREA" in result["PK50AEP"].equation
        assert result.citation is not None
        assert "Mastin" in result.citation.author

    def test_goat_creek_located_in_peak_region_2(self) -> None:
        """Addendum S5: the polygon picks GC1751, 100%, and only it is estimated."""
        watershed = delineate_and_get_characteristics("WA", 48.57426, -120.37893)
        assert watershed.polygon_geojson is not None

        located = locate_regression_regions("WA", watershed.polygon_geojson)
        assert [(r.code, r.percent_weight) for r in located] == [("GC1751", 100.0)]

        region_estimates, skipped = estimate_flow_statistics(
            "WA",
            watershed.characteristics,
            statistic_group_codes=["PFS"],
            watershed_polygon=watershed.polygon_geojson,
        )

        assert [(r.region_code, r.located) for r in region_estimates] == [("GC1751", True)]
        # 3,290 cfs live 2026-09-27; GC1750's 4,370 was never this basin's answer.
        assert region_estimates[0]["PK50AEP"].value == pytest.approx(3290.0, rel=0.02)
        assert {"PFS:GC1750", "PFS:GC1752", "PFS:GC1753"} <= set(skipped)

    def test_ogeechee_spans_three_regions_with_area_weights(self) -> None:
        """A live multi-region basin: SIR 2014-5030 regions 1/3/4 at 44/32/24%."""
        watershed = delineate_and_get_characteristics(
            "GA", 32.9676527079592, -82.3904040828971, characteristic_codes=["DRNAREA"]
        )
        assert watershed.polygon_geojson is not None

        weights = {
            r.code: r.percent_weight
            for r in locate_regression_regions("GA", watershed.polygon_geojson)
        }

        urban = {c: weights[c] for c in ("GC1539", "GC1541", "GC1542")}
        assert urban == pytest.approx({"GC1539": 44.0, "GC1541": 32.0, "GC1542": 24.0}, abs=2)
        assert sum(urban.values()) == pytest.approx(100.0, abs=0.5)
