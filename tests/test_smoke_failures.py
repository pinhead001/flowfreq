"""Each ungaged failure in the 2026-10-30 nationwide smoke run (docs/NATIONWIDE_SMOKE.md),
pinned as an offline test of the behavior it showed.

The smoke run picks one gage per jurisdiction and estimates at its coordinates. Where that
failed for a reason that is the library working as intended (a regulated basin, a basin
too large for the equations, a point off the stream network, a region StreamStats does
not delineate), the test below fixes that behavior, and ``tools/nationwide_smoke_points.csv``
moves the smoke point somewhere it should pass. WY (LONG_OUT from the pour point) and MN
(the Rainy River Basin) were defects, fixed with their own tests:
``test_ungaged.py::test_estimate_at_site_derives_long_out_from_the_pour_point`` and
``test_streamstats.py::...::test_minnesota_also_queries_the_rainy_river_basin``.
"""

from __future__ import annotations

import copy
import json
from typing import Any
from unittest.mock import Mock, patch

import pytest
import requests

from flowfreq.regression.equations import OutOfRangeError, evaluate
from flowfreq.regression.library import EquationsUnavailable, load_state
from flowfreq.streamstats import (
    DegenerateDelineationError,
    RegionNotDelineatedError,
    RegressionRegionLocation,
    RegulatedDelineationError,
    UnsnappablePointError,
    UnsupportedRegionError,
    delineate_watershed,
    snap_point,
)
from flowfreq.ungaged import assemble
from tests.fixtures.streamstats_responses import (
    DELINEATE_SSHYDRO_WITH_WARNING,
    SNAP_GOOD,
    SNAP_UNSNAPPABLE,
)

KW: dict = dict(aeps=None, allow_extrapolation=False, weighting_space="log")


def _response(json_data: Any = None, status_code: int = 200, text: str = "") -> Mock:
    r = Mock()
    r.status_code = status_code
    r.headers = {}
    r.url = "https://mocked.example/"
    r.text = text or (json.dumps(json_data) if json_data is not None else "")
    if json_data is not None:
        r.json.return_value = json_data
    else:
        r.json.side_effect = ValueError("no JSON")
    if status_code >= 400:
        r.raise_for_status.side_effect = requests.HTTPError(f"{status_code} error")
    else:
        r.raise_for_status.return_value = None
    return r


def _region(code: str, percent: float = 100.0) -> RegressionRegionLocation:
    return RegressionRegionLocation(
        region_id=0, code=code, name=code, percent_weight=percent, area_sq_mi=percent
    )


# FL 02315500, MI 04124000, NE 06784000: the snap service has no stream grid there.
@pytest.mark.parametrize(
    "state,body",
    [
        ("MI", "No data found for region='MI' ; Stream grid raster not found for region='MI'"),
        ("FL", "No data found for region='FL' ; Stream grid raster not found for region='FL'"),
        ("NE", "No data found for region='NE' ; Stream grid raster not found for region='NE'"),
    ],
)
def test_region_without_a_stream_grid(state: str, body: str) -> None:
    text = json.dumps({"detail": body})
    with patch(
        "flowfreq.streamstats.requests.get", return_value=_response(status_code=404, text=text)
    ):
        with pytest.raises(RegionNotDelineatedError, match="no stream grid") as info:
            snap_point(state, 44.4, -85.7)
    assert isinstance(info.value, UnsupportedRegionError)


def test_other_snap_404_is_still_an_http_error() -> None:
    text = json.dumps({"detail": "Not Found"})
    with patch(
        "flowfreq.streamstats.requests.get", return_value=_response(status_code=404, text=text)
    ):
        with pytest.raises(requests.HTTPError):
            snap_point("VT", 44.4, -72.7)


# ND 06354000 Cannonball River: StreamStats flags the basin REGULATED.
def _with_warning(message: str) -> Any:
    data = copy.deepcopy(DELINEATE_SSHYDRO_WITH_WARNING)

    def walk(obj: Any) -> None:
        if isinstance(obj, dict):
            for k, v in obj.items():
                if k.lower() == "warningmsg" and isinstance(v, str) and v.strip():
                    obj[k] = message
                else:
                    walk(v)
        elif isinstance(obj, list):
            for item in obj:
                walk(item)

    walk(data)
    return data


def test_regulated_basin_is_its_own_error() -> None:
    responses = [_response(SNAP_GOOD), _response(_with_warning("REGULATED: CANNONBALL RIVER"))]
    with patch("flowfreq.streamstats.requests.get", side_effect=responses):
        with pytest.raises(RegulatedDelineationError, match="CANNONBALL") as info:
            delineate_watershed("ND", 46.376, -100.934)
    assert isinstance(info.value, DegenerateDelineationError)


def test_other_delineation_warnings_stay_generic() -> None:
    responses = [_response(SNAP_GOOD), _response(DELINEATE_SSHYDRO_WITH_WARNING)]
    with patch("flowfreq.streamstats.requests.get", side_effect=responses):
        with pytest.raises(DegenerateDelineationError) as info:
            delineate_watershed("WA", 48.584, -120.370)
    assert not isinstance(info.value, RegulatedDelineationError)


# AK 15072000: the gage coordinates do not snap.
def test_unsnappable_alaska_gage() -> None:
    with patch(
        "flowfreq.streamstats.requests.get", return_value=_response(SNAP_UNSNAPPABLE)
    ) as get:
        with pytest.raises(UnsnappablePointError, match="move the point"):
            delineate_watershed("AK", 55.391604, -131.195583)
    assert get.call_count == 1


# MT 06191500 Yellowstone at Corwin Springs: CONTDA 2616 mi2 is past the SW limit.
def test_montana_southwest_refuses_a_basin_past_its_contda_limit() -> None:
    lib = load_state("MT")
    eq = lib.equation("SW", 0.667)
    limit = next(v for v in eq.variables if v.code == "CONTDA").maximum
    assert limit == pytest.approx(2472.17)
    with pytest.raises(OutOfRangeError, match="CONTDA"):
        evaluate(eq, {"CONTDA": 2616.0, "EL6000": 40.0})
    # The replacement smoke point, Willow Creek near Harrison (85.9 mi2), is inside it.
    assert evaluate(eq, {"CONTDA": 85.9, "EL6000": 40.0}).flow_cfs > 0


# OR 14137000: ss-hydro returned none of the western Oregon equations' climate and soil
# characteristics, so every matched region is refused with what it lacks.
def test_oregon_names_the_missing_characteristics() -> None:
    lib = load_state("OR")
    with pytest.raises(EquationsUnavailable, match="I24H2Y") as info:
        assemble(lib, [_region("GC729")], {"DRNAREA": 263.0}, **KW)
    for code in ("JANMAXT2K", "WATCAPORC", "SOILPERM"):
        assert code in str(info.value)


# AZ 09484000 Sabino Creek near Tucson (35.5 mi2), gaged 1% flood 14,403 cfs: the 2026-10-30
# smoke forced GC1618 (High Elevation), which NSS uses only at a mean basin elevation of
# 7,500 ft or more; Sabino's is about 5,980 ft, so region 5 (GC1622) applies.
def test_arizona_sabino_creek_belongs_to_region_5() -> None:
    lib = load_state("AZ")
    high = evaluate(lib.equation("1", 0.01), {"CONTDA": 35.5})
    se_basin = evaluate(lib.equation("5", 0.01), {"CONTDA": 35.5})
    assert high.flow_cfs == pytest.approx(1468, rel=0.01)  # the smoke run's 0.10 ratio
    assert se_basin.flow_cfs == pytest.approx(8604, rel=0.01)
    assert se_basin.interval is not None
    assert se_basin.interval[0] < 14403 < se_basin.interval[1]
