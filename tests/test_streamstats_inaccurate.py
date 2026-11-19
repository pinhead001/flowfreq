"""``allow_inaccurate``: a basin StreamStats warns "... MAY BE INACCURATE" about.

StreamStats says this of a basin that runs past its data, the Moyie River above Eastport
(Idaho, mostly in British Columbia) for one; the message below is the one it sent for that
basin on 2026-11-18. By default the basin is refused. With the option the delineation is
kept, the warning is an advisory, and the characteristics it names are withheld.
"""

from __future__ import annotations

import copy
from unittest.mock import patch

import pytest

from flowfreq.streamstats import (
    DegenerateDelineationError,
    RegulatedDelineationError,
    StreamStatsCache,
    delineate_and_get_characteristics,
    delineate_watershed,
    inaccurate_codes,
)
from flowfreq.ungaged import estimate_at_site
from tests.fixtures.streamstats_responses import (
    DELINEATE_SSHYDRO_GOOD,
    HYDRO_CHARACTERISTICS_GOOD,
    SNAP_GOOD,
)
from tests.test_streamstats import _mock_response

MOYIE = (
    "FOREST, SLOP30_10M, PRECIP, NFSL30_10M, ELEV, BSLDEM10M, NFSL30_30M, SLOP30_30M, AND "
    "BSLDEM30M MAY BE INACCURATE"
)
LON, LAT = SNAP_GOOD["input"]["coordinates"]


def _warned(message: str = MOYIE):
    data = copy.deepcopy(DELINEATE_SSHYDRO_GOOD)
    data["WarningMsg"] = message
    return data


def _hydro_with(*codes: str):
    """The good characteristics response with some of the named codes added."""
    hydro = copy.deepcopy(HYDRO_CHARACTERISTICS_GOOD)
    items = hydro["parameters"] if isinstance(hydro, dict) else hydro
    for code in codes:
        items.append({"code": code, "name": code, "description": code, "unit": "", "value": 12.5})
    return hydro


def _delineate(sshydro, hydro=HYDRO_CHARACTERISTICS_GOOD, **kwargs):
    with (
        patch(
            "flowfreq.streamstats.requests.get",
            side_effect=[_mock_response(SNAP_GOOD), _mock_response(sshydro)],
        ),
        patch("flowfreq.streamstats.requests.post", return_value=_mock_response(hydro)) as post,
    ):
        return delineate_and_get_characteristics("WA", LAT, LON, **kwargs), post


@pytest.mark.parametrize(
    "message,codes",
    [
        (
            MOYIE,
            [
                "FOREST",
                "SLOP30_10M",
                "PRECIP",
                "NFSL30_10M",
                "ELEV",
                "BSLDEM10M",
                "NFSL30_30M",
                "SLOP30_30M",
                "BSLDEM30M",
            ],
        ),  # fmt: skip
        ("PRECIP AND ELEV MAY BE INACCURATE.", ["PRECIP", "ELEV"]),
        ("elev may be inaccurate", ["ELEV"]),
        ("Delineated basin may contain karst", None),
        ("REGULATED: SOME RIVER; PRECIP MAY BE INACCURATE", None),
        ("PRECIP MAY BE INACCURATE and the basin is regulated", None),
    ],
)
def test_inaccurate_codes_reads_only_the_whole_message(message, codes):
    assert inaccurate_codes(message) == codes


def test_refused_by_default_and_says_how_to_keep_it():
    with pytest.raises(DegenerateDelineationError, match="MAY BE INACCURATE") as refusal:
        _delineate(_warned())
    assert "allow_inaccurate=True" in str(refusal.value)
    with pytest.raises(DegenerateDelineationError) as other:
        _delineate(_warned("Streamflow altered by upstream reservoirs."))
    assert "allow_inaccurate" not in str(other.value)


def test_kept_with_the_named_characteristics_withheld():
    result, post = _delineate(_warned(), _hydro_with("PRECIP", "ELEV"), allow_inaccurate=True)
    post.assert_called_once()
    assert result.advisories == [MOYIE]
    assert "PRECIP" not in result.characteristics and "ELEV" not in result.characteristics
    assert result.unavailable["PRECIP"].startswith("withheld: StreamStats says PRECIP may be")
    assert "12.5" in result.unavailable["ELEV"]
    # What the warning does not name is returned as computed.
    assert result["DRNAREA"].value == pytest.approx(412.0)
    assert result.polygon_geojson is not None


def test_any_other_warning_is_still_refused_with_the_option():
    with pytest.raises(DegenerateDelineationError, match="reservoirs"):
        _delineate(_warned("Streamflow altered by upstream reservoirs."), allow_inaccurate=True)
    both = _warned()
    both["elsewhere"] = {"WarningMsg": "REGULATED: SOME RIVER"}
    with pytest.raises(RegulatedDelineationError, match="SOME RIVER"):
        _delineate(both, allow_inaccurate=True)


def test_a_kept_result_never_answers_a_default_call(tmp_path):
    cache = StreamStatsCache(tmp_path / "cache.json")
    kept, _ = _delineate(_warned(), cache=cache, allow_inaccurate=True)
    with patch("flowfreq.streamstats.requests.get") as get:
        again = delineate_and_get_characteristics(
            "WA", LAT, LON, cache=cache, allow_inaccurate=True
        )
    get.assert_not_called()
    assert again.advisories == kept.advisories
    # The default call misses the cache, asks StreamStats, and is refused.
    with pytest.raises(DegenerateDelineationError, match="MAY BE INACCURATE"):
        _delineate(_warned(), cache=cache)


def test_the_polygon_alone_can_be_kept_too():
    responses = [_mock_response(SNAP_GOOD), _mock_response(_warned())]
    with patch("flowfreq.streamstats.requests.get", side_effect=responses):
        with pytest.raises(DegenerateDelineationError):
            delineate_watershed("WA", LAT, LON)
    responses = [_mock_response(SNAP_GOOD), _mock_response(_warned())]
    with patch("flowfreq.streamstats.requests.get", side_effect=responses):
        kept = delineate_watershed("WA", LAT, LON, allow_inaccurate=True)
    assert kept.advisories == [MOYIE] and kept.polygon_geojson is not None


def test_estimate_at_site_passes_it_on_and_records_it(monkeypatch):
    from tests.test_ungaged import WA_LAT, WA_LON, Fakes, _wa_site

    fakes = Fakes(_wa_site())
    monkeypatch.setattr("flowfreq.ungaged.delineate_and_get_characteristics", fakes.delineate)
    monkeypatch.setattr("flowfreq.ungaged.locate_regression_regions", fakes.locate)
    est = estimate_at_site(WA_LAT, WA_LON, state="WA", allow_inaccurate=True)
    assert fakes.delineate_calls[0]["allow_inaccurate"] is True
    assert est.provenance["options"]["allow_inaccurate"] is True
    plain = estimate_at_site(WA_LAT, WA_LON, state="WA")
    assert fakes.delineate_calls[1]["allow_inaccurate"] is False
    assert plain.provenance["options"]["allow_inaccurate"] is False


def test_a_regulated_warning_outranks_an_inaccurate_one_listed_first():
    both = _warned()
    both["later"] = {"WarningMsg": "REGULATED: SOME RIVER"}
    with pytest.raises(RegulatedDelineationError, match="SOME RIVER"):
        _delineate(both)
