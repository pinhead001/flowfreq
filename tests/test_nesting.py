"""Tests for flowfreq.nesting: which gages are flow-connected."""

from __future__ import annotations

from typing import Dict, List

import pandas as pd
import pytest

from flowfreq.nesting import (
    Basin,
    NLDISite,
    connections,
    flow_relation,
    nested_gage_pairs,
    nldi_upstream_sites,
)


def square(x0: float, y0: float, x1: float, y1: float) -> dict:
    return {"type": "Polygon", "coordinates": [[[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]]}


# A main stem draining south: "up" (small) inside "down" (large); "trib" beside
# "up", also inside "down"; "other" on a different stream.
BASINS = {
    "down": Basin("down", 100.0, 47.0, -120.0, square(-120.5, 47.0, -119.5, 48.0)),
    "up": Basin("up", 30.0, 47.6, -120.2, square(-120.4, 47.6, -120.0, 47.95)),
    "trib": Basin("trib", 20.0, 47.5, -119.8, square(-119.95, 47.5, -119.6, 47.9)),
    "other": Basin("other", 50.0, 46.5, -118.0, square(-118.5, 46.5, -117.5, 47.0)),
}


class TestFlowRelation:
    def test_upstream_and_downstream(self):
        assert flow_relation(BASINS["down"], BASINS["up"]) == "upstream"
        assert flow_relation(BASINS["up"], BASINS["down"]) == "downstream"

    def test_different_streams(self):
        assert flow_relation(BASINS["down"], BASINS["other"]) is None
        assert flow_relation(BASINS["up"], BASINS["trib"]) is None

    def test_colocated(self):
        twin = Basin("twin", 101.0, 47.0002, -120.0002, BASINS["down"].geometry)
        assert flow_relation(BASINS["down"], twin) == "colocated"

    def test_a_larger_basin_is_never_upstream(self):
        """A ragged boundary can put a big basin's outlet inside a small one."""
        big = Basin("big", 500.0, 47.7, -120.1, square(-121.0, 46.0, -119.0, 49.0))
        assert flow_relation(BASINS["up"], big) != "upstream"


def _site(site_no: str, comid: int, measure: float = 50.0) -> NLDISite:
    return NLDISite(site_no=site_no, name=site_no, comid=comid, reachcode="r", measure=measure)


# Upstream-with-tributaries navigation, each list including the start site.
NETWORK: Dict[str, List[NLDISite]] = {
    "down": [_site("down", 1), _site("up", 2), _site("trib", 3), _site("lake", 4)],
    "up": [_site("up", 2)],
    "trib": [_site("trib", 3)],
    "other": [_site("other", 9)],
}


def navigator(site_no: str) -> List[NLDISite]:
    return NETWORK[site_no]


class TestNestedGagePairs:
    def test_polygons_and_nldi_agree(self):
        out = nested_gage_pairs(list(BASINS), basins=BASINS, navigator=navigator)
        got = {(r.upstream, r.downstream) for r in out.itertuples()}
        assert got == {("up", "down"), ("trib", "down")}
        assert (out["polygon"] == "nested").all()
        assert (out["nldi"] == "nested").all()
        assert not out["conflict"].any()

    def test_sites_outside_the_list_are_ignored(self):
        out = nested_gage_pairs(["down", "up"], navigator=navigator)
        assert set(zip(out["upstream"], out["downstream"])) == {("up", "down")}

    def test_nldi_only(self):
        out = nested_gage_pairs(["down", "up", "trib", "other"], navigator=navigator)
        assert out["polygon"].isna().all()
        assert len(out) == 2

    def test_polygons_only(self):
        out = nested_gage_pairs(list(BASINS), basins=BASINS, navigator=None)
        assert len(out) == 2 and out["nldi"].isna().all()

    def test_disagreement_is_flagged_not_dropped(self):
        network = dict(NETWORK, down=[_site("down", 1), _site("up", 2)])  # NLDI misses trib
        out = nested_gage_pairs(list(BASINS), basins=BASINS, navigator=network.__getitem__)
        row = out.set_index("upstream").loc["trib"]
        assert row["polygon"] == "nested" and row["nldi"] == "not nested"
        assert row["conflict"]

    def test_same_comid_is_ordered_by_measure(self):
        # NHDPlus measure runs 0 at a flowline's downstream end to 100 upstream.
        a, b = _site("a", 7, measure=20.0), _site("b", 7, measure=80.0)
        net = {"a": [a, b], "b": [a, b]}
        out = nested_gage_pairs(["a", "b"], navigator=net.__getitem__)
        assert list(zip(out["upstream"], out["downstream"])) == [("b", "a")]

    def test_failed_navigation_leaves_the_verdict_unknown(self, caplog):
        def flaky(site_no: str) -> List[NLDISite]:
            if site_no == "up":
                raise ConnectionError("down for maintenance")
            return NETWORK[site_no]

        out = nested_gage_pairs(["down", "up"], navigator=flaky)
        # "down"'s own navigation still lists "up", so the pair is found.
        assert list(zip(out["upstream"], out["downstream"])) == [("up", "down")]
        assert "up" in caplog.text

    def test_needs_a_source(self):
        with pytest.raises(ValueError, match="basins or a navigator"):
            nested_gage_pairs(["a", "b"], basins=None, navigator=None)


class TestConnections:
    def test_upstream_is_the_donor_and_conflicts_are_dropped_on_request(self):
        pairs = pd.DataFrame(
            {
                "upstream": ["u1", "u2"],
                "downstream": ["d", "d"],
                "polygon": ["nested", "nested"],
                "nldi": ["nested", "not nested"],
                "relation": ["nested", "nested"],
                "conflict": [False, True],
            }
        )
        assert connections(pairs) == [("u1", "d"), ("u2", "d")]
        assert connections(pairs, drop_conflicts=True) == [("u1", "d")]

    def test_colocated_excluded_by_default(self):
        pairs = pd.DataFrame(
            {
                "upstream": ["a"],
                "downstream": ["b"],
                "polygon": ["colocated"],
                "nldi": [None],
                "relation": ["colocated"],
                "conflict": [False],
            }
        )
        assert connections(pairs) == []
        assert connections(pairs, include_colocated=True) == [("a", "b")]


class TestNldiParsing:
    def test_parses_and_strips_the_agency_prefix(self, monkeypatch):
        payload = {
            "features": [
                {
                    "properties": {
                        "identifier": "USGS-12447383",
                        "name": "METHOW RIVER ABOVE GOAT CREEK NEAR MAZAMA, WA",
                        "comid": 24383455,
                        "reachcode": "17020008000366",
                        "measure": 27.6592,
                    }
                }
            ]
        }

        class Response:
            status_code = 200

            def raise_for_status(self) -> None:
                pass

            def json(self) -> dict:
                return payload

        seen = {}

        def fake(method, url, **kw):
            seen["url"] = url
            seen["params"] = kw.get("params")
            return Response()

        monkeypatch.setattr("flowfreq.nesting._request_with_backoff", fake)
        sites = nldi_upstream_sites("12448500")
        assert sites == [
            NLDISite(
                "12447383",
                "METHOW RIVER ABOVE GOAT CREEK NEAR MAZAMA, WA",
                24383455,
                "17020008000366",
                27.6592,
            )
        ]
        assert seen["url"].endswith("/nwissite/USGS-12448500/navigation/UT/nwissite")
        assert float(seen["params"]["distance"]) > 1000


@pytest.mark.requires_network
def test_methow_nesting_live():
    """Goat Creek and the Chewuch are above Winthrop; the Twisp joins below it."""
    out = nested_gage_pairs(["12448500", "12447383", "12448000", "12448998"])
    got = set(zip(out["upstream"], out["downstream"]))
    assert ("12447383", "12448500") in got
    assert ("12448000", "12448500") in got
    assert not any("12448998" in pair for pair in got)
