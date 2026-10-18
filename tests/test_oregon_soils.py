"""Soil and geology helpers in tools/build_oregon_basin_characteristics.py.

The ``soils`` (STATSGO2 Ksat) and ``geology`` (SGMC young volcanics, NHD lakes)
stages; no network access.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import ModuleType

import numpy as np
import pandas as pd
import pytest

TOOL = Path(__file__).resolve().parents[1] / "tools" / "build_oregon_basin_characteristics.py"


@pytest.fixture(scope="module")
def tool() -> ModuleType:
    spec = importlib.util.spec_from_file_location("build_oregon_basin_characteristics", TOOL)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_wkt_polygon_with_hole(tool):
    g = tool.wkt_to_geojson("POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0), (1 1, 2 1, 2 2, 1 1))")
    assert g["type"] == "Polygon"
    assert len(g["coordinates"]) == 2
    assert g["coordinates"][0][1] == [4.0, 0.0]
    assert g["coordinates"][1][-1] == [1.0, 1.0]


def test_wkt_multipolygon(tool):
    g = tool.wkt_to_geojson(
        "MULTIPOLYGON (((0 0, 1 0, 1 1, 0 0)), ((5 5, 6 5, 6 6, 5 5), (5.2 5.1, 5.5 5.1, 5.5 5.3, 5.2 5.1)))"
    )
    assert g["type"] == "MultiPolygon"
    assert [len(p) for p in g["coordinates"]] == [1, 2]
    assert g["coordinates"][1][0][0] == [5.0, 5.0]


def test_wkt_rejects_other_geometries(tool):
    with pytest.raises(ValueError, match="unsupported"):
        tool.wkt_to_geojson("LINESTRING (0 0, 1 1)")


def test_mapunit_ksat_weights_by_thickness_then_component_percent(tool):
    rows = pd.DataFrame(
        [
            # component 1 (60 %): 10 cm at 100 um/s over 30 cm at 20 -> 40 um/s
            (1, 11, 60, 0, 10, 100.0),
            (1, 11, 60, 10, 40, 20.0),
            # component 2 (30 %): one horizon at 10
            (1, 12, 30, 0, 50, 10.0),
            # rock outcrop (10 %): no horizon data, left out of the weights
            (1, 13, 10, np.nan, np.nan, np.nan),
            # map unit 2: water only
            (2, 21, 100, np.nan, np.nan, np.nan),
        ],
        columns=["mukey", "cokey", "comppct_r", "hzdept_r", "hzdepb_r", "ksat_r"],
    )
    k = tool.mapunit_ksat(rows)
    assert k.loc[1] == pytest.approx((60 * 40 + 30 * 10) / 90)
    assert 2 not in k.index


def test_soilksat_raster_is_converted_to_inches_per_hour(tool):
    _, scale, offset = tool.RASTERS["SOILKSAT"]
    # 1 um/s = 3.6 mm/hr = 0.1417 in./hr
    assert 1.0 * scale + offset == pytest.approx(0.14173, rel=1e-4)
    assert "SOILKSAT" in tool.SOURCES


def test_arcgis_features_reads_cached_pages_until_a_short_one(tool, tmp_path):
    # A full page of 2, then a short page of 1, both cached: no network.
    feat = {"type": "Feature", "geometry": None, "properties": {}}
    (tmp_path / "page_000000.json").write_text(json.dumps([feat, feat]))
    (tmp_path / "page_000002.json").write_text(json.dumps([feat]))
    got = tool.arcgis_features("http://invalid", "1=1", (0, 0, 1, 1), "fid", tmp_path, 0.0, page=2)
    assert len(got) == 3


def test_geology_selects_young_volcanic_units_and_lakes(tool):
    assert "Igneous, volcanic" in tool.SGMC_WHERE
    assert "Pliocene" in tool.SGMC_WHERE and "Quaternary" in tool.SGMC_WHERE
    assert "390" in tool.NHD_WHERE and "436" in tool.NHD_WHERE
    for name in ("VOLCQP", "LAKEAREA"):
        assert tool.RASTERS[name][1] == 100.0  # mean of a 0/1 mask -> percent
        assert name in tool.SOURCES


def test_ogdc_selects_volcanic_units_no_older_than_pliocene(tool):
    assert set(tool.OGDC_VOLCANIC) == {
        "volcanic rocks",
        "vent and pyroclastic rocks",
        "volcaniclastic rocks",
    }
    assert "Quaternary" in tool.OGDC_YOUNG and "Pliocene/Quaternary" in tool.OGDC_YOUNG
    # ranges reaching into the Miocene (early High Cascade, Deschutes Formation) are out
    assert not any("Miocene" in a for a in tool.OGDC_YOUNG)
    for v in (*tool.OGDC_VOLCANIC, *tool.OGDC_YOUNG):
        assert f"'{v}'" in tool.OGDC_WHERE
    assert tool.OGDC_WHERE.startswith("ThematicRockType IN (")
    assert tool.RASTERS["VOLCQP_OG"][1] == 100.0
    assert "OGDC-8" in tool.SOURCES["VOLCQP_OG"]


def test_merge_state_mask_takes_the_state_map_inside_and_sgmc_outside(tool):
    inside = np.array([[1, 1], [0, 0]], dtype="uint8")
    state = np.array([[1, 0], [1, 1]], dtype="uint8")
    sgmc = np.array([[0, 1], [0, 1]], dtype="uint8")
    got = tool.merge_state_mask(inside, state, sgmc)
    assert got.tolist() == [[1, 0], [0, 1]]
    assert got.dtype == np.uint8
    with pytest.raises(ValueError, match="grids differ"):
        tool.merge_state_mask(inside, state, sgmc[:1])
