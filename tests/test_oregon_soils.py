"""The STATSGO2 Ksat helpers in tools/build_oregon_basin_characteristics.py (``soils`` stage)."""

from __future__ import annotations

import importlib.util
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
