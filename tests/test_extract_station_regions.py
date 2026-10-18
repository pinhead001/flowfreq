"""tools/extract_station_regions.py: parsers on small synthetic inputs in each report's layout."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def ext():
    spec = importlib.util.spec_from_file_location(
        "extract_station_regions", REPO_ROOT / "tools" / "extract_station_regions.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


WY_TEXT = """
   Table 10. Basin characteristics, selected streamflow-gaging stations, Wyoming and surrounding states
      Map        Station        Latitude      Longitude            (square
     number      number        (degrees)      (degrees)     Region miles)
          1     06036905        44.6203        110.8622         1        282         7,990
          2     06037000        44.6494        110.7839         1        118         7,910
         16     06207600        45.1622        108.8233         2          0.85      4,410
"""

NM_TEXT = """
Appendix 1. Selected basin, climatic, and flood characteristics; maximum peak discharge
                       Northeast plains flood region (1)
 07153500      Dry Cimarron River near Guy, New Mexico          611     0.1132      16.85     1
                                                                500
                       Northwest plateau flood region (2)
09367680       Chaco Wash at Chaco Canyon National Monument, New Mexico    611   0.0710   6,700    62
09367680       Catron Wash at Bisti Trading Post, New Mexico      26.5    0.1794     6,970     65
09367680          573       1,080        1,560        2,380       3,180        4,160     62
"""


def _w(x0, y0, x1, text):
    return f'<word xMin="{x0}" yMin="{y0}" xMax="{x1}" yMax="{y0 + 6}">{text}</word>'


NV_BBOX = "<page " + "".join(
    [
        _w(290, 123, 311, "Flood"),
        _w(290, 133, 316, "region"),
        _w(333, 103, 354, "Lati-"),
        # 10243700 Cleve Creek near Ely, Nev.   6   39.216
        _w(87, 200, 120, "10243700"),
        _w(125, 200, 146, "Cleve"),
        _w(151, 200, 172, "Creek"),
        _w(177, 200, 194, "near"),
        _w(199, 200, 215, "Ely,"),
        _w(303, 200, 307, "6"),
        _w(334, 200, 341, "39"),
        _w(343, 200, 358, ".216"),
        # 10172885 Great Salt Lake Desert Tributary No. 2   6  (the "2" is in the name)
        _w(87, 220, 120, "10172885"),
        _w(125, 220, 146, "Desert"),
        _w(260, 220, 270, "No."),
        _w(275, 220, 279, "2"),
        _w(303, 220, 307, "6"),
        # OCR comma in the station number
        _w(87, 240, 125, "10261,530"),
        _w(130, 240, 146, "Pine"),
        _w(300, 240, 309, "10"),
    ]
)


def test_parse_wy(ext):
    rows = ext.parse_wy(WY_TEXT)
    assert [(r["site_no"], r["report_region"], r["latitude"]) for r in rows] == [
        ("06036905", "1", "44.6203"),
        ("06037000", "1", "44.6494"),
        ("06207600", "2", "45.1622"),
    ]
    assert rows[2]["where"] == "Table 10, map number 16"


def test_parse_nm_keeps_a_number_printed_for_two_gages(ext):
    rows = ext.parse_nm(NM_TEXT)
    assert [(r["site_no"], r["report_region"], r["where"][-2:].strip()) for r in rows] == [
        ("07153500", "1", "1"),
        ("09367680", "2", "62"),
        ("09367680", "2", "65"),
    ]
    assert rows[2]["report_name"] == "Catron Wash at Bisti Trading Post, New Mexico"


def test_parse_nv_reads_the_region_column(ext):
    rows = ext.parse_nv(NV_BBOX, 155)
    got = [(r["site_no"], r["report_region"], r["latitude"]) for r in rows]
    assert got == [("10243700", "6", "39.216"), ("10172885", "6", ""), ("10261530", "10", "")]
    assert rows[0]["report_name"] == "Cleve Creek near Ely,"
    assert rows[0]["where"] == "data section, PDF p. 155"


def test_parse_nv_refuses_a_row_without_one_region(ext):
    bad = NV_BBOX.replace(_w(303, 200, 307, "6"), "")
    with pytest.raises(ValueError, match="no single region cell for 10243700"):
        ext.parse_nv(bad, 155)


def test_check_counts(ext):
    wy = [{"report_region": "1"}] * 139
    with pytest.raises(ValueError, match="WY: parsed"):
        ext.check_counts("WY", wy)
    nv = (
        [{"report_region": "6"}] * 80
        + [{"report_region": "10"}] * 104
        + [
            {"report_region": str(r)}
            for r, n in ((1, 165), (2, 108), (3, 35), (5, 37))
            for _ in range(n)
        ]
    )
    assert ext.check_counts("NV", nv)[6] == 80
    with pytest.raises(ValueError, match="NV region 6: 79 listed"):
        ext.check_counts("NV", nv[1:])
