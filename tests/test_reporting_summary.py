"""The `site_summary` report section: characteristics, a map, the flows by method and a plot.

Nothing here makes a request. The basin outlines are squares and the characteristics are
made up, since the section prints both as supplied.
"""

from __future__ import annotations

from typing import Any, Dict, List

import pytest

from flowfreq.reporting import (
    SECTION_KINDS,
    Figure,
    Paragraph,
    ReportError,
    ReportSpec,
    SectionSpec,
    Table,
    build_report,
    build_section,
)
from flowfreq.reporting.session import MAIN_TABLES
from flowfreq.reporting.summary import (
    SUMMARY_AEPS,
    SUMMARY_CHARACTERISTICS_CAPTION,
    SUMMARY_FLOWS_CAPTION,
)
from tests.test_reporting_comparison import ESTIMATES, SITE


def _square(lon: float, lat: float, half: float) -> Dict[str, Any]:
    ring = [
        [lon - half, lat - half],
        [lon + half, lat - half],
        [lon + half, lat + half],
        [lon - half, lat + half],
        [lon - half, lat - half],
    ]
    return {"type": "Polygon", "coordinates": [ring]}


LOCATIONS: List[Dict[str, Any]] = [
    {
        "name": "Gage",
        "lat": 36.04,
        "lon": -88.23,
        "characteristics": {"DRNAREA": 205.0, "PRECIP": 52.5},
        "basin": _square(-88.3, 36.0, 0.2),
    },
    {
        "name": "Site",
        "lat": 36.10,
        "lon": -88.30,
        "characteristics": {"DRNAREA": 150.0, "ELEV": 480.0},
        "basin": {
            "type": "MultiPolygon",
            "coordinates": [_square(-88.35, 36.0, 0.15)["coordinates"]],
        },
        "marker": "*",
    },
]
SUMMARY: Dict[str, Any] = {
    "site_name": SITE,
    "locations": LOCATIONS,
    "estimates": ESTIMATES,
    "characteristic_names": {"DRNAREA": "Drainage area (square miles)"},
    "characteristics_source": "made up for the tests.",
}


def _table(blocks: List[Any], caption: str) -> Table:
    (table,) = [b for b in blocks if isinstance(b, Table) and b.caption == caption]
    return table


@pytest.fixture(scope="module")
def built() -> Any:
    return build_section("site_summary", dict(SUMMARY))


def test_kind_is_registered():
    kind = SECTION_KINDS["site_summary"]
    assert kind.required == ("locations", "estimates")
    assert MAIN_TABLES["site_summary"] == SUMMARY_FLOWS_CAPTION
    assert SUMMARY_AEPS[0] == 0.667 and SUMMARY_AEPS[-1] == 0.002


def test_characteristics_are_printed_as_supplied(built):
    table = _table(built.blocks, SUMMARY_CHARACTERISTICS_CAPTION)
    assert table.header == ("Characteristic", "Code", "Gage", "Site")
    assert table.rows == (
        ("Drainage area (square miles)", "DRNAREA", "205", "150"),
        ("PRECIP", "PRECIP", "52.5", "-"),
        ("ELEV", "ELEV", "-", "480"),
    )
    text = "\n".join(b.text for b in built.blocks if isinstance(b, Paragraph))
    assert "printed as supplied; flowfreq did not compute them" in text
    assert "Source: made up for the tests." in text


def test_flows_are_the_comparison_sections_own(built):
    table = _table(built.blocks, SUMMARY_FLOWS_CAPTION)
    assert table.header[:2] == (
        "Annual exceedance probability (percent)",
        "Recurrence interval (years)",
    )
    assert table.header[2:] == tuple(f"{e['label']} (cfs)" for e in ESTIMATES)
    compared = build_section("comparison", {"estimates": ESTIMATES, "detail": False})
    theirs = {
        row[0]: row[1:4]
        for row in _table(compared.blocks, "Flood quantiles at the site by each method.").rows
    }
    for row in table.rows:
        assert row[2:] == theirs[row[0]]
    # Only the summary AEPs every estimate gives: Tennessee's set has no 1.5-year equation.
    assert [(row[0], row[1]) for row in table.rows] == [
        ("50", "2"),
        ("20", "5"),
        ("10", "10"),
        ("4", "25"),
        ("2", "50"),
        ("1", "100"),
        ("0.5", "200"),
        ("0.2", "500"),
    ]


def test_the_1_5_year_row_joins_0_67_and_0_667():
    """A Bulletin 17C fit tabulates 0.67 and Idaho's regression 0.667: one row, and said."""
    idaho = {
        "label": "Idaho regression",
        "kind": "ungaged",
        "options": {
            "state": "ID",
            "lat": 48.8,
            "lon": -116.2,
            "region": "1_2",
            "characteristics": {"DRNAREA": 205.0, "PRECPRIS10": 34.0},
        },
    }
    built = build_section(
        "site_summary", {"locations": LOCATIONS, "estimates": [ESTIMATES[0], idaho]}
    )
    table = _table(built.blocks, SUMMARY_FLOWS_CAPTION)
    assert (table.rows[0][0], table.rows[0][1]) == ("66.7", "1.5")
    text = "\n".join(b.text for b in built.blocks if isinstance(b, Paragraph))
    assert "slightly different probabilities" in text and "(0.667, 0.67)" in text


def test_a_map_and_a_curve_plot(built):
    figures = [b for b in built.blocks if isinstance(b, Figure)]
    assert len(figures) == 2
    assert figures[0].caption.startswith(
        "Locations used in this study, with the watershed above Gage, Site."
    )
    assert "no base map" in figures[0].caption
    assert "Regional regression, 90-percent prediction limits" in figures[1].caption
    assert all(f.png[:8] == b"\x89PNG\r\n\x1a\n" for f in figures)
    assert built.title == f"Summary: {SITE}"


def test_curves_can_be_chosen():
    built = build_section(
        "site_summary", {**SUMMARY, "curves": ["Regional regression", "Transposed from the gage"]}
    )
    caption = [b for b in built.blocks if isinstance(b, Figure)][1].caption
    assert "Gage, 90-percent" not in caption and "Transposed from the gage" in caption
    with pytest.raises(ReportError, match="curves must be a list of estimate labels"):
        build_section("site_summary", {**SUMMARY, "curves": ["Nope"]})


def test_locations_need_no_basin_or_characteristics():
    bare = [{"name": "Site", "lat": 36.1, "lon": -88.3}]
    built = build_section("site_summary", {"locations": bare, "estimates": ESTIMATES[:1]})
    assert not [b for b in built.blocks if isinstance(b, Table) and b.header[0] == "Characteristic"]
    assert [b for b in built.blocks if isinstance(b, Figure)][0].caption.startswith(
        "Locations used in this study. Outlines"
    )
    assert built.title == "Summary"


@pytest.mark.parametrize(
    "change,said",
    [
        ({"locations": []}, "locations must be a non-empty list"),
        ({"locations": [{"name": "A", "lat": 1.0}]}, "needs lat and lon"),
        ({"locations": [{"name": "", "lat": 1.0, "lon": 2.0}]}, "needs a name"),
        ({"locations": [{"name": "A", "lat": 1.0, "lon": 2.0, "area": 3}]}, "unknown key"),
        (
            {"locations": [{"name": "A", "lat": 1.0, "lon": 2.0, "basin": {"type": "Point"}}]},
            "GeoJSON Polygon or MultiPolygon",
        ),
        (
            {"locations": [{"name": "A", "lat": 1.0, "lon": 2.0, "characteristics": {"X": "wet"}}]},
            "a characteristic is not a number",
        ),
        ({"locations": [LOCATIONS[0], LOCATIONS[0]]}, "two locations have the same name"),
        ({"estimates": []}, "estimates must be a non-empty list"),
        ({"aeps": [0.995]}, "not every estimate has a quantile at AEP 0.995"),
    ],
)
def test_refusals(change, said):
    with pytest.raises(ReportError, match=said):
        build_section("site_summary", {**SUMMARY, **change})


def test_in_a_report_with_a_narrative(tmp_path):
    spec = ReportSpec(
        title="Big Sandy River",
        sections=[
            SectionSpec(
                "site_summary",
                title="Executive summary",
                options=dict(SUMMARY),
                narrative="Limitations.\n\n- The exponents are made up.",
            )
        ],
    )
    report = build_report(spec)
    text = report.to_markdown()
    assert "## 1. Executive summary" in text and SUMMARY_FLOWS_CAPTION in text
    assert "The exponents are made up." in text
    saved = report.spec.save(tmp_path / "report.json")
    assert build_report(ReportSpec.load(saved)).to_markdown() == text
    report.save(tmp_path / "report.docx")
