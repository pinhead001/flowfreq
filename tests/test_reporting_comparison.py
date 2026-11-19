"""The `comparison` report section: several estimates of a site's flood quantiles side by side.

Nothing here makes a request: the gage carries its peaks and the regression its
characteristics. The drainage-area exponents are made up, and cited as made up.
"""

from __future__ import annotations

from typing import Any, Dict, List

import pytest

from flowfreq.reporting import (
    SECTION_KINDS,
    Heading,
    Paragraph,
    ReportError,
    ReportSpec,
    SectionSpec,
    Table,
    build_report,
    build_section,
)
from flowfreq.reporting.sections import COMPARISON_CAPTION, cfs
from flowfreq.reporting.session import MAIN_TABLES, ReportSession
from tests.test_reporting_sites import GAGE, TENNESSEE
from tests.test_reporting_transposed import CITATION, EXPONENTS

TRANSPOSED: Dict[str, Any] = {
    "donor": dict(GAGE),
    "target_area_sqmi": 150.0,
    "exponents": dict(EXPONENTS),
    "detail": False,
}
ESTIMATES: List[Dict[str, Any]] = [
    {"label": "Gage", "kind": "gaged_ffa", "options": dict(GAGE)},
    {"label": "Transposed from the gage", "kind": "transposed", "options": dict(TRANSPOSED)},
    {"label": "Regional regression", "kind": "ungaged", "options": dict(TENNESSEE)},
]
SITE = "Big Sandy River above Bruceton"


def _table(blocks: List[Any], caption: str) -> Table:
    (table,) = [b for b in blocks if isinstance(b, Table) and b.caption == caption]
    return table


def _text(blocks: List[Any]) -> str:
    return "\n".join(b.text for b in blocks if isinstance(b, Paragraph))


def _number(cell: str) -> float:
    return float(cell.replace(",", ""))


@pytest.fixture(scope="module")
def built() -> Any:
    return build_section("comparison", {"site_name": SITE, "estimates": ESTIMATES})


def test_kind_is_registered():
    kind = SECTION_KINDS["comparison"]
    assert kind.required == ("estimates",)
    assert kind.optional == ("site_name", "aeps", "detail")
    assert MAIN_TABLES["comparison"] == COMPARISON_CAPTION


def test_each_column_is_that_methods_own_number(built):
    table = _table(built.blocks, COMPARISON_CAPTION)
    assert table.header == (
        "Annual exceedance probability (percent)",
        "Gage (cfs)",
        "Transposed from the gage (cfs)",
        "Regional regression (cfs)",
        "Largest over smallest",
    )
    by_aep = {row[0]: row for row in table.rows}
    # Only the AEPs every estimate gives: the regression and the exponents stop short of the fit.
    assert "99.5" not in by_aep and "1" in by_aep
    gage = build_section("gaged_ffa", dict(GAGE))
    own = _table(gage.blocks, "Flood quantiles with 90-percent confidence limits.")
    assert by_aep["1"][1] == {row[0]: row for row in own.rows}["1"][2]
    transposed = build_section("transposed", dict(TRANSPOSED))
    moved = _table(
        transposed.blocks, "Flood quantiles transposed from the gage by the drainage-area ratio."
    )
    assert by_aep["1"][2] == {row[0]: row for row in moved.rows}["1"][4]
    regression = build_section("ungaged", dict(TENNESSEE))
    fitted = _table(regression.blocks, "Regression flood quantiles with prediction limits.")
    assert by_aep["1"][3] in {cell for row in fitted.rows if row[0] == "1" for cell in row}
    flows = [_number(c) for c in by_aep["1"][1:4]]
    assert by_aep["1"][4] == f"{max(flows) / min(flows):.2f}"


def test_limits_are_tabulated_and_named(built):
    limits = _table(built.blocks, "Lower and upper limits of each estimate.")
    row = {r[0]: r for r in limits.rows}["1"]
    assert all(" to " in cell for cell in row[1:])
    low, high = (_number(part) for part in row[1].split(" to "))
    flows = _table(built.blocks, COMPARISON_CAPTION)
    assert low < _number({r[0]: r for r in flows.rows}["1"][1]) < high
    text = _text(built.blocks)
    assert "Gage, 90-percent confidence limits" in text
    assert "Regional regression, 90-percent prediction limits" in text
    assert "not the uncertainty of the transposition" in text
    assert "compared, not combined" in text


def test_detail_carries_each_estimate_under_its_label(built):
    headings = [b.text for b in built.blocks if isinstance(b, Heading) and b.level == 2]
    assert headings == ["Gage", "Transposed from the gage", "Regional regression"]
    assert built.title == f"Flood quantiles at {SITE} by 3 methods"
    assert CITATION in built.references and len(built.references) == len(set(built.references))
    brief = build_section("comparison", {"estimates": ESTIMATES[:2], "detail": False})
    assert not any(isinstance(b, Heading) for b in brief.blocks)
    assert brief.title == "Comparison of flood estimates"
    assert "in this report's other sections" in _text(brief.blocks)


def test_aeps_can_be_chosen_and_must_be_shared():
    built = build_section(
        "comparison", {"estimates": ESTIMATES, "aeps": [0.01, 0.1], "detail": False}
    )
    assert [row[0] for row in _table(built.blocks, COMPARISON_CAPTION).rows] == ["1", "10"]
    with pytest.raises(ReportError, match="not every estimate has a quantile at AEP 0.995"):
        build_section("comparison", {"estimates": ESTIMATES, "aeps": [0.995]})


def test_options_as_built_hold_each_resolved_estimate(built):
    kept = built.options["estimates"]
    assert [e["label"] for e in kept] == [e["label"] for e in ESTIMATES]
    assert kept[1]["options"]["donor_area_sqmi"] == 205.0  # resolved from the donor
    again = build_section("comparison", built.options)
    assert _table(again.blocks, COMPARISON_CAPTION) == _table(built.blocks, COMPARISON_CAPTION)


@pytest.mark.parametrize(
    "estimates,said",
    [
        (ESTIMATES[:1], "two or more estimates"),
        ("gage", "two or more estimates"),
        ([ESTIMATES[0], {**ESTIMATES[1], "label": "Gage"}], "the same label"),
        ([ESTIMATES[0], {**ESTIMATES[1], "label": " "}], "estimate 2 needs a label"),
        ([ESTIMATES[0], {**ESTIMATES[1], "kind": "weighted"}], "kind must be one of"),
        ([ESTIMATES[0], {**ESTIMATES[1], "note": "x"}], r"unknown key\(s\) note"),
        (
            [ESTIMATES[0], {"label": "R", "kind": "ungaged", "options": {"state": "TN"}}],
            r"missing option\(s\) lat, lon",
        ),
    ],
)
def test_refusals(estimates, said):
    with pytest.raises(ReportError, match=said):
        build_section("comparison", {"estimates": estimates})


def test_in_a_report_and_a_session(tmp_path):
    spec = ReportSpec(
        title="Big Sandy River",
        sections=[
            SectionSpec(
                "comparison", options={"site_name": SITE, "estimates": ESTIMATES, "detail": False}
            )
        ],
    )
    report = build_report(spec)
    text = report.to_markdown()
    assert COMPARISON_CAPTION in text and cfs(18069.0) in text
    saved = report.spec.save(tmp_path / "report.json")
    assert build_report(ReportSpec.load(saved)).to_markdown() == text
    session = ReportSession("Big Sandy River", output_dir=tmp_path)
    added = session.add_section(
        "comparison", {"site_name": SITE, "estimates": ESTIMATES, "detail": False}
    )
    assert added["result"].startswith(COMPARISON_CAPTION)


def test_the_1_5_year_row_is_kept_across_0_67_and_0_667():
    """A Bulletin 17C fit tabulates 0.67 and Idaho's regression 0.667: one row, not none."""
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
    built = build_section("comparison", {"estimates": [ESTIMATES[0], idaho], "detail": False})
    rows = {row[0]: row for row in _table(built.blocks, COMPARISON_CAPTION).rows}
    assert "66.7" in rows and "80" in rows
    assert "slightly different probabilities" in _text(built.blocks)
    assert "66.7 percent (0.667, 0.67)" in _text(built.blocks)


def test_nearby_aeps_are_never_taken_for_each_other():
    from flowfreq.reporting.sections import _matched_aeps, _row_at

    table = {0.005: {"aep": 0.005}, 0.01: {"aep": 0.01}, 0.67: {"aep": 0.67}}
    assert _row_at(table, 0.002) is None
    assert _row_at(table, 0.667)["aep"] == 0.67
    assert _row_at({}, 0.5) is None
    other = {0.002: {"aep": 0.002}, 0.01: {"aep": 0.01}}
    assert _matched_aeps("x", [table, other], None) == ([0.01], [])
    with pytest.raises(ReportError, match="no AEP in common"):
        _matched_aeps("x", [table, {0.5: {"aep": 0.5}}], None)
