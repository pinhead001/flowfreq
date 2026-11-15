"""The `transposed` report section: a gage's quantiles carried to a site on its stream.

Nothing here makes a request: the donor section carries its peaks. The drainage-area
exponents are made up, and cited as made up, so that no published report is credited
with values it does not hold.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np
import pytest

from flowfreq.regulation import RegulationMismatchError
from flowfreq.reporting import (
    SECTION_KINDS,
    Bullets,
    Heading,
    Paragraph,
    ReportError,
    ReportSpec,
    SectionSpec,
    Table,
    build_report,
    build_section,
)
from flowfreq.reporting.sections import BULLETIN_17C, cfs
from flowfreq.transpose import RegressionExponents, transpose_frequency
from flowfreq.workflow import run_ffa
from tests.test_reporting_sites import FLOWS, GAGE, SKEW, YEARS

CITATION = "Example, A.B., and Author, C.D., 2020, A made-up report for the tests, table 1."
#: Exponents by AEP, as a regional report's table of regression coefficients gives them.
EXPONENTS: Dict[str, Any] = {
    "aeps": [0.5, 0.1, 0.01],
    "exponents": [0.70, 0.66, 0.62],
    "citation": CITATION,
    "region": "Region 1",
}
#: Big Sandy's 205 square miles carried upstream to 150.
TRANSPOSED: Dict[str, Any] = {
    "donor": dict(GAGE),
    "target_area_sqmi": 150.0,
    "target_name": "Big Sandy River above Bruceton",
    "exponents": dict(EXPONENTS),
}
CAPTION = "Flood quantiles transposed from the gage by the drainage-area ratio."
DONOR_HEADING = "Analysis at the donor gage"


def _table(blocks: List[Any]) -> Table:
    (table,) = [b for b in blocks if isinstance(b, Table) and b.caption == CAPTION]
    return table


def _text(blocks: List[Any]) -> str:
    return "\n".join(b.text for b in blocks if isinstance(b, Paragraph))


def _bullets(blocks: List[Any]) -> List[str]:
    return [item for b in blocks if isinstance(b, Bullets) for item in b.items]


def _brief(**change: Any) -> Any:
    """The section without the donor's analysis below it, with some options changed."""
    return build_section("transposed", {**TRANSPOSED, "detail": False, **change})


@pytest.fixture(scope="module")
def transposed() -> Any:
    return build_section("transposed", dict(TRANSPOSED))


class TestTransposed:
    def test_kind_is_registered(self):
        kind = SECTION_KINDS["transposed"]
        assert kind.required == ("donor", "target_area_sqmi", "exponents")
        assert kind.optional == (
            "donor_area_sqmi",
            "target_name",
            "max_aep",
            "allow_out_of_range",
            "donor_site_class",
            "target_site_class",
            "allow_regulation_mismatch",
            "detail",
        )

    def test_table_is_transpose_frequency(self, transposed):
        fit = run_ffa(
            np.array(FLOWS), np.array(YEARS), site_no="03606500", stationarity=False, **SKEW
        )
        expected = transpose_frequency(
            fit["b17c"].results,
            205.0,
            150.0,
            RegressionExponents(
                aeps=EXPONENTS["aeps"], exponents=EXPONENTS["exponents"], citation=CITATION
            ),
        )
        limits = expected.confidence_limits.set_index("aep")
        table = _table(transposed.blocks)
        assert table.header == (
            "Annual exceedance probability (percent)",
            "Donor flow (cfs)",
            "Drainage-area exponent",
            "Factor applied",
            "Transposed flow (cfs)",
            "Lower 90-percent limit (cfs)",
            "Upper 90-percent limit (cfs)",
            "Exponent basis",
        )
        assert len(table.align) == 8 and table.align[-1] == "left"
        assert len(table.rows) == len(expected.quantiles) == 14
        for row, (_, quantile) in zip(table.rows, expected.quantiles.iterrows()):
            limit = limits.loc[quantile["aep"]]
            assert row == (
                f"{quantile['aep'] * 100:.3g}",
                cfs(quantile["donor_flow_cfs"]),
                f"{quantile['exponent']:.3f}",
                f"{quantile['area_ratio_factor']:.4f}",
                cfs(quantile["flow_cfs"]),
                cfs(limit["lower_5pct"]),
                cfs(limit["upper_5pct"]),
                quantile["source"],
            )

    def test_the_one_percent_flow(self, transposed):
        by_aep = {row[0]: row for row in _table(transposed.blocks).rows}
        # 18,069 cfs at the gage, times (150 / 205) ** 0.62.
        assert (150.0 / 205.0) ** 0.62 == pytest.approx(0.8239, abs=5e-5)
        assert by_aep["1"][1:] == (
            "18,069",
            "0.620",
            "0.8239",
            "14,888",
            "11,539",
            "22,339",
            "published",
        )
        assert by_aep["20"][2] == "0.674" and by_aep["20"][-1] == "interpolated"
        assert by_aep["0.2"][2] == "0.620" and by_aep["0.2"][-1] == "extrapolated"

    def test_title_method_and_references(self, transposed):
        assert transposed.title == "Transposed flood frequency at Big Sandy River above Bruceton"
        text = _text(transposed.blocks)
        assert (
            "Flood quantiles at Big Sandy River above Bruceton were transposed from the donor "
            "gage, USGS 03606500 Big Sandy River at Bruceton, TN, by the drainage-area ratio"
        ) in text
        assert "The exponents are from Example and Author (2020)." in text
        assert "No distribution was fitted at the site" in text
        assert "Q_target" not in text
        assert transposed.references == [BULLETIN_17C, CITATION]

    def test_inputs_and_caveats(self, transposed):
        bullets = _bullets(transposed.blocks)
        assert "Donor: USGS 03606500, 205.0 sq mi" in bullets
        assert "Target: Big Sandy River above Bruceton, 150.0 sq mi" in bullets
        assert "Area ratio: 0.7317" in bullets
        assert f"Exponent source: {CITATION}" in bullets and "Region: Region 1" in bullets
        assert any(
            b.startswith("**8 of 14 quantiles use an extrapolated exponent**") for b in bullets
        )
        assert any(b.startswith("**Regulation was not checked**") for b in bullets)
        assert any(b.startswith("Confidence limits are the donor's") for b in bullets)

    def test_detail_carries_the_donor_analysis_one_level_down(self):
        donor = {**GAGE, "peaks_table": True}  # for a heading inside the donor's analysis
        gage = build_section("gaged_ffa", dict(donor))
        transposed = build_section("transposed", {**TRANSPOSED, "donor": donor})
        headings = [(b.level, b.text) for b in transposed.blocks if isinstance(b, Heading)]
        assert [(b.level, b.text) for b in gage.blocks if isinstance(b, Heading)] == [
            (2, "Annual peaks")
        ]
        assert headings == [(2, "Caveats"), (2, DONOR_HEADING), (3, "Annual peaks")]
        assert len(transposed.blocks) == 6 + len(gage.blocks)
        captions = [b.caption for b in transposed.blocks if isinstance(b, Table) and b.caption]
        assert captions[0] == CAPTION
        assert "Flood quantiles with 90-percent confidence limits." in captions

    def test_without_detail_the_section_ends_with_its_caveats(self):
        content = _brief()
        assert [type(b) for b in content.blocks] == [Paragraph, Bullets, Table, Heading, Bullets]
        assert content.references == [BULLETIN_17C, CITATION]

    def test_the_options_as_built_rebuild_the_section(self, transposed):
        assert transposed.options["donor"] == build_section("gaged_ffa", dict(GAGE)).options
        assert transposed.options["donor_area_sqmi"] == 205.0
        assert transposed.options["exponents"] == EXPONENTS
        again = build_section("transposed", transposed.options)
        assert again.options == transposed.options
        assert _table(again.blocks) == _table(transposed.blocks)

    def test_the_callers_options_are_not_changed(self):
        options = {**TRANSPOSED, "donor": dict(GAGE), "detail": False}
        build_section("transposed", options)
        assert options["donor"] == GAGE and "donor_area_sqmi" not in options

    def test_one_exponent_for_every_aep(self):
        content = _brief(exponents={"exponent": 0.8, "citation": CITATION}, target_name=None)
        assert content.title == "Transposed flood frequency"
        assert "Flood quantiles at the ungaged site were transposed" in _text(content.blocks)
        rows = _table(content.blocks).rows
        assert {(row[2], row[-1]) for row in rows} == {("0.800", "constant")}
        by_aep = {row[0]: row for row in rows}
        assert by_aep["1"][3:5] == ("0.7789", "14,074")

    def test_the_donor_area_may_be_given_beside_the_donor(self):
        donor = {key: value for key, value in GAGE.items() if key != "drainage_area_sqmi"}
        with pytest.raises(ReportError, match="the donor has no drainage area; give donor_area"):
            _brief(donor=donor)
        content = _brief(donor=donor, donor_area_sqmi=205.0)
        assert _table(content.blocks) == _table(_brief().blocks)
        # Given both, the one beside the donor is the one used.
        other = _brief(donor_area_sqmi=180.0)
        assert "Area ratio: 0.8333" in _bullets(other.blocks)

    def test_max_aep_leaves_out_the_frequent_floods(self):
        rows = _table(_brief(max_aep=0.5).blocks).rows
        assert [row[0] for row in rows][:2] == ["50", "20"] and len(rows) == 8
        assert not any(row[-1] == "extrapolated" and float(row[0]) > 50 for row in rows)

    def test_an_area_ratio_outside_the_methods_range(self):
        with pytest.raises(ValueError, match="area ratio 0.2439 is outside the applicable range"):
            _brief(target_area_sqmi=50.0)
        content = _brief(target_area_sqmi=50.0, allow_out_of_range=True)
        assert any(
            b.startswith("**The area ratio 0.2439 is outside the applicable range (0.5, 1.5).**")
            for b in _bullets(content.blocks)
        )

    def test_a_regulated_donor_is_refused_unless_allowed(self):
        classes = {"donor_site_class": "regulated", "target_site_class": "reference"}
        with pytest.raises(RegulationMismatchError, match="donor is regulated"):
            _brief(**classes)
        allowed = _bullets(_brief(**classes, allow_regulation_mismatch=True).blocks)
        assert any(
            b.startswith("**Regulation mismatch, transferred by override**") for b in allowed
        )
        checked = _bullets(
            _brief(donor_site_class="reference", target_site_class="reference").blocks
        )
        assert not any("Regulation" in b for b in checked)

    def test_an_area_outside_the_regressions_range_is_a_caveat(self):
        exponents = {**EXPONENTS, "valid_area_range_sqmi": [200, 900]}
        content = _brief(exponents=exponents)
        caveats = content.blocks[-1]
        assert isinstance(caveats, Bullets)
        assert caveats.items[0] == (
            "**The site's drainage area, 150 square miles, is outside the range the cited "
            "regression was developed over** (200 to 900 square miles), so the exponent is "
            "applied beyond its basis."
        )
        assert not any("donor gage's drainage area" in item for item in caveats.items)
        inside = _brief(exponents={**EXPONENTS, "valid_area_range_sqmi": [10, 900]})
        assert not any("developed over" in item for item in _bullets(inside.blocks))

    @pytest.mark.parametrize(
        "change, message",
        [
            ({"donor": "03606500"}, "transposed: donor must be an object of gaged_ffa options"),
            ({"donor": {**GAGE, "regional_skw": 0.1}}, "gaged_ffa: unknown option.*regional_skw"),
            ({"target_area_sqmi": -5}, "target_area_sqmi must be a positive number, got -5"),
            ({"target_area_sqmi": True}, "target_area_sqmi must be a positive number, got True"),
            ({"target_area_sqmi": "small"}, "target_area_sqmi must be a positive number"),
            ({"donor_area_sqmi": 0}, "donor_area_sqmi must be a positive number, got 0"),
            ({"max_aep": "often"}, "max_aep must be a positive number, got 'often'"),
            ({"donor_site_class": "dammed"}, "donor_site_class must be one of reference, regul"),
            ({"target_site_class": 3}, "target_site_class must be one of reference, regulated"),
            ({"exponents": 0.8}, "exponents must be an object with a citation and either"),
            ({"exponents": {"exponent": 0.8}}, "exponents needs a citation naming the report"),
            ({"exponents": {"exponent": 0.8, "citation": " "}}, "exponents needs a citation"),
            ({"exponents": {"citation": CITATION}}, "not both and not neither"),
            ({"exponents": {**EXPONENTS, "exponent": 0.8}}, "not both and not neither"),
            (
                {"exponents": {"aeps": [0.5], "citation": CITATION}},
                "exponents by AEP needs both aeps and exponents",
            ),
            (
                {"exponents": {"exponent": "steep", "citation": CITATION}},
                "exponents exponent must be a positive number, got 'steep'",
            ),
            (
                {"exponents": {**EXPONENTS, "source": "a table"}},
                r"unknown exponents key\(s\) source; allowed: citation, exponent, aeps",
            ),
            (
                {"exponents": {**EXPONENTS, "valid_area_range_sqmi": 200}},
                "valid_area_range_sqmi must be the smallest and the largest drainage area",
            ),
            (
                {"exponents": {**EXPONENTS, "valid_area_range_sqmi": [0, 900]}},
                "exponents valid_area_range_sqmi must be a positive number, got 0",
            ),
        ],
    )
    def test_bad_options(self, change, message):
        with pytest.raises(ReportError, match=message):
            _brief(**change)

    def test_exponents_the_method_refuses_pass_through(self):
        exponents = {"aeps": [0.5, 0.1], "exponents": [0.7], "citation": CITATION}
        with pytest.raises(ValueError, match="aeps and exponents must be the same length"):
            _brief(exponents=exponents)

    def test_missing_part(self):
        with pytest.raises(ReportError, match="transposed: missing option.*exponents"):
            build_section("transposed", {"donor": dict(GAGE), "target_area_sqmi": 150.0})

    def test_in_a_report(self, tmp_path):
        spec = ReportSpec(
            "Big Sandy River",
            [SectionSpec("transposed", options=dict(TRANSPOSED))],
            date="2026-11-15",
        )
        report = build_report(spec)
        text = report.to_markdown()
        assert "## 1. Transposed flood frequency at Big Sandy River above Bruceton" in text
        assert "### Caveats" in text and f"### {DONOR_HEADING}" in text
        assert f"Table 1. {CAPTION}" in text
        assert "| 1 | 18,069 | 0.620 | 0.8239 | 14,888 | 11,539 | 22,339 | published |" in text
        assert text.count(CITATION) == 2  # the inputs, and the references
        path = tmp_path / "transposed.json"
        report.spec.save(path)
        assert build_report(ReportSpec.load(path)).to_markdown() == text

    def test_docx(self, tmp_path):
        docx = pytest.importorskip("docx")
        options = {**TRANSPOSED, "donor": {**GAGE, "peaks_table": True}}
        spec = ReportSpec("Big Sandy River", [SectionSpec("transposed", options=options)])
        path = build_report(spec).save(tmp_path / "transposed.docx")
        document = docx.Document(str(path))
        styles = {p.style.name for p in document.paragraphs}
        assert {"Heading 1", "Heading 2", "Heading 3"} <= styles
        cells = [c.text for t in document.tables for c in t.rows[0].cells]
        assert "Transposed flow (cfs)" in cells and "Exponent basis" in cells
