"""The `engine_check` report section: a gage's fit against a reference computation.

Nothing here needs the Fortran extension or makes a request. The reference is what
`emafitpr` returned for the fixture's record, stored as the section itself stores it, and
the live call is replaced wherever a test follows the path that makes one.
"""

from __future__ import annotations

import copy
import datetime as dt
from typing import Any, Dict, List

import pytest

from flowfreq.bulletin17c import Bulletin17C
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
from flowfreq.reporting.sections import BULLETIN_17C
from flowfreq.validation.reference import ReferenceResult
from tests.test_reporting_sites import FLOWS, GAGE

SOURCE = (
    "USGS peakfq 8.1.0, its Fortran routine emafitpr called on the same peaks, thresholds "
    "and regional skew by flowfreq 0.27.0.dev0"
)
#: `emafitpr` on the fixture's record (the vendored peakfq 8.1.0 through the f2py bridge,
#: 2026-11-15), as `engine_check` stores it, with the flows rounded to a thousandth of a cfs.
REFERENCE: Dict[str, Any] = {
    "source": SOURCE,
    "computed": "2026-11-15",
    "record": "00ec8592d3c6a928",
    "parameters": {
        "mean_log": 3.6909446183,
        "std_log": 0.267214073,
        "skew_at_site": -0.1874061748,
        "skew_weighted": -0.2809956404,
    },
    "low_outlier_count": 0,
    "low_outlier_threshold_cfs": 0.0,
    "quantiles": [
        {"aep": 0.995, "flow_cfs": 855.866, "lower_cfs": 399.375, "upper_cfs": 1236.947},
        {"aep": 0.99, "flow_cfs": 1034.51, "lower_cfs": 545.112, "upper_cfs": 1430.258},
        {"aep": 0.95, "flow_cfs": 1702.665, "lower_cfs": 1160.398, "upper_cfs": 2149.095},
        {"aep": 0.9, "flow_cfs": 2195.161, "lower_cfs": 1639.416, "upper_cfs": 2685.715},
        {"aep": 0.8, "flow_cfs": 2954.4, "lower_cfs": 2368.414, "upper_cfs": 3529.452},
        {"aep": 0.67, "flow_cfs": 3836.294, "lower_cfs": 3185.23, "upper_cfs": 4530.77},
        {"aep": 0.5, "flow_cfs": 5051.777, "lower_cfs": 4274.283, "upper_cfs": 5936.66},
        {"aep": 0.2, "flow_cfs": 8292.826, "lower_cfs": 7065.424, "upper_cfs": 9899.224},
        {"aep": 0.1, "flow_cfs": 10576.737, "lower_cfs": 8936.047, "upper_cfs": 13119.084},
        {"aep": 0.04, "flow_cfs": 13549.878, "lower_cfs": 11170.622, "upper_cfs": 18033.501},
        {"aep": 0.02, "flow_cfs": 15801.36, "lower_cfs": 12666.549, "upper_cfs": 22303.159},
        {"aep": 0.01, "flow_cfs": 18068.951, "lower_cfs": 14004.44, "upper_cfs": 27113.162},
        {"aep": 0.005, "flow_cfs": 20358.638, "lower_cfs": 15203.128, "upper_cfs": 32545.478},
        {"aep": 0.002, "flow_cfs": 23424.209, "lower_cfs": 16605.731, "upper_cfs": 40843.091},
    ],
}
CHECK: Dict[str, Any] = {"gage": dict(GAGE), "reference": REFERENCE}
QUANTILES = "Flood quantiles from flowfreq and from the reference."
PARAMETERS = "Fitted parameters from flowfreq and from the reference."
LIMITS = "90-percent confidence limits from flowfreq and from the reference."
GAGE_HEADING = "Analysis at the gage"
AGREE = "The two agree within the tolerances set."
DISAGREE = "**The two do not agree within the tolerances set.**"
ONE_PERCENT = ("1", "18,069", "18,069", "<0.001", "yes")


def _check(**changes: Any) -> Dict[str, Any]:
    """The check with a deep copy of its inputs, and `changes` laid over the options."""
    return {**copy.deepcopy(CHECK), **changes}


def _reference(**changes: Any) -> Dict[str, Any]:
    return {**copy.deepcopy(REFERENCE), **changes}


def _with_row(aep: float, **changes: Any) -> Dict[str, Any]:
    """The reference with one quantile row changed."""
    reference = copy.deepcopy(REFERENCE)
    (row,) = [r for r in reference["quantiles"] if r["aep"] == aep]
    row.update(changes)
    return reference


def _table(blocks: List[Any], caption: str) -> Table:
    (table,) = [b for b in blocks if isinstance(b, Table) and b.caption == caption]
    return table


def _captions(blocks: List[Any]) -> List[str]:
    return [b.caption for b in blocks if isinstance(b, Table)]


def _verdict(blocks: List[Any]) -> str:
    assert isinstance(blocks[0], Paragraph)
    return str(blocks[0].text)


def _as_emafitpr(reference: Dict[str, Any]) -> ReferenceResult:
    """What ``Bulletin17C.fortran_reference`` returns, from the numbers stored here."""
    rows = reference["quantiles"]
    return ReferenceResult(
        source="live emafitpr via flowfreq.peakfqr",
        low_outlier_count=reference["low_outlier_count"],
        low_outlier_threshold=reference["low_outlier_threshold_cfs"],
        parameters={**reference["parameters"], "mse_skew": 0.129},
        quantiles={r["aep"]: r["flow_cfs"] for r in rows},
        confidence_intervals={r["aep"]: (r["lower_cfs"], r["upper_cfs"]) for r in rows},
    )


@pytest.fixture
def no_extension(monkeypatch):
    """The Fortran extension is not built: asking for a live reference says so."""

    def absent(self, *args, **kwargs):
        raise ImportError("flowfreq.peakfqr is not built")

    monkeypatch.setattr(Bulletin17C, "fortran_reference", absent)


@pytest.fixture
def emafitpr(monkeypatch):
    """The live call, answered from the stored numbers; the list counts the calls."""
    calls: List[int] = []

    def live(self, *args, **kwargs):
        calls.append(1)
        return _as_emafitpr(REFERENCE)

    monkeypatch.setattr(Bulletin17C, "fortran_reference", live)
    return calls


class TestStoredReference:
    """A reference in the options is compared without the extension."""

    @pytest.fixture
    def content(self, no_extension):
        return build_section("engine_check", _check())

    def test_kind(self):
        kind = SECTION_KINDS["engine_check"]
        assert kind.required == ("gage",)
        assert set(kind.optional) == {
            "reference",
            "tolerance_pct",
            "parameter_tolerance_pct",
            "ci_tolerance_pct",
            "detail",
        }

    def test_title_and_verdict(self, content):
        assert content.title == (
            "Computation check at USGS 03606500 Big Sandy River at Bruceton, TN"
        )
        verdict = _verdict(content.blocks)
        assert verdict.startswith(
            "The Bulletin 17C fit at USGS 03606500 Big Sandy River at Bruceton, TN was "
            f"checked against a reference computation of the same record: {SOURCE} "
            f"(2026-11-15). {AGREE} The largest differences are <0.001 percent among the "
            "14 flood quantiles (tolerance 1 percent); 0.02"
        )
        assert " percent among their 90-percent confidence limits (tolerance 2 percent); " in (
            verdict
        )
        assert verdict.endswith(
            "<0.001 percent in the mean and standard deviation of log10 flow (tolerance 0.5 "
            "percent); 0.0000 in the skew coefficients (tolerance 0.05, in skew units). "
            "Neither censors a value as a low outlier."
        )
        assert "Outside tolerance" not in verdict

    def test_quantile_table(self, content):
        table = _table(content.blocks, QUANTILES)
        assert table.header == (
            "Annual exceedance probability (percent)",
            "flowfreq (cfs)",
            "Reference (cfs)",
            "Difference (percent)",
            "Within tolerance",
        )
        assert [row[0] for row in table.rows] == [
            "99.5",
            "99",
            "95",
            "90",
            "80",
            "67",
            "50",
            "20",
            "10",
            "4",
            "2",
            "1",
            "0.5",
            "0.2",
        ]
        assert ONE_PERCENT in table.rows
        assert {row[-1] for row in table.rows} == {"yes"}

    def test_parameter_table(self, content):
        table = _table(content.blocks, PARAMETERS)
        assert table.header == (
            "Parameter",
            "flowfreq",
            "Reference",
            "Difference",
            "Tolerance",
            "Within tolerance",
        )
        assert table.rows == (
            ("Mean of log10 flow", "3.6909", "3.6909", "<0.001 percent", "0.5 percent", "yes"),
            (
                "Standard deviation of log10 flow",
                "0.2672",
                "0.2672",
                "<0.001 percent",
                "0.5 percent",
                "yes",
            ),
            ("Station skew", "-0.1874", "-0.1874", "0.0000", "0.05", "yes"),
            ("Skew used", "-0.2810", "-0.2810", "0.0000", "0.05", "yes"),
            ("Values censored as low outliers", "0", "0", "0", "0", "yes"),
            ("Low-outlier threshold (cfs)", "none", "none", "", "", "not judged"),
        )

    def test_confidence_limit_table(self, content):
        table = _table(content.blocks, LIMITS)
        assert len(table.rows) == 14 and len(table.header) == 7
        assert ("10", "8,937", "13,117", "8,936", "13,119", "0.014", "yes") in table.rows

    def test_the_gage_follows_unless_detail_is_off(self, content, no_extension):
        headings = [b.text for b in content.blocks if isinstance(b, Heading)]
        assert headings == [GAGE_HEADING]
        after = content.blocks[content.blocks.index(Heading(2, GAGE_HEADING)) + 1 :]
        assert "Flood quantiles with 90-percent confidence limits." in _captions(after)
        brief = build_section("engine_check", _check(detail=False))
        assert _captions(brief.blocks) == [PARAMETERS, QUANTILES, LIMITS]
        assert not [b for b in brief.blocks if isinstance(b, Heading)]

    def test_options_and_references(self, content):
        assert content.options["reference"] == REFERENCE
        assert content.options["gage"]["peaks"]["peak_flow_cfs"] == FLOWS
        assert content.references[0] == BULLETIN_17C

    def test_in_a_report(self, no_extension, tmp_path):
        spec = ReportSpec(
            "Big Sandy River",
            [SectionSpec("engine_check", options=_check())],
            date="2026-11-15",
        )
        report = build_report(spec)
        text = report.to_markdown()
        assert "## 1. Computation check at USGS 03606500 Big Sandy River at Bruceton, TN" in text
        assert f"### {GAGE_HEADING}" in text
        assert f"Table 2. {QUANTILES}" in text
        assert "| 1 | 18,069 | 18,069 | <0.001 | yes |" in text
        assert AGREE in text
        assert report.spec.sections[0].options["reference"] == REFERENCE
        # The saved specification carries the reference, so the report rebuilds without one.
        path = tmp_path / "check.json"
        report.spec.save(path)
        assert build_report(ReportSpec.load(path)).to_markdown() == text

    def test_docx(self, no_extension, tmp_path):
        docx = pytest.importorskip("docx")
        spec = ReportSpec("Big Sandy River", [SectionSpec("engine_check", options=_check())])
        path = build_report(spec).save(tmp_path / "check.docx")
        document = docx.Document(str(path))
        cells = [[c.text for c in row.cells] for table in document.tables for row in table.rows]
        assert list(ONE_PERCENT) in cells
        assert any(AGREE in p.text for p in document.paragraphs)


class TestLiveReference:
    """Without a reference the record goes through peakfq, and the answer is stored."""

    def test_the_answer_is_stored(self, emafitpr):
        content = build_section("engine_check", {"gage": dict(GAGE)})
        assert emafitpr == [1]
        stored = content.options["reference"]
        assert stored["source"].startswith("USGS peakfq 8.1.0, its Fortran routine emafitpr")
        assert stored["computed"] == dt.date.today().isoformat()
        assert stored["record"] == REFERENCE["record"]
        # Only the four parameters the section prints are kept.
        assert stored["parameters"] == REFERENCE["parameters"]
        assert stored["quantiles"] == REFERENCE["quantiles"]
        assert stored["low_outlier_count"] == 0 and stored["low_outlier_threshold_cfs"] == 0.0
        assert AGREE in _verdict(content.blocks)
        assert ONE_PERCENT in _table(content.blocks, QUANTILES).rows

    def test_it_rebuilds_from_what_it_stored(self, monkeypatch, emafitpr):
        first = build_section("engine_check", {"gage": dict(GAGE), "detail": False})

        def absent(self, *args, **kwargs):
            raise ImportError("not built")

        monkeypatch.setattr(Bulletin17C, "fortran_reference", absent)
        again = build_section("engine_check", first.options)
        assert again.blocks == first.blocks and again.options == first.options

    def test_without_the_extension_it_says_how_to_get_one(self, no_extension):
        with pytest.raises(ReportError, match="which is not built in this environment") as error:
            build_section("engine_check", {"gage": dict(GAGE)})
        assert "python build_fortran/build.py" in str(error.value)
        assert "give the numbers of a reference computation in `reference`" in str(error.value)

    def test_a_stored_reference_belongs_to_its_record(self, no_extension):
        flows = list(FLOWS)
        flows[0] += 1.0
        other = {**GAGE, "peaks": {**GAGE["peaks"], "peak_flow_cfs": flows}}
        with pytest.raises(ReportError, match="computed for other peaks, thresholds or"):
            build_section("engine_check", {"gage": other, "reference": REFERENCE})
        # The skew decides the fit as well.
        other = {**GAGE, "regional_skew": GAGE["regional_skew"] + 0.1}
        with pytest.raises(ReportError, match="Remove reference to compute it again"):
            build_section("engine_check", {"gage": other, "reference": REFERENCE})
        # What does not decide it is free to change.
        renamed = {**GAGE, "site_name": "Big Sandy River", "peaks_table": True}
        build_section("engine_check", {"gage": renamed, "reference": REFERENCE})


class TestDisagreement:
    def test_a_quantile_outside_tolerance(self, no_extension):
        reference = _with_row(0.01, flow_cfs=19000.0)
        content = build_section("engine_check", _check(reference=reference, detail=False))
        verdict = _verdict(content.blocks)
        assert DISAGREE in verdict and AGREE not in verdict
        assert "4.900 percent among the 14 flood quantiles (tolerance 1 percent)" in verdict
        assert verdict.endswith("Outside tolerance: the flood quantiles.")
        rows = {row[0]: row for row in _table(content.blocks, QUANTILES).rows}
        assert rows["1"] == ("1", "18,069", "19,000", "4.900", "no")
        assert rows["2"][-1] == "yes"

    def test_a_wider_tolerance_accepts_it(self, no_extension):
        reference = _with_row(0.01, flow_cfs=19000.0)
        content = build_section("engine_check", _check(reference=reference, tolerance_pct=5))
        verdict = _verdict(content.blocks)
        assert AGREE in verdict and "(tolerance 5 percent)" in verdict
        assert content.options["tolerance_pct"] == 5

    def test_every_category_is_named(self, no_extension):
        reference = _with_row(0.01, flow_cfs=19000.0, lower_cfs=12000.0)
        reference["parameters"].update(std_log=0.30, skew_weighted=-0.10)
        reference["low_outlier_count"] = 2
        reference["low_outlier_threshold_cfs"] = 1500.0
        content = build_section("engine_check", _check(reference=reference, detail=False))
        verdict = _verdict(content.blocks)
        assert verdict.endswith(
            "Outside tolerance: the flood quantiles; the confidence limits; the mean or "
            "standard deviation; the skew; the number of values censored as low outliers "
            "(flowfreq 0, the reference 2)."
        )
        rows = {row[0]: row for row in _table(content.blocks, PARAMETERS).rows}
        assert rows["Mean of log10 flow"][-1] == "yes"
        assert rows["Standard deviation of log10 flow"][3:] == (
            "10.929 percent",
            "0.5 percent",
            "no",
        )
        assert rows["Station skew"][-1] == "yes"
        assert rows["Skew used"] == ("Skew used", "-0.2810", "-0.1000", "0.1810", "0.05", "no")
        assert rows["Values censored as low outliers"][1:] == ("0", "2", "2", "0", "no")
        assert rows["Low-outlier threshold (cfs)"][1:3] == ("none", "1,500")
        limits = {row[0]: row for row in _table(content.blocks, LIMITS).rows}
        assert limits["1"][-1] == "no" and limits["2"][-1] == "yes"

    def test_a_low_outlier_count_alone_is_a_disagreement(self, no_extension):
        content = build_section("engine_check", _check(reference=_reference(low_outlier_count=1)))
        verdict = _verdict(content.blocks)
        assert DISAGREE in verdict
        assert verdict.endswith(
            "Outside tolerance: the number of values censored as low outliers (flowfreq 0, "
            "the reference 1)."
        )

    def test_the_other_tolerances(self, no_extension):
        reference = _with_row(0.01, lower_cfs=13000.0)
        reference["parameters"]["std_log"] = 0.27
        tight = build_section("engine_check", _check(reference=reference))
        assert DISAGREE in _verdict(tight.blocks)
        loose = build_section(
            "engine_check",
            _check(reference=reference, ci_tolerance_pct=10, parameter_tolerance_pct=2.0),
        )
        verdict = _verdict(loose.blocks)
        assert AGREE in verdict
        assert "(tolerance 10 percent)" in verdict and "(tolerance 2 percent)" in verdict


class TestAnalystReference:
    """Numbers typed in from another program's output: only what was given is compared."""

    PRINTOUT: Dict[str, Any] = {
        "source": "A second program's printout for this gage, table of computed quantiles",
        "quantiles": [
            {"aep": 0.01, "flow_cfs": 18100.0},
            {"aep": 0.1, "flow_cfs": 10600.0},
            {"aep": 0.03, "flow_cfs": 14500.0},
        ],
    }

    def test_quantiles_only(self, no_extension):
        content = build_section("engine_check", {"gage": dict(GAGE), "reference": self.PRINTOUT})
        assert _captions(content.blocks)[:1] == [QUANTILES]
        assert PARAMETERS not in _captions(content.blocks)
        assert LIMITS not in _captions(content.blocks)
        table = _table(content.blocks, QUANTILES)
        # In order of AEP, and at an AEP the gage's own table does not have.
        assert [row[0] for row in table.rows] == ["10", "3", "1"]
        assert table.rows[2] == ("1", "18,069", "18,100", "0.172", "yes")
        assert table.rows[1][2] == "14,500"
        verdict = _verdict(content.blocks)
        assert verdict.startswith(
            "The Bulletin 17C fit at USGS 03606500 Big Sandy River at Bruceton, TN was "
            "checked against a reference computation of the same record: A second program's "
            f"printout for this gage, table of computed quantiles. {AGREE} The largest "
            "differences are "
        )
        assert verdict.endswith(" percent among the 3 flood quantiles (tolerance 1 percent).")
        assert "low outlier" not in verdict
        assert content.options["reference"] == self.PRINTOUT

    def test_limits_at_some_aeps_and_one_parameter(self, no_extension):
        reference = copy.deepcopy(self.PRINTOUT)
        reference["quantiles"][0].update(lower_cfs=14000.0, upper_cfs=27100.0)
        reference["parameters"] = {"skew_weighted": -0.28}
        content = build_section("engine_check", {"gage": dict(GAGE), "reference": reference})
        assert [row[0] for row in _table(content.blocks, LIMITS).rows] == ["1"]
        assert _table(content.blocks, PARAMETERS).rows == (
            ("Skew used", "-0.2810", "-0.2800", "0.0010", "0.05", "yes"),
        )

    def test_a_gage_without_a_name(self, no_extension):
        gage = {k: v for k, v in GAGE.items() if k not in ("site_no", "site_name")}
        content = build_section("engine_check", {"gage": gage, "reference": self.PRINTOUT})
        assert content.title == "Computation check"
        assert _verdict(content.blocks).startswith("The Bulletin 17C fit was checked against")


class TestRefusals:
    @pytest.mark.parametrize(
        "options, match",
        [
            ({"gage": "03606500"}, "engine_check: gage must be an object"),
            ({"gage": {**GAGE, "reference": {}}}, "unknown option.*reference"),
            ({"reference": "peakfq"}, "engine_check: reference must be an object"),
            (
                {"reference": _reference(station="03606500")},
                "unknown reference key.*station; allowed: source, computed",
            ),
            ({"reference": _reference(source=" ")}, "reference needs a source"),
            ({"reference": _reference(source=8.1)}, "reference needs a source"),
            ({"reference": _reference(parameters=[3.69])}, "reference parameters must be an"),
            (
                {"reference": _reference(parameters={"skew": -0.2})},
                "unknown reference parameter.*skew; allowed: mean_log, std_log, skew_at_site",
            ),
            (
                {"reference": _reference(parameters={"mean_log": "high"})},
                "the reference mean_log must be a number, got 'high'",
            ),
            ({"reference": _reference(quantiles=[])}, "reference needs quantiles"),
            ({"reference": _reference(quantiles={"0.01": 18069})}, "reference needs quantiles"),
            (
                {"reference": _reference(quantiles=[{"aep": 0.01}])},
                "each reference quantile needs aep and flow_cfs",
            ),
            (
                {"reference": _reference(quantiles=[{"aep": 0.01, "flow_cfs": 1.0, "ci": 2}])},
                "unknown reference quantile key.*ci; allowed: aep, flow_cfs",
            ),
            (
                {"reference": _reference(quantiles=[{"aep": 1, "flow_cfs": 18069.0}])},
                "a reference aep is a fraction between 0 and 1, got 1",
            ),
            (
                {"reference": _reference(quantiles=[{"aep": 0, "flow_cfs": 18069.0}])},
                "a reference aep must be a positive number, got 0",
            ),
            (
                {
                    "reference": _reference(
                        quantiles=[{"aep": 0.01, "flow_cfs": 1.0}, {"aep": 0.01, "flow_cfs": 2.0}]
                    )
                },
                "the reference gives AEP 0.01 more than once",
            ),
            (
                {"reference": _reference(quantiles=[{"aep": 0.01, "flow_cfs": -5.0}])},
                "the reference flow at AEP 0.01 must be a positive number, got -5.0",
            ),
            (
                {"reference": _with_row(0.01, upper_cfs=None)},
                "the reference at AEP 0.01 needs lower_cfs and upper_cfs together",
            ),
            (
                {"reference": _with_row(0.01, lower_cfs="low")},
                "the reference lower limit at AEP 0.01 must be a positive number, got 'low'",
            ),
            (
                {"reference": _reference(low_outlier_count=1.5)},
                "reference low_outlier_count must be a whole number, got 1.5",
            ),
            (
                {"reference": _reference(low_outlier_count=True)},
                "reference low_outlier_count must be a whole number, got True",
            ),
            (
                {"reference": _reference(low_outlier_threshold_cfs=-1.0)},
                "reference low_outlier_threshold_cfs must not be negative, got -1.0",
            ),
            ({"tolerance_pct": 0}, "engine_check: tolerance_pct must be a positive number"),
            ({"ci_tolerance_pct": "wide"}, "ci_tolerance_pct must be a positive number"),
            ({"parameter_tolerance_pct": -1}, "parameter_tolerance_pct must be a positive"),
            ({"tolerance": 5}, "unknown option.*tolerance"),
        ],
    )
    def test_refused(self, no_extension, options, match):
        with pytest.raises(ReportError, match=match):
            build_section("engine_check", _check(**options))

    def test_the_gage_is_required(self):
        with pytest.raises(ReportError, match="gage"):
            build_section("engine_check", {"reference": REFERENCE})
