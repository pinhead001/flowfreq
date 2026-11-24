"""The ``b17b_ffa`` report section: Bulletin 17B (superseded), computed only when asked for.

Everything here is offline: the peaks are inline, and the tests of a download replace
``gaged.download_peaks``.
"""

from __future__ import annotations

import copy
import json
import math
import warnings
from typing import Any, Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pytest

from flowfreq.bulletin17b import STANDARD_AEP, Bulletin17B, Bulletin17BError
from flowfreq.regional_skew import regional_skew_for
from flowfreq.reporting import (
    SECTION_KINDS,
    ReportError,
    ReportSpec,
    SectionSpec,
    build_report,
    build_section,
    gaged,
)
from flowfreq.reporting._base import section_title
from flowfreq.reporting.b17b import (
    B17B_CAPTION,
    BULLETIN_17B,
    DEFAULT_AEP,
    _b17b_ffa,
    _build_b17b_ffa,
    _frequency_figure,
)
from flowfreq.reporting.blocks import Bullets, Figure, Heading, Paragraph, Table
from flowfreq.reporting.sections import cfs
from flowfreq.reporting.session import MAIN_TABLES, ReportSession
from tests.fixtures.b17b import record

LABEL = "Bulletin 17B (superseded)"
REFUSED = "Bulletin 17B gives no answer for this record: "
#: A gage the GAGES-II screen classes regulated.
REGULATED_SITE = "01019000"


def peaks_spec(name: str) -> Dict[str, Any]:
    """An example's record as a spec stores it, its historic peaks under code 7."""
    rec = record(name)
    return {
        "water_year": list(rec.water_years) + [y for y, _ in rec.historic_peaks],
        "peak_flow_cfs": list(rec.peaks_cfs) + [q for _, q in rec.historic_peaks],
        "qualification_code": [None] * len(rec.water_years) + ["7"] * len(rec.historic_peaks),
    }


def with_peaks(extra: List[Tuple[int, float, Any]], name: str = "example1") -> Dict[str, Any]:
    """An example's peaks with further (year, flow, code) rows."""
    peaks = peaks_spec(name)
    for year, flow, code in extra:
        peaks["water_year"].append(year)
        peaks["peak_flow_cfs"].append(flow)
        peaks["qualification_code"].append(code)
    return peaks


def text_of(block: Any) -> str:
    if isinstance(block, Paragraph):
        return block.text
    if isinstance(block, Bullets):
        return "\n".join(block.items)
    return ""


def all_text(content: Any) -> str:
    return "\n".join(text_of(b) for b in content.blocks)


def tables(content: Any) -> Dict[Any, Table]:
    return {b.caption: b for b in content.blocks if isinstance(b, Table)}


def parameters(content: Any) -> Dict[str, str]:
    (caption,) = [c for c in tables(content) if c and "parameters" in c]
    return dict(tables(content)[caption].rows)


def figure_of(content: Any) -> Figure:
    (figure,) = [b for b in content.blocks if isinstance(b, Figure)]
    return figure


EX1: Dict[str, Any] = {"peaks": peaks_spec("example1"), "generalized_skew": 0.6}
SANDY: Dict[str, Any] = {"peaks": peaks_spec("appendix6_big_sandy"), "generalized_skew": -0.2}
SANDY_ROWS = [{"start_year": 1897, "end_year": 1929, "threshold_cfs": 18000.0}]


# -- labels -------------------------------------------------------------------------------


def test_section_builds_example1_with_superseded_labels():
    content = build_section("b17b_ffa", {"peaks": peaks_spec("example1"), "generalized_skew": 0.6})
    assert "Bulletin 17B (superseded)" in content.title
    assert any("superseded" in b.caption for b in content.blocks if getattr(b, "caption", None))


def test_every_caption_with_a_b17b_number_says_superseded():
    content = build_section("b17b_ffa", dict(EX1))
    captions = [b.caption for b in content.blocks if getattr(b, "caption", None)]
    assert len(captions) == 4
    assert [c for c in captions if LABEL not in c] == ["Station and peak-flow record."]
    assert B17B_CAPTION in captions and B17B_CAPTION.startswith(LABEL + " quantiles")
    assert content.title == "Flood frequency, " + LABEL
    named = build_section("b17b_ffa", {**EX1, "site_no": "01373500", "site_name": "Fishkill Creek"})
    assert named.title == "Flood frequency at USGS 01373500 Fishkill Creek, " + LABEL


def test_the_section_says_17b_is_superseded_and_why_it_is_there():
    said = all_text(build_section("b17b_ffa", dict(EX1)))
    assert "superseded by Bulletin 17C" in said and "comparison" in said


def _headings(report: Any) -> List[str]:
    return [b.text for b in report.blocks if isinstance(b, Heading) and b.level == 1]


def test_a_given_title_keeps_the_label_in_a_built_report():
    given = "Flood frequency at the gage"
    spec = ReportSpec(
        title="Report",
        sections=[
            SectionSpec("b17b_ffa", title=given, options=dict(EX1)),
            SectionSpec("b17b_ffa", title=f"Old study, {LABEL}", options=dict(EX1)),
            SectionSpec(
                "b17b_ffa", title="Old study, bulletin 17b (SUPERSEDED)", options=dict(EX1)
            ),
            SectionSpec("text", title=given, narrative="Words."),
        ],
    )
    assert _headings(build_report(spec))[:4] == [
        f"1. {given}, {LABEL}",
        f"2. Old study, {LABEL}",
        "3. Old study, bulletin 17b (SUPERSEDED)",
        f"4. {given}",
    ]


def test_a_given_title_keeps_the_label_in_a_session(tmp_path):
    session = ReportSession("Report", output_dir=tmp_path)
    added = session.add_section("b17b_ffa", dict(EX1), title="At the gage")
    assert added["title"] == f"At the gage, {LABEL}"
    assert session.add_text("At the gage", "Words.")["title"] == "At the gage"
    refused = session.add_section(
        "b17b_ffa", {**EX1, "low_outlier_threshold": 2190.0}, title="Refused"
    )
    assert refused["title"] == f"Refused, {LABEL}"
    assert f"At the gage, {LABEL}" in json.dumps(session.outline())
    session.render("report.md")
    assert f"1. At the gage, {LABEL}" in (tmp_path / "report.md").read_text(encoding="utf-8")


def test_only_the_b17b_kind_has_its_titles_labelled():
    for kind in SECTION_KINDS:
        expected = f"Given, {LABEL}" if kind == "b17b_ffa" else "Given"
        assert section_title(kind, "Given", "Own") == expected
        assert section_title(kind, None, f"Own, {LABEL}") == f"Own, {LABEL}"
        assert section_title(kind, "", "") == ""


# -- skew ---------------------------------------------------------------------------------


def test_section_uses_the_b17c_regional_skew_by_default():
    content, res = _b17b_ffa(
        {"peaks": peaks_spec("example1"), "regional_skew": 0.6, "regional_skew_se": 0.55}
    )
    assert res.skew_generalized == 0.6 and res.skew_generalized_mse == pytest.approx(0.3025)
    assert any("change in procedure" in text_of(b) for b in content.blocks)
    assert (
        "the regional skew a Bulletin 17C analysis with the same options uses, so the "
        "difference is the change in procedure" in all_text(content)
    )


def test_generalized_skew_overrides_the_regional_skew():
    content, res = _b17b_ffa(
        {
            "peaks": peaks_spec("example1"),
            "regional_skew": 0.1,
            "regional_skew_se": 0.55,
            "generalized_skew": 0.6,
        }
    )
    assert res.skew_generalized == 0.6 and res.skew_generalized_mse == 0.302
    said = all_text(content)
    assert "generalized skew from the earlier study" in said
    assert "in place of the regional skew" in said
    assert "change in procedure" not in said
    content, res = _b17b_ffa({**EX1, "generalized_skew_mse": 0.2})
    assert res.skew_generalized_mse == 0.2
    assert "in place of" not in all_text(content)


def test_a_generalized_skew_needs_no_lookup(monkeypatch):
    monkeypatch.setattr(
        "flowfreq.regional_skew.regional_skew_for", lambda *a, **k: pytest.fail("looked up")
    )
    content, res = _b17b_ffa({**EX1, "regional_skew_state": "CA"})
    assert res.skew_generalized == 0.6
    assert build_section("b17b_ffa", content.options).options == content.options


def test_no_skew_option_uses_station_skew_and_says_so():
    content, res = _b17b_ffa({"peaks": peaks_spec("example1")})
    assert res.skew_generalized is None and res.skew_used == res.skew_station
    assert "station skew was used alone" in all_text(content)
    assert "change in procedure" not in all_text(content)
    assert content.references == [BULLETIN_17B]


def test_skew_looked_up_by_state_is_frozen_and_cited():
    found = regional_skew_for("TN")
    content, res = _b17b_ffa({"peaks": peaks_spec("example1"), "regional_skew_state": "TN"})
    options = content.options
    assert "regional_skew_state" not in options
    assert options["regional_skew"] == pytest.approx(0.048)
    assert options["regional_skew_se"] == pytest.approx(math.sqrt(0.16))
    assert res.skew_generalized_mse == pytest.approx(0.16)
    assert content.references == [BULLETIN_17B, found.citation]
    assert "Messinger and others (2026): Skew Region 2" in all_text(content)
    again = build_section("b17b_ffa", options)
    assert again.options == options
    assert again.blocks == content.blocks


def test_a_lookup_without_a_published_skew_raises():
    with pytest.raises(ReportError, match="b17b_ffa: no regional skew: CA has 2 skew regions"):
        build_section("b17b_ffa", {"peaks": peaks_spec("example1"), "regional_skew_state": "CA"})


# -- bad options --------------------------------------------------------------------------


def _row(**change: Any) -> Dict[str, Any]:
    return {
        "perception_thresholds": [
            {"start_year": 1890, "end_year": 1944, "threshold_cfs": 9000.0, **change}
        ]
    }


@pytest.mark.parametrize(
    "extra, message",
    [
        ({"regional_skew": 0.1}, "regional_skew and regional_skew_se go together"),
        ({"regional_skew_se": 0.1}, "regional_skew and regional_skew_se go together"),
        ({"regional_skew": 0.1, "regional_skew_se": 0.0}, "regional_skew_se must be positive"),
        ({"regional_skew_state": "TN", "regional_skew_from_site": True}, "not both"),
        ({"regional_skew_region": "Skew Region 2"}, "needs regional_skew_state"),
        ({"regional_skew_from_site": True}, "needs site_no"),
        ({"generalized_skew": "0.6x"}, "generalized_skew must be a number"),
        ({"generalized_skew": True}, "generalized_skew must be a number"),
        ({"generalized_skew": math.inf}, "generalized_skew must be a number"),
        ({"generalized_skew": 0.6, "generalized_skew_mse": -1.0}, "must be positive"),
        ({"generalized_skew_mse": 0.2}, "generalized_skew_mse needs generalized_skew"),
        ({"low_outlier_threshold": -1.0}, "low_outlier_threshold must not be negative"),
        ({"low_outlier_threshold": "low"}, "low_outlier_threshold must be a number"),
        ({"station_skew_only": True, "generalized_skew": 0.6}, "station_skew_only"),
        ({"station_skew_only": True, "generalized_skew_mse": 0.3}, "station_skew_only"),
        ({"station_skew_only": True, "regional_skew": 0.1}, "station_skew_only"),
        ({"station_skew_only": True, "regional_skew_se": 0.1}, "station_skew_only"),
        ({"station_skew_only": True, "regional_skew_state": "TN"}, "station_skew_only"),
        ({"station_skew_only": True, "regional_skew_region": "x"}, "station_skew_only"),
        ({"station_skew_only": True, "regional_skew_from_site": True}, "station_skew_only"),
        ({"aeps": [0.01, 1e-12]}, "aeps must be between 1e-06 and 0.999999"),
        ({"aeps": [1 - 1e-12]}, "aeps must be between 1e-06 and 0.999999"),
        (_row(start_year="long ago"), "start_year must be a year"),
        ({"historic_period": [1900]}, r"historic_period must be \[first year, last year\]"),
        ({"historic_period": 1900}, r"historic_period must be \[first year, last year\]"),
        ({"historic_period": [1900.5, 1968]}, r"historic_period must be \[first year"),
        ({"historic_period": [1968, 1900]}, "historic_period ends before it begins"),
        ({"aeps": []}, "aeps must be a non-empty list"),
        ({"aeps": 0.01}, "aeps must be a non-empty list"),
        ({"aeps": ["x"]}, "aeps must be"),
        ({"aeps": [0.01, 1.5]}, "aeps must be"),
        ({"aeps": [0.67, 0.667]}, "same AEP given twice"),
        ({"perception_thresholds": [{"start_year": 1890}]}, "perception threshold needs"),
        (_row(threshold_cfs="much"), "a threshold must be a number"),
        (_row(end_year=1880), "ends before it begins"),
        ({"drainage_area_sqmi": "big"}, "drainage_area_sqmi must be a number"),
        ({"regional_skw": 0.1}, "unknown option"),
        ({"use_default_skew": True}, "unknown option"),
        ({"stationarity": False}, "unknown option"),
    ],
)
def test_bad_options_raise(extra, message):
    with pytest.raises(ReportError, match="b17b_ffa: .*" + message):
        build_section("b17b_ffa", {"peaks": peaks_spec("example1"), **extra})


def test_bad_peaks_raise():
    with pytest.raises(ReportError, match="b17b_ffa: give site_no"):
        build_section("b17b_ffa", {})
    with pytest.raises(ReportError, match="b17b_ffa: peak_flow_cfs holds a non-finite value"):
        build_section("b17b_ffa", {"peaks": with_peaks([(1980, math.nan, None)])})


@pytest.mark.parametrize("code", [{"a": 1}, ["7", "2"], True, 7.5, math.inf])
@pytest.mark.parametrize("apply", [True, False])
def test_an_unreadable_code_raises(code, apply):
    peaks = with_peaks([(1975, 2000.0, code)])
    with pytest.raises(ReportError, match="b17b_ffa: the qualification code of 1975"):
        build_section("b17b_ffa", {**EX1, "peaks": peaks, "apply_peak_codes": apply})


@pytest.mark.parametrize("code", [3, 3.0, "3"])
def test_a_code_may_be_written_as_a_number(code):
    peaks = with_peaks([(1975, 90000.0, code)])
    content, res = _b17b_ffa({**EX1, "peaks": peaks, "site_no": "99999999"})
    assert res.n_years == 24 and "1975 (90,000 cfs, code 3)" in all_text(content)


@pytest.mark.parametrize("zero", [0, 0.0])
def test_a_low_outlier_threshold_of_zero_means_none(zero):
    content, res = _b17b_ffa({**EX1, "low_outlier_threshold": zero})
    _, plain = _b17b_ffa(dict(EX1))
    assert res.low_outlier_threshold == plain.low_outlier_threshold
    assert "replaces eq. 8a" not in "\n".join(res.steps)
    assert build_section("b17b_ffa", content.options).blocks == content.blocks


def test_station_skew_only_uses_the_station_skew_and_says_so():
    content, res = _b17b_ffa({"peaks": peaks_spec("example1"), "station_skew_only": True})
    assert res.skew_generalized is None and res.skew_used == res.skew_station
    assert "station skew was used alone, as asked" in all_text(content)
    assert content.options["station_skew_only"] is True
    _, res = _b17b_ffa({**EX1, "station_skew_only": False})
    assert res.skew_generalized == 0.6


def test_threshold_years_may_be_text():
    rows = [{"start_year": "1897", "end_year": "1929", "threshold_cfs": "18000"}]
    content, res = _b17b_ffa({**SANDY, "perception_thresholds": rows})
    assert res.historic.h == 77 and "1897 to 1929, 18,000 cfs" in all_text(content)


@pytest.mark.parametrize("name", ["example1", "example4"])
def test_extreme_aeps_do_not_overflow(name):
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        content, res = _b17b_ffa(
            {"peaks": peaks_spec(name), "generalized_skew": 0.6, "aeps": [1 - 1e-6, 0.01, 1e-6]}
        )
    assert np.all(np.isfinite(res.quantiles[["flow_cfs", "lower", "upper"]].to_numpy()))
    assert _scatters(res)["Annual peaks"] > 20


@pytest.mark.parametrize("aep", [1e-7, 1e-12, 1 - 1e-7, 1 - 1e-12])
def test_the_class_refuses_aeps_it_cannot_evaluate(aep):
    rec = record("example1")
    with pytest.raises(Bulletin17BError, match="1e-06 and 0.999999"):
        Bulletin17B(rec.peaks_cfs, rec.water_years, aeps=[aep])


# -- refusals -----------------------------------------------------------------------------


def test_refusal_renders_a_sentence_not_an_error():
    content, res = _b17b_ffa({"peaks": peaks_spec("example1"), "low_outlier_threshold": 2190.0})
    assert res is None and "gives no answer" in text_of(content.blocks[-1])
    said = text_of(content.blocks[-1])
    assert said.startswith(REFUSED) and "of the 24 years are truncated" in said
    assert said.endswith(".") and not said.endswith("..")
    assert LABEL in content.title
    assert B17B_CAPTION not in tables(content)
    assert not any(isinstance(b, Figure) for b in content.blocks)
    assert content.references == [BULLETIN_17B]
    again = build_section("b17b_ffa", json.loads(json.dumps(content.options)))
    assert again.blocks == content.blocks


def test_a_regulated_record_is_a_sentence_and_the_override_is_stated():
    given = {**EX1, "site_no": REGULATED_SITE}
    content, res = _b17b_ffa(given)
    assert res is None
    said = text_of(content.blocks[-1])
    assert said.startswith(
        REFUSED + "USGS " + REGULATED_SITE + " is classed regulated (dam storage"
    )
    assert "17C" not in said and "allow_regulated=True" not in said and "England" not in said
    assert "Bulletin 17B" in said and "The section's allow_regulated option overrides" in said
    assert content.references == [BULLETIN_17B]
    content, res = _b17b_ffa({**given, "allow_regulated": True})
    assert res is not None
    assert "classified as regulated" in all_text(content)
    rows = dict(tables(content)["Station and peak-flow record."].rows)
    assert "regulated" in rows["Regulation screen"]


def test_a_code_6_peak_at_a_numbered_gage_is_the_regulated_sentence():
    peaks = with_peaks([(1975, 2000.0, "6")])
    content, res = _b17b_ffa({**EX1, "peaks": peaks, "site_no": "99999999"})
    assert res is None and "peak code 6" in text_of(content.blocks[-1])


def test_only_the_two_refusals_become_sentences(monkeypatch):
    def broken(self):
        raise RuntimeError("a fault")

    monkeypatch.setattr(Bulletin17B, "run_analysis", broken)
    with pytest.raises(RuntimeError, match="a fault"):
        build_section("b17b_ffa", dict(EX1))


# -- only when requested ------------------------------------------------------------------


def test_gaged_ffa_never_computes_b17b(monkeypatch):
    monkeypatch.setattr(
        "flowfreq.bulletin17b.Bulletin17B.run_analysis",
        lambda self: pytest.fail("B17B computed without being requested"),
    )
    gaged.clear_fit_cache()
    build_section(
        "gaged_ffa",
        {"peaks": peaks_spec("example1"), "regional_skew": 0.6, "regional_skew_se": 0.55},
    )


# -- the numbers --------------------------------------------------------------------------


def test_default_aeps_add_the_5_and_25_year_floods():
    assert DEFAULT_AEP == (0.99, 0.9, 0.5, 0.2, 0.1, 0.05, 0.04, 0.02, 0.01, 0.005, 0.002)
    assert 0.2 not in STANDARD_AEP and 0.04 not in STANDARD_AEP
    content, res = _b17b_ffa(dict(EX1))
    assert tuple(res.quantiles["aep"]) == DEFAULT_AEP
    assert [r[0] for r in tables(content)[B17B_CAPTION].rows][:4] == ["99", "90", "50", "20"]
    _, res = _b17b_ffa({**EX1, "aeps": [0.01, 0.1]})
    assert list(res.quantiles["aep"]) == [0.01, 0.1]


def test_quantile_table_carries_the_class_results():
    content, res = _b17b_ffa(dict(EX1))
    table = tables(content)[B17B_CAPTION]
    assert table.header == (
        "Annual exceedance probability (percent)",
        "Recurrence interval (years)",
        "Flow (cfs)",
        "Lower 90-percent limit (cfs)",
        "Upper 90-percent limit (cfs)",
        "Expected probability (percent)",
    )
    row = next(r for r in table.rows if r[0] == "1")
    q = res.quantiles.set_index("aep").loc[0.01]
    assert row == (
        "1",
        "100",
        cfs(q["flow_cfs"]),
        cfs(q["lower"]),
        cfs(q["upper"]),
        f"{100 * q['expected_probability']:.3g}",
    )
    rec = record("example1")
    alone = Bulletin17B(
        rec.peaks_cfs, rec.water_years, generalized_skew=0.6, aeps=[0.01]
    ).run_analysis()
    assert q["flow_cfs"] == pytest.approx(float(alone.quantiles["flow_cfs"][0]), rel=1e-12)
    assert "not applied to the flows" in all_text(content)


def test_hand_computation_rounds_the_skew():
    content, res = _b17b_ffa({**EX1, "hand_computation": True})
    assert res.skew_used == 0.7 and round(res.skew_weighted, 4) != 0.7
    assert "hand computation" in all_text(content)
    assert "hand computation" not in all_text(build_section("b17b_ffa", dict(EX1)))


def test_procedure_list_is_the_steps_and_parameters_follow_the_results():
    content, res = _b17b_ffa(dict(EX1))
    (steps,) = [b for b in content.blocks if isinstance(b, Bullets)]
    assert steps.items == tuple(res.steps)
    rows = parameters(content)
    assert rows["Station skew"] == f"{res.skew_station:.3f}"
    assert rows["Skew used"] == f"{res.skew_used:.3f}"
    assert rows["Mean of log10 flow"] == f"{res.mean_log:.4f}"
    assert rows["High outliers"] == "none" and rows["Low outliers"] == "none"


def test_big_sandy_from_codes_matches_appendix_6():
    content, res = _b17b_ffa({**SANDY, "perception_thresholds": SANDY_ROWS})
    assert res.historic.h == 77 and res.historic.z == 3
    flow = res.quantiles.set_index("aep").loc[0.01, "flow_cfs"]
    assert flow == pytest.approx(24391, rel=1e-3)
    rows = parameters(content)
    assert rows["Years in the historic period, H"] == "77"
    assert rows["Historic peaks, Z"] == "3"
    assert "Historically adjusted skew" in rows
    said = all_text(content)
    assert "3 historic (code 7)" in said and "1897 to 1929, 18,000 cfs" in said
    assert "historically weighted" in figure_of(content).caption


def test_a_conditional_adjustment_names_the_synthetic_statistics():
    content, res = _b17b_ffa({"peaks": peaks_spec("example3"), "generalized_skew": -0.1})
    assert res.conditional is not None
    rows = parameters(content)
    assert rows["Synthetic skew"] == f"{res.skew_station:.3f}"
    assert "Mean of log10 flow (synthetic)" in rows
    assert "P~" in "\n".join(rows)


# -- historic period ----------------------------------------------------------------------


def test_an_explicit_historic_period_wins_and_the_section_says_so():
    content, res = _b17b_ffa({**SANDY, "historic_period": [1890, 1973]})
    assert res.historic.h == 84
    said = all_text(content)
    assert "The historic period given, 1890 to 1973, was used" in said
    assert "1897 to 1973" in said
    assert content.options["historic_period"] == [1890, 1973]
    content, res = _b17b_ffa({**SANDY, "historic_period": [1897, 1973]})
    assert res.historic.h == 77 and "was used" not in all_text(content)
    content, res = _b17b_ffa({**SANDY, "historic_period": (1897, 1973)})
    assert content.options["historic_period"] == [1897, 1973]


def test_the_period_statement_survives_a_refusal():
    content, res = _b17b_ffa({**SANDY, "historic_period": [1900, 1973]})
    assert res is None and "outside the historic period" in text_of(content.blocks[-1])
    assert "The historic period given, 1900 to 1973, was used" in all_text(content)


def test_a_historic_period_that_misses_the_record_is_a_refusal():
    content, res = _b17b_ffa({**EX1, "historic_period": [1950, 1968]})
    assert res is None and "does not cover the systematic record" in all_text(content)


# -- codes --------------------------------------------------------------------------------


def test_removed_peaks_are_printed():
    content, res = _b17b_ffa({**EX1, "peaks": with_peaks([(1975, 90000.0, "3")])})
    assert res.n_years == 24
    assert "1975 (90,000 cfs, code 3)" in all_text(content)
    assert "emoved" not in all_text(build_section("b17b_ffa", dict(EX1)))


def test_removed_peaks_are_printed_when_17b_refuses_too():
    peaks = with_peaks([(1975, 90000.0, "3")])
    content, res = _b17b_ffa({**EX1, "peaks": peaks, "low_outlier_threshold": 2190.0})
    assert res is None and "1975 (90,000 cfs, code 3)" in all_text(content)


def test_codes_can_be_turned_off():
    given = {**EX1, "peaks": with_peaks([(1975, 2000.0, "3")]), "apply_peak_codes": False}
    content, res = _b17b_ffa(given)
    assert res.n_years == 25
    said = all_text(content)
    assert "codes were not applied" in said and "emoved" not in said


def test_a_threshold_over_a_historic_peak_refuses_without_codes():
    content, res = _b17b_ffa(
        {**SANDY, "apply_peak_codes": False, "perception_thresholds": SANDY_ROWS}
    )
    assert res is None
    said = text_of(content.blocks[-1])
    assert said.startswith(REFUSED) and "no measured peak" in said


def test_an_unobserved_threshold_period_is_named():
    content, res = _b17b_ffa(
        {
            **EX1,
            "perception_thresholds": [
                {"start_year": 1980, "end_year": 1985, "threshold_cfs": math.inf}
            ],
        }
    )
    assert res is not None and "1980 to 1985, not observed" in all_text(content)


# -- figure -------------------------------------------------------------------------------


def _scatters(res: Any) -> Dict[str, int]:
    """Points per legend entry of the figure, after checking all of them are visible."""
    fig = _frequency_figure(res, "somewhere")
    ax = fig.axes[0]
    labelled = [c for c in ax.collections if not c.get_label().startswith("_")]
    counts = {c.get_label(): len(c.get_offsets()) for c in labelled}
    points = np.vstack([np.asarray(c.get_offsets()) for c in labelled])
    left, right = ax.get_xlim()
    bottom, top = ax.get_ylim()
    title = ax.get_title()
    plt.close(fig)
    assert left < right and np.all((points[:, 0] >= left) & (points[:, 0] <= right))
    assert 0 < bottom and np.all((points[:, 1] >= bottom) & (points[:, 1] <= top))
    assert LABEL in title and "somewhere" in title
    return counts


def test_figure_drops_zero_and_below_base_years_and_counts_them():
    peaks = with_peaks([(1969, 0.0, None), (1970, 500.0, "4")])
    content, res = _b17b_ffa({**EX1, "peaks": peaks})
    assert res.n_zeros == 1 and res.n_below_base == 1
    caption = figure_of(content).caption
    assert "Not shown: 1 zero year and 1 year below the gage base." in caption
    assert "Not shown" not in figure_of(build_section("b17b_ffa", dict(EX1))).caption
    assert _scatters(res) == {"Annual peaks": 24}


def test_figure_marks_historic_peaks_and_low_outliers():
    _, sandy = _b17b_ffa(dict(SANDY))
    assert _scatters(sandy) == {"Annual peaks": 44, "Historic peaks": 3}
    _, low = _b17b_ffa({"peaks": peaks_spec("example3"), "generalized_skew": -0.1})
    assert low.low_outliers
    counts = _scatters(low)
    assert counts["Low outliers"] == len(low.low_outliers)
    assert sum(counts.values()) == low.n_years


def test_figure_shows_positions_beyond_the_usual_axis():
    """A long historic period puts the largest peak past the 0.2-percent flood."""
    _, res = _b17b_ffa({**EX1, "peaks": with_peaks([(1200, 50000.0, "7")])})
    assert res.historic.h == 769
    assert _scatters(res)["Historic peaks"] == 1


def test_figure_png_is_in_the_section():
    content = build_section("b17b_ffa", dict(EX1))
    assert figure_of(content).png.startswith(b"\x89PNG")
    assert "90-percent limits" in figure_of(content).caption


# -- options, offline rebuild -------------------------------------------------------------


def test_downloaded_peaks_are_written_into_the_options(monkeypatch):
    fetched = {
        "peaks": {**peaks_spec("example1"), "source": "test double", "retrieved": "2026-11-24"},
        "site_name": "Fishkill Creek",
        "drainage_area_sqmi": 190.0,
    }
    monkeypatch.setattr(gaged, "download_peaks", lambda site_no: copy.deepcopy(fetched))
    given = {"site_no": "99999999", "generalized_skew": 0.6}
    content = build_section("b17b_ffa", given)
    assert given == {"site_no": "99999999", "generalized_skew": 0.6}
    rows = dict(tables(content)["Station and peak-flow record."].rows)
    assert rows["Data source"] == "test double, retrieved 2026-11-24"
    assert rows["Station name"] == "Fishkill Creek"
    monkeypatch.setattr(gaged, "download_peaks", lambda site_no: pytest.fail("downloaded"))
    options = json.loads(json.dumps(content.options))
    again = build_section("b17b_ffa", options)
    assert again.options == content.options and again.blocks == content.blocks


def test_a_download_error_raises(monkeypatch):
    def down(site_no):
        raise ReportError("annual peaks for USGS 1 could not be downloaded: offline")

    monkeypatch.setattr(gaged, "download_peaks", down)
    with pytest.raises(ReportError, match="could not be downloaded"):
        build_section("b17b_ffa", {"site_no": "01373500"})


def test_given_options_are_not_changed():
    given = {
        "peaks": peaks_spec("appendix6_big_sandy"),
        "regional_skew_state": "TN",
        "historic_period": (1897, 1973),
        "perception_thresholds": copy.deepcopy(SANDY_ROWS),
    }
    before = copy.deepcopy(given)
    content = _build_b17b_ffa(given)
    assert given == before
    content.options["peaks"]["water_year"].append(1)
    assert given == before


# -- registry and session -----------------------------------------------------------------


def test_the_kind_is_registered_and_nowhere_else():
    kind = SECTION_KINDS["b17b_ffa"]
    assert kind.required == () and kind.build is _build_b17b_ffa
    assert list(SECTION_KINDS)[:3] == ["text", "gaged_ffa", "b17b_ffa"]
    assert kind.summary == (
        "Bulletin 17B (superseded) flood frequency at a gage, for comparison with Bulletin 17C. "
        "Computed only when this section is requested."
    )
    for estimates in ("site_summary", "comparison"):
        assert "b17b_ffa" not in SECTION_KINDS[estimates].optional


def test_site_summary_does_not_take_a_b17b_estimate():
    estimates = [{"label": "Old", "kind": "b17b_ffa", "options": dict(EX1)}]
    locations = [{"name": "Gage", "lat": 36.0, "lon": -88.2}]
    with pytest.raises(ReportError, match="b17b_ffa"):
        build_section("site_summary", {"estimates": estimates, "locations": locations})


def test_session_adds_the_section_and_returns_its_main_table(tmp_path):
    assert MAIN_TABLES["b17b_ffa"] == B17B_CAPTION
    session = ReportSession("Report", output_dir=tmp_path)
    added = session.add_section("b17b_ffa", dict(EX1))
    assert "Expected probability (percent)" in added["result"] and "cautions" not in added
    session.render("report.md")
    assert LABEL in (tmp_path / "report.md").read_text(encoding="utf-8")


def test_session_reply_carries_the_refusal(tmp_path):
    session = ReportSession("Report", output_dir=tmp_path)
    added = session.add_section("b17b_ffa", {**EX1, "low_outlier_threshold": 2190.0})
    assert added["result"] is None and LABEL in added["title"]
    (said,) = added["cautions"]
    assert said.startswith(REFUSED)
    session.render("report.md")
    assert REFUSED in (tmp_path / "report.md").read_text(encoding="utf-8")


def test_the_reference_gives_the_bulletin_and_its_address():
    content = build_section("b17b_ffa", dict(EX1))
    assert content.references == [BULLETIN_17B]
    assert "Interagency Advisory Committee on Water Data, 1982" in BULLETIN_17B
    assert "https://pubs.usgs.gov/unnumbered/70275162/report.pdf" in BULLETIN_17B
