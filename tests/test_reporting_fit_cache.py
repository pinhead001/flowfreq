"""One Bulletin 17C fit per gage per report, however many sections show it.

From the review of 2026-11-19 (TODO.md, item 4): a report with a gage section, a summary and
a comparison fitted the gage three times (and again on every render), about a second each.
``flowfreq.reporting.gaged`` now memoises the fit on the gage's options. Nothing here makes a
request: the gage carries its peaks and the regression its characteristics.
"""

from __future__ import annotations

import copy
from typing import Any, Dict, List

import pytest

import flowfreq.workflow as workflow
from flowfreq.reporting import build_section, gaged
from flowfreq.reporting.session import ReportSession
from tests.test_reporting_sites import GAGE, TENNESSEE

LOCATIONS: List[Dict[str, Any]] = [
    {"name": "Gage", "lat": 36.04, "lon": -88.23, "characteristics": {"DRNAREA": 205.0}},
]


@pytest.fixture(autouse=True)
def cold_cache():
    gaged.clear_fit_cache()
    yield
    gaged.clear_fit_cache()


@pytest.fixture
def fits(monkeypatch) -> List[Dict[str, Any]]:
    """Count calls of the fit behind every gaged section."""
    calls: List[Dict[str, Any]] = []
    real = workflow.run_ffa

    def counted(*args: Any, **kwargs: Any) -> Any:
        calls.append(kwargs)
        return real(*args, **kwargs)

    monkeypatch.setattr(workflow, "run_ffa", counted)
    return calls


def test_a_summary_and_a_comparison_reuse_the_gage_fit(tmp_path, fits):
    session = ReportSession("Big Sandy River", output_dir=tmp_path)
    session.add_section("gaged_ffa", dict(GAGE))  # s1
    session.add_section("ungaged", dict(TENNESSEE))  # s2
    session.summarize(["s1", "s2"], LOCATIONS, labels=["Gage", "Regression"], site_name="x")
    session.compare(["s1", "s2"], labels=["Gage", "Regression"], site_name="x")
    session.weight("s1", "s2", keep_parts=True)
    session.render("report.md")
    assert len(fits) == 1


def test_different_options_fit_again(fits):
    build_section("gaged_ffa", dict(GAGE))
    build_section("gaged_ffa", {**GAGE, "regional_skew": -0.3})
    build_section("gaged_ffa", {**GAGE, "low_outlier_threshold": 1500.0})
    assert len(fits) == 3
    build_section("gaged_ffa", {**GAGE, "regional_skew": -0.3})
    assert len(fits) == 3


def test_different_peaks_fit_again(fits):
    build_section("gaged_ffa", dict(GAGE))
    changed = copy.deepcopy(GAGE)
    changed["peaks"]["peak_flow_cfs"][0] *= 2.0
    build_section("gaged_ffa", changed)
    assert len(fits) == 2


def test_clearing_the_cache_fits_again(fits):
    build_section("gaged_ffa", dict(GAGE))
    gaged.clear_fit_cache()
    build_section("gaged_ffa", dict(GAGE))
    assert len(fits) == 2


def test_a_reused_fit_builds_the_same_section(fits):
    first = build_section("gaged_ffa", dict(GAGE))
    second = build_section("gaged_ffa", dict(GAGE))
    assert len(fits) == 1
    assert second.title == first.title
    assert second.options == first.options
    assert len(second.blocks) == len(first.blocks)
    assert [type(b) for b in second.blocks] == [type(b) for b in first.blocks]


def test_a_caller_changing_a_result_does_not_change_the_next(fits):
    first = build_section("gaged_ffa", copy.deepcopy(GAGE))
    first.blocks.clear()
    first.options["site_name"] = "changed by a caller"
    first.options["peaks"]["water_year"].append(9999)
    second = build_section("gaged_ffa", copy.deepcopy(GAGE))
    assert len(fits) == 1
    assert second.blocks
    assert second.options["site_name"] == GAGE["site_name"]
    assert 9999 not in second.options["peaks"]["water_year"]


def test_the_given_options_are_not_changed(fits):
    given = copy.deepcopy(GAGE)
    build_section("gaged_ffa", given)
    build_section("gaged_ffa", given)
    assert given == GAGE


def test_the_cache_is_bounded(fits):
    for k in range(gaged.FIT_CACHE_SIZE + 1):
        build_section("gaged_ffa", {**GAGE, "low_outlier_threshold": 1000.0 + k})
    assert len(fits) == gaged.FIT_CACHE_SIZE + 1
    build_section("gaged_ffa", {**GAGE, "low_outlier_threshold": 1000.0})
    assert len(fits) == gaged.FIT_CACHE_SIZE + 2


def test_a_failed_fit_is_not_cached(monkeypatch):
    calls: List[int] = []

    def failing(*args: Any, **kwargs: Any) -> Dict[str, Any]:
        calls.append(1)
        return {"error": "made to fail", "parameters": None, "b17c": None}

    monkeypatch.setattr(workflow, "run_ffa", failing)
    for _ in range(2):
        with pytest.raises(Exception, match="made to fail"):
            build_section("gaged_ffa", dict(GAGE))
    assert len(calls) == 2


def test_a_downloaded_gage_is_fitted_once_too(tmp_path, fits, monkeypatch):
    """Found by a live chat run, 2026-11-24: the download's own fit was not kept."""
    downloads: List[str] = []

    def downloaded(site_no: str) -> Dict[str, Any]:
        downloads.append(site_no)
        return {"peaks": copy.deepcopy(GAGE["peaks"]), "site_name": "Big Sandy River"}

    monkeypatch.setattr(gaged, "download_peaks", downloaded)
    given = {k: v for k, v in GAGE.items() if k not in ("peaks", "site_name")}
    session = ReportSession("Big Sandy River", output_dir=tmp_path)
    session.add_section("gaged_ffa", dict(given))
    session.render("report.md")
    assert len(fits) == 1 and downloads == ["03606500"]
    # Asked for again without peaks, it is downloaded and fitted again: the record may have grown.
    build_section("gaged_ffa", dict(given))
    assert len(fits) == 2 and len(downloads) == 2
