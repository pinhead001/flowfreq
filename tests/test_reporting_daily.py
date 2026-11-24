"""The `low_flow`, `flow_duration` and `future_flow` report sections.

Nothing here makes a request: the daily record is synthetic and carried inline, and the
download is a stub. The change factors are made up, and cited as made up.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np
import pandas as pd
import pytest

from flowfreq.future_flow import apply_change_factors, factor_set_from_dict
from flowfreq.lowflow import LowFlowFrequency
from flowfreq.regime import flow_duration_curve
from flowfreq.reporting import (
    SECTION_KINDS,
    Figure,
    Heading,
    Paragraph,
    ReportError,
    ReportSpec,
    SectionSpec,
    Table,
    build_report,
    build_section,
    daily,
)
from flowfreq.reporting.daily import (
    FLOW_DURATION_CAPTION,
    LOW_FLOW_CAPTION,
    RIGGS_1972,
    SEARCY_1959,
    download_daily,
)
from flowfreq.reporting.sections import FUTURE_FLOW_CAPTION, cfs
from flowfreq.reporting.session import MAIN_TABLES, ReportSession
from flowfreq.workflow import run_ffa
from tests.test_reporting_sites import FLOWS, GAGE, SKEW, TENNESSEE, YEARS

START = "1990-04-01"
DAYS = pd.date_range(START, "2020-03-31", freq="D")


def _flows() -> np.ndarray:
    rng = np.random.default_rng(1)
    doy = DAYS.dayofyear.to_numpy()
    wave = 200.0 + 150.0 * np.sin(2.0 * np.pi * (doy - 100) / 365.25)
    return np.round(wave + rng.gamma(2.0, 20.0, len(DAYS)), 1)


VALUES = _flows()
FRAME = pd.DataFrame({"flow_cfs": VALUES}, index=DAYS)
STATION: Dict[str, Any] = {
    "site_no": "01234567",
    "site_name": "Test Creek near Nowhere",
    "drainage_area_sqmi": 42.0,
    "daily": {"start": START, "flow_cfs": VALUES.tolist(), "source": "synthetic record"},
}
CITATION = "Example, A.B., 2024, A made-up climate guidance for the tests, table 2."
FACTORS: Dict[str, Any] = {
    "name": "made-up 2050",
    "aeps": [0.5, 0.1, 0.01],
    "factors": [1.05, 1.10, 1.20],
    "scenario": "SSP5-8.5",
    "horizon": "2050",
    "citation": CITATION,
    "level": "site",
}


def _table(blocks: List[Any], caption: str) -> Table:
    (table,) = [b for b in blocks if isinstance(b, Table) and b.caption == caption]
    return table


def _text(blocks: List[Any]) -> str:
    return "\n".join(b.text for b in blocks if isinstance(b, Paragraph))


class TestLowFlow:
    @pytest.fixture
    def built(self) -> Any:
        return build_section("low_flow", dict(STATION))

    def test_kind_is_registered(self):
        kind = SECTION_KINDS["low_flow"]
        assert kind.required == ()
        assert {"daily", "n_day", "non_exceedance", "year_type", "distribution"} <= set(
            kind.optional
        )
        assert MAIN_TABLES["low_flow"] == LOW_FLOW_CAPTION

    def test_quantiles_are_the_analysis_own(self, built):
        fit = LowFlowFrequency(FRAME)
        fit.run_analysis()
        expected = fit.compute_confidence_limits()
        table = _table(built.blocks, LOW_FLOW_CAPTION)
        assert [row[0] for row in table.rows] == ["50", "20", "10", "5", "2", "1"]
        assert [row[2] for row in table.rows] == ["7Q2", "7Q5", "7Q10", "7Q20", "7Q50", "7Q100"]
        for row, (_, want) in zip(table.rows, expected.iterrows()):
            assert row[3:] == (
                cfs(want["flow_cfs"]),
                cfs(want["lower_90pct"]),
                cfs(want["upper_90pct"]),
            )

    def test_title_method_reference_and_figure(self, built):
        assert built.title == "Low-flow frequency at USGS 01234567 Test Creek near Nowhere"
        text = _text(built.blocks)
        assert "annual minimum 7-day mean flow" in text
        assert "climatic year (April 1 to March 31)" in text and "Riggs, 1972" in text
        assert built.references == [RIGGS_1972]
        assert sum(isinstance(b, Figure) for b in built.blocks) == 1
        record = dict(_table(built.blocks, "Station and daily-flow record.").rows)
        assert record["Period used"] == "1990-04-01 to 2020-03-31"
        assert record["Days without one"] == "0"
        assert record["Data source"] == "synthetic record"

    def test_annual_minima_table_can_be_left_out(self, built):
        caption = "Annual minimum 7-day mean flows used in the analysis."
        assert len(_table(built.blocks, caption).rows) == 30
        brief = build_section("low_flow", {**STATION, "annual_table": False})
        assert not any(isinstance(b, Heading) for b in brief.blocks)

    def test_options_change_the_analysis(self):
        built = build_section(
            "low_flow",
            {
                **STATION,
                "n_day": 30,
                "non_exceedance": [0.1, 0.5, 0.3],
                "year_type": "water",
                "distribution": "lognormal",
                "annual_table": False,
            },
        )
        fit = LowFlowFrequency(FRAME, n_day=30, year_type="water", distribution="lognormal")
        fit.run_analysis()
        want = fit.compute_quantiles(np.array([0.5, 0.3, 0.1]))["flow_cfs"]
        table = _table(built.blocks, LOW_FLOW_CAPTION)
        assert [row[3] for row in table.rows] == [cfs(v) for v in want]
        # 1 / 0.3 is not a whole number of years, so the row has no nQy name.
        assert [row[2] for row in table.rows] == ["30Q2", "", "30Q10"]
        assert "lognormal distribution" in _text(built.blocks)
        assert "water year (October 1 to September 30)" in _text(built.blocks)

    def test_a_period_and_missing_days(self):
        values: List[Any] = VALUES.tolist()
        values[400:430] = [None] * 30
        built = build_section(
            "low_flow",
            {
                **STATION,
                "daily": {"start": START, "flow_cfs": values},
                "start_date": "1991-04-01",
                "annual_table": False,
            },
        )
        record = dict(_table(built.blocks, "Station and daily-flow record.").rows)
        assert record["Period used"] == "1991-04-01 to 2020-03-31"
        assert record["Days without one"] == "30"
        fitted = dict(
            _table(built.blocks, "Fitted distribution of the annual minimum 7-day mean flow.").rows
        )
        assert fitted["Years used"] == "28"
        assert "1 year(s) with too few days of record were left out" in _text(built.blocks)

    def test_zero_flow_years_are_said(self):
        values = VALUES.copy()
        for year in (1995, 2001, 2007):
            values[(DAYS.year == year) & (DAYS.month == 9)] = 0.0
        built = build_section(
            "low_flow",
            {**STATION, "daily": {"start": START, "flow_cfs": values.tolist()}},
        )
        text = _text(built.blocks)
        assert "3 of the 30 years had a minimum of zero" in text
        assert "conditional probability adjustment" in text

    @pytest.mark.parametrize(
        "change,said",
        [
            ({"n_day": 0}, "n_day must be a whole number of days"),
            ({"year_type": "fiscal"}, "year_type must be one of"),
            ({"distribution": "gev"}, "distribution must be lp3 or lognormal"),
            ({"non_exceedance": [0.1, 1.0]}, "strictly between 0 and 1"),
            ({"non_exceedance": [0.1, 0.1]}, "repeats a value"),
            ({"daily": {"start": START}}, "daily needs start"),
            ({"daily": {"start": "April 1", "flow_cfs": [1.0]}}, "must be a date as YYYY-MM-DD"),
            ({"daily": {"start": START, "flow_cfs": [1.0], "units": "cfs"}}, "unknown daily key"),
            ({"daily": {"start": START, "flow_cfs": ["wet"]}}, "is not numeric"),
            ({"end_date": "1980-01-01"}, "no value in the period asked for"),
            ({"end_date": "1995-03-31"}, "the low-flow analysis failed: Only 5 usable"),
        ],
    )
    def test_refusals(self, change, said):
        with pytest.raises(ReportError, match=said):
            build_section("low_flow", {**STATION, **change})

    def test_needs_a_station_or_a_record(self):
        with pytest.raises(ReportError, match="give site_no .* or daily"):
            build_section("low_flow", {"site_name": "Somewhere"})


class TestFlowDuration:
    @pytest.fixture
    def built(self) -> Any:
        return build_section("flow_duration", dict(STATION))

    def test_statistics_are_the_curve_own_and_lowercase(self, built):
        expected = flow_duration_curve(FRAME)
        table = _table(built.blocks, FLOW_DURATION_CAPTION)
        assert [row[1] for row in table.rows] == [
            "q1",
            "q5",
            "q10",
            "q20",
            "q50",
            "q80",
            "q90",
            "q95",
            "q99",
        ]
        assert [row[2] for row in table.rows] == [cfs(v) for v in expected["flow_cfs"]]
        assert MAIN_TABLES["flow_duration"] == FLOW_DURATION_CAPTION

    def test_title_method_reference_and_figure(self, built):
        assert built.title == "Flow duration at USGS 01234567 Test Creek near Nowhere"
        text = _text(built.blocks)
        assert "equaled or exceeded on N percent of those days" in text
        assert "not an annual probability" in text
        assert built.references == [SEARCY_1959]
        assert sum(isinstance(b, Figure) for b in built.blocks) == 1

    def test_percents_asked_for_and_missing_days(self):
        values: List[Any] = VALUES.tolist()
        values[100:150] = [None] * 50
        built = build_section(
            "flow_duration",
            {
                **STATION,
                "daily": {"start": START, "flow_cfs": values},
                "exceedance_pct": [90, 10, 50],
            },
        )
        table = _table(built.blocks, FLOW_DURATION_CAPTION)
        assert [row[0] for row in table.rows] == ["10", "50", "90"]
        assert "50 day(s) in the period have no value" in _text(built.blocks)

    @pytest.mark.parametrize(
        "change,said",
        [
            ({"exceedance_pct": [0, 50]}, "strictly between 0 and 100"),
            ({"exceedance_pct": []}, "non-empty list"),
            ({"n_day": 7}, r"unknown option\(s\) n_day"),
        ],
    )
    def test_refusals(self, change, said):
        with pytest.raises(ReportError, match=said):
            build_section("flow_duration", {**STATION, **change})


class TestDownload:
    @pytest.fixture
    def gage(self, monkeypatch) -> List[Any]:
        asked: List[Any] = []

        class FakeGage:
            def __init__(self, site_no: str) -> None:
                self.site_no = site_no
                self.site_name = "Test Creek near Nowhere"
                self.drainage_area = 42.0

            def download_daily_flow(self, start_date=None, end_date=None):
                asked.append((self.site_no, start_date, end_date))
                # Days with no published value are absent from the download.
                return FRAME.drop(FRAME.index[10:15])

        monkeypatch.setattr("flowfreq.usgs.USGSgage", FakeGage)
        return asked

    def test_the_record_is_stored_with_its_gaps(self, gage):
        fetched = download_daily("01234567", "1990-04-01", None)
        assert gage == [("01234567", "1990-04-01", None)]
        record = fetched["daily"]
        assert record["start"] == START and len(record["flow_cfs"]) == len(DAYS)
        assert record["flow_cfs"][10:15] == [None] * 5
        assert record["flow_cfs"][9] == VALUES[9]
        assert "site 01234567" in record["source"] and record["retrieved"]

    def test_a_section_downloads_once_and_rebuilds_offline(self, gage, monkeypatch):
        built = build_section("flow_duration", {"site_no": "01234567"})
        assert len(gage) == 1
        assert built.options["site_name"] == "Test Creek near Nowhere"
        assert built.options["drainage_area_sqmi"] == 42.0

        def no_network(*args: Any, **kwargs: Any) -> Any:
            raise AssertionError("the record is in the options")

        monkeypatch.setattr(daily, "download_daily", no_network)
        again = build_section("flow_duration", built.options)
        assert _table(again.blocks, FLOW_DURATION_CAPTION) == _table(
            built.blocks, FLOW_DURATION_CAPTION
        )

    def test_the_catalog_names_a_station_the_download_does_not(self, monkeypatch):
        class Unnamed:
            def __init__(self, site_no: str) -> None:
                self.site_no = site_no
                self.site_name = None
                self.drainage_area = None

            def download_daily_flow(self, **kwargs: Any) -> Any:
                return FRAME

        monkeypatch.setattr("flowfreq.usgs.USGSgage", Unnamed)
        known = download_daily("03606500")
        assert "BIG SANDY" in known["site_name"].upper()
        assert known["drainage_area_sqmi"] == pytest.approx(205.0, rel=0.02)
        unknown = download_daily("00000000")
        assert unknown["site_name"] is None and unknown["drainage_area_sqmi"] is None

    def test_a_failed_download_is_a_refusal(self, monkeypatch):
        class Broken:
            def __init__(self, site_no: str) -> None:
                self.site_no = site_no

            def download_daily_flow(self, **kwargs: Any) -> Any:
                raise ValueError("no daily data for the site")

        monkeypatch.setattr("flowfreq.usgs.USGSgage", Broken)
        with pytest.raises(ReportError, match="could not be downloaded: no daily data"):
            build_section("low_flow", {"site_no": "01234567"})


class TestFutureFlow:
    @pytest.fixture
    def built(self) -> Any:
        return build_section("future_flow", {"gage": dict(GAGE), "factor_set": dict(FACTORS)})

    def test_kind_is_registered(self):
        kind = SECTION_KINDS["future_flow"]
        assert kind.required == ("factor_set",)
        assert kind.optional == ("gage", "regression", "aeps", "detail")
        assert MAIN_TABLES["future_flow"] == FUTURE_FLOW_CAPTION

    def test_flows_are_the_fit_scaled(self, built):
        fit = run_ffa(np.array(FLOWS, dtype=float), np.array(YEARS), stationarity=False, **SKEW)
        present = fit["b17c"].results.quantiles
        present = present[(present["aep"] <= 0.5) & (present["aep"] >= 0.01)]
        expected = apply_change_factors(
            present.sort_values("aep", ascending=False), factor_set_from_dict(FACTORS)
        ).table
        table = _table(built.blocks, FUTURE_FLOW_CAPTION)
        assert len(table.rows) == len(expected) >= 5
        for row, (_, want) in zip(table.rows, expected.iterrows()):
            assert row[:4] == (
                f"{want['aep'] * 100:.3g}",
                cfs(want["regulatory_flow_cfs"]),
                f"{want['change_factor']:.3f}",
                cfs(want["future_flow_cfs"]),
            )
        by_aep = {row[0]: row for row in table.rows}
        assert by_aep["1"][2] == "1.200" and by_aep["1"][4] == "published"
        assert by_aep["2"][4] == "interpolated" and 1.10 < float(by_aep["2"][2]) < 1.20
        # Nothing outside the factors' range is scaled unless asked for.
        assert "0.2" not in by_aep and "99" not in by_aep

    def test_says_whose_factors_they_are(self, built):
        text = _text(built.blocks)
        assert "SSP5-8.5, 2050, from Example (2024)" in text
        assert "supplied by the analyst" in text and "ships no change factors" in text
        assert "limits do not carry over" in text
        assert built.references[-1] == CITATION
        assert built.title == (
            "Future-condition flood frequency at USGS 03606500 Big Sandy River at Bruceton, TN"
        )
        supplied = _table(built.blocks, "Change factors as supplied.")
        assert supplied.rows == (("50", "1.050"), ("10", "1.100"), ("1", "1.200"))

    def test_an_aep_beyond_the_factors_is_held_and_marked(self):
        built = build_section(
            "future_flow",
            {"gage": dict(GAGE), "factor_set": dict(FACTORS), "aeps": [0.01, 0.002]},
        )
        table = _table(built.blocks, FUTURE_FLOW_CAPTION)
        assert [(row[0], row[2], row[4]) for row in table.rows] == [
            ("1", "1.200", "published"),
            ("0.2", "1.200", "held beyond the published AEPs"),
        ]
        assert "nearest published factor" in _text(built.blocks)

    def test_detail_carries_the_present_analysis(self, built):
        headings = [b.text for b in built.blocks if isinstance(b, Heading)]
        assert headings[0] == "Present-condition analysis"
        brief = build_section(
            "future_flow", {"gage": dict(GAGE), "factor_set": dict(FACTORS), "detail": False}
        )
        assert not any(isinstance(b, Heading) for b in brief.blocks)

    def test_a_regression_estimate_is_scaled_too(self):
        built = build_section(
            "future_flow",
            {"regression": dict(TENNESSEE), "factor_set": dict(FACTORS), "detail": False},
        )
        assert "the regional regression estimate" in _text(built.blocks)
        table = _table(built.blocks, FUTURE_FLOW_CAPTION)
        by_aep = {row[0]: row for row in table.rows}
        present = float(by_aep["1"][1].replace(",", ""))
        future = float(by_aep["1"][3].replace(",", ""))
        assert future == pytest.approx(1.2 * present, rel=1e-3)
        assert built.title.startswith("Future-condition regression estimate")

    @pytest.mark.parametrize(
        "change,said",
        [
            ({"regression": dict(TENNESSEE)}, "one of the two"),
            ({"gage": None}, "one of the two"),
            ({"factor_set": {**FACTORS, "table": 2}}, "unknown factor_set key"),
            ({"factor_set": {k: v for k, v in FACTORS.items() if k != "citation"}}, "citation"),
            ({"factor_set": {**FACTORS, "factors": [1.0, -1.0, 1.0]}}, "factors must be positive"),
            ({"aeps": [0.03]}, "no quantile at AEP 0.03"),
            ({"factor_set": {**FACTORS, "aeps": [0.3, 0.03, 0.003]}}, None),
        ],
    )
    def test_refusals(self, change, said):
        given = {"gage": dict(GAGE), "factor_set": dict(FACTORS), "detail": False, **change}
        if said is None:
            # Factors at AEPs the fit does not tabulate still cover the AEPs between them.
            assert build_section("future_flow", given)
            return
        with pytest.raises(ReportError, match=said):
            build_section("future_flow", given)

    def test_no_aep_in_range_is_refused(self):
        narrow = {**FACTORS, "aeps": [0.03], "factors": [1.1]}
        with pytest.raises(ReportError, match="no AEP within the factor set's range"):
            build_section(
                "future_flow", {"gage": dict(GAGE), "factor_set": narrow, "detail": False}
            )


class TestInAReport:
    def test_all_three_build_save_and_rebuild(self, tmp_path):
        spec = ReportSpec(
            title="Test Creek",
            sections=[
                SectionSpec("low_flow", options={**STATION, "annual_table": False}),
                SectionSpec("flow_duration", options=dict(STATION)),
                SectionSpec(
                    "future_flow",
                    options={"gage": dict(GAGE), "factor_set": dict(FACTORS), "detail": False},
                ),
            ],
        )
        report = build_report(spec)
        text = report.to_markdown()
        assert LOW_FLOW_CAPTION in text and FLOW_DURATION_CAPTION in text
        assert FUTURE_FLOW_CAPTION in text
        assert RIGGS_1972 in text and SEARCY_1959 in text and CITATION in text
        saved = report.spec.save(tmp_path / "report.json")
        again = build_report(ReportSpec.load(saved))
        assert again.to_markdown() == text


class TestSession:
    @pytest.fixture
    def session(self, tmp_path) -> ReportSession:
        made = ReportSession("Big Sandy River", output_dir=tmp_path)
        made.add_section("gaged_ffa", GAGE)
        made.add_section("ungaged", TENNESSEE)
        made.add_text("Closing", "A paragraph.")
        return made

    def test_future_flow_goes_after_its_base_which_stays(self, session):
        added = session.future_flow("s1", FACTORS)
        assert added["section_id"] == "s4" and added["number"] == 2
        assert added["kind"] == "future_flow"
        assert added["result"].startswith(FUTURE_FLOW_CAPTION)
        assert "| published |" in added["result"]
        options = session.spec.sections[1].options
        assert options["gage"]["site_no"] == "03606500" and options["detail"] is False
        assert "regression" not in options

    def test_future_flow_of_a_regression_and_a_position(self, session):
        added = session.future_flow("s2", FACTORS, aeps=[0.01], position=4)
        assert added["number"] == 4
        assert "regression" in session.spec.sections[3].options
        assert added["result"].count("| published |") == 1

    def test_future_flow_refuses_another_kind(self, session):
        with pytest.raises(ReportError, match="not a gaged_ffa or an ungaged section"):
            session.future_flow("s3", FACTORS)
        assert len(session.outline()["sections"]) == 3

    def test_daily_sections_report_their_main_table(self, session):
        low = session.add_section("low_flow", {**STATION, "annual_table": False})
        assert low["result"].startswith(LOW_FLOW_CAPTION) and "| 7Q10 |" in low["result"]
        duration = session.add_section("flow_duration", dict(STATION))
        assert duration["result"].startswith(FLOW_DURATION_CAPTION)
        assert "| q50 |" in duration["result"]


@pytest.mark.parametrize("kind", ["low_flow", "flow_duration"])
def test_a_bad_drainage_area_is_refused_in_words(kind):
    with pytest.raises(ReportError, match="drainage_area_sqmi must be a positive number"):
        build_section(kind, {**STATION, "drainage_area_sqmi": "big"})


@pytest.mark.parametrize("kind", ["low_flow", "flow_duration"])
def test_a_record_with_no_usable_flow_is_refused_in_words(kind):
    bad = {**STATION, "daily": {"start": START, "flow_cfs": [-1.0] * 4000}}
    with pytest.raises(ReportError, match="flow"):
        build_section(kind, bad)
