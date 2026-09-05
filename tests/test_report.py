"""Tests for flowfreq.report.HydroReport.

The report is the artefact that leaves the library -- what gets attached to a
submittal and read by someone who did not run the analysis. So these tests are
mostly about whether the numbers on the page are the numbers that were fitted,
and whether sections that depend on optional data appear only when that data is
present.

No network: the gage is a stub carrying exactly the five attributes the report
reads (site_no, site_name, drainage_area, daily_data, peak_data). Building it
by hand rather than mocking USGSgage keeps the test honest about that surface --
if the report starts reading a sixth attribute, this fails rather than silently
picking up a Mock's auto-attribute.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
import pytest

from flowfreq.bulletin17c import Bulletin17C
from flowfreq.core import AnalysisMethod
from flowfreq.report import HydroReport
from tests.fixtures.big_sandy import REGIONAL_SKEW, REGIONAL_SKEW_SD, SYSTEMATIC_PEAKS


@dataclass
class StubGage:
    """Exactly the gage surface HydroReport touches. No more."""

    site_no: str = "03606500"
    site_name: Optional[str] = "Big Sandy River at Bruceton, TN"
    drainage_area: Optional[float] = 205.0
    daily_data: Optional[pd.DataFrame] = None
    peak_data: Optional[pd.DataFrame] = field(default=None)


def _peak_frame():
    years = sorted(SYSTEMATIC_PEAKS)
    return pd.DataFrame(
        {
            "water_year": years,
            "peak_date": pd.to_datetime([f"{y}-05-15" for y in years]),
            "peak_flow_cfs": [float(SYSTEMATIC_PEAKS[y]) for y in years],
        }
    )


def _daily_frame(days=400):
    """Minimal daily series: a DatetimeIndex and a flow_cfs column."""
    idx = pd.date_range("2000-10-01", periods=days, freq="D")
    flow = 500 + 400 * np.sin(np.arange(days) / 58.0) ** 2
    return pd.DataFrame({"flow_cfs": flow}, index=idx)


@pytest.fixture(scope="module")
def analysis_mom():
    years = np.array(sorted(SYSTEMATIC_PEAKS))
    flows = np.array([SYSTEMATIC_PEAKS[y] for y in years], dtype=float)
    b = Bulletin17C(peak_flows=flows, water_years=years)
    b.run_analysis(method="mom")
    return b


@pytest.fixture(scope="module")
def analysis_ema_weighted():
    years = np.array(sorted(SYSTEMATIC_PEAKS))
    flows = np.array([SYSTEMATIC_PEAKS[y] for y in years], dtype=float)
    b = Bulletin17C(
        peak_flows=flows,
        water_years=years,
        regional_skew=REGIONAL_SKEW,
        regional_skew_mse=REGIONAL_SKEW_SD**2,
    )
    b.run_analysis(method="ema")
    return b


@pytest.fixture
def report(analysis_mom):
    return HydroReport(StubGage(peak_data=_peak_frame()), analysis_mom)


class TestConstruction:
    def test_accepts_a_bulletin17c_facade(self, analysis_mom):
        """Bulletin17C is a facade; the report has to unwrap its analyzer."""
        rep = HydroReport(StubGage(), analysis_mom)
        assert rep._results is analysis_mom.results
        assert rep._analysis is analysis_mom._analyzer

    def test_accepts_the_underlying_analysis_directly(self, analysis_mom):
        inner = analysis_mom._analyzer
        rep = HydroReport(StubGage(), inner)
        assert rep._analysis is inner
        assert rep._results is inner.results

    def test_figures_starts_empty(self, report):
        assert report.figures == {}

    def test_figures_returns_a_copy(self, report):
        """A caller mutating the returned dict must not corrupt the report."""
        report.figures["injected"] = "/tmp/nope.png"
        assert report.figures == {}


class TestReportText:
    def test_names_the_site(self, report):
        text = report.generate_report_text()
        assert "03606500" in text
        assert "Big Sandy River at Bruceton, TN" in text

    def test_reports_the_fitted_moments_verbatim(self, report, analysis_mom):
        """The page must carry the numbers that were fitted, to the digit."""
        r = analysis_mom.results
        text = report.generate_report_text()
        assert f"{r.mean_log:.4f}" in text
        assert f"{r.std_log:.4f}" in text
        assert f"{r.skew_station:.4f}" in text

    def test_names_the_method(self, report):
        assert "Method of Moments" in report.generate_report_text()

    def test_ema_run_names_ema_and_reports_convergence(self, analysis_ema_weighted):
        text = HydroReport(StubGage(), analysis_ema_weighted).generate_report_text()
        assert "Expected Moments Algorithm" in text
        assert "EMA Converged" in text

    def test_mom_run_omits_the_ema_rows(self, report):
        """Reporting EMA iterations for a MOM fit would be a lie of omission."""
        text = report.generate_report_text()
        assert "EMA Iterations" not in text

    def test_regional_skew_rows_appear_only_when_supplied(self, report, analysis_ema_weighted):
        assert "Regional Skew Coefficient" not in report.generate_report_text()
        assert "Regional Skew Coefficient" in (
            HydroReport(StubGage(), analysis_ema_weighted).generate_report_text()
        )

    def test_drainage_area_omitted_when_unknown(self, analysis_mom):
        with_da = HydroReport(StubGage(), analysis_mom).generate_report_text()
        without = HydroReport(StubGage(drainage_area=None), analysis_mom).generate_report_text()
        assert "Drainage Area" in with_da
        assert "Drainage Area" not in without

    def test_unknown_site_name_falls_back_rather_than_printing_none(self, analysis_mom):
        text = HydroReport(StubGage(site_name=None), analysis_mom).generate_report_text()
        assert "Unknown" in text
        assert "None" not in text.splitlines()[2]

    def test_period_of_record_comes_from_the_peak_data(self, report):
        years = sorted(SYSTEMATIC_PEAKS)
        text = report.generate_report_text()
        assert f"WY {years[0]} - {years[-1]}" in text

    def test_period_of_record_omitted_without_peak_data(self, analysis_mom):
        text = HydroReport(StubGage(peak_data=None), analysis_mom).generate_report_text()
        assert "Period of Record" not in text

    def test_cites_bulletin_17c(self, report):
        text = report.generate_report_text()
        assert "Bulletin 17C" in text
        assert "10.3133/tm4B5" in text


class TestFrequencyTable:
    def test_one_row_per_confidence_limit(self, report, analysis_mom):
        text = report.generate_report_text()
        table = text.split("Table 1.")[1].split("## 5.")[0]
        rows = [ln for ln in table.splitlines() if ln.startswith("| ") and "---" not in ln]
        assert len(rows) == len(analysis_mom.results.confidence_limits) + 1  # + header

    def test_quantiles_match_the_analysis(self, report, analysis_mom):
        text = report.generate_report_text()
        for _, row in analysis_mom.results.confidence_limits.iterrows():
            assert f"{row['flow_cfs']:,.0f}" in text

    def test_bounds_are_ordered_on_every_row(self, report):
        text = report.generate_report_text()
        table = text.split("Table 1.")[1].split("## 5.")[0]
        for ln in table.splitlines():
            cells = [c.strip() for c in ln.split("|")[1:-1]]
            if len(cells) != 5 or not re.match(r"^[\d,]+$", cells[2]):
                continue
            flow, lo, hi = (float(cells[i].replace(",", "")) for i in (2, 3, 4))
            assert lo <= flow <= hi, ln


class TestAppendixPeakTable:
    def test_lists_every_peak(self, report):
        text = report.generate_report_text()
        appendix = text.split("Appendix B")[1]
        rows = [ln for ln in appendix.splitlines() if ln.startswith("| ") and "---" not in ln]
        assert len(rows) == len(SYSTEMATIC_PEAKS) + 1  # + header

    def test_handles_a_missing_peak_date(self, analysis_mom):
        """NWIS peak records can carry a null date; the appendix must not crash."""
        df = _peak_frame()
        df.loc[0, "peak_date"] = pd.NaT
        text = HydroReport(StubGage(peak_data=df), analysis_mom).generate_report_text()
        assert "N/A" in text.split("Appendix B")[1]

    def test_appendix_empty_without_peak_data(self, analysis_mom):
        text = HydroReport(StubGage(peak_data=None), analysis_mom).generate_report_text()
        appendix = text.split("Appendix B")[1]
        assert not [ln for ln in appendix.splitlines() if re.match(r"^\| \d{4} \|", ln)]


class TestSaveReport:
    def test_writes_what_generate_returns(self, report, tmp_path):
        """Identical but for the Report Date line, which is datetime.now()."""

        def without_date(text):
            return [ln for ln in text.splitlines() if not ln.startswith("**Report Date:**")]

        out = tmp_path / "report.md"
        report.save_report(str(out))
        assert without_date(out.read_text()) == without_date(report.generate_report_text())

    def test_creates_the_file(self, report, tmp_path):
        out = tmp_path / "nested" / "report.md"
        out.parent.mkdir()
        report.save_report(str(out))
        assert out.is_file() and out.stat().st_size > 0


class TestFigures:
    def test_skips_figures_whose_data_is_absent(self, report, tmp_path):
        """A stub with no daily data must not produce daily-data figures."""
        figs = report.generate_all_figures(str(tmp_path))
        assert "daily_timeseries" not in figs
        assert "summary_hydrograph" not in figs

    def test_creates_the_output_directory(self, report, tmp_path):
        out = tmp_path / "made" / "here"
        report.generate_all_figures(str(out))
        assert out.is_dir()

    def test_recorded_paths_exist_on_disk(self, report, tmp_path):
        for name, path in report.generate_all_figures(str(tmp_path)).items():
            assert (tmp_path / path.rsplit("/", 1)[-1]).is_file(), name

    def test_daily_data_produces_both_daily_figures(self, analysis_mom, tmp_path):
        """The two daily figures are separate blocks writing separate files."""
        gage = StubGage(peak_data=_peak_frame(), daily_data=_daily_frame())
        figs = HydroReport(gage, analysis_mom).generate_all_figures(str(tmp_path))
        assert "daily_timeseries" in figs
        assert "summary_hydrograph" in figs
        assert figs["daily_timeseries"] != figs["summary_hydrograph"]
        for path in figs.values():
            assert (tmp_path / path.rsplit("/", 1)[-1]).is_file()

    def test_frequency_curve_is_produced_regardless_of_gage_data(self, analysis_mom, tmp_path):
        """It comes from the analysis, not the gage, so it is never skipped."""
        figs = HydroReport(StubGage(), analysis_mom).generate_all_figures(str(tmp_path))
        assert "frequency_curve" in figs

    def test_generate_all_figures_returns_the_live_dict(self, report, tmp_path):
        """Documents an inconsistency rather than asserting it is right.

        The `figures` property returns a copy; this returns `self._figures`
        itself, so a caller mutating it corrupts the report's state. Harmless
        today because nothing downstream mutates it, but the two accessors
        disagree about ownership and only one of them is defensive.
        """
        returned = report.generate_all_figures(str(tmp_path))
        returned["injected"] = "/tmp/nope.png"
        assert "injected" in report.figures
