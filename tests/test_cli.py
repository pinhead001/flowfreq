"""Tests for flowfreq.cli.

The ``compare`` subcommand's happy path needs the built f2py extension and is
covered end-to-end in ``tests/fortran_parity/test_live_cli_compare.py``
instead. Here, ``compare`` is exercised only on the paths that don't need the
extension: CSV validation (before ``compare_engines`` is ever called) and the
``ImportError``-to-``ClickException`` mapping (via a monkeypatched
``compare_engines``, so this file runs the same with or without the extension
built).
"""

from __future__ import annotations

import pandas as pd
import pytest
from click.testing import CliRunner

from flowfreq.cli import cli
from flowfreq.validation.comparisons import ComparisonResult


class _FakeResult:
    def __init__(self, passed: bool) -> None:
        self.passed = passed


class TestValidate:
    def test_all_pass_exits_zero(self, monkeypatch):
        results = {
            "site_a": _FakeResult(passed=True),
            "site_b": _FakeResult(passed=True),
        }
        monkeypatch.setattr("flowfreq.validation.benchmarks.run_all_benchmarks", lambda: results)
        monkeypatch.setattr("flowfreq.validation.benchmarks.print_benchmark_report", lambda r: None)
        runner = CliRunner()
        result = runner.invoke(cli, ["validate"])
        assert result.exit_code == 0, result.output

    def test_any_failure_exits_nonzero(self, monkeypatch):
        results = {
            "site_a": _FakeResult(passed=True),
            "site_b": _FakeResult(passed=False),
        }
        monkeypatch.setattr("flowfreq.validation.benchmarks.run_all_benchmarks", lambda: results)
        monkeypatch.setattr("flowfreq.validation.benchmarks.print_benchmark_report", lambda r: None)
        runner = CliRunner()
        result = runner.invoke(cli, ["validate"])
        assert result.exit_code != 0


class TestBenchmark:
    def test_text_format_is_the_default(self, monkeypatch):
        monkeypatch.setattr("flowfreq.validation.benchmarks.run_all_benchmarks", lambda: {})
        monkeypatch.setattr(
            "flowfreq.validation.reports.generate_text_report", lambda r: "TEXT REPORT"
        )
        runner = CliRunner()
        result = runner.invoke(cli, ["benchmark"])
        assert result.exit_code == 0
        assert "TEXT REPORT" in result.output

    def test_json_format(self, monkeypatch):
        monkeypatch.setattr("flowfreq.validation.benchmarks.run_all_benchmarks", lambda: {})
        monkeypatch.setattr(
            "flowfreq.validation.reports.generate_json_report", lambda r: '{"ok": true}'
        )
        runner = CliRunner()
        result = runner.invoke(cli, ["benchmark", "--format", "json"])
        assert result.exit_code == 0
        assert '{"ok": true}' in result.output

    def test_invalid_format_rejected(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["benchmark", "--format", "xml"])
        assert result.exit_code != 0


class TestCompare:
    def test_rejects_a_csv_missing_the_expected_columns(self, tmp_path):
        bad_csv = tmp_path / "bad.csv"
        pd.DataFrame({"year": [2000, 2001, 2002], "flow": [100.0, 200.0, 150.0]}).to_csv(
            bad_csv, index=False
        )
        runner = CliRunner()
        result = runner.invoke(cli, ["compare", "--peaks", str(bad_csv), "--station-skew"])
        assert result.exit_code != 0
        assert "water_year" in result.output
        assert "peak_flow_cfs" in result.output

    def test_missing_extension_reports_as_a_click_exception(self, tmp_path, monkeypatch):
        """Without the built extension, compare_engines raises ImportError --
        the CLI must surface that as a clean CLI error, not a traceback."""
        good_csv = tmp_path / "peaks.csv"
        pd.DataFrame(
            {"water_year": [2000, 2001, 2002], "peak_flow_cfs": [100.0, 200.0, 150.0]}
        ).to_csv(good_csv, index=False)

        def _raise(*args, **kwargs):
            raise ImportError("flowfreq.peakfqr requires the f2py Fortran extension")

        monkeypatch.setattr("flowfreq.workflow.compare_engines", _raise)
        runner = CliRunner()
        result = runner.invoke(cli, ["compare", "--station-skew", "--peaks", str(good_csv)])
        assert result.exit_code != 0
        assert "f2py Fortran extension" in result.output
        assert "Traceback" not in result.output

    def test_writes_markdown_and_exits_zero_on_pass(self, tmp_path, monkeypatch):
        good_csv = tmp_path / "peaks.csv"
        pd.DataFrame(
            {"water_year": [2000, 2001, 2002], "peak_flow_cfs": [100.0, 200.0, 150.0]}
        ).to_csv(good_csv, index=False)

        class _FakeReport:
            comparison = ComparisonResult(passed=True)

            def to_markdown(self) -> str:
                return "# Engine comparison: fake\n\nPASS\n"

        monkeypatch.setattr("flowfreq.workflow.compare_engines", lambda **kw: _FakeReport())
        runner = CliRunner()
        result = runner.invoke(cli, ["compare", "--station-skew", "--peaks", str(good_csv)])
        assert result.exit_code == 0, result.output
        assert "# Engine comparison: fake" in result.output

    def test_exits_nonzero_when_comparison_fails(self, tmp_path, monkeypatch):
        good_csv = tmp_path / "peaks.csv"
        pd.DataFrame(
            {"water_year": [2000, 2001, 2002], "peak_flow_cfs": [100.0, 200.0, 150.0]}
        ).to_csv(good_csv, index=False)

        class _FakeReport:
            comparison = ComparisonResult(passed=False)

            def to_markdown(self) -> str:
                return "# Engine comparison: fake\n\nFAIL\n"

        monkeypatch.setattr("flowfreq.workflow.compare_engines", lambda **kw: _FakeReport())
        runner = CliRunner()
        result = runner.invoke(cli, ["compare", "--station-skew", "--peaks", str(good_csv)])
        assert result.exit_code != 0

    def test_no_skew_choice_is_a_usage_error(self, tmp_path, monkeypatch):
        good_csv = tmp_path / "peaks.csv"
        pd.DataFrame(
            {"water_year": [2000, 2001, 2002], "peak_flow_cfs": [100.0, 200.0, 150.0]}
        ).to_csv(good_csv, index=False)
        called = []
        monkeypatch.setattr("flowfreq.workflow.compare_engines", lambda **kw: called.append(kw))
        result = CliRunner().invoke(cli, ["compare", "--peaks", str(good_csv)])
        assert result.exit_code == 2, result.output
        assert "No regional skew chosen" in result.output
        assert "Traceback" not in result.output
        assert not called

    def test_skew_flags_reach_compare_engines(self, tmp_path, monkeypatch):
        good_csv = tmp_path / "peaks.csv"
        pd.DataFrame(
            {"water_year": [2000, 2001, 2002], "peak_flow_cfs": [100.0, 200.0, 150.0]}
        ).to_csv(good_csv, index=False)
        seen = {}

        class _FakeReport:
            comparison = ComparisonResult(passed=True)

            def to_markdown(self) -> str:
                return "ok"

        def _capture(**kw):
            seen.update(kw)
            return _FakeReport()

        monkeypatch.setattr("flowfreq.workflow.compare_engines", _capture)
        result = CliRunner().invoke(cli, ["compare", "--peaks", str(good_csv), "--default-skew"])
        assert result.exit_code == 0, result.output
        assert seen["use_default_skew"] is True
        assert seen["station_skew_only"] is False
        assert seen["regional_skew"] is None

    def test_historical_and_thresholds_reach_compare_engines(self, tmp_path, monkeypatch):
        peaks = tmp_path / "peaks.csv"
        pd.DataFrame({"water_year": [2000, 2001], "peak_flow_cfs": [100.0, 200.0]}).to_csv(
            peaks, index=False
        )
        hist = tmp_path / "hist.csv"
        pd.DataFrame({"water_year": [1897], "peak_flow_cfs": [25000.0]}).to_csv(hist, index=False)
        seen = {}

        class _FakeReport:
            comparison = ComparisonResult(passed=True)

            def to_markdown(self) -> str:
                return "ok"

        def _capture(**kw):
            seen.update(kw)
            return _FakeReport()

        monkeypatch.setattr("flowfreq.workflow.compare_engines", _capture)
        result = CliRunner().invoke(
            cli,
            [
                "compare",
                "--station-skew",
                "--peaks",
                str(peaks),
                "--historical",
                str(hist),
                "--threshold",
                "1890",
                "1929",
                "18000",
                "--threshold",
                "1974",
                "1987",
                "5000",
            ],
        )
        assert result.exit_code == 0, result.output
        assert seen["historical_peaks"] == [(1897, 25000.0)]
        assert seen["perception_thresholds"] == {(1890, 1929): 18000.0, (1974, 1987): 5000.0}

    def test_neither_flag_passes_none(self, tmp_path, monkeypatch):
        peaks = tmp_path / "peaks.csv"
        pd.DataFrame({"water_year": [2000], "peak_flow_cfs": [100.0]}).to_csv(peaks, index=False)
        seen = {}

        class _FakeReport:
            comparison = ComparisonResult(passed=True)

            def to_markdown(self) -> str:
                return "ok"

        monkeypatch.setattr(
            "flowfreq.workflow.compare_engines", lambda **kw: seen.update(kw) or _FakeReport()
        )
        result = CliRunner().invoke(cli, ["compare", "--station-skew", "--peaks", str(peaks)])
        assert result.exit_code == 0, result.output
        assert seen["historical_peaks"] is None
        assert seen["perception_thresholds"] is None

    @pytest.mark.parametrize(
        "threshold, message",
        [(["1930", "1920", "100"], "is after end"), (["1920", "1930", "0"], "must be positive")],
    )
    def test_bad_threshold_is_a_usage_error(self, tmp_path, threshold, message):
        peaks = tmp_path / "peaks.csv"
        pd.DataFrame({"water_year": [2000], "peak_flow_cfs": [100.0]}).to_csv(peaks, index=False)
        result = CliRunner().invoke(
            cli, ["compare", "--station-skew", "--peaks", str(peaks), "--threshold", *threshold]
        )
        assert result.exit_code == 2
        assert message in result.output

    def test_historical_csv_needs_the_columns(self, tmp_path):
        peaks = tmp_path / "peaks.csv"
        pd.DataFrame({"water_year": [2000], "peak_flow_cfs": [100.0]}).to_csv(peaks, index=False)
        hist = tmp_path / "hist.csv"
        pd.DataFrame({"year": [1897], "flow": [25000.0]}).to_csv(hist, index=False)
        result = CliRunner().invoke(
            cli, ["compare", "--station-skew", "--peaks", str(peaks), "--historical", str(hist)]
        )
        assert result.exit_code == 2
        assert "--historical CSV is missing" in result.output

    def test_output_option_writes_a_file_instead_of_stdout(self, tmp_path, monkeypatch):
        good_csv = tmp_path / "peaks.csv"
        pd.DataFrame(
            {"water_year": [2000, 2001, 2002], "peak_flow_cfs": [100.0, 200.0, 150.0]}
        ).to_csv(good_csv, index=False)
        out_path = tmp_path / "report.md"

        class _FakeReport:
            comparison = ComparisonResult(passed=True)

            def to_markdown(self) -> str:
                return "# Engine comparison: fake\n"

        monkeypatch.setattr("flowfreq.workflow.compare_engines", lambda **kw: _FakeReport())
        runner = CliRunner()
        result = runner.invoke(
            cli, ["compare", "--station-skew", "--peaks", str(good_csv), "--output", str(out_path)]
        )
        assert result.exit_code == 0, result.output
        assert out_path.is_file()
        assert "# Engine comparison: fake" in out_path.read_text()
        assert "# Engine comparison" not in result.output


class TestComparePeakCodes:
    """``--peaks`` honours a ``qualification_code`` column (#30)."""

    @staticmethod
    def _capture(monkeypatch):
        seen: dict = {}

        class _FakeReport:
            comparison = ComparisonResult(passed=True)

            def to_markdown(self) -> str:
                return "ok"

        monkeypatch.setattr(
            "flowfreq.workflow.compare_engines", lambda **kw: seen.update(kw) or _FakeReport()
        )
        return seen

    def test_codes_reach_compare_engines_as_strings(self, tmp_path, monkeypatch):
        peaks = tmp_path / "peaks.csv"
        pd.DataFrame(
            {
                "water_year": [2000, 2001, 2002],
                "peak_flow_cfs": [100.0, 200.0, 150.0],
                "qualification_code": ["7", "", "6"],
            }
        ).to_csv(peaks, index=False)
        seen = self._capture(monkeypatch)
        result = CliRunner().invoke(cli, ["compare", "--station-skew", "--peaks", str(peaks)])
        assert result.exit_code == 0, result.output
        # Read as strings: an all-numeric code column must not come back "7.0".
        assert seen["peak_codes"][0] == "7"
        assert seen["peak_codes"][2] == "6"
        assert seen["apply_peak_codes"] is True

    def test_ignore_peak_codes_turns_them_off(self, tmp_path, monkeypatch):
        peaks = tmp_path / "peaks.csv"
        pd.DataFrame(
            {"water_year": [2000], "peak_flow_cfs": [100.0], "qualification_code": ["6"]}
        ).to_csv(peaks, index=False)
        seen = self._capture(monkeypatch)
        result = CliRunner().invoke(
            cli, ["compare", "--station-skew", "--ignore-peak-codes", "--peaks", str(peaks)]
        )
        assert result.exit_code == 0, result.output
        assert seen["apply_peak_codes"] is False

    def test_no_code_column_passes_none(self, tmp_path, monkeypatch):
        peaks = tmp_path / "peaks.csv"
        pd.DataFrame({"water_year": [2000], "peak_flow_cfs": [100.0]}).to_csv(peaks, index=False)
        seen = self._capture(monkeypatch)
        result = CliRunner().invoke(cli, ["compare", "--station-skew", "--peaks", str(peaks)])
        assert result.exit_code == 0, result.output
        assert seen["peak_codes"] is None

    def test_an_uncodable_record_is_a_clean_error_naming_the_year(self, tmp_path):
        """Real compare_engines: the code step fails before the extension is needed.

        A code 4 peak alone is an interval peak now; a code 4 peak that is also
        historic (code 7) is still one neither engine can express."""
        peaks = tmp_path / "peaks.csv"
        years = list(range(2001, 2016))
        codes = [""] * len(years)
        codes[2] = "4,7"
        pd.DataFrame(
            {
                "water_year": years,
                "peak_flow_cfs": [100.0 + 10 * i for i in range(len(years))],
                "qualification_code": codes,
            }
        ).to_csv(peaks, index=False)
        result = CliRunner().invoke(cli, ["compare", "--station-skew", "--peaks", str(peaks)])
        assert result.exit_code != 0
        assert "2003" in result.output and "historic" in result.output
        assert "--ignore-peak-codes" in result.output
        assert "Traceback" not in result.output

    def test_codes_with_historical_csv_is_a_clean_error(self, tmp_path):
        peaks = tmp_path / "peaks.csv"
        pd.DataFrame(
            {
                "water_year": [2000, 2001, 2002],
                "peak_flow_cfs": [100.0, 200.0, 150.0],
                "qualification_code": ["7", "", ""],
            }
        ).to_csv(peaks, index=False)
        hist = tmp_path / "hist.csv"
        pd.DataFrame({"water_year": [1897], "peak_flow_cfs": [25000.0]}).to_csv(hist, index=False)
        result = CliRunner().invoke(
            cli,
            ["compare", "--station-skew", "--peaks", str(peaks), "--historical", str(hist)],
        )
        assert result.exit_code != 0
        assert "historical_peaks" in result.output
        assert "Traceback" not in result.output
