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
        result = runner.invoke(cli, ["compare", "--peaks", str(bad_csv)])
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
        result = runner.invoke(cli, ["compare", "--peaks", str(good_csv)])
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
        result = runner.invoke(cli, ["compare", "--peaks", str(good_csv)])
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
        result = runner.invoke(cli, ["compare", "--peaks", str(good_csv)])
        assert result.exit_code != 0

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
            cli, ["compare", "--peaks", str(good_csv), "--output", str(out_path)]
        )
        assert result.exit_code == 0, result.output
        assert out_path.is_file()
        assert "# Engine comparison: fake" in out_path.read_text()
        assert "# Engine comparison" not in result.output
