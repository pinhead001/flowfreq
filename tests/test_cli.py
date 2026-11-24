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

    def test_a_historic_interval_record_compares(self, tmp_path):
        """Real compare_engines: a code 4 peak that is also historic (code 7) is fitted.

        It used to be refused as a record neither engine could express; both
        now build siteQT's dtype = 1 interval row for it."""
        pytest.importorskip("flowfreq.peakfqr", exc_type=ImportError)
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
        assert result.exit_code == 0, result.output
        assert "**PASS**" in result.output
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


class _FakeEstimate:
    def to_markdown(self) -> str:
        return "# Ungaged estimate\nSYNTHETIC-MARKDOWN"


class TestUngaged:
    def test_prints_markdown_and_passes_options(self, monkeypatch):
        calls = {}

        def fake(lat, lon, **kw):
            calls.update(lat=lat, lon=lon, **kw)
            return _FakeEstimate()

        monkeypatch.setattr("flowfreq.ungaged.estimate_at_site", fake)
        result = CliRunner().invoke(
            cli,
            [
                "ungaged", "--lat", "46.5", "--lon", "-110.2", "--state", "MT",
                "--aep", "0.01", "--aep", "0.1",
                "--allow-nss-snapshot", "--allow-extrapolation",
            ],
        )  # fmt: skip
        assert result.exit_code == 0, result.output
        assert "SYNTHETIC-MARKDOWN" in result.output
        assert calls["lat"] == 46.5 and calls["lon"] == -110.2 and calls["state"] == "MT"
        assert calls["aeps"] == [0.01, 0.1]
        assert calls["allow_nss_snapshot"] is True
        assert calls["allow_extrapolation"] is True

    def test_defaults_leave_aeps_unset(self, monkeypatch):
        calls = {}

        def fake(lat, lon, **kw):
            calls.update(kw)
            return _FakeEstimate()

        monkeypatch.setattr("flowfreq.ungaged.estimate_at_site", fake)
        result = CliRunner().invoke(
            cli, ["ungaged", "--lat", "46.5", "--lon", "-110.2", "--state", "MT"]
        )
        assert result.exit_code == 0, result.output
        assert calls["aeps"] is None
        assert calls["allow_nss_snapshot"] is False
        assert calls["allow_extrapolation"] is False
        assert calls["region"] is None

    @pytest.mark.parametrize(
        "exc_name", ["EquationsUnavailable", "UnsnappablePointError", "UnsupportedRegionError"]
    )
    def test_known_errors_exit_one_with_message(self, monkeypatch, exc_name):
        from flowfreq.regression.library import EquationsUnavailable
        from flowfreq.streamstats import UnsnappablePointError, UnsupportedRegionError

        exc = {
            "EquationsUnavailable": EquationsUnavailable,
            "UnsnappablePointError": UnsnappablePointError,
            "UnsupportedRegionError": UnsupportedRegionError,
        }[exc_name]

        def fake(lat, lon, **kw):
            raise exc("no equations for the test")

        monkeypatch.setattr("flowfreq.ungaged.estimate_at_site", fake)
        result = CliRunner().invoke(cli, ["ungaged", "--lat", "1", "--lon", "2", "--state", "ZZ"])
        assert result.exit_code == 1
        assert "no equations for the test" in result.output
        assert "Traceback" not in result.output

    def test_cli_region_option_and_overlap_message(self, monkeypatch):
        from flowfreq.ungaged import OverlappingRegionsError

        calls = {}

        def fake(lat, lon, **kw):
            calls.update(kw)
            if kw["region"] is None:
                raise OverlappingRegionsError(
                    "TN",
                    [
                        ("GC344", "MultiVariable_Area_2", 2003),
                        ("GC2016", "TN_Peakflow_Region_2_SIR_2024-5130", 2025),
                    ],
                )
            return _FakeEstimate()

        monkeypatch.setattr("flowfreq.ungaged.estimate_at_site", fake)
        base = ["ungaged", "--lat", "36.3", "--lon", "-87.0", "--state", "TN"]
        result = CliRunner().invoke(cli, base)
        assert result.exit_code == 1
        assert "Traceback" not in result.output
        lines = result.output.splitlines()
        assert any("GC344" in ln and "2003" in ln and "GC2016" not in ln for ln in lines)
        assert any("GC2016" in ln and "2025" in ln and "GC344" not in ln for ln in lines)
        assert "--region" in result.output

        result = CliRunner().invoke(cli, base + ["--region", "GC2016"])
        assert result.exit_code == 0, result.output
        assert calls["region"] == "GC2016"
        assert "SYNTHETIC-MARKDOWN" in result.output

    def test_requires_state(self):
        result = CliRunner().invoke(cli, ["ungaged", "--lat", "1", "--lon", "2"])
        assert result.exit_code != 0

    def test_out_of_range_hints_at_extrapolation_flag(self, monkeypatch):
        from flowfreq.regression.equations import OutOfRangeError

        def fake(lat, lon, **kw):
            raise OutOfRangeError("DRNAREA 5 below calibrated minimum 10")

        monkeypatch.setattr("flowfreq.ungaged.estimate_at_site", fake)
        result = CliRunner().invoke(cli, ["ungaged", "--lat", "1", "--lon", "2", "--state", "MT"])
        assert result.exit_code == 1
        assert "Error: DRNAREA 5 below calibrated minimum 10" in result.output
        assert "--allow-extrapolation" in result.output
        assert "Traceback" not in result.output

    def test_transport_error_is_one_line_exit_one(self, monkeypatch):
        import requests

        def fake(lat, lon, **kw):
            raise requests.ConnectionError("StreamStats unreachable")

        monkeypatch.setattr("flowfreq.ungaged.estimate_at_site", fake)
        result = CliRunner().invoke(cli, ["ungaged", "--lat", "1", "--lon", "2", "--state", "MT"])
        assert result.exit_code == 1
        assert "Error: StreamStats unreachable" in result.output
        assert "Traceback" not in result.output

    def test_unknown_state_is_clean_error(self):
        # No monkeypatch: the unknown code is rejected before any network call.
        result = CliRunner().invoke(cli, ["ungaged", "--lat", "1", "--lon", "2", "--state", "ZZ"])
        assert result.exit_code == 1
        assert "Unknown jurisdiction code 'ZZ'" in result.output
        assert "Traceback" not in result.output

    def test_aep_out_of_range_rejected(self):
        result = CliRunner().invoke(
            cli, ["ungaged", "--lat", "1", "--lon", "2", "--state", "MT", "--aep", "1.5"]
        )
        assert result.exit_code == 2


class TestReport:
    @staticmethod
    def _spec(tmp_path, **options):
        import json

        years = list(range(1990, 2020))
        flows = [900.0 + 37.0 * ((7 * i) % 23) for i in range(len(years))]
        spec = {
            "title": "Example Creek",
            "date": "2026-11-14",
            "sections": [
                {"kind": "text", "title": "Purpose", "narrative": "Why."},
                {
                    "kind": "gaged_ffa",
                    "options": {
                        "peaks": {"water_year": years, "peak_flow_cfs": flows},
                        "regional_skew_state": "IA",
                        "stationarity": False,
                        **options,
                    },
                },
            ],
        }
        path = tmp_path / "spec.json"
        path.write_text(json.dumps(spec), encoding="utf-8")
        return path

    def test_writes_markdown_and_the_spec_as_built(self, tmp_path):
        from flowfreq.reporting import ReportSpec

        out, frozen = tmp_path / "report.md", tmp_path / "frozen.json"
        result = CliRunner().invoke(
            cli, ["report", str(self._spec(tmp_path)), "-o", str(out), "--save-spec", str(frozen)]
        )
        assert result.exit_code == 0, result.output
        assert f"Wrote {out} and {frozen}" in result.output
        text = out.read_text(encoding="utf-8")
        assert text.splitlines()[0] == "# Example Creek"
        assert "Flood quantiles with 90-percent confidence limits." in text
        assert (tmp_path / "report_files" / "figure_1.png").exists()
        # The looked-up skew is recorded as a value, so the frozen spec needs no lookup.
        options = ReportSpec.load(frozen).sections[1].options
        assert options["regional_skew"] == pytest.approx(-0.4)
        assert "regional_skew_state" not in options

    def test_writes_docx(self, tmp_path):
        pytest.importorskip("docx")
        out = tmp_path / "report.docx"
        result = CliRunner().invoke(cli, ["report", str(self._spec(tmp_path)), "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert result.output.strip() == f"Wrote {out}"
        assert out.read_bytes()[:2] == b"PK"

    def test_a_spec_error_is_one_line_exit_one(self, tmp_path):
        spec = self._spec(tmp_path, regional_skew=0.1)
        result = CliRunner().invoke(cli, ["report", str(spec), "-o", str(tmp_path / "r.md")])
        assert result.exit_code == 1
        assert "Error: section 2 (gaged_ffa)" in result.output
        assert "Traceback" not in result.output
        assert not (tmp_path / "r.md").exists()

    def test_unsupported_output_suffix(self, tmp_path):
        result = CliRunner().invoke(
            cli, ["report", str(self._spec(tmp_path)), "-o", str(tmp_path / "r.pdf")]
        )
        assert result.exit_code == 1
        assert "use .docx or .md" in result.output

    def test_missing_spec_file_is_a_usage_error(self, tmp_path):
        result = CliRunner().invoke(
            cli, ["report", str(tmp_path / "absent.json"), "-o", str(tmp_path / "r.md")]
        )
        assert result.exit_code == 2


class TestUngagedSupplied:
    BASE = ["ungaged", "--lat", "29.43", "--lon", "-100", "--state", "TX"]

    def test_characteristics_reach_the_estimate(self, monkeypatch):
        seen = {}

        class Done:
            def to_markdown(self):
                return "done"

        def fake(lat, lon, **kwargs):
            seen.update(kwargs)
            return Done()

        monkeypatch.setattr("flowfreq.ungaged.estimate_at_site", fake)
        args = ["--char", "MCSLOPE=0.00326", "--char", "CONTDA=737"]
        result = CliRunner().invoke(cli, self.BASE + args)
        assert result.exit_code == 0, result.output
        assert seen["supplied_characteristics"] == {"MCSLOPE": 0.00326, "CONTDA": 737.0}
        result = CliRunner().invoke(cli, self.BASE)
        assert result.exit_code == 0 and seen["supplied_characteristics"] is None

    @pytest.mark.parametrize("bad", ["MCSLOPE", "MCSLOPE=steep", "=3"])
    def test_a_bad_one_is_refused(self, bad, monkeypatch):
        monkeypatch.setattr("flowfreq.ungaged.estimate_at_site", lambda *a, **k: 1 / 0)
        result = CliRunner().invoke(cli, self.BASE + ["--char", bad])
        assert result.exit_code == 2 and "CODE=VALUE" in result.output
