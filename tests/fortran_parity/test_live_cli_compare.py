"""Live verification of the ``flowfreq compare`` CLI subcommand.

Requires the built f2py extension, same as the ``compare_engines`` function it
wraps. Uses Click's ``CliRunner`` rather than shelling out, matching how the
rest of this test suite avoids subprocess dependence on the installed console
script.
"""

from __future__ import annotations

import pandas as pd
import pytest

pytest.importorskip(
    "flowfreq.peakfqr",
    reason="Fortran extension not built; run python build_fortran/build.py "
    "(see docs/FORTRAN_UPLOAD.md)",
    # See tests/fortran_parity/test_fortran_oracles.py's exc_type comment: pytest
    # 9.1 narrowed importorskip's default exc_type to ModuleNotFoundError, but this
    # module deliberately re-raises as ImportError.
    exc_type=ImportError,
)

pytestmark = pytest.mark.requires_fortran


@pytest.fixture()
def powder_river_csv(tmp_path):
    from tests.fixtures.wymt_peaks import load_site

    site = load_site("06326500.00")
    years = sorted(site.peaks)
    path = tmp_path / "powder_river_peaks.csv"
    pd.DataFrame({"water_year": years, "peak_flow_cfs": [site.peaks[y] for y in years]}).to_csv(
        path, index=False
    )
    return path, site


def test_compare_prints_markdown_and_exits_zero_on_pass(powder_river_csv):
    from click.testing import CliRunner

    from flowfreq.cli import cli

    path, site = powder_river_csv
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "compare",
            "--peaks",
            str(path),
            "--site",
            "Powder River",
            "--regional-skew",
            str(site.regional_skew),
            "--regional-skew-se",
            str(site.regional_skew_mse**0.5),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "# Engine comparison: Powder River" in result.output
    assert "PASS" in result.output
    assert "## Quantiles (cfs)" in result.output


def test_compare_writes_output_file(powder_river_csv, tmp_path):
    from click.testing import CliRunner

    from flowfreq.cli import cli

    path, site = powder_river_csv
    out_path = tmp_path / "report.md"
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "compare",
            "--peaks",
            str(path),
            "--regional-skew",
            str(site.regional_skew),
            "--regional-skew-se",
            str(site.regional_skew_mse**0.5),
            "--output",
            str(out_path),
        ],
    )
    assert result.exit_code == 0, result.output
    assert out_path.is_file()
    assert "# Engine comparison" in out_path.read_text()


def test_compare_rejects_a_csv_missing_the_expected_columns(tmp_path):
    from click.testing import CliRunner

    from flowfreq.cli import cli

    bad_csv = tmp_path / "bad.csv"
    pd.DataFrame({"year": [2000, 2001, 2002], "flow": [100.0, 200.0, 150.0]}).to_csv(
        bad_csv, index=False
    )
    runner = CliRunner()
    result = runner.invoke(cli, ["compare", "--peaks", str(bad_csv)])
    assert result.exit_code != 0
    assert "water_year" in result.output


def test_compare_exits_nonzero_when_engines_disagree():
    """A comparison that fails its tolerance must exit nonzero.

    Cains Coulee used to fail on its own (a 0.058 skew_weighted gap, since
    fixed); the engines now agree to ~0.001%, so the failure is forced with a
    tolerance no real comparison can meet."""
    from click.testing import CliRunner

    from flowfreq.cli import cli
    from tests.fixtures.wymt_peaks import load_site

    site = load_site("06327450.00")
    years = sorted(site.peaks)

    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "cains_coulee_peaks.csv"
        pd.DataFrame({"water_year": years, "peak_flow_cfs": [site.peaks[y] for y in years]}).to_csv(
            path, index=False
        )

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "compare",
                "--peaks",
                str(path),
                "--regional-skew",
                str(site.regional_skew),
                "--regional-skew-se",
                str(site.regional_skew_mse**0.5),
                "--tolerance-pct",
                "1e-12",
            ],
        )
    assert result.exit_code != 0
    assert "FAIL" in result.output
