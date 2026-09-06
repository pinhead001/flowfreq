"""
FlowFreq command-line interface.

Provides CLI commands for validation and benchmarking of the
Bulletin 17C implementation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import click


@click.group()
def cli() -> None:
    """FlowFreq - Hydrologic frequency analysis tools."""
    pass


@cli.command()
def validate() -> None:
    """Run validation benchmarks against the reference expectations."""
    from flowfreq.validation.benchmarks import print_benchmark_report, run_all_benchmarks

    click.echo("Running validation benchmarks...")
    results = run_all_benchmarks()
    print_benchmark_report(results)

    n_pass = sum(1 for r in results.values() if r.passed)
    n_total = len(results)
    if n_pass < n_total:
        raise SystemExit(1)


@cli.command()
@click.option("--format", "fmt", type=click.Choice(["text", "json"]), default="text")
def benchmark(fmt: str) -> None:
    """Run benchmarks and generate a report.

    Parameters
    ----------
    fmt : str
        Output format: 'text' or 'json'.
    """
    from flowfreq.validation.benchmarks import run_all_benchmarks
    from flowfreq.validation.reports import generate_json_report, generate_text_report

    click.echo("Running benchmarks...")
    results = run_all_benchmarks()

    if fmt == "json":
        click.echo(generate_json_report(results))
    else:
        click.echo(generate_text_report(results))


@cli.command()
@click.option(
    "--peaks",
    "peaks_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help=(
        "CSV with 'water_year' and 'peak_flow_cfs' columns -- the same shape "
        "flowfreq.usgs.USGSgage.download_peak_flow produces, so its output can be "
        "saved straight to CSV and used here."
    ),
)
@click.option("--site", "site_name", default="", help="Site name/number for the report title.")
@click.option(
    "--regional-skew",
    type=float,
    default=None,
    help="Regional skew coefficient. Defaults to the nationwide B17C generalized skew (-0.302).",
)
@click.option("--regional-skew-se", type=float, default=0.55, help="Regional skew standard error.")
@click.option(
    "--low-outlier-threshold",
    type=float,
    default=None,
    help="User-supplied PILF threshold in cfs. Omit to let MGBT decide, as both engines do by default.",
)
@click.option(
    "--tolerance-pct",
    type=float,
    default=1.0,
    help="Quantile agreement tolerance, percent, for the pass/fail verdict.",
)
@click.option(
    "--output",
    "output_path",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Write the markdown report here instead of stdout.",
)
def compare(
    peaks_path: Path,
    site_name: str,
    regional_skew: Optional[float],
    regional_skew_se: float,
    low_outlier_threshold: Optional[float],
    tolerance_pct: float,
    output_path: Optional[Path],
) -> None:
    """Compare the native EMA against the vendored USGS Fortran (peakfq 8.1.0).

    Requires the f2py extension to be built (``python build_fortran/build.py``,
    needs gfortran and meson) -- this is the one thing this command cannot do
    without, and it says so rather than silently falling back to anything, per
    ``docs/FORTRAN_ENGINE_DESIGN.md`` section 9.

    Historical peaks and perception thresholds are not yet exposed here --
    a record that needs them should go through
    :func:`flowfreq.workflow.compare_engines` directly.
    """
    import pandas as pd

    from flowfreq.workflow import B17C_DEFAULT_SKEW, compare_engines

    peaks_df = pd.read_csv(peaks_path)
    missing = {"water_year", "peak_flow_cfs"} - set(peaks_df.columns)
    if missing:
        raise click.UsageError(
            f"--peaks CSV is missing column(s) {sorted(missing)}; expected 'water_year' and "
            "'peak_flow_cfs' (the shape USGSgage.download_peak_flow produces)."
        )

    skew = regional_skew if regional_skew is not None else B17C_DEFAULT_SKEW

    try:
        report = compare_engines(
            peak_flows=peaks_df["peak_flow_cfs"].to_numpy(dtype=float),
            water_years=peaks_df["water_year"].to_numpy(dtype=int),
            regional_skew=skew,
            regional_skew_se=regional_skew_se,
            user_low_outlier_threshold=low_outlier_threshold,
            site_name=site_name,
            tolerance_pct=tolerance_pct,
        )
    except ImportError as exc:
        raise click.ClickException(str(exc)) from exc

    markdown = report.to_markdown()
    if output_path is not None:
        output_path.write_text(markdown)
        click.echo(f"Wrote {output_path}")
    else:
        click.echo(markdown)

    if not report.comparison.passed:
        raise SystemExit(1)


if __name__ == "__main__":
    cli()
