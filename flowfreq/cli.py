"""
FlowFreq command-line interface.

Provides CLI commands for validation and benchmarking of the
Bulletin 17C implementation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

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
        "saved straight to CSV and used here. An optional 'qualification_code' column "
        "(NWIS peak codes) is applied as peakfq does unless --ignore-peak-codes."
    ),
)
@click.option("--site", "site_name", default="", help="Site name/number for the report title.")
@click.option(
    "--regional-skew",
    type=float,
    default=None,
    help=(
        "Regional skew coefficient from a published study; give --regional-skew-se with it. "
        "There is no default: pass this, --station-skew, or --default-skew."
    ),
)
@click.option("--regional-skew-se", type=float, default=None, help="Regional skew standard error.")
@click.option(
    "--station-skew",
    "station_skew_only",
    is_flag=True,
    default=False,
    help="Use the at-site skew alone, with no regional weighting.",
)
@click.option(
    "--default-skew",
    "use_default_skew",
    is_flag=True,
    default=False,
    help=(
        "Explicitly accept the unsourced fallback regional skew (-0.302, SE 0.55). "
        "Bulletin 17C gives no national default; prefer a published regional study."
    ),
)
@click.option(
    "--low-outlier-threshold",
    type=float,
    default=None,
    help="User-supplied PILF threshold in cfs. Omit to let MGBT decide, as both engines do by default.",
)
@click.option(
    "--historical",
    "historical_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "CSV of historical (non-systematic) peaks, with 'water_year' and 'peak_flow_cfs' "
        "columns. Use with --threshold to say what flood would have been recorded."
    ),
)
@click.option(
    "--threshold",
    "thresholds",
    type=(int, int, float),
    multiple=True,
    metavar="START END LOWER_CFS",
    help=(
        "Perception threshold: water years START..END (inclusive), where any peak above "
        "LOWER_CFS would have been recorded. Repeatable. Years in the period with no peak "
        "are censored below LOWER_CFS."
    ),
)
@click.option(
    "--ignore-peak-codes",
    is_flag=True,
    default=False,
    help=(
        "Fit every peak as an exact systematic value, ignoring the CSV's "
        "'qualification_code' column (the behaviour before codes were applied)."
    ),
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
    regional_skew_se: Optional[float],
    station_skew_only: bool,
    use_default_skew: bool,
    low_outlier_threshold: Optional[float],
    historical_path: Optional[Path],
    thresholds: Tuple[Tuple[int, int, float], ...],
    ignore_peak_codes: bool,
    tolerance_pct: float,
    output_path: Optional[Path],
) -> None:
    """Compare the native EMA against the vendored USGS Fortran (peakfq 8.1.0).

    Requires the f2py extension to be built (``python build_fortran/build.py``,
    needs gfortran and meson) -- this is the one thing this command cannot do
    without, and it says so rather than silently falling back to anything, per
    ``docs/FORTRAN_ENGINE_DESIGN.md`` section 9.

    Historical peaks (``--historical``) and perception thresholds
    (``--threshold``) are passed to both engines unchanged, as
    :func:`flowfreq.workflow.compare_engines` takes them.

    A ``qualification_code`` column in the ``--peaks`` CSV is applied to both
    engines alike, as peakfq 8.1.0's ``siteQT`` applies it (code 7 historic;
    3/O/6/C removed, 4/8 a less-than/greater-than interval). A record neither
    engine can express exactly is an error naming its years;
    ``--ignore-peak-codes`` fits every peak as exact.
    """
    import pandas as pd

    from flowfreq.workflow import compare_engines, resolve_regional_skew

    # Codes as strings: an all-numeric code column ("7") is otherwise read as
    # float and comes back as "7.0".
    peaks_df = pd.read_csv(peaks_path, dtype={"qualification_code": str})
    missing = {"water_year", "peak_flow_cfs"} - set(peaks_df.columns)
    if missing:
        raise click.UsageError(
            f"--peaks CSV is missing column(s) {sorted(missing)}; expected 'water_year' and "
            "'peak_flow_cfs' (the shape USGSgage.download_peak_flow produces)."
        )

    historical_peaks: Optional[List[Tuple[int, float]]] = None
    if historical_path is not None:
        hist_df = pd.read_csv(historical_path)
        missing = {"water_year", "peak_flow_cfs"} - set(hist_df.columns)
        if missing:
            raise click.UsageError(
                f"--historical CSV is missing column(s) {sorted(missing)}; expected "
                "'water_year' and 'peak_flow_cfs'."
            )
        historical_peaks = [
            (int(y), float(q)) for y, q in zip(hist_df["water_year"], hist_df["peak_flow_cfs"])
        ]

    perception_thresholds: Optional[Dict[Tuple[int, int], float]] = None
    if thresholds:
        perception_thresholds = {}
        for start, end, lower in thresholds:
            if start > end:
                raise click.UsageError(f"--threshold start {start} is after end {end}.")
            if lower <= 0:
                raise click.UsageError(f"--threshold lower bound must be positive, got {lower}.")
            perception_thresholds[(start, end)] = lower

    # Settle the skew before any fitting, so a missing choice is a usage error
    # rather than a traceback from deep inside the comparison.
    try:
        resolve_regional_skew(regional_skew, regional_skew_se, use_default_skew, station_skew_only)
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc

    peak_codes: Optional[List[object]] = None
    if "qualification_code" in peaks_df.columns:
        peak_codes = peaks_df["qualification_code"].tolist()

    from flowfreq.psf_convert import UnsupportedSpecError

    try:
        report = compare_engines(
            peak_flows=peaks_df["peak_flow_cfs"].to_numpy(dtype=float),
            water_years=peaks_df["water_year"].to_numpy(dtype=int),
            regional_skew=regional_skew,
            regional_skew_se=regional_skew_se,
            use_default_skew=use_default_skew,
            station_skew_only=station_skew_only,
            user_low_outlier_threshold=low_outlier_threshold,
            historical_peaks=historical_peaks,
            perception_thresholds=perception_thresholds,
            site_name=site_name,
            tolerance_pct=tolerance_pct,
            peak_codes=peak_codes,
            apply_peak_codes=not ignore_peak_codes,
        )
    except ImportError as exc:
        raise click.ClickException(str(exc)) from exc
    except UnsupportedSpecError as exc:
        raise click.ClickException(
            f"{exc} (pass --ignore-peak-codes to fit every peak as exact)"
        ) from exc
    except ValueError as exc:
        # Codes peakfq acts on combined with --historical/--threshold, or a
        # record the engines reject: a clean message, not a traceback.
        raise click.ClickException(str(exc)) from exc

    markdown = report.to_markdown()
    if output_path is not None:
        output_path.write_text(markdown, encoding="utf-8")
        click.echo(f"Wrote {output_path}")
    else:
        click.echo(markdown)

    if not report.comparison.passed:
        raise SystemExit(1)


if __name__ == "__main__":
    cli()
