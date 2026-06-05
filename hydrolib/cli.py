"""
HydroLib command-line interface.

Provides CLI commands for validation and benchmarking of the
Bulletin 17C implementation.
"""

from __future__ import annotations

import click


@click.group()
def cli() -> None:
    """HydroLib - Hydrologic frequency analysis tools."""
    pass


@cli.command()
def validate() -> None:
    """Run validation benchmarks against expected PeakfqSA results."""
    from hydrolib.validation.benchmarks import print_benchmark_report, run_all_benchmarks

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
    from hydrolib.validation.benchmarks import run_all_benchmarks
    from hydrolib.validation.reports import generate_json_report, generate_text_report

    click.echo("Running benchmarks...")
    results = run_all_benchmarks()

    if fmt == "json":
        click.echo(generate_json_report(results))
    else:
        click.echo(generate_text_report(results))


if __name__ == "__main__":
    cli()
