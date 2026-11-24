"""
FlowFreq command-line interface.

Provides CLI commands for validation and benchmarking of the
Bulletin 17C implementation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

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


def _code_values(chars: Tuple[str, ...]) -> Dict[str, float]:
    """``CODE=VALUE`` options as a dict; :class:`click.BadParameter` for anything else."""
    values: Dict[str, float] = {}
    for item in chars:
        name, _, raw = item.partition("=")
        try:
            values[name.strip()] = float(raw)
        except ValueError:
            raise click.BadParameter(f"{item!r} is not CODE=VALUE", param_hint="--char")
        if not name.strip():
            raise click.BadParameter(f"{item!r} is not CODE=VALUE", param_hint="--char")
    return values


@cli.command()
@click.option("--lat", type=float, required=True, help="Site latitude, decimal degrees.")
@click.option("--lon", type=float, required=True, help="Site longitude, decimal degrees.")
@click.option("--state", required=True, help="Two-letter jurisdiction code (e.g. MT).")
@click.option(
    "--aep",
    "aeps",
    type=click.FloatRange(0.0, 1.0, min_open=True, max_open=True),
    multiple=True,
    help=(
        "Annual exceedance probability (fraction); repeat for several. Default: every AEP "
        "that all the regression regions used publish."
    ),
)
@click.option(
    "--allow-nss-snapshot",
    is_flag=True,
    help="Accept equations from the unreviewed NSS snapshot tier.",
)
@click.option(
    "--allow-extrapolation",
    is_flag=True,
    help="Allow basin characteristics outside the equations' calibrated range.",
)
@click.option(
    "--region",
    default=None,
    help=(
        "Where NSS returns overlapping regions (alternative equation sets for the same "
        "area), the one to use: an NSS code (e.g. GC2016) or library region code. Give "
        'several, comma-separated ("GC1, GC2"), to use one set per side of a region '
        "boundary; together they must cover the whole basin. A set NSS does not serve "
        '("Statewide OmegaEM" in TX; "Urban 2026" in AL, GA, MS, NC, SC and TN) is used for '
        "the whole basin when named."
    ),
)
@click.option(
    "--char",
    "chars",
    multiple=True,
    help="A basin characteristic StreamStats does not return, as CODE=VALUE; repeat for each.",
)
def ungaged(
    lat: float,
    lon: float,
    state: str,
    aeps: Tuple[float, ...],
    allow_nss_snapshot: bool,
    allow_extrapolation: bool,
    region: Optional[str],
    chars: Tuple[str, ...],
) -> None:
    """Estimate flood-frequency quantiles at an ungaged site (regional regression).

    Delineates the basin with StreamStats, evaluates the state's offline
    regression equations and prints the estimate with its provenance as
    Markdown. Exits 1 with the message when no equations are available, the
    point cannot be snapped to a stream, or the region is unsupported. When NSS
    returns overlapping regions (alternative equation sets for one area), lists
    them one per line and exits 1 until --region chooses one.
    """
    import requests

    from flowfreq.regression.equations import OutOfRangeError
    from flowfreq.ungaged import estimate_at_site

    supplied = _code_values(chars)
    try:
        estimate = estimate_at_site(
            lat,
            lon,
            state=state,
            aeps=list(aeps) if aeps else None,
            allow_nss_snapshot=allow_nss_snapshot,
            allow_extrapolation=allow_extrapolation,
            region=(
                [c.strip() for c in region.split(",")]
                if region is not None and "," in region
                else region
            ),
            supplied_characteristics=supplied or None,
        )
    except OutOfRangeError as exc:
        raise click.ClickException(f"{exc} (pass --allow-extrapolation to proceed)") from exc
    except KeyError as exc:  # unknown --state: KeyError quotes str(); use the message
        raise click.ClickException(str(exc.args[0]) if exc.args else "unknown key") from exc
    except (ValueError, LookupError, requests.RequestException) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(estimate.to_markdown())


@cli.command("ungaged-chars")
@click.option("--state", required=True, help="Two-letter jurisdiction code (e.g. FL).")
@click.option(
    "--region",
    required=True,
    help=(
        "The region from the report's map: a library region code or NSS code. For a basin "
        'spanning regions, give each with its percent, comma-separated ("1=60,2=40").'
    ),
)
@click.option(
    "--char",
    "chars",
    multiple=True,
    required=True,
    help="A basin characteristic as CODE=VALUE (e.g. DRNAREA=105); repeat for each.",
)
@click.option("--lat", type=float, required=True, help="Site latitude, decimal degrees.")
@click.option("--lon", type=float, required=True, help="Site longitude, decimal degrees.")
@click.option(
    "--aep",
    "aeps",
    type=click.FloatRange(0.0, 1.0, min_open=True, max_open=True),
    multiple=True,
    help="Annual exceedance probability (fraction); repeat for several. Default: every AEP.",
)
@click.option("--allow-nss-snapshot", is_flag=True, help="Accept the unreviewed NSS snapshot tier.")
@click.option(
    "--allow-extrapolation",
    is_flag=True,
    help="Allow basin characteristics outside the equations' calibrated range.",
)
def ungaged_chars(
    state: str,
    region: str,
    chars: Tuple[str, ...],
    lat: float,
    lon: float,
    aeps: Tuple[float, ...],
    allow_nss_snapshot: bool,
    allow_extrapolation: bool,
) -> None:
    """Estimate flood quantiles from basin characteristics you supply (no StreamStats).

    For sites StreamStats can't delineate, such as Florida, Michigan and Nebraska. Name the
    region from the report's map and give each characteristic the region's equations use,
    measured as the report defines it. Prints the estimate as Markdown.
    """
    from flowfreq.regression.equations import OutOfRangeError
    from flowfreq.ungaged import estimate_from_characteristics

    values = _code_values(chars)
    if "=" in region:
        try:
            chosen: Union[str, Dict[str, float]] = {
                k.strip(): float(v) for k, v in (p.split("=") for p in region.split(","))
            }
        except ValueError:
            raise click.BadParameter(f"{region!r} is not REGION=PERCENT,...", param_hint="--region")
    else:
        chosen = region
    try:
        estimate = estimate_from_characteristics(
            state,
            values,
            region=chosen,
            lat=lat,
            lon=lon,
            aeps=list(aeps) if aeps else None,
            allow_nss_snapshot=allow_nss_snapshot,
            allow_extrapolation=allow_extrapolation,
        )
    except OutOfRangeError as exc:
        raise click.ClickException(f"{exc} (pass --allow-extrapolation to proceed)") from exc
    except KeyError as exc:
        raise click.ClickException(str(exc.args[0]) if exc.args else "unknown key") from exc
    except (ValueError, LookupError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(estimate.to_markdown())


@cli.command()
@click.argument("spec", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    required=True,
    help="The report to write: .docx (needs the 'report' extra) or .md.",
)
@click.option(
    "--save-spec",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Also write the specification as built, with downloaded peaks and looked-up skews "
        "recorded in it, so the report rebuilds offline."
    ),
)
def report(spec: Path, output: Path, save_spec: Optional[Path]) -> None:
    """Build a report from a JSON specification.

    SPEC lists the report's sections (flowfreq.reporting.SECTION_KINDS) and their inputs; every
    table and figure is computed from them. See flowfreq.reporting.
    """
    from flowfreq.reporting import ReportSpec, build_report

    # Before the build, which may download records: a suffix the writer will refuse.
    if output.suffix.lower() not in (".docx", ".md"):
        raise click.ClickException(
            f"cannot write a report as {output.suffix or 'a file with no suffix'}: use .docx or .md"
        )
    try:
        built = build_report(ReportSpec.load(spec))
        built.save(output)
    except ImportError as exc:
        raise click.ClickException(str(exc)) from exc
    except KeyError as exc:
        raise click.ClickException(str(exc.args[0]) if exc.args else "unknown key") from exc
    except (ValueError, LookupError) as exc:
        raise click.ClickException(str(exc)) from exc
    if save_spec is not None:
        built.spec.save(save_spec)
        click.echo(f"Wrote {output} and {save_spec}")
    else:
        click.echo(f"Wrote {output}")


def _parse_months(ctx: click.Context, param: click.Parameter, value: str) -> Tuple[int, ...]:
    """Click callback: ``"11,12,1"`` to a tuple of months 1-12."""
    try:
        months = tuple(int(m) for m in value.split(","))
    except ValueError as exc:
        raise click.BadParameter(f"{value!r} is not a comma list of integers") from exc
    bad = [m for m in months if not 1 <= m <= 12]
    if bad:
        raise click.BadParameter(f"months must be integers 1-12, got {bad}")
    return months


def _parse_day(ctx: click.Context, param: click.Parameter, value: str) -> str:
    """Click callback: a ``YYYY-MM-DD`` date, returned as given."""
    import pandas as pd

    try:
        pd.to_datetime(value, format="%Y-%m-%d")
    except ValueError as exc:
        raise click.BadParameter(f"{value!r} is not a YYYY-MM-DD date") from exc
    return value


@cli.command("ice-events")
@click.option("--gauge", "gauges", multiple=True, required=True, help="USGS site number (repeat).")
@click.option(
    "--air-station", required=True, help="NCEI air-temperature station, e.g. USC00123456."
)
@click.option("--start", required=True, callback=_parse_day, help="First day, YYYY-MM-DD.")
@click.option("--end", required=True, callback=_parse_day, help="Last day, YYYY-MM-DD.")
@click.option(
    "--months",
    required=True,
    callback=_parse_months,
    help="Analysis months, comma list of 1-12 (e.g. 11,12,1,2,3).",
)
@click.option("--air-c", required=True, type=float, help="Cold when TMIN is at or below this, C.")
@click.option(
    "--lookback-days",
    required=True,
    type=click.IntRange(min=0),
    help="Days before a day whose cold air still counts.",
)
@click.option(
    "--merge-gap-days",
    required=True,
    type=click.IntRange(min=0),
    help="Missing days allowed inside one event.",
)
@click.option(
    "--min-cold-spell-days",
    required=True,
    type=click.IntRange(min=1),
    help="Shortest cold spell without gauge ice reported in the QA/QC file.",
)
@click.option(
    "--out",
    "out_dir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="Output directory, created if missing.",
)
def ice_events(
    gauges: Tuple[str, ...],
    air_station: str,
    start: str,
    end: str,
    months: Tuple[int, ...],
    air_c: float,
    lookback_days: int,
    merge_gap_days: int,
    min_cold_spell_days: int,
    out_dir: Path,
) -> None:
    """Gauge ice events from daily qualifiers and cold air.

    A gauge day with the ICE qualifier AND cold air (NCEI TMIN at or below --air-c
    on the day or within --lookback-days before it) is an event day; days within
    --merge-gap-days of each other join one event, counted in the winter and month of
    its first day. Writes air_<station>.csv, daily_<site>.csv, gauge_ice_events.csv,
    gauge_ice_qaqc.csv and ice_events_run.json into --out. Every threshold is an option.
    """
    import json
    from datetime import datetime, timezone

    import pandas as pd
    import requests

    import flowfreq
    from flowfreq import ghcn, ice, waterdata

    try:
        reversed_dates = pd.Timestamp(start) > pd.Timestamp(end)
    except ValueError:
        reversed_dates = False  # a malformed date is reported by the fetch below
    if reversed_dates:
        raise click.ClickException(f"--start {start} is after --end {end}")
    try:
        air = ghcn.fetch_ghcn_daily(air_station, start, end)
    except (requests.RequestException, ValueError) as exc:
        raise click.ClickException(f"air station {air_station}: {exc}") from exc
    air_tmin_days = int(air["tmin_c"].notna().sum())
    if air_tmin_days == 0:
        raise click.ClickException(
            f"air station {air_station}: no daily minimum temperature between {start} and {end}"
        )
    cold = ice.cold_air_days(air["tmin_c"], air_c, lookback_days)
    cold = cold.reindex(air.index).fillna(False).astype(bool)

    # Every gauge is fetched before anything is written, so a failure leaves no partial output.
    dailies: Dict[str, "pd.DataFrame"] = {}
    for site in gauges:
        try:
            dailies[site] = waterdata.download_daily(
                site, start_date=start, end_date=end, keep_null=True
            )
        except (requests.RequestException, ValueError, LookupError) as exc:
            raise click.ClickException(f"gauge {site}: {exc}") from exc

    out_dir.mkdir(parents=True, exist_ok=True)
    air.to_csv(out_dir / f"air_{air_station}.csv")
    in_months = pd.DatetimeIndex(air.index).month.isin(months)
    event_rows: List["pd.DataFrame"] = []
    qaqc_rows: List["pd.DataFrame"] = []
    per_gauge: Dict[str, Dict[str, int]] = {}
    for site, daily in dailies.items():
        flagged = ice.ice_days_from_qualifiers(daily)
        day_ice = pd.Series(flagged["ice"].to_numpy(), index=pd.DatetimeIndex(flagged["date"]))
        day_ice = day_ice[~day_ice.index.duplicated(keep="last")]
        has_record = pd.Series(True, index=day_ice.index).reindex(air.index, fill_value=False)
        day_ice = day_ice.reindex(air.index).fillna(False).astype(bool)
        recorded = in_months & has_record.to_numpy(dtype=bool)

        flow = daily[["flow_cfs", "qualification_code"]].copy()
        flow.index = pd.DatetimeIndex(flagged["date"]).strftime("%Y-%m-%d")
        flow.index.name = "date"
        flow.to_csv(out_dir / f"daily_{site}.csv")

        both = pd.Series(air.index[(day_ice & cold).to_numpy()])
        event_rows.append(ice.events_in_months(both, months, merge_gap_days).assign(site=site))
        qa = ice.gauge_ice_qaqc(
            day_ice[recorded], cold[recorded], min_cold_spell_days, merge_gap_days
        )
        qaqc_rows.append(qa.assign(site=site))
        per_gauge[site] = {"daily_rows": int(len(daily)), "ice_days": int(flagged["ice"].sum())}

    event_cols = ["site", "winter", "month", "start", "end", "n_days"]
    qaqc_cols = ["site", "kind", "start", "end", "n_days"]
    events = pd.concat(event_rows, ignore_index=True).reindex(columns=event_cols)
    qaqc = pd.concat(qaqc_rows, ignore_index=True).reindex(columns=qaqc_cols)
    for frame in (events, qaqc):
        for col in ("start", "end"):
            frame[col] = pd.to_datetime(frame[col]).dt.strftime("%Y-%m-%d")
    events.to_csv(out_dir / "gauge_ice_events.csv", index=False)
    qaqc.to_csv(out_dir / "gauge_ice_qaqc.csv", index=False)

    run = {
        "arguments": {
            "gauge": list(gauges),
            "air_station": air_station,
            "start": start,
            "end": end,
            "months": list(months),
            "air_c": air_c,
            "lookback_days": lookback_days,
            "merge_gap_days": merge_gap_days,
            "min_cold_spell_days": min_cold_spell_days,
            "out": str(out_dir),
        },
        "flowfreq_version": flowfreq.__version__,
        "utc_time": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "air_tmin_days": air_tmin_days,
        "gauges": per_gauge,
    }
    (out_dir / "ice_events_run.json").write_text(json.dumps(run, indent=2), encoding="utf-8")
    click.echo(f"Wrote {len(events)} ice events and {len(qaqc)} QA/QC rows to {out_dir}")


if __name__ == "__main__":
    cli()
