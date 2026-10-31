"""Nationwide smoke run: one gaged fit and one ungaged estimate per jurisdiction.

For each of the 56 jurisdictions, picks one catalog gage (the longest record among
reference gages, else any unregulated one), fits Bulletin 17C at the gage with
station skew, then asks :func:`flowfreq.estimate_at_site` for the regression
estimate at the gage coordinates. Results go to a Markdown table. Failures are
recorded, never fixed here: each is a finding.

Per-state results are cached as JSON under ``--cache`` (default
``.cache/nationwide_smoke``, gitignored), so ``--resume`` skips finished states.

Usage::

    python tools/nationwide_smoke.py --out docs/NATIONWIDE_SMOKE.md --resume
    python tools/nationwide_smoke.py --points points.csv --states VT

``--points`` takes a CSV with header ``state,lat,lon,note``: a hand-picked pour point
per state for the ungaged estimate (StreamStats exposes no snap distance, so gage
coordinates that do not snap need a manual point). An optional ``region`` column
passes ``region=`` to :func:`flowfreq.estimate_at_site` (one code, or several separated
by commas for a boundary basin); a row with a region may leave ``lat`` and ``lon``
blank to keep the gage coordinates.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import sys
import time
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

logger = logging.getLogger("nationwide_smoke")

KNOWN_GAPS = (
    "DC, VI, GU, MP and AS have no regression-equation library, so their ungaged "
    "estimates are expected to fail with `EquationsUnavailable`.",
    "StreamStats does not delineate in the territories and only partly in AK and HI, "
    "so those ungaged estimates may fail at delineation.",
    "StreamStats' snap service has no stream grid for FL, MI and NE (HTTP 404, "
    "`RegionNotDelineatedError`), so no point there can be delineated through its API.",
    "A state with no unregulated catalog gage has no row to fit; it is reported as such.",
)

# Reference gages first, then unclassified, then urban. Regulated never qualifies.
_CLASS_ORDER = {"reference": 0, "unknown": 1, "urban": 2}


def select_site(catalog: pd.DataFrame, state: str) -> Optional[pd.Series]:
    """Pick the smoke-test gage for ``state``.

    The longest record (``n_peaks``) among ``reference`` gages; failing that,
    among unclassified ones, then urban ones. ``regulated`` gages are never
    chosen.

    Returns
    -------
    pandas.Series or None
        The catalog row, or ``None`` when the state has no such gage.
    """
    sub = catalog[catalog["state"] == state]
    sub = sub[sub["regulation_class"].isin(_CLASS_ORDER)]
    if sub.empty:
        return None
    rank = sub["regulation_class"].map(_CLASS_ORDER)
    sub = sub[rank == rank.min()]
    return sub.sort_values(["n_peaks", "site_no"], ascending=[False, True]).iloc[0]


def _err(exc: BaseException) -> str:
    return f"{type(exc).__name__}: {' '.join(str(exc).split())[:300]}"


def load_points(path: Path) -> Dict[str, Dict[str, Any]]:
    """Read hand-picked ungaged pour points from a CSV with header ``state,lat,lon,note``.

    Returns
    -------
    dict
        Uppercased state code to ``{"lat": float, "lon": float, "note": str}``, plus
        ``"region"`` (a code, or a list of codes) when the optional ``region`` column is
        set. With a region, blank ``lat`` and ``lon`` give ``None`` (keep the gage).

    Raises
    ------
    ValueError
        If a required column is missing, a coordinate is non-numeric or out of range,
        or a state appears twice.
    """
    points: Dict[str, Dict[str, Any]] = {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        fields = [f.strip() for f in (reader.fieldnames or [])]
        for col in ("state", "lat", "lon"):
            if col not in fields:
                raise ValueError(f"points file {path} is missing required column: {col}")
        for raw in reader:
            r = {(k or "").strip(): (v or "").strip() for k, v in raw.items() if k is not None}
            state = r["state"].upper()
            codes = [c.strip() for c in r.get("region", "").split(",") if c.strip()]
            entry: Dict[str, Any] = {"note": r.get("note", "")}
            if codes and not r["lat"] and not r["lon"]:
                entry.update(lat=None, lon=None)
            else:
                try:
                    lat = float(r["lat"])
                    lon = float(r["lon"])
                except ValueError as exc:
                    raise ValueError(f"{state}: non-numeric coordinate ({exc})") from exc
                if not (-90.0 <= lat <= 90.0) or not (-180.0 <= lon <= 180.0):
                    raise ValueError(f"{state}: coordinate out of range (lat={lat}, lon={lon})")
                entry.update(lat=lat, lon=lon)
            if codes:
                entry["region"] = codes[0] if len(codes) == 1 else codes
            if state in points:
                raise ValueError(f"duplicate state in points file: {state}")
            points[state] = {"lat": entry.pop("lat"), "lon": entry.pop("lon"), **entry}
    return points


def smoke_one(
    row: pd.Series,
    *,
    allow_nss_snapshot: bool = True,
    point: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Run the gaged fit and the ungaged estimate for one catalog gage.

    With ``point`` (``{"lat", "lon", "note"}``) the ungaged estimate uses that pour
    point instead of the gage coordinates; the gaged fit is unchanged. A ``"region"``
    entry is passed as ``region=``, and ``None`` coordinates keep the gage's.

    Never raises: each failure is recorded as ``"gaged: Type: message"`` or
    ``"ungaged: Type: message"`` under ``error``.
    """
    from flowfreq import estimate_at_site
    from flowfreq.usgs import USGSgage
    from flowfreq.workflow import run_ffa

    moved = bool(point) and point.get("lat") is not None
    region = point.get("region") if point else None
    where = f"override: {point.get('note', '')}" if moved and point else "gage"
    if region is not None:
        where += f", region={region if isinstance(region, str) else ', '.join(region)}"
    out: Dict[str, Any] = {
        "state": row["state"],
        "ungaged_point": where,
        "site_no": row["site_no"],
        "n_peaks": None if pd.isna(row.get("n_peaks")) else int(row["n_peaks"]),
        "ffa_ok": False,
        "ffa_q1pct": None,
        "ungaged_ok": False,
        "ungaged_q1pct": None,
        "library_status": None,
        "error": None,
    }
    errors: List[str] = []
    try:
        peaks = USGSgage(row["site_no"]).download_peak_flow()
        peaks = peaks.dropna(subset=["peak_flow_cfs"])
        res = run_ffa(
            peaks["peak_flow_cfs"].to_numpy(float),
            peaks["water_year"].to_numpy(int),
            station_skew_only=True,
            peak_codes=peaks["qualification_code"].tolist(),
            site_no=row["site_no"],
            allow_regulated=False,
        )
        if res.get("error"):
            raise RuntimeError(res["error"])
        q = res["quantile_df"]
        out["ffa_q1pct"] = float(q.loc[(q["AEP (%)"] - 0.01).abs().idxmin(), "Flow (cfs)"])
        out["ffa_ok"] = True
    except Exception as exc:  # noqa: BLE001 - a smoke run records every failure
        logger.warning("%s gaged %s: %s", row["state"], row["site_no"], _err(exc))
        errors.append("gaged: " + _err(exc))
    try:
        est = estimate_at_site(
            float(point["lat"]) if moved and point else float(row["latitude"]),
            float(point["lon"]) if moved and point else float(row["longitude"]),
            state=row["state"],
            allow_nss_snapshot=allow_nss_snapshot,
            region=region,
        )
        qd = est.quantiles
        out["ungaged_q1pct"] = float(qd.loc[(qd["aep"] - 0.01).abs().idxmin(), "flow_cfs"])
        out["library_status"] = est.library_status
        out["ungaged_ok"] = True
    except Exception as exc:  # noqa: BLE001
        logger.warning("%s ungaged %s: %s", row["state"], row["site_no"], _err(exc))
        errors.append("ungaged: " + _err(exc))
    out["error"] = "; ".join(errors) or None
    return out


def _fmt(v: Any) -> str:
    return "" if v is None or (isinstance(v, float) and v != v) else f"{v:,.0f}"


def _kinds(error: Optional[str]) -> List[str]:
    """Turn ``'gaged: A: x; ungaged: B: y'`` into ``['gaged A', 'ungaged B']``."""
    return [
        f"{side} {name}"
        for side, name in re.findall(r"(?:^|; )(gaged|ungaged): (\w+):", error or "")
    ]


def write_report(
    rows: Sequence[Dict[str, Any]], path: Path, *, run_date: Optional[str] = None
) -> None:
    """Write the per-jurisdiction table, known gaps and failure summary as Markdown."""
    run_date = run_date or date.today().isoformat()
    n = len(rows)
    n_ffa = sum(bool(r["ffa_ok"]) for r in rows)
    n_ung = sum(bool(r["ungaged_ok"]) for r in rows)
    lines = [
        "# Nationwide smoke run",
        "",
        f"Run date: {run_date}. One catalog gage per jurisdiction (longest record, "
        "reference class first, else unregulated); a Bulletin 17C fit with station skew "
        "at the gage, and a regression estimate (NSS snapshot allowed) at the gage "
        "coordinates. Q1% is the 1% AEP flow in cfs. Generated by "
        "`tools/nationwide_smoke.py`.",
        "",
        "## Known gaps",
        "",
        *[f"- {g}" for g in KNOWN_GAPS],
        "",
        "## Summary",
        "",
        f"- Gaged fits succeeded: {n_ffa} of {n}",
        f"- Ungaged estimates succeeded: {n_ung} of {n}",
        "",
        "Failures by kind:",
        "",
    ]
    groups: Dict[str, List[str]] = defaultdict(list)
    for r in rows:
        for k in _kinds(r.get("error")):
            groups[k].append(r["state"])
    if groups:
        lines += ["| Kind | Count | States |", "|---|---|---|"]
        for k, states in sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            lines.append(f"| {k} | {len(states)} | {', '.join(sorted(states))} |")
    else:
        lines.append("None.")
    lines += [
        "",
        "## Results",
        "",
        "| State | Gage | Years | Gaged | Q1% gaged | Ungaged | Q1% ungaged | Library | Point | Error |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        err = (r.get("error") or "").replace("|", "\\|")
        lines.append(
            f"| {r['state']} | {r.get('site_no') or ''} | {r.get('n_peaks') or ''} | "
            f"{'ok' if r['ffa_ok'] else 'FAIL'} | {_fmt(r['ffa_q1pct'])} | "
            f"{'ok' if r['ungaged_ok'] else 'FAIL'} | {_fmt(r['ungaged_q1pct'])} | "
            f"{r.get('library_status') or ''} | {r.get('ungaged_point') or 'gage'} | {err} |"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=REPO_ROOT / "docs" / "NATIONWIDE_SMOKE.md")
    ap.add_argument("--cache", type=Path, default=REPO_ROOT / ".cache" / "nationwide_smoke")
    ap.add_argument("--resume", action="store_true", help="reuse cached per-state results")
    ap.add_argument("--states", nargs="*", help="limit to these state codes")
    ap.add_argument(
        "--points",
        type=Path,
        default=None,
        help="CSV (state,lat,lon,note) of hand-picked ungaged pour points per state",
    )
    args = ap.parse_args(argv)
    points = load_points(args.points) if args.points else {}
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

    from flowfreq.catalog import load_catalog
    from flowfreq.regression.jurisdictions import JURISDICTIONS

    catalog = load_catalog()
    args.cache.mkdir(parents=True, exist_ok=True)
    rows: List[Dict[str, Any]] = []
    for j in JURISDICTIONS:
        if args.states and j.code not in args.states:
            continue
        cache = args.cache / f"{j.code}.json"
        if args.resume and cache.exists():
            rows.append(json.loads(cache.read_text(encoding="utf-8")))
            continue
        site = select_site(catalog, j.code)
        t0 = time.time()
        if site is None:
            row: Dict[str, Any] = {
                "state": j.code,
                "site_no": None,
                "n_peaks": None,
                "ffa_ok": False,
                "ffa_q1pct": None,
                "ungaged_ok": False,
                "ungaged_q1pct": None,
                "library_status": None,
                "error": "gaged: NoSite: no unregulated catalog gage; "
                "ungaged: NoSite: no gage coordinates",
            }
        else:
            row = smoke_one(site, point=points.get(j.code))
        logger.info("%s done in %.0f s: %s", j.code, time.time() - t0, row["error"] or "ok")
        cache.write_text(json.dumps(row), encoding="utf-8")
        rows.append(row)
    write_report(rows, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
