"""Snapshot USGS NSS regression equations to versioned JSON, with a diff report.

NSS scenario templates (``GET /nssservices/regions/{r}/Scenarios``) carry each
regression region's parameters and limits but **no equation**. Checked
2026-09-25: neither ``/nssservices/apiconfig`` nor the ``RegressionRegions``
resource (``/regressionregions/{id}``, ``/{id}/limitations``) exposes one
either. The literal equation string appears only in an Estimate result. So for
each state and statistic group this tool:

1. GETs the templates;
2. fills every parameter with a synthetic in-range value, the geometric
   midpoint of its limits, or the arithmetic midpoint when ``min <= 0``;
3. POSTs them to ``/nssservices/Scenarios/Estimate``;
4. harvests each result's ``equation``, ``sep``, ``errors`` (ASEp),
   ``equivalentYears``, statistic code, region id/code/name, citation and the
   parameter limits and units;
5. parses the equation with :func:`flowfreq.regression.nss.parse_equation` and
   round-trips it (see below).

**Gating parameters.** Some parameters have no limits and act as region
selectors. Oregon's ``ORREG2`` must equal the region's GC number (729) before
Estimate answers at all. A region that Estimate drops, or returns with no
results, is retried with each unbounded parameter set to the region's GC
number, 0 and 1. A region that still answers nothing is recorded as
``unresolved`` with the reason, never dropped.

**Round trip.** Each parsed equation is evaluated with
:func:`flowfreq.regression.equations.evaluate` at the submitted inputs and
compared with two things:

- a direct evaluation of the NSS string, which must agree to 1e-9 relative;
- NSS's own returned value. NSS rounds every value it returns to three
  significant figures, so the check is that the two agree within half a unit
  in the third significant figure.

**What a snapshot is not.** NSS supplies no report table number, covariance
matrix, model-error variance or site count, and its ``sep`` field does not
behave like the report's log10 SEP. Montana's NW region, for example, returns
0.043 at 4% AEP against 0.41 at 66.7%. It is therefore stored verbatim as
``nss_sep`` and never mapped to ``sep_log``. A snapshot is a cross-check and a
starting point for ``flowfreq/data/regression/<STATE>.json``, not a verified
library, and this tool never writes there.

Output goes to ``data/nss_snapshots/<STATE>_<YYYY-MM-DD>.json`` at the
repository root. That location is outside the package on purpose: the files
are not shipped, and are never loaded as equations. Each run prints a diff
against the newest earlier snapshot for the same state, so an upstream change
to an equation shows up in review.

Usage::

    python tools/snapshot_nss.py --state WA --state ID
    python tools/snapshot_nss.py --state OR --group all --out /tmp/nss

Roadmap: ``docs/MASTER_ROADMAP.md`` §3.1, issue #35.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import requests

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from flowfreq.regression.equations import Citation, OutOfRangeError, evaluate  # noqa: E402
from flowfreq.regression.nss import (  # noqa: E402
    NSSEquationError,
    evaluate_expression,
    parse_equation,
    statistic_code_to_aep,
    to_regression_equation,
)
from flowfreq.streamstats import (  # noqa: E402
    GENERIC_HOST,
    MAX_CONCURRENCY,
    StreamStatsResponseError,
    StreamStatsTransportError,
    _fetch_citations,
    _request_with_backoff,
    list_statistic_groups,
)

SNAPSHOT_SCHEMA = "flowfreq.nss_snapshot/1"
DEFAULT_OUT = REPO_ROOT / "data" / "nss_snapshots"
UNIT_SYSTEM = 2  # US customary
#: A parsed equation must reproduce the direct evaluation of its NSS string to this.
PARSE_RTOL = 1e-9
#: NSS rounds returned values to this many significant figures.
NSS_SIG_FIGS = 3
#: Values tried, per unbounded parameter, for a region Estimate did not answer.
GATING_FALLBACKS = (0.0, 1.0)
UNBOUNDED_DEFAULT = 1.0

NSS_SEP_NOTE = (
    "nss_sep is NSS's raw 'sep' field, stored verbatim. It is not confirmed to be the "
    "report's log10 standard error of prediction (MT NW region: 0.043 at 4% AEP vs "
    "0.41 at 66.7%), so it is never mapped to sep_log."
)
SNAPSHOT_NOTE = (
    "Harvested from NSS Scenarios/Estimate at synthetic in-range inputs. NSS supplies "
    "no report table, covariance, model-error variance or site count, so this is a "
    "cross-check and starting point, not a verified equation library."
)


# ---------------------------------------------------------------------------
# Synthetic inputs
# ---------------------------------------------------------------------------


def synthetic_value(limits: Optional[Dict[str, Any]]) -> Optional[float]:
    """An in-range value for a parameter, or None when it has no limits.

    The geometric midpoint when ``min > 0``, else the arithmetic midpoint.
    With only one bound, a value just inside it.
    """
    limits = limits or {}
    lo, hi = limits.get("min"), limits.get("max")
    if lo is not None and hi is not None:
        lo, hi = float(lo), float(hi)
        return math.sqrt(lo * hi) if lo > 0 else (lo + hi) / 2.0
    if lo is not None:
        lo = float(lo)
        return lo * 2.0 if lo > 0 else lo + 1.0
    if hi is not None:
        hi = float(hi)
        return hi / 2.0 if hi > 0 else hi - 1.0
    return None


def _gc_number(region_code: str) -> Optional[float]:
    m = re.fullmatch(r"GC(\d+)", region_code or "")
    return float(m.group(1)) if m else None


def _fill(region: Dict[str, Any], gating: Dict[str, float]) -> Dict[str, Any]:
    """Copy a template region, with every parameter given a value."""
    out = json.loads(json.dumps(region))
    for p in out.get("parameters", []):
        value = synthetic_value(p.get("limits"))
        p["value"] = gating.get(p["code"], UNBOUNDED_DEFAULT) if value is None else value
    return out


def _unbounded(region: Dict[str, Any]) -> List[str]:
    return [
        p["code"] for p in region.get("parameters", []) if synthetic_value(p.get("limits")) is None
    ]


# ---------------------------------------------------------------------------
# NSS calls
# ---------------------------------------------------------------------------


def _get_json(url: str, params: Dict[str, Any], timeout: float) -> Any:
    response = _request_with_backoff(requests.get, url, params=params, timeout=timeout)
    response.raise_for_status()
    try:
        return response.json()
    except ValueError as exc:
        raise StreamStatsResponseError(f"{url} returned invalid JSON") from exc


def _estimate(scenario: Dict[str, Any], timeout: float) -> List[Dict[str, Any]]:
    """POST one scenario; return the regression regions NSS answered with."""
    url = f"https://{GENERIC_HOST}/nssservices/Scenarios/Estimate"
    response = _request_with_backoff(requests.post, url, json_body=[scenario], timeout=timeout)
    if response.status_code >= 400:
        raise StreamStatsResponseError(
            f"Estimate rejected ({response.status_code}): {response.text[:300]}"
        )
    try:
        data = response.json()
    except ValueError as exc:
        raise StreamStatsResponseError("Estimate returned invalid JSON") from exc
    if not isinstance(data, list):
        raise StreamStatsResponseError(f"Estimate returned {type(data).__name__}, not a list")
    return [rr for sc in data for rr in sc.get("regressionRegions", []) or []]


def _gating_candidates(region: Dict[str, Any]) -> List[Dict[str, float]]:
    codes = _unbounded(region)
    if not codes:
        return []
    values: List[float] = []
    gc = _gc_number(str(region.get("code", "")))
    if gc is not None:
        values.append(gc)
    values += [v for v in GATING_FALLBACKS if v not in values]
    combos = [dict(zip(codes, combo)) for combo in itertools.product(values, repeat=len(codes))]
    first = {c: UNBOUNDED_DEFAULT for c in codes}
    return [c for c in combos if c != first]


# ---------------------------------------------------------------------------
# Record building
# ---------------------------------------------------------------------------


def _rounding_bound(value: float) -> float:
    """Half a unit in the last of NSS_SIG_FIGS significant figures of ``value``."""
    if value == 0:
        return 0.0
    exponent = math.floor(math.log10(abs(value))) - (NSS_SIG_FIGS - 1)
    return 0.5 * 10.0**exponent


def _rel(a: float, b: float) -> float:
    return abs(a - b) / abs(b) if b else abs(a - b)


def _equation_record(
    stat: Dict[str, Any],
    region: Dict[str, Any],
    inputs: Dict[str, float],
    citation: Citation,
) -> Dict[str, Any]:
    code = str(stat.get("code", ""))
    equation = str(stat.get("equation") or "")
    errors = [{k: e.get(k) for k in ("code", "name", "value")} for e in stat.get("errors") or []]
    asep = next((e["value"] for e in errors if e["code"] == "ASEp"), None)
    try:
        aep: Optional[float] = statistic_code_to_aep(code)
    except ValueError:
        aep = None
    nss_value = stat.get("value")
    rec: Dict[str, Any] = {
        "statistic_code": code,
        "statistic_name": str(stat.get("name", "")),
        "aep": aep,
        "unit": str((stat.get("unit") or {}).get("abbr", "")),
        "equation": equation,
        "nss_value": nss_value,
        "nss_sep": stat.get("sep"),
        "asep_pct": asep,
        "errors": errors,
        "equivalent_years": stat.get("equivalentYears"),
        "interval_bounds": stat.get("intervalBounds"),
        "parsed": None,
        "parse_error": None,
        "round_trip": None,
    }

    rt: Dict[str, Any] = {}
    try:
        direct = evaluate_expression(equation, inputs)
        rt["direct"] = direct
    except (NSSEquationError, KeyError, ZeroDivisionError, OverflowError) as exc:
        direct = None
        rt["direct_error"] = f"{type(exc).__name__}: {exc}"
    if direct is not None and nss_value is not None:
        bound = _rounding_bound(float(nss_value))
        rt["nss_rounding_bound"] = bound
        rt["direct_within_nss_rounding"] = abs(direct - float(nss_value)) <= bound * (1 + 1e-9)

    try:
        parsed = parse_equation(equation)
    except NSSEquationError as exc:
        rec["parse_error"] = str(exc)
        rec["round_trip"] = rt
        return rec
    rec["parsed"] = parsed.to_dict()

    limits = {
        p["code"]: ((p.get("limits") or {}).get("min"), (p.get("limits") or {}).get("max"))
        for p in region.get("parameters", [])
    }
    units = {
        p["code"]: str((p.get("unitType") or {}).get("abbr", ""))
        for p in region.get("parameters", [])
    }
    try:
        if aep is not None:
            eq = to_regression_equation(
                parsed,
                region_code=str(region.get("code", "")),
                region_name=str(region.get("name", "")),
                aep=aep,
                citation=citation,
                limits=limits,
                units=units,
            )
            evaluated = evaluate(eq, inputs).flow_cfs
            rt["method"] = "flowfreq.regression.equations.evaluate"
        else:
            evaluated = 10.0 ** parsed.log10_value(inputs)
            rt["method"] = "ParsedEquation.log10_value (statistic is not an AEP)"
    except (KeyError, OutOfRangeError, ValueError) as exc:
        rt["evaluate_error"] = f"{type(exc).__name__}: {exc}"
        rec["round_trip"] = rt
        return rec
    rt["evaluated"] = evaluated
    if direct is not None:
        rt["rel_vs_direct"] = _rel(evaluated, direct)
    if nss_value is not None:
        rt["rel_vs_nss"] = _rel(evaluated, float(nss_value))
        rt["within_nss_rounding"] = abs(evaluated - float(nss_value)) <= rt.get(
            "nss_rounding_bound", 0.0
        ) * (1 + 1e-9)
    rt["ok"] = bool(
        rt.get("rel_vs_direct", math.inf) <= PARSE_RTOL and rt.get("within_nss_rounding", False)
    )
    rec["round_trip"] = rt
    return rec


def _region_record(
    group: Dict[str, Any],
    template: Dict[str, Any],
    answered: Optional[Dict[str, Any]],
    submitted: Dict[str, Any],
    gating: Dict[str, float],
    tried: List[Dict[str, float]],
    citations: Dict[int, Any],
) -> Dict[str, Any]:
    params = []
    for p_tmpl, p_sub in zip(template.get("parameters", []), submitted.get("parameters", [])):
        lim = p_tmpl.get("limits") or {}
        params.append(
            {
                "code": p_tmpl.get("code"),
                "name": p_tmpl.get("name"),
                "units": (p_tmpl.get("unitType") or {}).get("abbr", ""),
                "min": lim.get("min"),
                "max": lim.get("max"),
                "submitted": p_sub.get("value"),
            }
        )
    citation_id = template.get("citationID")
    cit = citations.get(citation_id) if citation_id is not None else None
    rec: Dict[str, Any] = {
        "statistic_group_id": group["id"],
        "statistic_group_code": group["code"],
        "region_id": template.get("id"),
        "region_code": template.get("code"),
        "region_name": template.get("name"),
        "citation_id": citation_id,
        "status": "resolved",
        "parameters": params,
        "gating_values": gating,
        "equations": [],
    }
    results = (answered or {}).get("results") or []
    if not results:
        rec["status"] = "unresolved"
        why = (
            "Estimate returned the region with no results"
            if answered is not None
            else "Estimate omitted the region from its response"
        )
        if tried:
            why += f"; tried unbounded-parameter values {tried}"
        elif not _unbounded(template):
            why += " (all parameters bounded; submitted in-range midpoints)"
        rec["unresolved_reason"] = why
        return rec
    inputs = {p["code"]: float(p["submitted"]) for p in params}
    citation = Citation(
        publication=(cit.title if cit else f"NSS citation {citation_id}") or "NSS",
        table="not supplied by NSS",
    )
    rec["equations"] = [_equation_record(s, template, inputs, citation) for s in results]
    return rec


def snapshot_state(
    state: str, groups: Sequence[str] = ("PFS",), *, timeout: float = 90.0
) -> Dict[str, Any]:
    """Harvest every regression equation NSS has for one state.

    Parameters
    ----------
    state : str
        NSS region code, e.g. ``"WA"``.
    groups : sequence of str
        Statistic group codes, or ``["all"]`` for every group NSS lists for the state.
    """
    state = state.upper()
    available = list_statistic_groups(state, timeout=timeout)
    wanted = {g.upper() for g in groups}
    chosen = [g for g in available if "ALL" in wanted or str(g["code"]).upper() in wanted]
    missing = wanted - {"ALL"} - {str(g["code"]).upper() for g in available}
    templates_url = f"https://{GENERIC_HOST}/nssservices/regions/{state}/Scenarios"

    staged: List[Tuple[Dict[str, Any], Dict[str, Any], Any, Any, Any, Any]] = []
    for group in chosen:
        scenarios = _get_json(
            templates_url,
            {"statisticgroups": group["id"], "unitsystem": UNIT_SYSTEM},
            timeout,
        )
        for scenario in scenarios:
            regions = scenario.get("regressionRegions", []) or []
            first_gating = {r["id"]: {c: UNBOUNDED_DEFAULT for c in _unbounded(r)} for r in regions}
            filled = [_fill(r, first_gating[r["id"]]) for r in regions]
            body = {k: v for k, v in scenario.items() if k != "links"}
            body["regressionRegions"] = filled
            answered = {rr.get("id"): rr for rr in _estimate(body, timeout)}
            for tmpl, sub in zip(regions, filled):
                ans = answered.get(tmpl["id"])
                gating = first_gating[tmpl["id"]]
                tried: List[Dict[str, float]] = []
                if not (ans or {}).get("results") and gating:
                    tried.append(gating)
                    for cand in _gating_candidates(tmpl):
                        tried.append(cand)
                        retry = _fill(tmpl, cand)
                        one = dict(body, regressionRegions=[retry])
                        got = {rr.get("id"): rr for rr in _estimate(one, timeout)}
                        if (got.get(tmpl["id"]) or {}).get("results"):
                            ans, sub, gating = got[tmpl["id"]], retry, cand
                            break
                staged.append((group, tmpl, ans, sub, gating, tried))

    region_ids = sorted({t["id"] for _, t, *_ in staged})
    citations: Dict[int, Any] = {}
    if region_ids:
        try:
            citations = _fetch_citations(region_ids, timeout=timeout)
        except (StreamStatsResponseError, StreamStatsTransportError, requests.RequestException):
            citations = {}
    regions = [_region_record(g, t, a, s, gat, tr, citations) for g, t, a, s, gat, tr in staged]
    regions.sort(key=lambda r: (str(r["statistic_group_code"]), str(r["region_code"])))

    eqs = [e for r in regions for e in r["equations"]]
    summary = {
        "statistic_groups": len(chosen),
        "regions": len(regions),
        "regions_resolved": sum(r["status"] == "resolved" for r in regions),
        "regions_unresolved": sum(r["status"] == "unresolved" for r in regions),
        "equations": len(eqs),
        "equations_parsed": sum(e["parsed"] is not None for e in eqs),
        "equations_refused": sum(e["parse_error"] is not None for e in eqs),
        "round_trip_ok": sum(bool((e["round_trip"] or {}).get("ok")) for e in eqs),
        "round_trip_failed": sum(
            e["parsed"] is not None and not (e["round_trip"] or {}).get("ok") for e in eqs
        ),
    }
    return {
        "schema": SNAPSHOT_SCHEMA,
        "state": state,
        "retrieved_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "service": {
            "templates": templates_url,
            "estimate": f"https://{GENERIC_HOST}/nssservices/Scenarios/Estimate",
            "unit_system": UNIT_SYSTEM,
        },
        "tool": "tools/snapshot_nss.py",
        "note": SNAPSHOT_NOTE,
        "nss_sep_note": NSS_SEP_NOTE,
        "statistic_groups_requested": sorted(wanted),
        "statistic_groups_missing": sorted(missing),
        "summary": summary,
        "citations": {
            str(k): {
                "title": c.title,
                "author": c.author,
                "citation_url": c.citation_url,
                "last_year_of_data": c.last_year_of_data,
            }
            for k, c in sorted(citations.items())
        },
        "regions": regions,
    }


# ---------------------------------------------------------------------------
# Diff report
# ---------------------------------------------------------------------------


def _limits(region: Dict[str, Any]) -> Dict[str, Tuple[Any, Any, Any]]:
    return {p["code"]: (p.get("min"), p.get("max"), p.get("units")) for p in region["parameters"]}


def diff_snapshots(old: Dict[str, Any], new: Dict[str, Any]) -> List[str]:
    """Human-readable differences between two snapshots of one state.

    Compares regions, their status and parameter limits, and each equation's
    string and ASEp. Returns an empty list when nothing changed.
    """

    def key(r: Dict[str, Any]) -> Tuple[str, str]:
        return (str(r["statistic_group_code"]), str(r["region_code"]))

    a = {key(r): r for r in old.get("regions", [])}
    b = {key(r): r for r in new.get("regions", [])}
    out: List[str] = []
    for k in sorted(a.keys() - b.keys()):
        out.append(f"- region {k[0]}:{k[1]} ({a[k]['region_name']}) removed")
    for k in sorted(b.keys() - a.keys()):
        out.append(
            f"+ region {k[0]}:{k[1]} ({b[k]['region_name']}) added, "
            f"{len(b[k]['equations'])} equations"
        )
    for k in sorted(a.keys() & b.keys()):
        ra, rb, tag = a[k], b[k], f"{k[0]}:{k[1]}"
        if ra["status"] != rb["status"]:
            out.append(f"~ {tag} status {ra['status']} -> {rb['status']}")
        if ra.get("region_name") != rb.get("region_name"):
            out.append(f"~ {tag} name {ra.get('region_name')!r} -> {rb.get('region_name')!r}")
        la, lb = _limits(ra), _limits(rb)
        for code in sorted(la.keys() | lb.keys()):
            if la.get(code) != lb.get(code):
                out.append(f"~ {tag} parameter {code}: {la.get(code)} -> {lb.get(code)}")
        ea = {e["statistic_code"]: e for e in ra["equations"]}
        eb = {e["statistic_code"]: e for e in rb["equations"]}
        for s in sorted(ea.keys() - eb.keys()):
            out.append(f"- {tag} {s} removed")
        for s in sorted(eb.keys() - ea.keys()):
            out.append(f"+ {tag} {s} added: {eb[s]['equation']}")
        for s in sorted(ea.keys() & eb.keys()):
            if ea[s]["equation"] != eb[s]["equation"]:
                out.append(f"~ {tag} {s}: {ea[s]['equation']} -> {eb[s]['equation']}")
            if ea[s].get("asep_pct") != eb[s].get("asep_pct"):
                out.append(f"~ {tag} {s} ASEp {ea[s].get('asep_pct')} -> {eb[s].get('asep_pct')}")
    return out


def previous_snapshot(out_dir: Path, state: str) -> Optional[Path]:
    """The newest existing snapshot for ``state`` in ``out_dir``, if any."""
    files = sorted(out_dir.glob(f"{state.upper()}_????-??-??.json"))
    return files[-1] if files else None


def _report(snap: Dict[str, Any]) -> List[str]:
    s = snap["summary"]
    lines = [
        f"{snap['state']}: {s['regions']} regions ({s['regions_resolved']} resolved, "
        f"{s['regions_unresolved']} unresolved), {s['equations']} equations "
        f"({s['equations_parsed']} parsed, {s['equations_refused']} refused), "
        f"round trip {s['round_trip_ok']} ok / {s['round_trip_failed']} failed"
    ]
    for r in snap["regions"]:
        if r["status"] == "unresolved":
            lines.append(
                f"  unresolved {r['statistic_group_code']}:{r['region_code']} "
                f"{r['region_name']}: {r['unresolved_reason']}"
            )
    return lines


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--state", action="append", required=True, help="NSS region code; repeat")
    ap.add_argument(
        "--group",
        action="append",
        help="statistic group code (default PFS); repeat, or 'all' for every group",
    )
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT, help=f"default {DEFAULT_OUT}")
    ap.add_argument(
        "--concurrency",
        type=int,
        default=1,
        help=f"states fetched in parallel, capped at {MAX_CONCURRENCY}",
    )
    ap.add_argument("--timeout", type=float, default=90.0)
    args = ap.parse_args(argv)
    groups = args.group or ["PFS"]
    args.out.mkdir(parents=True, exist_ok=True)
    today = datetime.now(timezone.utc).date().isoformat()
    workers = max(1, min(args.concurrency, MAX_CONCURRENCY))

    states = [s.upper() for s in args.state]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        snaps = list(pool.map(lambda s: snapshot_state(s, groups, timeout=args.timeout), states))

    for state, snap in zip(states, snaps):
        prev_path = previous_snapshot(args.out, state)
        prev = json.loads(prev_path.read_text(encoding="utf-8")) if prev_path else None
        path = args.out / f"{state}_{today}.json"
        path.write_text(json.dumps(snap, indent=1) + "\n", encoding="utf-8")
        print("\n".join(_report(snap)))
        print(f"  wrote {path}")
        if prev_path is None or prev is None:
            print("  no previous snapshot to diff against")
            continue
        changes = diff_snapshots(prev, snap)
        label = f"{prev_path.name} (retrieved {prev.get('retrieved_utc')})"
        if changes:
            print(f"  {len(changes)} change(s) since {label}:")
            print("\n".join("    " + c for c in changes))
        else:
            print(f"  no equation changes since {label}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
