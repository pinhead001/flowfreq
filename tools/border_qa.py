"""Cross-border consistency of the verified regression equations, on every border.

For each pair of adjacent jurisdictions with verified equations, every catalog gage within
``--border-km`` of the neighbour's outline (20+ peaks, 1-3,000 mi2, not classed regulated or
urban) is delineated once in its own state's StreamStats region. Both states' equations are then
evaluated on that one basin, at AEP 50, 10 and 1 percent, and compared with each other and with
a Bulletin 17C fit at the gage (station skew). Results go to ``docs/BORDER_QA.md`` and
``docs/border_qa/rows.csv``. Nothing is changed: differences are documented, not smoothed.

How each state's estimate is made:

- **Home state:** NSS ``bylocation`` on the delineated basin, then
  :func:`flowfreq.ungaged.assemble`, so every library rule applies (limits, preferences, gates,
  Indiana's region at the point).
- **Neighbour:** the basin rarely lies in the neighbour's StreamStats grid, so the neighbour's
  equations are evaluated on the home basin's characteristics, as ``tools/wave_qa.py`` does.
  Its regions come from NSS ``bylocation`` on the basin where the basin reaches the neighbour,
  rescaled to cover the basin; otherwise from a circle around the outlet, grown until NSS
  answers, taking the region with the largest share. A characteristic the home region does not
  compute may be filled from an equivalent one (:data:`EXTRA_PROXIES`, beside
  :data:`flowfreq.validation.wave_qa.PROXIES`); each substitution is recorded on the row.
- Both sides allow extrapolation and record it. An overlap the library leaves to the caller
  takes the first set not named drainage-area-only, and records that.

Stages, cached under ``--cache`` (default ``.cache/border_qa``, gitignored), so a run resumes::

    python tools/border_qa.py candidates
    python tools/border_qa.py run --workers 3
    python tools/border_qa.py report

Network: StreamStats, NSS, TIGERweb and the Water Data API (set ``USGS_API_KEY``).
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys

from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import pandas as pd
import requests

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from flowfreq import streamstats as ss  # noqa: E402
from flowfreq.catalog import load_catalog  # noqa: E402
from flowfreq.regression.library import load_state  # noqa: E402
from flowfreq.ungaged import OverlappingRegionsError, assemble  # noqa: E402
from flowfreq.validation import wave_qa as qa  # noqa: E402

logger = logging.getLogger("border_qa")

AEPS: Tuple[float, ...] = (0.5, 0.1, 0.01)
MIN_PEAKS = 20
AREA_RANGE = (1.0, 3000.0)
#: A delineation whose area differs from the catalog's by more than this factor is not the
#: gage's basin (a mis-snapped point) and is not compared.
DA_TOL = 1.25
#: Median log10 ratio, over 3 or more gages, above which a pair's offset is called systematic.
SYSTEMATIC_LOG = 0.1
TIGER = (
    "https://tigerweb.geo.census.gov/arcgis/rest/services/TIGERweb/State_County/MapServer/0/query"
)

#: Adjacent jurisdictions, both with verified equations; each pair once.
ADJACENT: Tuple[Tuple[str, str], ...] = tuple(
    (a, b)
    for entry in (
        "AL:FL,GA,MS,TN AR:LA,MO,MS,OK,TN,TX AZ:CA,CO,NM,NV,UT CA:NV,OR CO:KS,NE,NM,OK,UT,WY "
        "CT:MA,NY,RI DE:MD,NJ,PA FL:GA GA:NC,SC,TN IA:IL,MN,MO,NE,SD,WI ID:MT,NV,OR,UT,WA,WY "
        "IL:IN,KY,MO,WI IN:KY,MI,OH KS:MO,NE,OK KY:MO,OH,TN,VA,WV LA:MS,TX MA:NH,NY,RI,VT "
        "MD:PA,VA,WV ME:NH MI:OH,WI MN:ND,SD,WI MO:NE,OK,TN MS:TN MT:ND,SD,WY NC:SC,TN,VA ND:SD "
        "NE:SD,WY NH:VT NJ:NY,PA NM:OK,TX,UT NV:OR,UT NY:PA,VT OH:PA,WV OK:TX OR:WA PA:WV SD:WY "
        "TN:VA UT:WY VA:WV"
    ).split()
    for a, bs in [entry.split(":")]
    for b in bs.split(",")
)

#: Equivalent characteristics beyond ``wave_qa.PROXIES``, in order of preference: the same
#: quantity by another method or data year. Slopes: the 10-85 channel slope in ft/mi along the
#: main channel or the longest flow path. Basin shape: by main channel or longest flow path.
#: Impervious area: NLCD years.
EXTRA_PROXIES: Dict[str, Tuple[str, ...]] = {
    "CSL10_85": ("CSL10_85fm", "CSL1085LFP"),
    "CSL10_85fm": ("CSL10_85", "CSL1085LFP"),
    "CSL1085LFP": ("CSL10_85", "CSL10_85fm"),
    "BSHAPE": ("BSHAPELFP",),
    "BSHAPELFP": ("BSHAPE",),
    "LC11IMP": ("LC16IMP", "LC19IMP", "LC06IMP"),
    "LC06IMP": ("LC11IMP", "LC01IMP", "LC16IMP"),
    "LC19IMP": ("LC21IMP", "LC16IMP", "LC11IMP"),
    "IMPNLCD01": ("LC01IMP", "LC06IMP", "LC11IMP"),
}


def neighbors(state: str) -> List[str]:
    """The states adjacent to ``state`` in :data:`ADJACENT`."""
    return sorted({b for a, b in ADJACENT if a == state} | {a for a, b in ADJACENT if b == state})


def fill_proxies(
    characteristics: Mapping[str, float], lat: float, lon: float
) -> Tuple[Dict[str, float], List[str]]:
    """``wave_qa.prepare_characteristics`` plus :data:`EXTRA_PROXIES`; substitutions listed."""
    values, subs = qa.prepare_characteristics(characteristics, lat, lon)
    have = dict(values)
    extra: List[str] = []
    for target, sources in EXTRA_PROXIES.items():
        if target in have:
            continue
        src = next((s for s in sources if s in have), None)
        if src is not None:
            values[target] = have[src]
            extra.append(f"{target}={src}")
    return values, list(subs) + extra


# ----------------------------------------------------------------------------
# candidates
# ----------------------------------------------------------------------------


def _read(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _write(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, default=str), encoding="utf-8")
    tmp.replace(path)


def state_outlines(cache: Path) -> Dict[str, Dict[str, Any]]:
    """Census TIGERweb outlines (about 200 m generalised), cached per state."""
    out = {}
    for st in sorted({s for pair in ADJACENT for s in pair}):
        path = cache / "states" / f"{st}.json"
        if not path.exists():
            r = requests.get(
                TIGER,
                params={
                    "where": f"STUSAB='{st}'",
                    "outFields": "STUSAB",
                    "f": "geojson",
                    "outSR": 4326,
                    "maxAllowableOffset": 0.002,
                    "returnGeometry": "true",
                },
                timeout=120,
            )
            r.raise_for_status()
            feats = r.json()["features"]
            if len(feats) != 1:
                raise RuntimeError(f"TIGERweb returned {len(feats)} features for {st}")
            _write(path, feats[0]["geometry"])
        out[st] = _read(path)
    return out


def candidates(cache: Path, border_km: float) -> pd.DataFrame:
    """One row per (gage, neighbour) within ``border_km`` of the neighbour's outline."""
    outlines = state_outlines(cache)
    cat = load_catalog()
    cat = cat[pd.to_numeric(cat["n_peaks"], errors="coerce") >= MIN_PEAKS]
    cat = cat[~cat["regulation_class"].isin(["regulated", "urban"])]
    cat = cat[cat["drainage_area_sqmi"].between(*AREA_RANGE)]
    rows = []
    for home in sorted(outlines):
        for _, g in cat[cat["state"] == home].iterrows():
            for other in neighbors(home):
                d = qa.distance_km(float(g["latitude"]), float(g["longitude"]), outlines[other])
                if d <= border_km:
                    rows.append(
                        {
                            "site_no": g["site_no"],
                            "site_name": g["site_name"],
                            "home": home,
                            "neighbor": other,
                            "distance_km": round(d, 2),
                            "latitude": float(g["latitude"]),
                            "longitude": float(g["longitude"]),
                            "drainage_area_sqmi": float(g["drainage_area_sqmi"]),
                            "n_peaks": int(g["n_peaks"]),
                        }
                    )
    df = pd.DataFrame(rows).sort_values(["home", "neighbor", "distance_km"])
    df.to_csv(cache / "candidates.csv", index=False)
    logger.info("candidates: %d rows, %d gages", len(df), df["site_no"].nunique())
    return df


# ----------------------------------------------------------------------------
# run
# ----------------------------------------------------------------------------


def _at_site(site_no: str) -> Dict[str, Any]:
    from flowfreq.usgs import USGSgage
    from flowfreq.workflow import run_ffa

    p = USGSgage(site_no).download_peak_flow().dropna(subset=["peak_flow_cfs"])
    res = run_ffa(
        p["peak_flow_cfs"].to_numpy(float),
        p["water_year"].to_numpy(int),
        station_skew_only=True,
        peak_codes=p["qualification_code"].tolist(),
        site_no=site_no,
        allow_regulated=False,
    )
    if res.get("error"):
        return {"error": str(res["error"])[:200]}
    q = res["quantile_df"]
    out: Dict[str, Any] = {"n": int(len(p))}
    for aep in AEPS:
        row = q.iloc[(q["AEP (%)"] - aep).abs().argsort()[:1]].iloc[0]
        out[f"{aep:g}"] = [
            float(row["Flow (cfs)"]),
            float(row["Lower 90% CI"]),
            float(row["Upper 90% CI"]),
        ]
    return out


def default_pick(candidates: Sequence[Sequence[Any]]) -> str:
    """The set to use where the library leaves an overlap to the caller: a rural full set.

    An urban set on a rural basin, or a drainage-area-only set where a full set applies, would
    make a border difference that is the tool's choice, not the equations'. Candidates are the
    library's ``(nss_code, region_name, year)``.
    """

    def rank(c: Sequence[Any]) -> Tuple[int, int]:
        name = str(c[1]).lower()
        return ("urban" in name, "drainage area only" in name)

    return str(sorted(candidates, key=rank)[0][0])


def _assemble(lib: Any, regions: Sequence[Any], values: Mapping[str, float]) -> Dict[str, Any]:
    pick = None
    try:
        a = assemble(
            lib, regions, values, aeps=None, allow_extrapolation=True, weighting_space="log"
        )
    except OverlappingRegionsError as exc:
        pick = default_pick(exc.candidates)
        a = assemble(
            lib,
            regions,
            values,
            aeps=None,
            allow_extrapolation=True,
            weighting_space="log",
            region=pick,
        )
    out: Dict[str, Any] = {
        "regions": [f"{r.code}:{a.region_codes.get(r.code, '?')}" for r in a.used],
        "picked": pick,
        "choice": a.region_choice_source,
    }
    # Every AEP the regions publish (Indiana publishes no 50 percent), then those compared.
    q = a.quantiles.set_index("aep")
    for aep in AEPS:
        hit = [x for x in q.index if abs(x - aep) < 1e-9]
        if hit:
            row = q.loc[hit[0]]
            out[f"{aep:g}"] = [float(row["flow_cfs"]), bool(row["extrapolated"])]
    return out


def _neighbor_regions(
    state: str, ws: ss.WatershedCharacteristics, distance_km: float
) -> Tuple[List[Any], str]:
    lib = load_state(state, warn=False)
    known = set(qa.nss_region_map(lib))
    located = [
        r
        for r in ss.locate_regression_regions(state, ws.polygon_geojson)
        if r.code in known and r.percent_weight > 0
    ]
    if located:
        total = sum(r.percent_weight for r in located)
        # Partitions (sum near the basin's share in the state) are rescaled to cover it;
        # overlapping sets keep their own weights for the library to resolve.
        if total <= 100.0 + 1.5:
            located = [replace(r, percent_weight=100.0 * r.percent_weight / total) for r in located]
        return located, "basin"
    lat, lon = ws.snap.snapped_lat, ws.snap.snapped_lon
    radius = max(3.0, 1.3 * distance_km)
    for _ in range(6):
        found = [
            r
            for r in ss.locate_regression_regions(state, qa.circle(lat, lon, radius))
            if r.code in known and r.percent_weight > 0
        ]
        if found:
            top = max(r.percent_weight for r in found)
            return [
                replace(r, percent_weight=100.0) for r in found if r.percent_weight >= top - 1
            ], f"nearest ({radius:.0f} km)"
        radius *= 1.6
    return [], "no region located"


def run_site(cache: Path, site_no: str, rows: pd.DataFrame) -> None:
    """Delineate one gage, fit it, and evaluate its home and every neighbour; cached."""
    path = cache / "sites" / f"{site_no}.json"
    if path.exists() and "TransportError" not in path.read_text(encoding="utf-8"):
        return  # a record with a service failure in it is redone
    g = rows.iloc[0]
    rec: Dict[str, Any] = {
        "site_no": site_no,
        "home": g["home"],
        "drainage_area_sqmi": float(g["drainage_area_sqmi"]),
    }
    try:
        rec["at_site"] = _at_site(site_no)
    except Exception as exc:  # noqa: BLE001 - a QA row records every failure
        rec["at_site"] = {"error": f"{type(exc).__name__}: {str(exc)[:160]}"}
    try:
        ws = ss.delineate_and_get_characteristics(
            g["home"], float(g["latitude"]), float(g["longitude"])
        )
    except Exception as exc:  # noqa: BLE001
        rec["error"] = f"{type(exc).__name__}: {str(exc)[:160]}"
        if "Transport" not in rec["error"]:
            _write(path, rec)
        return
    raw = {k: c.value for k, c in ws.characteristics.items() if isinstance(c.value, (int, float))}
    values, subs = fill_proxies(raw, ws.snap.snapped_lat, ws.snap.snapped_lon)
    rec["characteristics"] = raw  # so a later version can re-evaluate without StreamStats
    area = values.get("DRNAREA") or values.get("CONTDA")
    rec["delineated_sqmi"] = area
    rec["mis_snapped"] = not (area and 1 / DA_TOL <= area / rec["drainage_area_sqmi"] <= DA_TOL)
    try:
        home_regions = ss.locate_regression_regions(g["home"], ws.polygon_geojson)
        rec["home_est"] = _assemble(load_state(g["home"], warn=False), home_regions, values)
    except Exception as exc:  # noqa: BLE001
        rec["home_est"] = {"error": f"{type(exc).__name__}: {str(exc)[:160]}"}
    rec["neighbors"] = {}
    for _, row in rows.iterrows():
        nb = row["neighbor"]
        entry: Dict[str, Any] = {"distance_km": float(row["distance_km"])}
        try:
            regions, how = _neighbor_regions(nb, ws, float(row["distance_km"]))
            entry["how"] = how
            if regions:
                used = _assemble(load_state(nb, warn=False), regions, values)
                entry.update(used)
                entry["substitutions"] = [
                    s for s in subs if s.split("=")[0] in _codes_needed(nb, used)
                ]
        except Exception as exc:  # noqa: BLE001
            entry["error"] = f"{type(exc).__name__}: {str(exc)[:160]}"
        rec["neighbors"][nb] = entry
    _write(path, rec)
    logger.info(
        "%s %s: home %s; %s",
        site_no,
        g["home"],
        rec["home_est"].get("0.01", rec["home_est"].get("error")),
        {k: v.get("0.01", v.get("error", v.get("how"))) for k, v in rec["neighbors"].items()},
    )


def _codes_needed(state: str, used: Mapping[str, Any]) -> set:
    lib = load_state(state, warn=False)
    codes = {r.split(":", 1)[1] for r in used.get("regions", [])}
    return {v.code for e in lib.equations if e.region_code in codes for v in e.variables}


def run(cache: Path, workers: int, limit: Optional[int]) -> None:
    cand = pd.read_csv(cache / "candidates.csv", dtype={"site_no": str})
    groups = list(cand.groupby("site_no", sort=False))
    if limit:
        groups = groups[:limit]

    def one(item: Tuple[str, pd.DataFrame]) -> None:
        try:
            run_site(cache, item[0], item[1])
        except Exception:  # noqa: BLE001 - keep going; the site is retried on resume
            logger.exception("site %s failed", item[0])

    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(one, groups))


# ----------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------


def rows_from_cache(cache: Path) -> pd.DataFrame:
    """One row per (gage, neighbour, AEP) with both states' flows and the at-site fit."""
    out = []
    for path in sorted((cache / "sites").glob("*.json")):
        rec = _read(path)
        for nb, entry in (rec.get("neighbors") or {}).items():
            for aep in AEPS:
                key = f"{aep:g}"
                home = (rec.get("home_est") or {}).get(key)
                other = entry.get(key)
                site = (rec.get("at_site") or {}).get(key)
                out.append(
                    {
                        "site_no": rec["site_no"],
                        "home": rec["home"],
                        "neighbor": nb,
                        "aep": aep,
                        "distance_km": entry.get("distance_km"),
                        "catalog_sqmi": rec.get("drainage_area_sqmi"),
                        "delineated_sqmi": rec.get("delineated_sqmi"),
                        "mis_snapped": rec.get("mis_snapped", True),
                        "home_regions": " ".join((rec.get("home_est") or {}).get("regions", [])),
                        "home_q": home[0] if home else None,
                        "home_extrapolated": home[1] if home else None,
                        "neighbor_regions": " ".join(entry.get("regions", [])),
                        "neighbor_how": entry.get("how"),
                        "neighbor_q": other[0] if other else None,
                        "neighbor_extrapolated": other[1] if other else None,
                        "substitutions": " ".join(entry.get("substitutions", [])),
                        "picked": entry.get("picked"),
                        "at_site_q": site[0] if site else None,
                        "at_site_lo": site[1] if site else None,
                        "at_site_hi": site[2] if site else None,
                        "why_not": _why_not(rec, entry, key),
                        "clean": _clean(rec, entry, key),
                    }
                )
    return pd.DataFrame(out)


def _why_not(rec: Mapping[str, Any], entry: Mapping[str, Any], key: str) -> str:
    if rec.get("error"):
        return "delineation: " + str(rec["error"]).split(":")[0]
    if rec.get("mis_snapped"):
        return "delineated area differs from the catalog's"
    if key not in (rec.get("home_est") or {}):
        return "home: " + str((rec.get("home_est") or {}).get("error", "no estimate")).split(":")[0]
    if key not in entry:
        if entry.get("how") == "no region located":
            return "neighbour: no region located"
        err = str(entry.get("error", "no estimate"))
        if "needs a basin characteristic" in err:
            return "neighbour: characteristic not computed in the home region"
        return "neighbour: " + err.split(":")[0]
    return ""


def _clean(rec: Mapping[str, Any], entry: Mapping[str, Any], key: str) -> bool:
    """A comparison neither side reached by extrapolating or by a default overlap pick.

    Extrapolation here is usually a set used far outside its range because the set that
    applies needs a characteristic the home region doesn't compute (Georgia's 1 mi2 sets on
    Alabama basins); an urban default pick puts an urban set on a rural basin. Either makes a
    difference that is the check's, not the equations'.
    """
    home = (rec.get("home_est") or {}).get(key)
    other = entry.get(key)
    if not home or not other:
        return False
    urban_pick = bool(entry.get("picked")) and "urban" in " ".join(entry.get("regions", [])).lower()
    return not (home[1] or other[1] or urban_pick)


def summarize(rows: pd.DataFrame) -> Dict[str, Any]:
    """Counts and factors by AEP and, at 1 percent, by pair; from :func:`rows_from_cache`."""
    ok = rows[rows["why_not"] == ""].copy()
    ok["log_ratio"] = (ok["neighbor_q"] / ok["home_q"]).map(math.log10)
    by_aep = []
    for aep, g in ok.groupby("aep", sort=False):
        a = g["log_ratio"].abs()
        by_aep.append(
            {
                "aep": aep,
                "n": len(g),
                "median_factor": round(10 ** a.median(), 2),
                "within_1.5x": int((a < math.log10(1.5)).sum()),
                "1.5_to_2x": int(((a >= math.log10(1.5)) & (a < math.log10(2))).sum()),
                "2_to_3x": int(((a >= math.log10(2)) & (a < math.log10(3))).sum()),
                "3x_or_more": int((a >= math.log10(3)).sum()),
            }
        )
    one = ok[ok["aep"] == 0.01].copy()
    clean_col = one["clean"] if "clean" in one else pd.Series(True, index=one.index)
    one["clean"] = clean_col.fillna(False).astype(bool)
    pairs = []
    for (h, n), g in one.groupby(["home", "neighbor"]):
        med = float(g["log_ratio"].median())
        c = g[g["clean"]]
        cmed = float(c["log_ratio"].median()) if len(c) else float("nan")
        pairs.append(
            {
                "home": h,
                "neighbor": n,
                "n": len(g),
                "median_neighbor_over_home": round(10**med, 2),
                "systematic": len(g) >= 3 and abs(med) > SYSTEMATIC_LOG,
                "clean_n": len(c),
                "clean_median": round(10**cmed, 2) if len(c) else None,
                "clean_systematic": len(c) >= 3 and abs(cmed) > SYSTEMATIC_LOG,
            }
        )
    fit = one.dropna(subset=["at_site_q"])
    home_err = (fit["home_q"] / fit["at_site_q"]).map(math.log10).abs()
    nb_err = (fit["neighbor_q"] / fit["at_site_q"]).map(math.log10).abs()
    return {
        "by_aep": by_aep,
        "pairs": pairs,
        "home_nearer": int((home_err < nb_err).sum()),
        "with_fit": len(fit),
        "why_not": rows[(rows["aep"] == 0.01) & (rows["why_not"] != "")]["why_not"]
        .value_counts()
        .to_dict(),
        "gages": int(rows["site_no"].nunique()),
        "not_delineated": {},
        "compared_gages": int(ok["site_no"].nunique()),
        "clean_at_1pct": int(one["clean"].sum()),
    }


def report(cache: Path, out_md: Path, out_csv: Path) -> None:
    rows = rows_from_cache(cache)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(out_csv, index=False, float_format="%.6g")
    s = summarize(rows)
    failed: Dict[str, int] = {}
    for path in (cache / "sites").glob("*.json"):
        err = _read(path).get("error")
        if err:
            kind = str(err).split(":")[0]
            failed[kind] = failed.get(kind, 0) + 1
    s["not_delineated"] = failed
    total = len(pd.read_csv(cache / "candidates.csv", dtype={"site_no": str})["site_no"].unique())
    lines = [
        "# Cross-border consistency, every border",
        "",
        "Generated by `tools/border_qa.py` from its cache; rows in `docs/border_qa/rows.csv`. Method in the tool's "
        "docstring. Differences are documented, not smoothed. Values are the neighbour's flow over the home "
        "state's, on the same delineated basin.",
        "",
        f"Candidate gages: {total}. Delineated: {s['gages']}; compared at one AEP or more: "
        f"{s['compared_gages']}. Not delineated: "
        + (", ".join(f"{k} {v}" for k, v in sorted(failed.items(), key=lambda x: -x[1])) or "none")
        + f"; not yet run: {total - s['gages'] - sum(failed.values())}.",
        "",
        "| AEP (%) | Compared | Median factor | Within 1.5x | 1.5-2x | 2-3x | 3x or more |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for b in s["by_aep"]:
        lines.append(
            f"| {b['aep'] * 100:g} | {b['n']} | {b['median_factor']} | {b['within_1.5x']} | "
            f"{b['1.5_to_2x']} | {b['2_to_3x']} | {b['3x_or_more']} |"
        )
    lines += [
        "",
        f"At 1 percent, the home state's equation was nearer the at-site fit at {s['home_nearer']} of "
        f"{s['with_fit']} comparisons.",
        "",
        "## Pairs at 1 percent",
        "",
        "Systematic: 3 or more gages and a median factor beyond 1.26 (0.1 log10). Clean: neither side "
        "extrapolated and no urban set picked by default; "
        f"{s['clean_at_1pct']} of the comparisons at 1 percent are clean.",
        "",
        "| Home | Neighbour | Gages | Median neighbour/home | Systematic | Clean gages | Clean median "
        "| Clean systematic |",
        "|---|---|---:|---:|---|---:|---:|---|",
    ]
    for p in sorted(s["pairs"], key=lambda p: (p["home"], p["neighbor"])):
        lines.append(
            f"| {p['home']} | {p['neighbor']} | {p['n']} | {p['median_neighbor_over_home']} | "
            f"{'yes' if p['systematic'] else ''} | {p['clean_n']} | "
            f"{'' if p['clean_median'] is None else p['clean_median']} | "
            f"{'yes' if p['clean_systematic'] else ''} |"
        )
    lines += ["", "## Not compared, at 1 percent", "", "| Reason | Rows |", "|---|---:|"]
    for why, n in sorted(s["why_not"].items(), key=lambda x: -x[1]):
        lines.append(f"| {why} | {n} |")
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    logger.info("wrote %s and %s", out_md, out_csv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", choices=["candidates", "run", "report"])
    ap.add_argument("--cache", type=Path, default=REPO_ROOT / ".cache" / "border_qa")
    ap.add_argument("--border-km", type=float, default=20.0)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--limit", type=int, default=None, help="run at most this many gages")
    ap.add_argument("--out", type=Path, default=REPO_ROOT / "docs" / "BORDER_QA.md")
    ap.add_argument("--csv", type=Path, default=REPO_ROOT / "docs" / "border_qa" / "rows.csv")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
    for noisy in ("flowfreq", "urllib3"):
        logging.getLogger(noisy).setLevel(logging.ERROR)
    args.cache.mkdir(parents=True, exist_ok=True)
    if args.stage == "candidates":
        candidates(args.cache, args.border_km)
    elif args.stage == "run":
        run(args.cache, args.workers, args.limit)
    else:
        report(args.cache, args.out, args.csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
