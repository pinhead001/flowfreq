"""Wave 1 and 2 definition-of-done QA (``docs/MASTER_ROADMAP.md`` §3.2).

Runs the cross-border consistency check and the transposition LOOCV for the
Columbia (WA, OR, ID, MT) and Colorado (CO, UT, WY, NM, AZ, NV) waves and writes
the tables that ``docs/WAVE_QA_W1_W2.md`` reports. The computations live in
:mod:`flowfreq.validation.wave_qa`; this script fetches and caches the inputs.

Stages, each cached under ``--cache`` (default ``.cache/wave_qa``, gitignored),
so an interrupted run resumes and a finished one re-runs offline:

``states``      State outlines, Census TIGERweb (``State_County`` layer 0), ~200 m
                generalised.
``candidates``  Catalog gages (``flowfreq/data/gage_catalog.csv.gz``) of the ten
                states with >= 20 years of peaks, not classed regulated or urban,
                that lie within ``--border-km`` of a neighbouring wave state
                (``border``), or farther but close enough for their basin to reach
                the line (``straddle?``: distance < 1.2 sqrt(drainage area)).
``peaks``       Annual peaks, one bulk Water Data API ``peaks`` query per state
                (``tools/build_gage_catalog.py``'s cached fetch, paced).
``b17c``        Bulletin 17C (EMA, :func:`flowfreq.workflow.run_ffa`) at every
                eligible gage of the ten states, with peak codes applied. Regional
                skew where a verified study covers the gage's HUC8 (the PNW study for
                WA/OR/ID, advisory in western MT), else station skew -- no Wave 2
                state has a verified regional skew yet.
``delineate``   StreamStats delineation and basin characteristics for each
                candidate, in its own state's StreamStats region.
``locate``      Regression regions of each basin: NSS ``bylocation`` for the home
                state and the neighbour; for a basin wholly outside the neighbour's
                regions, the neighbour region nearest the outlet (a circle grown
                until NSS answers). Eastern Oregon's regions, not in NSS, come from
                the StreamStats region layer 31 polygons. The basin's share in the
                neighbour state is a 400-point grid sample against its outline.
``evaluate``    Both states' equations on each basin, by AEP, against each other
                and against the gage's B17C quantile (threshold: see
                :mod:`flowfreq.validation.wave_qa`).
``dar``         Leave-one-out drainage-area-ratio transposition over the B17C
                network; regions of targets outside Wave 1's catalog column come
                from NSS ``bylocation`` on a 0.5 km circle.
``report``      Markdown tables and CSVs under ``--out`` (default ``docs/wave_qa``).

Usage::

    python tools/wave_qa.py all
    python tools/wave_qa.py evaluate report     # re-run offline from the cache

Network: Water Data API through :func:`flowfreq.waterdata.request` (1000 requests
an hour per IP without ``USGS_API_KEY``; the bulk peaks query is about 15 requests
in all); StreamStats and NSS through :mod:`flowfreq.streamstats`, serially.
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import requests

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tools"))

from flowfreq import streamstats as ss  # noqa: E402
from flowfreq.catalog import load_catalog, peaks_by_site  # noqa: E402
from flowfreq.regression.library import load_state  # noqa: E402
from flowfreq.validation import wave_qa as qa  # noqa: E402

logger = logging.getLogger("wave_qa")

STATES: Tuple[str, ...] = qa.WAVE1 + qa.WAVE2
AEPS: Tuple[float, ...] = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002)
REPORT_AEPS: Tuple[float, ...] = (0.5, 0.1, 0.01)
MIN_YEARS = 20
#: A delineation whose DRNAREA differs from the catalog's drainage area by more than
#: this factor is not the gage's basin (a mis-snapped point) and is dropped.
DA_TOL = 1.2
#: Median log10 difference above which a consistent one-sided border offset is
#: reported as systematic (0.1 log10 = 26 %, about half of a typical SEP).
SYSTEMATIC_LOG = 0.1
TIGER = (
    "https://tigerweb.geo.census.gov/arcgis/rest/services/TIGERweb/State_County/MapServer/0/query"
)


# ----------------------------------------------------------------------------
# cache helpers
# ----------------------------------------------------------------------------


def _read(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _write(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data), encoding="utf-8")
    tmp.replace(path)


def _cached(path: Path, fetch: Any) -> Any:
    if path.exists():
        return _read(path)
    data = fetch()
    _write(path, data)
    return data


# ----------------------------------------------------------------------------
# states
# ----------------------------------------------------------------------------


def state_outlines(cache: Path) -> Dict[str, Dict[str, Any]]:
    def fetch(st: str) -> Dict[str, Any]:
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
        return dict(feats[0]["geometry"])

    return {st: _cached(cache / "states" / f"{st}.json", lambda st=st: fetch(st)) for st in STATES}


# ----------------------------------------------------------------------------
# candidates
# ----------------------------------------------------------------------------


def eligible_catalog() -> pd.DataFrame:
    cat = load_catalog()
    cat = cat[cat["state"].isin(STATES)]
    cat = cat[pd.to_numeric(cat["n_peaks"], errors="coerce") >= MIN_YEARS]
    cat = cat[~cat["regulation_class"].isin(["regulated", "urban"])]
    cat = cat.dropna(subset=["latitude", "longitude"])
    return cat.reset_index(drop=True)


def candidates(cache: Path, border_km: float) -> pd.DataFrame:
    path = cache / "candidates.csv"
    outlines = state_outlines(cache)
    cat = eligible_catalog()
    rows = []
    for a, b in qa.BORDER_PAIRS:
        for home, other in ((a, b), (b, a)):
            sub = cat[cat["state"] == home]
            for _, g in sub.iterrows():
                d = qa.distance_km(g["latitude"], g["longitude"], outlines[other])
                da = float(g["drainage_area_sqmi"]) if pd.notna(g["drainage_area_sqmi"]) else 0.0
                reach = 1.2 * math.sqrt(da * 2.59) if da > 0 else 0.0
                if d <= border_km:
                    kind = "border"
                elif d <= min(reach, 150.0):
                    kind = "straddle?"
                else:
                    continue
                rows.append(
                    {
                        "site_no": g["site_no"],
                        "site_name": g["site_name"],
                        "home": home,
                        "neighbor": other,
                        "kind": kind,
                        "distance_km": round(d, 2),
                        "latitude": g["latitude"],
                        "longitude": g["longitude"],
                        "drainage_area_sqmi": g["drainage_area_sqmi"],
                        "huc8": g["huc8"],
                        "n_peaks": g["n_peaks"],
                        "regulation_class": g["regulation_class"],
                    }
                )
    df = pd.DataFrame(rows).sort_values(["home", "neighbor", "distance_km"])
    df.to_csv(path, index=False)
    logger.info(
        "candidates: %d (site, neighbour) rows, %d sites; %s",
        len(df),
        df["site_no"].nunique(),
        df.groupby("kind").size().to_dict(),
    )
    return df


# ----------------------------------------------------------------------------
# peaks and B17C
# ----------------------------------------------------------------------------


def state_peaks(cache: Path, pace: float) -> Dict[str, pd.DataFrame]:
    import build_gage_catalog as bgc  # tools/, on sys.path above

    api = bgc.Paced(pace)
    out: Dict[str, pd.DataFrame] = {}
    for st in STATES:
        fips = bgc.POSTAL_FIPS[st]
        feats = bgc.state_peak_features(api, fips, cache / "peaks", False)
        out.update(peaks_by_site(feats))
        logger.info("peaks %s: %d features (%d API calls so far)", st, len(feats), api.calls)
    return out


def _skew_for(state: str, huc8: Optional[str], lat: float, lon: float) -> Tuple[Any, Any, str]:
    from flowfreq.regional_skew import regional_skew_at_huc

    if state in qa.WAVE1 and huc8 and str(huc8) != "nan":
        try:
            rs = regional_skew_at_huc(str(huc8), state, point=(lat, lon), timeout=60)
            return rs.skew, math.sqrt(rs.skew_mse), f"{rs.skew_region} ({rs.citation[:40]})"
        except Exception as exc:  # unavailable, ambiguous HU, service down
            return None, None, f"station (regional unavailable: {str(exc)[:60]})"
    return None, None, "station (no verified regional skew)"


def _fit_one(
    args: Tuple[str, str, Optional[str], float, float, List[Dict[str, Any]]]
) -> Dict[str, Any]:
    site, state, huc8, lat, lon, records = args
    from flowfreq.workflow import run_ffa

    frame = pd.DataFrame(records)
    frame = frame[pd.to_numeric(frame["peak_flow_cfs"], errors="coerce").notna()]
    frame = frame.sort_values("water_year")
    out: Dict[str, Any] = {"site_no": site, "state": state, "n": int(len(frame))}
    skew, skew_se, source = _skew_for(state, huc8, lat, lon)
    out["skew_source"] = source
    try:
        res = run_ffa(
            frame["peak_flow_cfs"].to_numpy(float),
            frame["water_year"].to_numpy(int),
            regional_skew=skew,
            regional_skew_se=skew_se,
            station_skew_only=skew is None,
            peak_codes=frame["qualification_code"].tolist(),
        )
    except Exception as exc:
        out["error"] = str(exc)[:200]
        return out
    if res["error"]:
        out["error"] = str(res["error"])[:200]
        return out
    q = res["quantile_df"]
    out["skew_used"] = res["parameters"]["skew_used"]
    for aep in AEPS:
        k = int(np.argmin(np.abs(q["AEP (%)"].to_numpy(float) - aep)))
        if abs(float(q["AEP (%)"].iloc[k]) - aep) > 1e-6:
            continue
        out[f"q_{aep:g}"] = float(q["Flow (cfs)"].iloc[k])
        out[f"lo_{aep:g}"] = float(q["Lower 90% CI"].iloc[k])
        out[f"hi_{aep:g}"] = float(q["Upper 90% CI"].iloc[k])
    return out


def b17c(cache: Path, pace: float, workers: int) -> pd.DataFrame:
    path = cache / "b17c.csv"
    done = pd.read_csv(path, dtype={"site_no": str}) if path.exists() else pd.DataFrame()
    have = set(done["site_no"]) if len(done) else set()
    cat = eligible_catalog()
    peaks = state_peaks(cache, pace)
    jobs = []
    for _, g in cat.iterrows():
        if g["site_no"] in have or g["site_no"] not in peaks:
            continue
        jobs.append(
            (
                g["site_no"],
                g["state"],
                g["huc8"] if pd.notna(g["huc8"]) else None,
                float(g["latitude"]),
                float(g["longitude"]),
                peaks[g["site_no"]].to_dict("records"),
            )
        )
    logger.info("b17c: %d to fit, %d cached", len(jobs), len(have))
    rows: List[Dict[str, Any]] = []
    t0 = time.monotonic()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_fit_one, j) for j in jobs]
        for i, f in enumerate(as_completed(futs), 1):
            rows.append(f.result())
            if i % 100 == 0 or i == len(futs):
                logger.info("b17c: %d/%d (%.0f s)", i, len(futs), time.monotonic() - t0)
                pd.concat([done, pd.DataFrame(rows)], ignore_index=True).to_csv(path, index=False)
    out = pd.concat([done, pd.DataFrame(rows)], ignore_index=True)
    out.to_csv(path, index=False)
    return out


# ----------------------------------------------------------------------------
# StreamStats
# ----------------------------------------------------------------------------


def _delineate_one(cache: Path, g: pd.Series, label: str) -> bool:
    path = cache / "streamstats" / f"{g['site_no']}.json"
    t0 = time.monotonic()
    errors = []
    # The home state's StreamStats region first. A gage whose basin lies mostly
    # across the line is sometimes refused there ("OUTSIDE BOUNDARY OF DATASET")
    # and served by the neighbour's region instead; the region used is recorded.
    for region in (g["home"], g["neighbor"]):
        try:
            w = ss.delineate_and_get_characteristics(
                region, float(g["latitude"]), float(g["longitude"]), timeout=180
            )
            break
        except Exception as exc:
            errors.append(f"{region}: {type(exc).__name__}: {str(exc)[:250]}")
    else:
        _write(path, {"error": " | ".join(errors)})
        logger.warning("%s %s: %s", label, g["site_no"], " | ".join(errors)[:300])
        return False
    _write(path, w.to_dict())
    logger.info(
        "%s %s %s: %.1f mi2 (catalog %s), %.0f s",
        label,
        g["site_no"],
        g["home"],
        w.characteristics["DRNAREA"].value if "DRNAREA" in w else float("nan"),
        g["drainage_area_sqmi"],
        time.monotonic() - t0,
    )
    return True


def delineate(cache: Path, limit: Optional[int], workers: int = 2) -> None:
    from concurrent.futures import ThreadPoolExecutor

    cand = pd.read_csv(cache / "candidates.csv", dtype={"site_no": str, "huc8": str})
    # Border gages first, nearest first; the straddle candidates after.
    cand = cand.sort_values(["kind", "distance_km"])
    sites = cand.drop_duplicates("site_no")
    if limit:
        sites = sites.head(limit)

    def pending(site: str) -> bool:
        p = cache / "streamstats" / f"{site}.json"
        # A failure recorded before the neighbour fallback existed is retried.
        return not p.exists() or ("error" in _read(p) and " | " not in _read(p)["error"])

    todo = [g for _, g in sites.iterrows() if pending(g["site_no"])]
    logger.info("delineate: %d to do of %d", len(todo), len(sites))
    workers = max(1, min(workers, ss.MAX_CONCURRENCY))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [
            ex.submit(_delineate_one, cache, g, f"[{i}/{len(todo)}]") for i, g in enumerate(todo, 1)
        ]
        ok = sum(f.result() for f in futs)
    logger.info("delineate: %d ok, %d failed this run", ok, len(todo) - ok)


SS_REGIONS = "https://gis.streamstats.usgs.gov/arcgis/rest/services/nss/regions/MapServer"


def region_layer(cache: Path, state: str) -> List[Tuple[List[str], Dict[str, Any]]]:
    """A state's peak-flow region polygons from the StreamStats region layer (cached)."""
    layer, grids = qa.REGION_LAYERS[state]

    def fetch() -> List[Dict[str, Any]]:
        out: List[Dict[str, Any]] = []
        offset = 0
        where = " OR ".join(f"grid_name='{g}'" for g in grids)
        while True:
            r = requests.get(
                f"{SS_REGIONS}/{layer}/query",
                params={
                    "where": where,
                    "outFields": "grid_name,GRIDCODE,Name",
                    "returnGeometry": "true",
                    "outSR": "4326",
                    "maxAllowableOffset": "0.0005",  # ~50 m, as tools/build_gage_catalog.py
                    "orderByFields": "OBJECTID",
                    "resultOffset": offset,
                    "resultRecordCount": 50,
                    "f": "json",
                },
                timeout=300,
            )
            r.raise_for_status()
            payload = r.json()
            if "error" in payload:
                raise RuntimeError(f"StreamStats regions layer {layer}: {payload['error']}")
            feats = payload.get("features") or []
            out.extend(
                {"gridcode": f["attributes"]["GRIDCODE"], "rings": f["geometry"]["rings"]}
                for f in feats
            )
            if not feats or not payload.get("exceededTransferLimit"):
                return out
            offset += len(feats)

    polys = _cached(cache / "regions" / f"{state}.json", fetch)
    return [
        (
            [c.strip().upper() for c in str(p["gridcode"]).split(",") if c.strip()],
            qa.esri_rings_to_geojson(p["rings"]),
        )
        for p in polys
    ]


def locate_state(
    cache: Path,
    state: str,
    w: ss.WatershedCharacteristics,
    pts: Sequence[Tuple[float, float]],
) -> Tuple[List[Tuple[str, float]], str]:
    """A state's peak-flow regions for a basin, as (NSS code, percent), and how found.

    Region-layer states: the basin's grid sample against the layer polygons; a basin
    in none of them takes the polygon nearest its outlet. New Mexico: ``HIGHREG``
    and Wyoming: ``WYPK_IND``, both computed only by that state's own StreamStats
    delineation. Nevada: not locatable (no region geometry anywhere).
    """
    lib = load_state(state)
    codes = set(qa.nss_region_map(lib))
    chars = {k: c.value for k, c in w.characteristics.items()}
    lat, lon = w.snap.snapped_lat, w.snap.snapped_lon
    if state in qa.REGION_LAYERS:
        layer = region_layer(cache, state)
        located: Dict[str, float] = {}
        for gcodes, geom in layer:
            gcodes = [c for c in gcodes if c in codes]
            if not gcodes:
                continue
            frac = qa.fraction_inside(pts, geom)
            if frac <= 0:
                continue
            for c in qa.resolve_overlay(state, gcodes, chars):
                located[c] = located.get(c, 0.0) + 100.0 * frac
        if located:
            return sorted(located.items()), "region layer"
        dists = [
            (qa.distance_km(lat, lon, geom), gcodes)
            for gcodes, geom in layer
            if any(c in codes for c in gcodes)
        ]
        d, gcodes = min(dists, key=lambda x: x[0])
        res = qa.resolve_overlay(state, [c for c in gcodes if c in codes], chars)
        return [(c, 100.0) for c in res], f"nearest region ({d:.0f} km)"
    if state == "NM":
        if w.region == "NM" and "HIGHREG" in chars:
            return [(qa.highreg_region(chars["HIGHREG"]), 100.0)], "HIGHREG"
        return [], "NM regions need HIGHREG, computed only on an NM delineation"
    if state == "WY":
        # NSS has no polygons for WY's WRIR 03-4107 regions (bylocation answers only
        # a different statewide set, GC2098-GC2115, everywhere). StreamStats' WY
        # delineation computes WYPK_IND, its peak-flow region index (1-6), which is
        # read here as WRIR 03-4107's region number.
        if w.region == "WY" and "WYPK_IND" in chars:
            region = str(int(chars["WYPK_IND"]))
            code = next((c for c, r in qa.nss_region_map(lib).items() if r == region), None)
            if code:
                return [(code, 100.0)], "WYPK_IND"
        return [], "WY regions need WYPK_IND, computed only on a WY delineation"
    # Nevada: no region layer, and NSS bylocation returns nothing anywhere in the state.
    return [], f"{state} regions have no geometry in NSS or the StreamStats region layers"


def locate(cache: Path) -> None:
    cand = pd.read_csv(cache / "candidates.csv", dtype={"site_no": str, "huc8": str})
    outlines = state_outlines(cache)
    for i, (_, g) in enumerate(cand.iterrows(), 1):
        src = cache / "streamstats" / f"{g['site_no']}.json"
        path = cache / "locate" / f"{g['site_no']}_{g['neighbor']}.json"
        if not src.exists() or path.exists():
            continue
        data = _read(src)
        if "error" in data:
            continue
        w = ss.WatershedCharacteristics.from_dict(data)
        if w.polygon_geojson is None:
            continue
        da_cat = g["drainage_area_sqmi"]
        if "DRNAREA" in w and pd.notna(da_cat) and da_cat > 0:
            r = w.characteristics["DRNAREA"].value / float(da_cat)
            if not 1 / DA_TOL <= r <= DA_TOL:
                continue  # mis-snapped; evaluate() lists it
        pts = qa.sample_points(w.polygon_geojson["geometry"], 400)
        out: Dict[str, Any] = {
            "frac_in_neighbor": qa.fraction_inside(pts, outlines[g["neighbor"]]),
            "frac_in_home": qa.fraction_inside(pts, outlines[g["home"]]),
        }
        try:
            for role, st in (("home", g["home"]), ("neighbor", g["neighbor"])):
                located, how = locate_state(cache, st, w, pts)
                out[role] = {"state": st, "located": located, "how": how}
        except Exception as exc:
            logger.warning("locate %s/%s: %s", g["site_no"], g["neighbor"], str(exc)[:200])
            continue
        _write(path, out)
        logger.info(
            "[%d/%d] %s %s->%s: %s / %s",
            i,
            len(cand),
            g["site_no"],
            g["home"],
            g["neighbor"],
            out["home"]["located"],
            out["neighbor"]["located"],
        )


# ----------------------------------------------------------------------------
# evaluate
# ----------------------------------------------------------------------------


def load_fits(cache: Path) -> pd.DataFrame:
    """The B17C fits, one row per gage, without failed or degenerate fits.

    ``degenerate`` (:func:`flowfreq.validation.wave_qa.b17c_degenerate`) rows are
    kept in ``b17c.csv`` and counted in the report, but never used as a reference.
    """
    fits = pd.read_csv(cache / "b17c.csv", dtype={"site_no": str}).drop_duplicates("site_no")
    ok = fits["error"].isna() if "error" in fits else pd.Series(True, index=fits.index)
    fits["degenerate"] = [
        (not k) or qa.b17c_degenerate(float(a), float(b))
        for k, a, b in zip(ok, fits["q_0.5"], fits["q_0.01"])
    ]
    return fits


def evaluate(cache: Path) -> pd.DataFrame:
    cand = pd.read_csv(cache / "candidates.csv", dtype={"site_no": str, "huc8": str})
    fits = load_fits(cache)
    fits = fits[~fits["degenerate"]].set_index("site_no")
    libs = {st: load_state(st) for st in STATES}
    rows: List[Dict[str, Any]] = []

    dropped: List[Dict[str, Any]] = []
    for _, g in cand.iterrows():
        ss_path = cache / "streamstats" / f"{g['site_no']}.json"
        if not ss_path.exists():
            continue
        data = _read(ss_path)
        if "error" in data:
            dropped.append(
                {
                    "site_no": g["site_no"],
                    "home": g["home"],
                    "reason": "delineation failed",
                    "detail": data["error"][:200],
                }
            )
            continue
        w = ss.WatershedCharacteristics.from_dict(data)
        if "DRNAREA" not in w:
            dropped.append(
                {
                    "site_no": g["site_no"],
                    "home": g["home"],
                    "reason": "no DRNAREA returned",
                    "detail": f"polygon {w.polygon_area_sq_mi or float('nan'):.0f} mi2",
                }
            )
            continue
        loc_path = cache / "locate" / f"{g['site_no']}_{g['neighbor']}.json"
        fit = fits.loc[g["site_no"]] if g["site_no"] in fits.index else None
        b17c_q: Optional[Dict[float, Tuple[float, float, float]]] = None
        if fit is not None and pd.isna(fit.get("error", np.nan)):
            b17c_q = {
                aep: (fit[f"q_{aep:g}"], fit[f"lo_{aep:g}"], fit[f"hi_{aep:g}"])
                for aep in AEPS
                if f"q_{aep:g}" in fit and pd.notna(fit[f"q_{aep:g}"])
            }
        da_ss = w.characteristics["DRNAREA"].value if "DRNAREA" in w else None
        da_cat = g["drainage_area_sqmi"]
        da_ratio = (
            float(da_ss) / float(da_cat) if (da_ss and pd.notna(da_cat) and da_cat > 0) else None
        )
        if da_ratio is None or not (1 / DA_TOL <= da_ratio <= DA_TOL):
            # The catalog coordinate snapped to the wrong stream (a 0.04 mi2 side
            # channel for an 856 mi2 river is typical): not this gage's basin.
            dropped.append(
                {
                    "site_no": g["site_no"],
                    "home": g["home"],
                    "reason": "drainage-area mismatch",
                    "detail": f"StreamStats {da_ss} vs catalog {da_cat} mi2",
                }
            )
            continue
        if not loc_path.exists():
            continue
        loc = _read(loc_path)
        frac = float(loc["frac_in_neighbor"])
        kind = "straddles" if frac > 0 else ("border" if g["kind"] == "border" else "inland")
        base = {
            "site_no": g["site_no"],
            "site_name": g["site_name"],
            "home": g["home"],
            "neighbor": g["neighbor"],
            "kind": kind,
            "distance_km": g["distance_km"],
            "frac_in_neighbor": round(frac, 3),
            "ss_region": w.region,
            "drnarea": w.characteristics["DRNAREA"].value if "DRNAREA" in w else None,
            "neighbor_how": loc["neighbor"]["how"],
            "n_peaks": None if fit is None else fit.get("n"),
            "skew_source": None if fit is None else fit.get("skew_source"),
        }
        for r in qa.compare_basin(
            libs,
            g["home"],
            g["neighbor"],
            {k: c.value for k, c in w.characteristics.items()},
            (w.snap.snapped_lat, w.snap.snapped_lon),
            [tuple(x) for x in loc["home"]["located"]],
            [tuple(x) for x in loc["neighbor"]["located"]],
            b17c_q,
            AEPS,
        ):
            rows.append({**base, **r})
    df = pd.DataFrame(rows)
    df.to_csv(cache / "evaluate.csv", index=False)
    dr = pd.DataFrame(dropped, columns=["site_no", "home", "reason", "detail"])
    dr = dr.drop_duplicates("site_no")
    dr.to_csv(cache / "dropped.csv", index=False)
    logger.info("evaluate: dropped %s", dr.groupby("reason").size().to_dict())
    logger.info(
        "evaluate: %d rows, %d site-pairs",
        len(df),
        df[["site_no", "neighbor"]].drop_duplicates().shape[0],
    )
    return df


# ----------------------------------------------------------------------------
# DAR LOOCV
# ----------------------------------------------------------------------------


def dar(cache: Path) -> pd.DataFrame:
    fits = load_fits(cache)
    fits = fits[~fits["degenerate"]]
    cat = eligible_catalog().set_index("site_no")
    net = fits.join(
        cat[["huc8", "latitude", "longitude", "drainage_area_sqmi", "regression_region"]],
        on="site_no",
    )
    net = net.dropna(subset=["huc8", "drainage_area_sqmi"])
    for aep in AEPS:
        net[f"log_q_{aep:g}"] = np.log10(net[f"q_{aep:g}"].where(net[f"q_{aep:g}"] > 0))
    # First pass with a unit exponent finds which targets have a donor at all.
    probe = qa.dar_loocv(net, {(s, AEPS[0]): 1.0 for s in net["site_no"]}, AEPS[:1])
    targets = set(probe["site_no"])
    libs = {st: load_state(st) for st in STATES}
    exps: Dict[Tuple[str, float], float] = {}
    regions_used: Dict[str, str] = {}
    for _, t in net[net["site_no"].isin(targets)].iterrows():
        lib = libs[t["state"]]
        # The catalog's gage-point region (WA, OR, ID, MT, CO, UT, AZ). WY, NM and NV
        # have no region geometry, so their targets get no exponent and are skipped.
        if not (isinstance(t["regression_region"], str) and t["regression_region"]):
            continue
        codes = [c.strip().upper() for c in t["regression_region"].split(";")]
        if t["state"] == "AZ" and len(codes) > 1:
            # 'GC1618;<2-5>': High Elevation region 1 needs the basin's mean
            # elevation, which the catalog lacks; the location's own region is used.
            codes = [c for c in codes if c != "GC1618"]
        codes = qa.resolve_overlay(t["state"], codes, {})
        weights, _ = qa.region_weights(lib, [(c, 100.0) for c in codes])
        if not weights:
            continue
        region = max(weights, key=lambda r: weights[r])
        regions_used[t["site_no"]] = region
        eq_region = "2B" if (t["state"] == "OR" and region == "2") else region
        for aep in AEPS:
            b = qa.area_exponent(lib, eq_region, aep)
            if b is not None:
                exps[(t["site_no"], aep)] = b
    out = qa.dar_loocv(net, exps, AEPS)
    out["region"] = out["site_no"].map(regions_used)
    out.to_csv(cache / "dar_loocv.csv", index=False)
    logger.info("dar: %d target-AEP rows, %d targets", len(out), out["site_no"].nunique())
    return out


# ----------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------


def _md(df: pd.DataFrame, floatfmt: str = ".3f") -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                cells.append("" if pd.isna(v) else format(v, floatfmt))
            else:
                cells.append("" if v is None or (not isinstance(v, str) and pd.isna(v)) else str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def _dominant(regions: Any) -> str:
    """The largest-weight region of a ``"r:w;r:w"`` string ("" when none)."""
    if not isinstance(regions, str) or not regions:
        return ""
    parts = [p.rsplit(":", 1) for p in regions.split(";")]
    return max(parts, key=lambda p: float(p[1]))[0]


def report(cache: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    ev = pd.read_csv(cache / "evaluate.csv", dtype={"site_no": str})
    for c in ("flag_hn", "flag_home_b17c", "flag_neighbor_b17c"):
        ev[c] = ev[c].map({True: 1, False: 0, "True": 1, "False": 0}).astype("Int64")
    ev["pair"] = ev["home"] + "-" + ev["neighbor"]
    ev["home_region"] = ev["home_regions"].map(_dominant)
    ev["neighbor_region"] = ev["neighbor_regions"].map(_dominant)
    parts: List[str] = []

    # 0. sample accounting
    sites = ev.drop_duplicates(["site_no", "neighbor"])
    acct = sites.groupby(["pair", "kind"]).size().unstack(fill_value=0).reset_index()
    acct.to_csv(out / "sample.csv", index=False)
    parts.append("## Sample: evaluated site-pairs by kind\n\n" + _md(acct))
    dr = pd.read_csv(cache / "dropped.csv", dtype={"site_no": str})
    dr.to_csv(out / "dropped_sites.csv", index=False)
    parts.append(
        "## Candidate gages dropped before evaluation\n\n"
        + _md(dr.groupby(["home", "reason"]).size().reset_index(name="n"))
    )
    fits = load_fits(cache)
    err = fits["error"].notna() if "error" in fits else pd.Series(False, index=fits.index)
    fa = (
        fits.assign(failed=err, degenerate_only=fits["degenerate"] & ~err)
        .groupby("state")
        .agg(
            gages=("site_no", "size"),
            failed=("failed", "sum"),
            degenerate=("degenerate_only", "sum"),
        )
        .reset_index()
    )
    fa.to_csv(out / "b17c_fits.csv", index=False)
    parts.append("## B17C fits on the network\n\n" + _md(fa))

    # 1. border discontinuity, border and straddling basins only
    cross = ev[ev["kind"].isin(["border", "straddles"])]
    both = cross.dropna(subset=["diff_log"]).copy()
    both["absd"] = both["diff_log"].abs()
    both["pos"] = (both["diff_log"] > 0).astype(int)
    by_pair = (
        both[both["aep"].isin(REPORT_AEPS)]
        .groupby(["pair", "aep"])
        .agg(
            n=("diff_log", "size"),
            median_diff_log=("diff_log", "median"),
            median_abs=("absd", "median"),
            max_abs=("absd", "max"),
            share_home_higher=("pos", "mean"),
            median_thr=("thr_log", "median"),
            flagged=("flag_hn", "sum"),
        )
        .reset_index()
    )
    by_pair["median_ratio"] = 10.0 ** by_pair["median_diff_log"]
    # A pair-level screen, descriptive rather than a test: one state's equations sit
    # on the same side of the other's at three or more sites, by more than 0.1 log10
    # (26 %) at the median. No single site need be flagged for this to be a
    # consistent discontinuity at the line.
    by_pair["systematic"] = (
        (by_pair["n"] >= 3)
        & (by_pair["median_diff_log"].abs() > SYSTEMATIC_LOG)
        & ((by_pair["share_home_higher"] >= 0.75) | (by_pair["share_home_higher"] <= 0.25))
    )
    by_pair.to_csv(out / "border_by_pair.csv", index=False)
    parts.append("## Border discontinuity by state pair and AEP\n\n" + _md(by_pair))

    # 2. flagged site-pairs (border and straddling only)
    fl = both[both["flag_hn"] == 1]
    flagged = (
        fl.groupby(["site_no", "home", "neighbor"])
        .agg(
            site_name=("site_name", "first"),
            kind=("kind", "first"),
            km=("distance_km", "first"),
            frac_nb=("frac_in_neighbor", "first"),
            da=("drnarea", "first"),
            home_regions=("home_regions", "first"),
            neighbor_regions=("neighbor_regions", "first"),
            neighbor_how=("neighbor_how", "first"),
            proxies=("neighbor_proxies", "first"),
            aeps=("aep", lambda s: ",".join(f"{a:g}" for a in sorted(s, reverse=True))),
            worst_ratio=("diff_log", lambda s: 10.0 ** s.loc[s.abs().idxmax()]),
            worst_thr_ratio=("thr_log", lambda s: 10.0 ** s.max()),
        )
        .reset_index()
        .sort_values(["home", "neighbor", "site_no"])
    )
    flagged.to_csv(out / "flagged_pairs.csv", index=False)
    parts.append(f"## Flagged site-pairs ({len(flagged)})\n\n" + _md(flagged, ".2f"))

    # 3. not evaluable, by side, state and reason
    unev = cross[cross["aep"] == 0.01]
    reasons = []
    for role in ("home", "neighbor"):
        r = unev[unev[f"{role}_why"].fillna("") != ""]
        for _, x in r.iterrows():
            why = str(x[f"{role}_why"])
            m = re.search(r"\[(.*)\]", why)
            kind = (
                f"missing {m.group(1)}"
                if m and why.startswith("missing")
                else ("Region 2A/2B needs ELEV" if "mean elevation" in why else why[:70])
            )
            reasons.append(
                {
                    "site_no": x["site_no"],
                    "pair": x["pair"],
                    "side": role,
                    "state": x[role],
                    "reason": kind,
                }
            )
    rs = pd.DataFrame(reasons, columns=["site_no", "pair", "side", "state", "reason"])
    rs.to_csv(out / "unevaluable.csv", index=False)
    parts.append(
        "## Border site-pairs not evaluable at AEP 0.01\n\n"
        + _md(rs.groupby(["side", "state", "reason"]).size().reset_index(name="n"))
    )

    # 4. equations vs B17C, every evaluated basin (inland included), one row per site
    rows = []
    for role in ("home", "neighbor"):
        sub = ev.dropna(subset=[f"{role}_vs_b17c_log"])
        if role == "home":
            sub = sub.drop_duplicates(["site_no", "aep"])
        rows.append(
            pd.DataFrame(
                {
                    "role": role,
                    "state": sub[role],
                    "region": sub[f"{role}_region"],
                    "site_no": sub["site_no"],
                    "kind": sub["kind"],
                    "aep": sub["aep"],
                    "log_est": sub[f"{role}_vs_b17c_log"],
                    "log_obs": 0.0,
                    "flag": sub[f"flag_{role}_b17c"],
                }
            )
        )
    vb = pd.concat(rows, ignore_index=True)
    vb.to_csv(out / "equation_vs_b17c.csv", index=False)
    rep = vb[vb["aep"].isin(REPORT_AEPS)]
    s = qa.summarize_log_errors(rep, ["role", "state", "aep"])
    s = s.merge(
        rep.groupby(["role", "state", "aep"])["flag"].sum().reset_index(name="flagged"),
        on=["role", "state", "aep"],
    )
    s.to_csv(out / "equation_vs_b17c_by_state.csv", index=False)
    parts.append("## Equations vs. gage B17C\n\n" + _md(s))
    fb = (
        vb[vb["flag"] == 1]
        .merge(ev[["site_no", "site_name"]].drop_duplicates("site_no"), on="site_no")
        .groupby(["role", "state", "region", "site_no"])
        .agg(
            site_name=("site_name", "first"),
            kind=("kind", "first"),
            aeps=("aep", lambda a: ",".join(f"{x:g}" for x in sorted(set(a), reverse=True))),
            worst_ratio=("log_est", lambda e: 10.0 ** e.loc[e.abs().idxmax()]),
        )
        .reset_index()
        .sort_values(["role", "state", "site_no"])
    )
    fb.to_csv(out / "flagged_vs_b17c.csv", index=False)
    parts.append(
        f"## Equation vs. B17C disagreements beyond the threshold ({len(fb)})\n\n"
        "`worst_ratio` is equation / B17C at the AEP with the largest log difference.\n\n"
        + _md(fb, ".2f")
    )
    home = rep[rep["role"] == "home"]
    s2 = qa.summarize_log_errors(home, ["state", "region", "aep"])
    s2.to_csv(out / "regression_by_region.csv", index=False)
    parts.append(
        "## Home-state equations vs. B17C by dominant region (AEP 0.01)\n\n"
        + _md(s2[s2["aep"] == 0.01])
    )

    # 5. DAR LOOCV
    if (cache / "dar_loocv.csv").exists():
        d = pd.read_csv(cache / "dar_loocv.csv", dtype={"site_no": str})
        d = d[np.isfinite(d["log_obs"]) & (d["log_obs"] >= 0) & np.isfinite(d["log_est"])]
        s3 = qa.summarize_log_errors(d[d["aep"].isin(REPORT_AEPS)], ["state", "aep"])
        s3.to_csv(out / "dar_loocv_by_state.csv", index=False)
        s4 = qa.summarize_log_errors(d, ["state", "region", "aep"])
        s4.to_csv(out / "dar_loocv_by_region.csv", index=False)
        s5 = qa.summarize_log_errors(d, ["aep"])
        parts.append(
            "## DAR transposition LOOCV, all targets by AEP\n\n"
            + _md(s5)
            + "\n\n## DAR transposition LOOCV by state\n\n"
            + _md(s3)
        )
    (out / "tables.md").write_text("\n\n".join(parts) + "\n", encoding="utf-8")
    logger.info("report written to %s", out)


def vignette(cache: Path, path: Path) -> None:
    """Capture one evaluated basin per state for ``docs/vignettes/wave1_wave2_ungaged.py``.

    Per state, the home basin that evaluates at AEP 0.01 with every input in range,
    preferring one that spans two regions (to show the area weighting), then the
    longest record. Only the characteristics its equations use are kept.
    """
    ev = pd.read_csv(cache / "evaluate.csv", dtype={"site_no": str})
    ev = ev[(ev["aep"] == 0.01) & ev["home_q"].notna() & ev["home_oor"].isna()]
    libs = {st: load_state(st) for st in STATES}
    sites: Dict[str, Any] = {}
    for st in STATES:
        sub = ev[ev["home"] == st].drop_duplicates("site_no")
        if sub.empty:
            continue
        sub = sub.assign(n_reg=sub["home_regions"].str.count(":"))
        g = sub.sort_values(["n_reg", "n_peaks"], ascending=False).iloc[0]
        w = ss.WatershedCharacteristics.from_dict(
            _read(cache / "streamstats" / f"{g['site_no']}.json")
        )
        loc = _read(next((cache / "locate").glob(f"{g['site_no']}_*.json")))
        weights, _ = qa.region_weights(libs[st], [tuple(x) for x in loc["home"]["located"]])
        chars_all, _ = qa.prepare_characteristics(
            {k: c.value for k, c in w.characteristics.items()},
            w.snap.snapped_lat,
            w.snap.snapped_lon,
        )
        need = set()
        for region in weights:
            for r in ["2A", "2B"] if (st == "OR" and region == "2") else [region]:
                need |= {
                    v.code for e in libs[st].equations if e.region_code == r for v in e.variables
                }
        if st == "OR" and "2" in weights:
            need.add("ELEV")
        fit = (
            pd.read_csv(cache / "b17c.csv", dtype={"site_no": str})
            .set_index("site_no")
            .loc[g["site_no"]]
        )
        sites[st] = {
            "site_no": g["site_no"],
            "site_name": g["site_name"],
            "outlet": [round(w.snap.snapped_lat, 6), round(w.snap.snapped_lon, 6)],
            "streamstats_region": w.region,
            "captured": (w.provenance.requested_at_utc[:10] if w.provenance else ""),
            "region_weights": {k: round(v, 4) for k, v in weights.items()},
            "characteristics": {k: chars_all[k] for k in sorted(need) if k in chars_all},
            "b17c_q1pct": round(float(fit["q_0.01"]), 1),
            "b17c_n": int(fit["n"]),
        }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(sites, indent=2) + "\n", encoding="utf-8")
    logger.info("vignette sites written to %s (%d states)", path, len(sites))


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "stages",
        nargs="+",
        choices=[
            "all",
            "states",
            "candidates",
            "peaks",
            "b17c",
            "delineate",
            "locate",
            "evaluate",
            "dar",
            "report",
            "vignette",
        ],
    )
    ap.add_argument("--cache", type=Path, default=REPO_ROOT / ".cache" / "wave_qa")
    ap.add_argument("--out", type=Path, default=REPO_ROOT / "docs" / "wave_qa")
    ap.add_argument("--border-km", type=float, default=20.0)
    ap.add_argument("--pace", type=float, default=4.0, help="seconds between Water Data calls")
    ap.add_argument("--workers", type=int, default=4, help="processes for the B17C fits")
    ap.add_argument("--limit", type=int, default=None, help="delineate at most this many sites")
    ap.add_argument("--ss-workers", type=int, default=2, help="parallel StreamStats delineations")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logging.getLogger("flowfreq").setLevel(logging.ERROR)
    stages = args.stages
    if "all" in stages:
        stages = [
            "states",
            "candidates",
            "b17c",
            "delineate",
            "locate",
            "evaluate",
            "dar",
            "report",
            "vignette",
        ]
    for st in stages:
        if st == "states":
            state_outlines(args.cache)
        elif st == "candidates":
            candidates(args.cache, args.border_km)
        elif st == "peaks":
            state_peaks(args.cache, args.pace)
        elif st == "b17c":
            b17c(args.cache, args.pace, args.workers)
        elif st == "delineate":
            delineate(args.cache, args.limit, args.ss_workers)
        elif st == "locate":
            locate(args.cache)
        elif st == "evaluate":
            evaluate(args.cache)
        elif st == "dar":
            dar(args.cache)
        elif st == "report":
            report(args.cache, args.out)
        elif st == "vignette":
            vignette(args.cache, REPO_ROOT / "docs" / "vignettes" / "wave1_wave2_sites.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
