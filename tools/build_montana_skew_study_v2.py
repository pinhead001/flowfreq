"""Montana regional skew study v2 (B-WLS/B-GLS) -- inputs to USGS standard.

NOT a USGS study. Northwest studies plan, tasks 5-7; see
``docs/MONTANA_REGIONAL_SKEW_V2.md``. v1 (``build_montana_skew_study.py``,
``docs/MONTANA_REGIONAL_SKEW_PROVISIONAL.md``) stays reproducible; v2 reuses its
download, screening, EMA and redundancy helpers and closes the gaps v1 lists
under Limitations, the way the Pacific Northwest study (SIR 2016-5083 app. B)
and the Arizona study (SIR 2014-5211 app. 1) were built:

- **Neighbouring gages in shared basins.** Wyoming and Idaho gages in HUC4
  1002 (Missouri headwaters), 1007 (Yellowstone), 1008 (Bighorn) and 1701
  (Kootenai-Pend Oreille-Spokane), from the same bulk OGC queries.
- **Regulation screen.** v1's peak-code cut (record before the first peak
  coded 5, 6 or C), plus the packaged GAGES-II screen
  (:func:`flowfreq.regulation.load_screen`): a gage classed ``regulated`` or
  ``urban`` is excluded unless its record was already cut at a coded
  regulation peak -- uncoded regulation has no known start year.
- **Basin centroids and characteristics** from StreamStats
  (:func:`flowfreq.streamstats.batch_get_characteristics`, cached): mean basin
  elevation ``ELEV`` and mean annual precipitation ``PRECIP``; the centroid is
  the area-weighted centroid of the delineated polygon. Redundancy screening
  and cross-correlation distances use centroids.
- **Distance to the Continental Divide** (``div_dist_mi``, positive east): from
  the basin centroid to the boundary between WBD HUC2 region 17 and regions
  10/09, the Rocky Mountain Front proxy v1's residuals point to.

Stages are cached under ``--cache``; an interrupted run resumes. OGC requests
are paced (1000/hour without ``USGS_API_KEY``); StreamStats runs at
``--concurrency`` (<= 3).

Usage::

    python tools/build_montana_skew_study_v2.py --cache .cache/mt_skew_v2 inputs
    python tools/build_montana_skew_study_v2.py --cache .cache/mt_skew_v2 fit
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_montana_skew_study as v1  # noqa: E402

from flowfreq.regulation import load_screen  # noqa: E402
from flowfreq.skew_study import distance_miles  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "data" / "skew_study" / "montana_v2_inputs.csv"
RESULTS = ROOT / "data" / "skew_study" / "montana_v2_results.json"
VERSION = "montana-provisional-v2"

STATES = {"MT": "30", "WY": "56", "ID": "16"}
#: Neighbour-state gages are kept only in these HUC4s (basins shared with MT).
NEIGHBOUR_HUC4 = ("1002", "1007", "1008", "1701")
WBD_HUC2 = "https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer/1/query"
SS_CODES = ("DRNAREA", "ELEV", "PRECIP")

log = logging.getLogger("mt_skew_v2")


# ----------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------


def fetch_state(cache: Path, state: str, pause: float) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """(peaks, sites) for one state from the bulk OGC queries, cached per state."""
    code = STATES[state]
    feats = v1._paged(
        cache,
        f"peaks_{state}",
        f"{v1.API}/peaks/items",
        {
            "f": "json",
            "state_code": code,
            "parameter_code": "00060",
            "skipGeometry": "true",
            "limit": 20000,
            "properties": "monitoring_location_id,water_year,value,qualifier,time",
        },
        pause,
    )
    rows = []
    for f in feats:
        p = f["properties"]
        site = str(p["monitoring_location_id"]).replace("USGS-", "")
        rows.append(
            {
                "site_no": site,
                "water_year": p.get("water_year"),
                "peak_flow_cfs": pd.to_numeric(p.get("value"), errors="coerce"),
                "qualification_code": v1.qualifiers_to_codes(
                    p.get("qualifier"), site, p.get("water_year")
                ),
            }
        )
    peaks = pd.DataFrame(rows).dropna(subset=["water_year", "peak_flow_cfs"])
    peaks["water_year"] = peaks["water_year"].astype(int)
    sfeats = v1._paged(
        cache,
        f"sites_{state}",
        f"{v1.API}/monitoring-locations/items",
        {
            "f": "json",
            "state_code": code,
            "site_type_code": "ST",
            "limit": 10000,
            "properties": "monitoring_location_number,monitoring_location_name,"
            "hydrologic_unit_code,drainage_area,altitude",
        },
        pause,
    )
    srows = []
    for f in sfeats:
        p = f["properties"]
        lon, lat = ((f.get("geometry") or {}).get("coordinates") or [None, None])[:2]
        srows.append(
            {
                "site_no": str(p.get("monitoring_location_number")),
                "state": state,
                "name": p.get("monitoring_location_name"),
                "huc": str(p.get("hydrologic_unit_code") or ""),
                "drainage_area": p.get("drainage_area"),
                "altitude": p.get("altitude"),
                "lat": lat,
                "lon": lon,
            }
        )
    sites = pd.DataFrame(srows).drop_duplicates("site_no")
    if state != "MT":
        sites = sites[sites["huc"].str[:4].isin(NEIGHBOUR_HUC4)]
    return peaks[peaks["site_no"].isin(sites["site_no"])], sites


# ----------------------------------------------------------------------------
# StreamStats basins and the Continental Divide
# ----------------------------------------------------------------------------


def polygon_centroid(geometry: Dict[str, Any]) -> Tuple[float, float]:
    """Area-weighted (lat, lon) centroid of a GeoJSON Polygon/MultiPolygon (outer rings)."""
    polys = (
        geometry["coordinates"] if geometry["type"] == "MultiPolygon" else [geometry["coordinates"]]
    )
    a_sum = cx = cy = 0.0
    for poly in polys:
        ring = np.asarray(poly[0], dtype=float)
        x, y = ring[:, 0], ring[:, 1]
        cross = x[:-1] * y[1:] - x[1:] * y[:-1]
        a = cross.sum() / 2.0
        if a == 0:
            continue
        cx += ((x[:-1] + x[1:]) * cross).sum() / 6.0
        cy += ((y[:-1] + y[1:]) * cross).sum() / 6.0
        a_sum += a
    if a_sum == 0:
        raise ValueError("degenerate polygon")
    return cy / a_sum, cx / a_sum


def streamstats_basins(
    cache: Path, gages: pd.DataFrame, concurrency: int
) -> Dict[str, Dict[str, Any]]:
    """Per gage: centroid, ELEV, PRECIP, DRNAREA from StreamStats (cached JSON per gage)."""
    from flowfreq.streamstats import StreamStatsCache, batch_get_characteristics

    d = cache / "streamstats"
    d.mkdir(parents=True, exist_ok=True)
    out: Dict[str, Dict[str, Any]] = {}
    todo = []
    for r in gages.itertuples():
        f = d / f"{r.site_no}.json"
        if f.exists():
            out[r.site_no] = json.loads(f.read_text(encoding="utf-8"))
        else:
            todo.append((r.site_no, r.state, float(r.lat), float(r.lon)))
    log.info("StreamStats: %d cached, %d to delineate", len(out), len(todo))
    ss_cache = StreamStatsCache(cache / "streamstats_raw")
    for i in range(0, len(todo), 12):
        chunk = todo[i : i + 12]
        res, errs = batch_get_characteristics(
            chunk,
            concurrency=concurrency,
            cache=ss_cache,
            characteristic_codes=list(SS_CODES),
            validate_regions=False,
            timeout=120.0,
        )
        for site, _, _, _ in chunk:
            rec: Dict[str, Any]
            if site in res:
                w = res[site]
                lat_c, lon_c = polygon_centroid(w.polygon_geojson["geometry"])
                rec = {"lat_cent": lat_c, "lon_cent": lon_c}
                for code in SS_CODES:
                    c = w.characteristics.get(code)
                    rec[code] = None if c is None else float(c.value)
            else:
                rec = {"error": errs.get(site, "unknown")[:300]}
            (d / f"{site}.json").write_text(json.dumps(rec), encoding="utf-8")
            out[site] = rec
        log.info("StreamStats %d/%d", min(i + 12, len(todo)), len(todo))
    return out


def divide_points(cache: Path) -> np.ndarray:
    """(lat, lon) points along the HUC2 17 / (10, 09) boundary, from the WBD (cached)."""
    f = cache / "divide_points.json"
    if f.exists():
        return np.asarray(json.loads(f.read_text()), dtype=float)
    rings: Dict[str, List[np.ndarray]] = {}
    for huc2 in ("17", "10", "09"):
        r = requests.get(
            WBD_HUC2,
            params={
                "where": f"huc2='{huc2}'",
                "outFields": "huc2",
                "returnGeometry": "true",
                "outSR": "4326",
                "maxAllowableOffset": "0.005",
                "f": "json",
            },
            timeout=600,
        )
        r.raise_for_status()
        geo = r.json()["features"][0]["geometry"]["rings"]
        rings[huc2] = [np.asarray(x, dtype=float) for x in geo]
    west = np.concatenate(rings["17"])
    east = np.concatenate(rings["10"] + rings["09"])
    # Region 17 vertices within ~3 km of a region 10/09 vertex, inside the study box.
    box = (west[:, 1] > 41.5) & (west[:, 1] < 49.5) & (west[:, 0] > -117.5) & (west[:, 0] < -108)
    west = west[box]
    ebox = (east[:, 1] > 41.0) & (east[:, 1] < 50.0) & (east[:, 0] > -118) & (east[:, 0] < -107.5)
    east = east[ebox]
    keep = []
    for i in range(0, len(west), 2000):
        w = west[i : i + 2000]
        dd = np.hypot(
            (w[:, None, 0] - east[None, :, 0]) * math.cos(math.radians(46)),
            w[:, None, 1] - east[None, :, 1],
        )
        keep.append(w[dd.min(axis=1) < 0.03])
    pts = np.concatenate(keep)[:, ::-1]  # (lat, lon)
    f.write_text(json.dumps(pts.tolist()))
    log.info("Continental Divide: %d boundary points", len(pts))
    return pts


def divide_distance(
    lat: np.ndarray, lon: np.ndarray, huc2: Sequence[str], pts: np.ndarray
) -> np.ndarray:
    """Signed great-circle miles to the nearest divide point; positive east (HUC2 10/09)."""
    out = np.empty(len(lat))
    for i, (la, lo) in enumerate(zip(lat, lon)):
        d = distance_miles(np.r_[la, pts[:, 0]], np.r_[lo, pts[:, 1]])[0, 1:]
        out[i] = float(d.min()) * (-1.0 if str(huc2[i]) == "17" else 1.0)
    return out


# ----------------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------------


def build_inputs(args: argparse.Namespace) -> pd.DataFrame:
    peaks_l, sites_l = [], []
    for st in STATES:
        p, s = fetch_state(args.cache, st, args.pause)
        peaks_l.append(p)
        sites_l.append(s)
    peaks = pd.concat(peaks_l)
    peaks = peaks[peaks["water_year"] <= args.last_wy]
    sites = pd.concat(sites_l).drop_duplicates("site_no").set_index("site_no")
    log.info("%d peaks at %d sites", len(peaks), peaks["site_no"].nunique())

    screened = v1.screen(peaks, args.min_years)
    ema = v1.run_ema(args.cache, screened)
    reg = load_screen()
    rows = []
    for site, v in sorted(screened.items()):
        meta = sites.loc[site]
        row: Dict[str, Any] = {
            "site_no": site,
            "state": meta["state"],
            "name": meta["name"],
            "huc": meta["huc"],
            "lat": meta["lat"],
            "lon": meta["lon"],
            "drainage_area": meta["drainage_area"],
            "altitude": meta["altitude"],
            "regulated_from_wy": v["regulated_from"],
            "n_pre_regulation": v["n"],
            "reg_class": reg["site_class"].get(site, "not in GAGES-II"),
        }
        e = ema.get(site, {})
        for k in ("skew", "skew_mse", "prl", "n_systematic", "n_peaks", "n_low_outliers",
                  "yb", "ye", "ph", "at_site_option"):  # fmt: skip
            row[k] = e.get(k)
        reason = v["reason"] or e.get("error", "")
        if not reason and not (pd.notna(meta["drainage_area"]) and meta["drainage_area"] > 0):
            reason = "no drainage area"
        if not reason and str(meta["huc"])[:2] not in ("09", "10", "17"):
            reason = f"HUC {meta['huc']} outside the study basins"
        if (
            not reason
            and row["reg_class"] in ("regulated", "urban")
            and v["regulated_from"] is None
        ):
            reason = f"GAGES-II {row['reg_class']}, no coded regulation start"
        row["excluded"] = reason
        rows.append(row)
    df = pd.DataFrame(rows)

    cand = df[df["excluded"] == ""]
    ss = streamstats_basins(args.cache, cand, args.concurrency)
    for col in ("lat_cent", "lon_cent", *SS_CODES):
        df[col] = df["site_no"].map(lambda s, c=col: ss.get(s, {}).get(c))
    bad = (df["excluded"] == "") & (df["lat_cent"].isna() | df["ELEV"].isna())
    df.loc[bad, "excluded"] = (
        "no StreamStats basin ("
        + df.loc[bad, "site_no"].map(
            lambda s: str(ss.get(s, {}).get("error", "ELEV unavailable"))[:80]
        )
        + ")"
    )

    ok = df[df["excluded"] == ""].reset_index(drop=True)
    pts = divide_points(args.cache)
    df["div_dist_mi"] = np.nan
    idx = df.index[df["excluded"] == ""]
    df.loc[idx, "div_dist_mi"] = divide_distance(
        ok["lat_cent"].to_numpy(float),
        ok["lon_cent"].to_numpy(float),
        ok["huc"].str[:2].tolist(),
        pts,
    )
    ok = df[df["excluded"] == ""].reset_index(drop=True)
    centroids = ok.rename(
        columns={"lat": "lat_gage", "lon": "lon_gage", "lat_cent": "lat", "lon_cent": "lon"}
    )
    for s in v1.redundancy(centroids):
        df.loc[df["site_no"] == s, "excluded"] = "redundant (nested, DAR <= 5)"
    INPUTS.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(INPUTS, index=False, lineterminator="\n")
    used = df[df["excluded"] == ""]
    log.info("wrote %s: %d sites, %d used (%s)", INPUTS.name, len(df), len(used),
             used["state"].value_counts().to_dict())  # fmt: skip
    return df


#: Rocky Mountain Front indicator: basin centroid east of the Divide, within this many miles.
FRONT_MILES = 60.0


def fit(args: argparse.Namespace) -> Dict[str, Any]:
    """B-WLS/B-GLS models on the v2 inputs; writes RESULTS."""
    df = pd.read_csv(INPUTS, dtype={"site_no": str, "huc": str})
    used = df[df["excluded"].isna() | (df["excluded"] == "")].reset_index(drop=True)
    peaks = pd.concat([fetch_state(args.cache, st, args.pause)[0] for st in STATES])
    peaks = peaks[peaks["water_year"] <= args.last_wy]
    screened = v1.screen(peaks[peaks["site_no"].isin(used["site_no"])], args.min_years)
    geo = used.rename(
        columns={"lat": "lat_gage", "lon": "lon_gage", "lat_cent": "lat", "lon_cent": "lon"}
    )
    cmodel = v1.correlation_model(screened, geo, args.min_concurrent)
    corr = v1.skew_correlation_matrix(
        geo["lat"], geo["lon"], geo["prl"], geo["yb"], geo["ye"], geo["ph"], cmodel, args.kappa
    )
    g, p = used["skew"].to_numpy(float), used["prl"].to_numpy(float)
    const = v1.fit_regional_skew(g, p, correlation=corr, kappa=args.kappa)
    used["div_100mi"] = used["div_dist_mi"].astype(float) / 100.0
    used["elev_kft"] = used["ELEV"].astype(float) / 1000.0
    used["precip_10in"] = used["PRECIP"].astype(float) / 10.0
    used["front"] = ((used["div_dist_mi"] > 0) & (used["div_dist_mi"] <= FRONT_MILES)).astype(int)
    used["east"] = (used["huc"].str[:2] != "17").astype(int)
    candidates = {
        "DIVIDE": ["div_100mi"],
        "ELEV": ["elev_kft"],
        "PRECIP": ["precip_10in"],
        "FRONT": ["front"],
        "EAST": ["east"],
        "FRONT+ELEV": ["front", "elev_kft"],
    }
    models: Dict[str, Any] = {"CONSTANT": v1._model_dict(const)}
    eligible = {"CONSTANT": const.avp_new}
    for name, cols in candidates.items():
        if used[cols].isna().any().any():
            models[name] = {"skipped": "missing values"}
            continue
        m = v1.fit_regional_skew(
            g, p, used[cols].to_numpy(float), cols, correlation=corr, kappa=args.kappa,
            _constant_sigma2=const.sigma2,
        )  # fmt: skip
        md = v1._model_dict(m)
        t = np.abs(np.asarray(m.beta[1:]) / np.asarray(m.beta_sd[1:]))
        md["t_added"] = [round(float(x), 2) for x in t]
        models[name] = md
        if np.all(t >= 2.0):
            eligible[name] = m.avp_new
    selected = min(eligible, key=lambda k: eligible[k])
    # The table row covers MT east of the Divide (HUC2 10/09); west of it the
    # verified PNW study applies. So its constant is fitted on those gages only.
    east = np.flatnonzero(used["east"].to_numpy() == 1)
    for label, cols in (("CONSTANT_EAST_ONLY", []), ("PRECIP_EAST_ONLY", ["precip_10in"])):
        x = used.loc[east, cols].to_numpy(float) if cols else None
        m = v1.fit_regional_skew(
            g[east], p[east], x, cols or None, correlation=corr[np.ix_(east, east)],
            kappa=args.kappa,
        )  # fmt: skip
        models[label] = v1._model_dict(m)
    top = np.argsort(-const.influence)[:15]
    results: Dict[str, Any] = {
        "version": VERSION,
        "status": "provisional -- NOT a USGS study; returned only with allow_provisional=True",
        "last_water_year": args.last_wy,
        "n_sites": len(used),
        "n_by_state": used["state"].value_counts().to_dict(),
        "n_by_basin": used["huc"].str[:2].value_counts().to_dict(),
        "screen": {"min_years": args.min_years, "codes_removed": sorted(v1.SCREEN_CODES)},
        "front_miles": FRONT_MILES,
        "correlation_model": cmodel.__dict__,
        "kappa": args.kappa,
        "models": models,
        "selected": selected,
        "selection_rule": "smallest AVPnew among models whose added parameters have |beta/SD| >= 2",
        "high_influence": [
            {"site_no": used.loc[i, "site_no"], "skew": float(g[i]), "prl": float(p[i]),
             "div_dist_mi": float(used.loc[i, "div_dist_mi"])}
            for i in top
        ],  # fmt: skip
        "exclusions": df["excluded"]
        .fillna("used")
        .replace("", "used")
        .str.split(" \\(")
        .str[0]
        .value_counts()
        .to_dict(),
    }
    RESULTS.write_text(json.dumps(results, indent=1, default=float) + "\n", encoding="utf-8")
    log.info("selected %s; wrote %s", selected, RESULTS.name)
    return results


def main(argv: Sequence[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("stage", choices=["inputs", "fit"])
    ap.add_argument("--cache", type=Path, default=ROOT / ".cache" / "mt_skew_v2")
    ap.add_argument("--pause", type=float, default=4.0)
    ap.add_argument("--min-years", type=int, default=25)
    ap.add_argument("--min-concurrent", type=int, default=50)
    ap.add_argument("--last-wy", type=int, default=2025)
    ap.add_argument("--concurrency", type=int, default=3)
    ap.add_argument("--kappa", type=float, default=v1.DEFAULT_KAPPA)
    args = ap.parse_args(list(argv))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    logging.getLogger("flowfreq").setLevel(logging.ERROR)
    args.cache.mkdir(parents=True, exist_ok=True)
    if args.stage == "inputs":
        build_inputs(args)
        return 0
    res = fit(args)
    print(json.dumps({k: res["models"][k] for k in res["models"]}, indent=1, default=float))
    print("selected:", res["selected"])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
