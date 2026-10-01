"""Provisional Montana regional skew study (B-WLS/B-GLS) -- inputs and fit.

NOT a USGS study. See ``docs/MONTANA_REGIONAL_SKEW_PROVISIONAL.md``.

Stages, each cached under ``--cache`` so an interrupted run resumes:

1. ``peaks``: every Montana discharge peak from the USGS Water Data OGC API
   ``peaks`` collection (``state_code=30``), a handful of paged requests.
2. ``sites``: Montana stream monitoring locations (HUC, drainage area,
   altitude, location) from ``monitoring-locations``.
3. ``screen``: per site, the record before the first peak coded 5, 6 or C
   (regulation/diversion/urbanization evidence,
   :data:`flowfreq.peak_codes.REGULATION_EVIDENCE_CODES` and
   ``ALTERATION_EVIDENCE_CODES``) -- the "pre-regulation part of the record"
   rule of SIR 2016-5083 p. 5 -- and at least ``--min-years`` peaks in it.
4. ``ema``: station skew, its MSE and pseudo record length from flowfreq's own
   EMA with MGBT (:func:`flowfreq.skew_study.station_skew`), cached per site.
5. ``redundancy``: of gage pairs with drainage-area ratio <= 5 and
   standardized gage distance ``D / sqrt((DA_i + DA_j)/2) <= 0.5`` (SIR
   2016-5083 eqs. B3-B4, with gage locations standing in for basin centroids),
   the one with the shorter PRL is dropped.
6. ``correlation``: a Fisher-Z correlation-distance model (eq. B11) fitted to
   concurrent log peaks of site pairs with at least ``--min-concurrent`` common
   years.
7. ``fit``: B-WLS/B-GLS models (CONSTANT, and with candidate explanatory
   variables), written to ``--results``.

The OGC API allows 1000 requests/hour per IP without a key. Requests go through
:func:`flowfreq.waterdata.request` (``USGS_API_KEY``, 429/503 backoff) and are
paced (``--pause``); a full run makes three.

Usage::

    python tools/build_montana_skew_study.py --cache .cache/mt_skew
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from flowfreq import waterdata
from flowfreq.peak_codes import ALTERATION_EVIDENCE_CODES, REGULATION_EVIDENCE_CODES, parse_codes
from flowfreq.peak_sources import qualifiers_to_codes
from flowfreq.skew_study import (
    DEFAULT_KAPPA,
    CorrelationDistanceModel,
    distance_miles,
    fit_regional_skew,
    skew_correlation_matrix,
    station_skew,
)

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "data" / "skew_study" / "montana_provisional_v1_inputs.csv"
RESULTS = ROOT / "data" / "skew_study" / "montana_provisional_v1_results.json"
API = "https://api.waterdata.usgs.gov/ogcapi/v1/collections"
VERSION = "montana-provisional-v1"

log = logging.getLogger("mt_skew")
SCREEN_CODES = REGULATION_EVIDENCE_CODES | ALTERATION_EVIDENCE_CODES


# ----------------------------------------------------------------------------
# Paced, resumable OGC download
# ----------------------------------------------------------------------------


def _get(url: str, params: Optional[Dict[str, Any]], pause: float) -> Dict[str, Any]:
    """One paced request through :func:`flowfreq.waterdata.request`.

    That carries the ``USGS_API_KEY`` and the 429/503 backoff; the pause on
    top keeps a keyless run well inside 1000 requests/hour.
    """
    time.sleep(pause)
    return dict(waterdata.request(url, params, timeout=600).json())


def _paged(
    cache: Path, name: str, url: str, params: Dict[str, Any], pause: float
) -> List[Dict[str, Any]]:
    """All features of a collection query; each page cached, resume from the last."""
    d = cache / name
    d.mkdir(parents=True, exist_ok=True)
    pages = sorted(d.glob("page_*.json"))
    feats: List[Dict[str, Any]] = []
    nxt: Optional[str] = url
    q: Optional[Dict[str, Any]] = params
    for p in pages:
        payload = json.loads(p.read_text(encoding="utf-8"))
        feats.extend(payload["features"])
        nxt, q = payload.get("next"), None
    i = len(pages)
    while nxt:
        payload = _get(nxt, q, pause)
        link = next(
            (lk["href"] for lk in payload.get("links", []) if lk.get("rel") == "next"), None
        )
        page = {"features": payload.get("features") or [], "next": link}
        (d / f"page_{i:04d}.json").write_text(json.dumps(page), encoding="utf-8")
        feats.extend(page["features"])
        log.info("%s: page %d, %d features so far", name, i, len(feats))
        nxt, q, i = link, None, i + 1
    return feats


def fetch_peaks(cache: Path, pause: float) -> pd.DataFrame:
    feats = _paged(
        cache,
        "peaks",
        f"{API}/peaks/items",
        {
            "f": "json",
            "state_code": "30",
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
                "qualification_code": qualifiers_to_codes(
                    p.get("qualifier"), site, p.get("water_year")
                ),
            }
        )
    df = pd.DataFrame(rows).dropna(subset=["water_year", "peak_flow_cfs"])
    df["water_year"] = df["water_year"].astype(int)
    return df


def fetch_sites(cache: Path, pause: float) -> pd.DataFrame:
    feats = _paged(
        cache,
        "sites",
        f"{API}/monitoring-locations/items",
        {
            "f": "json",
            "state_code": "30",
            "site_type_code": "ST",
            "limit": 10000,
            "properties": "monitoring_location_number,monitoring_location_name,"
            "hydrologic_unit_code,drainage_area,altitude",
        },
        pause,
    )
    rows = []
    for f in feats:
        p = f["properties"]
        geom = f.get("geometry") or {}
        lon, lat = (geom.get("coordinates") or [None, None])[:2]
        rows.append(
            {
                "site_no": str(p.get("monitoring_location_number")),
                "name": p.get("monitoring_location_name"),
                "huc": str(p.get("hydrologic_unit_code") or ""),
                "drainage_area": p.get("drainage_area"),
                "altitude": p.get("altitude"),
                "lat": lat,
                "lon": lon,
            }
        )
    return pd.DataFrame(rows).drop_duplicates("site_no")


# ----------------------------------------------------------------------------
# Screening and EMA
# ----------------------------------------------------------------------------


def screen(peaks: pd.DataFrame, min_years: int) -> Dict[str, Dict[str, Any]]:
    """Per site: the pre-regulation record, and why a site is out."""
    out: Dict[str, Dict[str, Any]] = {}
    for site, g in peaks.groupby("site_no"):
        g = g.sort_values("water_year").drop_duplicates("water_year")
        coded = g["qualification_code"].map(lambda c: bool(parse_codes(c) & SCREEN_CODES))
        cut = int(g.loc[coded, "water_year"].min()) if coded.any() else None
        keep = g[g["water_year"] < cut] if cut is not None else g
        n = int(len(keep))
        reason = "" if n >= min_years else f"< {min_years} pre-regulation peaks ({n})"
        out[str(site)] = {"record": keep, "regulated_from": cut, "n": n, "reason": reason}
    return out


def run_ema(cache: Path, screened: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    d = cache / "ema"
    d.mkdir(parents=True, exist_ok=True)
    out: Dict[str, Dict[str, Any]] = {}
    todo = [s for s, v in screened.items() if not v["reason"]]
    for i, site in enumerate(sorted(todo)):
        f = d / f"{site}.json"
        if f.exists():
            out[site] = json.loads(f.read_text(encoding="utf-8"))
            continue
        try:
            s = station_skew(site, screened[site]["record"])
            rec: Dict[str, Any] = dict(s.__dict__)
        except Exception as exc:  # noqa: BLE001 -- record and move on
            log.warning("%s: EMA failed: %s", site, exc)
            rec = {"site_no": site, "error": str(exc)}
        f.write_text(json.dumps(rec), encoding="utf-8")
        out[site] = rec
        log.info("EMA %d/%d %s", i + 1, len(todo), site)
    return out


def redundancy(df: pd.DataFrame) -> List[str]:
    """Sites to drop as redundant (the shorter-PRL member of each nested pair)."""
    da = df["drainage_area"].to_numpy(float)
    d = distance_miles(df["lat"], df["lon"])
    sd = d / np.sqrt(0.5 * (da[:, None] + da[None, :]))
    dar = np.maximum(da[:, None] / da[None, :], da[None, :] / da[:, None])
    drop: set = set()
    prl = df["prl"].to_numpy(float)
    sites = df["site_no"].tolist()
    order = np.argsort(-prl)
    for a in order:
        if sites[a] in drop:
            continue
        for b in range(len(sites)):
            if b != a and sites[b] not in drop and sd[a, b] <= 0.5 and dar[a, b] <= 5:
                if prl[b] <= prl[a]:
                    drop.add(sites[b])
    return sorted(drop)


def correlation_model(
    screened: Dict[str, Dict[str, Any]], df: pd.DataFrame, min_concurrent: int
) -> CorrelationDistanceModel:
    recs = {
        s: screened[s]["record"].set_index("water_year")["peak_flow_cfs"] for s in df["site_no"]
    }
    dist = distance_miles(df["lat"], df["lon"])
    ds, rs = [], []
    sites = df["site_no"].tolist()
    for i in range(len(sites)):
        for j in range(i + 1, len(sites)):
            a, b = recs[sites[i]], recs[sites[j]]
            common = a.index.intersection(b.index)
            x, y = a.loc[common], b.loc[common]
            ok = (x > 0) & (y > 0)
            if ok.sum() >= min_concurrent:
                ds.append(dist[i, j])
                rs.append(np.corrcoef(np.log10(x[ok]), np.log10(y[ok]))[0, 1])
    log.info("correlation model: %d pairs with >= %d concurrent years", len(ds), min_concurrent)
    return CorrelationDistanceModel.fit(ds, rs, min_concurrent=min_concurrent)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------


def _model_dict(m: Any) -> Dict[str, Any]:
    return {k: (round(v, 5) if isinstance(v, float) else v) for k, v in m.summary().items()}


def main(argv: Sequence[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cache", type=Path, default=ROOT / ".cache" / "mt_skew")
    ap.add_argument("--pause", type=float, default=4.0, help="seconds between OGC requests")
    ap.add_argument("--min-years", type=int, default=25)
    ap.add_argument("--min-concurrent", type=int, default=50)
    ap.add_argument("--kappa", type=float, default=DEFAULT_KAPPA)
    ap.add_argument("--last-wy", type=int, default=2025, help="last water year used")
    args = ap.parse_args(list(argv))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    logging.getLogger("flowfreq").setLevel(logging.ERROR)
    args.cache.mkdir(parents=True, exist_ok=True)

    peaks = fetch_peaks(args.cache, args.pause)
    peaks = peaks[peaks["water_year"] <= args.last_wy]
    sites = fetch_sites(args.cache, args.pause).set_index("site_no")
    log.info(
        "%d peaks at %d sites; %d stream sites", len(peaks), peaks.site_no.nunique(), len(sites)
    )
    screened = screen(peaks[peaks["site_no"].isin(sites.index)], args.min_years)
    ema = run_ema(args.cache, screened)

    rows = []
    for site, v in sorted(screened.items()):
        meta = sites.loc[site]
        row: Dict[str, Any] = {
            "site_no": site,
            "name": meta["name"],
            "huc": meta["huc"],
            "lat": meta["lat"],
            "lon": meta["lon"],
            "drainage_area": meta["drainage_area"],
            "altitude": meta["altitude"],
            "regulated_from_wy": v["regulated_from"],
            "n_pre_regulation": v["n"],
        }
        e = ema.get(site, {})
        for k in ("skew", "skew_mse", "prl", "n_systematic", "n_peaks", "n_low_outliers",
                  "yb", "ye", "ph", "at_site_option"):  # fmt: skip
            row[k] = e.get(k)
        reason = v["reason"] or e.get("error", "")
        if not reason and not (pd.notna(meta["drainage_area"]) and meta["drainage_area"] > 0):
            reason = "no drainage area"
        if not reason and not str(meta["huc"])[:2] in ("09", "10", "17"):
            reason = f"HUC {meta['huc']} outside MT basins"
        row["excluded"] = reason
        rows.append(row)
    df = pd.DataFrame(rows)
    ok = df[df["excluded"] == ""].reset_index(drop=True)
    for s in redundancy(ok):
        df.loc[df["site_no"] == s, "excluded"] = "redundant (nested, DAR <= 5)"
    used = df[df["excluded"] == ""].reset_index(drop=True)
    used["east_of_divide"] = (used["huc"].str[:2] != "17").astype(int)
    log.info("%d sites used", len(used))

    INPUTS.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(INPUTS, index=False, lineterminator="\n")

    cmodel = correlation_model(screened, used, args.min_concurrent)
    corr = skew_correlation_matrix(
        used["lat"],
        used["lon"],
        used["prl"],
        used["yb"],
        used["ye"],
        used["ph"],
        cmodel,
        args.kappa,
    )
    g, p = used["skew"].to_numpy(float), used["prl"].to_numpy(float)
    const = fit_regional_skew(g, p, correlation=corr, kappa=args.kappa)
    results: Dict[str, Any] = {
        "version": VERSION,
        "status": "provisional -- NOT a USGS study; do not use as a default",
        "last_water_year": args.last_wy,
        "n_sites": len(used),
        "screen": {"min_years": args.min_years, "codes_removed": sorted(SCREEN_CODES)},
        "correlation_model": cmodel.__dict__,
        "kappa": args.kappa,
        "models": {"CONSTANT": _model_dict(const)},
    }
    s2_0 = const.sigma2
    cands = {
        "EAST": ["east_of_divide"],
        "LOGDA": ["log_da"],
        "ALT": ["alt_kft"],
        "LAT": ["lat"],
    }
    used["log_da"] = np.log10(used["drainage_area"].astype(float))
    used["alt_kft"] = used["altitude"].astype(float) / 1000.0
    for name, cols in cands.items():
        sub = used.dropna(subset=cols)
        if len(sub) < len(used):
            continue
        m = fit_regional_skew(
            g, p, used[cols].to_numpy(float), cols, correlation=corr, kappa=args.kappa,
            _constant_sigma2=s2_0,
        )  # fmt: skip
        results["models"][name] = _model_dict(m)
    # Sensitivity: east of the Divide only (Missouri + Hudson Bay), and PRL >= 35.
    for label, mask in (
        ("CONSTANT_EAST_ONLY", used["east_of_divide"] == 1),
        ("CONSTANT_PRL35", used["prl"] >= 35),
        # Large basins carry regulation the peak codes do not flag (e.g. the
        # Bighorn above Yellowstone River near Sidney).
        ("CONSTANT_DA_LE_3000", used["drainage_area"].astype(float) <= 3000),
        # SIR 2016-5083's screen: records extending through at least 1980.
        ("CONSTANT_THROUGH_1980", used["ye"].astype(float) >= 1980),
    ):
        idx = np.flatnonzero(mask.to_numpy())
        m = fit_regional_skew(g[idx], p[idx], correlation=corr[np.ix_(idx, idx)], kappa=args.kappa)
        results["models"][label] = _model_dict(m)
    for k in (2.8, 3.3):
        c = skew_correlation_matrix(
            used["lat"], used["lon"], used["prl"], used["yb"], used["ye"], used["ph"], cmodel, k
        )
        results["models"][f"CONSTANT_kappa{k}"] = _model_dict(
            fit_regional_skew(g, p, correlation=c, kappa=k)
        )
    top = np.argsort(-const.influence)[:15]
    results["high_influence"] = [
        {"site_no": used.loc[i, "site_no"], "skew": float(g[i]), "prl": float(p[i]),
         "residual": round(float(const.residuals[i]), 3)}
        for i in top
    ]  # fmt: skip
    results["exclusions"] = (
        df["excluded"].replace("", "used").str.split(" \\(").str[0].value_counts().to_dict()
    )
    RESULTS.write_text(json.dumps(results, indent=1, default=float) + "\n", encoding="utf-8")
    print(json.dumps(results["models"], indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
