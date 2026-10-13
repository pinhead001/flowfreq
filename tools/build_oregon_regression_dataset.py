"""Oregon peak-flow regression update, task 2: the station dataset.

NOT a USGS study. Oregon regression update plan
(``docs/superpowers/plans/2026-10-12-oregon-regression-update.md``), task 2.
Writes ``data/regression_study/oregon_stations.csv``: one row per candidate
gage, with its screening outcome, its skew source and, for the gages kept, the
flowfreq EMA quantiles at the eight AEPs the Pacific Northwest reports use.

Gages: every Oregon stream gage with annual peaks, plus Washington, Idaho,
Nevada and California gages within ``--border-miles`` (25) of the Oregon
boundary (U.S. Census TIGERweb). Peaks and sites come from bulk OGC queries
(:mod:`flowfreq.waterdata`, with ``USGS_API_KEY`` when set), cached and
resumable.

Screens, as the Montana v2 skew study
(``tools/build_montana_skew_study_v2.py``):

- the record before the first peak coded for regulation or alteration
  (codes 5, 6, C), and at least ``--min-years`` (10) such peaks, the threshold
  SIR 2016-5118 and SIR 2016-5083 use for regression gages;
- the packaged GAGES-II screen: a gage classed ``regulated`` or ``urban`` is
  out unless its record was already cut at a coded regulation peak;
- peaks through ``--last-wy`` (2025);
- at least ``--min-years`` nonzero peaks after the fit: a record of mostly
  zero-flow years leaves EMA two or three values to fit and an SD near 10 log
  units;
- drainage area at most ``MAX_AREA`` (20,000 mi2): the Columbia mainstem drains
  mostly outside the study area.

Nested gages (DAR <= 5 within half a basin length) are flagged in
``redundant_with``, not removed; GLS carries their correlated errors and task 5
decides.

EMA with MGBT (:func:`flowfreq.workflow.run_ffa`), the skew weighted with the
regional skew :func:`flowfreq.regional_skew.regional_skew_at_huc` gives for the
gage's own HUC12 and state, provisional rows allowed (the Oregon closed-basin
extension, HUC4 1712). Where no study covers the gage (the Klamath, 1801;
NV/CA closed basins) the fit is station-only and ``skew_source`` says so.

``validate`` refits the gages this dataset shares with SIR 2016-5083 Table A2
(Idaho; Wood and others, 2016) on that report's terms -- peaks through WY2013,
the Pacific Northwest skew -0.07 (MSE 0.18) -- and writes the ratios to
``oregon_sir2016_5083_check.csv``. Table A3's perception thresholds are not
applied; ``table_a3`` marks the gages that have them.

Usage::

    python tools/build_oregon_regression_dataset.py --cache <dir> inputs
    python tools/build_oregon_regression_dataset.py --cache <dir> validate
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
from _xlsx import read_sheet  # noqa: E402

from flowfreq.regional_skew import RegionalSkewUnavailable, regional_skew_at_huc  # noqa: E402
from flowfreq.regulation import load_screen  # noqa: E402
from flowfreq.skew_study import distance_miles  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "regression_study" / "oregon_stations.csv"
CHECK = ROOT / "data" / "regression_study" / "oregon_sir2016_5083_check.csv"
AEPS = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)
STATE_FIPS = {"41": "OR", "53": "WA", "16": "ID", "32": "NV", "06": "CA"}
#: Oregon's extent padded by ~35 miles, for the bulk site and peak queries.
BBOX = (-125.2, 41.4, -115.9, 46.8)
TIGER = (
    "https://tigerweb.geo.census.gov/arcgis/rest/services/TIGERweb/State_County/"
    "MapServer/0/query"
)
SIR_2016_5083_A = "https://pubs.usgs.gov/sir/2016/5083/sir20165083_appendixa.xlsx"
#: SIR 2016-5083's terms: peaks through WY2013, the PNW regional skew (app. B).
PNW_SKEW, PNW_MSE, IDAHO_LAST_WY = -0.07, 0.18, 2013
#: Larger basins are the Columbia mainstem, draining mostly outside the study area.
MAX_AREA = 20000.0

log = logging.getLogger("or_regression")


# ----------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------


def oregon_boundary(cache: Path) -> np.ndarray:
    """(lat, lon) points along the Oregon boundary, ~1 km apart (rings cached)."""
    f = cache / "oregon_boundary.json"
    if not f.exists():
        r = requests.get(
            TIGER,
            params={
                "where": "STUSAB='OR'",
                "outFields": "STUSAB",
                "returnGeometry": "true",
                "outSR": "4326",
                "maxAllowableOffset": "0.01",
                "f": "json",
            },
            timeout=300,
        )
        r.raise_for_status()
        rings = r.json()["features"][0]["geometry"]["rings"]
        f.write_text(json.dumps(rings))
    rings = json.loads(f.read_text())
    # The generalized rings keep only corners -- the 42nd-parallel border with
    # NV and CA is a few vertices hundreds of miles apart -- so densify each
    # edge to <= 0.01 degree before measuring distance to vertices.
    pts = []
    for ring in rings:
        r = np.asarray(ring, dtype=float)
        for p0, p1 in zip(r[:-1], r[1:]):
            k = max(1, int(np.ceil(np.abs(p1 - p0).max() / 0.01)))
            pts.append(p0 + (p1 - p0) * np.arange(k)[:, None] / k)
    return np.asarray(np.concatenate(pts)[:, ::-1])


def miles_to_boundary(lat: np.ndarray, lon: np.ndarray, boundary: np.ndarray) -> np.ndarray:
    """Great-circle (haversine) miles from each point to the nearest boundary point."""
    p1, l1 = np.radians(boundary[:, 0])[None, :], np.radians(boundary[:, 1])[None, :]
    p0, l0 = np.radians(lat)[:, None], np.radians(lon)[:, None]
    h = np.sin((p1 - p0) / 2) ** 2 + np.cos(p0) * np.cos(p1) * np.sin((l1 - l0) / 2) ** 2
    return np.asarray(2 * 3958.8 * np.arcsin(np.sqrt(h)).min(axis=1))


def fetch_sites(cache: Path, pause: float) -> pd.DataFrame:
    """Stream gages in :data:`BBOX` in the five states, from one paged OGC query."""
    feats = v1._paged(
        cache,
        "sites",
        f"{v1.API}/monitoring-locations/items",
        {
            "f": "json",
            "bbox": ",".join(str(x) for x in BBOX),
            "site_type_code": "ST",
            "limit": 10000,
            "properties": "monitoring_location_number,monitoring_location_name,state_code,"
            "hydrologic_unit_code,drainage_area,contributing_drainage_area,altitude",
        },
        pause,
    )
    rows = []
    for f in feats:
        p = f["properties"]
        lon, lat = ((f.get("geometry") or {}).get("coordinates") or [None, None])[:2]
        rows.append(
            {
                "site_no": str(p.get("monitoring_location_number")),
                "state": STATE_FIPS.get(str(p.get("state_code") or "").zfill(2), ""),
                "name": p.get("monitoring_location_name"),
                "huc": str(p.get("hydrologic_unit_code") or ""),
                "drainage_area": p.get("drainage_area"),
                "contributing_area": p.get("contributing_drainage_area"),
                "altitude": p.get("altitude"),
                "lat": lat,
                "lon": lon,
            }
        )
    sites = pd.DataFrame(rows).drop_duplicates("site_no")
    return sites[sites["state"] != ""].reset_index(drop=True)


def fetch_peaks(cache: Path, pause: float) -> pd.DataFrame:
    """Annual peaks in :data:`BBOX`, from one paged OGC query."""
    feats = v1._paged(
        cache,
        "peaks",
        f"{v1.API}/peaks/items",
        {
            "f": "json",
            "bbox": ",".join(str(x) for x in BBOX),
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
    return peaks


def study_sites(cache: Path, pause: float, border_miles: float) -> pd.DataFrame:
    """Oregon gages, and neighbouring-state gages within ``border_miles`` of Oregon."""
    sites = fetch_sites(cache, pause)
    boundary = oregon_boundary(cache)
    sites["border_miles"] = 0.0
    other = sites["state"] != "OR"
    ok = other & sites["lat"].notna()
    sites.loc[ok, "border_miles"] = miles_to_boundary(
        sites.loc[ok, "lat"].to_numpy(float), sites.loc[ok, "lon"].to_numpy(float), boundary
    )
    keep = (sites["state"] == "OR") | (ok & (sites["border_miles"] <= border_miles))
    return sites[keep].reset_index(drop=True)


# ----------------------------------------------------------------------------
# Skew and EMA
# ----------------------------------------------------------------------------


def skew_for(huc: str, state: str, lat: float, lon: float) -> Dict[str, Any]:
    """The regional skew a gage's EMA is weighted with, or why there is none."""
    try:
        rs = regional_skew_at_huc(huc, state, point=(lat, lon), allow_provisional=True)
    except RegionalSkewUnavailable as exc:
        return {"skew_source": f"at-site only ({str(exc)[:120]})"}
    src = rs.skew_region + (" [provisional]" if rs.provisional else "")
    return {"skew_source": src, "regional_skew": rs.skew, "regional_skew_mse": rs.skew_mse}


def fit(
    record: pd.DataFrame, site: str, skew: Optional[float], skew_mse: Optional[float]
) -> Dict[str, Any]:
    """EMA/MGBT on one record: moments, counts and the quantiles at :data:`AEPS`."""
    from flowfreq.workflow import run_ffa

    kw: Dict[str, Any] = (
        {"regional_skew": skew, "regional_skew_se": math.sqrt(float(skew_mse))}
        if skew is not None and skew_mse is not None
        else {"station_skew_only": True}
    )
    r = run_ffa(
        record["peak_flow_cfs"].to_numpy(float),
        record["water_year"].to_numpy(int),
        peak_codes=record["qualification_code"].tolist(),
        site_no=site,
        allow_regulated=True,
        **kw,
    )
    if r.get("error"):
        raise RuntimeError(r["error"])
    res = r["b17c"].results
    qd = r["quantile_df"]
    # Despite its header, "AEP (%)" holds fractions (0.5, 0.01, ...), largest first.
    aep = qd["AEP (%)"].to_numpy(float)[::-1]
    flow = qd["Flow (cfs)"].to_numpy(float)[::-1]
    out: Dict[str, Any] = {f"q_{a}": float(np.interp(a, aep, flow)) for a in AEPS}
    out.update(
        mean_log=float(res.mean_log),
        sd_log=float(res.std_log),
        skew_station=float(res.skew_station),
        skew_weighted=None if res.skew_weighted is None else float(res.skew_weighted),
        n_peaks=int(res.n_peaks),
        n_systematic=int(res.n_systematic),
        n_historical=int(res.n_historical),
        n_low_outliers=int(res.n_low_outliers),
        n_zeros=int(res.n_zeros),
        first_wy=int(record["water_year"].min()),
        last_wy=int(record["water_year"].max()),
        converged=bool(r.get("converged", True)),
    )
    return out


def run_fits(cache: Path, cand: pd.DataFrame, screened: Dict[str, Dict[str, Any]]) -> None:
    """EMA per candidate gage (cached as JSON); results land in ``cand`` in place."""
    d = cache / "ema"
    d.mkdir(parents=True, exist_ok=True)
    for i, r in enumerate(cand.itertuples()):
        f = d / f"{r.site_no}.json"
        if f.exists():
            rec = json.loads(f.read_text(encoding="utf-8"))
        else:
            rec = skew_for(r.huc, r.state, float(r.lat), float(r.lon))
            try:
                rec.update(
                    fit(
                        screened[r.site_no]["record"],
                        r.site_no,
                        rec.get("regional_skew"),
                        rec.get("regional_skew_mse"),
                    )
                )
            except Exception as exc:  # noqa: BLE001 -- recorded per gage
                rec["error"] = f"{type(exc).__name__}: {exc}"[:200]
            f.write_text(json.dumps(rec), encoding="utf-8")
            log.info("EMA %d/%d %s %s", i + 1, len(cand), r.site_no, rec.get("skew_source"))
        for k, v in rec.items():
            cand.loc[r.Index, k] = v


# ----------------------------------------------------------------------------
# Stages
# ----------------------------------------------------------------------------


def build_inputs(args: argparse.Namespace) -> pd.DataFrame:
    sites = study_sites(args.cache, args.pause, args.border_miles)
    peaks = fetch_peaks(args.cache, args.pause)
    peaks = peaks[peaks["site_no"].isin(sites["site_no"]) & (peaks["water_year"] <= args.last_wy)]
    log.info("%d peaks at %d gages", len(peaks), peaks["site_no"].nunique())
    screened = v1.screen(peaks, args.min_years)
    reg = load_screen()["site_class"]

    df = sites[sites["site_no"].isin(screened)].copy().reset_index(drop=True)
    df["regulated_from_wy"] = df["site_no"].map(lambda s: screened[s]["regulated_from"])
    df["n_pre_regulation"] = df["site_no"].map(lambda s: screened[s]["n"])
    df["reg_class"] = df["site_no"].map(lambda s: reg.get(s, "not in GAGES-II"))
    reasons = []
    for r in df.itertuples():
        reason = screened[r.site_no]["reason"]
        if not reason and not (pd.notna(r.drainage_area) and r.drainage_area > 0):
            reason = "no drainage area"
        if not reason and r.drainage_area > MAX_AREA:
            reason = f"drainage area > {MAX_AREA:,.0f} mi2 (mainstem)"
        if not reason and r.reg_class in ("regulated", "urban") and pd.isna(r.regulated_from_wy):
            reason = f"GAGES-II {r.reg_class}, no coded regulation start"
        reasons.append(reason)
    df["excluded"] = reasons

    cand = df[df["excluded"] == ""].copy()
    run_fits(args.cache, cand, screened)
    df = df.join(cand.drop(columns=df.columns.intersection(cand.columns)), how="left")
    err = (df["excluded"] == "") & df.get("error", pd.Series(index=df.index)).notna()
    df.loc[err, "excluded"] = "EMA failed"
    nonzero = df["n_peaks"] - df["n_zeros"]
    few = (df["excluded"] == "") & (nonzero < args.min_years)
    df.loc[few, "excluded"] = f"< {args.min_years} nonzero peaks"

    used = df["excluded"] == ""
    red = {s: "" for s in df["site_no"]}
    u = df[used].reset_index(drop=True)
    if len(u) > 1:
        da = u["drainage_area"].to_numpy(float)
        dist = distance_miles(u["lat"], u["lon"])
        sd = dist / np.sqrt(0.5 * (da[:, None] + da[None, :]))
        dar = np.maximum(da[:, None] / da[None, :], da[None, :] / da[:, None])
        np.fill_diagonal(sd, np.inf)
        for i, j in zip(*np.nonzero((sd <= 0.5) & (dar <= 5))):
            red[u.loc[i, "site_no"]] = ";".join(
                x for x in (red[u.loc[i, "site_no"]], u.loc[j, "site_no"]) if x
            )
    df["redundant_with"] = df["site_no"].map(red)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.sort_values(["state", "site_no"]).to_csv(
        OUT, index=False, lineterminator="\n", float_format="%.6g"
    )
    log.info(
        "wrote %s: %d gages, %d used (%s)",
        OUT.name,
        len(df),
        int(used.sum()),
        df.loc[used, "state"].value_counts().to_dict(),
    )
    return df


def table_a2(path: Path) -> Dict[str, Dict[str, float]]:
    """SIR 2016-5083 Table A2: site -> {AEP: quantile, cfs}."""
    rows = read_sheet(path, 1)
    head = next(r for r in rows if r and r[0].startswith("Map identification"))
    cols: Dict[float, int] = {}
    for j, h in enumerate(head):
        h = h.strip()
        if h.startswith("Q") and h[1:].replace(".", "").isdigit():
            cols[float(h[1:]) / 100.0] = j
    out: Dict[str, Dict[str, float]] = {}
    for r in rows:
        site = r[1].strip() if len(r) > 1 else ""
        if not site.isdigit():
            continue
        vals = {}
        for a in AEPS:
            if a in cols:
                s = r[cols[a]].replace(",", "").strip()
                vals[str(a)] = float(s) if s else float("nan")
        out[site] = vals
    return out


def validate(args: argparse.Namespace) -> pd.DataFrame:
    """flowfreq EMA vs SIR 2016-5083 Table A2 on the gages both datasets hold."""
    path = args.cache / "sir20165083_appendixa.xlsx"
    if not path.exists():
        r = requests.get(SIR_2016_5083_A, timeout=300)
        r.raise_for_status()
        path.write_bytes(r.content)
    a2 = table_a2(path)
    # Table A3: gages with perception-threshold settings (historic periods,
    # missing years), which this check does not apply.
    a3 = {r[1].strip() for r in read_sheet(path, 2) if len(r) > 1 and r[1].strip().isdigit()}
    df = pd.read_csv(OUT, dtype={"site_no": str, "huc": str})
    shared = sorted(set(df["site_no"]) & set(a2))
    peaks = fetch_peaks(args.cache, args.pause)
    rows: List[Dict[str, Any]] = []
    for site in shared:
        rec = peaks[(peaks["site_no"] == site) & (peaks["water_year"] <= IDAHO_LAST_WY)]
        rec = rec.sort_values("water_year").drop_duplicates("water_year")
        row: Dict[str, Any] = {
            "site_no": site,
            "state": df.loc[df["site_no"] == site, "state"].iloc[0],
            "n_wy2013": len(rec),
            "table_a3": site in a3,
        }
        try:
            got = fit(rec, site, PNW_SKEW, PNW_MSE)
        except Exception as exc:  # noqa: BLE001 -- recorded per gage
            row["error"] = f"{type(exc).__name__}: {exc}"[:200]
            rows.append(row)
            continue
        for a in AEPS:
            ref = a2[site].get(str(a), float("nan"))
            row[f"usgs_{a}"] = ref
            row[f"ratio_{a}"] = got[f"q_{a}"] / ref if ref and ref > 0 else float("nan")
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(CHECK, index=False, lineterminator="\n", float_format="%.6g")
    r = out[[f"ratio_{a}" for a in (0.5, 0.01, 0.002)]]
    log.info("wrote %s: %d shared gages; median ratios %s", CHECK.name, len(out),
             r.median().round(4).to_dict())  # fmt: skip
    return out


def main(argv: Sequence[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("stage", choices=["inputs", "validate"])
    ap.add_argument("--cache", type=Path, default=ROOT / ".cache" / "or_regression")
    ap.add_argument("--pause", type=float, default=1.0)
    ap.add_argument("--min-years", type=int, default=10)
    ap.add_argument("--last-wy", type=int, default=2025)
    ap.add_argument("--border-miles", type=float, default=25.0)
    args = ap.parse_args(list(argv))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    logging.getLogger("flowfreq").setLevel(logging.ERROR)
    args.cache.mkdir(parents=True, exist_ok=True)
    if args.stage == "inputs":
        build_inputs(args)
    else:
        validate(args)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
