"""Does the Pacific Northwest skew apply in Oregon's closed basins and the Klamath?

Northwest studies plan, task 9. NOT a USGS study. ``regional_skew_hucs.csv``
leaves HUC4 1712 (Oregon closed basins) and 1801 (Klamath) ``unresolved``:
SIR 2016-5118 (pp. 1, 23) defines the Pacific Northwest (PNW) study region
without them, yet 5 of the study's 290 gages are there (SIR 2016-5083 Table B1).

This tests the PNW value against the data. Every Oregon peak-flow gage in those
two basins is screened as in the Montana v2 study (pre-regulation record, at
least ``--min-years`` peaks, GAGES-II regulated/urban without a coded start
out, redundant nested pairs out); a CONSTANT B-WLS/B-GLS skew is fitted with
the PNW correlation-distance model (SIR 2016-5083 eq. B11) and compared with
G = -0.07 by t = (beta - (-0.07)) / SD(beta).

**Decision rule** (fixed in the plan before the data were looked at):
|t| < 2 -> "consistent"; otherwise "inconsistent". Also reported per basin.

Writes ``data/skew_study/pnw_edge_inputs.csv`` and ``pnw_edge_results.json``.

Usage::

    python tools/pnw_edge_skew_test.py --cache .cache/pnw_edge
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, Sequence

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_montana_skew_study as v1  # noqa: E402

from flowfreq.regulation import load_screen  # noqa: E402
from flowfreq.skew_study import CorrelationDistanceModel  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "data" / "skew_study" / "pnw_edge_inputs.csv"
RESULTS = ROOT / "data" / "skew_study" / "pnw_edge_results.json"
PNW_SKEW = -0.07
PNW_CORRELATION = CorrelationDistanceModel(a=0.21, b=-0.17, c=-0.0058)
BASINS = {"1712": "Oregon closed basins", "1801": "Klamath"}

log = logging.getLogger("pnw_edge")


def fetch_oregon(cache: Path, pause: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    feats = v1._paged(
        cache,
        "peaks_OR",
        f"{v1.API}/peaks/items",
        {
            "f": "json",
            "state_code": "41",
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
        "sites_OR",
        f"{v1.API}/monitoring-locations/items",
        {
            "f": "json",
            "state_code": "41",
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
                "name": p.get("monitoring_location_name"),
                "huc": str(p.get("hydrologic_unit_code") or ""),
                "drainage_area": p.get("drainage_area"),
                "lat": lat,
                "lon": lon,
            }
        )
    sites = pd.DataFrame(srows).drop_duplicates("site_no")
    sites = sites[sites["huc"].str[:4].isin(BASINS)].set_index("site_no")
    return peaks[peaks["site_no"].isin(sites.index)], sites


def fit_constant(used: pd.DataFrame) -> Dict[str, Any]:
    corr = v1.skew_correlation_matrix(
        used["lat"], used["lon"], used["prl"], used["yb"], used["ye"], used["ph"],
        PNW_CORRELATION, v1.DEFAULT_KAPPA,
    )  # fmt: skip
    m = v1.fit_regional_skew(
        used["skew"].to_numpy(float), used["prl"].to_numpy(float), correlation=corr
    )
    t = (float(m.beta[0]) - PNW_SKEW) / float(m.beta_sd[0])
    return {
        "n": int(len(used)),
        "constant": v1._model_dict(m),
        "t_vs_pnw": round(t, 2),
        "decision": "consistent" if abs(t) < 2 else "inconsistent",
    }


def main(argv: Sequence[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cache", type=Path, default=ROOT / ".cache" / "pnw_edge")
    ap.add_argument("--pause", type=float, default=2.0)
    ap.add_argument("--min-years", type=int, default=25)
    ap.add_argument("--last-wy", type=int, default=2025)
    args = ap.parse_args(list(argv))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    logging.getLogger("flowfreq").setLevel(logging.ERROR)
    args.cache.mkdir(parents=True, exist_ok=True)

    peaks, sites = fetch_oregon(args.cache, args.pause)
    peaks = peaks[peaks["water_year"] <= args.last_wy]
    screened = v1.screen(peaks, args.min_years)
    ema = v1.run_ema(args.cache, screened)
    reg = load_screen()
    rows = []
    for site, v in sorted(screened.items()):
        meta = sites.loc[site]
        e = ema.get(site, {})
        row: Dict[str, Any] = {
            "site_no": site,
            "name": meta["name"],
            "huc": meta["huc"],
            "basin": BASINS[str(meta["huc"])[:4]],
            "lat": meta["lat"],
            "lon": meta["lon"],
            "drainage_area": meta["drainage_area"],
            "n_pre_regulation": v["n"],
            "regulated_from_wy": v["regulated_from"],
            "reg_class": reg["site_class"].get(site, "not in GAGES-II"),
        }
        for k in ("skew", "skew_mse", "prl", "yb", "ye", "ph"):
            row[k] = e.get(k)
        reason = v["reason"] or e.get("error", "")
        if not reason and not (pd.notna(meta["drainage_area"]) and meta["drainage_area"] > 0):
            reason = "no drainage area"
        if (
            not reason
            and row["reg_class"] in ("regulated", "urban")
            and v["regulated_from"] is None
        ):
            reason = f"GAGES-II {row['reg_class']}, no coded regulation start"
        row["excluded"] = reason
        rows.append(row)
    df = pd.DataFrame(rows)
    ok = df[df["excluded"] == ""].reset_index(drop=True)
    for s in v1.redundancy(ok):
        df.loc[df["site_no"] == s, "excluded"] = "redundant (nested, DAR <= 5)"
    used = df[df["excluded"] == ""].reset_index(drop=True)
    INPUTS.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(INPUTS, index=False, lineterminator="\n")

    results: Dict[str, Any] = {
        "status": "provisional analysis -- NOT a USGS study",
        "pnw_skew": PNW_SKEW,
        "rule": "|t| < 2 -> consistent",
        "combined": fit_constant(used),
        "by_basin": {},
        "exclusions": df["excluded"].replace("", "used").str.split(" \\(").str[0]
        .value_counts().to_dict(),
        "mean_station_skew": round(float(np.mean(used["skew"])), 3),
    }  # fmt: skip
    for huc4, name in BASINS.items():
        sub = used[used["huc"].str[:4] == huc4].reset_index(drop=True)
        results["by_basin"][name] = fit_constant(sub) if len(sub) >= 5 else {"n": len(sub)}
    RESULTS.write_text(json.dumps(results, indent=1, default=float) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
