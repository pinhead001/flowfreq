"""Snake River Plain regional skew: feasibility of a B-WLS/B-GLS study (Idaho).

Northwest studies plan, task 8. NOT a USGS study. SIR 2016-5083 (p. 52) says the
Pacific Northwest regional skew is not valid on the Snake River Plain (its
"region 0") because of regulation, groundwater-surface water interaction and
infiltration, and gives no alternative. This asks whether enough unregulated
gages with adequate records exist there for a study of its own, following the
same procedure as the Montana v2 study (``build_montana_skew_study_v2.py``):

1. Idaho peaks and stream sites from the Water Data OGC API (the v2 cache is
   reused; bulk queries, a handful of requests).
2. Region-0 gages: the gage's 12-digit HUC is an ``excluded`` Snake River Plain
   unit in ``flowfreq/data/regional_skew_hucs.csv``; for a ``boundary`` unit the
   gage location decides (StreamStats region polygons, as
   :func:`flowfreq.regional_skew.regional_skew_at` does).
3. Screening as Montana v2: record before the first peak coded 5, 6 or C; at
   least ``--min-years`` peaks in it; GAGES-II ``regulated``/``urban`` gages out
   unless the record was cut at a coded regulation start; redundant nested pairs
   (gage locations) dropped.
4. **Decision rule** (fixed before looking at the result): feasible only with at
   least ``MIN_SITES`` = 20 screened gages. Then a CONSTANT B-WLS/B-GLS model is
   fitted and compared with the PNW value by t = (beta - (-0.07)) / SD(beta).

Writes ``data/skew_study/srp_inputs.csv`` and ``srp_results.json``.

Usage::

    python tools/snake_river_plain_skew.py --cache .cache/mt_skew_v2
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
import build_montana_skew_study_v2 as v2  # noqa: E402

from flowfreq.regional_skew import _on_snake_river_plain, load_huc_table  # noqa: E402
from flowfreq.regulation import load_screen  # noqa: E402
from flowfreq.skew_study import CorrelationDistanceModel  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "data" / "skew_study" / "srp_inputs.csv"
RESULTS = ROOT / "data" / "skew_study" / "srp_results.json"
MIN_SITES = 20
PNW_SKEW = -0.07
#: SIR 2016-5083 eq. B11, used when too few region-0 pairs share 50 years.
PNW_CORRELATION = CorrelationDistanceModel(a=0.21, b=-0.17, c=-0.0058)

log = logging.getLogger("srp_skew")


def idaho_sites(cache: Path, pause: float) -> pd.DataFrame:
    """Every Idaho stream site (the v2 sites cache holds the unfiltered pages)."""
    feats = v1._paged(cache, "sites_ID", f"{v1.API}/monitoring-locations/items", {}, pause)
    rows = []
    for f in feats:
        p = f["properties"]
        lon, lat = ((f.get("geometry") or {}).get("coordinates") or [None, None])[:2]
        rows.append(
            {
                "site_no": str(p.get("monitoring_location_number")),
                "name": p.get("monitoring_location_name"),
                "huc": str(p.get("hydrologic_unit_code") or ""),
                "drainage_area": p.get("drainage_area"),
                "lat": lat,
                "lon": lon,
            }
        )
    return pd.DataFrame(rows).drop_duplicates("site_no").set_index("site_no")


def in_region_0(sites: pd.DataFrame) -> pd.Series:
    """'excluded', 'boundary-in', 'boundary-out' or '' per site, from its HU12."""
    hucs = load_huc_table()
    srp = hucs[hucs["skew_region"] == "Snake River Plain"].set_index("huc")["rule"]
    out = {}
    for site, r in sites.iterrows():
        rule = srp.get(str(r["huc"])[:12], "")
        if rule == "boundary":
            on = _on_snake_river_plain(float(r["lat"]), float(r["lon"]))
            rule = "boundary-in" if on else "boundary-out"
        out[site] = rule
    return pd.Series(out)


def main(argv: Sequence[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cache", type=Path, default=ROOT / ".cache" / "mt_skew_v2")
    ap.add_argument("--pause", type=float, default=4.0)
    ap.add_argument("--min-years", type=int, default=25)
    ap.add_argument("--last-wy", type=int, default=2025)
    args = ap.parse_args(list(argv))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    logging.getLogger("flowfreq").setLevel(logging.ERROR)

    v2.fetch_state(args.cache, "ID", args.pause)  # fills the cache if a v2 run has not
    sites = idaho_sites(args.cache, args.pause)
    peaks = _all_id_peaks(args)
    peaks = peaks[peaks["water_year"] <= args.last_wy]
    with_peaks = sites[sites.index.isin(peaks["site_no"].unique())]
    loc = in_region_0(with_peaks)
    plain = with_peaks[loc.reindex(with_peaks.index).isin(["excluded", "boundary-in"])]
    log.info("%d Idaho peak sites; %d in region 0", len(with_peaks), len(plain))

    screened = v1.screen(peaks[peaks["site_no"].isin(plain.index)], args.min_years)
    reg = load_screen()
    ema = v1.run_ema(args.cache / "srp", screened)
    rows = []
    for site, v in sorted(screened.items()):
        meta = plain.loc[site]
        e = ema.get(site, {})
        row: Dict[str, Any] = {
            "site_no": site,
            "name": meta["name"],
            "huc": meta["huc"],
            "lat": meta["lat"],
            "lon": meta["lon"],
            "drainage_area": meta["drainage_area"],
            "region0": loc[site],
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
    if len(ok) > 1:
        for s in v1.redundancy(ok):
            df.loc[df["site_no"] == s, "excluded"] = "redundant (nested, DAR <= 5)"
    used = df[df["excluded"] == ""].reset_index(drop=True)
    INPUTS.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(INPUTS, index=False, lineterminator="\n")

    results: Dict[str, Any] = {
        "status": "provisional analysis -- NOT a USGS study",
        "min_sites": MIN_SITES,
        "region0_peak_sites": int(len(plain)),
        "screened_in": int(len(used)),
        "exclusions": df["excluded"].replace("", "used").str.split(" \\(").str[0]
        .value_counts().to_dict(),
    }  # fmt: skip
    if len(used) < MIN_SITES:
        results["decision"] = "not feasible"
    else:
        try:
            cmodel = v1.correlation_model(screened, used, 50)
        except ValueError:
            cmodel = PNW_CORRELATION
        corr = v1.skew_correlation_matrix(
            used["lat"], used["lon"], used["prl"], used["yb"], used["ye"], used["ph"],
            cmodel, v1.DEFAULT_KAPPA,
        )  # fmt: skip
        m = v1.fit_regional_skew(
            used["skew"].to_numpy(float), used["prl"].to_numpy(float), correlation=corr
        )
        t = (float(m.beta[0]) - PNW_SKEW) / float(m.beta_sd[0])
        results.update(
            decision="feasible",
            constant=v1._model_dict(m),
            t_vs_pnw=round(t, 2),
            correlation_model=cmodel.__dict__,
        )
    RESULTS.write_text(json.dumps(results, indent=1, default=float) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=1, default=float))
    return 0


def _all_id_peaks(args: argparse.Namespace) -> pd.DataFrame:
    """Idaho peaks for every site (fetch_state filters sites, not the cached pages)."""
    feats = v1._paged(args.cache, "peaks_ID", f"{v1.API}/peaks/items", {}, args.pause)
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
    out = pd.DataFrame(rows).dropna(subset=["water_year", "peak_flow_cfs"])
    out["water_year"] = out["water_year"].astype(int)
    return out


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
