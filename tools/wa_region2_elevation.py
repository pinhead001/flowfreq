"""Why the SIR 2016-5118 Peak Region 2 equation under-predicts the Methow (methow item 5).

The Region 2 equation (Table 6: ``log Q = b0 + b1 log A + b2 P + b3 CANOPY``) is 1.4-2.5x low
at the five Methow gages. This tool refits it with ``flowfreq.regression.gls`` on flowfreq EMA
quantiles through WY2025, adding candidate basin characteristics, and reports for each model
and AEP: the model error variance (MEV), the leave-one-out mean |dlog10 Q|, the t statistic of
each added term, and the five Methow at-site / regression ratios with each target gage left
out of the fit that predicts it. It also refits without any gage of the Methow HUC8
(17020008), an out-of-basin test, and summarizes base-model residuals by distance from the
Methow. Findings: ``docs/WA_REGION2_METHOW.md``.

GLS setup as ``tools/wa_region4_refit.py`` (WREG LP3 sampling variance with the station SD
smoothed on the design, the weighted skew, skew MSE 0.18; the Pacific Northwest
concurrent-peak correlation between basin centroids), with concurrent years counted from the
peak records themselves.

Usage::

    python tools/wa_region2_elevation.py --build --cache <dir>   # inputs CSV (network)
    python tools/wa_region2_elevation.py                         # results JSON (offline)

``--build`` fetches each Region 2 gage's StreamStats GageStats characteristics and computes
its EMA quantiles through WY2025 (Water Data API peaks, PNW skew -0.07 with MSE 0.18, peak
codes applied), caching responses under ``--cache``.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tools"))

import wa_region4_refit as wr  # noqa: E402

from flowfreq.regression.gls import fit_gls, lp3_quantile_sampling_variance  # noqa: E402

DATA = REPO_ROOT / "data" / "regression_study"
INPUTS = DATA / "wa_2016_5118_r2_inputs.csv"
OUT = DATA / "wa_2016_5118_r2_elevation.json"
LAST_WY = 2025
AEPS = wr.AEPS
SCORE_AEPS = (0.5, 0.1, 0.01)
METHOW_HUC8 = "17020008"
#: The five Methow gages methow-flow transfers from; SIR 2016-5118 calibrated on two.
METHOW = {
    "12447383": "Methow ab Goat Cr",
    "12448998": "Twisp R",
    "12449500": "Methow at Twisp",
    "12448500": "Methow at Winthrop",
    "12448000": "Chewuch at Winthrop",
}
#: GageStats characteristics carried into the inputs (SIR 2016-5118 citations).
CHARS = ("ELEV", "MINBELEV", "ELEVMAX", "RELIEF", "JANMINTMP", "FOREST", "I24H2Y")
MODELS: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("base [1, log A, P, CAN]", ()),
    ("+ELEV", ("ELEV",)),
    ("+MINBELEV", ("MINBELEV",)),
    ("+JANMINTMP", ("JANMINTMP",)),
    ("+I24H2Y", ("I24H2Y",)),
    ("+FOREST", ("FOREST",)),
    ("+METHOW_HUC8", ("METHOW_HUC8",)),
    ("+ELEV+METHOW_HUC8", ("ELEV", "METHOW_HUC8")),
)


# ----------------------------------------------------------------------------- build


def _gagestats(site: str, cache: Path) -> Dict[str, Any]:
    import requests

    path = cache / f"gs_{site}.json"
    if not path.exists():
        for _ in range(3):
            try:
                r = requests.get(
                    f"https://streamstats.usgs.gov/gagestatsservices/stations/{site}", timeout=90
                )
                if r.status_code == 200:
                    path.write_text(r.text, encoding="utf-8")
                    break
            except requests.RequestException:
                time.sleep(3)
    if not path.exists():
        return {}
    d = json.loads(path.read_text(encoding="utf-8"))
    out: Dict[str, Any] = {}
    preferred: Dict[str, bool] = {}
    for c in d.get("characteristics", []):
        code = (c.get("variableType") or {}).get("code")
        citation = c.get("citation") if isinstance(c.get("citation"), dict) else {}
        is_2016 = str(citation.get("title") or "").startswith("2016")
        if code in CHARS and (code not in out or (is_2016 and not preferred.get(code))):
            out[code] = c.get("value")
            preferred[code] = is_2016
    return out


def _ema(site: str, cache: Path) -> Dict[str, Any]:
    path = cache / f"ema_{site}_{LAST_WY}.json"
    peaks_path = cache / f"peaks_{site}.csv"
    if path.exists():
        rec = json.loads(path.read_text())
        if "error" in rec or not peaks_path.exists():
            return rec if "error" in rec else {"error": "cached EMA without its peaks file"}
    else:
        from flowfreq.usgs import USGSgage
        from flowfreq.workflow import run_ffa

        try:
            pk = USGSgage(site).download_peak_flow()
        except ValueError as exc:
            rec = {"error": f"no peaks: {exc}"[:160]}
            path.write_text(json.dumps(rec))
            return rec
        pk.to_csv(peaks_path, index=False)
        pk = pk[pk["water_year"] <= LAST_WY].dropna(subset=["peak_flow_cfs"])
        r = run_ffa(
            pk["peak_flow_cfs"].to_numpy(),
            pk["water_year"].to_numpy(),
            regional_skew=-0.07,
            regional_skew_se=0.18**0.5,
            peak_codes=pk["qualification_code"].fillna("").astype(str).tolist(),
            site_no=site,
            allow_regulated=True,
        )
        if r.get("error"):
            rec = {"error": str(r["error"])}
        else:
            qd = r["quantile_df"]
            aep, flow = qd["AEP (%)"].to_numpy()[::-1], qd["Flow (cfs)"].to_numpy()[::-1]
            rec = {f"q_{a}": float(np.interp(a, aep, flow)) for a in AEPS}
            rec.update(sd=float(r["parameters"]["std_log"]),
                       skew_w=float(r["parameters"]["skew_weighted"]), n=int(len(pk)))  # fmt: skip
        path.write_text(json.dumps(rec))
    if "error" in rec:
        return rec
    pk = pd.read_csv(peaks_path)
    pk = pk[(pk["water_year"] <= LAST_WY) & pk["peak_flow_cfs"].notna()]
    rec["peak_years"] = ";".join(str(y) for y in sorted(set(pk["water_year"].astype(int))))
    return rec


def build(cache: Path) -> pd.DataFrame:
    """Region 2 calibration gages and the Methow targets, with characteristics and EMA."""
    from flowfreq.catalog import load_catalog

    cache.mkdir(parents=True, exist_ok=True)
    sites = pd.read_csv(wr.SITES, dtype={"site_no": str})
    r2 = sites[sites["region"] == 2]
    r2 = r2[(r2["used"] == True) | r2["site_no"].isin(METHOW)]  # noqa: E712
    huc8 = load_catalog().set_index("site_no")["huc8"].astype(str)
    rows = []
    for _, s in r2.iterrows():
        rec = {"site_no": s["site_no"], "used": bool(s["used"]),
               "huc8": huc8.get(s["site_no"], "")}  # fmt: skip
        rec.update(_gagestats(s["site_no"], cache))
        ema = _ema(s["site_no"], cache)
        if "error" in ema:
            logging.warning("%s: no EMA (%s)", s["site_no"], ema["error"])
            continue
        rec.update(ema)
        rows.append(rec)
    out = pd.DataFrame(rows)
    out.to_csv(INPUTS, index=False)
    print(f"wrote {INPUTS.relative_to(REPO_ROOT)}: {len(out)} gages")
    return out


# --------------------------------------------------------------------------- analysis


def load_frame() -> pd.DataFrame:
    sites = pd.read_csv(wr.SITES, dtype={"site_no": str})
    keep = ["site_no", "lat_cent", "lon_cent", "DRNAREA", "PRECIP", "CANOPY_PCT"]
    inputs = pd.read_csv(INPUTS, dtype={"site_no": str, "huc8": str, "peak_years": str})
    d = inputs.merge(sites[keep], on="site_no", how="left")
    d["METHOW_HUC8"] = (d["huc8"] == METHOW_HUC8).astype(float)
    return d


def design(df: pd.DataFrame, extra: Sequence[str]) -> np.ndarray:
    cols = [np.ones(len(df)), np.log10(df["DRNAREA"].to_numpy(float)),
            df["PRECIP"].to_numpy(float), df["CANOPY_PCT"].to_numpy(float)]  # fmt: skip
    return np.column_stack(cols + [df[t].to_numpy(float) for t in extra])


def concurrent_years(df: pd.DataFrame) -> np.ndarray:
    years = [set(int(y) for y in str(s).split(";") if y) for s in df["peak_years"]]
    return np.array([[len(a & b) for b in years] for a in years], dtype=float)


def sigma(df: pd.DataFrame, X: np.ndarray, aep: float) -> np.ndarray:
    coef, *_ = np.linalg.lstsq(X, df["sd"].to_numpy(float), rcond=None)
    n = df["n"].to_numpy(float)
    var = lp3_quantile_sampling_variance(
        X @ coef, df["skew_w"].to_numpy(float), n, aep, wr.PNW_SKEW_MSE
    )
    rho = wr.PNW_CORRELATION.rho(
        wr.distance_miles(df["lat_cent"].to_numpy(float), df["lon_cent"].to_numpy(float))
    )
    s = rho * np.sqrt(np.outer(var, var)) * concurrent_years(df) / np.sqrt(np.outer(n, n))
    np.fill_diagonal(s, var)
    return np.asarray(s, dtype=float)


def evaluate_model(
    cal: pd.DataFrame, targets: pd.DataFrame, extra: Sequence[str], aep: float
) -> Dict[str, Any]:
    X = design(cal, extra)
    y = np.log10(cal[f"q_{aep}"].to_numpy(float))
    S = sigma(cal, X, aep)
    res = fit_gls(y, X, S)
    se = np.sqrt(np.diag(res.cov_beta))
    loo = []
    for i in range(len(cal)):
        k = np.r_[0:i, i + 1 : len(cal)]
        loo.append(abs(float(X[i] @ fit_gls(y[k], X[k], S[np.ix_(k, k)]).beta) - y[i]))
    Xt = design(targets, extra)
    ratios = {}
    for j, site in enumerate(targets["site_no"]):
        k = np.flatnonzero(cal["site_no"].to_numpy() != site)
        b = fit_gls(y[k], X[k], S[np.ix_(k, k)]).beta
        ratios[site] = float(targets[f"q_{aep}"].iloc[j] / 10 ** float(Xt[j] @ b))
    return {
        "m": len(cal),
        "beta": [float(b) for b in res.beta],
        "mev": float(res.model_error_variance),
        "loo_mean_abs_dlog": float(np.mean(loo)),
        "t": {t: float(res.beta[4 + i] / se[4 + i]) for i, t in enumerate(extra)},
        "methow_ratio": ratios,
        "methow_geomean": float(np.exp(np.mean(np.log(list(ratios.values()))))),
    }


def analyse() -> Dict[str, Any]:
    d = load_frame()
    cal_all = d[d["used"]].reset_index(drop=True)
    targets = d[d["site_no"].isin(METHOW)].set_index("site_no").loc[list(METHOW)].reset_index()
    out: Dict[str, Any] = {"last_water_year": LAST_WY, "methow": METHOW, "models": {}}
    for label, extra in MODELS:
        cal = cal_all.dropna(subset=list(extra)).reset_index(drop=True)
        out["models"][label] = {
            str(aep): evaluate_model(cal, targets, extra, aep) for aep in SCORE_AEPS
        }
    outside = cal_all[(cal_all["METHOW_HUC8"] == 0) & cal_all["ELEV"].notna()].reset_index(
        drop=True
    )
    out["outside_methow_huc8"] = {
        label: {str(aep): evaluate_model(outside, targets, extra, aep) for aep in SCORE_AEPS}
        for label, extra in MODELS[:2]
    }
    cal = cal_all
    dist = wr.distance_miles(np.r_[48.5, cal["lat_cent"]], np.r_[-120.3, cal["lon_cent"]])[0, 1:]
    bands: Dict[str, Any] = {}
    for aep in SCORE_AEPS:
        X = design(cal, ())
        y = np.log10(cal[f"q_{aep}"].to_numpy(float))
        resid = y - X @ fit_gls(y, X, sigma(cal, X, aep)).beta
        bands[str(aep)] = [
            {
                "miles": [lo, hi],
                "n": int(sel.sum()),
                "mean": float(resid[sel].mean()),
                "share_positive": float(np.mean(resid[sel] > 0)),
            }  # fmt: skip
            for lo, hi in ((0, 40), (40, 80), (80, 150), (150, 1000))
            for sel in [(dist >= lo) & (dist < hi)]
            if sel.any()
        ]
    out["base_residuals_by_distance"] = bands
    return out


def report(out: Dict[str, Any]) -> None:
    for aep in map(str, SCORE_AEPS):
        print(f"\nAEP {aep}")
        print(f"  {'model':<22} {'m':>3} {'MEV':>7} {'LOO':>7} {'t':>12}  Methow ratios (gm)")
        for label, by_aep in out["models"].items():
            r = by_aep[aep]
            t = ", ".join(f"{v:.2f}" for v in r["t"].values())
            ratios = " ".join(f"{v:.2f}" for v in r["methow_ratio"].values())
            print(f"  {label:<22} {r['m']:>3} {r['mev']:>7.4f} {r['loo_mean_abs_dlog']:>7.4f} "
                  f"{t:>12}  {ratios} ({r['methow_geomean']:.2f})")  # fmt: skip
    print("\nOutside the Methow HUC8:")
    for label, by_aep in out["outside_methow_huc8"].items():
        print("  " + label + ": " + "; ".join(
            f"AEP {a} gm {r['methow_geomean']:.2f}" + (f" t {list(r['t'].values())[0]:.2f}" if r["t"] else "")
            for a, r in by_aep.items()))  # fmt: skip


def main(argv: Sequence[str] = ()) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--build", action="store_true", help="rebuild the inputs CSV (network)")
    ap.add_argument("--cache", type=Path, default=REPO_ROOT / ".cache" / "wa_region2")
    args = ap.parse_args(list(argv))
    logging.disable(logging.WARNING)
    if args.build:
        build(args.cache)
    out = analyse()
    OUT.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    report(out)
    print(f"\nwrote {OUT.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
