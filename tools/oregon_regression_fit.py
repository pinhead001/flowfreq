"""Oregon peak-flow regression update, tasks 4-5: regions, GLS fits, diagnostics.

NOT a USGS study. Oregon regression update plan
(``docs/superpowers/plans/2026-10-12-oregon-regression-update.md``), tasks 4 and 5.
Fits GLS regressions (:mod:`flowfreq.regression.gls`, the WREG procedure) of the
station log10 quantiles in ``oregon_stations.csv`` on the basin characteristics
in ``oregon_basin_characteristics.csv``.

Regions (task 4). NSS cannot place a basin in Oregon's current regions: every
western basin intersects all of Regions 1, 2A and 2B at 100 %, and OWRD's six
eastern regions have no geometry. So regions here are sets of hydrologic units,
assignable from a basin's outlet HUC alone. Start: **West** (HUC4 1708, 1709,
1710 and the north coastal Klamath HUC6 180101) and **East** (everything else).
Each fitted region's HUC6 residual clusters are then tested; a cluster whose
mean residual differs from zero by more than two standard errors, with enough
gages to stand alone, becomes its own region (``SPLITS``).

GLS specification, as the Washington refit (``tools/wa_region4_refit.py``):

- y: log10 station quantile at each of the eight AEPs;
- sampling variance: :func:`~flowfreq.regression.gls.lp3_quantile_sampling_variance`
  with the station SD smoothed by OLS on the region's design, the skew the
  quantile used, the skew's MSE (0.18 for the PNW rows, 0.302 -- the Bulletin
  17B map's -- for station-only fits), and n = peaks in the record;
- cross-correlation: the PNW concurrent-peak model (SIR 2016-5083 eq. B11)
  between basin centroids; concurrent years from first and last water years.

Variables (task 5): log10 DRNAREA is forced; others enter by forward selection
at AEP 0.01 while every added term has |t| >= 2 and every VIF < 5, and the
selected set is used at all eight AEPs. Leave-one-out cross-validation gives a
prediction SEP to set beside the AVP-based one and the current equations'.

Writes ``data/regression_study/oregon_regression_results.json`` and, with
``--store``, ``flowfreq/data/regression/provisional/OR.json``.

Usage::

    python tools/oregon_regression_fit.py [--store]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from flowfreq.regression.gls import GLSResult, fit_gls, lp3_quantile_sampling_variance
from flowfreq.skew_study import CorrelationDistanceModel, distance_miles

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "regression_study"
STATIONS = DATA / "oregon_stations.csv"
BASINS = DATA / "oregon_basin_characteristics.csv"
OUT = DATA / "oregon_regression_results.json"
PROVISIONAL = ROOT / "flowfreq" / "data" / "regression" / "provisional" / "OR.json"
OR_JSON = ROOT / "flowfreq" / "data" / "regression" / "OR.json"

AEPS = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)
PNW_CORRELATION = CorrelationDistanceModel(a=0.21, b=-0.17, c=-0.0058)
PNW_SKEW_MSE = 0.18
#: Skew MSE for a station-only fit: the Bulletin 17B skew map's.
AT_SITE_SKEW_MSE = 0.302
SELECT_AEP = 0.01

#: Explanatory variables: design-column name -> (CSV column, stored code, transform,
#: scale, units). Codes name their source so a StreamStats value of a same-named
#: characteristic (computed from other grids) is never silently substituted.
CANDIDATES: Dict[str, Tuple[str, str, str, float, str]] = {
    "log10 DRNAREA": ("DRNAREA", "DRNAREA", "log10", 1.0, "mi2"),
    "log10 PRECIP": ("PRECIP", "PRECIP_PRISM9120", "log10", 1.0, "in"),
    "log10 I24H2Y": ("I24H2Y", "I24H2Y_NA2", "log10", 1.0, "in"),
    "JANMAXT": ("JANMAXT", "JANMAXT_PRISM9120", "identity", 1.0, "degF"),
    "ELEV/1000": ("ELEV", "ELEV_3DEP", "identity", 0.001, "ft"),
    "log10 BSLOPD": ("BSLOPD", "BSLOPD_3DEP", "log10", 1.0, "degrees"),
    "log10 SOILKSAT": ("SOILKSAT", "SOILKSAT_STATSGO2", "log10", 1.0, "in/hr"),
}
FORCED = ("log10 DRNAREA",)

#: Base regions by HUC prefix (longest prefix wins); anything unlisted is East.
BASE_REGIONS: Dict[str, str] = {"1708": "W", "1709": "W", "1710": "W", "180101": "W"}
#: Region splits, filled by :func:`split_regions`: HUC6 -> region code.
SPLITS: Dict[str, str] = {}
REGION_NAMES = {
    "W": "West of the Cascade crest (HUC4 1708, 1709, 1710; HUC6 180101), less split HUC6s",
    "E": "East of the Cascade crest (all other HUCs), less split HUC6s",
}
MIN_REGION_SITES = 25
#: WBD names of the HUC6s a split can take.
HUC6_NAMES = {
    "170701": "Middle Columbia",
    "170800": "Lower Columbia",
    "170900": "Willamette",
    "171001": "Washington Coastal",
    "171002": "Northern Oregon Coastal",
    "171003": "Southern Oregon Coastal",
    "180101": "Northern California Coastal",
}


def region_of(huc: str) -> str:
    table = {**BASE_REGIONS, **SPLITS}
    for k in sorted(table, key=len, reverse=True):
        if str(huc).startswith(k):
            return table[k]
    return "E"


def load() -> pd.DataFrame:
    st = pd.read_csv(STATIONS, dtype={"site_no": str, "huc": str})
    st = st[st["excluded"].isna()]
    bc = pd.read_csv(BASINS, dtype={"site_no": str})
    df = st.merge(bc, on="site_no", how="inner", suffixes=("", "_basin"))
    need = [c[0] for c in CANDIDATES.values()]
    df = df[df["basin_source"].notna()].dropna(subset=need).reset_index(drop=True)
    df["region"] = df["huc"].map(region_of)
    return df


def column(df: pd.DataFrame, name: str) -> np.ndarray:
    col, _, transform, scale, _ = CANDIDATES[name]
    x = df[col].to_numpy(dtype=float) * scale
    return np.log10(x) if transform == "log10" else x


def design(df: pd.DataFrame, names: Sequence[str]) -> np.ndarray:
    return np.column_stack([np.ones(len(df))] + [column(df, n) for n in names])


def sigma(df: pd.DataFrame, X: np.ndarray, aep: float) -> np.ndarray:
    coef, *_ = np.linalg.lstsq(X, df["sd_log"].to_numpy(dtype=float), rcond=None)
    sd = np.clip(X @ coef, 0.02, None)
    n = df["n_peaks"].to_numpy(dtype=float)
    skew = df["skew_weighted"].fillna(df["skew_station"]).to_numpy(dtype=float)
    mse = np.where(df["skew_weighted"].notna(), PNW_SKEW_MSE, AT_SITE_SKEW_MSE)
    var = np.concatenate(
        [lp3_quantile_sampling_variance(sd[i : i + 1], skew[i : i + 1], n[i : i + 1], aep, mse[i])
         for i in range(len(df))]
    )  # fmt: skip
    y0 = df["first_wy"].to_numpy(dtype=float)
    y1 = df["last_wy"].to_numpy(dtype=float)
    # Overlap of the record spans, capped at each record's peak count: a record
    # with gaps cannot share more years than it has.
    conc = np.clip(np.minimum.outer(y1, y1) - np.maximum.outer(y0, y0) + 1.0, 0.0, None)
    conc = np.minimum(conc, np.minimum.outer(n, n))
    rho = PNW_CORRELATION.rho(
        distance_miles(df["lat_cent"].to_numpy(dtype=float), df["lon_cent"].to_numpy(dtype=float))
    )
    s = rho * np.sqrt(np.outer(var, var)) * conc / np.sqrt(np.outer(n, n))
    np.fill_diagonal(s, var)
    s = (s + s.T) / 2.0
    # The estimator need not be positive definite; clip to the nearest matrix that is.
    w, v = np.linalg.eigh(s)
    if w.min() < 1e-10:
        s = (v * np.clip(w, 1e-10, None)) @ v.T
    return np.asarray(s, dtype=float)


def y_of(df: pd.DataFrame, aep: float) -> np.ndarray:
    return np.log10(df[f"q_{aep}"].to_numpy(dtype=float))


def vif(X: np.ndarray) -> List[float]:
    """Variance inflation factors of the non-intercept columns."""
    Z = X[:, 1:]
    out = []
    for j in range(Z.shape[1]):
        others = np.column_stack([np.ones(len(Z)), np.delete(Z, j, axis=1)])
        coef, *_ = np.linalg.lstsq(others, Z[:, j], rcond=None)
        resid = Z[:, j] - others @ coef
        r2 = 1.0 - resid.var() / Z[:, j].var() if Z[:, j].var() > 0 else 1.0
        out.append(float(1.0 / max(1.0 - r2, 1e-12)))
    return out


def fit(df: pd.DataFrame, names: Sequence[str], aep: float) -> GLSResult:
    X = design(df, names)
    return fit_gls(y_of(df, aep), X, sigma(df, X, aep), names=["const", *names])


def t_values(res: GLSResult) -> np.ndarray:
    return np.asarray(res.beta / np.sqrt(np.diag(res.cov_beta)))


def select(df: pd.DataFrame) -> Tuple[List[str], List[Dict[str, Any]]]:
    """Forward selection at :data:`SELECT_AEP`: |t| >= 2 and VIF < 5 for every term."""
    names = list(FORCED)
    log: List[Dict[str, Any]] = []
    while True:
        best: Optional[Tuple[float, str]] = None
        for cand in CANDIDATES:
            if cand in names:
                continue
            trial = names + [cand]
            X = design(df, trial)
            if max(vif(X)) >= 5.0:
                log.append({"add": cand, "rejected": "VIF >= 5"})
                continue
            res = fit(df, trial, SELECT_AEP)
            t = t_values(res)[1:]
            if np.any(np.abs(t) < 2.0):
                log.append({"add": cand, "rejected": "|t| < 2", "t": t.round(2).tolist()})
                continue
            if best is None or res.avp < best[0]:
                best = (res.avp, cand)
        if best is None:
            return names, log
        names.append(best[1])
        log.append({"added": best[1], "avp": best[0]})


def huc6_clusters(
    df: pd.DataFrame, res: GLSResult, X: np.ndarray, aep: float
) -> List[Dict[str, Any]]:
    """Mean residual per HUC6 with its standard error (independent-sites approximation)."""
    lam = np.diag(sigma(df, X, aep)) + res.model_error_variance
    out = []
    for h6, idx in df.groupby(df["huc"].str[:6]).indices.items():
        if len(idx) < 5:
            continue
        r = float(res.residuals[idx].mean())
        se = float(np.sqrt(lam[idx].mean() / len(idx)))
        out.append({"huc6": h6, "n": len(idx), "mean_resid": r, "t": r / se})
    return sorted(out, key=lambda c: -abs(c["t"]))


def loocv_sep(df: pd.DataFrame, names: Sequence[str], aep: float) -> float:
    """Leave-one-out SEP (percent): each site predicted by the GLS fitted without it."""
    X = design(df, names)
    y = y_of(df, aep)
    S = sigma(df, X, aep)
    err = []
    for i in range(len(df)):
        keep = np.arange(len(df)) != i
        res = fit_gls(y[keep], X[keep], S[np.ix_(keep, keep)])
        err.append(y[i] - float(X[i] @ res.beta))
    mse = float(np.mean(np.square(err)) - np.mean(np.diag(S)))
    mse = max(mse, 1e-12)
    return 100.0 * math.sqrt(math.exp(math.log(10.0) ** 2 * mse) - 1.0)


def current_seps() -> Dict[str, Dict[str, float]]:
    """Published SEPs (percent) of the current Oregon equations, by region and AEP."""
    d = json.loads(OR_JSON.read_text(encoding="utf-8"))
    out: Dict[str, Dict[str, float]] = {}
    for e in d["equations"]:
        if e.get("sep_log") is None:
            continue
        sep = 100.0 * math.sqrt(math.exp(math.log(10.0) ** 2 * e["sep_log"] ** 2) - 1.0)
        out.setdefault(e["region_code"], {})[str(e["aep"])] = round(sep, 1)
    return out


def split_regions(df: pd.DataFrame) -> List[Dict[str, Any]]:
    """Task 4: split off HUC6 residual clusters until none qualifies.

    Each pass fits every region (variables re-selected) and takes the HUC6
    cluster with the largest |t| at AEP 0.01 among those with |t| >= 2, at least
    :data:`MIN_REGION_SITES` gages, and leaving at least as many behind. It
    becomes region ``<parent>-<HUC6>``. Adopted splits go into :data:`SPLITS`
    and ``df["region"]`` in place.
    """
    steps: List[Dict[str, Any]] = []
    while True:
        best: Optional[Tuple[float, str, str, int]] = None
        for region, g in df.groupby("region"):
            g = g.reset_index(drop=True)
            names, _ = select(g)
            X = design(g, names)
            res = fit(g, names, SELECT_AEP)
            for c in huc6_clusters(g, res, X, SELECT_AEP):
                ok = (
                    abs(c["t"]) >= 2.0
                    and c["n"] >= MIN_REGION_SITES
                    and len(g) - c["n"] >= MIN_REGION_SITES
                )
                if ok and (best is None or abs(c["t"]) > best[0]):
                    best = (abs(c["t"]), str(region), c["huc6"], c["n"])
        if best is None:
            return steps
        t, parent, huc6, n = best
        code = f"{parent}-{huc6}"
        SPLITS[huc6] = code
        REGION_NAMES[code] = f"{HUC6_NAMES.get(huc6, 'HUC6')} (HUC6 {huc6}), split from {parent}"
        df["region"] = df["huc"].map(region_of)
        steps.append({"split": huc6, "from": parent, "region": code, "n": n, "t": round(t, 2)})


def fixed_splits(df: pd.DataFrame, huc6s: Sequence[str]) -> List[Dict[str, Any]]:
    """Adopt the given HUC6 splits, in order, without re-deriving them (``--splits``)."""
    steps: List[Dict[str, Any]] = []
    for huc6 in huc6s:
        parent = region_of(huc6)
        code = f"{parent}-{huc6}"
        SPLITS[huc6] = code
        REGION_NAMES[code] = f"{HUC6_NAMES.get(huc6, 'HUC6')} (HUC6 {huc6}), split from {parent}"
        df["region"] = df["huc"].map(region_of)
        n = int((df["region"] == code).sum())
        steps.append({"split": huc6, "from": parent, "region": code, "n": n, "t": None})
    return steps


def run(store: bool, splits_given: Optional[Sequence[str]] = None) -> Dict[str, Any]:
    df = load()
    splits = split_regions(df) if splits_given is None else fixed_splits(df, splits_given)
    out: Dict[str, Any] = {
        "status": "provisional -- NOT a USGS study",
        "n_sites": len(df),
        "n_by_region": df["region"].value_counts().to_dict(),
        "region_splits": splits,
        "regions": {},
        "current_equation_sep_percent": current_seps(),
    }
    equations: List[Dict[str, Any]] = []
    for region, g in sorted(df.groupby("region")):
        g = g.reset_index(drop=True)
        names, sel_log = select(g)
        X = design(g, names)
        res01 = fit(g, names, SELECT_AEP)
        reg: Dict[str, Any] = {
            "name": REGION_NAMES.get(region, region),
            "n": len(g),
            "variables": names,
            "selection": sel_log,
            "vif": vif(X),
            "huc6_clusters_aep_0.01": huc6_clusters(g, res01, X, SELECT_AEP),
            "fits": {},
        }
        for aep in AEPS:
            res = fit(g, names, aep)
            reg["fits"][str(aep)] = {
                "beta": res.beta.tolist(),
                "t": t_values(res).round(2).tolist(),
                "cov_beta": res.cov_beta.tolist(),
                "mev": res.model_error_variance,
                "avp": res.avp,
                "sep_percent": res.sep_percent,
                "pseudo_r2": res.pseudo_r2,
                "max_leverage": float(res.leverage.max()),
                "loocv_sep_percent": loocv_sep(g, names, aep) if aep in (0.5, 0.01) else None,
            }
            equations.append(equation_json(region, g, names, aep, res))
        out["regions"][region] = reg
    OUT.write_text(json.dumps(out, indent=1, default=float) + "\n", encoding="utf-8")
    if store:
        write_provisional(equations, df)
    return out


def equation_json(
    region: str, g: pd.DataFrame, names: Sequence[str], aep: float, res: GLSResult
) -> Dict[str, Any]:
    variables = []
    for n in names:
        col, code, transform, scale, units = CANDIDATES[n]
        v: Dict[str, Any] = {
            "code": code,
            "transform": transform,
            "units": units,
            "minimum": float(g[col].min()),
            "maximum": float(g[col].max()),
        }
        if scale != 1.0:
            v["scale"] = scale
        variables.append(v)
    return {
        "region_code": region,
        "region_name": REGION_NAMES.get(region, region),
        "aep": aep,
        "intercept": float(res.beta[0]),
        "variables": variables,
        "coefficients": [float(b) for b in res.beta[1:]],
        "citation": {
            "publication": "flowfreq provisional Oregon regression (NOT a USGS study), "
            "tools/oregon_regression_fit.py",
            "table": "data/regression_study/oregon_regression_results.json",
            "year": 2026,
        },
        "sep_log": math.sqrt(res.avp),
        "avp": res.avp,
        "model_error_variance": res.model_error_variance,
        "covariance": res.cov_beta.tolist(),
        "n_sites": len(g),
    }


def write_provisional(equations: List[Dict[str, Any]], df: pd.DataFrame) -> None:
    hucs = {**BASE_REGIONS, **SPLITS}
    data = {
        "schema_version": 1,
        "state": "OR",
        "status": "provisional",
        "source_reports": [
            "flowfreq provisional Oregon regression update (2026), NOT a USGS study; "
            "docs/OREGON_REGRESSION_UPDATE.md"
        ],
        "equations": equations,
        "notes": (
            "Developed by flowfreq following the USGS GLS procedure (WREG). Regions are "
            f"hydrologic-unit sets: {json.dumps(hucs)}; any other HUC is region E. "
            "Variable codes name their sources (PRISM 1991-2020 normals, NOAA Atlas 2, "
            "3DEP, STATSGO2 Ksat): values of same-named StreamStats characteristics are not "
            "interchangeable. See docs/OREGON_REGRESSION_UPDATE.md."
        ),
    }
    PROVISIONAL.parent.mkdir(parents=True, exist_ok=True)
    PROVISIONAL.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")


def main(argv: Sequence[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--store", action="store_true")
    ap.add_argument(
        "--splits",
        default=None,
        help="comma-separated HUC6 splits to adopt in order instead of re-deriving them",
    )
    args = ap.parse_args(list(argv))
    out = run(args.store, None if args.splits is None else args.splits.split(","))
    for region, reg in out["regions"].items():
        f = reg["fits"]["0.01"]
        print(region, reg["n"], reg["variables"], "SEP %.1f%%" % f["sep_percent"],
              "LOOCV %.1f%%" % f["loocv_sep_percent"], "pseudo-R2 %.3f" % f["pseudo_r2"])  # fmt: skip
        for c in reg["huc6_clusters_aep_0.01"][:6]:
            print("   ", c)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
