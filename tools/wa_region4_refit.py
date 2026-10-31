"""Refit Washington SIR 2016-5118's regional regressions with flowfreq's GLS.

Northwest studies plan, task 4. Two purposes:

1. Validation of :mod:`flowfreq.regression.gls` against a published USGS
   regression. Regions 1-3 are refit at all eight AEPs and compared with the
   stored ``WA.json`` (Tables 6 and 7).
2. The Region 4 covariance basis (#65). Table 7 prints Region 4's
   ``[X' Lambda^-1 X]^-1`` with no word on the scale of P. The refit is done
   with raw P and with P/10, and each is compared with the printed matrix.

Inputs (both under ``data/regression_study/``):

- ``wa_2016_5118_sites.csv`` (tools/build_wa_gls_dataset.py): Table 5
  characteristics, Table 2 record spans, Table 8 station quantiles.
- ``wa_2016_5118_idor_ema.csv``: station quantiles for the 51 Idaho and Oregon
  calibration gages Table 8 does not print. They are computed with flowfreq EMA
  (peaks through WY2014, Pacific Northwest skew -0.07, MSE 0.18), the report's
  procedure (p. 13). Built by ``tools/build_wa_gls_dataset.py --ema``.

GLS specification (WREG, Eng and others 2009). It follows the report (p. 31)
where the report says what it did:

- y: log10 station quantile (Table 8 line 1).
- Sampling variance: :func:`~flowfreq.regression.gls.lp3_quantile_sampling_variance`.
  It uses the station SD smoothed by OLS on the region's design, the station
  weighted skew, a skew MSE of 0.18 (the PNW regional skew), and n = the
  historical period length (Table 2).
- Cross-correlation: the PNW concurrent-peak model (SIR 2016-5083 eq. B11,
  ``Z = 0.21 + exp(-0.17 - 0.0058 D)``), between basin centroids. Concurrent
  years come from each gage's first and last year in Table 2.

Usage::

    python tools/wa_region4_refit.py    # writes data/regression_study/wa_2016_5118_refit.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

from flowfreq.regression.gls import fit_gls, lp3_quantile_sampling_variance
from flowfreq.skew_study import CorrelationDistanceModel, distance_miles

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "regression_study"
SITES = DATA / "wa_2016_5118_sites.csv"
IDOR = DATA / "wa_2016_5118_idor_ema.csv"
OUT = DATA / "wa_2016_5118_refit.json"
WA_JSON = ROOT / "flowfreq" / "data" / "regression" / "WA.json"

AEPS = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)
PNW_SKEW_MSE = 0.18
PNW_CORRELATION = CorrelationDistanceModel(a=0.21, b=-0.17, c=-0.0058)

#: Table 7, Region 4, upper triangle (a, b, c, d, e, f) of the symmetric 3x3
#: [1, log10 A, P] matrix, and its model error variance, as printed (p. 42).
#: Transcribed by PR #65 for these three AEPs; tests/test_regression_washington.py
#: pins them against the report's workbook and Table 8.
T7_REGION4: Dict[float, Tuple[float, Tuple[float, ...]]] = {
    0.5: (0.045, (4.633e-3, -5.756e-4, -3.314e-4, 4.384e-4, 1.945e-6, 3.270e-5)),
    0.01: (0.047, (6.467e-3, -6.873e-4, -4.498e-4, 5.221e-4, 8.927e-7, 4.273e-5)),
    0.002: (0.053, (7.663e-3, -7.937e-4, -5.319e-4, 6.073e-4, 2.921e-7, 5.033e-5)),
}


def sym3(up: Tuple[float, ...]) -> np.ndarray:
    a, b, c, d, e, f = up
    return np.array([[a, b, c], [b, d, e], [c, e, f]])


def load_sites() -> pd.DataFrame:
    d = pd.read_csv(SITES, dtype={"site_no": str})
    d = d[d["used"]].copy()
    if IDOR.exists():
        e = pd.read_csv(IDOR, dtype={"site_no": str}).set_index("site_no")
        for col in [f"q_{a}" for a in AEPS] + ["sd", "skew_w"]:
            d[col] = d[col].fillna(d["site_no"].map(e[col]))
    d["source"] = np.where(d["site_no"].isin(_idor_sites()), "flowfreq_ema", "table8")
    return d


def _idor_sites() -> List[str]:
    if not IDOR.exists():
        return []
    return list(pd.read_csv(IDOR, dtype={"site_no": str})["site_no"])


def design(df: pd.DataFrame, region: int, p_scale: float = 1.0) -> Tuple[np.ndarray, List[str]]:
    la = np.log10(df["DRNAREA"].to_numpy(dtype=float))
    P = df["PRECIP"].to_numpy(dtype=float)
    one = np.ones(len(df))
    if region in (1, 2):
        return np.column_stack([one, la, P, df["CANOPY_PCT"].to_numpy(dtype=float)]), [
            "const",
            "log10 DRNAREA",
            "PRECIP",
            "CANOPY_PCT",
        ]
    if region == 3:
        return np.column_stack([one, la, np.log10(P)]), ["const", "log10 DRNAREA", "log10 PRECIP"]
    name = "PRECIP" if p_scale == 1.0 else f"PRECIP/{p_scale:g}"
    return np.column_stack([one, la, P / p_scale]), ["const", "log10 DRNAREA", name]


def sigma(df: pd.DataFrame, X: np.ndarray, aep: float) -> np.ndarray:
    coef, *_ = np.linalg.lstsq(X, df["sd"].to_numpy(dtype=float), rcond=None)
    sd = X @ coef
    n = df["hist_len"].to_numpy(dtype=float)
    var = lp3_quantile_sampling_variance(
        sd, df["skew_w"].to_numpy(dtype=float), n, aep, PNW_SKEW_MSE
    )
    y0 = df["first_year"].to_numpy(dtype=float)
    y1 = df["last_year"].to_numpy(dtype=float)
    conc = np.clip(np.minimum.outer(y1, y1) - np.maximum.outer(y0, y0) + 1.0, 0.0, None)
    # A record span can exceed the record length where a record has gaps; concurrent
    # years cannot exceed either record, and more would make Sigma indefinite.
    conc = np.minimum(conc, np.minimum.outer(n, n))
    rho = PNW_CORRELATION.rho(
        distance_miles(df["lat_cent"].to_numpy(dtype=float), df["lon_cent"].to_numpy(dtype=float))
    )
    s = rho * np.sqrt(np.outer(var, var)) * conc / np.sqrt(np.outer(n, n))
    np.fill_diagonal(s, var)
    return np.asarray(s, dtype=float)


def raw_beta(eq: Dict[str, Any]) -> List[float]:
    """``[b0, b1, ...]`` of a stored equation in its raw variables' basis.

    ``WA.json`` stores some variables rescaled (Region 4's PRECPRIS10 as ``scale = 0.1``
    with its coefficient x10, the P/10 basis); the refit's design uses the raw values, so
    an identity term's scale and offset fold into its coefficient and the intercept, and a
    log10 term's into the intercept.
    """
    b0 = float(eq["intercept"])
    out: List[float] = []
    for v, c in zip(eq["variables"], eq["coefficients"]):
        scale, offset = float(v.get("scale", 1.0)), float(v.get("offset", 0.0))
        transform = v.get("transform", "log10")
        if transform not in ("identity", "log10") or (transform == "log10" and scale <= 0):
            raise ValueError(f"{v['code']}: no raw-variable form for {transform} scale {scale}")
        if transform == "identity":
            b0 += c * offset
            out.append(c * scale)
        else:
            if offset:
                raise ValueError(f"{v['code']}: a log10 offset has no raw-variable form")
            b0 += c * float(np.log10(scale))
            out.append(float(c))
    return [b0, *out]


def published(region: int, aep: float) -> Dict[str, Any]:
    wa = json.loads(WA_JSON.read_text(encoding="utf-8"))
    for e in wa["equations"]:
        if e["region_code"] == str(region) and abs(e["aep"] - aep) < 1e-12:
            return dict(e)
    raise KeyError((region, aep))


def refit() -> Dict[str, Any]:
    d = load_sites()
    out: Dict[str, Any] = {"aeps": list(AEPS), "regions": {}}
    for region in (1, 2, 3, 4):
        df = d[(d["region"] == region) & d["q_0.01"].notna()].reset_index(drop=True)
        bases = (1.0, 10.0) if region == 4 else (1.0,)
        reg: Dict[str, Any] = {
            "m": len(df),
            "m_flowfreq_ema": int((df["source"] == "flowfreq_ema").sum()),
        }
        for p_scale in bases:
            X, names = design(df, region, p_scale)
            per_aep = {}
            for aep in AEPS:
                y = np.log10(df[f"q_{aep}"].to_numpy(dtype=float))
                res = fit_gls(y, X, sigma(df, X, aep), names=names)
                pub = published(region, aep)
                entry: Dict[str, Any] = {
                    "beta": res.beta.tolist(),
                    "cov_beta": res.cov_beta.tolist(),
                    "mev": res.model_error_variance,
                    "avp": res.avp,
                    "pub_beta": raw_beta(pub),
                    "pub_mev": pub["model_error_variance"],
                }
                if "covariance" in pub:
                    entry["diag_ratio"] = (
                        np.diag(res.cov_beta) / np.diag(np.asarray(pub["covariance"]))
                    ).tolist()
                if region == 4 and aep in T7_REGION4:
                    t7 = sym3(T7_REGION4[aep][1])
                    entry["diag_ratio_t7"] = (np.diag(res.cov_beta) / np.diag(t7)).tolist()
                per_aep[str(aep)] = entry
            reg["basis_P" if p_scale == 1.0 else f"basis_P/{p_scale:g}"] = per_aep
        out["regions"][str(region)] = reg
    out["region4"] = region4_decision(out["regions"]["4"], d)
    return out


def region4_decision(reg4: Dict[str, Any], d: pd.DataFrame) -> Dict[str, Any]:
    """Which P basis Table 7's Region 4 matrix is in, by the refit and by Table 6's Sp."""
    scores = {}
    for basis in ("basis_P", "basis_P/10"):
        logs = [abs(np.log(reg4[basis][str(a)]["diag_ratio_t7"][2])) for a in T7_REGION4]
        scores[basis] = float(np.median(logs))
    selected = min(scores, key=lambda b: scores[b])
    ratio_p = float(np.median([reg4[selected][str(a)]["diag_ratio_t7"][2] for a in T7_REGION4]))
    df = d[(d["region"] == 4)].reset_index(drop=True)
    la = np.log10(df["DRNAREA"].to_numpy(dtype=float))
    P = df["PRECIP"].to_numpy(dtype=float)
    sp = {}
    for basis, scale in (("basis_P", 1.0), ("basis_P/10", 10.0)):
        X = np.column_stack([np.ones(len(df)), la, P / scale])
        rows = {}
        for aep, (mev, up) in T7_REGION4.items():
            u = sym3(up)
            avp = mev + float(np.mean(np.einsum("ij,jk,ik->i", X, u, X)))
            rows[str(aep)] = {
                "sqrt_avp": float(np.sqrt(avp)),
                "table6_sep_log": published(4, aep)["sep_log"],
            }
            rows[str(aep)]["ratio"] = rows[str(aep)]["sqrt_avp"] / rows[str(aep)]["table6_sep_log"]
        sp[basis] = rows
    return {
        "median_abs_log_ratio_P_entry": scores,
        "selected_basis": selected,
        "selected_P_entry_ratio": ratio_p,
        "table6_sp_check": sp,
    }


def main() -> int:
    out = refit()
    OUT.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    r4 = out["region4"]
    print(f"wrote {OUT.relative_to(ROOT)}")
    print(
        "Region 4 basis:",
        r4["selected_basis"],
        "P-entry ratio",
        round(r4["selected_P_entry_ratio"], 3),
    )
    for region in ("1", "2", "3", "4"):
        reg = out["regions"][region]
        e = reg["basis_P"]["0.01"]
        print(
            f"R{region} m={reg['m']} (ema {reg['m_flowfreq_ema']}) AEP 0.01 beta",
            np.round(e["beta"], 4),
            "pub",
            np.round(e["pub_beta"], 4),
            f"mev {e['mev']:.4f} pub {e['pub_mev']}",
            "diag",
            np.round(e.get("diag_ratio", e.get("diag_ratio_t7", [])), 3),
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
