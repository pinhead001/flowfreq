"""Store the 2026 urban flood-frequency equations: Alabama, Mississippi and Tennessee, and
Tennessee's report in Georgia, North Carolina and South Carolina as well.

Neither report prints a coefficient: the equations are in their ScienceBase data releases,
as WREG output. This reads the unpacked releases, checks them against the releases' own
per-gage results, and writes the equations into the state files as region ``Urban 2026``.

- SIR 2026-5044 (Wagner and Anderson, 2026), Alabama and Mississippi; data release
  doi:10.5066/P13GPMON (``WREG.zip``: ``RegressionEquations.csv``, ``WLSresults.csv``,
  ``WLS_output/AEP_*/``). Weighted least squares on 42 streamgages. The release has no
  covariance matrix, so it is computed here as (X' W^-1 X)^-1 from ``modelInput.txt`` and
  ``Weighting.txt``, and accepted only if it reproduces the printed coefficients, their
  standard errors and every gage's variance of prediction.
- SIR 2025-5104 (Wagner and Ladd, 2026), Tennessee; data release doi:10.5066/P1TE9KTH
  (``GLS_results.zip``: ``GLSresults.csv``, ``AEP_*/coefs.txt``, ``CoVarMat.txt``,
  ``PerformanceMetrics.txt``). Generalized least squares on 136 streamgages, 86 of them in
  Georgia (32), North Carolina (38) and South Carolina (16), where the same set is stored by
  the owner's decision of 2026-11-23: the report doesn't say it applies there in place of
  those states' own equations.

Usage::

    python tools/build_urban_2026.py --alms <dir with RegressionEquations.csv> \
        --tn <GLS_results dir>

Every check failure stops the run before a file is written.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA = REPO_ROOT / "flowfreq" / "data" / "regression"
REGION = "Urban 2026"
AEPS: Tuple[str, ...] = ("0.5", "0.2", "0.1", "0.04", "0.02", "0.01", "0.005", "0.002")

ALMS_REPORT = (
    "Wagner, D.M., and Anderson, B.T., 2026, Magnitude and frequency of urban floods in "
    "Alabama and Mississippi using data through the 2024 water year: U.S. Geological Survey "
    "Scientific Investigations Report 2026-5044, 13 p., https://doi.org/10.3133/sir20265044."
)
TN_REPORT = (
    "Wagner, D.M., and Ladd, D.E., 2026, Estimating the magnitude and frequency of floods at "
    "ungaged locations on urban streams in Tennessee and parts of Alabama, Georgia, "
    "Mississippi, North Carolina, and South Carolina, using data through the 2022 water year: "
    "U.S. Geological Survey Scientific Investigations Report 2025-5104, 17 p., "
    "https://doi.org/10.3133/sir20255104."
)


def _metrics(path: Path) -> Dict[str, float]:
    frame = pd.read_csv(path, sep="\t", skiprows=1, header=None, index_col=0)
    return {str(k): float(v) for k, v in frame[1].items()}


def _coefficients(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def _design(path: Path) -> Tuple[List[str], np.ndarray]:
    """Site numbers and the X matrix (constant first) of a WREG ``modelInput.txt``."""
    frame = pd.read_csv(path, sep="\t", dtype={"selected_sites": str})
    columns = [c for c in frame.columns if c not in ("selected_sites", "Y")]
    return list(frame["selected_sites"]), frame[columns].to_numpy(dtype=float)


def sep_log(sp_percent: float) -> float:
    """A standard error of prediction in percent as log10 units (the library's convention)."""
    return math.sqrt(math.log(1.0 + (sp_percent / 100.0) ** 2)) / math.log(10.0)


def _same_3sf(a: float, b: float) -> bool:
    """True when ``a`` rounds to the three significant figures printed as ``b``."""
    if b == 0:
        return abs(a) < 0.5
    step = 10.0 ** (math.floor(math.log10(abs(b))) - 2)
    return abs(a - b) <= 0.5 * step * (1 + 1e-9)


def _check_gages(
    label: str,
    sites: Sequence[str],
    x: np.ndarray,
    beta: np.ndarray,
    cov: np.ndarray,
    mev: float,
    printed: pd.DataFrame,
    flow_column: str,
    variance_column: str,
) -> None:
    flows = 10.0 ** (x @ beta)
    variances = mev + np.einsum("ij,jk,ik->i", x, cov, x)
    table = printed.set_index("site_no")
    for site, flow, variance in zip(sites, flows, variances):
        row = table.loc[site]
        if not _same_3sf(float(flow), float(row[flow_column])):
            raise SystemExit(f"{label} {site}: flow {flow:.6g} is not {row[flow_column]}")
        if abs(variance - float(row[variance_column])) > 5.1e-5:
            raise SystemExit(
                f"{label} {site}: variance of prediction {variance:.6f} is not "
                f"{row[variance_column]}"
            )


def alms_equations(root: Path) -> List[Dict[str, Any]]:
    """The eight SIR 2026-5044 equations, checked against the release's 42 gages."""
    printed = pd.read_csv(root / "WLSresults.csv", dtype={"site_no": str})
    summary = pd.read_csv(root / "RegressionEquations.csv", dtype=str).dropna(subset=["AEP"])
    summary = summary.set_index("AEP")
    out = []
    for aep in AEPS:
        folder = root / "WLS_output" / f"AEP_{aep}"
        coefs = _coefficients(folder / "coefs.txt")
        metrics = _metrics(folder / "PerformanceMetrics.txt")
        sites, x = _design(folder / "modelInput.txt")
        weights = pd.read_csv(folder / "Weighting.txt", sep="\t", index_col=0).to_numpy(float)
        if not np.allclose(weights, np.diag(np.diag(weights))):
            raise SystemExit(f"AL/MS AEP {aep}: the weighting matrix is not diagonal")
        inverse = np.diag(1.0 / np.diag(weights))
        cov = np.linalg.inv(x.T @ inverse @ x)
        y = pd.read_csv(folder / "modelInput.txt", sep="\t")["Y"].to_numpy(float)
        beta = cov @ x.T @ inverse @ y
        printed_beta = coefs["Coefficient"].to_numpy(float)
        if not np.allclose(beta, printed_beta, rtol=0, atol=1e-10):
            raise SystemExit(f"AL/MS AEP {aep}: refit {beta} is not the release's {printed_beta}")
        if not np.allclose(np.sqrt(np.diag(cov)), coefs["Standard Error"], rtol=1e-8):
            raise SystemExit(f"AL/MS AEP {aep}: standard errors are not the release's")
        mev = metrics["ModErrVar"]
        _check_gages(
            f"AL/MS AEP {aep}",
            sites,
            x,
            printed_beta,
            cov,
            mev,
            printed,
            f"WLS_{aep}",
            f"VarPred_{aep}",
        )
        variances = mev + np.einsum("ij,jk,ik->i", x, cov, x)
        if abs(float(variances.mean()) - metrics["AVP"]) > 1e-9:
            raise SystemExit(f"AL/MS AEP {aep}: mean variance of prediction is not the AVP")
        rounded = [float(summary.loc[aep, k]) for k in ("b", "a", "c")]
        if not np.allclose(np.round(printed_beta, 4), rounded, atol=1e-12):
            raise SystemExit(f"AL/MS AEP {aep}: coefs.txt does not round to RegressionEquations")
        out.append(
            {
                "region_code": REGION,
                "region_name": "Urban, Alabama and Mississippi (SIR 2026-5044; not in NSS)",
                "aep": float(aep),
                "intercept": float(printed_beta[0]),
                "variables": [
                    {
                        "code": "DRNAREA",
                        "transform": "log10",
                        "units": "mi2",
                        "minimum": 1.06,
                        "maximum": 53.6,
                    },
                    {
                        "code": "LC11DEV",
                        "transform": "identity",
                        "units": "percent",
                        "minimum": 20.0,
                        "maximum": 100.0,
                    },
                ],
                "coefficients": [float(printed_beta[1]), float(printed_beta[2])],
                "citation": {
                    "publication": ALMS_REPORT,
                    "table": (
                        "Data release doi:10.5066/P13GPMON, WLS_output/AEP_*/coefs.txt and "
                        "PerformanceMetrics.txt; ranges from the report, p. 8"
                    ),
                    "year": 2026,
                    "url": "https://doi.org/10.3133/sir20265044",
                },
                "sep_log": round(sep_log(metrics["Sp"]), 6),
                "avp": metrics["AVP"],
                "model_error_variance": mev,
                "covariance": [[float(v) for v in row] for row in cov],
                "n_sites": len(sites),
            }
        )
    return out


def tn_equations(
    root: Path, region_name: str = "Urban, Tennessee (SIR 2025-5104; not in NSS)"
) -> List[Dict[str, Any]]:
    """The eight SIR 2025-5104 equations, checked against the release's 136 gages."""
    printed = pd.read_csv(root / "GLSresults.csv", dtype={"site_no": str})
    summary = pd.read_csv(root / "RegressionEquations.csv", dtype=str).dropna(subset=["AEP"])
    summary = summary.set_index("AEP")
    names = ["constant", "log10_DRNAREA", "LC11DEV", "L3_Piedmont", "L3_Ridge_Valley"]
    out = []
    for aep in AEPS:
        folder = root / f"AEP_{aep}"
        coefs = _coefficients(folder / "coefs.txt")
        if list(coefs.index) != names:
            raise SystemExit(f"TN AEP {aep}: unexpected terms {list(coefs.index)}")
        metrics = _metrics(folder / "PerformanceMetrics.txt")
        sites, x = _design(folder / "modelInput.txt")
        cov_frame = pd.read_csv(folder / "CoVarMat.txt", sep="\t", index_col=0)
        if list(cov_frame.index) != names or list(cov_frame.columns) != names:
            raise SystemExit(f"TN AEP {aep}: covariance order is {list(cov_frame.index)}")
        cov = cov_frame.to_numpy(float)
        beta = coefs["Coefficient"].to_numpy(float)
        if not np.allclose(np.sqrt(np.diag(cov)), coefs["Standard Error"], rtol=1e-5):
            raise SystemExit(f"TN AEP {aep}: covariance diagonal is not the standard errors")
        mev = metrics["ModErrVar"]
        _check_gages(
            f"TN AEP {aep}", sites, x, beta, cov, mev, printed, f"GLS_{aep}", f"GLSvar_{aep}"
        )
        rounded = [float(summary.loc[aep, k]) for k in ("b", "a", "c", "d", "e")]
        if not np.allclose(np.round(beta, 4), rounded, atol=1e-12):
            raise SystemExit(f"TN AEP {aep}: coefs.txt does not round to RegressionEquations")
        percent = {"transform": "identity", "units": "percent"}
        out.append(
            {
                "region_code": REGION,
                "region_name": region_name,
                "aep": float(aep),
                "intercept": float(beta[0]),
                "variables": [
                    {
                        "code": "DRNAREA",
                        "transform": "log10",
                        "units": "mi2",
                        "minimum": 0.164,
                        "maximum": 93.4,
                    },
                    {"code": "LC11DEV", **percent, "minimum": 26.0, "maximum": 100.0},
                    {"code": "L3_PIEDMNT", **percent, "minimum": 0.0, "maximum": 100.0},
                    {"code": "L3_RDGVLY", **percent, "minimum": 0.0, "maximum": 100.0},
                ],
                "coefficients": [float(v) for v in beta[1:]],
                "citation": {
                    "publication": TN_REPORT,
                    "table": (
                        "Data release doi:10.5066/P1TE9KTH (ver. 1.1), GLS_results/AEP_*/coefs.txt, "
                        "CoVarMat.txt and PerformanceMetrics.txt; ranges from the report, p. 12"
                    ),
                    "year": 2026,
                    "url": "https://doi.org/10.3133/sir20255104",
                },
                "sep_log": round(sep_log(metrics["Sp"]), 6),
                "avp": metrics["AVP"],
                "model_error_variance": mev,
                "covariance": [[float(v) for v in row] for row in cov],
                "n_sites": len(sites),
                "applicability": [
                    {"code": "LC11IMP", **percent, "minimum": 10.0, "maximum": 100.0}
                ],
            }
        )
    return out


ALMS_SOURCE = {
    "citation": ALMS_REPORT,
    "doi": "10.3133/sir20265044",
    "url": "https://pubs.usgs.gov/publication/sir20265044",
    "pdf": "https://pubs.usgs.gov/sir/2026/5044/sir20265044.pdf",
    "data_release": "https://doi.org/10.5066/P13GPMON",
    "tables": ["Data release: RegressionEquations.csv, WLSresults.csv, WLS_output"],
    "nss": {"regression_regions": [], "last_year_of_data": 2024, "checked": "2026-11-18"},
}
TN_SOURCE = {
    "citation": TN_REPORT,
    "doi": "10.3133/sir20255104",
    "url": "https://pubs.usgs.gov/publication/sir20255104",
    "pdf": "https://pubs.usgs.gov/sir/2025/5104/sir20255104.pdf",
    "data_release": "https://doi.org/10.5066/P1TE9KTH",
    "tables": ["Data release: RegressionEquations.csv, GLSresults.csv, AEP_*/CoVarMat.txt"],
    "nss": {"regression_regions": [], "last_year_of_data": 2022, "checked": "2026-11-18"},
}

ALMS_NOTE = (
    " Urban 2026: Wagner and Anderson, SIR 2026-5044, one set for the urban areas of Alabama "
    "and Mississippi at 8 AEPs, Q = DRNAREA^a * 10^(b + c*LC11DEV), by weighted least squares "
    "on 42 streamgages (29 in Alabama, 13 in Mississippi). The report prints no equation, "
    "coefficient or worked example; everything is from its data release "
    "(doi:10.5066/P13GPMON), written by tools/build_urban_2026.py. Coefficients are "
    "coefs.txt's at full precision (RegressionEquations.csv rounds them to 4 decimals, "
    "which moves a flow by up to 1.1 percent); model_error_variance is ModErrVar, avp is AVP "
    "and sep_log is Sp (percent) converted as sqrt(ln(1 + (Sp/100)^2)) / ln(10). The release "
    "has no covariance matrix. The one stored is (X' W^-1 X)^-1 computed from the release's "
    "modelInput.txt and Weighting.txt (order constant, log10 DRNAREA, LC11DEV), and is not a "
    "USGS-printed value: it reproduces the release's coefficients to 1e-10, their standard "
    "errors, and with the model error variance every gage's VarPred in WLSresults.csv at all "
    "8 AEPs to its 4 decimals; its mean equals the AVP. Flows reproduce at all 42 gages and "
    "8 AEPs to the 3 significant figures printed. Ranges are the report's (p. 8): DRNAREA "
    "1.06 to 53.6 mi2 and LC11DEV 20 to 100 percent (developed land, NLCD 2011 classes "
    "21-24); not for streams substantially affected by regulation. The LC11DEV coefficient "
    "falls with AEP and is -0.0002 at 0.2 percent. NSS served no region for this report on "
    "2026-11-18, so nothing locates this set: name it with region='Urban 2026'."
)
AL_NOTE = (
    " The report says it supersedes SIR 2010-5012 (cover and p. 8), the stored region "
    "'Urban', which NSS and StreamStats still serve as GC1504; that set stays, with a caution, "
    "until NSS changes. Unlike it, the 2026 set states no exclusion of flood region 3."
)
MS_NOTE = (
    " The report says the set supersedes Mississippi's 1991 urban equations (Landers and "
    "Wilson), which were never stored here. Regions A, B and C stop at 20 percent LC11DEV, "
    "where this set starts. The report does not say how it relates to the Mississippi "
    "Alluvial Plain set."
)
TN_NOTE = (
    " Urban 2026: Wagner and Ladd, SIR 2025-5104 (published 2026-02-05), one set at 8 AEPs "
    "for urban streams in Tennessee and the parts of AL, GA, MS, NC and SC in six Level 3 "
    "ecoregions, Q = DRNAREA^a * 10^(b + c*LC11DEV + d*L3_PIEDMNT + e*L3_RDGVLY), by "
    "generalized least squares on 136 streamgages (20 in Tennessee). The report names "
    "only Tennessee for StreamStats and says nothing of replacing another state's equations; "
    "the set is stored in GA, NC and SC too, by the owner's decision of 2026-11-23. The report prints no equation or coefficient; "
    "everything is from its data release (doi:10.5066/P1TE9KTH, ver. 1.1), written by "
    "tools/build_urban_2026.py. Coefficients are coefs.txt's at full precision (the 4 "
    "decimals of RegressionEquations.csv move a flow by up to 1.9 percent, since the "
    "ecoregion terms multiply by up to 100); covariance is CoVarMat.txt's (order constant, "
    "log10 DRNAREA, LC11DEV, L3_PIEDMNT, L3_RDGVLY), model_error_variance is ModErrVar, avp is "
    "AVP and sep_log is Sp converted as elsewhere. Flows and variances of prediction "
    "reproduce GLSresults.csv at all 136 gages and 8 AEPs to the figures printed. L3_PIEDMNT "
    "and L3_RDGVLY are the percent of the basin in the Piedmont and in the Ridge and Valley "
    "Level 3 ecoregions (EPA, 2013); a Tennessee basin has no Piedmont. Ranges (report p. "
    "12): DRNAREA 0.164 to 93.4 mi2, LC11DEV 26 to 100 percent. Applicability: 10 percent or "
    "more developed imperviousness (LC11IMP, NLCD 2011), below which the report sends "
    "Tennessee sites to the rural equations of SIR 2024-5130; not for regulated or tidal "
    "streams, nor the Middle Atlantic or Southern Coastal Plain. The report does not say it "
    "supersedes WRI 84-4182 or WRI 84-4110 ('Urban statewide' and 'Memphis urban'), which NSS "
    "still serves as GC351 and GC350 and which stay, with a caution. Its abstract says the "
    "equations were incorporated into StreamStats and its summary says they will be; NSS "
    "served no region for it on 2026-11-18, so nothing locates this set: name it with "
    "region='Urban 2026'. Whether StreamStats computes the two ecoregion percentages for "
    "Tennessee is not known."
)

AL_CAUTION = (
    "SIR 2026-5044 (Wagner and Anderson, 2026) says it supersedes these urban equations "
    "(SIR 2010-5012). NSS and StreamStats still serve this set, so it is what a StreamStats "
    "region lookup returns; the 2026 set is stored as region 'Urban 2026' and takes DRNAREA "
    "and LC11DEV."
)
TN_CAUTION = (
    "SIR 2025-5104 (Wagner and Ladd, 2026) updated the methods for urban streams in Tennessee "
    "with data through water year 2022, and says that where the basin is 10 percent or more "
    "developed imperviousness its equations 'should be used to predict streamflows' (p. 12). "
    "They are stored as region 'Urban 2026'. The report does not name this set or say it "
    "supersedes it, and NSS and StreamStats still serve it."
)

#: What any estimate with SIR 2025-5104's equations should carry, whatever the state.
LIMITS_CAUTION = (
    " It is not valid in the Middle Atlantic or Southern Coastal Plain Level 3 ecoregions or for "
    "regulated or tidal streams; its accuracy is not known outside the study area, which is "
    "six ecoregions and had no urban gage in the Mississippi Alluvial Plain, Central "
    "Appalachians or Blue Ridge (pp. 5, 12). Nothing here checks "
    "which ecoregion the basin is in. L3_PIEDMNT and L3_RDGVLY are percent of the basin, 0 to "
    "100, not fractions. At a gage used to develop the equations the report calls weighting "
    "them with the at-site estimate inappropriate, because the two are not independent (p. 13)."
)
TN_2026_CAUTION = (
    "SIR 2025-5104 (Wagner and Ladd, 2026) is an urban set fitted to gages in six southeastern "
    "states, 20 of 136 in Tennessee." + LIMITS_CAUTION
)


#: For Georgia, North Carolina and South Carolina. ``{gages}`` is the state's count of the 136.
SE_NOTE = (
    " Urban 2026: Wagner and Ladd, SIR 2025-5104 (published 2026-02-05), the set described in "
    "TN.json, identical here: 8 AEPs, Q = DRNAREA^a * 10^(b + c*LC11DEV + d*L3_PIEDMNT + "
    "e*L3_RDGVLY), by generalized least squares on 136 streamgages, {gages} of them in this "
    "state, from its data release (doi:10.5066/P1TE9KTH, ver. 1.1) by "
    "tools/build_urban_2026.py. It is here by the owner's decision of 2026-11-23, for urban "
    "basins above 1 mi2 and 10 percent impervious, where the other sets stored here are SIR "
    "2014-5030's (Urban 2014); the report does not say it applies in place of this state's equations, and names only "
    "Tennessee for StreamStats. It is in no region_preference tier and NSS locates no region "
    "for it, so it is used only when named: region='Urban 2026', which named_regions lets "
    "stand in place of what NSS located, for the whole basin. L3_PIEDMNT and L3_RDGVLY are "
    "the percent of the basin in the Piedmont and in the Ridge and Valley Level 3 ecoregions "
    "(EPA, 2013), which the caller supplies (supplied_characteristics) unless StreamStats has "
    "them; basin_shares refuses a sum over 100. Ranges (report p. "
    "12): DRNAREA 0.164 to 93.4 mi2, LC11DEV 26 to 100 percent; applicability 10 percent or "
    "more developed imperviousness (LC11IMP). Not valid in the Middle Atlantic or Southern "
    "Coastal Plain Level 3 ecoregions, nor for regulated or tidal streams; no urban gage was "
    "in the Blue Ridge (p. 5)."
)
SC_NOTE = (
    " At the report's 16 South Carolina gages the geometric mean of the gage's flow over the "
    "equation's is 0.69 at 10 percent AEP and 0.76 at 1 percent (GLSresults.csv), and at the "
    "8 of them outside the Piedmont 0.65 and 0.70."
)
SE_CAUTION = (
    "SIR 2025-5104 (Wagner and Ladd, 2026) is an urban set fitted to gages in six southeastern "
    "states; it does not say it applies in place of this state's equations. The USGS also "
    "publishes equations for urban basins here by hydrologic region, SIR 2014-5030 (Feaster, "
    "Gotvald and Weaver, 2014, Table 7), to 436 mi2 in HR1, 459 in HR3 and 53.5 in HR4, "
    "though the largest urban basin it used was 378, 5.67 and 43.9 mi2 (Table 11); "
    "flowfreq has them as the 'Urban 2014' regions, so compare with those." + LIMITS_CAUTION
)
SC_CAUTION = (
    " It estimates high in this state: at-site flow over the equations' flow at its 16 South "
    "Carolina gages has a geometric mean of 0.69 at 10 percent AEP (0.36 to 1.20, 14 of 16 "
    "below 1) and 0.76 at 1 percent AEP (0.33 to 2.71, 12 of 16 below 1)."
)
SOUTHEAST = {"GA": 32, "NC": 38, "SC": 16}


def _store(
    state: str,
    equations: List[Dict[str, Any]],
    source: Dict[str, Any],
    note: str,
    cautions: Dict[str, str],
    shares: bool = True,
) -> None:
    path = DATA / f"{state}.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    data["equations"] = [e for e in data["equations"] if e["region_code"] != REGION] + equations
    reports = [r for r in data["source_reports"] if r.get("doi") != source["doi"]]
    data["source_reports"] = reports + [source]
    marker = f" {REGION}: "
    base = data.get("notes", "")
    data["notes"] = (base.split(marker)[0] if marker in base else base) + note
    if cautions:
        data["region_cautions"] = {**data.get("region_cautions", {}), **cautions}
    others = [r for r in data.get("named_regions") or [] if r != REGION]
    data["named_regions"] = [REGION] + others
    if shares:
        data["basin_shares"] = [["L3_PIEDMNT", "L3_RDGVLY"]]
    # The state files are indent-1 ASCII JSON with CRLF line ends.
    text = json.dumps(data, indent=1) + "\n"
    path.write_bytes(text.replace("\n", "\r\n").encode("ascii"))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--alms", type=Path, required=True, help="unpacked WREG.zip (AL/MS)")
    parser.add_argument("--tn", type=Path, required=True, help="unpacked GLS_results (TN)")
    args = parser.parse_args(argv)
    alms = alms_equations(args.alms)
    tennessee = tn_equations(args.tn)
    _store("AL", alms, ALMS_SOURCE, ALMS_NOTE + AL_NOTE, {"Urban": AL_CAUTION}, False)
    _store("MS", alms, ALMS_SOURCE, ALMS_NOTE + MS_NOTE, {}, False)
    _store(
        "TN",
        tennessee,
        TN_SOURCE,
        TN_NOTE,
        {
            "Urban statewide": TN_CAUTION,
            "Memphis urban": TN_CAUTION,
            REGION: TN_2026_CAUTION,
        },
    )
    southeast = tn_equations(args.tn, "Urban, southeastern states (SIR 2025-5104; not in NSS)")
    for state, gages in SOUTHEAST.items():
        south_carolina = state == "SC"
        _store(
            state,
            southeast,
            TN_SOURCE,
            SE_NOTE.format(gages=gages) + (SC_NOTE if south_carolina else ""),
            {REGION: SE_CAUTION + (SC_CAUTION if south_carolina else "")},
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
