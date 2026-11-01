"""Build ``flowfreq/data/regression/provisional/WA.json``: Washington's published set with the
Peak Region 2 equations refit with mean basin elevation.

``docs/WA_REGION2_METHOW.md`` found that mean basin elevation (``ELEV``, feet), which the
SIR 2016-5118 Region 2 design leaves out, explains about half of the equation's bias at the
Methow at common floods (t = 4.1 at the 50% AEP, regional rather than Methow-specific) and
fades by the 10% AEP (t = 1.8) and 1% AEP (t = 0.4). The owner chose (2026-10-31) to ship it
as a provisional set for methow-flow: Region 2 becomes
``log Q = b0 + b1 log A + b2 PRECPRIS10 + b3 CANOPY_PCT + b4 ELEV``; every other region is the
published set, unchanged.

**Every Region 2 AEP is refit, not only the 50% and 20%.** The first plan (refit equations at
the common floods, published ones at the rarer) gave a frequency curve that decreased between
the 20% and 10% AEPs at 7 of 90 Region 2 gages, all above about 5,300 ft (Methow above Goat
Creek and the Chewuch among them): the refit equations raise high basins, the published ones
do not. Refit at every AEP, the elevation coefficient tapers smoothly to zero (t from 4.1 at
50% to 0.1 at 0.2%), the curve increases at all 90 gages, and at the 10% AEP and rarer the
refit estimates are a median 0.95-0.99 of the published ones.

The fit is ``tools/wa_region2_elevation.py``'s ``+ELEV`` model (flowfreq EMA quantiles through
WY2025, WREG sampling covariance, ``flowfreq.regression.gls.fit_gls``) on the Region 2
calibration gages with an ``ELEV`` value. Each refit equation carries its GLS covariance, model
error variance and average variance of prediction, so prediction intervals are computed as for
the published equations. Variable limits are the calibration gages' ranges.

Usage::

    python tools/wa_region2_provisional.py          # writes the provisional file (offline)
"""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tools"))

import wa_region2_elevation as r2  # noqa: E402

from flowfreq.regression.gls import fit_gls  # noqa: E402

PUBLISHED = REPO_ROOT / "flowfreq" / "data" / "regression" / "WA.json"
OUT = REPO_ROOT / "flowfreq" / "data" / "regression" / "provisional" / "WA.json"
REFIT_AEPS = r2.AEPS
REGION = "2"
NAMES = ["intercept", "log10 DRNAREA", "PRECPRIS10", "CANOPY_PCT", "ELEV"]
CITATION = {
    "publication": (
        "flowfreq provisional Washington Peak Region 2 refit with mean basin elevation "
        "(NOT a USGS study), tools/wa_region2_provisional.py; design and calibration gages of "
        "SIR 2016-5118"
    ),
    "table": "docs/WA_REGION2_METHOW.md; data/regression_study/wa_2016_5118_r2_inputs.csv",
    "year": 2026,
}


def _limits(values: np.ndarray) -> Dict[str, float]:
    return {"minimum": float(np.nanmin(values)), "maximum": float(np.nanmax(values))}


def refit_equations(published: Dict[str, Any]) -> List[Dict[str, Any]]:
    """The two refit Region 2 equations, in the library schema."""
    d = r2.load_frame()
    cal = d[d["used"] & d["ELEV"].notna()].reset_index(drop=True)
    X = r2.design(cal, ("ELEV",))
    out = []
    for aep in REFIT_AEPS:
        base = next(
            e for e in published["equations"] if e["region_code"] == REGION and e["aep"] == aep
        )
        y = np.log10(cal[f"q_{aep}"].to_numpy(float))
        res = fit_gls(y, X, r2.sigma(cal, X, aep), names=NAMES)
        variables = [
            {**copy.deepcopy(v), **_limits(cal[col].to_numpy(float))}
            for v, col in zip(base["variables"], ("DRNAREA", "PRECIP", "CANOPY_PCT"))
        ]
        variables.append(
            {"code": "ELEV", "transform": "identity", "units": "ft", **_limits(cal["ELEV"])}
        )
        out.append(
            {
                "region_code": REGION,
                "region_name": base["region_name"],
                "aep": aep,
                "intercept": float(res.beta[0]),
                "variables": variables,
                "coefficients": [float(b) for b in res.beta[1:]],
                "citation": CITATION,
                "sep_log": float(np.sqrt(res.avp)),
                "avp": float(res.avp),
                "model_error_variance": float(res.model_error_variance),
                "covariance": [[float(c) for c in row] for row in res.cov_beta],
                "n_sites": int(len(cal)),
            }
        )
    return out


def build() -> Dict[str, Any]:
    published = json.loads(PUBLISHED.read_text(encoding="utf-8"))
    lib = copy.deepcopy(published)
    refit = {e["aep"]: e for e in refit_equations(published)}
    lib["status"] = "provisional"
    lib["equations"] = [
        refit.get(e["aep"], e) if e["region_code"] == REGION else e for e in lib["equations"]
    ]
    lib["source_reports"] = list(published["source_reports"]) + [
        {"citation": CITATION["publication"], "version": "2026-10-31"}
    ]
    lib["notes"] = (
        "PROVISIONAL (flowfreq, not USGS). Washington's published SIR 2016-5118 set "
        "(load_state('WA')) with one change: Peak Region 2 (NSS GC1751) adds mean basin "
        "elevation, ELEV in feet (identity), refit at every AEP by GLS on the Region 2 "
        "calibration gages with flowfreq EMA quantiles through WY2025 "
        "(tools/wa_region2_provisional.py). ELEV is significant at the 50% AEP (t = 4.1) and "
        "explains about half of the published equation's bias at the Methow; its coefficient "
        "tapers to zero by the 0.2% AEP. Every AEP is refit, not only the common floods, "
        "because mixing refit and published equations gave a decreasing frequency curve at high "
        "basins (docs/WA_REGION2_METHOW.md). Every other region is the published set. "
        "Published notes follow. " + str(published.get("notes", ""))
    )
    return lib


def main() -> int:
    lib = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(lib, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    for e in lib["equations"]:
        if e["region_code"] == REGION and e["aep"] in REFIT_AEPS:
            print(
                f"AEP {e['aep']}: b0 {e['intercept']:.4f}, b {[round(c, 5) for c in e['coefficients']]}, "
                f"MEV {e['model_error_variance']:.4f}, SEP {e['sep_log']:.4f} log, n {e['n_sites']}"
            )
    print(f"wrote {OUT.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
