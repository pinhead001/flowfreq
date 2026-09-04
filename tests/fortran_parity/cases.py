"""Input construction for Fortran parity cases.

Single source of truth: the golden-file generator and both parity tests build
their inputs from here, so a golden file can never drift from the inputs the
tests believe produced it.

Interval construction follows ``siteQT`` in ``vendor/peakfqr/R/readInputs.R``:

* every water year in ``[min(threshold start), max(threshold end)]`` gets a row,
  including years with no observed peak;
* an observed peak is exact, ``ql == qu``;
* a year with no peak is censored against the *lower* perception threshold,
  ``ql = Qmin`` and ``qu = tl``;
* ``dtype`` is 1 only for peaks carrying the USGS historic flag (peak code 7),
  not for every year inside the historical period.

All flows are in real space here; ``build_emafit_inputs`` returns log10, which is
what ``emafitpr`` expects.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Tuple

import numpy as np

QMIN = 1e-20  # peakfqr/R/main.R
QMAX = 1e20


@dataclass(frozen=True)
class ParityCase:
    """One site's worth of EMA inputs, in real space."""

    site_no: str
    description: str
    systematic: Dict[int, float]
    historical: Dict[int, float]
    threshold_start: int
    threshold_end: int
    threshold_lower: float
    regional_skew: float
    regional_skew_mse: float
    aeps: Tuple[float, ...]
    eps: float = 0.90
    weight_opt: int = 1  # 1=HWN, 2=ERL, 3=INV
    gbthrsh0: float = -99.0  # <= -6 runs MGBT
    extra: Dict[str, Any] = field(default_factory=dict)


def build_emafit_inputs(case: ParityCase) -> Dict[str, Any]:
    """Expand a case into the exact argument set ``emafitpr`` takes (log10)."""
    ql: List[float] = []
    qu: List[float] = []
    tl: List[float] = []
    tu: List[float] = []
    dtype: List[int] = []

    for _year, flow in sorted(case.systematic.items()):
        ql.append(flow)
        qu.append(flow)
        tl.append(QMIN)
        tu.append(QMAX)
        dtype.append(0)

    for _year, flow in sorted(case.historical.items()):
        ql.append(flow)
        qu.append(flow)
        tl.append(case.threshold_lower)
        tu.append(QMAX)
        dtype.append(1)  # carries the historic peak flag

    for year in range(case.threshold_start, case.threshold_end + 1):
        if year in case.historical or year in case.systematic:
            continue
        ql.append(QMIN)
        qu.append(case.threshold_lower)  # censored against the lower threshold
        tl.append(case.threshold_lower)
        tu.append(QMAX)
        dtype.append(0)  # no historic flag -> 0, per siteQT

    log = np.log10
    return {
        "ql": log(np.asarray(ql, dtype=float)),
        "qu": log(np.asarray(qu, dtype=float)),
        "tl": log(np.asarray(tl, dtype=float)),
        "tu": log(np.asarray(tu, dtype=float)),
        "dtype": np.asarray(dtype, dtype=np.int32),
        "reg_m": 0.0,
        "reg_m_mse": -1e99,
        "reg_sd": 1.0,
        "reg_sd_mse": 1e99,
        "r_g": case.regional_skew,
        "r_g_mse": case.regional_skew_mse,
        "gbthrsh0": case.gbthrsh0,
        "pq": 1.0 - np.asarray(case.aeps, dtype=float),
        "eps": case.eps,
        "wght_opt_n": case.weight_opt,
    }


def call_emafitpr(case: ParityCase) -> Dict[str, Any]:
    """Run the vendored Fortran for *case* and return its outputs verbatim (log10)."""
    from flowfreq.peakfqr import emafitpr

    a = build_emafit_inputs(case)
    out = emafitpr(
        a["ql"],
        a["qu"],
        a["tl"],
        a["tu"],
        a["dtype"],
        a["reg_m"],
        a["reg_m_mse"],
        a["reg_sd"],
        a["reg_sd_mse"],
        a["r_g"],
        a["r_g_mse"],
        a["gbthrsh0"],
        a["pq"],
        a["eps"],
        a["wght_opt_n"],
    )
    (
        gbval,
        gbns,
        gbnzero,
        gbnlow,
        _gbp,
        _gbqs,
        as_g_mse,
        as_g_mse_syst,
        as_g_prl,
        cmoms,
        yp,
        ci_low,
        ci_high,
        var_est,
        wdout,
        *_rest,
    ) = out
    n = len(a["ql"])
    return {
        "n": n,
        "mgbt": {
            "gbval": float(gbval),
            "gbns": int(gbns),
            "gbnzero": int(gbnzero),
            "gbnlow": int(gbnlow),
        },
        "skew": {
            "as_G_mse_o": float(as_g_mse),
            "as_G_mse_Syst_o": float(as_g_mse_syst),
            "as_G_PRL_o": float(as_g_prl),
            "Wdout": float(wdout),
        },
        "cmoms": np.asarray(cmoms, dtype=float).tolist(),
        "quantiles": {
            "yp": np.asarray(yp).tolist(),
            "ci_low": np.asarray(ci_low).tolist(),
            "ci_high": np.asarray(ci_high).tolist(),
            "var_est": np.asarray(var_est).tolist(),
        },
    }


def big_sandy_case() -> ParityCase:
    """USGS 03606500, Big Sandy River at Bruceton TN — the primary validation site."""
    from tests.fixtures.big_sandy import (
        HISTORICAL_PEAKS,
        REGIONAL_SKEW,
        REGIONAL_SKEW_SD,
        STATION_NAME,
        SYSTEMATIC_PEAKS,
        THRESHOLDS,
    )

    historical_threshold = min(t["lower"] for t in THRESHOLDS if t["lower"] > 0)
    return ParityCase(
        site_no="03606500",
        description=STATION_NAME,
        systematic=dict(SYSTEMATIC_PEAKS),
        historical=dict(HISTORICAL_PEAKS),
        threshold_start=min(t["start"] for t in THRESHOLDS),
        threshold_end=max(t["end"] for t in THRESHOLDS),
        threshold_lower=historical_threshold,
        regional_skew=REGIONAL_SKEW,
        regional_skew_mse=REGIONAL_SKEW_SD**2,
        aeps=tuple(
            sorted(
                (
                    0.995,
                    0.99,
                    0.95,
                    0.9,
                    0.8,
                    0.6667,
                    0.5,
                    0.2,
                    0.1,
                    0.04,
                    0.02,
                    0.01,
                    0.005,
                    0.002,
                ),
                reverse=True,
            )
        ),
    )


# Standard AEP set, shared by every case so goldens stay comparable.
_AEPS = tuple(
    sorted(
        (0.995, 0.99, 0.95, 0.9, 0.8, 0.6667, 0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002),
        reverse=True,
    )
)


def wymt_case(site_no: str) -> ParityCase:
    """A Wyoming/Montana site from the vendored peakfqr test data.

    These sites are plain systematic records -- contiguous water years, no
    historic peaks, no zero flows -- so their EMA inputs are just the peaks,
    with the systematic perception thresholds ``(Qmin, Qmax)`` throughout.
    That keeps the case construction unambiguous, which is the point: they
    exist to test the parity machinery against a *second* site, not to
    exercise censoring reconstruction.

    Two things Big Sandy cannot test that these can. Its at-site skew is
    0.0066, below the 0.04 floor at ``emafit.f:763``, so the Halloween
    determinant ratio is never reached there; both sites here are well above
    it. And its conditioning is the only conditioning the live-vs-golden
    tolerances have ever been calibrated against.

    Parameters
    ----------
    site_no : str
        Station number as it appears in the peakfqr CSVs.

    Returns
    -------
    ParityCase

    Raises
    ------
    ValueError
        If the site is not the plain systematic record this assumes.
    """
    from tests.fixtures.wymt_peaks import load_site

    site = load_site(site_no)
    if site.n_historical or not site.is_contiguous:
        raise ValueError(
            f"site {site_no} is not a contiguous systematic record "
            f"(historic peaks: {site.n_historical}, contiguous: {site.is_contiguous}); "
            "it needs perception thresholds built explicitly"
        )
    if site.regional_skew is None:
        raise ValueError(f"site {site_no} has no regional skew, so weighting is not exercised")

    years = sorted(site.peaks)
    return ParityCase(
        site_no=site.site_no,
        description=f"{site.name}, {years[0]}-{years[-1]}",
        systematic=dict(site.peaks),
        historical={},
        threshold_start=years[0],
        threshold_end=years[-1],
        # Unused: with contiguous years and no historic peaks the censored-fill
        # loop in build_emafit_inputs adds no rows, so no year is ever measured
        # against this. QMIN keeps it harmless if that ever changes.
        threshold_lower=QMIN,
        regional_skew=site.regional_skew,
        regional_skew_mse=site.regional_skew_mse,
        aeps=_AEPS,
    )


CASES = {
    "big_sandy_03606500": big_sandy_case,
    # Powder River near Locate MT: 85 contiguous years, no PILFs. The plain
    # counterpart to Big Sandy -- same machinery, different conditioning.
    "powder_river_06326500": lambda: wymt_case("06326500.00"),
    # Cains Coulee at Glendive MT: 32 years and 11 PILFs under MGBT, so the
    # censoring here is produced by the fit rather than supplied as input.
    "cains_coulee_06327450": lambda: wymt_case("06327450.00"),
}
