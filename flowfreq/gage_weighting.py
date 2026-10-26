"""
Weighting an at-site flood frequency estimate with a regional regression.

Bulletin 17C Appendix 9: at a gage, the at-site log quantile and the
regression's log estimate are combined by the inverse of their variances::

    log Qw = (V_r log Q_s + V_s log Q_r) / (V_s + V_r)
    V_w    = V_s V_r / (V_s + V_r)

``V_s`` is the variance of the at-site estimate, peakfq's ``var_est``
(``FrequencyResults.confidence_limits["variance_log"]``). ``V_r`` is the
regression's variance of prediction at the gage
(:attr:`flowfreq.regression.RegressionEstimate.prediction_variance`).

What this module cannot do
--------------------------
- **The two estimates are treated as independent.** A gage that calibrated
  the regression is correlated with it, which makes ``V_w`` slightly too
  small. State reports (e.g. SIR 2016-5118) make the same assumption.
- **The regression is used as given.** Its basin characteristics are the
  caller's, so a regression evaluated at a different delineation of the gage
  than the report's gives a different weight.
- **Years-of-record weighting is a different method.** Reports that weight by
  record length and equivalent years (Sauer, 1974) are served by
  :func:`flowfreq.regression.nevada.weight_gaged_estimate`.

References
----------
England, J.F., Jr., and others, 2019, Guidelines for determining flood flow
frequency, Bulletin 17C (ver. 1.1): U.S. Geological Survey Techniques and
Methods, book 4, chap. B5, Appendix 9.
"""

from __future__ import annotations

import logging
import math
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd
from scipy.special import ndtri

from flowfreq.core import FrequencyResults
from flowfreq.regression.equations import RegressionEstimate

logger = logging.getLogger(__name__)

#: Tolerance for matching a regression's AEP to an at-site AEP.
_AEP_TOL = 1e-9


def _by_aep(regression: Sequence[RegressionEstimate]) -> Dict[float, RegressionEstimate]:
    out: Dict[float, RegressionEstimate] = {}
    for est in regression:
        aep = float(est.equation.aep)
        if any(abs(aep - a) <= _AEP_TOL for a in out):
            raise ValueError(f"more than one regression estimate at AEP {aep:g}")
        out[aep] = est
    return out


def require_variance_log(results: FrequencyResults) -> pd.DataFrame:
    """Return ``results.confidence_limits``, refusing one without ``variance_log``.

    Raises
    ------
    ValueError
        ``confidence_limits`` is missing or has no ``variance_log`` column (an analysis
        run before the column existed).
    """
    cl = results.confidence_limits
    if cl is None or "variance_log" not in cl.columns:
        raise ValueError(
            "results.confidence_limits has no 'variance_log' column; rerun the analysis "
            "with this version of flowfreq"
        )
    return cl


def weight_with_regression(
    results: FrequencyResults,
    regression: Sequence[RegressionEstimate],
    *,
    level: float = 0.90,
    allow_extrapolated: bool = False,
) -> pd.DataFrame:
    """Weight at-site quantiles with regression estimates (Bulletin 17C Appendix 9).

    Parameters
    ----------
    results : FrequencyResults
        The gage's analysis. Its ``confidence_limits`` must carry
        ``variance_log``.
    regression : sequence of RegressionEstimate
        At most one per AEP, evaluated at the gage.
    level : float, default 0.90
        Two-sided level of the weighted interval.
    allow_extrapolated : bool, default False
        Weight with a regression evaluated outside its calibrated range.
        Without it, those AEPs are left unweighted.

    Returns
    -------
    pd.DataFrame
        One row per at-site AEP: ``aep``, ``flow_cfs`` and ``variance_log``
        (at-site), ``regression_flow_cfs``, ``regression_variance_log``,
        ``weight_at_site``, ``weighted_flow_cfs``, ``weighted_variance_log``,
        ``weighted_lower_cfs``, ``weighted_upper_cfs`` and ``note``. The
        weighted columns are NaN, and ``note`` says why, at an AEP with no
        regression estimate, one without a variance of prediction, or one
        outside its calibrated range.

    Raises
    ------
    ValueError
        ``confidence_limits`` has no ``variance_log`` column, two estimates
        share an AEP, or ``level`` is not in (0, 1).

    Notes
    -----
    The regression enters as ``log_flow``, before any bias correction the
    equation applies to ``flow_cfs``: the weighting is in log space.
    """
    if not 0.0 < level < 1.0:
        raise ValueError(f"level must be in (0, 1), got {level}")
    cl = require_variance_log(results)
    regs = _by_aep(regression)
    z = float(ndtri(0.5 + level / 2.0))

    rows: List[Dict[str, object]] = []
    for aep, q_s, v_s in zip(cl["aep"], cl["flow_cfs"], cl["variance_log"]):
        aep = float(aep)
        row: Dict[str, object] = {
            "aep": aep,
            "flow_cfs": float(q_s),
            "variance_log": float(v_s),
            "regression_flow_cfs": np.nan,
            "regression_variance_log": np.nan,
            "weight_at_site": np.nan,
            "weighted_flow_cfs": np.nan,
            "weighted_variance_log": np.nan,
            "weighted_lower_cfs": np.nan,
            "weighted_upper_cfs": np.nan,
            "note": "",
        }
        rows.append(row)
        est = next((e for a, e in regs.items() if abs(a - aep) <= _AEP_TOL), None)
        if est is None:
            row["note"] = "no regression estimate at this AEP"
            continue
        row["regression_flow_cfs"] = est.flow_cfs
        v_r = est.prediction_variance
        if v_r is None:
            row["note"] = f"regression has no variance of prediction ({est.interval_method})"
            continue
        row["regression_variance_log"] = v_r
        if est.extrapolated and not allow_extrapolated:
            row["note"] = "regression outside its calibrated range: " + "; ".join(est.out_of_range)
            continue
        if not (math.isfinite(est.log_flow) and q_s > 0 and v_s > 0 and v_r > 0):
            row["note"] = "a flow or variance is zero or not finite"
            continue
        w = v_r / (v_s + v_r)
        log_w = w * math.log10(q_s) + (1.0 - w) * est.log_flow
        v_w = v_s * v_r / (v_s + v_r)
        half = z * math.sqrt(v_w)
        row.update(
            weight_at_site=w,
            weighted_flow_cfs=10.0**log_w,
            weighted_variance_log=v_w,
            weighted_lower_cfs=10.0 ** (log_w - half),
            weighted_upper_cfs=10.0 ** (log_w + half),
            note="extrapolated regression" if est.extrapolated else "",
        )
    unused = [a for a in regs if not any(abs(a - float(x)) <= _AEP_TOL for x in cl["aep"])]
    if unused:
        logger.warning(
            "regression AEP(s) %s have no at-site quantile and were not used",
            ", ".join(f"{a:g}" for a in unused),
        )
    return pd.DataFrame(rows)
