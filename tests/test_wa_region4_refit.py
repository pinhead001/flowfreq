"""The Washington GLS refit (tools/wa_region4_refit.py; docs/WA_REGION4_REFIT.md).

Reads the committed results JSON. Two things are pinned:

- Validation of flowfreq.regression.gls on a published USGS regression: the refit
  predicts the published equations' flows at the calibration gages, and its model
  error variance is close to Table 7's.
- The Region 4 decision: Table 7's printed matrix is in the P/10 basis.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

DATA = Path(__file__).resolve().parents[1] / "data" / "regression_study"
AEPS = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)


@pytest.fixture(scope="module")
def refit():
    return json.loads((DATA / "wa_2016_5118_refit.json").read_text(encoding="utf-8"))


def test_calibration_sets_complete(refit):
    m = {r: refit["regions"][r]["m"] for r in ("1", "2", "3", "4")}
    assert m == {"1": 93, "2": 89, "3": 142, "4": 139}
    # Idaho/Oregon gages Table 8 omits are flowfreq EMA estimates.
    ema = {r: refit["regions"][r]["m_flowfreq_ema"] for r in ("1", "2", "3", "4")}
    assert ema == {"1": 22, "2": 2, "3": 0, "4": 27}


@pytest.mark.parametrize("region", ["1", "2", "3", "4"])
def test_model_error_variance_close_to_table_7(refit, region):
    for aep in AEPS:
        e = refit["regions"][region]["basis_P"][str(aep)]
        assert 0.8 <= e["mev"] / e["pub_mev"] <= 1.2, (region, aep)


def _prediction_gap(refit, region: str, aep: float) -> float:
    """Median |log10 Q_refit - log10 Q_published| over the region's calibration gages."""
    sites = pd.read_csv(DATA / "wa_2016_5118_sites.csv", dtype={"site_no": str})
    cal = sites[sites["used"] & (sites["region"] == int(region))]
    la = np.log10(cal["DRNAREA"].to_numpy())
    p = cal["PRECIP"].to_numpy()
    if region in ("1", "2"):
        X = np.column_stack([np.ones(len(cal)), la, p, cal["CANOPY_PCT"].to_numpy()])
    elif region == "3":
        X = np.column_stack([np.ones(len(cal)), la, np.log10(p)])
    else:
        X = np.column_stack([np.ones(len(cal)), la, p])
    e = refit["regions"][region]["basis_P"][str(aep)]
    return float(np.median(np.abs(X @ np.array(e["beta"]) - X @ np.array(e["pub_beta"]))))


@pytest.mark.parametrize("region", ["3", "4"])
def test_refit_reproduces_published_predictions(refit, region):
    # Region 3 is built entirely from printed Table 8 values; Region 4 mostly.
    for aep in AEPS:
        assert _prediction_gap(refit, region, aep) <= 0.045, aep


@pytest.mark.parametrize("region", ["1", "2"])
def test_refit_regions_1_2_close_at_frequent_aeps(refit, region):
    # Regions 1 and 2 drift at rare AEPs (median gap up to 0.12 log units at AEP
    # 0.002 in Region 1, whose 22 Idaho/Oregon gages are flowfreq EMA estimates
    # without Table 2's perception thresholds). At AEP >= 0.04 they agree.
    for aep in (0.5, 0.2, 0.1, 0.04):
        assert _prediction_gap(refit, region, aep) <= 0.03, aep


def test_region_3_covariance_matches_table_7(refit):
    for aep in AEPS:
        ratios = refit["regions"]["3"]["basis_P"][str(aep)]["diag_ratio"]
        assert all(0.67 <= r <= 1.5 for r in ratios), (aep, ratios)


def test_region_4_basis_is_p_over_10(refit):
    r4 = refit["region4"]
    assert r4["selected_basis"] == "basis_P/10"
    assert 0.67 <= r4["selected_P_entry_ratio"] <= 1.5
    scores = r4["median_abs_log_ratio_P_entry"]
    # Raw P is off by about two orders of magnitude in the P variance.
    assert scores["basis_P"] > 4.0 and scores["basis_P/10"] < 0.3
    for aep, row in r4["table6_sp_check"]["basis_P/10"].items():
        assert row["ratio"] == pytest.approx(1.0, abs=0.01), aep
    for aep, row in r4["table6_sp_check"]["basis_P"].items():
        assert row["ratio"] > 2.5, aep
