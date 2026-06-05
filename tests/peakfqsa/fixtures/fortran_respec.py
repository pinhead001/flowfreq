"""
Test fixtures from peakfqr test-fortran.R (Respec EMA tests).

These test the low-level Fortran routines: moms_p3, p3est_ema, var_mom, and qP3.
Data is provided in real space (ql/qu) as used in the R tests, which operate
on real-valued flow intervals (not log10).
"""

from __future__ import annotations

# --- Input data for moms_p3 and p3est_ema tests ---
# 34 flow intervals (ql, qu) from the Respec test suite
QL: list[float] = [
    68.0,
    20.0,
    345.0,
    139.0,
    680.0,
    970.0,
    355.0,
    76.0,
    110.0,
    12.0,
    117.0,
    102.0,
    74.0,
    82.0,
    221.0,
    84.0,
    55.0,
    117.0,
    225.0,
    94.0,
    0.0,
    320.0,
    0.0,
    0.0,
    41.0,
    27.0,
    205.0,
    174.0,
    100.0,
    108.0,
    387.0,
    29.0,
    111.0,
    178.0,
]

QU: list[float] = [
    68.0,
    20.0,
    345.0,
    139.0,
    680.0,
    970.0,
    355.0,
    76.0,
    110.0,
    12.0,
    117.0,
    102.0,
    74.0,
    82.0,
    221.0,
    84.0,
    55.0,
    117.0,
    225.0,
    94.0,
    40.0,
    320.0,
    40.0,
    40.0,
    41.0,
    27.0,
    205.0,
    174.0,
    100.0,
    108.0,
    387.0,
    29.0,
    111.0,
    178.0,
]

# --- "truth" data: 200 standardized flow intervals (used for truth_test_*) ---
QL_TRUTH: list[float] = [
    -0.69,
    2.04,
    0.36,
    0.17,
    -0.16,
    0.59,
    -0.84,
    -0.63,
    1.03,
    -0.35,
    0.86,
    0.69,
    0.81,
    -0.57,
    -0.47,
    1.24,
    0.29,
    -0.33,
    0.22,
    0.39,
    0.08,
    -0.12,
    1.13,
    -0.25,
    0.67,
    -0.28,
    0.17,
    -0.23,
    -0.59,
    -0.7,
    -0.45,
    -0.26,
    -0.56,
    -0.47,
    -0.26,
    0.26,
    0.35,
    0.64,
    0.11,
    -0.54,
    0,
    0.22,
    0.03,
    -0.12,
    -0.06,
    -0.08,
    0.51,
    -0.58,
    0.17,
    0.01,
    -0.76,
    -1.1,
    0.25,
    0.13,
    -0.63,
    -0.53,
    1.14,
    1.05,
    -0.42,
    0.36,
    0.52,
    -0.54,
    0.27,
    0.86,
    -0.3,
    -0.37,
    0.55,
    0.24,
    -0.05,
    -0.38,
    -0.04,
    0.55,
    -0.62,
    0.26,
    -0.28,
    -0.15,
    -0.13,
    -0.93,
    0.46,
    0.52,
    -0.34,
    -0.15,
    -0.08,
    0.04,
    0.42,
    0.63,
    0.66,
    -0.33,
    0.05,
    -0.14,
    -0.34,
    0.76,
    0.36,
    0.93,
    0.27,
    0.08,
    0.22,
    0.93,
    -0.67,
    0.3,
    0.08,
    0.65,
    -0.28,
    0.38,
    -0.05,
    0.17,
    0.65,
    -0.07,
    -0.29,
    -0.33,
    -0.29,
    0.4,
    0.61,
    -0.42,
    -0.14,
    -0.65,
    0.41,
    0.4,
    0.19,
    -0.67,
    0.2,
    -0.41,
    -0.28,
    1.12,
    1.3,
    -0.61,
    0.44,
    1.07,
    0.14,
    0.45,
    -0.77,
    -0.18,
    -0.03,
    -0.05,
    -0.31,
    -0.54,
    0.45,
    0.54,
    0.8,
    -0.24,
    0.73,
    0.45,
    -0.63,
    -0.33,
    0.09,
    -0.25,
    0.37,
    0.87,
    0.43,
    -0.86,
    0.18,
    -0.3,
    -0.56,
    -0.15,
    -0.73,
    0.41,
    -0.47,
    0.3,
    0.34,
    1.4,
    -0.95,
    -0.49,
    1.27,
    0.34,
    0.2,
    -0.74,
    0.36,
    0.57,
    -0.15,
    -0.53,
    0.73,
    0.27,
    0.06,
    0.34,
    1.46,
    -0.48,
    1.26,
    0.33,
    -0.08,
    -0.45,
    0.52,
    -0.84,
    0.55,
    -0.25,
    -0.36,
    0.02,
    -0.15,
    -0.53,
    0.46,
    0,
    0.38,
    -0.57,
    0.23,
    -0.59,
    0.34,
    -0.76,
    0.09,
    -0.3,
    0.24,
    0.92,
]

QU_TRUTH: list[float] = list(QL_TRUTH)  # ql == qu for truth data (exact obs)


# ========================
# test_moms_p3 expected
# ========================
MOMS_P3_EXPECTED: dict[str, list[float]] = {
    # weight_option: [mean, variance, skew]
    "orig": [165.5409898, 40060.40734, 1.423055596],
    "ERL": [165.5409898, 40060.40734, 1.385585543],
    "HWN": [165.5409898, 40060.40734, 1.449544227],
}

MOMS_P3_INPUTS: dict[str, dict[str, float]] = {
    "orig": {"n": 34, "Wd": 1.0, "nG_factor": 34 * 1.0 * 0.160 / 0.206},
    "ERL": {"n": 34, "Wd": 1.056, "nG_factor": 34 * 1.056 * 0.160 / 0.206},
    "HWN": {"n": 34, "Wd": 0.962, "nG_factor": 34 * 0.962 * 0.160 / 0.206},
}

MOMS_P3_RG: float = -0.145
MOMS_P3_MC_OLD: list[float] = [0.0, 1.0, 0.0]


# ========================
# test_p3est_ema expected
# ========================
P3EST_EMA_EXPECTED: dict[str, list[float]] = {
    # weight_option (Wd value): [mean, variance, skew]
    "orig": [167.236599216544, 39541.5780590362, 1.46729880561388],
    "ERL": [167.238920642166, 39540.8926394831, 1.42946503056195],
    "HWN": [167.234871013990, 39542.0884205924, 1.49402326308862],
}

P3EST_EMA_INPUTS: dict[str, float] = {
    "r_G": -0.145,
    "as_M_mse": 0.006306109073542849,
    "as_S2_mse": 0.002770797581276051,
    "as_G_mse": 0.155132529907670,
    "r_M": 0.0,
    "r_S2": 1.0,
    "r_M_mse": -99.0,
    "r_S2_mse": -3996.0,
    "r_G_mse": 0.206115990877151,
}

P3EST_EMA_WD: dict[str, float] = {
    "orig": 1.0,
    "ERL": 1.056,
    "HWN": 0.962,
}


# ========================
# truth_test_moms_p3
# ========================
TRUTH_MOMS_P3_EXPECTED: list[float] = [0.07415, 0.307708822, 0.452403784]
TRUTH_MOMS_P3_NG: float = 0.0
TRUTH_MOMS_P3_RG: float = 0.0
TRUTH_MOMS_P3_MC_OLD: list[float] = [0.0, 1.0, 0.0]


# ========================
# truth_test_var_mom
# ========================
TRUTH_VAR_MOM_MC: list[float] = [0.07415, 0.307708822, 0.452403784]
TRUTH_VAR_MOM_N: int = 200
TRUTH_VAR_MOM_NTHRESH: int = 1
TRUTH_VAR_MOM_TL: float = -20.0
TRUTH_VAR_MOM_TU: float = 20.0

# Expected diagonal values of variance-covariance matrix
# var_m_x = [mc[2]/n, 2*mc[2]^2/(n-1), formula_for_skew_var]
# Note: skew variance test is commented out in peakfqr (known issue)


# ========================
# truth_test_qP3
# ========================
QP3_EXCEEDANCE_PROBS: list[float] = [0.9980, 0.9900, 0.9000, 0.5000, 0.1000, 0.0100, 0.0020]
QP3_MU: float = 0.07415
QP3_S2: float = 0.307708822
# Tested with g=+2 and g=-2; the test verifies round-trip (quantile -> CDF -> prob)
