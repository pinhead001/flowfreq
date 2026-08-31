"""
Test fixtures from peakfqr test-moments.R (Wyoming/Montana multi-site tests).

These tests run the full peakfq() function on the wymt_ffa_2022A test dataset
and compare site info, moments, MGBT results, quantiles, and plotting positions
against PeakFQ v7.4 expected output CSVs.

Since these tests require running the full peakfq pipeline (not individual
Fortran routines), the fixture stores the expected CSV filenames and column
mappings for validation after the Fortran bridge is available.
"""

from __future__ import annotations

# Test dataset identifier
DATASET_NAME: str = "wymt_ffa_2022A"
DATASET_DESCRIPTION: str = (
    "Wyoming/Montana FFA 2022A multi-site test — "
    "validates full peakfq pipeline against PeakFQ v7.4 output"
)

# Input files (relative to peakfqr/inst/testdata/)
INPUT_PSF: str = "wymt_ffa_2022A.psf"
INPUT_DATA: str = "wymt_ffa_2022A_WATSTORE.TXT"

# Expected output CSVs (relative to peakfqr/inst/testdata/)
EXPECTED_INFO_CSV: str = "wymt_ffa_2022A_EXPinfo_7_4.csv"
EXPECTED_DATA_CSV: str = "wymt_ffa_2022A_EXPdata_7_4.csv"
EXPECTED_EMP_CSV: str = "wymt_ffa_2022A_EMPdata_7_4.csv"
EXPECTED_MGBT_CSV: str = "wymt_ffa_2022A_MGBT_7_5_1.csv"

# Test column sets (matching R test expectations)
SITE_INFO_COLS: list[str] = [
    "site_no",
    "station_nm",
    "BegYear",
    "EndYear",
    "HistPeaks",
    "SkewOption",
]

MOMENTS_COLS: list[str] = ["site_no", "Mean", "StandDev", "AtSiteSkew"]
MOMENTS_TOLERANCE: float = 0.03

MGBT_COLS: list[str] = ["site_no", "PILF_Method", "PILF_Thresh", "PILFs", "PILF_0s"]

QUANTILE_COLS: list[str] = [
    "site_no",
    "EXC_Prob",
    "Estimate",
    "Variance",
    "Conf_Low",
    "Conf_Up",
]
QUANTILE_TOLERANCE: float = 0.001
# Note: Quantile test only uses Station skew sites (not Weighted/Generalized)

PLOTTING_POS_COLS: list[str] = ["site_no", "peak_WY", "plot_pos"]
PLOTTING_POS_TOLERANCE: float = 0.001
# Note: Sites 06326960 and 06328100 excluded due to v7.4 rounding bug

# Sites to exclude from specific tests
MGBT_EXCLUDE_SITES: list[str] = ["06326960.00"]  # Code 4 treatment difference
PLOTTING_POS_EXCLUDE_SITES: list[str] = ["06326960.00", "06328100.00"]  # Rounding bug
