"""
HU02 (Northeast US) test fixtures from peakfqr inst/testdata/extra_tests/HU02.

Reference: PeakfqSA v7.5.1 expected outputs.
Selected representative stations covering different EMA scenarios.
"""

from typing import Any

# ---------------------------------------------------------------------------
# Standard exceedance probabilities
# ---------------------------------------------------------------------------
EXCEEDANCE_PROBS: list[float] = [
    0.995,
    0.99,
    0.98,
    0.975,
    0.96,
    0.95,
    0.9,
    0.8,
    0.7,
    0.6667,
    0.6,
    0.5704,
    0.5,
    0.4292,
    0.4,
    0.3,
    0.2,
    0.1,
    0.05,
    0.04,
    0.025,
    0.02,
    0.01,
    0.005,
    0.002,
]


# ---------------------------------------------------------------------------
# Station: 01199050 - Salmon Creek at Lime Rock, CT
# EMA, station skew, historical period with 2 hist peaks, 65-year record
# ---------------------------------------------------------------------------
STATION_01199050: dict[str, Any] = {
    "site_no": "01199050",
    "name": "Salmon Creek at Lime Rock, CT",
    "analysis": "EMA",
    "beg_year": 1949,
    "end_year": 2013,
    "hist_length": 65,
    "skew_option": "Station",
    "expected_params": {
        "skew": 0.548,
        "mean": 2.844,
        "std_dev": 0.301,
        "at_site_skew": 0.548,
        "at_site_mseg": 0.133,
        "at_site_mseg_gaged_only": 0.139,
        "reg_skew": -999,
        "reg_mseg": 0,
        "gaged_peaks": 52,
        "hist_peaks": 2,
        "pilf_method": "MGBT",
        "pilf_thresh": 0,
        "pilfs": 0,
        "ema_num_iter": 24,
    },
    "psf_config": {
        "pcpt_thresh": [
            {"beg": 1949, "end": 2013, "low": 0, "high": 1e20, "comment": "DEFAULT"},
            {"beg": 1950, "end": 1954, "low": 6300, "high": 1e20, "comment": "HISTORIC 1"},
            {"beg": 1956, "end": 1961, "low": 6300, "high": 1e20, "comment": "HISTORIC 1"},
        ],
        "peaks": [
            {"year": 1949, "value": 2100, "type": "H"},
            {"year": 1955, "value": 6300, "type": "H"},
        ],
        "skew_opt": "Station",
        "lo_type": "MGBT",
    },
    "expected_quantiles": {
        0.995: 166.6,
        0.99: 184.3,
        0.98: 207.2,
        0.975: 216.2,
        0.96: 238.3,
        0.95: 250.8,
        0.9: 301.8,
        0.8: 385.2,
        0.7: 466.0,
        0.6667: 494.0,
        0.6: 553.6,
        0.5704: 581.9,
        0.5: 655.6,
        0.4292: 742.2,
        0.4: 782.9,
        0.3: 955.2,
        0.2: 1220,
        0.1: 1753,
        0.05: 2411,
        0.04: 2655,
        0.025: 3227,
        0.02: 3529,
        0.01: 4609,
        0.005: 5941,
        0.002: 8180,
    },
    "expected_ci": {
        0.5: (549.9, 780.7),
        0.1: (1399, 2412),
        0.05: (1840, 3781),
        0.02: (2510, 6996),
        0.01: (3094, 11290),
        0.002: (4746, 35360),
    },
}


# ---------------------------------------------------------------------------
# Station: 01200000 - Tenmile River near Gaylordsville, CT
# EMA, station skew, 84-year record with missing years, no low outliers
# ---------------------------------------------------------------------------
STATION_01200000: dict[str, Any] = {
    "site_no": "01200000",
    "name": "Tenmile River near Gaylordsville, CT",
    "analysis": "EMA",
    "beg_year": 1930,
    "end_year": 2013,
    "hist_length": 84,
    "skew_option": "Station",
    "expected_params": {
        "skew": 0.389,
        "mean": 3.508,
        "std_dev": 0.261,
        "at_site_skew": 0.389,
        "at_site_mseg": 0.087,
        "at_site_mseg_gaged_only": 0.087,
        "reg_skew": -999,
        "reg_mseg": 0,
        "gaged_peaks": 81,
        "hist_peaks": 0,
        "pilf_method": "MGBT",
        "pilf_thresh": 0,
        "pilfs": 0,
        "ema_num_iter": 4,
    },
    "psf_config": {
        "pcpt_thresh": [
            {"beg": 1930, "end": 2013, "low": 0, "high": 1e20, "comment": "DEFAULT"},
            {"beg": 1989, "end": 1991, "low": 1e20, "high": 1e20, "comment": "MISSING1"},
        ],
        "skew_opt": "Station",
        "lo_type": "MGBT",
    },
    "expected_quantiles": {
        0.995: 852.2,
        0.99: 945.9,
        0.98: 1065,
        0.975: 1111,
        0.96: 1223,
        0.95: 1285,
        0.9: 1534,
        0.8: 1925,
        0.7: 2288,
        0.6667: 2411,
        0.6: 2667,
        0.5704: 2788,
        0.5: 3095,
        0.4292: 3447,
        0.4: 3609,
        0.3: 4278,
        0.2: 5258,
        0.1: 7096,
        0.05: 9200,
        0.04: 9943,
        0.025: 11630,
        0.02: 12490,
        0.01: 15430,
        0.005: 18820,
        0.002: 24140,
    },
    "expected_ci": {
        0.5: (2742, 3494),
        0.1: (6073, 8785),
        0.05: (7632, 12360),
        0.02: (9862, 19200),
        0.01: (11690, 26640),
        0.002: (16500, 56000),
    },
}


# ---------------------------------------------------------------------------
# Convenience collection
# ---------------------------------------------------------------------------
ALL_HU02_STATIONS: dict[str, dict[str, Any]] = {
    "01199050": STATION_01199050,
    "01200000": STATION_01200000,
}
