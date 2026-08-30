"""
Big Sandy River at Bruceton, TN (USGS gage 03606500)
Test fixture data sourced from PeakfqSA User Manual (Tim Cohn, USGS, 2012).
Used as the primary validation case for hybrid Bulletin 17C implementation.
"""

# Systematic annual peaks 1930-1973 (cfs)
SYSTEMATIC_PEAKS: dict[int, float] = {
    1930: 9100,
    1931: 2060,
    1932: 7820,
    1933: 3220,
    1934: 5580,
    1935: 17000,
    1936: 6740,
    1937: 13800,
    1938: 4270,
    1939: 5940,
    1940: 1680,
    1941: 1200,
    1942: 10100,
    1943: 3780,
    1944: 5340,
    1945: 5630,
    1946: 12000,
    1947: 3980,
    1948: 6130,
    1949: 4740,
    1950: 9880,
    1951: 5230,
    1952: 4260,
    1953: 5000,
    1954: 3320,
    1955: 5480,
    1956: 11800,
    1957: 5150,
    1958: 3350,
    1959: 2400,
    1960: 1460,
    1961: 3770,
    1962: 7480,
    1963: 2740,
    1964: 3100,
    1965: 7180,
    1966: 1920,
    1967: 9060,
    1968: 3080,
    1969: 2800,
    1970: 4330,
    1971: 5080,
    1972: 12000,
    1973: 7640,
}

# Historical floods (known to exceed 18,000 cfs threshold)
HISTORICAL_PEAKS: dict[int, float] = {
    1897: 25000,
    1919: 21000,
    1927: 18500,
}

# Perception thresholds
THRESHOLDS: list[dict[str, float]] = [
    {"start": 1890, "end": 1929, "lower": 18000.0, "upper": 1e50},
    {"start": 1930, "end": 1973, "lower": 0.0, "upper": 1e50},
]

BEGYEAR: int = 1890
ENDYEAR: int = 1973
REGIONAL_SKEW: float = -0.5
REGIONAL_SKEW_SD: float = 0.55
STATION_NAME: str = "BIG SANDY RIVER AT BRUCETON, TN, 1890-1973"

# Expected results from PeakfqSA manual (page 26-27)
EXPECTED_PARAMETERS: dict[str, float] = {
    "mean_log": 3.717272,
    "std_log": 0.289200,
    "skew_weighted": -0.118702,
}

EXPECTED_QUANTILES: dict[float, float] = {
    # AEP: discharge (cfs) — from PeakfqSA manual output
    0.9950: 871.25,
    0.9900: 1045.59,
    0.9500: 1706.18,
    0.9000: 2203.77,
    0.8000: 2990.15,
    0.6667: 3957.50,
    0.5000: 5284.36,
    0.2000: 9166.15,
    0.1000: 12134.65,
    0.0400: 16276.60,
    0.0200: 19617.73,
    0.0100: 23158.65,
    0.0050: 26912.12,
    0.0020: 32217.14,
}

EXPECTED_CONFIDENCE_INTERVALS: dict[float, tuple[float, float]] = {
    # AEP: (lower_95, upper_95)
    0.1000: (9766.00, 15218.32),
    0.0200: (15154.99, 29124.18),
    0.0100: (17388.03, 37986.08),
}

TOLERANCE_PERCENT: float = 1.0  # Results must match within 1%


# ---------------------------------------------------------------------------
# peakfq 8.1.0 reference values
# ---------------------------------------------------------------------------
# Regenerated from the vendored USGS Fortran (vendor/peakfqr, peakfq 8.1.0) via
# tools/gen_fortran_golden.py. Use THESE for parity work.
#
# The EXPECTED_* values above come from the 2012 PeakfqSA manual and are kept as
# historical record only -- they are NOT reproducible by peakfq 8.1.0. The skew
# weighting changed: 8.1.0 uses HWN, "a generalization of the PeakFQ 7.4.1
# algorithm using an optimized adjustment factor when censored data are present"
# (vendor/peakfqr/R/fortranWrappers.R), and Big Sandy has 37 censored intervals.
# Feeding 8.1.0's own reported at-site skew MSE (0.09437) through the standard
# Bulletin 17C weighting reproduces the manual to within 4%; 8.1.0's internal
# weighting differs by 32%. See docs/FORTRAN_UPLOAD.md section 6.0b.

PEAKFQ_810_PARAMETERS: dict[str, float] = {
    "mean_log": 3.717508,
    "std_log": 0.291043,
    "skew_weighted": -0.156306,
    "skew_at_site": 0.006601,
}

PEAKFQ_810_QUANTILES: dict[float, float] = {
    0.995: 841.75,
    0.99: 1016.70,
    0.95: 1683.48,
    0.9: 2187.72,
    0.8: 2985.60,
    0.6667: 3966.78,
    0.5: 5309.91,
    0.2: 9213.05,
    0.1: 12169.84,
    0.04: 16257.05,
    0.02: 19523.97,
    0.01: 22959.39,
    0.005: 26573.14,
    0.002: 31635.64,
}

# 90% interval: the 5th and 95th percentiles (eps=0.90)
PEAKFQ_810_CONFIDENCE_INTERVALS: dict[float, tuple[float, float]] = {
    0.1: (10170.10, 14646.38),
    0.02: (15650.11, 26088.55),
    0.01: (17903.45, 32590.70),
}
