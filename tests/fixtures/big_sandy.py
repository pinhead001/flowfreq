"""
Big Sandy River at Bruceton, TN (USGS gage 03606500).

The primary validation case. The numbers are not defined here: they load from
``flowfreq/validation/data/big_sandy_03606500.json``, which ships with the
package because ``flowfreq benchmark`` needs them and ``tests/`` is not
packaged. Duplicating them here would give the benchmark and the test suite two
sources of truth for the same site, which is how they drift.

Two sets of expectations, deliberately:

``PEAKFQ_810_*``
    Regenerated from the vendored USGS Fortran (peakfq 8.1.0) by
    ``tools/gen_fortran_golden.py``. Use these for parity work.

``EXPECTED_*``
    From the 2012 PeakfqSA User Manual (Cohn, USGS), pages 26-27, kept as
    historical record. They are **not reproducible** by peakfq 8.1.0: its skew
    weighting is HWN, "a generalization of the PeakFQ 7.4.1 algorithm using an
    optimized adjustment factor when censored data are present"
    (``vendor/peakfqr/R/fortranWrappers.R``), and Big Sandy has 37 censored
    intervals. Feeding 8.1.0's own reported at-site skew MSE (0.09437) through
    the standard Bulletin 17C weighting reproduces the manual to within 4%;
    8.1.0's internal weighting differs by 32%. See ``docs/FORTRAN_UPLOAD.md``
    section 6.0b.
"""

from __future__ import annotations

from flowfreq.validation.benchmarks import load_case

_CASE = load_case("big_sandy_03606500")
_REF = _CASE["reference"]
_MANUAL = _CASE["historical_reference"]

# Systematic annual peaks 1930-1973 (cfs)
SYSTEMATIC_PEAKS: dict[int, float] = {
    int(y): float(q) for y, q in _CASE["systematic_peaks"].items()
}

# Historical floods (known to exceed the 18,000 cfs threshold)
HISTORICAL_PEAKS: dict[int, float] = {
    int(y): float(q) for y, q in _CASE["historical_peaks"].items()
}

# Perception thresholds
THRESHOLDS: list[dict[str, float]] = [dict(t) for t in _CASE["thresholds"]]

BEGYEAR: int = int(_CASE["begyear"])
ENDYEAR: int = int(_CASE["endyear"])
REGIONAL_SKEW: float = float(_CASE["regional_skew"])
REGIONAL_SKEW_SD: float = float(_CASE["regional_skew_sd"])
STATION_NAME: str = _CASE["station_name"]

# --- 2012 PeakfqSA manual, historical record only (see the module docstring) --
EXPECTED_PARAMETERS: dict[str, float] = dict(_MANUAL["parameters"])
EXPECTED_QUANTILES: dict[float, float] = {
    float(a): float(q) for a, q in _MANUAL["quantiles"].items()
}
EXPECTED_CONFIDENCE_INTERVALS: dict[float, tuple[float, float]] = {
    float(a): (float(lo), float(hi)) for a, (lo, hi) in _MANUAL["confidence_intervals"].items()
}

TOLERANCE_PERCENT: float = 1.0  # Results must match within 1%

# --- peakfq 8.1.0, the reference for parity work ------------------------------
PEAKFQ_810_PARAMETERS: dict[str, float] = dict(_REF["parameters"])
PEAKFQ_810_QUANTILES: dict[float, float] = {
    float(a): float(q) for a, q in _REF["quantiles"].items()
}
# 90% interval: the 5th and 95th percentiles (eps=0.90)
PEAKFQ_810_CONFIDENCE_INTERVALS: dict[float, tuple[float, float]] = {
    float(a): (float(lo), float(hi)) for a, (lo, hi) in _REF["confidence_intervals"].items()
}
