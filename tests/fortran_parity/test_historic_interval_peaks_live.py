"""Historic interval peaks (code 7 with 4 or 8) against a live ``emafitpr`` call.

``siteQT`` (``vendor/peakfqr/R/readInputs.R``) makes a historic less-than or
greater-than peak an interval row with ``dtype = 1``. ``gbtest`` leaves
``dtype = 1`` rows out of MGBT's sample but still recodes them below the
cutoff and counts them in ``gbnlow``. Two records exercise that, each fitted
by both engines from the same arguments:

* a synthetic record: 50 systematic peaks, a 1900-1959 historic period
  with a 4000 cfs perception threshold, one exact historic flood, a historic
  greater-than peak and a historic less-than peak below the low-outlier
  cutoff, under both MGBT and a fixed threshold;
* the same peaks through a ``.psf``: the peak codes ``"7"``, ``"4,7"`` and
  ``"7,8"`` in a WATSTORE-style frame, converted by ``convert_station`` and
  checked row for row against ``siteQT``'s rows before fitting.

No vendored WATSTORE file has a peak coded both 7 and 4/8 (checked:
``wymt_ffa_2022A_WATSTORE.TXT``, ``HU02_WATSTORE.txt``,
``HU0405_WATSTORE.txt``; 35,759 peaks), so there is no real record here.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip(
    "flowfreq.peakfqr",
    reason="Fortran extension not built; run python build_fortran/build.py "
    "(see docs/FORTRAN_UPLOAD.md)",
    # See tests/fortran_parity/test_fortran_oracles.py's exc_type comment.
    exc_type=ImportError,
)

pytestmark = [pytest.mark.requires_fortran]


def _systematic():
    rng = np.random.default_rng(11)
    years = np.arange(1960, 2010)
    flows = np.round(10 ** rng.normal(3.0, 0.3, len(years)), 0)
    return flows, years


def _kwargs(**extra):
    flows, years = _systematic()
    kw = dict(
        peak_flows=flows,
        water_years=years,
        historical_peaks=[(1925, 9000.0)],
        historical_interval_peaks=[(1938, 6000.0, np.inf), (1950, 0.0, 150.0)],
        perception_thresholds={(1900, 1959): 4000.0, (1960, 2009): 0.0},
        regional_skew=-0.1,
        regional_skew_mse=0.3,
    )
    kw.update(extra)
    return kw


def _compare(kw):
    from flowfreq.bulletin17c import Bulletin17C

    native = Bulletin17C(**kw).run_analysis(engine="native")
    fortran = Bulletin17C(**kw).run_analysis(engine="fortran")
    return native, fortran


def _assert_parity(native, fortran):
    assert native.n_peaks == fortran.n_peaks
    assert native.n_historical == fortran.n_historical == 3
    assert native.n_low_outliers == fortran.n_low_outliers
    assert native.n_mgbt_outliers == fortran.n_mgbt_outliers
    assert native.low_outlier_threshold == pytest.approx(fortran.low_outlier_threshold)
    assert abs(native.mean_log - fortran.mean_log) < 1e-6
    assert abs(native.std_log - fortran.std_log) < 1e-6
    assert abs(native.skew_station - fortran.skew_station) < 1e-5
    assert abs(native.skew_used - fortran.skew_used) < 1e-5
    qn = native.quantiles.set_index("aep")["flow_cfs"]
    qf = fortran.quantiles.set_index("aep")["flow_cfs"]
    assert float((qn / qf - 1).abs().max()) < 1e-5
    cn = native.confidence_limits.set_index("aep")
    cf = fortran.confidence_limits.set_index("aep")
    for col in ("lower_5pct", "upper_5pct"):
        assert float((cn[col] / cf[col] - 1).abs().max()) < 1e-3


def test_fixed_threshold():
    """The 150 cfs historic less-than peak is below a 300 cfs cutoff: recoded, counted."""
    flows, _ = _systematic()
    native, fortran = _compare(_kwargs(user_low_outlier_threshold=300.0))
    _assert_parity(native, fortran)
    below = int(np.sum(flows < 300.0))
    assert native.n_mgbt_outliers == below
    assert native.n_low_outliers == below + 1


def test_mgbt():
    native, fortran = _compare(_kwargs())
    _assert_parity(native, fortran)


def test_psf_codes_build_sitetq_rows_and_fit():
    from flowfreq.psf import parse_psf
    from flowfreq.psf_convert import convert_station

    flows, years = _systematic()
    frame = pd.DataFrame(
        {
            "water_year": np.concatenate([years, [1925, 1938, 1950]]),
            "peak_flow_cfs": np.concatenate([flows, [9000.0, 6000.0, 150.0]]),
            "qualification_code": [""] * len(years) + ["7", "7,8", "4,7"],
        }
    )
    spec = parse_psf(
        "Station S1\n"
        "PCPT_Thresh 1900 1959 4000 1E+20\nPCPT_Thresh 1960 2009 0 1E+20\n"
        "SkewOpt Weighted\nGenSkew -0.1\nSkewSE 0.547722557505\nLOType MGBT\n"
    ).stations["S1"]
    station = convert_station(spec, frame)
    assert station.historical_interval_peaks == ((1938, 6000.0, 1e20), (1950, 1e-20, 150.0))
    for engine in ("native", "fortran"):
        assert station.unsupported_reasons(engine) == []

    native = station.run("native")
    reference = station.fortran_reference()
    p = reference.parameters
    assert native.n_low_outliers == reference.low_outlier_count
    assert abs(native.mean_log - p["mean_log"]) < 1e-6
    assert abs(native.skew_used - p["skew_weighted"]) < 1e-5
    for aep, flow in reference.quantiles.items():
        q = native.quantiles
        got = float(q.loc[np.isclose(q["aep"], aep), "flow_cfs"].iloc[0])
        assert got == pytest.approx(flow, rel=1e-5), aep
