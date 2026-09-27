"""The .psf converter against live ``emafitpr`` calls, with the ``.psf``'s own values.

For each WY/MT station, :meth:`StationInputs.fortran_reference` runs the
vendored Fortran on the ``siteQT`` rows with ``main.R``'s arguments -- what
peakfq 8.1.0 itself does with this ``.psf``. Two claims are checked against
it:

* Every station ``Bulletin17C``'s Fortran engine can express (all but
  06328100, whose 2021 upper threshold of 407 cfs has no ``Bulletin17C``
  argument) gets bit-identical moments and quantiles through
  ``Bulletin17C(**kwargs).run_analysis(engine="fortran")``. That pins the
  argument mapping -- including the station-skew sentinel, which
  ``fortran_engine`` sends as ``r_G_mse = 1e10`` where ``main.R`` sends
  ``-1e99`` -- rather than just the rows.
* The native engine on Powder River, uncensored, matches it to 1e-6 with the
  ``.psf``'s unrounded regional skew. Nothing censored, so the local-gfortran
  drift ``CLAUDE.md`` describes does not reach this comparison.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.fixtures.paths import SKIP_REASON, TESTDATA_AVAILABLE
from tests.fixtures.paths import testdata_path as _testdata

pytest.importorskip(
    "flowfreq.peakfqr",
    reason="Fortran extension not built; run python build_fortran/build.py "
    "(see docs/FORTRAN_UPLOAD.md)",
    # See tests/fortran_parity/test_fortran_oracles.py's exc_type comment.
    exc_type=ImportError,
)

pytestmark = [
    pytest.mark.requires_fortran,
    pytest.mark.requires_peakfqr_testdata,
    pytest.mark.skipif(not TESTDATA_AVAILABLE, reason=SKIP_REASON),
]


@pytest.fixture(scope="module")
def stations():
    from flowfreq.psf_convert import convert_psf

    return convert_psf(_testdata("wymt_ffa_2022A.psf"))


def test_fortran_engine_via_bulletin17c_is_bit_identical(stations):
    checked = []
    for sid, s in stations.items():
        if s.unsupported_reasons("fortran"):
            continue
        ref = s.fortran_reference()
        results = s.run("fortran")
        assert results.mean_log == ref.parameters["mean_log"], sid
        assert results.skew_used == ref.parameters["skew_weighted"], sid
        q = results.quantiles
        for aep, flow in ref.quantiles.items():
            got = float(q.loc[np.isclose(q["aep"], aep), "flow_cfs"].iloc[0])
            assert got == flow, (sid, aep)
        checked.append(sid)
    assert len(checked) == len(stations) - 1


def test_every_station_runs_directly(stations):
    for sid, s in stations.items():
        ref = s.fortran_reference()
        assert np.isfinite(ref.parameters["mean_log"]), sid


def test_native_powder_river_matches_live_fortran(stations):
    s = stations["06326500.00"]
    ref = s.fortran_reference()
    results = s.run("native")
    p = ref.parameters
    assert abs(results.mean_log - p["mean_log"]) < 1e-6
    assert abs(results.std_log - p["std_log"]) < 1e-6
    assert abs(results.skew_weighted - p["skew_weighted"]) < 1e-6
    assert abs(results.skew_station - p["skew_at_site"]) < 1e-6


def _zero_flow_stations(stations):
    return {
        sid: s for sid, s in stations.items() if any(r.is_exact and r.ql <= 1e-20 for r in s.rows)
    }


def test_every_station_but_06328100_now_runs_natively(stations):
    """Zero-flow rows closed the last structural gap; 06328100's 407 cfs upper
    perception threshold is the one thing Bulletin17C still cannot express."""
    refused = sorted(sid for sid, s in stations.items() if s.unsupported_reasons("native"))
    assert refused == ["06328100.00"]


def test_native_zero_flow_stations_match_live_fortran(stations):
    """12 WY/MT stations with zero flows, refused natively until zero years got rows.

    Measured worst: weighted skew 8.2e-5 and quantiles 0.012% (06329350);
    every other station under 1e-5 and 0.004%. Took siteQT's Qmin row for each
    zero year, gbtest's gbtmin (1e-6) lower bound on censored rows (06329570
    was 4.6% off without it), and an exact MGBT cutoff (06328900 and 06326960
    were 75% and 54% off when 10**log10 rounding censored a peak on the cutoff).
    """
    zero = {sid: s for sid, s in _zero_flow_stations(stations).items() if sid != "06328100.00"}
    assert len(zero) == 12
    for sid, s in zero.items():
        ref = s.fortran_reference()
        results = s.run("native")
        p = ref.parameters
        assert abs(results.mean_log - p["mean_log"]) < 1e-5, sid
        assert abs(results.skew_used - p["skew_weighted"]) < 2e-4, sid
        q = results.quantiles
        for aep, flow in ref.quantiles.items():
            got = float(q.loc[np.isclose(q["aep"], aep), "flow_cfs"].iloc[0])
            assert abs(got / flow - 1) < 5e-4, (sid, aep, got, flow)
