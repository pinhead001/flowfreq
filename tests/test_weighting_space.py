"""A state's `weighting_space`: how its report combines the regions of a basin that spans two.

Fifteen reports average the regions' flows by drainage area, so their files set "linear".
The rest keep the library's default, an average in log: their reports say nothing on such a
basin, or were not read (TODO.md).
"""

from __future__ import annotations

import pytest

from flowfreq import estimate_from_characteristics
from flowfreq.regression.library import available_states, load_state

#: States whose report says to weight the flows. FL and NJ say so in words only.
LINEAR = ["AZ", "CA", "CO", "FL", "GA", "IA", "KY", "MT", "NC", "NJ", "NV", "NY", "SC", "SD", "WY"]


def _flow(state, characteristics, region, aep, **kwargs):
    estimate = estimate_from_characteristics(
        state, characteristics, region=region, lat=40.0, lon=-100.0, aeps=[aep], **kwargs
    )
    assert len(estimate.quantiles) == 1
    return estimate, float(estimate.quantiles.loc[0, "flow_cfs"])


@pytest.mark.parametrize("state", LINEAR)
def test_these_states_weight_in_flow(state):
    lib = load_state(state)
    assert lib.weighting_space == "linear"
    assert "weighting_space" in lib.notes  # the file says which page of the report


def test_no_other_state_sets_a_weighting_space():
    others = [s for s in available_states() if s not in LINEAR]
    assert others
    assert [s for s in others if load_state(s).weighting_space is not None] == []


def test_new_york_example_2_is_the_default():
    """SIR 2006-5112, Genesee River: 53.5 percent region 5, 46.5 percent region 6, 79,900."""
    characteristics = {
        "DRNAREA": 2467,
        "CSL10_85": 8.05,
        "PRECIP": 33.92,
        "STORAGE": 1.08,
        "MAR": 14.64,
        "EL1200": 58.8,
        "SLOPERATIO": 0.019,
    }
    regions = {"5": 53.5, "6": 46.5}
    estimate, flow = _flow("NY", characteristics, regions, 0.02)
    assert flow == pytest.approx(79900, rel=5e-4)
    assert estimate.provenance["options"]["weighting_space"] == "linear"
    logged, in_log = _flow("NY", characteristics, regions, 0.02, weighting_space="log")
    assert in_log == pytest.approx(74800, rel=1e-3)
    assert logged.provenance["options"]["weighting_space"] == "log"


@pytest.mark.parametrize(("aep", "printed", "in_log"), [(0.1, 1680, 1210), (0.01, 9500, 6930)])
def test_nevada_equation_6_is_the_default(aep, printed, in_log):
    """WSP 2433, p. 41: 57 mi2 at 6,500 ft, 21 mi2 in region 6 and 36 in region 10."""
    characteristics = {"DRNAREA": 57, "ELEV": 6500}
    regions = {"6": 21 / 57 * 100, "10": 36 / 57 * 100}
    _, flow = _flow("NV", characteristics, regions, aep)
    assert flow == pytest.approx(printed, rel=3e-3)
    _, logged = _flow("NV", characteristics, regions, aep, weighting_space="log")
    assert logged == pytest.approx(in_log, rel=3e-3)


def test_nevada_region_6_zero_flow_weights_in_flow_only():
    """Region 6 prints Q = 0 at the 2-year flood, which has no logarithm to average."""
    characteristics = {"DRNAREA": 57, "ELEV": 6500}
    regions = {"6": 21 / 57 * 100, "10": 36 / 57 * 100}
    _, flow = _flow("NV", characteristics, regions, 0.5)
    assert flow > 0
    with pytest.raises(ValueError, match="cannot weight in log space"):
        _flow("NV", characteristics, regions, 0.5, weighting_space="log")


def test_montana_two_region_example_is_the_default():
    """SIR 2015-5019-F, p. 23: Q10 = 0.36 (339) + 0.64 (420) = 391 ft3/s."""
    characteristics = {"CONTDA": 17.27, "EL6000": 0.0, "SLOP30_30M": 0.05, "ET0306MOD": 1.15}
    upstream = 6.21 / 17.27 * 100
    regions = {"UYCM": upstream, "ECP": 100 - upstream}
    _, flow = _flow("MT", characteristics, regions, 0.1)
    assert round(flow) == 391
    _, logged = _flow("MT", characteristics, regions, 0.1, weighting_space="log")
    assert round(logged) == 389


def test_an_unknown_weighting_space_is_an_error():
    with pytest.raises(ValueError, match="weighting_space must be 'log' or 'linear'"):
        _flow("NY", {"DRNAREA": 2467}, {"5": 53.5, "6": 46.5}, 0.02, weighting_space="area")


def test_colorado_set_awaiting_nss_weights_in_flow_too():
    """SIR 2025-5047, p. 32, repeats SIR 2016-5099's rule."""
    lib = load_state("CO", allow_awaiting_nss=True, warn=False)
    assert lib.status == "awaiting_nss"
    assert lib.weighting_space == "linear"
