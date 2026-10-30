"""Interpolation between two gages on one stream (methow-flow priority item 4).

Upstreamed from methow-flow's ``reach_flows`` bracketing: a site between an upstream and
a downstream gage takes both gages' floods, each transposed to the site by the
drainage-area ratio, blended in log space with a weight linear in log drainage area. At
either gage's area the result is that gage's own transposed value, so a reach estimate
has no step where a single nearest donor would switch from one gage to the other.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from flowfreq.bulletin17c import Bulletin17C
from flowfreq.transpose import (
    RegressionExponents,
    bracket_weight,
    interpolate_between_gages,
    transpose_frequency,
)

EXPONENTS = RegressionExponents(
    aeps=[0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002],
    exponents=[0.83, 0.81, 0.80, 0.78, 0.77, 0.76, 0.75, 0.74],
    citation="Example SIR 0000-0000 table 6, hydrologic region 3",
)
WIDE = dict(area_ratio_range=(0.2, 5.0), max_aep=0.5)
UP_AREA, DOWN_AREA = 100.0, 300.0


@pytest.fixture(scope="module")
def fits():
    from tests.fixtures.big_sandy import REGIONAL_SKEW, REGIONAL_SKEW_SD, SYSTEMATIC_PEAKS

    years = np.array(sorted(SYSTEMATIC_PEAKS))
    flows = np.array([SYSTEMATIC_PEAKS[y] for y in years], dtype=float)

    def fit(scale: float):
        return Bulletin17C(
            peak_flows=flows * scale,
            water_years=years,
            regional_skew=REGIONAL_SKEW,
            regional_skew_mse=REGIONAL_SKEW_SD**2,
        ).run_analysis(method="ema")

    return fit(0.4), fit(1.0)


class TestBracketWeight:
    def test_ends_and_middle(self) -> None:
        assert bracket_weight(UP_AREA, DOWN_AREA, UP_AREA) == 1.0
        assert bracket_weight(UP_AREA, DOWN_AREA, DOWN_AREA) == 0.0
        mid = math.sqrt(UP_AREA * DOWN_AREA)
        assert bracket_weight(UP_AREA, DOWN_AREA, mid) == pytest.approx(0.5)

    @pytest.mark.parametrize(
        "up,down,target",
        [(100.0, 300.0, 50.0), (100.0, 300.0, 301.0), (100.0, 100.0, 100.0), (0.0, 300.0, 10.0)],
    )
    def test_refuses_a_target_outside_the_bracket(self, up, down, target) -> None:
        with pytest.raises(ValueError, match="between"):
            bracket_weight(up, down, target)


class TestInterpolateBetweenGages:
    def test_columns_and_aeps(self, fits) -> None:
        up, down = fits
        out = interpolate_between_gages(up, UP_AREA, down, DOWN_AREA, 150.0, EXPONENTS, **WIDE)
        assert list(out.columns) == [
            "aep",
            "flow_cfs",
            "upstream_flow_cfs",
            "downstream_flow_cfs",
            "weight_upstream",
        ]
        t = transpose_frequency(up, UP_AREA, 150.0, EXPONENTS, **WIDE).quantiles
        assert list(out["aep"]) == pytest.approx(list(t["aep"]))
        assert (out["aep"] <= 0.5).all()

    def test_log_blend_of_the_two_transfers(self, fits) -> None:
        up, down = fits
        target = 150.0
        out = interpolate_between_gages(up, UP_AREA, down, DOWN_AREA, target, EXPONENTS, **WIDE)
        w = bracket_weight(UP_AREA, DOWN_AREA, target)
        tu = transpose_frequency(up, UP_AREA, target, EXPONENTS, **WIDE).quantiles
        td = transpose_frequency(down, DOWN_AREA, target, EXPONENTS, **WIDE).quantiles
        np.testing.assert_allclose(out["upstream_flow_cfs"], tu["flow_cfs"])
        np.testing.assert_allclose(out["downstream_flow_cfs"], td["flow_cfs"])
        expected = np.exp(w * np.log(tu["flow_cfs"]) + (1 - w) * np.log(td["flow_cfs"]))
        np.testing.assert_allclose(out["flow_cfs"], expected)
        assert (out["weight_upstream"] == w).all()

    def test_passes_through_each_gage_at_its_own_area(self, fits) -> None:
        up, down = fits
        at_up = interpolate_between_gages(up, UP_AREA, down, DOWN_AREA, UP_AREA, EXPONENTS, **WIDE)
        own = transpose_frequency(up, UP_AREA, UP_AREA, EXPONENTS, **WIDE).quantiles
        np.testing.assert_allclose(at_up["flow_cfs"], own["flow_cfs"])
        at_down = interpolate_between_gages(
            up, UP_AREA, down, DOWN_AREA, DOWN_AREA, EXPONENTS, **WIDE
        )
        own = transpose_frequency(down, DOWN_AREA, DOWN_AREA, EXPONENTS, **WIDE).quantiles
        np.testing.assert_allclose(at_down["flow_cfs"], own["flow_cfs"])

    def test_estimate_lies_between_the_transfers(self, fits) -> None:
        up, down = fits
        out = interpolate_between_gages(up, UP_AREA, down, DOWN_AREA, 200.0, EXPONENTS, **WIDE)
        lo = np.minimum(out["upstream_flow_cfs"], out["downstream_flow_cfs"])
        hi = np.maximum(out["upstream_flow_cfs"], out["downstream_flow_cfs"])
        assert ((out["flow_cfs"] >= lo - 1e-9) & (out["flow_cfs"] <= hi + 1e-9)).all()

    def test_refuses_a_target_outside_the_bracket(self, fits) -> None:
        up, down = fits
        with pytest.raises(ValueError, match="between"):
            interpolate_between_gages(up, UP_AREA, down, DOWN_AREA, 400.0, EXPONENTS, **WIDE)

    def test_transpose_options_reach_both_transfers(self, fits) -> None:
        up, down = fits
        # The default area-ratio band refuses the 1/3 transfer from the downstream gage.
        with pytest.raises(ValueError):
            interpolate_between_gages(up, UP_AREA, down, DOWN_AREA, UP_AREA, EXPONENTS)
