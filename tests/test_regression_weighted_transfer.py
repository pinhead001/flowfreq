"""Tests for the SIR 2016-5118 (Ries, 2007) regression-weighted area-ratio transfer."""

from __future__ import annotations

import numpy as np
import pytest

from flowfreq.transpose import RIES_2007_RANGE, regression_weighted_transfer

# SIR 2016-5118 p. 45: Mineral Creek near Mineral, WA (12083000), 2 % AEP.
A_G, A_U, B = 76.75, 56.96, 0.897
Q_G_WTD, Q_U_REG = 10_600.0, 6_990.0


class TestWorkedExample:
    def test_equation_11_to_the_printed_precision(self):
        out = regression_weighted_transfer(Q_G_WTD, Q_U_REG, A_G, A_U, B)
        assert round(float(out.transferred_cfs), -1) == 8_110.0

    def test_equation_12_from_the_printed_intermediate(self):
        """The report carries 8,110 into eq. 12 and prints 7,532.4."""
        w = 2 * abs(A_G - A_U) / A_G
        assert w * Q_U_REG + (1 - w) * 8_110.0 == pytest.approx(7_532.4, abs=0.05)

    def test_equation_12_to_the_printed_answer(self):
        out = regression_weighted_transfer(Q_G_WTD, Q_U_REG, A_G, A_U, B)
        assert float(out.regression_weight) == pytest.approx(39.58 / 76.75, abs=1e-4)
        assert round(float(out.flow_cfs), -1) == 7_530.0


class TestWeights:
    def test_ratio_one_is_the_transfer(self):
        out = regression_weighted_transfer(1000.0, 400.0, 50.0, 50.0, 0.8)
        assert float(out.regression_weight) == 0.0
        assert float(out.flow_cfs) == pytest.approx(1000.0)

    @pytest.mark.parametrize("ratio", RIES_2007_RANGE)
    def test_band_edges_are_the_regression(self, ratio):
        out = regression_weighted_transfer(1000.0, 400.0, 50.0, 50.0 * ratio, 0.8)
        assert float(out.regression_weight) == pytest.approx(1.0)
        assert float(out.flow_cfs) == pytest.approx(400.0)

    def test_vectorized_over_aeps(self):
        out = regression_weighted_transfer(
            np.array([1000.0, 2000.0]), np.array([800.0, 1500.0]), 50.0, 60.0, np.array([0.8, 0.7])
        )
        assert out.flow_cfs.shape == (2,)
        expected_t = np.array([1000.0, 2000.0]) * 1.2 ** np.array([0.8, 0.7])
        np.testing.assert_allclose(out.transferred_cfs, expected_t)
        w = 2 * 10.0 / 50.0
        np.testing.assert_allclose(
            out.flow_cfs, w * np.array([800.0, 1500.0]) + (1 - w) * expected_t
        )


class TestGuardrails:
    def test_outside_the_band_raises(self):
        with pytest.raises(ValueError, match="0.5"):
            regression_weighted_transfer(1000.0, 400.0, 50.0, 80.0, 0.8)

    def test_non_positive_area_raises(self):
        with pytest.raises(ValueError, match="positive"):
            regression_weighted_transfer(1000.0, 400.0, 0.0, 10.0, 0.8)

    def test_negative_flow_raises(self):
        with pytest.raises(ValueError, match="flow"):
            regression_weighted_transfer(-1.0, 400.0, 50.0, 50.0, 0.8)

    def test_shapes_must_broadcast(self):
        with pytest.raises(ValueError):
            regression_weighted_transfer(np.ones(2), np.ones(3), 50.0, 50.0, 0.8)
