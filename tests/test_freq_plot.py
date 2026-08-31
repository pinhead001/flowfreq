"""Tests for hydrolib.freq_plot."""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest

from hydrolib.bulletin17c import Bulletin17C
from hydrolib.freq_plot import plot_frequency_curve_streamlit

matplotlib.use("Agg")


@pytest.fixture
def big_sandy_b17c():
    """Create a Bulletin17C instance from Big Sandy fixture data."""
    from tests.fixtures.big_sandy import (
        HISTORICAL_PEAKS,
        REGIONAL_SKEW,
        REGIONAL_SKEW_SD,
        SYSTEMATIC_PEAKS,
    )

    years = np.array(sorted(SYSTEMATIC_PEAKS.keys()))
    flows = np.array([SYSTEMATIC_PEAKS[y] for y in years])
    historical = [(y, q) for y, q in HISTORICAL_PEAKS.items()]

    b17c = Bulletin17C(
        peak_flows=flows,
        water_years=years,
        regional_skew=REGIONAL_SKEW,
        regional_skew_mse=REGIONAL_SKEW_SD**2,
        historical_peaks=historical,
    )
    b17c.run_analysis(method="ema")
    return b17c


def test_returns_figure(big_sandy_b17c):
    """Verify function returns a matplotlib Figure."""
    fig = plot_frequency_curve_streamlit(big_sandy_b17c)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_axes_labels(big_sandy_b17c):
    """Verify y-axis label contains discharge info."""
    fig = plot_frequency_curve_streamlit(big_sandy_b17c)
    ax = fig.axes[0]
    ylabel = ax.get_ylabel()
    assert "cfs" in ylabel.lower() or "discharge" in ylabel.lower()
    plt.close(fig)


def test_with_big_sandy_data(big_sandy_b17c):
    """Run with Big Sandy fixture; verify no exception and plot elements exist."""
    fig = plot_frequency_curve_streamlit(
        big_sandy_b17c,
        site_name="Big Sandy River at Bruceton, TN",
        site_no="03606500",
    )
    ax = fig.axes[0]
    # Should have scatter + line + fill + CI dashes at minimum
    assert len(ax.collections) >= 1  # scatter + fill
    assert len(ax.lines) >= 1  # fitted curve
    assert ax.get_yscale() == "log"
    plt.close(fig)


class TestExtraCurves:
    """Overlay curves carrying their own moments, not just their own skew.

    Ported from the `dev` branch, the one library delta on it that main had not
    superseded. `skew_curves` varies only the skew of the fit already on the
    axes; showing the effect of a *refit* -- a PILF-threshold override, a
    perception-threshold-aware EMA -- needs its own mean and standard
    deviation too.
    """

    OVERLAY = {"Threshold-aware EMA": (3.72, 0.29, -0.16, 44)}

    def test_overlay_adds_a_curve_and_a_band(self, big_sandy_b17c):
        base = plot_frequency_curve_streamlit(big_sandy_b17c)
        n_lines, n_collections = len(base.axes[0].lines), len(base.axes[0].collections)
        plt.close(base)

        fig = plot_frequency_curve_streamlit(big_sandy_b17c, extra_curves=self.OVERLAY)
        ax = fig.axes[0]
        # the curve itself plus two dotted bounds, and one filled band
        assert len(ax.lines) == n_lines + 3
        assert len(ax.collections) == n_collections + 1
        plt.close(fig)

    def test_overlay_is_labelled_for_the_legend(self, big_sandy_b17c):
        fig = plot_frequency_curve_streamlit(big_sandy_b17c, extra_curves=self.OVERLAY)
        labels = [ln.get_label() for ln in fig.axes[0].lines]
        assert "Threshold-aware EMA" in labels
        plt.close(fig)

    def test_overlay_uses_its_own_moments(self, big_sandy_b17c):
        """A different mean must move the curve; otherwise the parameter is decorative."""
        low = plot_frequency_curve_streamlit(
            big_sandy_b17c, extra_curves={"low": (3.0, 0.29, -0.16, 44)}
        )
        high = plot_frequency_curve_streamlit(
            big_sandy_b17c, extra_curves={"high": (4.0, 0.29, -0.16, 44)}
        )
        y_low = [ln for ln in low.axes[0].lines if ln.get_label() == "low"][0].get_ydata()
        y_high = [ln for ln in high.axes[0].lines if ln.get_label() == "high"][0].get_ydata()
        assert np.all(y_high > y_low)
        plt.close(low)
        plt.close(high)

    def test_none_and_empty_are_no_ops(self, big_sandy_b17c):
        base = plot_frequency_curve_streamlit(big_sandy_b17c)
        n_lines = len(base.axes[0].lines)
        plt.close(base)
        for value in (None, {}):
            fig = plot_frequency_curve_streamlit(big_sandy_b17c, extra_curves=value)
            assert len(fig.axes[0].lines) == n_lines
            plt.close(fig)
