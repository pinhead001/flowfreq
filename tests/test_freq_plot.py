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
    from tests.peakfqsa.fixtures.big_sandy import (
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
