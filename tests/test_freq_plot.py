"""Tests for flowfreq.freq_plot."""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from flowfreq.bulletin17c import Bulletin17C
from flowfreq.freq_plot import _lp3_quantiles, plot_frequency_curve, plot_peak_flows_with_thresholds

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
    fig = plot_frequency_curve(big_sandy_b17c)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_axes_labels(big_sandy_b17c):
    """Verify y-axis label contains discharge info."""
    fig = plot_frequency_curve(big_sandy_b17c)
    ax = fig.axes[0]
    ylabel = ax.get_ylabel()
    assert "cfs" in ylabel.lower() or "discharge" in ylabel.lower()
    plt.close(fig)


def test_with_big_sandy_data(big_sandy_b17c):
    """Run with Big Sandy fixture; verify no exception and plot elements exist."""
    fig = plot_frequency_curve(
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
        base = plot_frequency_curve(big_sandy_b17c)
        n_lines, n_collections = len(base.axes[0].lines), len(base.axes[0].collections)
        plt.close(base)

        fig = plot_frequency_curve(big_sandy_b17c, extra_curves=self.OVERLAY)
        ax = fig.axes[0]
        # the curve itself plus two dotted bounds, and one filled band
        assert len(ax.lines) == n_lines + 3
        assert len(ax.collections) == n_collections + 1
        plt.close(fig)

    def test_overlay_is_labelled_for_the_legend(self, big_sandy_b17c):
        fig = plot_frequency_curve(big_sandy_b17c, extra_curves=self.OVERLAY)
        labels = [ln.get_label() for ln in fig.axes[0].lines]
        assert "Threshold-aware EMA" in labels
        plt.close(fig)

    def test_overlay_uses_its_own_moments(self, big_sandy_b17c):
        """A different mean must move the curve; otherwise the parameter is decorative."""
        low = plot_frequency_curve(big_sandy_b17c, extra_curves={"low": (3.0, 0.29, -0.16, 44)})
        high = plot_frequency_curve(big_sandy_b17c, extra_curves={"high": (4.0, 0.29, -0.16, 44)})
        y_low = [ln for ln in low.axes[0].lines if ln.get_label() == "low"][0].get_ydata()
        y_high = [ln for ln in high.axes[0].lines if ln.get_label() == "high"][0].get_ydata()
        assert np.all(y_high > y_low)
        plt.close(low)
        plt.close(high)

    def test_none_and_empty_are_no_ops(self, big_sandy_b17c):
        base = plot_frequency_curve(big_sandy_b17c)
        n_lines = len(base.axes[0].lines)
        plt.close(base)
        for value in (None, {}):
            fig = plot_frequency_curve(big_sandy_b17c, extra_curves=value)
            assert len(fig.axes[0].lines) == n_lines
            plt.close(fig)


class TestBackwardCompatibleAlias:
    """The app pins an older tag and imports the pre-rename name."""

    def test_alias_is_the_same_function(self):
        from flowfreq import freq_plot

        assert freq_plot.plot_frequency_curve_streamlit is freq_plot.plot_frequency_curve

    def test_alias_still_importable_by_name(self):
        """How flowfreq-app imports it today; breaking this breaks the app."""
        from flowfreq.freq_plot import plot_frequency_curve_streamlit

        assert callable(plot_frequency_curve_streamlit)


@pytest.fixture
def peak_df():
    """Simple synthetic annual-peak record, unrelated to any LP3 fit."""
    years = np.arange(1980, 2020)
    rng = np.random.default_rng(0)
    flows = 1000.0 * 10.0 ** rng.normal(0, 0.25, size=len(years))
    return pd.DataFrame({"water_year": years, "peak_flow_cfs": flows})


class TestPlotPeakFlowsWithThresholds:
    """Basic behavior of the bar-chart plotter, previously untested."""

    def test_returns_figure(self, peak_df):
        fig = plot_peak_flows_with_thresholds(peak_df)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_one_bar_per_year(self, peak_df):
        fig = plot_peak_flows_with_thresholds(peak_df)
        ax = fig.axes[0]
        assert len(ax.patches) == len(peak_df)
        plt.close(fig)

    def test_mgbt_threshold_drawn(self, peak_df):
        fig = plot_peak_flows_with_thresholds(peak_df, mgbt_threshold=500.0)
        ax = fig.axes[0]
        assert any(abs(line.get_ydata()[0] - 500.0) < 1e-6 for line in ax.lines)
        plt.close(fig)

    def test_peaks_below_the_threshold_are_censored_hollow(self, peak_df):
        """A PILF/MGBT cut hollows out the bars below it -- the third feature
        moved in from the app's plot_peak_timeseries (TODO.md), so an applied
        override is visible on the record itself."""
        threshold = float(peak_df["peak_flow_cfs"].median())
        n_below = int((peak_df["peak_flow_cfs"] < threshold).sum())
        assert 0 < n_below < len(peak_df)  # fixture actually exercises both sides

        fig = plot_peak_flows_with_thresholds(peak_df, mgbt_threshold=threshold)
        ax = fig.axes[0]
        # Total bars drawn still equals one per year, split solid/hollow.
        assert len(ax.patches) == len(peak_df)
        hollow_patches = [p for p in ax.patches if p.get_facecolor()[3] == 0]
        assert len(hollow_patches) == n_below
        plt.close(fig)

    def test_no_threshold_leaves_every_bar_solid(self, peak_df):
        fig = plot_peak_flows_with_thresholds(peak_df)
        ax = fig.axes[0]
        assert all(p.get_facecolor()[3] != 0 for p in ax.patches)
        plt.close(fig)

    def test_threshold_source_appears_in_the_legend_label(self, peak_df):
        fig = plot_peak_flows_with_thresholds(
            peak_df, mgbt_threshold=500.0, mgbt_threshold_source="override"
        )
        ax = fig.axes[0]
        labels = [line.get_label() for line in ax.lines]
        assert any("override" in lbl for lbl in labels)
        plt.close(fig)


class TestReturnPeriodLinesAndMaxPeakAnnotation:
    """The two features ported in from the app's ``plot_peak_timeseries`` (TODO.md).

    ``plot_peak_flows_with_thresholds`` previously had neither return-period
    reference lines nor a max-peak recurrence annotation, so switching the app
    to the library function would have lost both. ``lp3_params`` is the new
    opt-in that supplies them.
    """

    # Roughly Big Sandy's weighted-skew fit; exact values don't matter here,
    # only that the same triple is used to compute expectations and to plot.
    LP3 = (3.7175, 0.2910, -0.1563)

    def test_no_lp3_params_is_a_no_op(self, peak_df):
        base = plot_peak_flows_with_thresholds(peak_df)
        n_lines, n_texts = len(base.axes[0].lines), len(base.axes[0].texts)
        plt.close(base)

        fig = plot_peak_flows_with_thresholds(peak_df, lp3_params=None)
        assert len(fig.axes[0].lines) == n_lines
        assert len(fig.axes[0].texts) == n_texts
        plt.close(fig)

    def test_one_reference_line_per_return_period(self, peak_df):
        fig = plot_peak_flows_with_thresholds(
            peak_df, lp3_params=self.LP3, return_periods=(2, 10, 100)
        )
        ax = fig.axes[0]
        assert len(ax.lines) == 3
        plt.close(fig)

    def test_reference_lines_match_the_lp3_quantiles(self, peak_df):
        return_periods = (2, 10, 100)
        fig = plot_peak_flows_with_thresholds(
            peak_df, lp3_params=self.LP3, return_periods=return_periods
        )
        ax = fig.axes[0]
        drawn = sorted(line.get_ydata()[0] for line in ax.lines)
        mean_log, std_log, skew = self.LP3
        expected = sorted(_lp3_quantiles(mean_log, std_log, skew, 1.0 / np.array(return_periods)))
        assert drawn == pytest.approx(expected, rel=1e-9)
        plt.close(fig)

    def test_reference_lines_get_a_single_legend_entry(self, peak_df):
        fig = plot_peak_flows_with_thresholds(
            peak_df, lp3_params=self.LP3, return_periods=(2, 10, 100)
        )
        ax = fig.axes[0]
        labels = [ln.get_label() for ln in ax.lines]
        assert labels.count("Return Period (LP3 fit)") == 1
        plt.close(fig)

    def test_return_period_lines_are_labelled_on_the_plot(self, peak_df):
        fig = plot_peak_flows_with_thresholds(
            peak_df, lp3_params=self.LP3, return_periods=(2, 10, 100)
        )
        texts = {t.get_text() for t in fig.axes[0].texts}
        assert {"2-yr", "10-yr", "100-yr"} <= texts
        plt.close(fig)

    def test_max_peak_annotation_reports_correct_recurrence_interval(self):
        """A peak sized to exactly the 100-yr LP3 quantile should read ~100-yr."""
        mean_log, std_log, skew = self.LP3
        years = np.arange(1980, 2000)
        flows = np.full(len(years), 500.0)
        hundred_yr_flow = float(_lp3_quantiles(mean_log, std_log, skew, np.array([0.01]))[0])
        flows[-1] = hundred_yr_flow
        df = pd.DataFrame({"water_year": years, "peak_flow_cfs": flows})

        fig = plot_peak_flows_with_thresholds(df, lp3_params=self.LP3, return_periods=(2, 10, 100))
        texts = [t.get_text() for t in fig.axes[0].texts]
        assert any("cfs" in t and "100-yr" in t for t in texts)
        plt.close(fig)

    def test_annotate_max_peak_false_suppresses_the_annotation(self, peak_df):
        fig = plot_peak_flows_with_thresholds(peak_df, lp3_params=self.LP3, annotate_max_peak=False)
        texts = [t.get_text() for t in fig.axes[0].texts]
        assert not any("cfs" in t for t in texts)
        plt.close(fig)

    def test_a_peak_far_beyond_the_largest_return_period_reports_a_bound(self):
        """The LP3 CDF saturates; the annotation must not divide by zero."""
        mean_log, std_log, skew = self.LP3
        years = np.arange(1980, 2000)
        flows = np.full(len(years), 500.0)
        flows[-1] = 10.0 ** (mean_log + 20 * std_log)  # absurdly far out on the tail
        df = pd.DataFrame({"water_year": years, "peak_flow_cfs": flows})

        fig = plot_peak_flows_with_thresholds(df, lp3_params=self.LP3, return_periods=(2, 10, 100))
        texts = [t.get_text() for t in fig.axes[0].texts]
        assert any("> 100-yr" in t for t in texts)
        plt.close(fig)
