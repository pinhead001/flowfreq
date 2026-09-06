"""Tests for flowfreq.plots.

Note this is a different plotting module from `flowfreq.freq_plot` (covered by
`tests/test_freq_plot.py`): `freq_plot.plot_frequency_curve` takes a fitted
`Bulletin17C`, while `plots.plot_frequency_curve` here takes a fitted
`B17CEngine` plus the records used to fit it.
"""

from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest

from flowfreq.batch import run_multi_site
from flowfreq.core import PeakRecord
from flowfreq.engine import B17CEngine
from flowfreq.plots import (
    apply_b17c_style,
    plot_frequency_curve,
    plot_multi_site_comparison,
    plotting_positions,
)
from tests.fixtures.big_sandy import SYSTEMATIC_PEAKS

matplotlib.use("Agg")


@pytest.fixture(scope="module")
def records():
    years = sorted(SYSTEMATIC_PEAKS)
    return [PeakRecord(year=y, flow=SYSTEMATIC_PEAKS[y]) for y in years]


@pytest.fixture
def fitted_engine(records):
    engine = B17CEngine()
    engine.fit(records)
    return engine


class TestApplyB17cStyle:
    def test_sets_expected_rcparams(self):
        apply_b17c_style()
        assert plt.rcParams["figure.dpi"] == 140
        assert plt.rcParams["axes.grid"] is True
        assert plt.rcParams["font.size"] == 10

    def test_idempotent(self):
        """Calling it twice must not raise or change the outcome."""
        apply_b17c_style()
        apply_b17c_style()
        assert plt.rcParams["figure.facecolor"] == "white"
        assert plt.rcParams["axes.facecolor"] == "white"


class TestPlottingPositions:
    def test_weibull_formula(self):
        result = plotting_positions([40.0, 30.0, 20.0, 10.0])
        assert result == pytest.approx([5.0, 2.5, 5 / 3, 1.25])

    def test_length_only_ignores_values(self):
        """Documented as sorting the input, but it only ever consults len():
        the values themselves never enter the computation."""
        a = plotting_positions([1.0, 2.0, 3.0])
        b = plotting_positions([300.0, 1.0, 42.0])
        assert list(a) == list(b)

    def test_single_value(self):
        assert plotting_positions([5.0]) == pytest.approx([2.0])

    def test_object_without_len_raises(self):
        with pytest.raises(TypeError):
            plotting_positions(5)


class TestPlotFrequencyCurve:
    def test_returns_figure(self, fitted_engine, records):
        fig = plot_frequency_curve(fitted_engine, records)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_default_title_cites_bulletin_17c(self, fitted_engine, records):
        fig = plot_frequency_curve(fitted_engine, records)
        assert "Bulletin 17C" in fig.axes[0].get_title()
        plt.close(fig)

    def test_custom_title_overrides_default(self, fitted_engine, records):
        fig = plot_frequency_curve(fitted_engine, records, title="Big Sandy")
        assert fig.axes[0].get_title() == "Big Sandy"
        plt.close(fig)

    def test_stats_annotation_reports_fitted_params(self, fitted_engine, records):
        fig = plot_frequency_curve(fitted_engine, records)
        texts = [t.get_text() for t in fig.axes[0].texts]
        joined = "\n".join(texts)
        assert f"n = {fitted_engine.n}" in joined
        assert f"{fitted_engine.skew:.3f}" in joined
        plt.close(fig)

    def test_bankfull_annotation_can_be_disabled(self, fitted_engine, records):
        with_band = plot_frequency_curve(fitted_engine, records, show_bankfull=True)
        without_band = plot_frequency_curve(fitted_engine, records, show_bankfull=False)
        assert len(with_band.axes[0].patches) > len(without_band.axes[0].patches)
        plt.close(with_band)
        plt.close(without_band)

    def test_yscale_is_honoured(self, fitted_engine, records):
        fig = plot_frequency_curve(fitted_engine, records, yscale="linear")
        assert fig.axes[0].get_yscale() == "linear"
        plt.close(fig)

    def test_save_path_writes_a_file(self, fitted_engine, records, tmp_path):
        out = tmp_path / "freq.png"
        fig = plot_frequency_curve(fitted_engine, records, save_path=str(out))
        assert out.is_file()
        plt.close(fig)

    def test_unfitted_engine_raises(self, records):
        """engine.params is None until fit() runs; quantiles_with_ci refuses."""
        with pytest.raises(RuntimeError, match="fit"):
            plot_frequency_curve(B17CEngine(), records)


class TestPlotMultiSiteComparison:
    @pytest.fixture
    def two_site_results(self, records):
        return run_multi_site({"siteA": records, "siteB": records[:20]})

    def test_returns_figure(self, two_site_results):
        fig = plot_multi_site_comparison(two_site_results)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_one_bar_per_site(self, two_site_results):
        fig = plot_multi_site_comparison(two_site_results)
        assert len(fig.axes[0].get_yticklabels()) == len(two_site_results)
        plt.close(fig)

    def test_sites_sorted_by_estimate(self, two_site_results):
        fig = plot_multi_site_comparison(two_site_results, return_period=100)
        labels = [t.get_text() for t in fig.axes[0].get_yticklabels()]
        bars = sorted(fig.axes[0].patches, key=lambda p: p.get_y())
        widths = [b.get_width() for b in bars]
        assert widths == sorted(widths)
        assert set(labels) == set(two_site_results)
        plt.close(fig)

    def test_default_title_names_the_return_period(self, two_site_results):
        fig = plot_multi_site_comparison(two_site_results, return_period=50)
        assert "50-Year" in fig.axes[0].get_title()
        plt.close(fig)

    def test_error_entries_are_excluded(self, records):
        """A site whose fit failed must not appear as a bar."""
        results = run_multi_site({"good": records, "bad": records[:1]})
        assert "error" in results["bad"]
        fig = plot_multi_site_comparison(results)
        labels = [t.get_text() for t in fig.axes[0].get_yticklabels()]
        assert labels == ["good"]
        plt.close(fig)

    def test_no_valid_sites_raises(self):
        with pytest.raises(ValueError, match="No valid results"):
            plot_multi_site_comparison({})

    def test_missing_return_period_raises(self, two_site_results):
        """Every site's ci dict is keyed by the requested return periods only."""
        with pytest.raises(ValueError, match="No valid results"):
            plot_multi_site_comparison(two_site_results, return_period=99999)
