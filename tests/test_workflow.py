"""Tests for flowfreq.workflow -- the one-call analysis entry points."""

import numpy as np
import pytest

from flowfreq.workflow import (
    B17C_DEFAULT_SKEW,
    B17C_DEFAULT_SKEW_SE,
    DEFAULT_RETURN_INTERVALS,
    SKEW_OPTIONS,
    _low_outlier_source,
    build_skew_curves_dict,
    compute_skew_tables,
    resolve_regional_skew,
    run_ffa,
)
from tests.fixtures.big_sandy import REGIONAL_SKEW, REGIONAL_SKEW_SD, SYSTEMATIC_PEAKS


def _big_sandy_arrays():
    """Extract peak_flows and water_years arrays from big_sandy fixture."""
    years = np.array(sorted(SYSTEMATIC_PEAKS.keys()))
    flows = np.array([SYSTEMATIC_PEAKS[y] for y in years])
    return flows, years


class TestRunFFA:
    def test_run_ffa_returns_expected_keys(self):
        flows, years = _big_sandy_arrays()
        result = run_ffa(flows, years, use_default_skew=True)
        expected_keys = {"b17c", "converged", "method", "parameters", "quantile_df", "error"}
        assert set(result.keys()) == expected_keys

    def test_run_ffa_with_big_sandy_data(self):
        flows, years = _big_sandy_arrays()
        result = run_ffa(
            flows, years, regional_skew=REGIONAL_SKEW, regional_skew_se=REGIONAL_SKEW_SD
        )
        assert result["error"] is None
        assert result["converged"] is True
        assert result["method"] in ("ema", "mom")
        assert len(result["quantile_df"]) == 9
        assert result["b17c"] is not None

    def test_run_ffa_handles_bad_data(self):
        result = run_ffa(np.array([]), np.array([]), use_default_skew=True)
        # Should return an error for empty input
        assert result["error"] is not None

    def test_default_skew_opt_in_matches_passing_the_value(self):
        """use_default_skew=True is exactly the old silent default, now labelled."""
        flows, years = _big_sandy_arrays()
        opted = run_ffa(flows, years, use_default_skew=True)
        explicit = run_ffa(
            flows,
            years,
            regional_skew=B17C_DEFAULT_SKEW,
            regional_skew_se=B17C_DEFAULT_SKEW_SE,
        )
        assert opted["parameters"]["mean_log"] == explicit["parameters"]["mean_log"]
        assert opted["parameters"]["skew_weighted"] == explicit["parameters"]["skew_weighted"]
        assert opted["parameters"]["regional_skew_source"] == "default"
        assert explicit["parameters"]["regional_skew_source"] == "user"


class TestSkewChoiceIsRequired:
    """No silent regional skew: Bulletin 17C has no national default (p. 31)."""

    def test_no_choice_raises_before_fitting(self):
        flows, years = _big_sandy_arrays()
        with pytest.raises(ValueError, match="No regional skew chosen"):
            run_ffa(flows, years)

    def test_the_error_names_every_way_out(self):
        with pytest.raises(ValueError) as exc:
            resolve_regional_skew(None, None)
        for option in ("regional_skew_se", "station_skew_only", "use_default_skew"):
            assert option in str(exc.value)

    def test_default_opt_in_warns(self, caplog):
        with caplog.at_level("WARNING", logger="flowfreq.workflow"):
            skew, mse, source = resolve_regional_skew(None, None, use_default_skew=True)
        assert (skew, source) == (B17C_DEFAULT_SKEW, "default")
        assert mse == pytest.approx(B17C_DEFAULT_SKEW_SE**2)
        assert "unsourced" in caplog.text

    def test_station_only_passes_no_regional_skew(self):
        assert resolve_regional_skew(None, None, station_skew_only=True) == (
            None,
            None,
            "station",
        )

    def test_station_only_fit_uses_the_station_skew(self):
        flows, years = _big_sandy_arrays()
        result = run_ffa(flows, years, station_skew_only=True)
        assert result["error"] is None
        p = result["parameters"]
        assert p["regional_skew"] is None
        assert p["regional_skew_source"] == "station"
        assert p["skew_used"] == pytest.approx(p["skew_station"])

    def test_supplied_skew_is_squared_into_an_mse(self):
        assert resolve_regional_skew(-0.07, 0.6) == (-0.07, pytest.approx(0.36), "user")

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"regional_skew": -0.07, "regional_skew_se": 0.6, "use_default_skew": True},
            {"regional_skew": -0.07, "regional_skew_se": 0.6, "station_skew_only": True},
            {"use_default_skew": True, "station_skew_only": True},
        ],
    )
    def test_more_than_one_choice_raises(self, kwargs):
        args = {"regional_skew": None, "regional_skew_se": None, **kwargs}
        with pytest.raises(ValueError, match="Choose one skew source"):
            resolve_regional_skew(**args)

    @pytest.mark.parametrize("skew, se", [(-0.07, None), (None, 0.6)])
    def test_skew_and_se_go_together(self, skew, se):
        with pytest.raises(ValueError, match="go together"):
            resolve_regional_skew(skew, se)

    def test_non_positive_se_raises(self):
        with pytest.raises(ValueError, match="must be positive"):
            resolve_regional_skew(-0.07, 0.0)


class TestPilfOverride:
    """The override has to change the fit, or the control is decorative."""

    OVERRIDE = 2000.0  # censors four Big Sandy peaks; EMA still converges

    def test_override_censors_more_peaks_than_mgbt(self):
        flows, years = _big_sandy_arrays()
        base = run_ffa(peak_flows=flows, water_years=years, use_default_skew=True)
        forced = run_ffa(
            peak_flows=flows,
            water_years=years,
            low_outlier_threshold_override=self.OVERRIDE,
            use_default_skew=True,
        )
        assert forced["parameters"]["n_low_outliers"] > base["parameters"]["n_low_outliers"]
        assert forced["parameters"]["low_outlier_threshold"] == pytest.approx(self.OVERRIDE)
        assert forced["parameters"]["low_outlier_source"] == "override"
        assert base["parameters"]["low_outlier_source"] == "MGBT"

    def test_override_moves_the_fitted_moments(self):
        flows, years = _big_sandy_arrays()
        base = run_ffa(peak_flows=flows, water_years=years, use_default_skew=True)
        forced = run_ffa(
            peak_flows=flows,
            water_years=years,
            low_outlier_threshold_override=self.OVERRIDE,
            use_default_skew=True,
        )
        assert base["parameters"]["mean_log"] != forced["parameters"]["mean_log"]

    def test_zero_and_none_both_mean_use_mgbt(self):
        flows, years = _big_sandy_arrays()
        zero = run_ffa(
            peak_flows=flows,
            water_years=years,
            use_default_skew=True,
            low_outlier_threshold_override=0.0,
        )
        none = run_ffa(
            peak_flows=flows,
            water_years=years,
            use_default_skew=True,
            low_outlier_threshold_override=None,
        )
        assert zero["parameters"]["low_outlier_source"] == "MGBT"
        assert zero["parameters"]["mean_log"] == none["parameters"]["mean_log"]

    def test_mom_fallback_censors_on_the_override(self):
        """A high override stops EMA converging, and MOM now censors on it too.

        MethodOfMoments applies the Bulletin 17B conditional-probability
        adjustment, so the number reported as the threshold is also the one
        that shaped the fit.
        """
        flows, years = _big_sandy_arrays()
        result = run_ffa(
            peak_flows=flows,
            water_years=years,
            use_default_skew=True,
            low_outlier_threshold_override=6000.0,
        )
        assert result["method"] == "mom"
        assert result["parameters"]["low_outlier_threshold"] == pytest.approx(6000.0)
        assert result["parameters"]["low_outlier_source"] == "override"

    def test_mom_fallback_moments_shift_with_the_override(self):
        """The claim the label makes must actually hold: MOM censors now."""
        flows, years = _big_sandy_arrays()
        forced = run_ffa(
            peak_flows=flows,
            water_years=years,
            use_default_skew=True,
            low_outlier_threshold_override=6000.0,
        )
        assert forced["method"] == "mom"
        plain = run_ffa(peak_flows=flows, water_years=years, use_default_skew=True)
        assert forced["parameters"]["mean_log"] != pytest.approx(plain["parameters"]["mean_log"])

    def test_source_helper_covers_every_combination(self):
        assert _low_outlier_source(None) == "MGBT"
        assert _low_outlier_source(1000.0) == "override"


class TestComputeSkewTables:
    """Public API since the split; it had no direct coverage while it lived in the app."""

    def test_one_table_per_selected_label(self):
        flows, years = _big_sandy_arrays()
        result = run_ffa(flows, years, use_default_skew=True)
        tables = compute_skew_tables(result, SKEW_OPTIONS)
        assert set(tables) == set(SKEW_OPTIONS)
        for df in tables.values():
            assert len(df) == len(DEFAULT_RETURN_INTERVALS)
            assert list(df.columns) == [
                "Return Interval (yr)",
                "AEP (%)",
                "Flow (cfs)",
                "Lower 90% CI",
                "Upper 90% CI",
            ]

    def test_bounds_bracket_the_estimate(self):
        flows, years = _big_sandy_arrays()
        tables = compute_skew_tables(
            run_ffa(flows, years, use_default_skew=True), ["Weighted Skew"]
        )
        df = tables["Weighted Skew"]
        assert (df["Lower 90% CI"] < df["Flow (cfs)"]).all()
        assert (df["Flow (cfs)"] < df["Upper 90% CI"]).all()

    def test_differing_skews_give_differing_curves(self):
        flows, years = _big_sandy_arrays()
        tables = compute_skew_tables(
            run_ffa(flows, years, use_default_skew=True), ["Station Skew", "Regional Skew"]
        )
        station = tables["Station Skew"]["Flow (cfs)"].to_numpy()
        regional = tables["Regional Skew"]["Flow (cfs)"].to_numpy()
        assert not np.allclose(station, regional)

    def test_failed_analysis_yields_no_tables(self):
        assert compute_skew_tables({"error": "boom", "b17c": None}, SKEW_OPTIONS) == {}
        assert compute_skew_tables({"error": None, "b17c": None}, SKEW_OPTIONS) == {}

    def test_unknown_label_is_skipped_not_raised(self):
        flows, years = _big_sandy_arrays()
        tables = compute_skew_tables(
            run_ffa(flows, years, use_default_skew=True), ["Weighted Skew", "Nonsense Skew"]
        )
        assert set(tables) == {"Weighted Skew"}


class TestBuildSkewCurvesDict:
    def test_returns_selected_skews_only(self):
        flows, years = _big_sandy_arrays()
        result = run_ffa(flows, years, use_default_skew=True)
        curves = build_skew_curves_dict(result, ["Station Skew", "Regional Skew"])
        assert set(curves) == {"Station Skew", "Regional Skew"}
        assert curves["Regional Skew"] == pytest.approx(-0.302)

    def test_empty_when_nothing_resolves(self):
        """An empty dict is the documented signal to fall back to the default."""
        assert build_skew_curves_dict({"parameters": {}}, SKEW_OPTIONS) == {}
