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
    compare_engines,
    compute_skew_tables,
    peak_code_kwargs,
    resolve_regional_skew,
    run_ffa,
)
from tests.fixtures.big_sandy import REGIONAL_SKEW, REGIONAL_SKEW_SD, SYSTEMATIC_PEAKS
from tests.fixtures.paths import SKIP_REASON, TESTDATA_AVAILABLE
from tests.fixtures.paths import testdata_path as _testdata


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

    @pytest.fixture
    def ema_never_converges(self, monkeypatch):
        """Force the MOM fallback.

        A 6000 cfs override used to trip it on its own: the fixed point then
        stopped after 100 iterations, short of convergence. Following
        ``p3est_ema``'s own iteration limit it converges (827 iterations), as
        it does in peakfq, so the fallback is exercised directly instead.
        """
        from flowfreq.bulletin17c import ExpectedMomentsAlgorithm

        real = ExpectedMomentsAlgorithm._ema_fixed_point

        def not_converged(self, *args, **kwargs):
            mean, std, skew, _converged, iterations = real(self, *args, **kwargs)
            return mean, std, skew, False, iterations

        monkeypatch.setattr(ExpectedMomentsAlgorithm, "_ema_fixed_point", not_converged)

    def test_mom_fallback_censors_on_the_override(self, ema_never_converges):
        """When EMA does not converge, MOM censors on the override too.

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

    def test_mom_fallback_moments_shift_with_the_override(self, ema_never_converges):
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


# --------------------------------------------------------------------------- #
# Peak qualification codes (#30)
# --------------------------------------------------------------------------- #

_CODE_YEARS = np.arange(2001, 2016)
_CODE_FLOWS = np.array([float(100 + 37 * i % 211) for i in range(len(_CODE_YEARS))])


def _codes(at: dict) -> list:
    """Blank codes, with ``at[index] = code`` set."""
    out = [""] * len(_CODE_YEARS)
    for index, code in at.items():
        out[index] = code
    return out


class TestPeakCodeKwargs:
    def test_nothing_to_apply_returns_none(self):
        assert peak_code_kwargs(_CODE_FLOWS, _CODE_YEARS, None) is None
        assert peak_code_kwargs(_CODE_FLOWS, _CODE_YEARS, ["2"] * len(_CODE_YEARS)) is None

    def test_opting_out_returns_none_and_logs(self, caplog):
        kw = peak_code_kwargs(_CODE_FLOWS, _CODE_YEARS, _codes({3: "6"}), apply_peak_codes=False)
        assert kw is None
        assert "code 6: 1" in caplog.text

    def test_code_7_is_historic_and_codes_6_and_c_removed(self):
        kw = peak_code_kwargs(_CODE_FLOWS, _CODE_YEARS, _codes({0: "7", 3: "6", 5: "2,C"}))
        assert kw["historical_peaks"] == [(2001, _CODE_FLOWS[0])]
        years = list(kw["water_years"])
        assert 2001 not in years and 2004 not in years and 2006 not in years
        assert len(kw["peak_flows"]) == len(_CODE_YEARS) - 3

    def test_both_engines_are_checked(self):
        kw = peak_code_kwargs(
            _CODE_FLOWS, _CODE_YEARS, _codes({0: "7"}), engines=("native", "fortran")
        )
        assert kw is not None

    def test_skew_and_low_outlier_settings_carry_through(self):
        kw = peak_code_kwargs(
            _CODE_FLOWS,
            _CODE_YEARS,
            _codes({3: "3"}),
            regional_skew=-0.07,
            regional_skew_mse=0.18,
            user_low_outlier_threshold=120.0,
        )
        assert kw["regional_skew"] == -0.07
        assert kw["regional_skew_mse"] == 0.18
        assert kw["user_low_outlier_threshold"] == 120.0

    def test_censored_codes_become_interval_peaks(self):
        """Codes 4 and 8 are siteQT's (Qmin, q) and (q, Qmax) rows, for both engines."""
        kw = peak_code_kwargs(
            _CODE_FLOWS, _CODE_YEARS, _codes({2: "4", 4: "8"}), engines=("native", "fortran")
        )
        assert kw["interval_peaks"] == [(2003, 1e-20, _CODE_FLOWS[2]), (2005, _CODE_FLOWS[4], 1e20)]
        years = list(kw["water_years"])
        assert 2003 not in years and 2005 not in years

    def test_a_historic_interval_becomes_a_historical_interval_peak(self):
        """Code 7 with 4 is siteQT's (Qmin, q) row with dtype = 1, for both engines."""
        kw = peak_code_kwargs(
            _CODE_FLOWS, _CODE_YEARS, _codes({2: "4,7"}), engines=("native", "fortran")
        )
        assert kw["historical_interval_peaks"] == [(2003, 1e-20, _CODE_FLOWS[2])]
        assert kw["interval_peaks"] is None
        assert 2003 not in list(kw["water_years"])

    def test_a_removed_code_wins_over_a_censoring_one(self):
        """siteQT applies 6/C after 4/8, so "4,6" is removed, not censored."""
        kw = peak_code_kwargs(_CODE_FLOWS, _CODE_YEARS, _codes({2: "4,6"}))
        assert 2003 not in list(kw["water_years"])

    @pytest.mark.parametrize(
        "kwargs, match",
        [
            ({"historical_peaks": [(1990, 5.0)]}, "historical_peaks"),
            ({"perception_thresholds": {(1990, 2000): 50.0}}, "perception_thresholds"),
            ({"engines": ()}, "at least one"),
        ],
    )
    def test_calling_errors(self, kwargs, match):
        with pytest.raises(ValueError, match=match):
            peak_code_kwargs(_CODE_FLOWS, _CODE_YEARS, _codes({0: "7"}), **kwargs)

    def test_misaligned_codes_raise(self):
        with pytest.raises(ValueError, match="aligned"):
            peak_code_kwargs(_CODE_FLOWS, _CODE_YEARS, ["7"])

    def test_codes_need_water_years(self):
        with pytest.raises(ValueError, match="water_years"):
            peak_code_kwargs(_CODE_FLOWS, None, _codes({0: "7"}))


def _q100(result: dict) -> float:
    df = result["quantile_df"].set_index("Return Interval (yr)")
    return float(df.loc[100, "Flow (cfs)"])


class TestRunFFAPeakCodes:
    def test_codes_are_applied_when_supplied(self):
        plain = run_ffa(_CODE_FLOWS, _CODE_YEARS, station_skew_only=True)
        coded = run_ffa(
            _CODE_FLOWS, _CODE_YEARS, station_skew_only=True, peak_codes=_codes({3: "6"})
        )
        assert coded["error"] is None
        assert coded["parameters"]["peak_codes_applied"] == {"6": 1}
        assert coded["b17c"].results.n_peaks == len(_CODE_YEARS) - 1
        assert coded["parameters"]["mean_log"] != plain["parameters"]["mean_log"]

    def test_opting_out_matches_no_codes(self):
        plain = run_ffa(_CODE_FLOWS, _CODE_YEARS, station_skew_only=True)
        off = run_ffa(
            _CODE_FLOWS,
            _CODE_YEARS,
            station_skew_only=True,
            peak_codes=_codes({3: "6"}),
            apply_peak_codes=False,
        )
        assert off["parameters"]["peak_codes_applied"] == {}
        assert off["parameters"]["mean_log"] == plain["parameters"]["mean_log"]

    def test_codes_peakfq_ignores_change_nothing(self):
        plain = run_ffa(_CODE_FLOWS, _CODE_YEARS, station_skew_only=True)
        coded = run_ffa(
            _CODE_FLOWS, _CODE_YEARS, station_skew_only=True, peak_codes=["2"] * len(_CODE_YEARS)
        )
        assert coded["parameters"]["mean_log"] == plain["parameters"]["mean_log"]
        assert coded["parameters"]["skew_station"] == plain["parameters"]["skew_station"]

    def test_codes_may_be_an_array(self):
        coded = run_ffa(
            _CODE_FLOWS,
            _CODE_YEARS,
            station_skew_only=True,
            peak_codes=np.array(_codes({3: "6"}), dtype=object),
        )
        assert coded["parameters"]["peak_codes_applied"] == {"6": 1}

    def test_censored_code_is_fitted_as_an_interval(self):
        result = run_ffa(
            _CODE_FLOWS, _CODE_YEARS, station_skew_only=True, peak_codes=_codes({2: "4"})
        )
        assert result["error"] is None
        assert result["parameters"]["peak_codes_applied"] == {"4": 1}
        assert result["b17c"].results.n_peaks == len(_CODE_YEARS)
        assert result["b17c"].results.n_censored == 1

    def test_historic_interval_is_fitted_as_a_historic_row(self):
        result = run_ffa(
            _CODE_FLOWS, _CODE_YEARS, station_skew_only=True, peak_codes=_codes({2: "4,7"})
        )
        assert result["error"] is None
        r = result["b17c"].results
        assert r.n_peaks == len(_CODE_YEARS)
        assert r.n_historical == 1
        assert r.n_censored == 1

    def test_misaligned_codes_raise(self):
        with pytest.raises(ValueError, match="aligned"):
            run_ffa(_CODE_FLOWS, _CODE_YEARS, station_skew_only=True, peak_codes=["7"])

    def test_codes_and_thresholds_together_raise(self):
        with pytest.raises(ValueError, match="perception_thresholds"):
            run_ffa(
                _CODE_FLOWS,
                _CODE_YEARS,
                station_skew_only=True,
                peak_codes=_codes({0: "7"}),
                perception_thresholds=[
                    {"start_year": 1990, "end_year": 2000, "threshold_cfs": 50.0}
                ],
            )

    @pytest.mark.skipif(not TESTDATA_AVAILABLE, reason=SKIP_REASON)
    def test_real_regulated_record(self):
        """USGS 01426500 (HU02 WATSTORE): 50 of 102 peaks carry code 6 and are
        removed, as peakfq removes them. The CHANGELOG quotes these numbers."""
        from flowfreq.watstore import read_watstore

        frame = read_watstore(_testdata("extra_tests/HU02_WATSTORE.txt"))
        site = frame[frame["site_no"].astype(str).str.strip() == "01426500"]
        flows, years = site["peak_flow_cfs"].to_numpy(), site["water_year"].to_numpy()
        before = run_ffa(flows, years, station_skew_only=True)
        after = run_ffa(
            flows, years, station_skew_only=True, peak_codes=site["qualification_code"].tolist()
        )
        assert before["b17c"].results.n_peaks == 102
        assert after["b17c"].results.n_peaks == 52
        assert after["parameters"]["peak_codes_applied"] == {"6": 50, "7": 1}
        assert _q100(before) == pytest.approx(36336, rel=1e-3)
        assert _q100(after) == pytest.approx(38772, rel=1e-3)

    @pytest.mark.skipif(not TESTDATA_AVAILABLE, reason=SKIP_REASON)
    def test_real_less_than_record(self):
        """USGS 01362100 (HU02 WATSTORE): WY1985's code 4 peak (< 475 cfs) is an
        interval now, not a refusal. The CHANGELOG quotes these numbers."""
        from flowfreq.watstore import read_watstore

        frame = read_watstore(_testdata("extra_tests/HU02_WATSTORE.txt"))
        site = frame[frame["site_no"].astype(str).str.strip() == "01362100"]
        flows, years = site["peak_flow_cfs"].to_numpy(), site["water_year"].to_numpy()
        before = run_ffa(flows, years, station_skew_only=True, apply_peak_codes=False)
        after = run_ffa(
            flows, years, station_skew_only=True, peak_codes=site["qualification_code"].tolist()
        )
        assert after["error"] is None
        assert after["parameters"]["peak_codes_applied"] == {"4": 1}
        assert after["b17c"]._interval_peaks == [(1985, 1e-20, 475.0)]
        assert _q100(before) == pytest.approx(3153.8, rel=1e-3)
        assert _q100(after) == pytest.approx(3139.4, rel=1e-3)


class TestCompareEnginesPeakCodes:
    """The code step runs before the extension is needed, so these run anywhere."""

    def test_historic_interval_reaches_both_engines(self):
        pytest.importorskip("flowfreq.peakfqr", exc_type=ImportError)
        report = compare_engines(
            _CODE_FLOWS, _CODE_YEARS, station_skew_only=True, peak_codes=_codes({2: "4,7"})
        )
        assert report.native.n_historical == report.reference.n_historical == 1
        assert report.comparison.passed

    def test_censored_codes_reach_both_engines(self):
        pytest.importorskip("flowfreq.peakfqr", exc_type=ImportError)
        report = compare_engines(
            _CODE_FLOWS,
            _CODE_YEARS,
            station_skew_only=True,
            peak_codes=_codes({2: "4", 4: "8"}),
        )
        assert report.native.n_censored == 2
        assert report.comparison.passed

    def test_codes_and_historical_peaks_together_raise(self):
        with pytest.raises(ValueError, match="historical_peaks"):
            compare_engines(
                _CODE_FLOWS,
                _CODE_YEARS,
                station_skew_only=True,
                peak_codes=_codes({0: "7"}),
                historical_peaks=[(1990, 500.0)],
            )
