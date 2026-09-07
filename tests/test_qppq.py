"""Tests for flowfreq.qppq.

The one worth reading is ``TestMeltTimingPhaseError``. It pins a measured
finding that contradicts the obvious intuition: grouping the duration curves
by coarse season does not fix a melt-timing offset between donor and target,
it makes it *worse*, because the offset happens inside a four-month season
block. Aligning the timing with a lag is what works. If a future change makes
seasonal grouping look like the fix again, that test should fail.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from flowfreq.qppq import (
    MONTHLY_SEASONS,
    SNOWMELT_SEASONS,
    FlowDurationCurve,
    apply_lag,
    center_of_timing,
    estimate_donor_lag,
    loocv_qppq,
    performance,
    qppq,
    rank_donors,
    seasonal_curves,
)


def synthetic_snowmelt(seed: int, scale: float, phase_days: int = 0, years: int = 8):
    """A snowmelt hydrograph: a spring freshet on a low winter baseflow.

    ``phase_days`` shifts the melt later, which is what a higher, colder
    basin does relative to a lower one.
    """
    index = pd.date_range("2000-10-01", periods=365 * years, freq="D")
    doy = index.dayofyear.to_numpy() + phase_days
    melt = np.exp(-(((doy - 150) % 365) ** 2) / (2 * 45.0**2))
    rng = np.random.default_rng(seed)
    flows = scale * (0.15 + 3.0 * melt) * np.exp(rng.normal(0, 0.25, len(index)))
    return pd.DataFrame({"flow_cfs": flows}, index=index)


@pytest.fixture(scope="module")
def donor() -> pd.DataFrame:
    return synthetic_snowmelt(1, 100.0)


@pytest.fixture(scope="module")
def target_aligned() -> pd.DataFrame:
    return synthetic_snowmelt(2, 60.0)


@pytest.fixture(scope="module")
def target_late() -> pd.DataFrame:
    """Melts 25 days later than the donor."""
    return synthetic_snowmelt(3, 60.0, phase_days=-25)


class TestFlowDurationCurve:
    def test_round_trip_is_the_identity(self, donor):
        """A flow mapped to its probability and back must come back
        unchanged. Everything else in this module rests on it."""
        curve = FlowDurationCurve.from_daily(donor, "donor")
        recovered = curve.flow_at(curve.probability_of(curve.flows))
        np.testing.assert_allclose(recovered, curve.flows, rtol=1e-9)

    def test_curve_decreases_with_exceedance(self, donor):
        curve = FlowDurationCurve.from_daily(donor)
        flows = curve.flow_at([0.01, 0.1, 0.5, 0.9, 0.99])
        assert np.all(np.diff(flows) < 0)

    def test_negative_discharge_rejected(self):
        frame = pd.DataFrame(
            {"flow_cfs": [10.0, -1.0, 5.0]}, index=pd.date_range("2000-01-01", periods=3)
        )
        with pytest.raises(ValueError, match="negative discharge"):
            FlowDurationCurve.from_daily(frame)

    def test_all_zero_record_rejected(self):
        frame = pd.DataFrame({"flow_cfs": [0.0, 0.0]}, index=pd.date_range("2000-01-01", periods=2))
        with pytest.raises(ValueError, match="entirely zero flow"):
            FlowDurationCurve.from_daily(frame)

    def test_intermittent_record_reports_zero_fraction(self):
        flows = [0.0] * 20 + list(np.linspace(1.0, 50.0, 80))
        frame = pd.DataFrame({"flow_cfs": flows}, index=pd.date_range("2000-01-01", periods=100))
        curve = FlowDurationCurve.from_daily(frame)
        assert curve.is_intermittent
        assert curve.p_zero == pytest.approx(0.20)

    def test_zero_block_returns_zero_not_a_small_number(self):
        flows = [0.0] * 20 + list(np.linspace(1.0, 50.0, 80))
        frame = pd.DataFrame({"flow_cfs": flows}, index=pd.date_range("2000-01-01", periods=100))
        curve = FlowDurationCurve.from_daily(frame)
        assert curve.flow_at([0.95])[0] == 0.0

    def test_probability_outside_unit_interval_rejected(self, donor):
        curve = FlowDurationCurve.from_daily(donor)
        with pytest.raises(ValueError, match="must lie in"):
            curve.flow_at([1.5])

    def test_missing_column_raises(self):
        with pytest.raises(KeyError, match="flow_cfs"):
            FlowDurationCurve.from_daily(pd.DataFrame({"q": [1.0]}))


class TestSeasonalCurves:
    def test_one_curve_per_season(self, donor):
        curves = seasonal_curves(donor, label="donor")
        assert set(curves) == set(SNOWMELT_SEASONS)

    def test_monthly_grouping(self, donor):
        curves = seasonal_curves(donor, MONTHLY_SEASONS)
        assert len(curves) == 12

    def test_freshet_is_wetter_than_recession(self, donor):
        curves = seasonal_curves(donor)
        assert curves["freshet"].flow_at([0.5])[0] > curves["recession"].flow_at([0.5])[0]

    def test_incomplete_grouping_rejected(self, donor):
        with pytest.raises(ValueError, match="do not cover month"):
            seasonal_curves(donor, {"half": (1, 2, 3, 4, 5, 6)})

    def test_overlapping_grouping_rejected(self, donor):
        with pytest.raises(ValueError, match="in both"):
            seasonal_curves(donor, {"a": tuple(range(1, 13)), "b": (6,)})


class TestQppq:
    def test_transfers_a_series(self, donor, target_aligned):
        result = qppq(
            donor,
            FlowDurationCurve.from_daily(donor),
            FlowDurationCurve.from_daily(target_aligned),
        )
        assert len(result.series) == len(donor)
        assert np.all(np.isfinite(result.series["flow_cfs"]))

    def test_transferring_through_the_donors_own_curve_is_the_identity(self, donor):
        """Donor -> probability -> donor must return the donor."""
        curve = FlowDurationCurve.from_daily(donor)
        result = qppq(donor, curve, curve)
        np.testing.assert_allclose(
            result.series["flow_cfs"], result.series["donor_flow_cfs"], rtol=1e-9
        )

    def test_seasonal_and_annual_curves_cannot_be_mixed(self, donor, target_aligned):
        with pytest.raises(ValueError, match="both be single curves"):
            qppq(
                donor,
                seasonal_curves(donor),
                FlowDurationCurve.from_daily(target_aligned),
            )

    def test_seasonal_transfer_labels_each_day(self, donor, target_aligned):
        result = qppq(donor, seasonal_curves(donor), seasonal_curves(target_aligned))
        assert set(result.series["season"]) == set(SNOWMELT_SEASONS)
        assert result.seasonal

    def test_annual_transfer_carries_the_phase_error_note(self, donor, target_aligned):
        result = qppq(
            donor,
            FlowDurationCurve.from_daily(donor),
            FlowDurationCurve.from_daily(target_aligned),
        )
        assert any("phase error" in note for note in result.notes)


class TestMeltTimingPhaseError:
    """The measured finding, pinned.

    A 25-day melt offset between donor and target is the failure mode QPPQ
    has in a snowmelt basin. These tests record what does and does not fix
    it, because the intuitive answer (group by season) is wrong.
    """

    @staticmethod
    def _score(donor, target, donor_curve, target_curve, seasons=None):
        result = qppq(donor, donor_curve, target_curve, seasons=seasons)
        paired = pd.concat(
            [
                target["flow_cfs"].rename("observed"),
                result.series["flow_cfs"].rename("estimated"),
            ],
            axis=1,
            join="inner",
        ).dropna()
        return performance(paired["observed"], paired["estimated"])

    def test_phase_error_degrades_the_transfer(self, donor, target_aligned, target_late):
        aligned = self._score(
            donor,
            target_aligned,
            FlowDurationCurve.from_daily(donor),
            FlowDurationCurve.from_daily(target_aligned),
        )
        late = self._score(
            donor,
            target_late,
            FlowDurationCurve.from_daily(donor),
            FlowDurationCurve.from_daily(target_late),
        )
        assert aligned["log_nse"] > 0.8
        assert late["log_nse"] < 0.5, "a 25-day offset should visibly break the transfer"

    def test_coarse_seasons_do_not_help_and_can_hurt(self, donor, target_late):
        """The counterintuitive one. Three four-month seasons are longer
        than the offset being corrected, so both sites stay labelled
        'freshet' and the narrower within-season curve amplifies the
        mismatch instead of absorbing it."""
        annual = self._score(
            donor,
            target_late,
            FlowDurationCurve.from_daily(donor),
            FlowDurationCurve.from_daily(target_late),
        )
        three_season = self._score(
            donor, target_late, seasonal_curves(donor), seasonal_curves(target_late)
        )
        assert three_season["log_nse"] <= annual["log_nse"]

    def test_monthly_grouping_helps(self, donor, target_late):
        annual = self._score(
            donor,
            target_late,
            FlowDurationCurve.from_daily(donor),
            FlowDurationCurve.from_daily(target_late),
        )
        monthly = self._score(
            donor,
            target_late,
            seasonal_curves(donor, MONTHLY_SEASONS),
            seasonal_curves(target_late, MONTHLY_SEASONS),
            seasons=MONTHLY_SEASONS,
        )
        assert monthly["log_nse"] > annual["log_nse"]

    def test_lag_alignment_is_what_actually_fixes_it(self, donor, target_aligned, target_late):
        """Shifting the donor by the estimated lag should recover close to
        the skill the no-offset pair gets."""
        best_lag = int(estimate_donor_lag(target_late, donor).iloc[0]["lag_days"])
        assert 20 <= best_lag <= 30, f"should recover the imposed 25-day offset, got {best_lag}"

        lagged = apply_lag(donor, best_lag)
        corrected = self._score(
            lagged,
            target_late,
            FlowDurationCurve.from_daily(lagged),
            FlowDurationCurve.from_daily(target_late),
        )
        uncorrected = self._score(
            donor,
            target_late,
            FlowDurationCurve.from_daily(donor),
            FlowDurationCurve.from_daily(target_late),
        )
        assert corrected["log_nse"] > uncorrected["log_nse"] + 0.4
        assert corrected["log_nse"] > 0.8


class TestTimingTools:
    def test_center_of_timing_is_in_the_melt_season(self, donor):
        table = center_of_timing(donor)
        mean_ct = table.loc[table["complete"], "center_of_timing_dowy"].mean()
        # Water year starts Oct 1; a late-May freshet is around day 230-280.
        assert 200 < mean_ct < 320

    def test_center_of_timing_detects_the_offset(self, donor, target_late):
        donor_ct = center_of_timing(donor)
        target_ct = center_of_timing(target_late)
        difference = (
            target_ct.loc[target_ct["complete"], "center_of_timing_dowy"].mean()
            - donor_ct.loc[donor_ct["complete"], "center_of_timing_dowy"].mean()
        )
        assert 15 < difference < 35, "should see roughly the imposed 25-day offset"

    def test_incomplete_years_are_flagged_not_computed(self):
        index = pd.date_range("2000-10-01", periods=100, freq="D")
        frame = pd.DataFrame({"flow_cfs": np.ones(100)}, index=index)
        table = center_of_timing(frame)
        assert not table["complete"].any()
        assert table["center_of_timing_dowy"].isna().all()

    def test_estimate_donor_lag_returns_the_whole_profile(self, donor, target_late):
        """Not just the argmax: a flat or double-peaked profile means the
        lag is not identified, and that has to be visible."""
        profile = estimate_donor_lag(target_late, donor, max_lag_days=30)
        assert len(profile) == 61
        assert set(profile.columns) == {"lag_days", "spearman", "n_concurrent"}

    def test_apply_lag_shifts_the_index(self, donor):
        shifted = apply_lag(donor, 10)
        assert (shifted.index[0] - donor.index[0]).days == 10


class TestDonorRankingAndScoring:
    def test_ranks_the_better_correlated_donor_first(self, donor, target_aligned, target_late):
        ranked = rank_donors(
            target_aligned, {"aligned_twin": donor, "offset": apply_lag(donor, 60)}
        )
        assert ranked.iloc[0]["donor"] == "aligned_twin"

    def test_thin_overlap_is_flagged_not_dropped(self, donor):
        short = donor.iloc[:10]
        ranked = rank_donors(short, {"donor": donor}, min_concurrent_days=365)
        assert not bool(ranked.iloc[0]["sufficient"])
        assert ranked.iloc[0]["n_concurrent"] == 10

    def test_performance_of_a_perfect_estimate(self):
        observed = np.array([1.0, 5.0, 10.0, 50.0, 100.0])
        metrics = performance(observed, observed)
        assert metrics["nse"] == pytest.approx(1.0)
        assert metrics["log_nse"] == pytest.approx(1.0)
        assert metrics["kge"] == pytest.approx(1.0)
        assert metrics["pbias"] == pytest.approx(0.0)

    def test_performance_shape_mismatch_raises(self):
        with pytest.raises(ValueError, match="shape mismatch"):
            performance([1.0, 2.0], [1.0])

    def test_log_nse_and_nse_disagree_on_low_flow_error(self):
        """Why all of them are reported. This estimate is exact at high flow
        and wrong by up to a factor of ten at low flow -- the regime that
        limits habitat. NSE calls it near-perfect; log-NSE is materially
        worse, and the dry-end bias is the one that actually shows the
        problem for what it is."""
        observed = np.array([1.0, 2.0, 5.0, 100.0, 500.0])
        estimated = np.array([10.0, 12.0, 14.0, 100.0, 500.0])
        metrics = performance(observed, estimated)

        assert metrics["nse"] > 0.99, "NSE is blind to this"
        assert metrics["nse"] - metrics["log_nse"] > 0.3, "log-NSE should see it"
        assert metrics["pbias_dry"] > 100.0, "dry-end bias should be unmissable"


class TestLoocv:
    def test_scores_every_site(self, donor, target_aligned):
        table = loocv_qppq({"a": donor, "b": target_aligned}, seasonal=False)
        assert set(table["site"]) == {"a", "b"}
        assert {"log_nse", "kge", "pbias_dry"} <= set(table.columns)

    def test_needs_two_sites(self, donor):
        with pytest.raises(ValueError, match="at least two sites"):
            loocv_qppq({"only": donor})

    def test_records_which_donor_was_used(self, donor, target_aligned):
        table = loocv_qppq({"a": donor, "b": target_aligned}, seasonal=False)
        assert set(table["donor"]) <= {"a", "b"}
        assert not (table["site"] == table["donor"]).any(), "a site cannot donate to itself"
