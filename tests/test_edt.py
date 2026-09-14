"""Tests for flowfreq.edt (EDT Level 2 hydrology attribute statistics).

The load-bearing test here is
``TestQ2yrRating::test_reproduces_every_published_table3_example``. The rating
transform was calibrated against D&L Table 3, so that test is the only thing
standing between the anchors and a silent regression -- and the anchors encode
a non-obvious fact (index 4 sits at the printed +110, not at its band
midpoint) that a well-meaning edit would "tidy" away.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from flowfreq.edt import (
    DIEL_NO_VARIATION_IN_PER_HR,
    DIEL_STAGE_BANDS,
    EDT_HYDROLOGY_ATTRIBUTES,
    LOW_FLOW_CHANGE_BANDS,
    LOW_FLOW_WINDOWS,
    Q2YR_CHANGE_BANDS,
    Q2YR_RATING_ANCHORS,
    Q2YR_RATING_MAX_RESIDUAL,
    TABLE3_Q2YR_EXAMPLES,
    TQMEAN_FLASHINESS_BANDS,
    TQMEAN_HIGH_FLOW_BANDS,
    audit_record,
    index_band,
    low_flow_change,
    q2yr_rating,
    split_periods,
    tqmean_change,
)

P1 = (pd.Timestamp("1960-01-01"), pd.Timestamp("1984-12-31"))
P2 = (pd.Timestamp("1985-01-01"), pd.Timestamp("2009-12-31"))


def _steady(years: int, start: str, level: float = 100.0, summer_low: float = 100.0):
    """A deterministic series: constant `level`, dropping to `summer_low` for
    a 90-day summer window each year. No noise, so the n-day minimum is
    exactly `summer_low` and TQmean is high (most days sit at one value)."""
    idx = pd.date_range(start, periods=int(365.25 * years), freq="D")
    doy = idx.dayofyear.to_numpy()
    flow = np.where((doy >= 190) & (doy < 280), summer_low, level)
    return pd.DataFrame({"flow_cfs": flow.astype(float)}, index=idx)


def _seasonal(years: int, start: str, level: float = 100.0, low: float = 10.0):
    """A non-flashy series with a real distribution: a high plateau for most of
    the year dropping to `low` for 65 days. Most days sit above the annual
    mean, so TQmean is high -- unlike `_steady`, which is constant and whose
    TQmean is legitimately 0 (no day exceeds the mean when every day equals
    it). Deterministic, so the direction of any change is by construction."""
    idx = pd.date_range(start, periods=int(365.25 * years), freq="D")
    doy = idx.dayofyear.to_numpy()
    flow = np.where((doy >= 200) & (doy < 265), low, level)
    return pd.DataFrame({"flow_cfs": flow.astype(float)}, index=idx)


def _spiky(years: int, start: str, base: float = 10.0, spike: float = 2000.0):
    """A flashy series: a low base with a spike every 14th day. Few days sit
    above the annual mean, so TQmean is low -- the opposite end of the
    flashiness scale from `_steady`, by construction rather than by chance."""
    idx = pd.date_range(start, periods=int(365.25 * years), freq="D")
    flow = np.full(len(idx), base)
    flow[::14] = spike
    return pd.DataFrame({"flow_cfs": flow}, index=idx)


class TestIndexBand:
    """Tests for index_band."""

    def test_bands_are_contiguous_and_disjoint(self) -> None:
        """A gap or an overlap in a ladder would silently misclassify."""
        for bands in (Q2YR_CHANGE_BANDS, TQMEAN_FLASHINESS_BANDS, LOW_FLOW_CHANGE_BANDS):
            ordered = sorted(bands, key=lambda b: b.low)
            for lower, upper in zip(ordered, ordered[1:]):
                assert lower.high == upper.low, f"{bands} is not contiguous"

    def test_published_band_edges(self) -> None:
        """Spot-check each ladder against the Categorical Conclusions tables."""
        assert index_band(0.0, Q2YR_CHANGE_BANDS).index == 2
        assert index_band(25.0, Q2YR_CHANGE_BANDS).index == 3
        assert index_band(-45.0, Q2YR_CHANGE_BANDS).index == 0
        assert index_band(-30.0, LOW_FLOW_CHANGE_BANDS).index == 3
        assert index_band(-60.0, LOW_FLOW_CHANGE_BANDS).index == 4
        assert index_band(80.0, LOW_FLOW_CHANGE_BANDS).index == 0
        assert index_band(-10.0, TQMEAN_FLASHINESS_BANDS).index == 3
        assert index_band(20.0, TQMEAN_FLASHINESS_BANDS).index == 0

    def test_lower_edge_is_inclusive(self) -> None:
        """Bands are half-open, so the printed boundary lands in the upper band."""
        assert index_band(20.0, Q2YR_CHANGE_BANDS).index == 3
        assert index_band(19.999, Q2YR_CHANGE_BANDS).index == 2

    def test_diel_ladder_is_in_stage_units_not_percent(self) -> None:
        """Attribute 17's ladder is inches of stage per hour, and index 0 -- not
        index 2 -- is its pristine reference."""
        assert index_band(1.0, DIEL_STAGE_BANDS).index == 1
        assert index_band(4.0, DIEL_STAGE_BANDS).index == 2
        assert index_band(18.0, DIEL_STAGE_BANDS).index == 4

    def test_provisional_edge_is_flagged(self) -> None:
        """The diel 0/1 boundary is a local convention; saying so is the point."""
        result = index_band(DIEL_NO_VARIATION_IN_PER_HR / 2, DIEL_STAGE_BANDS)
        assert result.index == 0
        assert any("local convention" in w for w in result.warnings)

    def test_beyond_printed_edge_is_flagged(self) -> None:
        """A value past what D&L printed is still classified, but not silently."""
        result = index_band(300.0, Q2YR_CHANGE_BANDS)
        assert result.index == 4
        assert any("beyond the outer edge" in w for w in result.warnings)

    def test_within_printed_edge_is_not_flagged(self) -> None:
        assert index_band(50.0, Q2YR_CHANGE_BANDS).warnings == ()

    def test_tqmean_high_flow_ladder_cannot_reach_index_0_or_1(self) -> None:
        """Attribute 15's TQmean option bottoms out at index 2.

        Indices 0 and 1 are defined by a decrease in Q2yr or by known
        regulation, not by TQmean, so a reach whose peaks are suppressed
        cannot be rated from TQmean alone. An increase in TQmean is index 2,
        because the published cell reads "<5% reduction" and an increase
        satisfies that as written.
        """
        assert {b.index for b in TQMEAN_HIGH_FLOW_BANDS} == {2, 3, 4}
        assert index_band(-8.0, TQMEAN_HIGH_FLOW_BANDS).index == 3
        assert index_band(40.0, TQMEAN_HIGH_FLOW_BANDS).index == 2

    def test_non_finite_raises(self) -> None:
        for bad in (np.nan, np.inf, None):
            with pytest.raises(ValueError, match="non-finite"):
                index_band(bad, Q2YR_CHANGE_BANDS)


class TestQ2yrRating:
    """Tests for q2yr_rating, the one calibrated transform in this module."""

    def test_reproduces_every_published_table3_example(self) -> None:
        """The calibration test. D&L Table 3 is the only published set of
        percent-change/rating pairs in the hydrology sections, so it is both
        where the anchors came from and the only check on them."""
        for name, pct, published in TABLE3_Q2YR_EXAMPLES:
            got = q2yr_rating(pct)
            assert abs(got - published) <= Q2YR_RATING_MAX_RESIDUAL, (
                f"{name}: {pct:+.1f}% gave {got:.3f}, D&L Table 3 publishes " f"{published}"
            )

    def test_index_4_anchors_at_the_printed_edge_not_the_midpoint(self) -> None:
        """The non-obvious anchor, asserted so a tidy-up cannot remove it.

        Every other index sits at its band's midpoint, per the footnote to
        D&L 36.6/37.6. Index 4 sits at the +110 its cell prints. Using the
        midpoint (+75) instead misses Table 3 by more than twenty times the
        tolerance, so this is a fact about the source, not a preference.
        """
        assert (110.0, 4.0) in Q2YR_RATING_ANCHORS
        midpoint_anchors = [(0.0, 2.0), (30.0, 3.0), (75.0, 4.0)]
        xs = [a[0] for a in midpoint_anchors]
        ys = [a[1] for a in midpoint_anchors]
        worst = max(
            abs(float(np.interp(pct, xs, ys)) - published)
            for _, pct, published in TABLE3_Q2YR_EXAMPLES
            if pct > 0
        )
        assert worst > 10 * Q2YR_RATING_MAX_RESIDUAL

    def test_zero_change_is_the_pristine_rating(self) -> None:
        assert q2yr_rating(0.0) == pytest.approx(2.0)

    def test_monotonic_in_percent_change(self) -> None:
        values = [q2yr_rating(p) for p in range(-100, 200, 5)]
        assert all(b >= a for a, b in zip(values, values[1:]))

    def test_clamped_to_the_index_range(self) -> None:
        """Past the outermost anchors it clamps; an index cannot exceed 4."""
        assert q2yr_rating(10_000.0) == 4.0
        assert q2yr_rating(-99.9) == 0.0

    def test_band_floor_disagrees_with_the_rating(self) -> None:
        """The trap this module documents: flooring the band is not the rating.

        Six of the eight published examples disagree, so anything that reports
        a band where a rating is expected is wrong by about one index unit.
        """
        disagreements = sum(
            int(np.floor(published)) != index_band(pct, Q2YR_CHANGE_BANDS).index
            for _, pct, published in TABLE3_Q2YR_EXAMPLES
        )
        assert disagreements == 6

    def test_non_finite_raises(self) -> None:
        with pytest.raises(ValueError, match="finite"):
            q2yr_rating(np.nan)


class TestSplitPeriods:
    """Tests for split_periods."""

    def test_default_splits_at_the_midpoint(self) -> None:
        data = _steady(40, "1960-01-01")
        (s1, e1), (s2, e2) = split_periods(data)
        assert s1 == data.index.min()
        assert e2 == data.index.max()
        assert e1 == s2 - pd.Timedelta(days=1)
        assert abs((e1 - s1).days - (e2 - s2).days) <= 2

    def test_explicit_split_honoured(self) -> None:
        data = _steady(40, "1960-01-01")
        (_, e1), (s2, _) = split_periods(data, "1980-01-01")
        assert s2 == pd.Timestamp("1980-01-01")
        assert e1 == pd.Timestamp("1979-12-31")

    def test_split_outside_record_raises(self) -> None:
        data = _steady(10, "1960-01-01")
        with pytest.raises(ValueError, match="outside the record"):
            split_periods(data, "2050-01-01")

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="empty record"):
            split_periods(pd.DataFrame({"flow_cfs": []}, index=pd.DatetimeIndex([])))


class TestAuditRecord:
    """Tests for audit_record."""

    def test_long_record_with_split_is_sufficient_for_all_three(self) -> None:
        audit = audit_record(_steady(50, "1960-01-01"), split="1985-01-01")
        assert set(audit["attribute_id"]) == {"15", "16", "18"}
        assert audit["sufficient"].all()

    def test_short_record_fails_the_forty_year_statistics_only(self) -> None:
        """A 25-year record cannot support Q2yr or the low-flow change (both
        need ~40 yr with 20 per state) but is long enough for TQmean, which
        needs ~10 per state. This asymmetry is the whole reason TQmean is the
        viable route on a thin record."""
        audit = audit_record(_steady(25, "1960-01-01"), split="1972-07-01").set_index(
            "attribute_id"
        )
        assert not audit.loc["15", "sufficient"]
        assert not audit.loc["16", "sufficient"]
        assert audit.loc["18", "sufficient"]

    def test_without_a_split_the_per_state_clause_is_unchecked(self) -> None:
        """No development date means the 20-years-per-state condition cannot be
        tested, so nothing is reported sufficient and the note says why."""
        audit = audit_record(_steady(60, "1950-01-01"))
        assert not audit["sufficient"].any()
        assert audit["years_period1"].isna().all()
        assert all("unchecked" in note for note in audit["note"])

    def test_years_measured_from_valid_days_only(self) -> None:
        """Negative values are missing, so they must not extend the record."""
        data = _steady(20, "1960-01-01")
        data.loc[data.index[-int(365.25 * 5) :], "flow_cfs"] = -999.0
        audit = audit_record(data)
        assert audit["years_available"].iloc[0] == pytest.approx(15.0, abs=0.2)

    def test_missing_column_raises(self) -> None:
        with pytest.raises(KeyError, match="flow_cfs"):
            audit_record(pd.DataFrame({"q": [1.0]}, index=pd.date_range("2000-01-01", periods=1)))

    def test_all_invalid_raises(self) -> None:
        data = _steady(5, "1960-01-01")
        data["flow_cfs"] = -1.0
        with pytest.raises(ValueError, match="no valid"):
            audit_record(data)


class TestTqmeanChange:
    """Tests for tqmean_change."""

    def test_steady_to_spiky_is_a_tqmean_reduction(self) -> None:
        """Constructed so the direction cannot be accidental: a steady series
        has most days above its own mean (high TQmean), a spiked one has few
        (low TQmean). Lower TQmean is flashier, which is index 3 or 4."""
        data = pd.concat([_seasonal(25, "1960-01-01"), _spiky(25, "1985-01-01")])
        change = tqmean_change(data, periods=(P1, P2))
        assert change.value1 > change.value2
        assert change.pct_change < -15.0
        assert change.band is not None and change.band.index == 4
        assert change.attribute_id == "18"

    def test_unchanged_record_rates_as_pristine(self) -> None:
        data = pd.concat([_seasonal(25, "1960-01-01"), _seasonal(25, "1985-01-01")])
        change = tqmean_change(data, periods=(P1, P2))
        assert change.pct_change == pytest.approx(0.0, abs=2.0)
        assert change.band is not None and change.band.index == 2

    def test_high_flow_ladder_selects_attribute_15(self) -> None:
        data = pd.concat([_seasonal(25, "1960-01-01"), _spiky(25, "1985-01-01")])
        change = tqmean_change(data, periods=(P1, P2), ladder="high_flow")
        assert change.attribute_id == "15"
        assert change.band is not None and change.band.index == 4

    def test_three_window_constructions_warn_differently(self) -> None:
        """A halving, a stated split, and explicit windows are three different
        things, and an earlier cut of this told a caller who had supplied a
        development-state date that their periods "were derived by halving the
        record". Each case gets its own message, or none."""
        data = pd.concat([_seasonal(25, "1960-01-01"), _spiky(25, "1985-01-01")])

        halved = tqmean_change(data).warnings
        assert any("halving the record" in w for w in halved)

        split_at = tqmean_change(data, split="1985-01-01").warnings
        assert any("splitting at 1985-01-01" in w for w in split_at)
        assert not any("halving" in w for w in split_at)

        explicit = tqmean_change(data, periods=(P1, P2)).warnings
        assert not any("halving" in w or "splitting at" in w for w in explicit)

    def test_short_period_warns_but_still_computes(self) -> None:
        data = pd.concat([_seasonal(4, "1960-01-01"), _spiky(4, "1964-01-01")])
        change = tqmean_change(
            data,
            periods=(
                (pd.Timestamp("1960-01-01"), pd.Timestamp("1963-12-31")),
                (pd.Timestamp("1964-01-01"), pd.Timestamp("1967-12-31")),
            ),
        )
        assert np.isfinite(change.pct_change)
        assert any("complete years in the shorter period" in w for w in change.warnings)

    def test_empty_period_gives_nan_and_says_so(self) -> None:
        data = _seasonal(25, "1960-01-01")
        change = tqmean_change(
            data,
            periods=(P1, (pd.Timestamp("2030-01-01"), pd.Timestamp("2035-01-01"))),
        )
        assert np.isnan(change.pct_change)
        assert change.band is None
        assert any("record-length failure" in w for w in change.warnings)

    def test_bad_ladder_raises(self) -> None:
        with pytest.raises(ValueError, match="ladder must be one of"):
            tqmean_change(_steady(25, "1960-01-01"), ladder="peak")

    def test_missing_column_raises(self) -> None:
        with pytest.raises(KeyError, match="flow_cfs"):
            tqmean_change(pd.DataFrame({"q": [1.0]}, index=pd.date_range("2000-01-01", periods=1)))

    def test_to_row_is_flat_and_carries_the_band(self) -> None:
        data = pd.concat([_seasonal(25, "1960-01-01"), _spiky(25, "1985-01-01")])
        row = tqmean_change(data, periods=(P1, P2)).to_row()
        assert row["attribute_id"] == "18"
        assert row["index_band"] == 4
        assert isinstance(row["pct_change"], float)
        assert not any(isinstance(v, (dict, list, tuple)) for v in row.values())


class TestLowFlowChange:
    """Tests for low_flow_change."""

    def test_known_reduction_recovered_exactly(self) -> None:
        """The summer low is a constant by construction, so a 100 -> 60 cfs
        change is exactly -40% and must land in band 3 (-20 to -50%)."""
        data = pd.concat(
            [
                _steady(25, "1960-01-01", summer_low=100.0),
                _steady(25, "1985-01-01", summer_low=60.0),
            ]
        )
        change = low_flow_change(data, n_day=45, periods=(P1, P2))
        assert change.value1 == pytest.approx(100.0, abs=0.5)
        assert change.value2 == pytest.approx(60.0, abs=0.5)
        assert change.pct_change == pytest.approx(-40.0, abs=1.0)
        assert change.band is not None and change.band.index == 3
        assert change.attribute_id == "16"

    def test_severe_reduction_reaches_band_4(self) -> None:
        data = pd.concat(
            [
                _steady(25, "1960-01-01", summer_low=100.0),
                _steady(25, "1985-01-01", summer_low=30.0),
            ]
        )
        change = low_flow_change(data, n_day=45, periods=(P1, P2))
        assert change.pct_change == pytest.approx(-70.0, abs=1.0)
        assert change.band is not None and change.band.index == 4

    def test_increase_rates_at_the_low_index_end(self) -> None:
        """Regulation that stabilises or raises low flows is index 0 or 1, not
        a better version of index 2 -- the ladder runs in both directions."""
        data = pd.concat(
            [
                _steady(25, "1960-01-01", summer_low=50.0),
                _steady(25, "1985-01-01", summer_low=100.0),
            ]
        )
        change = low_flow_change(data, n_day=45, periods=(P1, P2))
        assert change.pct_change == pytest.approx(100.0, abs=2.0)
        assert change.band is not None and change.band.index == 0

    def test_both_published_windows_run(self) -> None:
        data = pd.concat(
            [
                _steady(25, "1960-01-01", summer_low=100.0),
                _steady(25, "1985-01-01", summer_low=60.0),
            ]
        )
        for window in LOW_FLOW_WINDOWS:
            change = low_flow_change(data, n_day=window, periods=(P1, P2))
            assert change.pct_change == pytest.approx(-40.0, abs=2.0)
            assert not any("not one of" in w for w in change.warnings)

    def test_unpublished_window_is_flagged_as_a_sensitivity_check(self) -> None:
        data = pd.concat(
            [
                _steady(25, "1960-01-01", summer_low=100.0),
                _steady(25, "1985-01-01", summer_low=60.0),
            ]
        )
        change = low_flow_change(data, n_day=7, periods=(P1, P2))
        assert any("sensitivity check" in w for w in change.warnings)

    def test_non_positive_window_raises(self) -> None:
        with pytest.raises(ValueError, match="n_day must be >= 1"):
            low_flow_change(_steady(10, "1960-01-01"), n_day=0)

    def test_missing_column_raises(self) -> None:
        with pytest.raises(KeyError, match="flow_cfs"):
            low_flow_change(
                pd.DataFrame({"q": [1.0]}, index=pd.date_range("2000-01-01", periods=1))
            )


class TestAttributeRegistry:
    """Tests for EDT_HYDROLOGY_ATTRIBUTES."""

    def test_retired_attribute_is_present_and_says_so(self) -> None:
        """37 is kept so a caller reaching for it learns it is retired rather
        than finding nothing and assuming an omission."""
        spec = EDT_HYDROLOGY_ATTRIBUTES["37"]
        assert "RETIRED" in spec.name
        assert "retired" in spec.statistic.lower()

    def test_diel_and_natural_regime_need_no_reference_period(self) -> None:
        """17's pristine value is fixed at index 0 a priori and 36 is the
        regime itself, so neither needs a paired template run."""
        assert not EDT_HYDROLOGY_ATTRIBUTES["17"].needs_reference_period
        assert not EDT_HYDROLOGY_ATTRIBUTES["36"].needs_reference_period
        for attribute_id in ("15", "16", "18"):
            assert EDT_HYDROLOGY_ATTRIBUTES[attribute_id].needs_reference_period

    def test_tqmean_record_requirement_is_shorter_than_q2yr(self) -> None:
        """The fact that decides which statistic a thin record can support."""
        assert EDT_HYDROLOGY_ATTRIBUTES["18"].min_years_total == 10.0
        assert EDT_HYDROLOGY_ATTRIBUTES["15"].min_years_total == 40.0

    def test_contested_rating_month_is_recorded_not_resolved(self) -> None:
        """18.8 says the highest-runoff month; Appendix A says the lowest-flow
        month. That is a contradiction in the source, and the registry must
        surface it rather than silently pick one."""
        assert "CONTESTED" in EDT_HYDROLOGY_ATTRIBUTES["18"].rating_month

    def test_every_flow_attribute_is_shaped(self) -> None:
        for attribute_id in ("15", "16", "17", "18"):
            assert EDT_HYDROLOGY_ATTRIBUTES[attribute_id].shaped
        assert not EDT_HYDROLOGY_ATTRIBUTES["36"].shaped


class TestSplitApplicability:
    """A single split date does not suit every record in a batch.

    Regression tests for what a first 13-gage run actually produced: a split at
    the very start of a short record came back as "0.0 years" and a NaN change,
    with nothing naming the cause. The fix is in the audit rather than in the
    slicing -- how short is too short is a sufficiency judgement.
    """

    def test_degenerate_period_is_reported_not_raised(self) -> None:
        """A record starting days before the split slices fine and answers
        nothing, so it is flagged rather than erroring."""
        data = _steady(20, "2000-01-05")
        audit = audit_record(data, split="2000-01-08")
        assert not audit["split_applicable"].any()
        assert not audit["sufficient"].any()
        assert all("does not suit this record" in note for note in audit["note"])

    def test_the_note_names_the_records_own_range(self) -> None:
        """So a batch run can tell which gage the split does not suit, and pick
        another date for it."""
        audit = audit_record(_steady(20, "2000-01-05"), split="2000-01-08")
        assert all("2000-01-05" in note for note in audit["note"])

    def test_applicable_split_is_marked_applicable(self) -> None:
        audit = audit_record(_steady(50, "1960-01-01"), split="1985-01-01")
        assert audit["split_applicable"].all()
        assert audit["sufficient"].all()

    def test_no_split_is_not_applicable_either(self) -> None:
        """Absent a development date the clause is unchecked, which is not the
        same as satisfied."""
        audit = audit_record(_steady(50, "1960-01-01"))
        assert not audit["split_applicable"].any()

    def test_applicable_but_insufficient_is_distinguishable(self) -> None:
        """The two failure modes need different remedies, so they are separate
        columns: a different split date versus more data."""
        audit = audit_record(_steady(25, "1960-01-01"), split="1972-07-01").set_index(
            "attribute_id"
        )
        assert audit.loc["15", "split_applicable"]
        assert not audit.loc["15", "sufficient"]

    def test_split_outside_the_record_still_raises(self) -> None:
        with pytest.raises(ValueError, match="outside the record"):
            split_periods(_steady(10, "1960-01-01"), "2050-01-01")
