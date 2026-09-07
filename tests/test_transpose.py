"""Tests for flowfreq.transpose (drainage-area ratio transposition).

The correctness plan these follow is in ``docs/TRANSPOSITION_DESIGN.md`` §8.
Two of them are worth naming here because they catch whole classes of bug
that no eyeball review would:

- ``test_unit_ratio_is_the_identity``: an area ratio of 1.0 must return the
  donor's own numbers bitwise, for every exponent. Any mistake in the
  exponent plumbing -- wrong column, transposed arrays, an off-by-one in the
  interpolation -- breaks this while still producing plausible output
  everywhere else.
- ``test_interpolation_is_the_identity_at_published_points``: interpolating
  at an AEP the regression actually publishes must return the published
  exponent unchanged, not something a hair off it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from flowfreq.bulletin17c import Bulletin17C
from flowfreq.transpose import (
    DEFAULT_AREA_RATIO_RANGE,
    RegressionExponents,
    TransposedResults,
    transpose_frequency,
)

# A realistically shaped published set: the eight AEPs a flood regression
# normally reports, with the exponent easing down as the event gets rarer.
PUBLISHED_AEPS = [0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002]
PUBLISHED_B = [0.83, 0.81, 0.80, 0.78, 0.77, 0.76, 0.75, 0.74]

CITATION = "Example SIR 0000-0000 table 6, hydrologic region 3"


@pytest.fixture()
def exponents() -> RegressionExponents:
    return RegressionExponents(
        aeps=PUBLISHED_AEPS,
        exponents=PUBLISHED_B,
        citation=CITATION,
        region="EX-3",
        valid_area_range_sqmi=(1.0, 3000.0),
    )


@pytest.fixture(scope="module")
def donor() -> "object":
    """A real Big Sandy fit, so quantiles and limits have real shape."""
    from tests.fixtures.big_sandy import REGIONAL_SKEW, REGIONAL_SKEW_SD, SYSTEMATIC_PEAKS

    years = np.array(sorted(SYSTEMATIC_PEAKS))
    flows = np.array([SYSTEMATIC_PEAKS[y] for y in years], dtype=float)
    b17c = Bulletin17C(
        peak_flows=flows,
        water_years=years,
        regional_skew=REGIONAL_SKEW,
        regional_skew_mse=REGIONAL_SKEW_SD**2,
    )
    return b17c.run_analysis(method="ema")


class TestRegressionExponents:
    def test_requires_a_citation(self):
        with pytest.raises(ValueError, match="citation"):
            RegressionExponents(aeps=PUBLISHED_AEPS, exponents=PUBLISHED_B, citation="")

    def test_whitespace_citation_is_not_a_citation(self):
        with pytest.raises(ValueError, match="citation"):
            RegressionExponents(aeps=PUBLISHED_AEPS, exponents=PUBLISHED_B, citation="   ")

    def test_length_mismatch_rejected(self):
        with pytest.raises(ValueError, match="same length"):
            RegressionExponents(aeps=[0.5, 0.1], exponents=[0.8], citation=CITATION)

    def test_aep_outside_unit_interval_rejected(self):
        with pytest.raises(ValueError, match="strictly between 0 and 1"):
            RegressionExponents(aeps=[0.5, 1.0], exponents=[0.8, 0.8], citation=CITATION)

    def test_nonpositive_exponent_rejected(self):
        with pytest.raises(ValueError, match="finite and positive"):
            RegressionExponents(aeps=[0.5, 0.1], exponents=[0.8, 0.0], citation=CITATION)

    def test_repeated_aep_rejected(self):
        with pytest.raises(ValueError, match="repeated"):
            RegressionExponents(aeps=[0.5, 0.5], exponents=[0.8, 0.7], citation=CITATION)

    def test_stored_sorted_frequent_end_first(self, exponents):
        """Index 0 is the frequent end wherever the class indexes itself."""
        assert exponents.aeps[0] == pytest.approx(0.50)
        assert exponents.aeps[-1] == pytest.approx(0.002)
        assert np.all(np.diff(exponents.aeps) < 0)

    def test_input_order_does_not_matter(self, exponents):
        reversed_set = RegressionExponents(
            aeps=PUBLISHED_AEPS[::-1], exponents=PUBLISHED_B[::-1], citation=CITATION
        )
        np.testing.assert_array_equal(reversed_set.aeps, exponents.aeps)
        np.testing.assert_array_equal(reversed_set.exponents, exponents.exponents)

    def test_constant_requires_a_citation_too(self):
        with pytest.raises(ValueError, match="citation"):
            RegressionExponents.constant(1.0, citation="")

    def test_constant_applies_everywhere(self):
        const = RegressionExponents.constant(
            1.0, citation="assumed 1.0 for a sensitivity check, not from a regression"
        )
        assert const.is_constant
        frame = const.exponent_at([0.995, 0.5, 0.002])
        assert np.all(frame["exponent"] == 1.0)
        assert set(frame["source"]) == {"constant"}


class TestExponentInterpolation:
    def test_interpolation_is_the_identity_at_published_points(self, exponents):
        frame = exponents.exponent_at(PUBLISHED_AEPS)
        np.testing.assert_allclose(frame["exponent"].to_numpy(), PUBLISHED_B, rtol=0, atol=1e-12)
        assert set(frame["source"]) == {"published"}

    def test_pchip_is_also_the_identity_at_published_points(self, exponents):
        frame = exponents.exponent_at(PUBLISHED_AEPS, interpolation="pchip")
        np.testing.assert_allclose(frame["exponent"].to_numpy(), PUBLISHED_B, rtol=0, atol=1e-12)

    def test_interior_point_is_interpolated_and_bracketed(self, exponents):
        frame = exponents.exponent_at([0.03])  # between 0.04 and 0.02
        assert frame["source"].iloc[0] == "interpolated"
        assert 0.77 <= frame["exponent"].iloc[0] <= 0.78

    def test_linear_and_pchip_differ_immaterially(self, exponents):
        """Measured rather than asserted: the docstring claims the choice
        between them does not matter over a normal published set."""
        query = np.linspace(0.003, 0.49, 60)
        linear = exponents.exponent_at(query, interpolation="linear")["exponent"].to_numpy()
        pchip = exponents.exponent_at(query, interpolation="pchip")["exponent"].to_numpy()
        max_diff = float(np.max(np.abs(linear - pchip)))
        assert max_diff < 0.005, f"interpolants differ by {max_diff}, no longer immaterial"

    def test_unknown_interpolation_rejected(self, exponents):
        with pytest.raises(ValueError, match="unknown interpolation"):
            exponents.exponent_at([0.01], interpolation="cubic")

    def test_unknown_extrapolation_rejected(self, exponents):
        with pytest.raises(ValueError, match="unknown extrapolation"):
            exponents.exponent_at([0.01], extrapolation="wild")


class TestExponentExtrapolation:
    def test_frequent_end_is_flagged_extrapolated(self, exponents):
        frame = exponents.exponent_at([0.995, 0.90, 0.67])
        assert set(frame["source"]) == {"extrapolated"}

    def test_clamp_holds_the_endpoint(self, exponents):
        frame = exponents.exponent_at([0.995, 0.67], extrapolation="clamp")
        assert np.all(frame["exponent"] == PUBLISHED_B[0])

    def test_clamp_holds_the_rare_endpoint_too(self, exponents):
        frame = exponents.exponent_at([0.001], extrapolation="clamp")
        assert frame["exponent"].iloc[0] == pytest.approx(PUBLISHED_B[-1])

    def test_linear_extrapolation_continues_the_trend(self, exponents):
        """b rises toward the frequent end in this set, so extrapolating
        below AEP 0.5 must exceed the endpoint rather than hold it."""
        frame = exponents.exponent_at([0.90], extrapolation="linear")
        assert frame["exponent"].iloc[0] > PUBLISHED_B[0]

    def test_error_mode_refuses(self, exponents):
        with pytest.raises(ValueError, match="outside the published range"):
            exponents.exponent_at([0.995], extrapolation="error")

    def test_error_mode_allows_in_range_queries(self, exponents):
        frame = exponents.exponent_at([0.03], extrapolation="error")
        assert frame["source"].iloc[0] == "interpolated"

    def test_linear_extrapolation_to_0995_stays_under_the_default_bound(self, exponents):
        """Measured, not assumed: on a normally-sloped published set even a
        full linear extrapolation from AEP 0.5 out to 0.995 only reaches
        ~0.89, so the default (0.5, 1.0) bound does not bite here. The bound
        is a backstop against a pathological set, not a routine correction."""
        frame = exponents.exponent_at([0.995], extrapolation="linear")
        assert 0.85 < frame["exponent"].iloc[0] < 1.0
        assert not bool(frame["clamped"].iloc[0])

    def test_b_bounds_clamp_is_applied_and_flagged(self, exponents):
        """With a bound that does bite, the value is clamped and says so."""
        frame = exponents.exponent_at([0.995], extrapolation="linear", b_bounds=(0.5, 0.85))
        assert frame["exponent"].iloc[0] == pytest.approx(0.85)
        assert bool(frame["clamped"].iloc[0])

    def test_inverted_bounds_rejected(self, exponents):
        with pytest.raises(ValueError, match="inverted"):
            exponents.exponent_at([0.01], b_bounds=(1.0, 0.5))


class TestTransposeFrequency:
    def test_transposes_the_whole_standard_range_down_to_0995(self, donor, exponents):
        """The full STANDARD_AEP set, 0.995 through 0.002, comes through."""
        result = transpose_frequency(donor, 178.0, 131.0, exponents)
        assert result.quantiles["aep"].max() == pytest.approx(0.995)
        assert result.quantiles["aep"].min() == pytest.approx(0.002)
        assert len(result.quantiles) == len(donor.quantiles)
        assert np.all(np.isfinite(result.quantiles["flow_cfs"]))

    def test_arithmetic_matches_the_closed_form(self, donor, exponents):
        result = transpose_frequency(donor, 178.0, 131.0, exponents)
        ratio = 131.0 / 178.0
        expected = result.quantiles["donor_flow_cfs"] * ratio ** result.quantiles["exponent"]
        np.testing.assert_allclose(result.quantiles["flow_cfs"], expected, rtol=0, atol=0)

    def test_unit_ratio_is_the_identity(self, donor, exponents):
        """Ratio 1.0 must return the donor's numbers bitwise, whatever b is."""
        result = transpose_frequency(donor, 150.0, 150.0, exponents)
        np.testing.assert_array_equal(
            result.quantiles["flow_cfs"].to_numpy(),
            result.quantiles["donor_flow_cfs"].to_numpy(),
        )

    def test_smaller_target_gives_smaller_flows(self, donor, exponents):
        result = transpose_frequency(donor, 178.0, 131.0, exponents)
        assert np.all(result.quantiles["flow_cfs"] < result.quantiles["donor_flow_cfs"])

    def test_quantiles_stay_monotone_in_aep(self, donor, exponents):
        result = transpose_frequency(donor, 178.0, 131.0, exponents)
        ordered = result.quantiles.sort_values("aep", ascending=False)
        assert np.all(np.diff(ordered["flow_cfs"].to_numpy()) > 0)

    def test_exponent_of_one_is_the_plain_area_ratio(self, donor):
        naive = RegressionExponents.constant(1.0, citation="assumed 1.0, sensitivity check")
        result = transpose_frequency(donor, 178.0, 89.0, naive)
        np.testing.assert_allclose(
            result.quantiles["flow_cfs"],
            result.quantiles["donor_flow_cfs"] * 0.5,
            rtol=1e-12,
        )

    def test_confidence_limits_scale_by_the_same_factor(self, donor, exponents):
        result = transpose_frequency(donor, 178.0, 131.0, exponents)
        assert not result.confidence_limits.empty
        merged = result.confidence_limits.merge(
            result.quantiles[["aep", "area_ratio_factor"]], on="aep"
        )
        donor_limits = donor.confidence_limits.set_index("aep")
        for _, row in merged.iterrows():
            expected = donor_limits.loc[row["aep"], "lower_5pct"] * row["area_ratio_factor"]
            assert row["lower_5pct"] == pytest.approx(expected, rel=1e-12)

    def test_limits_are_empty_when_the_donor_has_none(self, donor, exponents):
        stripped = donor
        original = stripped.confidence_limits
        try:
            stripped.confidence_limits = pd.DataFrame()
            result = transpose_frequency(stripped, 178.0, 131.0, exponents)
            assert result.confidence_limits.empty
        finally:
            stripped.confidence_limits = original

    def test_max_aep_trims_the_frequent_end(self, donor, exponents):
        result = transpose_frequency(donor, 178.0, 131.0, exponents, max_aep=0.67)
        assert result.quantiles["aep"].max() == pytest.approx(0.67)


class TestGuardrails:
    def test_area_ratio_below_range_raises(self, donor, exponents):
        with pytest.raises(ValueError, match="outside the applicable range"):
            transpose_frequency(donor, 178.0, 70.0, exponents)  # ratio 0.39

    def test_area_ratio_above_range_raises(self, donor, exponents):
        with pytest.raises(ValueError, match="outside the applicable range"):
            transpose_frequency(donor, 100.0, 160.0, exponents)  # ratio 1.6

    @pytest.mark.parametrize("target", [50.0, 150.0])
    def test_band_endpoints_are_allowed(self, donor, exponents, target):
        result = transpose_frequency(donor, 100.0, target, exponents)
        assert result.provenance.area_ratio_in_range

    def test_override_transposes_and_records_the_violation(self, donor, exponents):
        result = transpose_frequency(donor, 178.0, 70.0, exponents, allow_out_of_range=True)
        assert not result.provenance.area_ratio_in_range
        assert "outside the applicable range" in result.to_markdown()

    def test_zero_area_rejected(self, donor, exponents):
        with pytest.raises(ValueError, match="donor_area_sqmi must be positive"):
            transpose_frequency(donor, 0.0, 100.0, exponents)

    def test_negative_target_area_rejected(self, donor, exponents):
        with pytest.raises(ValueError, match="target_area_sqmi must be positive"):
            transpose_frequency(donor, 178.0, -5.0, exponents)

    def test_donor_without_quantiles_rejected(self, exponents):
        empty = Bulletin17C(
            peak_flows=np.array([100.0, 200.0, 150.0, 300.0, 250.0]),
            water_years=np.array([2000, 2001, 2002, 2003, 2004]),
        ).run_analysis(method="mom")
        empty.quantiles = pd.DataFrame()
        with pytest.raises(ValueError, match="carries no quantiles"):
            transpose_frequency(empty, 178.0, 131.0, exponents)

    def test_max_aep_excluding_everything_rejected(self, donor, exponents):
        """Below the donor's rarest AEP (0.002), nothing survives the filter."""
        with pytest.raises(ValueError, match="excluded every quantile"):
            transpose_frequency(donor, 178.0, 131.0, exponents, max_aep=0.001)


class TestProvenanceAndReport:
    def test_extrapolated_aeps_are_listed(self, donor, exponents):
        result = transpose_frequency(donor, 178.0, 131.0, exponents)
        extrapolated = result.provenance.extrapolated_aeps
        # Everything more frequent than the published 0.50 endpoint.
        assert set(np.round(extrapolated, 3)) >= {0.995, 0.99, 0.95, 0.9, 0.8, 0.67}

    def test_citation_survives_into_provenance(self, donor, exponents):
        result = transpose_frequency(donor, 178.0, 131.0, exponents)
        assert result.provenance.citation == CITATION
        assert result.provenance.region == "EX-3"

    def test_markdown_reports_the_arithmetic_and_the_caveat(self, donor, exponents):
        result = transpose_frequency(
            donor, 178.0, 131.0, exponents, target_name="Test Creek", donor_site_no="03606500"
        )
        md = result.to_markdown()
        assert "Test Creek" in md
        assert "03606500" in md
        assert CITATION in md
        assert "extrapolated exponent" in md
        assert "understates the true uncertainty" in md
        # One row per transposed quantile, plus header and separator.
        assert md.count("\n|") >= len(result.quantiles) + 2

    def test_markdown_says_so_when_nothing_was_extrapolated(self, donor, exponents):
        result = transpose_frequency(donor, 178.0, 131.0, exponents, max_aep=0.50)
        assert "within the cited regression's published range" in result.to_markdown()

    def test_result_is_not_a_frequency_results(self, donor, exponents):
        """No fit was performed at the target, so there are no moments to
        report and the type must not pretend otherwise."""
        result = transpose_frequency(donor, 178.0, 131.0, exponents)
        assert isinstance(result, TransposedResults)
        assert not hasattr(result, "mean_log")
        assert not hasattr(result, "skew_used")


def test_default_area_ratio_range_is_the_usgs_band():
    assert DEFAULT_AREA_RATIO_RANGE == (0.5, 1.5)
