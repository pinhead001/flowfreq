"""Regulation guardrail on transfers (roadmap §5.2, TODO tier 2 #6).

A regulated donor carries its regulation into every transfer, and an
unregulated donor carries none into a regulated target. ``transpose_*`` and
``qppq`` refuse either mismatch unless forced, and record what they were told
in provenance. ``known_diversions`` marks a site regulated when GAGES-II does
not see the diversion (Beaver Creek in the Methow).
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import pytest

from flowfreq.core import LowFlowResults
from flowfreq.qppq import FlowDurationCurve, qppq
from flowfreq.regime import flow_duration_curve
from flowfreq.regulation import (
    RegulatedRecordError,
    RegulationMismatchError,
    SiteClass,
    SiteClassification,
    check_transfer_regulation,
    classify_site,
    require_unregulated,
)
from flowfreq.transpose import (
    RegressionExponents,
    regression_weighted_transfer,
    transpose_duration,
    transpose_frequency,
    transpose_low_flow,
)

REFERENCE_SITE = "12048000"  # Dungeness River near Sequim, WA (GAGES-II Ref)
UNKNOWN_SITE = "03606500"  # Big Sandy at Bruceton: Non-ref, below both thresholds
ABSENT_SITE = "99999999"
DIVERSIONS = "irrigation diversions above the gage (2002 canal inventory)"


def _cls(site_class: SiteClass, site_no: str = "00000000") -> SiteClassification:
    return SiteClassification(site_no=site_no, site_class=site_class, basis=("test",))


REGULATED = _cls(SiteClass.REGULATED, "11111111")
REFERENCE = _cls(SiteClass.REFERENCE, "22222222")


class TestKnownDiversions:
    def test_reference_gage_becomes_regulated(self):
        out = classify_site(REFERENCE_SITE, known_diversions=DIVERSIONS)
        assert out.site_class is SiteClass.REGULATED
        assert out.caller_regulated
        assert f"caller: {DIVERSIONS}" in out.basis

    def test_gage_outside_gagesii_becomes_regulated(self):
        out = classify_site(ABSENT_SITE, known_diversions=DIVERSIONS)
        assert out.site_class is SiteClass.REGULATED
        assert out.caller_regulated
        assert not out.in_table

    def test_default_leaves_the_class_alone(self):
        out = classify_site(REFERENCE_SITE)
        assert out.site_class is SiteClass.REFERENCE
        assert not out.caller_regulated

    def test_blank_description_rejected(self):
        with pytest.raises(ValueError, match="known_diversions"):
            classify_site(REFERENCE_SITE, known_diversions="  ")

    def test_b17c_refuses_it(self):
        out = classify_site(REFERENCE_SITE, known_diversions=DIVERSIONS)
        with pytest.raises(RegulatedRecordError, match="caller"):
            require_unregulated(out)

    def test_as_dict_records_it(self):
        out = classify_site(REFERENCE_SITE, known_diversions=DIVERSIONS)
        assert out.as_dict()["caller_regulated"] is True


class TestCheckTransferRegulation:
    @pytest.mark.parametrize(
        "donor, target",
        [(REGULATED, REFERENCE), (REFERENCE, REGULATED), (SiteClass.REGULATED, SiteClass.UNKNOWN)],
    )
    def test_mismatch_refused(self, donor, target):
        with pytest.raises(RegulationMismatchError, match="allow_regulation_mismatch"):
            check_transfer_regulation(donor, target)

    def test_mismatch_message_names_both_sides(self):
        with pytest.raises(RegulationMismatchError) as info:
            check_transfer_regulation(REGULATED, REFERENCE)
        assert "11111111" in str(info.value) and "22222222" in str(info.value)

    def test_mismatch_is_a_value_error(self):
        assert issubclass(RegulationMismatchError, ValueError)

    def test_override_proceeds_and_is_recorded(self, caplog):
        with caplog.at_level(logging.WARNING, logger="flowfreq.regulation"):
            prov = check_transfer_regulation(REGULATED, REFERENCE, allow_mismatch=True)
        assert prov["regulation_checked"] is True
        assert prov["regulation_override"] is True
        assert prov["donor"]["site_class"] == "regulated"
        assert prov["target"]["site_class"] == "reference"
        assert "allow_regulation_mismatch" in caplog.text

    @pytest.mark.parametrize(
        "donor, target",
        [
            (SiteClass.REFERENCE, SiteClass.REFERENCE),
            (SiteClass.REFERENCE, SiteClass.URBAN),
            (SiteClass.UNKNOWN, SiteClass.URBAN),
            (SiteClass.URBAN, SiteClass.UNKNOWN),
        ],
    )
    def test_unregulated_pairs_pass_silently(self, donor, target, caplog):
        with caplog.at_level(logging.WARNING, logger="flowfreq.regulation"):
            prov = check_transfer_regulation(donor, target)
        assert prov["regulation_checked"] is True
        assert prov["regulation_override"] is False
        assert caplog.text == ""

    def test_both_regulated_proceeds_with_warning(self, caplog):
        with caplog.at_level(logging.WARNING, logger="flowfreq.regulation"):
            prov = check_transfer_regulation(REGULATED, SiteClass.REGULATED)
        assert prov["regulation_checked"] is True
        assert prov["regulation_override"] is False
        assert "both regulated" in caplog.text

    @pytest.mark.parametrize("donor, target", [(None, REFERENCE), (REFERENCE, None), (None, None)])
    def test_missing_side_warns_and_is_recorded(self, donor, target, caplog):
        with caplog.at_level(logging.WARNING, logger="flowfreq.regulation"):
            prov = check_transfer_regulation(donor, target)
        assert prov["regulation_checked"] is False
        assert prov["regulation_override"] is False
        assert "regulation" in caplog.text

    def test_site_class_value_recorded_without_site_no(self):
        prov = check_transfer_regulation(SiteClass.REFERENCE, SiteClass.URBAN)
        assert prov["donor"] == {"site_class": "reference"}
        assert prov["target"] == {"site_class": "urban"}

    def test_bad_type_rejected(self):
        with pytest.raises(TypeError, match="SiteClassification"):
            check_transfer_regulation("regulated", SiteClass.REFERENCE)  # type: ignore[arg-type]


# --- the four entry points ----------------------------------------------------


@pytest.fixture(scope="module")
def flood_donor():
    from flowfreq.bulletin17c import Bulletin17C
    from tests.fixtures.big_sandy import REGIONAL_SKEW, REGIONAL_SKEW_SD, SYSTEMATIC_PEAKS

    years = np.array(sorted(SYSTEMATIC_PEAKS))
    flows = np.array([SYSTEMATIC_PEAKS[y] for y in years], dtype=float)
    return Bulletin17C(
        peak_flows=flows,
        water_years=years,
        regional_skew=REGIONAL_SKEW,
        regional_skew_mse=REGIONAL_SKEW_SD**2,
    ).run_analysis(method="ema")


@pytest.fixture(scope="module")
def daily() -> pd.DataFrame:
    rng = np.random.default_rng(0)
    index = pd.date_range("2000-10-01", periods=3653, freq="D")
    flows = 200.0 * 10.0 ** rng.normal(0.0, 0.35, size=len(index))
    return pd.DataFrame({"flow_cfs": flows}, index=index)


def _aep_exponents() -> RegressionExponents:
    return RegressionExponents(
        aeps=[0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002],
        exponents=[0.83, 0.81, 0.80, 0.78, 0.77, 0.76, 0.75, 0.74],
        citation="Example SIR 0000-0000 table 6",
    )


def _low_flow() -> LowFlowResults:
    probs = (0.5, 0.2, 0.1, 0.04, 0.02)
    return LowFlowResults(
        n_years=40,
        n_zero_years=0,
        p_zero=0.0,
        n_day=7,
        year_type="climatic",
        distribution="lp3",
        mean_log=1.5,
        std_log=0.3,
        skew_station=-0.1,
        skew_used=-0.1,
        quantiles=pd.DataFrame(
            {
                "non_exceedance_prob": list(probs),
                "return_period": [1.0 / p for p in probs],
                "flow_cfs": [40.0, 30.0, 25.0, 20.0, 18.0],
            }
        ),
    )


def _run(entry: str, flood_donor, daily, **kwargs):
    if entry == "frequency":
        return transpose_frequency(flood_donor, 100.0, 80.0, _aep_exponents(), **kwargs)
    if entry == "duration":
        exps = RegressionExponents(
            aeps=[0.01, 0.10, 0.50, 0.90, 0.99],
            exponents=[0.95, 0.92, 0.87, 0.78, 0.70],
            citation="Example duration regression",
            probability_kind="exceedance",
        )
        return transpose_duration(flow_duration_curve(daily), 100.0, 80.0, exps, **kwargs)
    if entry == "low_flow":
        exps = RegressionExponents(
            aeps=[0.50, 0.20, 0.10, 0.04, 0.02],
            exponents=[0.88, 0.85, 0.83, 0.81, 0.80],
            citation="Example low-flow regression",
            probability_kind="non_exceedance",
        )
        return transpose_low_flow(
            _low_flow(), 100.0, 90.0, exps, hydrogeologic_setting="same aquifer", **kwargs
        )
    curve = FlowDurationCurve.from_daily(daily)
    return qppq(daily, curve, curve, **kwargs)


def _regulation(result) -> dict:
    holder = getattr(result, "provenance", result)
    return holder.regulation


ENTRIES = ["frequency", "duration", "low_flow", "qppq"]


@pytest.mark.parametrize("entry", ENTRIES)
class TestEntryPoints:
    def test_mismatch_refused(self, entry, flood_donor, daily):
        with pytest.raises(RegulationMismatchError):
            _run(entry, flood_donor, daily, donor_regulation=REGULATED, target_regulation=REFERENCE)

    def test_override_recorded(self, entry, flood_donor, daily):
        result = _run(
            entry,
            flood_donor,
            daily,
            donor_regulation=REGULATED,
            target_regulation=REFERENCE,
            allow_regulation_mismatch=True,
        )
        reg = _regulation(result)
        assert reg["regulation_override"] is True
        assert reg["donor"]["site_no"] == "11111111"

    def test_matching_pair_recorded(self, entry, flood_donor, daily):
        result = _run(
            entry,
            flood_donor,
            daily,
            donor_regulation=REFERENCE,
            target_regulation=SiteClass.REFERENCE,
        )
        reg = _regulation(result)
        assert reg["regulation_checked"] is True
        assert reg["regulation_override"] is False

    def test_default_warns_and_leaves_numbers_unchanged(self, entry, flood_donor, daily, caplog):
        checked = _run(
            entry,
            flood_donor,
            daily,
            donor_regulation=REFERENCE,
            target_regulation=SiteClass.REFERENCE,
        )
        with caplog.at_level(logging.WARNING, logger="flowfreq.regulation"):
            unchecked = _run(entry, flood_donor, daily)
        assert "regulation" in caplog.text
        assert _regulation(unchecked)["regulation_checked"] is False
        frame = "series" if entry == "qppq" else "quantiles"
        pd.testing.assert_frame_equal(getattr(checked, frame), getattr(unchecked, frame))


class TestReport:
    def test_markdown_names_the_override(self, flood_donor, daily):
        result = _run(
            "frequency",
            flood_donor,
            daily,
            donor_regulation=REGULATED,
            target_regulation=REFERENCE,
            allow_regulation_mismatch=True,
        )
        assert "regulat" in result.to_markdown().lower()
        assert "override" in result.to_markdown().lower()

    def test_markdown_names_an_unchecked_transfer(self, flood_donor, daily):
        result = _run("frequency", flood_donor, daily)
        assert "not checked" in result.to_markdown()


class TestRegressionWeightedTransfer:
    def test_mismatch_refused(self):
        with pytest.raises(RegulationMismatchError):
            regression_weighted_transfer(
                7000.0,
                6000.0,
                100.0,
                110.0,
                0.8,
                donor_regulation=REGULATED,
                target_regulation=REFERENCE,
            )

    def test_recorded(self):
        out = regression_weighted_transfer(
            7000.0,
            6000.0,
            100.0,
            110.0,
            0.8,
            donor_regulation=REFERENCE,
            target_regulation=REFERENCE,
        )
        assert out.regulation["regulation_checked"] is True

    def test_default_unchanged(self):
        out = regression_weighted_transfer(7000.0, 6000.0, 100.0, 110.0, 0.8)
        assert out.regulation["regulation_checked"] is False
        assert out.area_ratio == pytest.approx(1.1)
