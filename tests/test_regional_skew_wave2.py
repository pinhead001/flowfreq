"""Wave 2 regional skew: Arizona (SIR 2014-5211 app. 1).

The 12-digit HUs are the ``hydrologic_unit_code`` the Water Data API records
for each named gage (queried 2026-10-11).
"""

from __future__ import annotations

import pytest

from flowfreq.regional_skew import (
    RegionalSkewUnavailable,
    load_huc_table,
    load_table,
    regional_skew_at_huc,
    regional_skew_for,
)

AZ_SKEW, AZ_MSE, AZ_ERL = -0.09, 0.08, 85


def test_arizona_row_is_verified():
    row = load_table().set_index(["state", "skew_region"]).loc[("AZ", "Arizona")]
    assert row["status"] == "verified"
    assert (row["skew"], row["skew_mse"], row["effective_record_years"]) == (
        AZ_SKEW,
        AZ_MSE,
        AZ_ERL,
    )
    assert "SIR 2014-5211" in row["citation"] or "2014-5211" in row["citation"]
    assert "Table 1-2 (p. 55)" in row["table"]


def test_regional_skew_for_arizona():
    rs = regional_skew_for("AZ")
    assert (rs.skew, rs.skew_mse) == (AZ_SKEW, AZ_MSE)
    assert not rs.provisional and rs.advisory == ""


def test_arizona_members_are_huc8_valid_in_az_only():
    hucs = load_huc_table()
    az = hucs[hucs["skew_region"] == "Arizona"]
    assert len(az) == 84
    assert (az["rule"] == "member").all() and (az["valid_states"] == "AZ").all()
    assert az["huc"].str.len().eq(8).all()
    assert az["huc"].str[:2].isin(["14", "15", "16", "18"]).all()


@pytest.mark.parametrize(
    "site, huc",
    [
        ("09498500", "150601030905"),  # Salt River near Roosevelt, AZ
        ("09505800", "150602030108"),  # West Clear Creek near Camp Verde, AZ
    ],
)
def test_arizona_only_huc_infers_the_state(site, huc):
    rs = regional_skew_at_huc(huc)
    assert rs.state == "AZ" and rs.skew_region == "Arizona"
    assert (rs.skew, rs.skew_mse, rs.effective_record_years) == (AZ_SKEW, AZ_MSE, AZ_ERL)


@pytest.mark.parametrize(
    "huc, touches",
    [
        ("150502020706", "AZ, MX"),  # 09471000 San Pedro River at Charleston, AZ
        ("150400040905", "AZ, NM"),  # 09444500 San Francisco River at Clifton, AZ
    ],
)
def test_shared_huc8_needs_the_state(huc, touches):
    with pytest.raises(RegionalSkewUnavailable, match=f"touches {touches}.*Pass state="):
        regional_skew_at_huc(huc)
    assert regional_skew_at_huc(huc, "AZ").skew == AZ_SKEW


def test_new_mexico_side_of_a_shared_huc8_is_outside():
    # 09442680 San Francisco River near Reserve, NM: same HUC8 as Clifton.
    with pytest.raises(RegionalSkewUnavailable, match="NM is outside the Arizona study"):
        regional_skew_at_huc("150400040311", "NM")


def test_new_mexico_only_huc8_has_no_study():
    # 09430500 Gila River near Gila, NM.
    with pytest.raises(RegionalSkewUnavailable, match="No verified regional skew study"):
        regional_skew_at_huc("150400010907", "NM")


def test_new_mexico_is_pending_until_tillery_2026_is_published():
    row = load_table().set_index(["state", "skew_region"]).loc[("NM", "New Mexico")]
    assert row["status"] == "pending"
    assert "SIR 2008-5119, p. 5" in row["notes"] and "10.5066/P9QWEXCY" in row["notes"]
    with pytest.raises(
        RegionalSkewUnavailable,
        match=r"NM \(pending: No verified regional skew for New Mexico; see",
    ):
        regional_skew_for("NM")


@pytest.mark.parametrize(
    "state, report",
    [
        ("CO", "SIR 2025-5047"),
        ("UT", "SIR 2007-5158"),
        ("WY", "SIR 2025-5019"),
        ("NV", "WSP 2433"),
    ],
)
def test_wave2_states_without_a_bayesian_study_are_pending(state, report):
    row = load_table().set_index("state").loc[state]
    assert row["status"] == "pending" and report in row["notes"]
    with pytest.raises(
        RegionalSkewUnavailable,
        match=rf"{state} \(pending: No B-WLS/B-GLS regional skew study covers",
    ):
        regional_skew_for(state)
