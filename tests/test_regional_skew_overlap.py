"""Overlapping skew studies in the HUC membership table (#204).

Several B-WLS/B-GLS studies cover the same hydrologic units, and each state adopts one:
HU 05 is the Great Lakes and Ohio study in Ohio but SIR 2025-5110's Skew Region 2 in West
Virginia. The membership table therefore allows one row per study for the same HUC, so
long as their ``valid_states`` are disjoint, and ``regional_skew_at_huc`` picks the row
for the site's state. A state none of the matched rows covers keeps the old behavior
(raise), except that a state whose single usable row belongs to a study with no
membership rows at all (Maine's statewide study inside HU 01) falls back to it.

These tests use synthetic membership tables over the packaged value table.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from flowfreq.regional_skew import RegionalSkewUnavailable, load_huc_table, regional_skew_at_huc

SRC = "synthetic test row"


#: Studies whose packaged rows #204 added; the synthetic tables start without them.
ADDED = {
    "Great Lakes and Ohio River Basins",
    "Eastern New York and Pennsylvania",
    "Skew Region 1 (Atlantic drainage)",
    "Skew Region 2 (Gulf of America drainage)",
    "New England",
    "Arkansas-Louisiana",
}


def _table(tmp_path: Path, rows: list) -> Path:
    base = load_huc_table()
    base = base[~base["skew_region"].isin(ADDED)]
    extra = pd.DataFrame(
        [
            {
                "huc": h,
                "skew_region": r,
                "rule": rule,
                "valid_states": v,
                "source": SRC,
                "notes": "",
            }
            for h, r, rule, v in rows
        ]
    )
    path = tmp_path / "hucs.csv"
    pd.concat([base, extra], ignore_index=True).to_csv(path, index=False)
    return path


GLO = "Great Lakes and Ohio River Basins"
SR2 = "Skew Region 2 (Gulf of America drainage)"
NE = "New England"


@pytest.fixture
def ohio(tmp_path: Path) -> Path:
    return _table(tmp_path, [("05", GLO, "member", "OH;PA"), ("05", SR2, "member", "KY;WV")])


def test_duplicate_huc_with_disjoint_states_loads(ohio: Path) -> None:
    hucs = load_huc_table(ohio)
    assert sorted(hucs.loc[hucs["huc"] == "05", "skew_region"]) == [GLO, SR2]


def test_duplicate_huc_with_overlapping_states_is_refused(tmp_path: Path) -> None:
    path = _table(tmp_path, [("05", GLO, "member", "OH;WV"), ("05", SR2, "member", "KY;WV")])
    with pytest.raises(ValueError, match="duplicate HUC.*05"):
        load_huc_table(path)


@pytest.mark.parametrize(
    "huc, state, region",
    [
        ("050200010101", "WV", SR2),
        ("050100010101", "PA", GLO),
        ("050400010101", "OH", GLO),
        ("051002050101", "KY", SR2),
    ],
)
def test_the_site_state_picks_the_study(ohio: Path, huc: str, state: str, region: str) -> None:
    rs = regional_skew_at_huc(huc, state, huc_path=ohio)
    assert rs.skew_region == region and rs.state == state


def test_state_inferred_from_a_single_state_huc8(ohio: Path) -> None:
    # 05040001 (Tuscarawas) lies wholly in Ohio.
    assert regional_skew_at_huc("050400010101", huc_path=ohio).skew_region == GLO


def test_a_border_huc8_needs_state(ohio: Path) -> None:
    # 05030201 (Little Muskingum-Middle Island) touches OH and WV.
    with pytest.raises(RegionalSkewUnavailable, match="Pass state="):
        regional_skew_at_huc("050302010101", huc_path=ohio)


def test_a_state_no_row_covers_still_raises(ohio: Path) -> None:
    # Illinois (pending rows only) inside HU 05: no fallback is possible, so the old
    # "outside the study" refusal stands.
    with pytest.raises(RegionalSkewUnavailable, match="IL is outside"):
        regional_skew_at_huc("051201130101", "IL", huc_path=ohio)


def test_state_bounded_study_falls_back_inside_another_study(tmp_path: Path) -> None:
    # Maine's single verified row is its statewide study, which has no membership rows;
    # inside New England's HU 01 a Maine site gets it, labelled as the state fallback.
    path = _table(tmp_path, [("01", NE, "member", "CT;MA;NH;RI")])
    rs = regional_skew_at_huc("010500020101", "ME", huc_path=path)
    assert rs.skew_region == "Maine (statewide)"
    assert rs.advisory.startswith("STATE-LEVEL FALLBACK")
    assert regional_skew_at_huc("010802010101", "MA", huc_path=path).skew_region == NE


def test_an_exclusion_without_states_applies_to_every_state(tmp_path: Path) -> None:
    path = _table(tmp_path, [("01", NE, "member", "CT;MA;NH;RI"), ("01090002", NE, "excluded", "")])
    with pytest.raises(RegionalSkewUnavailable, match="no regional skew applies"):
        regional_skew_at_huc("010900020101", "MA", huc_path=path)


def test_a_state_specific_exclusion_spares_other_states(tmp_path: Path) -> None:
    path = _table(
        tmp_path, [("01", NE, "member", "CT;MA;NH;RI"), ("01090002", NE, "excluded", "MA")]
    )
    with pytest.raises(RegionalSkewUnavailable, match="no regional skew applies"):
        regional_skew_at_huc("010900020101", "MA", huc_path=path)
    assert regional_skew_at_huc("010900020101", "RI", huc_path=path).skew_region == NE


def test_exclusion_from_a_usable_study_gives_its_own_reason(tmp_path: Path) -> None:
    path = _table(
        tmp_path, [("01", NE, "member", "CT;MA;NH;RI"), ("01090002", NE, "excluded", "MA")]
    )
    with pytest.raises(RegionalSkewUnavailable, match="excluded from the New England study.*"
                       "synthetic test row"):  # fmt: skip
        regional_skew_at_huc("010900020101", "MA", huc_path=path)


def test_multistate_unit_no_row_covers_asks_for_state(tmp_path: Path) -> None:
    # 03050102 touches NC and SC; the only HU 03 row is Virginia's and West Virginia's.
    path = _table(tmp_path, [("03", SR2, "member", "VA;WV")])
    with pytest.raises(RegionalSkewUnavailable, match="touches NC, SC; pass state="):
        regional_skew_at_huc("030501020101", huc_path=path)
    assert regional_skew_at_huc("030501020101", "NC", huc_path=path).advisory.startswith(
        "STATE-LEVEL FALLBACK"
    )


# ----------------------------------------------------------------------------
# Packaged rows (#204), at real gages. HU12s are the Water Data API's
# hydrologic_unit_code for each site (retrieved 2026-10-30).
# ----------------------------------------------------------------------------

ENYPA = "Eastern New York and Pennsylvania"
SR1 = "Skew Region 1 (Atlantic drainage)"
ARLA = "Arkansas-Louisiana"


@pytest.mark.parametrize(
    "site, huc, state, region",
    [
        ("03237500", "050902010506", "OH", GLO),  # Ohio Brush Creek near West Union
        ("03069500", "050200040503", "WV", SR2),  # Cheat River near Parsons
        ("01611500", "020700030708", "WV", SR1),  # Cacapon River near Great Cacapon
        ("01543000", "020502020207", "PA", ENYPA),  # Driftwood Br Sinnemahoning Cr
        ("03015500", "050100011108", "PA", GLO),  # Brokenstraw Creek at Youngsville
        ("01350000", "020200050205", "NY", ENYPA),  # Schoharie Creek at Prattsville
        ("04221000", "041300020503", "NY", GLO),  # Genesee River at Wellsville
        ("01664000", "020801030602", "VA", SR1),  # Rappahannock River at Remington
        ("03170000", "050500011704", "VA", SR2),  # Little River at Graysontown
        ("03488000", "060101010201", "VA", SR2),  # N F Holston River near Saltville
        ("07030500", "080102100302", "TN", SR2),  # Wolf River at Rossville
        ("01193500", "010802050805", "CT", NE),  # Salmon River near East Hampton
        ("01096000", "010700040302", "MA", NE),  # Squannacook River near West Groton
        ("07351500", "111402060107", "LA", ARLA),  # Cypress Bayou near Keithville
        ("07072000", "110100110405", "AR", ARLA),  # Eleven Point River near Ravenden Springs
    ],
)
def test_packaged_rows_pick_the_adopted_study(site: str, huc: str, state: str, region: str) -> None:
    rs = regional_skew_at_huc(huc, state)
    assert rs.skew_region == region, site
    assert not rs.advisory, site


@pytest.mark.parametrize(
    "site, huc, state, region",
    [
        ("01030500", "010200030704", "ME", "Maine (statewide)"),  # Mattawamkeag River
        ("02143000", "030501020000", "NC", "Southeastern United States"),  # Henry Fork
    ],
)
def test_packaged_rows_leave_state_bounded_studies_to_the_fallback(
    site: str, huc: str, state: str, region: str
) -> None:
    rs = regional_skew_at_huc(huc, state)
    assert rs.skew_region == region and rs.advisory.startswith("STATE-LEVEL FALLBACK"), site


@pytest.mark.parametrize(
    "site, huc, state, match",
    [
        ("01374781", "020301010301", "NY", "excluded from the Eastern New York.*Lower Hudson"),
        ("011058837", "010900020203", "MA", "excluded from the New England study.*Cape Cod"),
        ("07381440", "080703000302", "LA", "excluded from the Arkansas-Louisiana study"),
        ("01639500", "020700090405", "MD", "MD is outside"),  # Atlantic MD rows are pending
    ],
)
def test_packaged_rows_refuse_where_no_skew_applies(
    site: str, huc: str, state: str, match: str
) -> None:
    with pytest.raises(RegionalSkewUnavailable, match=match):
        regional_skew_at_huc(huc, state)


def test_packaged_new_rows_need_no_state_inside_one_state() -> None:
    # 04130002 (Upper Genesee) lies wholly in New York, 02080103 (Rapidan-Upper
    # Rappahannock) wholly in Virginia.
    assert regional_skew_at_huc("041300020503").skew_region == GLO
    assert regional_skew_at_huc("020801030602").skew_region == SR1


# ----------------------------------------------------------------------------
# Review fixes (2026-10-31)
# ----------------------------------------------------------------------------


def test_multistate_unit_names_the_study_covering_a_touched_state() -> None:
    # 08020203 touches AR and MO. HU 08 carries Skew Region 2 (KY, TN) and
    # Arkansas-Louisiana (AR, LA) rows; the refusal must not name Skew Region 2.
    with pytest.raises(RegionalSkewUnavailable) as info:
        regional_skew_at_huc("080202030101")
    assert "Skew Region 2" not in str(info.value)
    assert regional_skew_at_huc("080202030101", "AR").skew_region == ARLA
