"""`flowfreq.regression.border`: the systematic border differences the library ships."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

from flowfreq.regression import (
    BORDER_KM,
    BorderOffset,
    border_offsets,
    border_offsets_between,
    border_offsets_for,
)
from flowfreq.regression.border import OFFSETS_TABLE

REPO_ROOT = Path(__file__).resolve().parent.parent
ROWS = REPO_ROOT / "docs" / "border_qa" / "rows.csv"


def test_the_shipped_table_has_the_systematic_pairs():
    offsets = border_offsets()
    assert len(offsets) == 15  # docs/VERIFICATION_PLAN.md
    assert all(isinstance(o, BorderOffset) for o in offsets)
    assert [(o.home, o.neighbor) for o in offsets] == sorted((o.home, o.neighbor) for o in offsets)
    # Systematic means 3 or more gages and a median beyond 0.1 in log10.
    assert all(o.gages >= 3 and not 0.795 < o.neighbor_over_home < 1.26 for o in offsets)
    assert all(0 <= o.home_nearer <= o.with_fit <= o.gages for o in offsets)


def test_an_offset_carries_the_comparison_with_the_at_site_fits():
    (ri_ct,) = border_offsets_between("RI", "CT")
    assert ri_ct == BorderOffset(
        home="RI",
        neighbor="CT",
        gages=14,
        neighbor_over_home=2.39,
        with_fit=14,
        home_over_at_site=0.87,
        neighbor_over_at_site=1.9,
        home_nearer=10,
    )


def test_between_finds_a_pair_from_either_side_and_ignores_case():
    assert border_offsets_between("ct", "ri") == border_offsets_between("RI", "CT")
    # Idaho and Montana are systematic in both directions: the first state's gages first.
    assert [(o.home, o.neighbor) for o in border_offsets_between("MT", "ID")] == [
        ("MT", "ID"),
        ("ID", "MT"),
    ]


def test_between_is_empty_for_a_pair_that_is_not_systematic_or_not_adjacent():
    assert border_offsets_between("KY", "TN") == ()  # compared, median 0.85
    assert border_offsets_between("WA", "FL") == ()


def test_for_a_state_lists_its_own_gages_first():
    pairs = [(o.home, o.neighbor) for o in border_offsets_for("IA")]
    assert pairs == [("IA", "MN"), ("IA", "WI"), ("SD", "IA")]
    assert border_offsets_for("HI") == ()


@pytest.mark.parametrize("call", [border_offsets_for, lambda c: border_offsets_between("TN", c)])
def test_an_unknown_code_is_an_error(call):
    with pytest.raises(KeyError, match="Unknown jurisdiction code 'ZZ'"):
        call("ZZ")


def test_describe_names_both_states_and_the_at_site_comparison():
    (ar_mo,) = border_offsets_between("AR", "MO")
    assert ar_mo.describe() == (
        "At 11 Arkansas gages within 20 km of Missouri, the Missouri equations gave a median "
        "0.73 times the Arkansas equations at the 1 percent AEP. The Arkansas equations were "
        "nearer the at-site Bulletin 17C fit at 9 of 11 gages."
    )
    assert BORDER_KM == 20.0


def test_describe_says_when_only_some_gages_have_a_fit():
    (ks_mo,) = border_offsets_between("KS", "MO")
    assert (ks_mo.gages, ks_mo.with_fit) == (4, 2)
    assert ks_mo.describe().endswith("fit at 2 of 2 gages with a fit.")


def test_describe_without_any_fit_stops_at_the_ratio(tmp_path):
    table = tmp_path / "offsets.csv"
    table.write_text(
        "home,neighbor,gages,neighbor_over_home,with_fit,home_over_at_site,"
        "neighbor_over_at_site,home_nearer\nTN,KY,3,1.5,0,,,0\n",
        encoding="utf-8",
    )
    (offset,) = border_offsets(table)
    assert offset.home_over_at_site is None and offset.neighbor_over_at_site is None
    assert offset.describe() == (
        "At 3 Tennessee gages within 20 km of Kentucky, the Kentucky equations gave a median "
        "1.5 times the Tennessee equations at the 1 percent AEP."
    )


def test_a_missing_table_is_an_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        border_offsets(tmp_path / "absent.csv")


def test_the_shipped_table_is_what_the_border_check_measured():
    # flowfreq/data/border_offsets.csv is generated (`python tools/border_qa.py offsets`)
    # from docs/border_qa/rows.csv. Regenerating one without the other fails here.
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    import border_qa as bq

    expected = bq.offsets(pd.read_csv(ROWS, dtype={"site_no": str}))
    shipped = pd.read_csv(OFFSETS_TABLE)
    pd.testing.assert_frame_equal(shipped, expected.reset_index(drop=True), check_dtype=False)
