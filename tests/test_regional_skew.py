"""Tests for flowfreq.regional_skew."""

from __future__ import annotations

import pytest

from flowfreq.regional_skew import (
    RegionalSkew,
    RegionalSkewUnavailable,
    available,
    load_table,
    regional_skew_for,
)

HEADER = "state,skew_region,skew,skew_mse,effective_record_years,citation,table,status,notes\n"


def _table(tmp_path, body):
    p = tmp_path / "skew.csv"
    p.write_text(HEADER + body)
    return p


def test_packaged_table_loads_and_wave1_pending_rows_remain():
    # WA/OR/ID were transcribed in Wave 1 (tests/test_regional_skew_wave1.py);
    # MT and Idaho's Snake River Plain have no representable B-GLS skew.
    df = load_table()
    assert set(df.loc[df["status"] == "pending", "state"]) >= {"ID", "MT"}


@pytest.mark.parametrize("state, region", [("MT", None), ("ID", "Snake River Plain")])
def test_pending_rows_are_never_returned(state, region):
    with pytest.raises(RegionalSkewUnavailable, match="pending"):
        regional_skew_for(state, region)


def test_no_entry():
    with pytest.raises(RegionalSkewUnavailable, match="no entry"):
        regional_skew_for("FL")


def test_verified_row(tmp_path):
    p = _table(tmp_path, "XX,R1,-0.1,0.12,50,Report 1,Table 2,verified,\n")
    # XX is not a jurisdiction, but the skew table does not check codes; the catalog does.
    rs = regional_skew_for("xx", path=p)
    assert rs == RegionalSkew("XX", "R1", -0.1, 0.12, "Report 1", "Table 2", 50.0)
    assert available("XX", path=p) == ["R1"]


def test_ambiguous_region(tmp_path):
    p = _table(tmp_path, "XX,R1,-0.1,0.12,,C,,verified,\nXX,R2,0.1,0.1,,C,,verified,\n")
    with pytest.raises(RegionalSkewUnavailable, match="Name one"):
        regional_skew_for("XX", path=p)
    assert regional_skew_for("XX", "R2", path=p).skew == 0.1


@pytest.mark.parametrize(
    "body, match",
    [
        ("XX,R1,,,,,,bogus,\n", "unknown status"),
        ("XX,R1,-0.1,,,C,,verified,\n", "missing skew_mse"),
        ("XX,R1,,,,,,pending,\nXX,R1,,,,,,pending,\n", "duplicate"),
    ],
)
def test_table_validation(tmp_path, body, match):
    with pytest.raises(ValueError, match=match):
        load_table(_table(tmp_path, body))


def test_regional_skew_requires_citation_and_positive_mse():
    with pytest.raises(ValueError, match="citation"):
        RegionalSkew("WA", "R", 0.0, 0.1, " ")
    with pytest.raises(ValueError, match="skew_mse"):
        RegionalSkew("WA", "R", 0.0, 0.0, "C")


def test_a_pending_region_makes_an_unnamed_lookup_ambiguous(tmp_path):
    """One verified and one pending row: still ambiguous without a region name."""
    p = _table(tmp_path, "XX,R1,-0.1,0.12,,C,,verified,\nXX,R2,,,,,,pending,\n")
    with pytest.raises(RegionalSkewUnavailable, match=r"'R2' \(pending\)"):
        regional_skew_for("XX", path=p)
    assert regional_skew_for("XX", "R1", path=p).skew == -0.1
