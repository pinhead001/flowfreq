"""Wave 1 (#34) regional skew rows: the transcribed values, and the gaps.

Every expected number below was read from the primary source cited in the
packaged table, not from memory:

- Wood and others (2016), SIR 2016-5083, appendix B (A.G. Veilleux), Table B2
  p. 53 and text p. 52: Pacific Northwest CONSTANT model, G = -0.07,
  AVPnew (the Bulletin 17 MSE) = 0.18, effective record length 41 years.
- Mastin and others (2016), SIR 2016-5118, appendix A, Table A2 p. 66 and
  text p. 23: the same study, G = -0.07, MSE 0.18, SE 0.4243.
"""

from __future__ import annotations

import math

import pytest

from flowfreq.regional_skew import RegionalSkewUnavailable, available, load_table, regional_skew_for

PNW_SKEW = -0.07
PNW_MSE = 0.18
PNW_ERL = 41.0


@pytest.mark.parametrize("state", ["WA", "OR", "ID"])
def test_pacific_northwest_rows(state):
    rs = regional_skew_for(state, "Pacific Northwest")
    assert rs.skew_region == "Pacific Northwest"
    assert rs.skew == PNW_SKEW
    assert rs.skew_mse == PNW_MSE
    assert rs.effective_record_years == PNW_ERL
    assert "https://doi.org/10.3133/sir2016" in rs.citation
    assert "p. " in rs.table


def test_pnw_mse_matches_published_standard_error():
    # SIR 2016-5118 p. 23 states the SE as 0.4243; MSE = SE**2.
    assert math.sqrt(regional_skew_for("WA").skew_mse) == pytest.approx(0.4243, abs=5e-5)


@pytest.mark.parametrize(
    "state, report", [("WA", "2016-5118"), ("OR", "2016-5083"), ("ID", "2016-5083")]
)
def test_citations_name_the_report(state, report):
    rs = regional_skew_for(state, "Pacific Northwest")
    assert f"Scientific Investigations Report {report}" in rs.citation


def test_idaho_snake_river_plain_is_excluded():
    # SIR 2016-5083 p. 52: the PNW model is not valid in the Snake River Plain.
    assert available("ID") == ["Pacific Northwest"]
    with pytest.raises(RegionalSkewUnavailable, match="pending"):
        regional_skew_for("ID", "Snake River Plain")
    df = load_table()
    row = df[(df["state"] == "ID") & (df["skew_region"] == "Snake River Plain")].iloc[0]
    assert "p. 52" in row["notes"]


def test_montana_stays_pending():
    # SIR 2025-5019 pp. 9-10: no B-GLS study covers MT; practice is the
    # spatially varying B17B Plate I map, which this table never returns.
    assert available("MT") == []
    with pytest.raises(RegionalSkewUnavailable, match="pending"):
        regional_skew_for("MT")
    df = load_table()
    notes = df.loc[df["state"] == "MT", "notes"].iloc[0]
    assert "SIR 2025-5019" in notes and "Plate I" in notes


def test_wave1_rows_all_have_notes():
    df = load_table()
    wave1 = df[df["state"].isin(["WA", "OR", "ID", "MT"])]
    assert wave1["notes"].notna().all()
    verified = wave1[wave1["status"] == "verified"]
    assert set(verified["state"]) == {"WA", "OR", "ID"}
    assert verified["table"].notna().all()


def test_idaho_needs_a_named_region():
    """ID has a verified PNW row and a pending Snake River Plain row: an unnamed
    lookup must not return the PNW value, which that study says is invalid on
    the Plain (SIR 2016-5083 p. 52)."""
    with pytest.raises(RegionalSkewUnavailable, match="Name one"):
        regional_skew_for("ID")
