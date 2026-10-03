"""Tests for flowfreq.future_flow.

Factor values are synthetic except in ``test_reproduces_nchrp_15_61_table_7_1``, which
uses a published worked example built in the test, never shipped as data.
"""

from __future__ import annotations

import json

import pandas as pd
import pytest

from flowfreq.future_flow import (
    ChangeFactorSet,
    apply_change_factors,
    available_factor_sets,
    factor_set_from_dict,
    national_source_review,
    select_factor_set,
)


def _set(**kw):
    base = dict(
        name="synthetic-national",
        aeps=(0.5, 0.1, 0.01),
        factors=(1.05, 1.10, 1.20),
        scenario="test",
        horizon="2050",
        citation="Synthetic test source",
    )
    base.update(kw)
    return ChangeFactorSet(**base)


QUANTILES = pd.DataFrame({"aep": [0.5, 0.1, 0.01, 0.002], "flow_cfs": [100.0, 200.0, 400.0, 600.0]})


def test_apply_keeps_regulatory_and_flags_extrapolation():
    fq = apply_change_factors(QUANTILES, _set())
    t = fq.table
    assert list(t["regulatory_flow_cfs"]) == [100.0, 200.0, 400.0, 600.0]
    assert list(t["change_factor"][:3]) == pytest.approx([1.05, 1.10, 1.20])
    assert t["future_flow_cfs"].iloc[2] == pytest.approx(480.0)
    assert list(t["extrapolated"]) == [False, False, False, True]
    assert t["change_factor"].iloc[3] == pytest.approx(1.20)  # endpoint held
    assert fq.provenance["n_extrapolated"] == 1 and fq.notes


def test_interpolates_between_published_points():
    f, extrap = _set().factor_at([0.05])
    assert 1.10 < f[0] < 1.20 and not extrap[0]


def test_missing_column():
    with pytest.raises(KeyError):
        apply_change_factors(QUANTILES.rename(columns={"flow_cfs": "q"}), _set())


@pytest.mark.parametrize(
    "kw",
    [
        {"citation": ""},
        {"scenario": " "},
        {"factors": (1.0, 1.0)},
        {"factors": (1.0, -1.0, 1.0)},
        {"aeps": (0.5, 0.5, 0.01)},
        {"legal_status": "mandatory"},
        {"level": "state"},  # no geography
    ],
)
def test_invalid_sets(kw):
    with pytest.raises(ValueError):
        _set(**kw)


def test_state_overrides_national_only_where_it_applies():
    nat = _set()
    wa = _set(name="synthetic-wa", level="state", geography=("WA",), legal_status="recommended")
    chosen, overrode = select_factor_set(nat, [wa], state="wa")
    assert chosen is wa and overrode == "synthetic-national"
    chosen, overrode = select_factor_set(nat, [wa], state="OR")
    assert chosen is nat and overrode is None
    fq = apply_change_factors(QUANTILES, chosen, overrode=overrode)
    assert fq.provenance["overrode"] is None


def test_select_rejects_bad_default_and_ambiguity():
    wa1 = _set(name="a", level="state", geography=("WA",))
    wa2 = _set(name="b", level="state", geography=("WA",))
    with pytest.raises(ValueError, match="national"):
        select_factor_set(wa1)
    with pytest.raises(ValueError, match="2 state sets"):
        select_factor_set(_set(), [wa1, wa2], state="WA")


# NCHRP Project 15-61 Design Practices Guide (Kilgore et al. 2019), Table 7.1, p. 80:
# a USGS regression equation evaluated with historical and projected mean annual
# precipitation for one hypothetical watershed.
NCHRP_7_1_AEP = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01)
NCHRP_7_1_HIST = (3290, 5060, 6800, 8920, 10700, 12600)
NCHRP_7_1_PROJ = (4400, 6320, 8220, 10500, 12400, 14300)


def test_reproduces_nchrp_15_61_table_7_1():
    """Per-AEP ratios from the published example reproduce its projected flows."""
    fs = ChangeFactorSet(
        name="nchrp-15-61-table-7.1-example",
        aeps=NCHRP_7_1_AEP,
        factors=tuple(p / h for p, h in zip(NCHRP_7_1_PROJ, NCHRP_7_1_HIST)),
        scenario="hypothetical MAP ratio 1.25",
        horizon="unspecified",
        citation="NCHRP Project 15-61 Design Practices (2019), Table 7.1, p. 80",
    )
    q = pd.DataFrame({"aep": NCHRP_7_1_AEP, "flow_cfs": NCHRP_7_1_HIST})
    t = apply_change_factors(q, fs).table
    assert list(t["future_flow_cfs"]) == pytest.approx(NCHRP_7_1_PROJ)
    assert list(t["regulatory_flow_cfs"]) == list(NCHRP_7_1_HIST)
    # The text (p. 80): "an estimated 16-percent increase in the 0.02 AEP (50-year) discharge".
    assert t.loc[t["aep"] == 0.02, "change_factor"].iloc[0] == pytest.approx(1.16, abs=0.005)


def test_no_national_set_ships_and_review_says_why():
    assert available_factor_sets() == ()
    review = {r.id: r for r in national_source_review()}
    assert {"fhwa-hec17-2016", "nchrp-15-61-2019", "noaa-atlas-15"} <= set(review)
    for r in review.values():
        assert not r.tabulated_national_factors
        assert r.citation and r.url.startswith("https://")
        for t in r.tables_examined:
            assert t["page"] and t["why_not_a_factor_set"]
    assert review["nchrp-15-61-2019"].provides == "method"


def test_available_factor_sets_reads_files(tmp_path):
    nat = {
        "name": "n",
        "aeps": [0.5, 0.01],
        "factors": [1.1, 1.2],
        "scenario": "s",
        "horizon": "2050",
        "citation": "c, Table 1, p. 2",
    }
    st = dict(nat, name="s1", level="state", geography=["WA"])
    (tmp_path / "a.factors.json").write_text(json.dumps(nat))
    (tmp_path / "b.factors.json").write_text(json.dumps([st]))
    (tmp_path / "notes.json").write_text("{}")  # ignored: wrong suffix
    sets = available_factor_sets(tmp_path)
    assert [s.name for s in sets] == ["n", "s1"]
    assert [s.name for s in available_factor_sets(tmp_path, level="state")] == ["s1"]
    assert sets[1].geography == ("WA",)
    (tmp_path / "c.factors.json").write_text(json.dumps(nat))
    with pytest.raises(ValueError, match="duplicate"):
        available_factor_sets(tmp_path)
    with pytest.raises(ValueError):
        available_factor_sets(tmp_path, level="county")


def test_factor_set_from_dict_requires_citation_fields():
    with pytest.raises(KeyError, match="citation"):
        factor_set_from_dict(
            {"name": "x", "aeps": [0.5], "factors": [1.0], "scenario": "s", "horizon": "h"}
        )
