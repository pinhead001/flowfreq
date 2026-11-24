"""Sets no NSS code reaches, used where NSS does locate a region, and supplied characteristics.

``named_regions`` in a state's file lists the sets NSS does not serve and a caller may put in
place of what NSS located: Texas's OmegaEM set and the 2026 urban sets. Before this they were
refused wherever NSS located one of the state's regions, which is nearly everywhere.
"""

from __future__ import annotations

import json
from dataclasses import replace
from typing import Any, Dict, Set

import pytest

import flowfreq.ungaged as ungaged
from flowfreq import OverlappingRegionsError, estimate_at_site
from flowfreq.regression import evaluate
from flowfreq.regression.library import DATA_DIR, EquationsUnavailable, load_state
from flowfreq.streamstats import RegressionRegionLocation
from flowfreq.ungaged import (
    PartialCoverageError,
    assemble,
    estimate_from_characteristics,
    region_for_nss_code,
)

from .test_ungaged import CHARS, KW, LIB, Fakes, _region, _wa_site, _watershed

NAMED = replace(LIB, named_regions=["R2"])
URBAN = {"DRNAREA": 20.0, "LC11DEV": 60.0, "L3_PIEDMNT": 100.0, "L3_RDGVLY": 0.0, "LC11IMP": 30.0}
TEXAS = {"CONTDA": 737.0, "PRECIPCENT": 24.5, "OmegaEM": 0.33}
STATES = ["TX", "AL", "MS", "TN", "GA", "NC", "SC"]


def _set(state: str) -> str:
    return "Statewide OmegaEM" if state == "TX" else "Urban 2026"


def _other(state: str) -> RegressionRegionLocation:
    """A region NSS serves in the state, as NSS would locate it."""
    lib = load_state(state)
    for number in range(1, 3000):
        code = f"GC{number}"
        if region_for_nss_code(lib, code) not in (None, _set(state)):
            return RegressionRegionLocation(1, code, "located", 100.0, 20.0)
    raise AssertionError(f"no NSS region in {state}")


def _codes(lib: Any, name: str) -> Set[str]:
    return {
        v.code
        for e in lib.equations
        if e.region_code == name
        for v in (*e.variables, *e.applicability)
    }


@pytest.mark.parametrize("state", STATES)
def test_the_state_files_name_their_sets_nss_does_not_serve(state: str) -> None:
    lib = load_state(state)
    assert lib.named_regions == [_set(state)]
    tiers = lib.region_preference.tiers if lib.region_preference else []
    assert all(_set(state) not in tier for tier in tiers)


@pytest.mark.parametrize("state", STATES)
@pytest.mark.parametrize("only_its_own", [True, False])
def test_a_named_set_is_used_where_nss_locates_a_region(state: str, only_its_own: bool) -> None:
    name, lib = _set(state), load_state(state)
    inputs: Dict[str, float] = {**URBAN, **TEXAS, "MCSLOPE": 0.00326}
    if state != "TX":
        del inputs["CONTDA"]
    if only_its_own:
        # Nothing the located set needs but what the named one takes.
        inputs = {k: v for k, v in inputs.items() if k in _codes(lib, name)}
    located = _other(state)
    a = assemble(
        lib,
        [located],
        inputs,
        aeps=[0.01],
        allow_extrapolation=False,
        weighting_space="log",
        region=name.lower(),
    )
    assert a.region_codes[a.used[0].code] == name and len(a.used) == 1
    flow = evaluate(lib.equation(name, 0.01), inputs).flow_cfs
    assert a.quantiles.loc[0, "flow_cfs"] == pytest.approx(flow)
    assert a.quantiles.loc[0, "lower"] < flow < a.quantiles.loc[0, "upper"]
    assert a.alternatives == [located] and a.region_codes[located.code] != name
    assert a.region_choice_source is not None
    assert "NSS does not serve it" in a.region_choice_source
    assert located.code in a.region_choice_source and "whole basin" in a.region_choice_source


def test_a_named_set_takes_the_place_of_overlapping_and_partial_regions() -> None:
    for located in (
        [_region("GC1", 100.0), _region("GC1", 100.0)],
        [_region("GC1", 60.0), _region("GC9", 40.0)],
    ):
        a = assemble(NAMED, located, CHARS, region="R2", **KW)
        assert [r.code for r in a.used] == ["R2"]
        assert [r.code for r in a.unmatched] == [r.code for r in located if r.code == "GC9"]


def test_a_named_set_nss_does_locate_goes_through_the_usual_rules() -> None:
    both = [_region("GC1", 100.0), _region("GC2", 100.0)]
    assert assemble(NAMED, both, CHARS, region="R2", **KW).used == [both[1]]
    part = [_region("GC1", 60.0), _region("GC2", 40.0)]
    assert len(assemble(NAMED, part, CHARS, region="R2", **KW).used) == 2


def test_a_set_not_listed_is_still_refused() -> None:
    with pytest.raises(EquationsUnavailable, match="not in use"):
        assemble(LIB, [_region("GC1", 100.0)], CHARS, region="R2", **KW)


def test_only_the_region_code_reaches_a_named_set() -> None:
    """A word of its name does not: 'Urban' must not mean 'Urban 2026' in Georgia."""
    with pytest.raises((EquationsUnavailable, OverlappingRegionsError)):
        assemble(
            load_state("GA"),
            [_other("GA")],
            URBAN,
            aeps=None,
            allow_extrapolation=False,
            weighting_space="log",
            region="Urban",
        )


def test_one_of_several_codes_cannot_be_a_named_set() -> None:
    with pytest.raises(EquationsUnavailable, match="'R2' is not in use"):
        assemble(NAMED, [_region("GC1", 100.0)], CHARS, region=["GC1", "R2"], **KW)


def test_a_named_set_lacking_a_characteristic_says_which() -> None:
    with pytest.raises(EquationsUnavailable, match="MCSLOPE"):
        assemble(
            load_state("TX"),
            [_other("TX")],
            TEXAS,
            aeps=None,
            allow_extrapolation=False,
            weighting_space="log",
            region="Statewide OmegaEM",
        )


def test_a_named_region_with_no_equations_is_refused(tmp_path) -> None:
    data = json.loads((DATA_DIR / "TX.json").read_text(encoding="utf-8"))
    data["named_regions"] = ["Nowhere"]
    (tmp_path / "TX.json").write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError, match="named_regions for regions with no equations.*Nowhere"):
        load_state("TX", data_dir=tmp_path)


# -- estimate_at_site: characteristics StreamStats does not return -------------------------


@pytest.fixture
def texas(monkeypatch: pytest.MonkeyPatch) -> Fakes:
    site = {
        **_wa_site(),
        "streamstats_region": "TX",
        "characteristics": {"DRNAREA": 760.0, "PRECIPCENT": 24.5, "OmegaEM": 0.33},
    }
    f = Fakes(site)
    f.regions = [_other("TX")]
    monkeypatch.setattr(ungaged, "delineate_and_get_characteristics", f.delineate)
    monkeypatch.setattr(ungaged, "locate_regression_regions", f.locate)
    return f


SUPPLIED = {"MCSLOPE": 0.00326, "CONTDA": 737.0}


def _texas(**kwargs: Any) -> Any:
    return estimate_at_site(29.43, -100.0, state="TX", region="Statewide OmegaEM", **kwargs)


def test_supplied_characteristics_complete_what_streamstats_returns(texas: Fakes) -> None:
    est = _texas(supplied_characteristics=SUPPLIED)
    direct = estimate_from_characteristics(
        "TX", {**TEXAS, "MCSLOPE": 0.00326}, region="Statewide OmegaEM", lat=29.43, lon=-100.0
    )
    assert list(est.quantiles["flow_cfs"]) == pytest.approx(list(direct.quantiles["flow_cfs"]))
    assert est.characteristics["MCSLOPE"] == 0.00326 and est.characteristics["DRNAREA"] == 760.0
    assert est.provenance["supplied_characteristics"] == SUPPLIED
    text = est.to_markdown()
    assert "Supplied by the caller, not StreamStats: CONTDA, MCSLOPE." in text
    assert "NSS does not serve it" in text
    assert est.provenance["alternative_regions"][0]["code"] == texas.regions[0].code


def test_nothing_supplied_is_not_mentioned(texas: Fakes) -> None:
    texas.site["characteristics"].update(SUPPLIED)
    est = _texas()
    assert "Supplied by the caller" not in est.to_markdown()
    assert est.provenance["supplied_characteristics"] == {}


@pytest.mark.parametrize(
    "supplied,message",
    [
        ({"MCSLOP": 0.003}, "MCSLOP.*no TX equation uses"),
        ({"MCSLOPE": float("nan")}, "MCSLOPE.*finite"),
        ({"MCSLOPE": "steep"}, "MCSLOPE.*finite"),
        ({"MCSLOPE": True}, "MCSLOPE.*finite"),
    ],
)
def test_a_bad_supplied_characteristic_is_refused_before_streamstats(
    texas: Fakes, supplied: Dict[str, Any], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        _texas(supplied_characteristics=supplied)
    assert texas.delineate_calls == []


def test_a_value_streamstats_returned_is_not_replaced(texas: Fakes) -> None:
    with pytest.raises(ValueError, match="StreamStats returned OmegaEM"):
        _texas(supplied_characteristics={**SUPPLIED, "OmegaEM": 0.2})


def test_a_withheld_value_may_be_supplied(texas: Fakes, monkeypatch: pytest.MonkeyPatch) -> None:
    def withheld(region: str, lat: float, lon: float, **kw: Any) -> Any:
        ws = _watershed(texas.site, lat, lon)
        kept = {k: c for k, c in ws.characteristics.items() if k != "OmegaEM"}
        return replace(ws, characteristics=kept, unavailable={"OmegaEM": "withheld: inaccurate"})

    monkeypatch.setattr(ungaged, "delineate_and_get_characteristics", withheld)
    est = _texas(supplied_characteristics={**SUPPLIED, "OmegaEM": 0.33})
    assert est.characteristics["OmegaEM"] == 0.33


@pytest.mark.parametrize("state", ["GA", "NC", "SC", "MS"])
def test_a_word_of_a_named_sets_name_is_not_its_code(state: str) -> None:
    """'Urban' took the 2026 set silently: it is the only region with the word in its name."""
    with pytest.raises(EquationsUnavailable, match="'Urban' is not one of its regions"):
        estimate_from_characteristics(state, URBAN, region="Urban", lat=33.7, lon=-84.4)
    with pytest.raises(EquationsUnavailable, match="not a library region"):
        assemble(
            load_state(state),
            [],
            URBAN,
            aeps=None,
            allow_extrapolation=False,
            weighting_space="log",
            region="Urban",
        )


# -- a basin NSS's regions only partly cover -----------------------------------------------


def _located(code: str, percent: float) -> RegressionRegionLocation:
    return RegressionRegionLocation(1, code, "located", percent, 20.0)


def _nss(state: str, region: str) -> str:
    """The NSS code that reaches a library region."""
    lib = load_state(state)
    return next(f"GC{n}" for n in range(1, 3000) if region_for_nss_code(lib, f"GC{n}") == region)


def _covered(state: str, *located: RegressionRegionLocation) -> float:
    a = assemble(
        load_state(state),
        list(located),
        {**URBAN, **TEXAS, "MCSLOPE": 0.00326},
        aeps=[0.01],
        allow_extrapolation=False,
        weighting_space="log",
        region=_set(state),
    )
    return a.covered_percent


def test_a_named_set_covers_what_nss_located_of_the_state() -> None:
    assert _covered("TX", _located(_nss("TX", "Statewide"), 45.0)) == 45.0
    assert _covered("TX", _located(_nss("TX", "Statewide"), 100.0)) == 100.0


def test_overlapping_reports_do_not_add_up_to_the_whole_basin() -> None:
    both = (_located(_nss("TX", "Statewide"), 60.0), _located(_nss("TX", "Brazos 3"), 60.0))
    assert _covered("TX", *both) == 60.0


def test_the_parts_of_one_report_do_add_up() -> None:
    parts = [_located(_nss("AL", r), p) for r, p in (("1", 70.0), ("2", 30.0))]
    assert _covered("AL", *parts) == 100.0


def test_a_partly_covered_basin_is_cautioned_and_may_be_refused(texas: Fakes) -> None:
    texas.regions = [_located(_nss("TX", "Statewide"), 45.0)]
    est = _texas(supplied_characteristics=SUPPLIED)
    assert est.provenance["covered_percent"] == 45.0
    assert "at least 45 percent" in est.advisory and "'Statewide OmegaEM'" in est.advisory
    assert "scaled" not in est.advisory  # nothing was: the set is used for the whole basin
    assert "45 %" in est.provenance["region_choice_source"]
    with pytest.raises(PartialCoverageError) as info:
        _texas(supplied_characteristics=SUPPLIED, min_coverage_percent=90)
    assert info.value.covered_percent == 45.0 and texas.regions[0].code in str(info.value)
    texas.regions = [_located(_nss("TX", "Statewide"), 100.0)]
    whole = _texas(supplied_characteristics=SUPPLIED, min_coverage_percent=90)
    assert "percent of the basin" not in whole.advisory.split("SIR 2009-5087")[0]


def test_the_choice_of_a_region_selector_may_be_supplied(monkeypatch: pytest.MonkeyPatch) -> None:
    lib = load_state("IN")
    assert lib.region_selector is not None
    code = lib.region_selector.variable
    assert ungaged._supplied(lib, {code: 2}) == {code: 2.0}


def test_a_limit_or_applicability_code_may_be_supplied() -> None:
    assert ungaged._supplied(load_state("GA"), {"LC11IMP": 30}) == {"LC11IMP": 30.0}


def test_a_named_set_first_in_a_sequence_is_refused_too() -> None:
    with pytest.raises(EquationsUnavailable, match="not in use"):
        assemble(NAMED, [_region("GC1", 100.0)], CHARS, region=["R2", "GC1"], **KW)


def test_sets_a_report_splits_by_size_are_not_parts_of_the_state() -> None:
    """Colorado small and large share a table and the same ground (32 mi2 divides them)."""
    both = [
        _located(_nss("TX", r), 45.0) for r in ("Colorado small", "Colorado large", "Statewide")
    ]
    assert _covered("TX", *both) == 45.0


def test_alternatives_of_one_report_in_different_tables_do_not_add_up() -> None:
    both = [_located(_nss("TN", r), 60.0) for r in ("HA1", "HA1 DA only")]
    assert _covered("TN", *both) == 60.0


def test_coverage_is_never_over_100_percent() -> None:
    parts = [_located(_nss("AL", r), p) for r, p in (("1", 70.0), ("2", 31.0))]
    assert _covered("AL", *parts) == 100.0


@pytest.mark.parametrize("located", [[], [_located("GC9999", 30.0)]])
def test_a_named_set_where_nss_locates_none_of_the_states_regions(
    texas: Fakes, located: list
) -> None:
    """Nothing says the basin is in Texas: a required coverage is refused, else a caution."""
    texas.regions = located
    est = _texas(supplied_characteristics=SUPPLIED)
    assert est.provenance["covered_percent"] == 0.0
    assert "NSS located none of the state's regression regions" in est.advisory
    assert "at least" not in est.advisory
    with pytest.raises(PartialCoverageError):
        _texas(supplied_characteristics=SUPPLIED, min_coverage_percent=1)


def test_a_region_of_the_states_map_is_not_held_to_that() -> None:
    """Wyoming and Nevada: NSS locates nothing anywhere, so it says nothing of the basin."""
    a = assemble(LIB, [], CHARS, region="R2", **KW)
    assert a.covered_percent == 100.0


def test_a_named_set_over_a_sliver_says_so(texas: Fakes) -> None:
    texas.regions = [_located(_nss("TX", "Statewide"), 0.0)]
    est = _texas(supplied_characteristics=SUPPLIED)
    assert "NSS located none of the state's regression regions" in est.advisory
