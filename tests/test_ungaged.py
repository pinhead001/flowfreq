"""Tests for ``flowfreq.ungaged``: ``SiteEstimate`` and the pure region assembly.

The library is synthetic: ``GC1`` matches a located NSS region by ``region_code``,
``R2`` matches ``GC2`` through its ``region_name`` (as the verified files write
"... (NSS GC1239)"), and ``R4`` needs a characteristic (``PRECIP``) the basin may
not carry.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import pytest

from flowfreq.bulletin17c import Bulletin17C
from flowfreq.regression.equations import (
    Citation,
    OutOfRangeError,
    RegressionEquation,
    RegressionEstimate,
    Variable,
)
from flowfreq.regression.library import (
    NSS_SNAPSHOT_SUBDIR,
    EquationsUnavailable,
    RegionGate,
    StateLibrary,
)
from flowfreq.streamstats import (
    Characteristic,
    Provenance,
    RegressionRegionLocation,
    SnapResult,
    UnsnappablePointError,
    WatershedCharacteristics,
)
from flowfreq.ungaged import (
    Assembly,
    OverlappingRegionsError,
    SiteEstimate,
    assemble,
    estimate_at_site,
    outlet_coordinates,
    region_for_nss_code,
)

CIT1 = Citation(publication="Synthetic report one", table="T1", year=2020)
CIT2 = Citation(publication="Synthetic report two", table="T2", year=2021)
AEPS = (0.5, 0.1, 0.01)
AREA = Variable("DRNAREA", "log10", "mi2", 1.0, 1000.0)
PRECIP = Variable("PRECIP", "identity", "in", 10.0, 60.0)
KW = dict(aeps=None, allow_extrapolation=False, weighting_space="log")


def _eq(region_code: str, aep: float, intercept: float, **kw) -> RegressionEquation:
    base = dict(
        region_code=region_code,
        aep=aep,
        intercept=intercept,
        variables=(AREA,),
        coefficients=(0.7,),
        citation=CIT1,
        sep_log=0.2,
        n_sites=30,
    )
    base.update(kw)
    return RegressionEquation(**base)


def _library() -> StateLibrary:
    eqs = []
    for i, aep in enumerate(AEPS):
        eqs.append(_eq("GC1", aep, 2.0 + 0.3 * i, region_name="Region 1 (NSS GC1)"))
        eqs.append(
            _eq(
                "R2",
                aep,
                1.5 + 0.3 * i,
                region_name="Region 2, Plains (WRIR 03-4107; NSS GC2)",
                citation=CIT2,
            )
        )
        eqs.append(
            _eq(
                "R4",
                aep,
                1.0 + 0.3 * i,
                variables=(AREA, PRECIP),
                coefficients=(0.7, 0.01),
                region_name="Region 4 (NSS GC4)",
            )
        )
    # GC1 alone publishes a 0.002 AEP: it must drop out of the default intersection.
    eqs.append(_eq("GC1", 0.002, 3.2, region_name="Region 1 (NSS GC1)"))
    return StateLibrary(state="XX", status="verified", equations=eqs, advisory="")


def _region(code: str, percent: float) -> RegressionRegionLocation:
    return RegressionRegionLocation(
        region_id=int(code[-1]),
        code=code,
        name=f"Located {code}",
        percent_weight=percent,
        area_sq_mi=percent,
    )


LIB = _library()
CHARS = {"DRNAREA": 100.0}


def _flow(region_code: str, aep: float, chars=CHARS) -> float:
    eq = LIB.equation(region_code, aep)
    assert eq.intercept is not None
    return 10.0 ** (eq.intercept + 0.7 * math.log10(chars["DRNAREA"]))


class TestRegionForNssCode:
    def test_equality_with_region_code_wins(self):
        assert region_for_nss_code(LIB, "GC1") == "GC1"

    def test_falls_back_to_a_region_name_containing_the_code(self):
        assert region_for_nss_code(LIB, "GC2") == "R2"

    def test_no_match_is_none(self):
        assert region_for_nss_code(LIB, "GC9") is None

    def test_the_code_must_be_a_whole_token_of_the_name(self):
        # "GC1618" must not be found by "GC161", nor "GC1" by "GC1618".
        lib = StateLibrary(
            state="AZ",
            status="verified",
            equations=[
                _eq("1", 0.01, 2.0, region_name="Flood region 1 (SIR 2014-5211; NSS GC1618)")
            ],
        )
        assert region_for_nss_code(lib, "GC161") is None
        assert region_for_nss_code(lib, "GC1") is None
        assert region_for_nss_code(lib, "GC1618") == "1"


class TestAssemble:
    def test_single_region_keeps_its_interval(self):
        a = assemble(LIB, [_region("GC1", 100.0)], CHARS, **KW)
        assert isinstance(a, Assembly)
        assert a.unmatched == [] and a.zero_weight == [] and a.missing == {}
        assert [r.code for r in a.used] == ["GC1"]
        assert a.region_codes == {"GC1": "GC1"}
        row = a.quantiles.set_index("aep").loc[0.01]
        ((fraction, est),) = a.estimates[0.01]
        assert fraction == 1.0
        assert isinstance(est, RegressionEstimate)
        assert row["flow_cfs"] == pytest.approx(_flow("GC1", 0.01))
        assert est.interval is not None
        assert (row["lower"], row["upper"]) == pytest.approx(est.interval)
        assert row["interval_level"] == est.interval_level == 0.90
        assert row["interval_method"] == est.interval_method
        assert bool(row["extrapolated"]) is False
        assert row["regions"] == "GC1"

    def test_two_regions_area_weighted_in_log_space(self):
        regions = [_region("GC1", 60.0), _region("GC2", 40.0)]
        a = assemble(LIB, regions, CHARS, **KW)
        assert a.unmatched == [] and a.missing == {}
        assert a.region_codes == {"GC1": "GC1", "GC2": "R2"}
        row = a.quantiles.set_index("aep").loc[0.1]
        q1, q2 = _flow("GC1", 0.1), _flow("R2", 0.1)
        expected = 10 ** (0.6 * math.log10(q1) + 0.4 * math.log10(q2))
        assert row["flow_cfs"] == pytest.approx(expected)
        assert pd.isna(row["lower"]) and pd.isna(row["upper"])
        assert pd.isna(row["interval_level"])
        assert row["interval_method"] == "none: area-weighted across 2 regions"
        assert row["regions"] == "GC1, GC2"
        assert [(f, e.equation.region_code) for f, e in a.estimates[0.1]] == [
            (0.6, "GC1"),
            (0.4, "R2"),
        ]

    def test_linear_weighting_space(self):
        regions = [_region("GC1", 60.0), _region("GC2", 40.0)]
        a = assemble(LIB, regions, CHARS, **{**KW, "aeps": [0.1], "weighting_space": "linear"})
        expected = 0.6 * _flow("GC1", 0.1) + 0.4 * _flow("R2", 0.1)
        assert a.quantiles.loc[0, "flow_cfs"] == pytest.approx(expected)

    def test_fractions_renormalise_after_an_unmatched_region(self):
        regions = [_region("GC1", 60.0), _region("GC2", 30.0), _region("GC3", 10.0)]
        a = assemble(LIB, regions, CHARS, **{**KW, "aeps": [0.01]})
        assert a.unmatched == [regions[2]]
        assert a.missing == {}
        assert a.used == regions[:2]
        fractions = [f for f, _ in a.estimates[0.01]]
        assert fractions == pytest.approx([2 / 3, 1 / 3])
        assert a.quantiles.loc[0, "regions"] == "GC1, GC2"

    def test_no_matching_region_raises_naming_the_codes(self):
        regions = [_region("GC9", 70.0), _region("GC8", 30.0)]
        with pytest.raises(EquationsUnavailable, match="GC9"):
            assemble(LIB, regions, CHARS, **KW)
        with pytest.raises(EquationsUnavailable, match="GC8"):
            assemble(LIB, regions, CHARS, **KW)

    def test_no_located_region_says_nss_located_none(self):
        with pytest.raises(
            EquationsUnavailable,
            match=r"XX: NSS located no regression region for this basin.*check the state "
            r"and that the point is on a stream NSS covers",
        ):
            assemble(LIB, [], CHARS, **KW)

    def test_zero_weight_region_is_dropped_before_evaluation(self):
        # NSS rounds percent_weight to a whole percent, so a sliver arrives as 0. GC4
        # would need PRECIP (missing) and, with DRNAREA out of range, would raise
        # OutOfRangeError if it were evaluated; dropped first, it does neither, and the
        # basin keeps GC1's own interval instead of taking the multi-region path.
        regions = [_region("GC1", 100.0), _region("GC4", 0.0)]
        a = assemble(LIB, regions, CHARS, **KW)
        assert a.zero_weight == [regions[1]]
        assert a.missing == {} and a.unmatched == []
        assert a.used == [regions[0]]
        assert a.region_codes == {"GC1": "GC1", "GC4": "R4"}
        row = a.quantiles.set_index("aep").loc[0.01]
        assert [(f, e.equation.region_code) for f, e in a.estimates[0.01]] == [(1.0, "GC1")]
        assert row["interval_method"] == a.estimates[0.01][0][1].interval_method
        assert row["lower"] == pytest.approx(a.estimates[0.01][0][1].interval[0])
        assert row["regions"] == "GC1"
        big = {"DRNAREA": 5000.0, "PRECIP": 30.0}
        a = assemble(LIB, regions, big, **{**KW, "aeps": [0.01], "allow_extrapolation": True})
        assert [e.equation.region_code for _, e in a.estimates[0.01]] == ["GC1"]

    def test_every_region_zero_weight_raises(self):
        with pytest.raises(EquationsUnavailable, match="GC1.*zero weight"):
            assemble(LIB, [_region("GC1", 0.0)], CHARS, **KW)

    def test_missing_characteristic_skips_the_region_and_records_it(self):
        regions = [_region("GC1", 50.0), _region("GC4", 50.0)]
        a = assemble(LIB, regions, CHARS, **KW)
        assert a.unmatched == []
        assert a.missing == {"GC4": ["PRECIP"]}
        assert a.used == [regions[0]]
        # GC1 carries the whole basin once GC4 is skipped, interval included.
        assert [(f, e.equation.region_code) for f, e in a.estimates[0.5]] == [(1.0, "GC1")]
        assert a.quantiles.loc[0, "regions"] == "GC1"
        assert a.quantiles.loc[0, "interval_method"] == a.estimates[0.5][0][1].interval_method
        assert a.quantiles.loc[0, "lower"] == pytest.approx(a.estimates[0.5][0][1].interval[0])

    def test_skipped_region_does_not_narrow_the_default_aeps(self):
        # GC4 publishes no 0.002; once it is skipped the default AEPs are GC1's alone.
        a = assemble(LIB, [_region("GC1", 50.0), _region("GC4", 50.0)], CHARS, **KW)
        assert list(a.quantiles["aep"]) == [0.5, 0.1, 0.01, 0.002]

    def test_skip_happens_before_an_unpublished_requested_aep_can_raise(self):
        # GC4 neither publishes 0.002 nor has PRECIP; it is skipped, not a failure.
        regions = [_region("GC1", 50.0), _region("GC4", 50.0)]
        a = assemble(LIB, regions, CHARS, **{**KW, "aeps": [0.002]})
        assert a.missing == {"GC4": ["PRECIP"]}
        assert list(a.quantiles["aep"]) == [0.002]
        assert a.quantiles.loc[0, "flow_cfs"] == pytest.approx(_flow("GC1", 0.002))

    def test_missing_characteristic_checks_only_the_requested_aeps(self):
        # R4 needs PRECIP at 0.01 alone: requesting the others leaves it evaluable, and
        # it is skipped only when the request reaches the AEP whose equation needs it.
        eqs = [e for e in LIB.equations if e.region_code != "R4"]
        for i, aep in enumerate(AEPS):
            kw: Dict[str, Any] = dict(region_name="Region 4 (NSS GC4)")
            if aep == 0.01:
                kw.update(variables=(AREA, PRECIP), coefficients=(0.7, 0.01))
            eqs.append(_eq("R4", aep, 1.0 + 0.3 * i, **kw))
        lib = StateLibrary(state="XX", status="verified", equations=eqs, advisory="")
        regions = [_region("GC1", 50.0), _region("GC4", 50.0)]
        a = assemble(lib, regions, CHARS, **{**KW, "aeps": [0.5, 0.1]})
        assert a.missing == {}
        assert a.used == regions
        assert [e.equation.region_code for _, e in a.estimates[0.1]] == ["GC1", "R4"]
        a = assemble(lib, regions, CHARS, **{**KW, "aeps": [0.5, 0.01]})
        assert a.missing == {"GC4": ["PRECIP"]}
        assert a.used == [regions[0]]
        # No request: every published AEP is checked, as before.
        a = assemble(lib, regions, CHARS, **KW)
        assert a.missing == {"GC4": ["PRECIP"]}

    def test_every_region_missing_a_characteristic_raises(self):
        with pytest.raises(EquationsUnavailable, match="GC4 needs PRECIP"):
            assemble(LIB, [_region("GC4", 100.0)], CHARS, **KW)

    def test_extrapolated_flag_follows_out_of_range(self):
        big = {"DRNAREA": 5000.0}
        with pytest.raises(OutOfRangeError):
            assemble(LIB, [_region("GC1", 100.0)], big, **{**KW, "aeps": [0.01]})
        regions = [_region("GC1", 60.0), _region("GC2", 40.0)]
        a = assemble(
            LIB, regions, big, aeps=[0.01], allow_extrapolation=True, weighting_space="log"
        )
        assert bool(a.quantiles.loc[0, "extrapolated"]) is True
        assert all(e.extrapolated for _, e in a.estimates[0.01])
        in_range = assemble(
            LIB,
            [_region("GC1", 100.0)],
            CHARS,
            aeps=[0.01],
            allow_extrapolation=True,
            weighting_space="log",
        )
        assert bool(in_range.quantiles.loc[0, "extrapolated"]) is False

    def test_default_aeps_are_the_intersection_largest_first(self):
        regions = [_region("GC1", 60.0), _region("GC2", 40.0)]
        a = assemble(LIB, regions, CHARS, **KW)
        assert list(a.quantiles["aep"]) == [0.5, 0.1, 0.01]
        assert list(a.estimates) == [0.5, 0.1, 0.01]
        alone = assemble(LIB, [_region("GC1", 100.0)], CHARS, **KW)
        assert list(alone.quantiles["aep"]) == [0.5, 0.1, 0.01, 0.002]

    def test_disjoint_default_aeps_raise_naming_the_regions(self):
        lib = StateLibrary(
            state="XX",
            status="verified",
            equations=[_eq("GC1", 0.5, 2.0), _eq("GC2", 0.01, 2.0)],
        )
        regions = [_region("GC1", 60.0), _region("GC2", 40.0)]
        with pytest.raises(
            EquationsUnavailable, match=r"no AEP in common: GC1 at 0.5; GC2 at 0.01"
        ):
            assemble(lib, regions, CHARS, **KW)

    def test_empty_aeps_are_refused(self):
        with pytest.raises(ValueError, match="aeps must be non-empty"):
            assemble(LIB, [_region("GC1", 100.0)], CHARS, **{**KW, "aeps": []})

    def test_requested_aeps_keep_the_caller_order(self):
        a = assemble(LIB, [_region("GC1", 100.0)], CHARS, **{**KW, "aeps": [0.01, 0.5]})
        assert list(a.quantiles["aep"]) == [0.01, 0.5]

    def test_unpublished_requested_aep_raises(self):
        regions = [_region("GC1", 60.0), _region("GC2", 40.0)]
        with pytest.raises(EquationsUnavailable, match="0.002"):
            assemble(LIB, regions, CHARS, **{**KW, "aeps": [0.002]})

    def test_frame_columns(self):
        a = assemble(LIB, [_region("GC1", 100.0)], CHARS, **KW)
        assert list(a.quantiles.columns) == [
            "aep",
            "flow_cfs",
            "lower",
            "upper",
            "interval_level",
            "interval_method",
            "extrapolated",
            "regions",
        ]


def _estimate(regions, assembly: Optional[Assembly] = None, **kw) -> SiteEstimate:
    a = assembly if assembly is not None else assemble(LIB, regions, CHARS, **KW)
    base = dict(
        lat=45.0,
        lon=-110.0,
        state="XX",
        method="regression",
        quantiles=a.quantiles,
        regions=a.used,
        unmatched_regions=a.unmatched,
        characteristics=dict(CHARS),
        library_status=LIB.status,
        advisory="Synthetic advisory text.",
        estimates=a.estimates,
        region_codes=a.region_codes,
        provenance={
            "missing_characteristics": a.missing,
            "zero_weight_regions": [r.code for r in a.zero_weight],
        },
    )
    base.update(kw)
    return SiteEstimate(**base)


class TestSiteEstimate:
    def test_to_markdown_names_regions_fractions_and_citations(self):
        regions = [_region("GC1", 60.0), _region("GC2", 30.0), _region("GC3", 10.0)]
        text = _estimate(regions).to_markdown()
        assert "45.0" in text and "-110.0" in text and "XX" in text
        assert "| GC1 | Located GC1 | 0.667 | Synthetic report one (2020), table T1 |" in text
        assert "| GC2 | Located GC2 | 0.333 | Synthetic report two (2021), table T2 |" in text
        assert "GC3" in text  # the unmatched region is reported, not silently dropped
        assert "Synthetic advisory text." in text
        assert "DRNAREA" in text and "100" in text
        assert "| AEP (%) | Flow (cfs) |" in text
        assert "| 1 |" in text
        assert "none: area-weighted across 2 regions" in text

    def test_to_markdown_two_nss_codes_on_one_library_region(self):
        # GC2 reaches R2 through its name and R2 by equality: two rows, one citation.
        regions = [_region("GC2", 75.0), _region("R2", 25.0)]
        text = _estimate(regions).to_markdown()
        assert "| GC2 | Located GC2 | 0.750 | Synthetic report two (2021), table T2 |" in text
        assert "| R2 | Located R2 | 0.250 | Synthetic report two (2021), table T2 |" in text

    def test_to_markdown_single_region_shows_interval(self):
        est = _estimate([_region("GC1", 100.0)])
        text = est.to_markdown()
        lower = est.quantiles.set_index("aep").loc[0.01, "lower"]
        assert f"{lower:,.0f}" in text
        assert "| GC1 | Located GC1 | 1.000 |" in text

    def test_to_markdown_caveats_name_skipped_and_dropped_regions(self):
        regions = [_region("GC1", 50.0), _region("GC4", 50.0), _region("GC2", 0.0)]
        text = _estimate(regions).to_markdown()
        assert "GC4 needs PRECIP" in text
        assert "dropped: GC2." in text
        assert "| GC1 | Located GC1 | 1.000 |" in text
        assert "| GC4 |" not in text and "| GC2 |" not in text

    def test_to_markdown_renders_a_weighted_frame(self):
        weighted = pd.DataFrame(
            {"aep": [0.01], "at_site": [12345.6], "weighted": [11000.0], "extrapolated": [False]}
        )
        text = _estimate([_region("GC1", 100.0)], weighted=weighted).to_markdown()
        assert "## At-site weighted with regression" in text
        assert "| aep | at_site | weighted | extrapolated |" in text
        assert "| 0.01 | 12,345.6 | 11,000 | False |" in text

    def test_weighted_and_provenance_default_empty(self):
        est = SiteEstimate(
            lat=0.0,
            lon=0.0,
            state="XX",
            method="regression",
            quantiles=pd.DataFrame(),
            regions=[],
            unmatched_regions=[],
            characteristics={},
            library_status="verified",
            advisory="",
            estimates={},
        )
        assert est.weighted is None and est.provenance == {} and est.region_codes == {}
        assert est.basin_geojson is None
        assert "XX" in est.to_markdown()


class TestOverlappingRegions:
    """NSS serves superseded and current equation sets for one area side by side, each
    at 100 %: they are alternatives, and only a caller can choose between them."""

    def test_overlapping_regions_raise_naming_candidates(self):
        regions = [_region("GC1", 100.0), _region("GC2", 100.0)]
        with pytest.raises(OverlappingRegionsError) as info:
            assemble(LIB, regions, CHARS, **KW)
        exc = info.value
        assert isinstance(exc, EquationsUnavailable)
        assert exc.candidates == [
            ("GC1", "Region 1 (NSS GC1)", 2020),
            ("GC2", "Region 2, Plains (WRIR 03-4107; NSS GC2)", 2021),
        ]
        text = str(exc)
        assert "GC1" in text and "GC2" in text and "2020" in text and "2021" in text
        assert "region=" in text and "--region" in text

    def test_overlap_with_region_chosen_uses_only_that_region(self):
        regions = [_region("GC1", 100.0), _region("GC2", 100.0)]
        a = assemble(LIB, regions, CHARS, region="GC2", **KW)
        assert [r.code for r in a.used] == ["GC2"]
        assert a.alternatives == [regions[0]]
        ((fraction, est),) = a.estimates[0.01]
        assert fraction == 1.0
        row = a.quantiles.set_index("aep").loc[0.01]
        assert row["flow_cfs"] == pytest.approx(_flow("R2", 0.01))
        assert est.interval is not None
        assert (row["lower"], row["upper"]) == pytest.approx(est.interval)
        assert row["regions"] == "GC2"

    def test_region_may_be_given_as_library_code(self):
        regions = [_region("GC1", 100.0), _region("GC2", 100.0)]
        a = assemble(LIB, regions, CHARS, region="R2", **KW)
        assert [r.code for r in a.used] == ["GC2"]
        assert a.alternatives == [regions[0]]

    def test_partition_is_unchanged_by_the_overlap_check(self):
        regions = [_region("GC1", 60.0), _region("GC2", 40.0)]
        a = assemble(LIB, regions, CHARS, **KW)
        assert a.alternatives == []
        assert [f for f, _ in a.estimates[0.1]] == pytest.approx([0.6, 0.4])

    def test_overlap_tolerance_matches_nss_rounding(self):
        a = assemble(LIB, [_region("GC1", 51.0), _region("GC2", 50.0)], CHARS, **KW)
        assert len(a.used) == 2 and a.alternatives == []
        with pytest.raises(OverlappingRegionsError):
            assemble(LIB, [_region("GC1", 100.0), _region("GC2", 2.0)], CHARS, **KW)

    def test_four_region_partition_within_rounding_is_averaged(self):
        # 25.5/25.5/25.5/23.5 rounds to 26/26/26/24 = 102: a partition, not an overlap.
        regions = [
            _region("GC1", 26.0),
            _region("GC2", 26.0),
            _region("GC4", 26.0),
            _region("R2", 24.0),
        ]
        chars = {"DRNAREA": 100.0, "PRECIP": 30.0}
        a = assemble(LIB, regions, chars, **{**KW, "aeps": [0.1]})
        assert a.alternatives == [] and a.used == regions
        assert [f for f, _ in a.estimates[0.1]] == pytest.approx(
            [26 / 102, 26 / 102, 26 / 102, 24 / 102]
        )
        assert a.quantiles.loc[0, "interval_method"] == "none: area-weighted across 4 regions"

    def test_region_choice_is_case_insensitive(self):
        regions = [_region("GC1", 100.0), _region("GC2", 100.0)]
        for choice in ("gc2", "r2"):
            a = assemble(LIB, regions, CHARS, region=choice, **KW)
            assert [r.code for r in a.used] == ["GC2"]

    def test_lower_case_choice_still_gets_the_specific_reason(self):
        regions = [_region("GC1", 100.0), _region("GC4", 100.0)]
        with pytest.raises(EquationsUnavailable, match="GC4") as info:
            assemble(LIB, regions, CHARS, region="gc4", **KW)
        assert "PRECIP" in str(info.value) and "skipped" in str(info.value)

    def test_library_code_shared_by_two_nss_codes_is_ambiguous(self):
        # GC3 and GC5 both reach library region R3 through its name: "R3" cannot pick one.
        lib = StateLibrary(
            state="XX",
            status="verified",
            equations=[_eq("R3", 0.01, 2.0, region_name="Region 3 (NSS GC3; NSS GC5)")],
        )
        regions = [_region("GC3", 100.0), _region("GC5", 100.0)]
        with pytest.raises(OverlappingRegionsError, match="ambiguous"):
            assemble(lib, regions, CHARS, region="R3", **KW)
        a = assemble(lib, regions, CHARS, region="GC5", **KW)
        assert [r.code for r in a.used] == ["GC5"]

    def test_unknown_region_choice_raises_naming_candidates(self):
        regions = [_region("GC1", 100.0), _region("GC2", 100.0)]
        with pytest.raises(OverlappingRegionsError, match="GC9") as info:
            assemble(LIB, regions, CHARS, region="GC9", **KW)
        assert [c[0] for c in info.value.candidates] == ["GC1", "GC2"]
        assert "GC1" in str(info.value) and "GC2" in str(info.value)

    def test_skipped_region_is_not_a_candidate(self):
        # GC4 lacks PRECIP and is skipped first; GC1 alone remains, so no overlap.
        regions = [_region("GC1", 100.0), _region("GC4", 100.0)]
        a = assemble(LIB, regions, CHARS, **KW)
        assert [r.code for r in a.used] == ["GC1"] and a.alternatives == []

    # Fix round 1, finding 1: a chosen set that covers only part of the basin.
    def test_chosen_region_covering_part_of_the_basin_is_refused(self):
        # A boundary basin: GC1 60 / GC2 40 on one side of the line, GC4 100 a set
        # overlapping both. GC1 alone would drop 40 % of the basin.
        chars = {"DRNAREA": 100.0, "PRECIP": 30.0}
        regions = [_region("GC1", 60.0), _region("GC2", 40.0), _region("GC4", 100.0)]
        with pytest.raises(OverlappingRegionsError, match="only 60 % of the basin") as info:
            assemble(LIB, regions, chars, region="GC1", **KW)
        assert isinstance(info.value, LookupError)
        assert "boundary" in str(info.value)
        assert [c[0] for c in info.value.candidates] == ["GC1", "GC2", "GC4"]

    def test_chosen_region_within_rounding_of_full_coverage_is_accepted(self):
        regions = [_region("GC1", 99.0), _region("GC2", 100.0)]
        a = assemble(LIB, regions, CHARS, region="GC1", **KW)
        assert [r.code for r in a.used] == ["GC1"]

    # Fix round 1, finding 2: a choice that names no region in use always raises.
    def test_choice_skipped_for_a_characteristic_raises_without_overlap(self):
        # GC4 lacks PRECIP; GC1 alone remains, so there is no overlap, but the
        # estimate would rest on a set the caller did not choose.
        regions = [_region("GC1", 100.0), _region("GC4", 100.0)]
        with pytest.raises(EquationsUnavailable, match="GC4") as info:
            assemble(LIB, regions, CHARS, region="GC4", **KW)
        assert "PRECIP" in str(info.value) and "skipped" in str(info.value)

    def test_choice_skipped_for_a_characteristic_raises_with_overlap(self):
        regions = [_region("GC1", 100.0), _region("GC2", 100.0), _region("GC4", 100.0)]
        with pytest.raises(OverlappingRegionsError) as info:
            assemble(LIB, regions, CHARS, region="GC4", **KW)
        text = str(info.value)
        assert "PRECIP" in text and "skipped" in text
        assert "is not one of the overlapping" not in text
        assert [c[0] for c in info.value.candidates] == ["GC1", "GC2"]

    def test_choice_of_a_zero_weight_region_raises(self):
        regions = [_region("GC1", 100.0), _region("GC2", 0.0)]
        with pytest.raises(EquationsUnavailable, match="GC2.*0 %"):
            assemble(LIB, regions, CHARS, region="GC2", **KW)

    def test_choice_of_an_unmatched_region_raises(self):
        regions = [_region("GC1", 100.0), _region("GC9", 100.0)]
        with pytest.raises(EquationsUnavailable, match="GC9.*no library region"):
            assemble(LIB, regions, CHARS, region="GC9", **KW)

    def test_unknown_choice_on_a_partition_raises_listing_located(self):
        regions = [_region("GC1", 60.0), _region("GC2", 40.0)]
        with pytest.raises(EquationsUnavailable, match="GC7") as info:
            assemble(LIB, regions, CHARS, region="GC7", **KW)
        assert "GC1" in str(info.value) and "GC2" in str(info.value)
        assert not isinstance(info.value, OverlappingRegionsError)

    def test_choice_in_use_on_a_partition_is_a_silent_no_op(self, caplog):
        regions = [_region("GC1", 60.0), _region("GC2", 40.0)]
        with caplog.at_level("WARNING", logger="flowfreq.ungaged"):
            a = assemble(LIB, regions, CHARS, region="R2", **KW)
        assert [r.code for r in a.used] == ["GC1", "GC2"] and a.alternatives == []
        assert [f for f, _ in a.estimates[0.1]] == pytest.approx([0.6, 0.4])
        assert not [r for r in caplog.records if "region=" in r.getMessage()]

    def test_to_markdown_lists_alternatives(self):
        regions = [_region("GC1", 100.0), _region("GC2", 100.0)]
        a = assemble(LIB, regions, CHARS, region="GC1", **KW)
        est = _estimate(
            regions,
            assembly=a,
            provenance={
                "alternative_regions": [
                    {
                        "code": "GC2",
                        "name": "Located GC2",
                        "region_code": "R2",
                        "citation": "Synthetic report two (2021), table T2",
                        "year": 2021,
                    }
                ]
            },
        )
        text = est.to_markdown()
        assert "| GC1 | Located GC1 | 1.000 |" in text
        assert "Not used (alternative NSS equation sets)" in text
        assert "GC2 (Located GC2): Synthetic report two (2021), table T2" in text
        assert "area-weighted" not in text


# --------------------------------------------------------------- estimate_at_site
#
# StreamStats is faked at the two names ``estimate_at_site`` calls, with the WA entry of
# the wave 1/2 vignette capture (Little Klickitat River near Goldendale, 14112000): a
# basin NSS places 0.92% in region 1 (GC1750) and 99.08% in region 4 (GC1753) of the
# verified WA library.

VIGNETTE = Path(__file__).parent.parent / "docs" / "vignettes" / "wave1_wave2_sites.json"
WA_NSS_CODES = {"1": "GC1750", "4": "GC1753"}


def _wa_site() -> Dict[str, Any]:
    site: Dict[str, Any] = json.loads(VIGNETTE.read_text(encoding="utf-8"))["WA"]
    return site


def _watershed(site: Dict[str, Any], lat: float, lon: float) -> WatershedCharacteristics:
    return WatershedCharacteristics(
        region=site["streamstats_region"],
        snap=SnapResult(lat, lon, lat + 1e-4, lon - 1e-4, True, 13.2),
        polygon_geojson={
            "type": "Feature",
            "geometry": {"type": "Polygon", "coordinates": [[[lon, lat], [lon + 0.1, lat]]]},
        },
        characteristics={
            code: Characteristic(code, code, "", value, "unit")
            for code, value in site["characteristics"].items()
        },
        provenance=Provenance(
            service_versions={"ss-hydro": "1.0"},
            request_urls=["https://example.invalid/ss-hydro"],
            requested_at_utc="2026-10-04T00:00:00+00:00",
            server_used="prodweba",
        ),
        polygon_area_sq_mi=site["characteristics"]["DRNAREA"],
        unavailable={"ELEV": "not found"},
    )


class Fakes:
    """Recording stand-ins for the two StreamStats calls."""

    def __init__(self, site: Dict[str, Any]) -> None:
        self.site = site
        self.regions: List[RegressionRegionLocation] = [
            RegressionRegionLocation(
                region_id=i,
                code=WA_NSS_CODES[code],
                name=f"Regression Region {code}",
                percent_weight=100.0 * weight,
                area_sq_mi=weight * site["characteristics"]["DRNAREA"],
            )
            for i, (code, weight) in enumerate(site["region_weights"].items(), start=717)
        ]
        self.delineate_calls: List[Dict[str, Any]] = []
        self.locate_calls: List[Dict[str, Any]] = []
        self.watershed: Optional[WatershedCharacteristics] = None

    def delineate(self, region: str, lat: float, lon: float, **kw: Any) -> WatershedCharacteristics:
        self.delineate_calls.append({"region": region, "lat": lat, "lon": lon, **kw})
        self.watershed = _watershed(self.site, lat, lon)
        return self.watershed

    def locate(
        self, region: str, polygon: Dict[str, Any], **kw: Any
    ) -> List[RegressionRegionLocation]:
        self.locate_calls.append({"region": region, "polygon": polygon, **kw})
        return list(self.regions)


@pytest.fixture
def fakes(monkeypatch) -> Fakes:
    f = Fakes(_wa_site())
    monkeypatch.setattr("flowfreq.ungaged.delineate_and_get_characteristics", f.delineate)
    monkeypatch.setattr("flowfreq.ungaged.locate_regression_regions", f.locate)
    return f


@pytest.fixture
def tmp_data_dir(tmp_path: Path) -> Path:
    """A library directory holding only an ``nss_snapshot`` Tennessee set (one region)."""
    path = tmp_path / NSS_SNAPSHOT_SUBDIR / "TN.json"
    path.parent.mkdir(parents=True)
    data = {
        "schema_version": 1,
        "state": "TN",
        "status": "nss_snapshot",
        "source_reports": [{"citation": "Synthetic report"}],
        "retrieved_utc": "2026-10-25T17:20:40+00:00",
        "equations": [
            {
                "region_code": "GC2015",
                "aep": 0.01,
                "intercept": 2.0,
                "variables": [{"code": "DRNAREA"}],
                "coefficients": [0.7],
                "citation": {"publication": "Synthetic", "table": "T1"},
            }
        ],
    }
    path.write_text(json.dumps(data), encoding="utf-8")
    return tmp_path


@pytest.fixture(scope="module")
def big_sandy_results():
    """A real Big Sandy EMA fit, whose confidence_limits carry ``variance_log``."""
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


WA_LAT, WA_LON = 45.844124, -120.796392


class TestEstimateAtSite:
    def test_washington_site_uses_the_verified_library_and_two_regions(self, fakes):
        from flowfreq.regression.library import load_state

        est = estimate_at_site(WA_LAT, WA_LON, state="WA", cache=None, timeout=12.5)
        assert isinstance(est, SiteEstimate)
        assert (est.lat, est.lon, est.state) == (WA_LAT, WA_LON, "WA")
        assert est.method == "regression"
        assert est.library_status == "verified" and est.advisory == ""
        assert est.weighted is None
        # The two StreamStats calls receive exactly what the caller passed.
        assert fakes.delineate_calls == [
            {"region": "WA", "lat": WA_LAT, "lon": WA_LON, "cache": None, "timeout": 12.5}
        ]
        assert fakes.locate_calls == [
            {"region": "WA", "polygon": fakes.watershed.polygon_geojson, "timeout": 12.5}
        ]
        # Both located GC codes map to library regions 1 and 4 through their names.
        assert est.regions == fakes.regions
        assert est.unmatched_regions == []
        assert est.region_codes == {"GC1750": "1", "GC1753": "4"}
        assert [est.region_codes[r.code] for r in est.regions] == ["1", "4"]
        assert est.characteristics == fakes.site["characteristics"]
        # Default AEPs: everything the WA library publishes, largest first.
        lib = load_state("WA")
        assert list(est.quantiles["aep"]) == sorted(lib.aeps("4"), reverse=True)
        assert list(est.estimates) == list(est.quantiles["aep"])
        fractions = [f for f, _ in est.estimates[0.01]]
        assert fractions == pytest.approx([0.0092, 0.9908])
        assert est.quantiles.set_index("aep").loc[0.01, "regions"] == "GC1750, GC1753"
        assert (
            est.quantiles.set_index("aep").loc[0.01, "interval_method"]
            == "none: area-weighted across 2 regions"
        )
        assert not est.quantiles["extrapolated"].any()
        text = est.to_markdown()
        assert "| GC1753 | Regression Region 4 | 0.991 |" in text

    def test_lower_case_state_is_accepted(self, fakes):
        est = estimate_at_site(WA_LAT, WA_LON, state="wa")
        assert est.provenance["library"]["status"] == "verified"

    def test_estimate_carries_the_watershed_polygon(self, fakes):
        est = estimate_at_site(WA_LAT, WA_LON, state="WA")
        # The same validated Feature StreamStats returned, for callers that draw it.
        assert est.basin_geojson is fakes.watershed.polygon_geojson
        assert est.basin_geojson["type"] == "Feature"
        # The report stays free of the polygon's vertices.
        assert "coordinates" not in est.to_markdown()

    def test_requested_aeps_and_options_reach_the_assembly(self, fakes):
        est = estimate_at_site(
            WA_LAT,
            WA_LON,
            state="WA",
            aeps=[0.01, 0.5],
            weighting_space="linear",
            allow_extrapolation=True,
        )
        assert list(est.quantiles["aep"]) == [0.01, 0.5]
        q1 = est.estimates[0.01][0][1].flow_cfs
        q4 = est.estimates[0.01][1][1].flow_cfs
        assert est.quantiles.loc[0, "flow_cfs"] == pytest.approx(0.0092 * q1 + 0.9908 * q4)
        assert est.provenance["options"] == {
            "allow_nss_snapshot": False,
            "allow_provisional": False,
            "allow_extrapolation": True,
            "weighting_space": "linear",
            "region": None,
        }

    def test_snapshot_state_needs_the_flag(self, fakes, tmp_data_dir):
        fakes.regions = [
            RegressionRegionLocation(
                region_id=1, code="GC2015", name="Region 1", percent_weight=100.0, area_sq_mi=82.2
            )
        ]
        with pytest.raises(
            EquationsUnavailable,
            match=r"^TN: no regression library file\. NSS serves peak-flow regions for this "
            r"jurisdiction",
        ) as info:
            estimate_at_site(36.0, -86.0, state="TN", data_dir=tmp_data_dir)
        assert "allow_nss_snapshot=True" in str(info.value)
        assert fakes.delineate_calls == []  # refused before any network call

        est = estimate_at_site(
            36.0, -86.0, state="TN", data_dir=tmp_data_dir, allow_nss_snapshot=True
        )
        assert est.library_status == "nss_snapshot"
        assert est.advisory.startswith("NSS SNAPSHOT (not report-verified): the TN equations")
        assert est.provenance["library"] == {
            "status": "nss_snapshot",
            "file": "nss_snapshot/TN.json",
            "retrieved_utc": "2026-10-25T17:20:40+00:00",
        }
        assert list(est.quantiles["aep"]) == [0.01]
        assert est.quantiles.loc[0, "flow_cfs"] == pytest.approx(
            10 ** (2.0 + 0.7 * math.log10(82.2))
        )
        assert "Advisory: NSS SNAPSHOT" in est.to_markdown()

    def test_no_library_file_names_nss_coverage(self, fakes):
        with pytest.raises(
            EquationsUnavailable,
            match=r"^DC: no regression library file\. NSS serves no peak-flow regions for "
            r"this jurisdiction$",
        ):
            estimate_at_site(38.9, -77.0, state="DC")
        assert fakes.delineate_calls == [] and fakes.locate_calls == []

    def test_no_library_file_in_an_override_directory(self, fakes, tmp_data_dir):
        # With data_dir the coverage check reads that directory: TN has a snapshot file
        # there and WA has nothing, however many packaged files say otherwise.
        with pytest.raises(EquationsUnavailable, match="WA: .*NSS serves no peak-flow"):
            estimate_at_site(WA_LAT, WA_LON, state="WA", data_dir=tmp_data_dir)

    def test_at_site_results_are_weighted_in_a_single_region(self, fakes, big_sandy_results):
        fakes.regions = [fakes.regions[1]]  # GC1753 alone
        est = estimate_at_site(WA_LAT, WA_LON, state="WA", at_site=big_sandy_results)
        assert est.method == "at-site weighted with regression"
        assert [r.code for r in est.regions] == ["GC1753"]
        assert est.weighted is not None
        w = est.weighted.set_index("aep")
        assert "variance_log" in w.columns and "weighted_flow_cfs" in w.columns
        # One row per at-site AEP; the regression's AEPs are weighted, the rest noted.
        assert list(est.weighted["aep"]) == list(big_sandy_results.confidence_limits["aep"])
        assert np.isfinite(w.loc[0.01, "weighted_flow_cfs"])
        assert w.loc[0.01, "regression_flow_cfs"] == pytest.approx(
            est.quantiles.set_index("aep").loc[0.01, "flow_cfs"]
        )
        assert w.loc[0.995, "note"] == "no regression estimate at this AEP"
        assert est.provenance["weighting"] == "at-site fit weighted with regression region GC1753"
        assert "## At-site weighted with regression" in est.to_markdown()

    def test_at_site_in_two_regions_is_not_weighted(self, fakes, big_sandy_results):
        est = estimate_at_site(WA_LAT, WA_LON, state="WA", at_site=big_sandy_results)
        assert est.method == "regression"
        assert est.weighted is None
        assert est.provenance["weighting"] == (
            "at-site fit not weighted: the basin spans 2 regression regions (GC1750, GC1753), "
            "whose area-weighted flow carries no variance of prediction"
        )

    def test_no_at_site_fit_says_so(self, fakes):
        est = estimate_at_site(WA_LAT, WA_LON, state="WA")
        assert est.provenance["weighting"] == "no at-site fit supplied"

    def test_provenance_carries_snap_delineation_and_status(self, fakes):
        est = estimate_at_site(WA_LAT, WA_LON, state="WA")
        p = est.provenance
        ws = fakes.watershed
        assert p["snap"] == ws.snap.to_dict()
        assert p["snap"]["snapped_lat"] == pytest.approx(WA_LAT + 1e-4)
        assert p["delineation"] == {
            "polygon_area_sq_mi": 82.2,
            "unavailable": {"ELEV": "not found"},
            "streamstats": ws.provenance.to_dict(),
        }
        assert p["library"] == {"status": "verified", "file": "WA.json", "retrieved_utc": ""}
        assert p["regions_located"] == [r.to_dict() for r in fakes.regions]
        assert p["missing_characteristics"] == {} and p["zero_weight_regions"] == []
        assert p["out_of_range"] == {}
        assert p["estimated_at_utc"].endswith("+00:00")
        assert p["weighting"] == "no at-site fit supplied"

    def test_provenance_records_out_of_range_inputs(self, fakes):
        fakes.site["characteristics"]["DRNAREA"] = 50000.0
        est = estimate_at_site(
            WA_LAT, WA_LON, state="WA", aeps=[0.01, 0.5], allow_extrapolation=True
        )
        assert est.quantiles["extrapolated"].all()
        out = est.provenance["out_of_range"]
        assert set(out) == {0.01, 0.5}
        assert set(out[0.01]) == {"GC1750", "GC1753"}
        assert out[0.01]["GC1753"] == list(est.estimates[0.01][1][1].out_of_range)
        assert any("DRNAREA=50000.0 > max" in s for s in out[0.01]["GC1753"])

    def test_streamstats_errors_propagate_unchanged(self, fakes, monkeypatch):
        def unsnappable(region, lat, lon, **kw):
            raise UnsnappablePointError("no stream within tolerance")

        monkeypatch.setattr("flowfreq.ungaged.delineate_and_get_characteristics", unsnappable)
        with pytest.raises(UnsnappablePointError, match="no stream within tolerance"):
            estimate_at_site(WA_LAT, WA_LON, state="WA")

    def test_a_delineation_without_a_polygon_is_refused_before_locating(self, fakes, monkeypatch):
        def no_polygon(region, lat, lon, **kw):
            ws = fakes.delineate(region, lat, lon, **kw)
            ws.polygon_geojson = None
            return ws

        monkeypatch.setattr("flowfreq.ungaged.delineate_and_get_characteristics", no_polygon)
        with pytest.raises(ValueError, match="WA: StreamStats returned no watershed polygon"):
            estimate_at_site(WA_LAT, WA_LON, state="WA")
        assert fakes.locate_calls == []

    def test_unknown_weighting_space_is_refused_before_any_call(self, fakes, monkeypatch):
        def no_library(*args: Any, **kw: Any) -> None:
            raise AssertionError("the library must not be loaded")

        monkeypatch.setattr("flowfreq.ungaged._load_library", no_library)
        with pytest.raises(ValueError, match="weighting_space"):
            estimate_at_site(WA_LAT, WA_LON, state="WA", weighting_space="geometric")
        assert fakes.delineate_calls == [] and fakes.locate_calls == []

    def test_empty_aeps_are_refused_before_any_call(self, fakes):
        with pytest.raises(ValueError, match="aeps must be non-empty"):
            estimate_at_site(WA_LAT, WA_LON, state="WA", aeps=[])
        assert fakes.delineate_calls == [] and fakes.locate_calls == []

    def test_at_site_without_variance_log_is_refused_before_streamstats(
        self, fakes, big_sandy_results
    ):
        import copy

        stale = copy.copy(big_sandy_results)
        stale.confidence_limits = big_sandy_results.confidence_limits.drop(columns=["variance_log"])
        with pytest.raises(ValueError, match="no 'variance_log' column"):
            estimate_at_site(WA_LAT, WA_LON, state="WA", at_site=stale)
        assert fakes.delineate_calls == [] and fakes.locate_calls == []

    def test_overlapping_regions_need_a_choice(self, fakes):
        for r in fakes.regions:
            r.percent_weight = 100.0
        with pytest.raises(OverlappingRegionsError, match="GC1750"):
            estimate_at_site(WA_LAT, WA_LON, state="WA")

    def test_region_choice_reaches_the_assembly_and_provenance(self, fakes, big_sandy_results):
        for r in fakes.regions:
            r.percent_weight = 100.0
        est = estimate_at_site(
            WA_LAT, WA_LON, state="WA", region="GC1753", at_site=big_sandy_results
        )
        assert [r.code for r in est.regions] == ["GC1753"]
        assert est.method == "at-site weighted with regression"
        (alt,) = est.provenance["alternative_regions"]
        assert alt["code"] == "GC1750" and alt["region_code"] == "1"
        assert alt["citation"] and alt["year"]
        assert est.provenance["options"]["region"] == "GC1753"
        assert "Not used (alternative NSS equation sets)" in est.to_markdown()

    def test_region_in_use_on_a_partition_is_recorded_as_given(self, fakes):
        est = estimate_at_site(WA_LAT, WA_LON, state="WA", region="GC1753")
        assert est.provenance["options"]["region"] == "GC1753"
        assert est.provenance["alternative_regions"] == []

    def test_region_not_in_use_never_reaches_provenance(self, fakes):
        with pytest.raises(EquationsUnavailable, match="GC9999"):
            estimate_at_site(WA_LAT, WA_LON, state="WA", region="GC9999")

    def test_no_overlap_leaves_alternatives_empty(self, fakes):
        est = estimate_at_site(WA_LAT, WA_LON, state="WA")
        assert est.provenance["alternative_regions"] == []

    def test_exported_from_the_package(self):
        import flowfreq

        assert flowfreq.estimate_at_site is estimate_at_site
        assert flowfreq.SiteEstimate is SiteEstimate
        assert "estimate_at_site" in flowfreq.__all__ and "SiteEstimate" in flowfreq.__all__


def test_overlapping_regions_error_is_exported_from_package() -> None:
    from flowfreq import OverlappingRegionsError as exported

    assert exported is OverlappingRegionsError


# ------------------------------------------------------------ outlet coordinates
#
# StreamStats does not return the outlet coordinates some equations use (WY region 1's
# LONG_OUT, regions 4 and 6's LAT_OUT; MN's LAT_OUT/LONG_OUT), so estimate_at_site takes
# them from the snapped pour point. Reports differ on the sign of a west longitude: WY's
# range is 105.21-111.34 (degrees west, positive), MN's -96.12 to -91.66 (signed).


class TestOutletCoordinates:
    def test_wyoming_takes_degrees_west_as_positive(self):
        from flowfreq.regression.library import load_state

        lib = load_state("WY", warn=False)
        out = outlet_coordinates(lib, 43.1, -108.4, {})
        assert out == {"LAT_OUT": pytest.approx(43.1), "LONG_OUT": pytest.approx(108.4)}

    def test_minnesota_keeps_the_sign(self):
        from flowfreq.regression.library import load_state

        lib = load_state("MN", allow_nss_snapshot=True, warn=False)
        out = outlet_coordinates(lib, 45.0, -93.5, {})
        assert out == {"LAT_OUT": pytest.approx(45.0), "LONG_OUT": pytest.approx(-93.5)}

    def test_a_value_streamstats_returned_is_kept(self):
        from flowfreq.regression.library import load_state

        lib = load_state("WY", warn=False)
        assert outlet_coordinates(lib, 43.1, -108.4, {"LONG_OUT": 108.0}) == {
            "LAT_OUT": pytest.approx(43.1)
        }

    def test_a_library_without_outlet_codes_gets_nothing(self):
        from flowfreq.regression.library import load_state

        assert outlet_coordinates(load_state("WA", warn=False), 46.0, -120.0, {}) == {}


def _outlet_library(tmp_path: Path) -> Path:
    """A WA ``nss_snapshot`` set whose equations use LONG_OUT as degrees west."""
    path = tmp_path / NSS_SNAPSHOT_SUBDIR / "WA.json"
    path.parent.mkdir(parents=True)
    eq = {
        "aep": 0.01,
        "intercept": 1.0,
        "variables": [
            {"code": "DRNAREA", "transform": "log10"},
            {"code": "LONG_OUT", "transform": "log10", "minimum": 100.0, "maximum": 130.0,
             "offset": -100.0},
        ],
        "coefficients": [0.7, 0.5],
        "citation": {"publication": "Synthetic", "table": "T1"},
    }  # fmt: skip
    data = {
        "schema_version": 1,
        "state": "WA",
        "status": "nss_snapshot",
        "source_reports": [{"citation": "Synthetic report"}],
        "retrieved_utc": "2026-10-30T00:00:00+00:00",
        "equations": [{**eq, "region_code": code} for code in ("GC1750", "GC1753")],
    }
    path.write_text(json.dumps(data), encoding="utf-8")
    return tmp_path


def test_estimate_at_site_derives_long_out_from_the_pour_point(fakes, tmp_path):
    data_dir = _outlet_library(tmp_path)
    est = estimate_at_site(
        WA_LAT, WA_LON, state="WA", data_dir=data_dir, allow_nss_snapshot=True, cache=None
    )
    snapped_lon = WA_LON - 1e-4  # the fake's snap
    assert est.provenance["derived_characteristics"] == {
        "LONG_OUT": {"value": pytest.approx(-snapped_lon), "source": "snapped pour point"}
    }
    area = fakes.site["characteristics"]["DRNAREA"]
    expected = 10 ** (1.0 + 0.7 * math.log10(area) + 0.5 * math.log10(-snapped_lon - 100.0))
    q = est.quantiles
    assert float(q.loc[q["aep"] == 0.01, "flow_cfs"].iloc[0]) == pytest.approx(expected)


def _area_class_library(small_max: float = 1000.0, large_min: float = 1000.0) -> StateLibrary:
    small = Variable("DRNAREA", "log10", "mi2", 0.4, small_max)
    large = Variable("DRNAREA", "log10", "mi2", large_min, 31100.0)
    eqs = []
    for aep in AEPS:
        eqs.append(_eq("GC56", aep, 2.0, variables=(small,), region_name="Areas under 1000"))
        eqs.append(_eq("GC57", aep, 2.5, variables=(large,), region_name="Areas over 1000"))
    return StateLibrary(state="AK", status="verified", equations=eqs, advisory="")


class TestDrainageAreaClasses:
    """Overlapping regions that split the state by drainage area (AK SIR 2016-5024's
    GC1656 under 1,000 mi2 and GC1657 over it, both located statewide) are chosen by
    the basin's DRNAREA rather than refused."""

    REGIONS = [_region("GC56", 100.0), _region("GC57", 100.0)]

    @pytest.mark.parametrize("area,code", [(62.0, "GC56"), (5000.0, "GC57")])
    def test_the_basin_area_picks_the_class(self, area: float, code: str) -> None:
        a = assemble(_area_class_library(), self.REGIONS, {"DRNAREA": area}, **KW)
        assert [r.code for r in a.used] == [code]
        assert [r.code for r in a.alternatives] == [{"GC56": "GC57", "GC57": "GC56"}[code]]
        assert a.region_choice_source is not None
        assert "drainage-area classes" in a.region_choice_source and code in a.region_choice_source
        assert not a.quantiles["interval_method"].str.startswith("none").any()

    def test_a_shared_endpoint_is_left_to_the_caller(self) -> None:
        with pytest.raises(OverlappingRegionsError):
            assemble(_area_class_library(), self.REGIONS, {"DRNAREA": 1000.0}, **KW)
        a = assemble(_area_class_library(), self.REGIONS, {"DRNAREA": 1000.0}, region="GC57", **KW)
        assert [r.code for r in a.used] == ["GC57"]

    def test_overlapping_ranges_are_alternatives_not_classes(self) -> None:
        lib = _area_class_library(small_max=2000.0, large_min=500.0)
        with pytest.raises(OverlappingRegionsError):
            assemble(lib, self.REGIONS, {"DRNAREA": 62.0}, **KW)

    def test_a_partial_cover_is_not_decided_by_area(self) -> None:
        regions = [_region("GC56", 60.0), _region("GC57", 100.0)]
        with pytest.raises(OverlappingRegionsError):
            assemble(_area_class_library(), regions, {"DRNAREA": 62.0}, **KW)

    def test_an_explicit_choice_still_wins(self) -> None:
        a = assemble(_area_class_library(), self.REGIONS, {"DRNAREA": 62.0}, region="GC56", **KW)
        assert [r.code for r in a.used] == ["GC56"]
        assert a.region_choice_source is None


class TestRegionGates:
    """Overlapping regions chosen by a region-selection rule on a characteristic that is
    not a regression variable (Arizona: High Elevation at a mean elevation of 7,500 ft or
    more, SIR 2014-5211)."""

    REGIONS = [_region("GC18", 100.0), _region("GC22", 100.0)]

    @staticmethod
    def _lib(gates=None) -> StateLibrary:
        eqs = []
        for aep in AEPS:
            eqs.append(_eq("1", aep, 2.0, region_name="Flood region 1 (High Elevation; NSS GC18)"))
            eqs.append(_eq("5", aep, 3.0, region_name="Flood region 5 (NSS GC22)"))
        lib = StateLibrary(state="AZ", status="verified", equations=eqs, advisory="")
        lib.region_gates = (
            {
                "1": RegionGate("ELEV", minimum=7500.0, units="ft", source="report"),
                "5": RegionGate("ELEV", maximum=7500.0, units="ft", source="report"),
            }
            if gates is None
            else gates
        )
        return lib

    @pytest.mark.parametrize("elev,code", [(5976.0, "GC22"), (7500.0, "GC18"), (9000.0, "GC18")])
    def test_the_gate_picks_the_region(self, elev, code):
        a = assemble(self._lib(), self.REGIONS, {"DRNAREA": 35.5, "ELEV": elev}, **KW)
        assert [r.code for r in a.used] == [code]
        assert a.region_choice_source.startswith("region gate: ELEV")
        assert len(a.alternatives) == 1

    def test_without_the_characteristic_the_caller_chooses(self):
        with pytest.raises(OverlappingRegionsError):
            assemble(self._lib(), self.REGIONS, {"DRNAREA": 35.5}, **KW)

    def test_an_ungated_candidate_leaves_it_to_the_caller(self):
        lib = self._lib({"1": RegionGate("ELEV", minimum=7500.0)})
        with pytest.raises(OverlappingRegionsError):
            assemble(lib, self.REGIONS, {"DRNAREA": 35.5, "ELEV": 5000.0}, **KW)

    def test_an_explicit_choice_still_wins(self):
        a = assemble(
            self._lib(), self.REGIONS, {"DRNAREA": 35.5, "ELEV": 5000.0}, region="GC18", **KW
        )
        assert [r.code for r in a.used] == ["GC18"] and a.region_choice_source is None

    def test_arizona_sabino_creek(self):
        from flowfreq.regression.library import load_state

        lib = load_state("AZ")
        regions = [
            RegressionRegionLocation(0, "GC1618", "high", 100.0, 35.5),
            RegressionRegionLocation(0, "GC1622", "se", 100.0, 35.5),
        ]
        a = assemble(lib, regions, {"CONTDA": 35.5, "ELEV": 5976.0}, **{**KW, "aeps": [0.01]})
        assert [r.code for r in a.used] == ["GC1622"]
        assert a.quantiles["flow_cfs"].iloc[0] == pytest.approx(8604, rel=0.01)

    OREGON = {
        "DRNAREA": 50.0,
        "BSLOPD": 15.0,
        "I24H2Y": 3.0,
        "JANMINT2K": 27.0,
        "JANMAXT2K": 40.0,
    }
    OREGON_REGIONS = [
        RegressionRegionLocation(0, "GC730", "2A", 100.0, 50.0),
        RegressionRegionLocation(0, "GC731", "2B", 100.0, 50.0),
    ]

    @pytest.mark.parametrize(
        "elev,code,lib_code", [(2500.0, "GC731", "2B"), (3500.0, "GC730", "2A")]
    )
    def test_oregon_region2_outside_the_transition_zone(self, elev, code, lib_code):
        """Outside 2,875-3,125 ft one region applies alone, as eq. 8's blend gives."""
        from flowfreq.regression.library import load_state
        from flowfreq.regression.oregon import estimate_region2

        chars = {**self.OREGON, "ELEV": elev}
        a = assemble(load_state("OR"), self.OREGON_REGIONS, chars, **{**KW, "aeps": [0.01]})
        assert [r.code for r in a.used] == [code]
        blend = estimate_region2(0.01, chars)
        assert blend.weights[lib_code] == 1.0
        assert a.quantiles["flow_cfs"].iloc[0] == pytest.approx(blend.flow_cfs, rel=1e-9)

    def test_oregon_region2_transition_zone_names_the_blend(self):
        from flowfreq.regression.library import load_state

        chars = {**self.OREGON, "ELEV": 3000.0}
        with pytest.raises(
            OverlappingRegionsError, match="no region gate admits.*estimate_region2"
        ):
            assemble(load_state("OR"), self.OREGON_REGIONS, chars, **KW)


def test_region_gate_rules():
    gate = RegionGate("ELEV", minimum=7500.0, units="ft")
    assert gate.admits(7500.0) and not gate.admits(7499.9)
    assert gate.describe() == "ELEV >= 7500 ft"
    below = RegionGate("ELEV", maximum=7500.0)
    assert below.admits(7499.9) and not below.admits(7500.0)
    with pytest.raises(ValueError, match="no minimum or maximum"):
        RegionGate.from_dict({"variable": "ELEV"})


def test_parse_state_refuses_gates_for_unknown_regions():
    from flowfreq.regression.library import parse_state

    data = json.loads(
        (Path(__file__).resolve().parent.parent / "flowfreq/data/regression/AZ.json").read_text(
            encoding="utf-8"
        )
    )
    assert set(parse_state(data).region_gates) == {"1", "2", "3", "4", "5"}
    data["region_gates"]["9"] = {"variable": "ELEV", "minimum": 1.0}
    with pytest.raises(ValueError, match="no equations"):
        parse_state(data)
