"""Tests for ``flowfreq.ungaged``: ``SiteEstimate`` and the pure region assembly.

The library is synthetic: ``GC1`` matches a located NSS region by ``region_code``,
``R2`` matches ``GC2`` through its ``region_name`` (as the verified files write
"... (NSS GC1239)"), and ``R4`` needs a characteristic (``PRECIP``) the basin may
not carry.
"""

from __future__ import annotations

import math

import pandas as pd
import pytest

from flowfreq.regression.equations import (
    Citation,
    OutOfRangeError,
    RegressionEquation,
    RegressionEstimate,
    Variable,
)
from flowfreq.regression.library import EquationsUnavailable, StateLibrary
from flowfreq.streamstats import RegressionRegionLocation
from flowfreq.ungaged import SiteEstimate, assemble, region_for_nss_code

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
        frame, estimates, unmatched, missing = assemble(LIB, [_region("GC1", 100.0)], CHARS, **KW)
        assert unmatched == [] and missing == {}
        row = frame.set_index("aep").loc[0.01]
        ((fraction, est),) = estimates[0.01]
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
        frame, estimates, unmatched, missing = assemble(LIB, regions, CHARS, **KW)
        assert unmatched == [] and missing == {}
        row = frame.set_index("aep").loc[0.1]
        q1, q2 = _flow("GC1", 0.1), _flow("R2", 0.1)
        expected = 10 ** (0.6 * math.log10(q1) + 0.4 * math.log10(q2))
        assert row["flow_cfs"] == pytest.approx(expected)
        assert pd.isna(row["lower"]) and pd.isna(row["upper"])
        assert pd.isna(row["interval_level"])
        assert row["interval_method"] == "none: area-weighted across 2 regions"
        assert row["regions"] == "GC1, GC2"
        assert [(f, e.equation.region_code) for f, e in estimates[0.1]] == [
            (0.6, "GC1"),
            (0.4, "R2"),
        ]

    def test_linear_weighting_space(self):
        regions = [_region("GC1", 60.0), _region("GC2", 40.0)]
        frame, *_ = assemble(
            LIB, regions, CHARS, **{**KW, "aeps": [0.1], "weighting_space": "linear"}
        )
        expected = 0.6 * _flow("GC1", 0.1) + 0.4 * _flow("R2", 0.1)
        assert frame.loc[0, "flow_cfs"] == pytest.approx(expected)

    def test_fractions_renormalise_after_an_unmatched_region(self):
        regions = [_region("GC1", 60.0), _region("GC2", 30.0), _region("GC3", 10.0)]
        frame, estimates, unmatched, missing = assemble(
            LIB, regions, CHARS, **{**KW, "aeps": [0.01]}
        )
        assert unmatched == [regions[2]]
        assert missing == {}
        fractions = [f for f, _ in estimates[0.01]]
        assert fractions == pytest.approx([2 / 3, 1 / 3])
        assert frame.loc[0, "regions"] == "GC1, GC2"

    def test_no_matching_region_raises_naming_the_codes(self):
        regions = [_region("GC9", 70.0), _region("GC8", 30.0)]
        with pytest.raises(EquationsUnavailable, match="GC9"):
            assemble(LIB, regions, CHARS, **KW)
        with pytest.raises(EquationsUnavailable, match="GC8"):
            assemble(LIB, regions, CHARS, **KW)

    def test_missing_characteristic_skips_the_region_and_records_it(self):
        regions = [_region("GC1", 50.0), _region("GC4", 50.0)]
        frame, estimates, unmatched, missing = assemble(LIB, regions, CHARS, **KW)
        assert unmatched == []
        assert missing == {"GC4": ["PRECIP"]}
        # GC1 carries the whole basin once GC4 is skipped, interval included.
        assert [(f, e.equation.region_code) for f, e in estimates[0.5]] == [(1.0, "GC1")]
        assert frame.loc[0, "regions"] == "GC1"
        assert frame.loc[0, "interval_method"] == estimates[0.5][0][1].interval_method
        assert frame.loc[0, "lower"] == pytest.approx(estimates[0.5][0][1].interval[0])

    def test_every_region_missing_a_characteristic_raises(self):
        with pytest.raises(KeyError, match="PRECIP"):
            assemble(LIB, [_region("GC4", 100.0)], CHARS, **KW)

    def test_extrapolated_flag_follows_out_of_range(self):
        big = {"DRNAREA": 5000.0}
        with pytest.raises(OutOfRangeError):
            assemble(LIB, [_region("GC1", 100.0)], big, **{**KW, "aeps": [0.01]})
        regions = [_region("GC1", 60.0), _region("GC2", 40.0)]
        frame, estimates, _, _ = assemble(
            LIB, regions, big, aeps=[0.01], allow_extrapolation=True, weighting_space="log"
        )
        assert bool(frame.loc[0, "extrapolated"]) is True
        assert all(e.extrapolated for _, e in estimates[0.01])
        in_range, *_ = assemble(
            LIB,
            [_region("GC1", 100.0)],
            CHARS,
            aeps=[0.01],
            allow_extrapolation=True,
            weighting_space="log",
        )
        assert bool(in_range.loc[0, "extrapolated"]) is False

    def test_default_aeps_are_the_intersection_largest_first(self):
        regions = [_region("GC1", 60.0), _region("GC2", 40.0)]
        frame, estimates, _, _ = assemble(LIB, regions, CHARS, **KW)
        assert list(frame["aep"]) == [0.5, 0.1, 0.01]
        assert list(estimates) == [0.5, 0.1, 0.01]
        alone, *_ = assemble(LIB, [_region("GC1", 100.0)], CHARS, **KW)
        assert list(alone["aep"]) == [0.5, 0.1, 0.01, 0.002]

    def test_requested_aeps_keep_the_caller_order(self):
        frame, *_ = assemble(LIB, [_region("GC1", 100.0)], CHARS, **{**KW, "aeps": [0.01, 0.5]})
        assert list(frame["aep"]) == [0.01, 0.5]

    def test_unpublished_requested_aep_raises(self):
        regions = [_region("GC1", 60.0), _region("GC2", 40.0)]
        with pytest.raises(EquationsUnavailable, match="0.002"):
            assemble(LIB, regions, CHARS, **{**KW, "aeps": [0.002]})

    def test_frame_columns(self):
        frame, *_ = assemble(LIB, [_region("GC1", 100.0)], CHARS, **KW)
        assert list(frame.columns) == [
            "aep",
            "flow_cfs",
            "lower",
            "upper",
            "interval_level",
            "interval_method",
            "extrapolated",
            "regions",
        ]


def _estimate(regions, **kw) -> SiteEstimate:
    frame, estimates, unmatched, missing = assemble(LIB, regions, CHARS, **KW)
    base = dict(
        lat=45.0,
        lon=-110.0,
        state="XX",
        method="regression",
        quantiles=frame,
        regions=[r for r in regions if r not in unmatched and r.code not in missing],
        unmatched_regions=unmatched,
        characteristics=dict(CHARS),
        library_status=LIB.status,
        advisory="Synthetic advisory text.",
        estimates=estimates,
        provenance={"missing_characteristics": missing},
    )
    base.update(kw)
    return SiteEstimate(**base)


class TestSiteEstimate:
    def test_to_markdown_names_regions_fractions_and_citations(self):
        regions = [_region("GC1", 60.0), _region("GC2", 30.0), _region("GC3", 10.0)]
        text = _estimate(regions).to_markdown()
        assert "45.0" in text and "-110.0" in text and "XX" in text
        assert "GC1" in text and "Located GC1" in text
        assert "0.667" in text and "0.333" in text
        assert "Synthetic report one" in text and "Synthetic report two" in text
        assert "GC3" in text  # the unmatched region is reported, not silently dropped
        assert "Synthetic advisory text." in text
        assert "DRNAREA" in text and "100" in text
        assert "| 0.01 |" in text
        assert "none: area-weighted across 2 regions" in text

    def test_to_markdown_single_region_shows_interval(self):
        est = _estimate([_region("GC1", 100.0)])
        text = est.to_markdown()
        lower = est.quantiles.set_index("aep").loc[0.01, "lower"]
        assert f"{lower:,.0f}" in text
        assert "1.000" in text  # the single region's fraction

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
        assert est.weighted is None and est.provenance == {}
        assert "XX" in est.to_markdown()
