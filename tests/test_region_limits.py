"""Default region rules for overlapping sets: stated limits, then preferences.

A limit (:class:`~flowfreq.regression.library.RegionLimit`) takes a set out of the running where
it cannot apply, so a truly urban or regulated basin still goes to the caller. A preference orders
sets that all apply. ``region=`` overrides both.
"""

import pytest

from flowfreq.regression.equations import OutOfRangeError
from flowfreq.regression.library import RegionLimit, load_state, parse_state
from flowfreq.streamstats import RegressionRegionLocation
from flowfreq.ungaged import OverlappingRegionsError, assemble, region_for_nss_code

KW = dict(aeps=[0.01], allow_extrapolation=False, weighting_space="log")


def _located(*weights):
    return [
        RegressionRegionLocation(
            region_id=i, code=code, name=code, percent_weight=weight, area_sq_mi=100.0
        )
        for i, (code, weight) in enumerate(weights)
    ]


def _used(assembly):
    return [r.code for r in assembly.used]


class TestRegionLimit:
    def test_bounds_are_inclusive(self):
        limit = RegionLimit("DRNAREA", minimum=0.1, maximum=1.0)
        assert not limit.excludes({"DRNAREA": 1.0})
        assert not limit.excludes({"DRNAREA": 0.1})
        assert limit.excludes({"DRNAREA": 1.01})
        assert limit.excludes({"DRNAREA": 0.09})

    def test_an_unknown_characteristic_excludes_nothing(self):
        assert not RegionLimit("DRNAREA", maximum=1.0).excludes({})
        assert not RegionLimit("DRNAREA", maximum=1.0).excludes({"DRNAREA": float("nan")})
        assert not RegionLimit("DAUNREG", below="CONTDA").excludes({"DAUNREG": 5.0})

    def test_below_compares_two_characteristics(self):
        limit = RegionLimit("DAUNREG", below="CONTDA")
        assert limit.excludes({"DAUNREG": 52.464, "CONTDA": 52.464})
        assert not limit.excludes({"DAUNREG": 40.0, "CONTDA": 52.464})

    def test_a_limit_needs_a_bound(self):
        with pytest.raises(ValueError, match="DRNAREA"):
            RegionLimit.from_dict({"variable": "DRNAREA"})

    def test_a_limit_for_an_unknown_region_is_refused(self):
        lib = load_state("GA")
        assert set(lib.region_limits) == {"HR1-5", "HR1", "HR3", "HR4", "HR5"}
        import json
        from pathlib import Path

        import flowfreq.regression.library as library

        path = Path(library.__file__).parent.parent / "data" / "regression" / "GA.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        data["region_limits"]["NOPE"] = [{"variable": "DRNAREA", "maximum": 1.0}]
        with pytest.raises(ValueError, match="region_limits"):
            parse_state(data)


class TestSmallBasinSets:
    """GA, NC and SC: the small rural sets stop at 1 mi2, where they cannot be an alternative."""

    CASES = [("GA", "GC1572", "GC1934"), ("NC", "GC1580", "GC1936"), ("SC", "GC1587", "GC1935")]

    @staticmethod
    def _rural(area):
        return {
            "DRNAREA": area,
            "PCTREG1": 100.0,
            "PCTREG2": 0.0,
            "PCTREG3": 0.0,
            "PCTREG4": 0.0,
            "PCTREG5": 0.0,
            "LC06IMP": 1.0,
            "LC19IMP": 1.0,
        }

    @pytest.mark.parametrize("state,small,rural", CASES)
    def test_above_1_mi2_the_rural_set_is_used(self, state, small, rural):
        a = assemble(
            load_state(state),
            _located((small, 100.0), (rural, 100.0)),
            self._rural(131.0),
            **KW,
        )
        assert _used(a) == [rural]
        assert [r.code for r in a.alternatives] == [small]
        assert "region limit" in a.region_choice_source
        assert small in a.region_choice_source

    @pytest.mark.parametrize("state,small,rural", CASES)
    def test_a_small_set_on_part_of_the_basin_no_longer_blocks_it(self, state, small, rural):
        a = assemble(
            load_state(state),
            _located((small, 78.0), (rural, 100.0)),
            self._rural(131.0),
            **KW,
        )
        assert _used(a) == [rural]

    @pytest.mark.parametrize("state,small,rural", CASES)
    def test_at_or_under_1_mi2_a_rural_basin_gets_the_newer_set(self, state, small, rural):
        a = assemble(
            load_state(state),
            _located((small, 100.0), (rural, 100.0)),
            self._rural(0.5),
            **KW,
        )
        assert _used(a) == [rural]
        assert "region preference" in a.region_choice_source

    @pytest.mark.parametrize("state,small,rural", CASES)
    def test_at_or_under_1_mi2_an_urban_basin_gets_the_small_set(self, state, small, rural):
        a = assemble(
            load_state(state),
            _located((small, 100.0), (rural, 100.0)),
            {**self._rural(0.5), "LC19IMP": 25.0, "LC06IMP": 25.0, "LC06DEV": 50.0, "I24H50Y": 8.0},
            **KW,
        )
        assert _used(a) == [small]
        assert [r.code for r in a.alternatives] == [rural]

    @pytest.mark.parametrize("state,small,rural", CASES)
    def test_above_1_mi2_an_urban_basin_is_refused(self, state, small, rural):
        """Nothing else applies, and the rural set's own LC19IMP range refuses it."""
        for located in (_located((small, 100.0), (rural, 100.0)), _located((rural, 100.0))):
            with pytest.raises(OutOfRangeError, match="LC19IMP=25"):
                assemble(
                    load_state(state),
                    located,
                    {**self._rural(131.0), "LC19IMP": 25.0},
                    **KW,
                )

    @pytest.mark.parametrize("state,small,rural", CASES)
    def test_exactly_10_percent_is_still_rural(self, state, small, rural):
        a = assemble(
            load_state(state),
            _located((rural, 100.0)),
            {**self._rural(131.0), "LC19IMP": 10.0},
            **KW,
        )
        assert _used(a) == [rural]

    @pytest.mark.parametrize("state,small,rural", CASES)
    def test_naming_the_small_set_still_reaches_it(self, state, small, rural):
        with pytest.raises(OutOfRangeError, match="DRNAREA"):
            assemble(
                load_state(state),
                _located((small, 100.0), (rural, 100.0)),
                self._rural(131.0),
                region=small,
                **KW,
            )


class TestRuralOverUrban:
    def test_al_rural_basin_gets_its_flood_region(self):
        a = assemble(
            load_state("AL"),
            _located(("GC1499", 100.0), ("GC1503", 100.0), ("GC1504", 100.0)),
            {"CONTDA": 194.0, "URBAN": 5.49},
            **KW,
        )
        assert _used(a) == ["GC1499"]
        assert sorted(r.code for r in a.alternatives) == ["GC1503", "GC1504"]

    def test_al_small_streams_stop_at_the_reports_5_mi2(self):
        a = assemble(
            load_state("AL"),
            _located(("GC1502", 100.0), ("GC1503", 100.0)),
            {"CONTDA": 9.0, "URBAN": 5.0},
            **KW,
        )
        assert _used(a) == ["GC1502"]

    def test_al_under_5_mi2_gets_small_streams(self):
        a = assemble(
            load_state("AL"),
            _located(("GC1502", 100.0), ("GC1503", 100.0)),
            {"CONTDA": 2.0, "URBAN": 5.0},
            **KW,
        )
        assert _used(a) == ["GC1503"]
        assert "region preference" in a.region_choice_source

    def test_al_urban_basin_stays_the_callers(self):
        """Not silently rural: the report says not to use those over 10 percent impervious."""
        with pytest.raises(OverlappingRegionsError):
            assemble(
                load_state("AL"),
                _located(("GC1499", 100.0), ("GC1504", 100.0)),
                {"CONTDA": 30.0, "URBAN": 25.0},
                **KW,
            )

    def test_al_boundary_basin_weights_its_flood_regions(self):
        a = assemble(
            load_state("AL"),
            _located(("GC1499", 60.0), ("GC1500", 40.0), ("GC1503", 100.0), ("GC1504", 100.0)),
            {"CONTDA": 194.0, "URBAN": 5.49},
            **KW,
        )
        assert _used(a) == ["GC1499", "GC1500"]

    def test_mo_rural_basin_gets_its_region(self):
        a = assemble(
            load_state("MO"),
            _located(("GC1483", 100.0), ("GC1486", 100.0)),
            {"DRNAREA": 100.0, "BSHAPE": 10.0, "IMPNLCD01": 0.768},
            **KW,
        )
        assert _used(a) == ["GC1483"]

    def test_mo_urban_basin_stays_the_callers(self):
        with pytest.raises(OverlappingRegionsError):
            assemble(
                load_state("MO"),
                _located(("GC1483", 100.0), ("GC1486", 100.0)),
                {"DRNAREA": 100.0, "BSHAPE": 10.0, "IMPNLCD01": 20.0},
                **KW,
            )

    def test_mo_basin_past_the_urban_range_gets_its_region(self):
        a = assemble(
            load_state("MO"),
            _located(("GC1483", 100.0), ("GC1486", 100.0)),
            {"DRNAREA": 300.0, "BSHAPE": 10.0, "IMPNLCD01": 20.0},
            **KW,
        )
        assert _used(a) == ["GC1483"]

    def test_the_callers_choice_reaches_the_urban_set(self):
        a = assemble(
            load_state("MO"),
            _located(("GC1483", 100.0), ("GC1486", 100.0)),
            {"DRNAREA": 100.0, "BSHAPE": 10.0, "IMPNLCD01": 3.0},
            region="GC1486",
            **KW,
        )
        assert _used(a) == ["GC1486"]


class TestIllinoisAdjusted:
    """The report's procedure outside region 2 is the adjusted one, so it comes first."""

    FULL = {"DRNAREA": 148.0, "STATSSLPWT": 5.0, "URBTHE2010": 0.021}

    def test_adjusted_first(self):
        a = assemble(
            load_state("IL"), _located(("GC2004", 100.0), ("GC2010", 100.0)), self.FULL, **KW
        )
        assert _used(a) == ["GC2010"]
        assert "region preference" in a.region_choice_source

    def test_unadjusted_when_urbthe_is_missing(self):
        values = {k: v for k, v in self.FULL.items() if k != "URBTHE2010"}
        a = assemble(load_state("IL"), _located(("GC2004", 100.0), ("GC2010", 100.0)), values, **KW)
        assert _used(a) == ["GC2004"]

    def test_region_2_has_one_set_and_joins_either_side(self):
        values = {
            **self.FULL,
            "FLC16DVLMH": 0.1,
            "FSSURGDC78": 0.1,
            "RELRELF": 10.0,
        }
        a = assemble(
            load_state("IL"),
            _located(("GC2004", 60.0), ("GC2010", 60.0), ("GC2002", 40.0)),
            values,
            **KW,
        )
        assert _used(a) == ["GC2010", "GC2002"]


class TestOklahomaNrcs:
    def test_no_structures_uses_the_unregulated_set(self):
        a = assemble(
            load_state("OK"),
            _located(("GC1487", 100.0), ("GC1488", 100.0)),
            {"CONTDA": 52.464, "DAUNREG": 52.464, "PRECIP": 30.0},
            **KW,
        )
        assert _used(a) == ["GC1487"]

    def test_structures_stay_the_callers(self):
        """Unregulated would overstate a basin whose structures cut its contributing area."""
        with pytest.raises(OverlappingRegionsError):
            assemble(
                load_state("OK"),
                _located(("GC1487", 100.0), ("GC1488", 100.0)),
                {"CONTDA": 52.464, "DAUNREG": 30.0, "PRECIP": 30.0},
                **KW,
            )


class TestFullSetsFirst:
    @pytest.mark.parametrize(
        "state,full,da_only",
        [
            ("CT", "GC1801", "GC1800"),
            ("ME", "GC1798", "GC1799"),
            ("VT", "GC2067", "GC2068"),
            ("TN", "GC2016", "GC2020"),
        ],
    )
    def test_full_then_drainage_area_only(self, state, full, da_only):
        lib = load_state(state)
        located = _located((full, 100.0), (da_only, 100.0))
        code = region_for_nss_code(lib, full)
        values = {v.code: (v.minimum + v.maximum) / 2 for v in lib.equation(code, 0.01).variables}
        a = assemble(lib, located, values, **KW)
        assert _used(a) == [full]
        assert [r.code for r in a.alternatives] == [da_only]
        assert "method choice" in a.region_choice_source

        rest = {"DRNAREA": values["DRNAREA"]}
        assert _used(assemble(lib, located, rest, **KW)) == [da_only]
        assert _used(assemble(lib, located, values, region=da_only, **KW)) == [da_only]

    def test_tn_urban_small_basin_stays_the_callers(self):
        with pytest.raises(OverlappingRegionsError):
            assemble(
                load_state("TN"),
                _located(("GC2016", 100.0), ("GC2020", 100.0), ("GC351", 100.0)),
                {
                    "DRNAREA": 10.0,
                    "CSL10_85": 20.0,
                    "CONTDA": 10.0,
                    "IMPERV": 30.0,
                    "I24H2Y": 3.5,
                },
                **KW,
            )


def test_ar_prefers_the_alluvial_plain_where_it_covers_the_basin():
    values = {"DRNAREA": 374.0, "CSLBlue": 0.82, "BSHAPELFP": 13.0}
    a = assemble(load_state("AR"), _located(("GC1917", 100.0), ("GC1719", 100.0)), values, **KW)
    assert _used(a) == ["GC1917"]
    assert [r.code for r in a.alternatives] == ["GC1719"]
    a = assemble(
        load_state("AR"),
        _located(("GC1917", 85.0), ("GC1719", 100.0)),
        values,
        **KW,
    )
    assert _used(a) == ["GC1719"]


def test_a_selector_for_an_unknown_region_is_refused():
    import json
    from pathlib import Path

    data = json.loads(
        (Path(__file__).parent.parent / "flowfreq/data/regression/IN.json").read_text("utf-8")
    )
    data["region_selector"]["values"]["1005"] = "NOPE"
    with pytest.raises(ValueError, match="region_selector"):
        parse_state(data)


def test_a_limit_never_leaves_part_of_the_basin_uncovered():
    """Both limits exclude the basin, but dropping them leaves 60 percent: nothing is dropped."""
    with pytest.raises(OverlappingRegionsError):
        assemble(
            load_state("AL"),
            _located(("GC1499", 60.0), ("GC1503", 100.0), ("GC1504", 100.0)),
            {"CONTDA": 194.0, "URBAN": 5.0},
            **KW,
        )


def test_limits_change_nothing_where_sets_do_not_overlap():
    a = assemble(load_state("AL"), _located(("GC1499", 100.0)), {"CONTDA": 194.0}, **KW)
    assert _used(a) == ["GC1499"]
    assert a.region_choice_source is None
