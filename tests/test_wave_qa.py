"""Tests for flowfreq.validation.wave_qa (roadmap §3.2 wave definition of done).

Offline: synthetic geometry, the packaged state libraries, and trimmed inputs for one
real border basin captured by ``tools/wave_qa.py`` (``tests/fixtures/wave_qa``). The
one ``requires_network`` test asks NSS for Washington's regions around a point.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from flowfreq.regression import evaluate, load_state
from flowfreq.validation import wave_qa as qa

FIXTURES = Path(__file__).parent / "fixtures" / "wave_qa"

SQUARE = {"type": "Polygon", "coordinates": [[[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]]}


@pytest.fixture(scope="module")
def libs():
    return {st: load_state(st) for st in qa.WAVE1 + qa.WAVE2}


class TestGeometry:
    def test_distance_inside_is_zero(self):
        assert qa.distance_km(0.5, 0.5, SQUARE) == 0.0

    def test_distance_outside(self):
        # One tenth of a degree east of the square at the equator: ~11.1 km.
        assert qa.distance_km(0.5, 1.1, SQUARE) == pytest.approx(11.13, rel=0.01)

    def test_feature_accepted(self):
        assert qa.point_in(0.5, 0.5, {"type": "Feature", "geometry": SQUARE})

    def test_feature_without_geometry_raises(self):
        with pytest.raises(ValueError):
            qa.point_in(0.5, 0.5, {"type": "Feature"})

    def test_circle_is_closed_with_radius(self):
        c = qa.circle(45.0, -110.0, 10.0, n=32)
        ring = c["coordinates"][0]
        assert ring[0] == ring[-1] and len(ring) == 33
        # area of a 10 km circle is 314 km2 = 121 mi2
        assert qa.area_sq_mi(c) == pytest.approx(121.3, rel=0.02)

    def test_circle_rejects_bad_radius(self):
        with pytest.raises(ValueError):
            qa.circle(45.0, -110.0, 0.0)

    def test_sample_points_fraction(self):
        pts = qa.sample_points(SQUARE, 400)
        assert len(pts) >= 240
        west_half = {
            "type": "Polygon",
            "coordinates": [[[0, 0], [0.5, 0], [0.5, 1], [0, 1], [0, 0]]],
        }
        assert qa.fraction_inside(pts, west_half) == pytest.approx(0.5, abs=0.05)

    def test_sample_points_rejects_small_n(self):
        with pytest.raises(ValueError):
            qa.sample_points(SQUARE, 5)

    def test_fraction_inside_needs_points(self):
        with pytest.raises(ValueError):
            qa.fraction_inside([], SQUARE)

    def test_esri_rings_keep_holes(self):
        outer = [[0, 0], [0, 4], [4, 4], [4, 0], [0, 0]]  # clockwise: exterior
        hole = [[1, 1], [2, 1], [2, 2], [1, 2], [1, 1]]  # counter-clockwise: hole
        g = qa.esri_rings_to_geojson([outer, hole])
        assert g["type"] == "MultiPolygon" and len(g["coordinates"]) == 1
        assert len(g["coordinates"][0]) == 2
        assert not qa.point_in(1.5, 1.5, g)
        assert qa.point_in(3.0, 3.0, g)


class TestRegionMapping:
    def test_washington_codes(self, libs):
        m = qa.nss_region_map(libs["WA"])
        assert m == {"GC1750": "1", "GC1751": "2", "GC1752": "3", "GC1753": "4"}

    def test_oregon_eastern_regions_added(self, libs):
        m = qa.nss_region_map(libs["OR"])
        assert m["GC1754"] == "E1" and m["GC1759"] == "E6" and m["GC729"] == "1"

    def test_navajo_and_channel_width_excluded(self, libs):
        assert not any(r.startswith("Navajo") for r in qa.nss_region_map(libs["CO"]).values())
        assert "W-AC" not in qa.nss_region_map(libs["MT"]).values()
        assert qa.excluded_reason("NWF-BF") and qa.excluded_reason("Navajo8")
        assert qa.excluded_reason("W") is None

    def test_weights_ignore_other_codes_and_normalise(self, libs):
        w, cov = qa.region_weights(libs["WA"], [("GC1750", 60), ("GC1752", 30), ("GC1294", 100)])
        assert w == pytest.approx({"1": 2 / 3, "3": 1 / 3})
        assert cov == pytest.approx(0.9)

    def test_oregon_2a_2b_collapse(self, libs):
        w, _ = qa.region_weights(libs["OR"], [("GC730", 100), ("GC731", 100)])
        assert w == {"2": 1.0}

    def test_no_region(self, libs):
        assert qa.region_weights(libs["WA"], [("GC1", 100)]) == ({}, 0.0)

    def test_arizona_high_elevation_overlay(self):
        codes = ["gc1619", "gc1618"]
        assert qa.resolve_overlay("AZ", codes, {"ELEV": 8000.0}) == ["GC1618"]
        assert qa.resolve_overlay("AZ", codes, {"ELEV": 6000.0}) == ["GC1619"]
        assert qa.resolve_overlay("AZ", ["GC1620"], {}) == ["GC1620"]
        with pytest.raises(KeyError):
            qa.resolve_overlay("AZ", codes, {})

    def test_oregon_western_polygon_is_region_2(self):
        got = qa.resolve_overlay("OR", ["gc10001", "gc730", "gc731", "gc729"], {})
        assert got == ["GC730", "GC731"]
        assert qa.resolve_overlay("OR", ["gc1756"], {}) == ["GC1756"]

    def test_highreg(self):
        assert qa.highreg_region(1096.0) == "GC1096"
        with pytest.raises(ValueError):
            qa.highreg_region(1096.5)


class TestEstimates:
    def test_prepare_characteristics(self):
        v, subs = qa.prepare_characteristics({"DRNAREA": 50.0}, 44.5, -108.25)
        assert v["LONG_OUT"] == 108.25 and v["LAT_OUT"] == 44.5 and v["CONTDA"] == 50.0
        assert subs == ("CONTDA=DRNAREA",)
        v, subs = qa.prepare_characteristics({"CONTDA": 50.0, "LAT_OUT": 1.0}, 44.5, -108.25)
        assert v["DRNAREA"] == 50.0 and v["LAT_OUT"] == 1.0 and subs == ("DRNAREA=CONTDA",)

    def test_vintage_proxies_prefer_first_source(self):
        v, subs = qa.prepare_characteristics(
            {"DRNAREA": 5.0, "CONTDA": 5.0, "PRECIP": 20.0, "PRECPRIS20": 21.0}, 40.0, -110.0
        )
        assert v["PRECPRIS10"] == 21.0 and "PRECPRIS10=PRECPRIS20" in subs
        assert v["PRECIP"] == 20.0  # present: never replaced

    def test_estimate_reports_only_proxies_it_used(self, libs):
        chars, subs = qa.prepare_characteristics(
            {"DRNAREA": 100.0, "PRECIP": 30.0, "LC16FOREST": 40.0}, 47.0, -116.0
        )
        est = qa.estimate_state(libs["ID"], {"5": 1.0}, chars, 0.01, substitutions=subs)
        assert est.substitutions == ("PRECPRIS10=PRECIP",)

    def test_single_region_matches_evaluate(self, libs):
        chars = {"DRNAREA": 79.22, "PRECPRIS10": 67.6, "CANOPY_PCT": 68.1}
        est = qa.estimate_state(libs["WA"], {"2": 1.0}, chars, 0.01)
        ref = evaluate(libs["WA"].equation("2", 0.01), chars)
        assert est.flow_cfs == pytest.approx(ref.flow_cfs, rel=1e-12)
        assert est.sep_log == libs["WA"].equation("2", 0.01).sep_log
        assert est.log_flow == pytest.approx(math.log10(ref.flow_cfs))

    def test_two_regions_weighted_linearly(self, libs):
        chars = {"DRNAREA": 100.0, "PRECPRIS10": 40.0, "CANOPY_PCT": 50.0}
        a = evaluate(libs["WA"].equation("1", 0.1), chars).flow_cfs
        b = evaluate(libs["WA"].equation("3", 0.1), chars).flow_cfs
        est = qa.estimate_state(libs["WA"], {"1": 0.25, "3": 0.75}, chars, 0.1)
        assert est.flow_cfs == pytest.approx(0.25 * a + 0.75 * b)

    def test_out_of_range_is_recorded_not_raised(self, libs):
        est = qa.estimate_state(libs["ID"], {"3": 1.0}, {"DRNAREA": 99999.0}, 0.01)
        assert est.out_of_range and "DRNAREA" in est.out_of_range[0]

    def test_oregon_region_2_blend(self, libs):
        chars = {"DRNAREA": 42.1, "BSLOPD": 23.8, "I24H2Y": 2.84, "JANMINT2K": 31.0}
        chars.update({"JANMAXT2K": 44.3, "ELEV": 2970.0})
        est = qa.estimate_state(libs["OR"], {"2": 1.0}, chars, 0.01)
        # Report pp. 45-46: the eq. 8 blend lies between the 2B and 2A estimates.
        assert 7700 < est.flow_cfs < 8430

    def test_missing_sep_propagates_none(self, libs):
        est = qa.estimate_state(libs["NV"], {"10": 1.0}, {"DRNAREA": 57.0}, 0.01)
        assert est.sep_log is None

    def test_try_estimate_reasons(self, libs):
        est, why = qa.try_estimate_state(libs["WA"], {"1": 1.0}, {"DRNAREA": 10.0}, 0.01)
        assert est is None and why.startswith("missing characteristic")
        est, why = qa.try_estimate_state(libs["WA"], {}, {"DRNAREA": 10.0}, 0.01)
        assert est is None and why.startswith("no equation")
        est, why = qa.try_estimate_state(libs["NV"], {"10": 1.0}, {"DRNAREA": 10.0}, 0.002)
        assert est is None and why.startswith("no equation")


class TestDiscrepancy:
    def test_threshold_and_flag(self):
        d = qa.discrepancy(math.log10(200), 0.2, math.log10(100), 0.2)
        assert d.threshold_log == pytest.approx(qa.Z90 * math.sqrt(0.08))
        assert not d.flagged and d.ratio == pytest.approx(2.0)
        assert qa.discrepancy(math.log10(400), 0.1, math.log10(100), 0.1).flagged

    def test_one_or_no_se(self):
        d = qa.discrepancy(1.0, None, 0.0, 0.3)
        assert d.threshold_log == pytest.approx(qa.Z90 * 0.3) and d.note == "one SE unknown"
        d = qa.discrepancy(1.0, None, 0.0, None)
        assert d.threshold_log is None and not d.flagged and d.note == "no SE"

    def test_rejects_bad_inputs(self):
        with pytest.raises(ValueError):
            qa.discrepancy(1.0, 0.1, 0.0, 0.1, z=0)
        with pytest.raises(ValueError):
            qa.discrepancy(1.0, -0.1, 0.0, 0.1)

    def test_b17c_degenerate(self):
        assert not qa.b17c_degenerate(500.0, 5000.0)
        assert qa.b17c_degenerate(1e-11, 1e15)  # mostly zero-flow short record
        assert qa.b17c_degenerate(2.0, 5000.0)  # ratio above 1,000
        assert qa.b17c_degenerate(float("nan"), 10.0)

    def test_b17c_sd(self):
        assert qa.b17c_log_sd(100.0, 100.0 * 10 ** (2 * qa.Z90 * 0.05)) == pytest.approx(0.05)
        with pytest.raises(ValueError):
            qa.b17c_log_sd(200.0, 100.0)


class TestSummaries:
    def test_bias_and_rmse(self):
        df = pd.DataFrame(
            {"state": ["A", "A", "B"], "log_est": [1.1, 0.9, 2.0], "log_obs": [1.0, 1.0, 1.7]}
        )
        s = qa.summarize_log_errors(df, ["state"]).set_index("state")
        assert s.loc["A", "n"] == 2 and s.loc["A", "bias_log"] == pytest.approx(0.0)
        assert s.loc["A", "rmse_log"] == pytest.approx(0.1)
        assert s.loc["B", "bias_pct"] == pytest.approx((10**0.3 - 1) * 100)

    def test_missing_column(self):
        with pytest.raises(KeyError):
            qa.summarize_log_errors(pd.DataFrame({"x": [1]}), ["state"])

    def _network(self):
        # Q100 = 100 * A**0.7 exactly, so the DAR transfer with b = 0.7 is exact.
        areas = [100.0, 120.0, 90.0, 1000.0]
        return pd.DataFrame(
            {
                "site_no": ["1", "2", "3", "4"],
                "state": ["WA"] * 4,
                "huc8": ["17020001", "17020001", "17020002", "17020001"],
                "latitude": [47.0, 47.1, 47.0, 47.05],
                "longitude": [-120.0, -120.0, -120.1, -120.0],
                "drainage_area_sqmi": areas,
                "log_q_0.01": [2 + 0.7 * math.log10(a) for a in areas],
            }
        )

    def test_dar_loocv_exact_power_law(self):
        net = self._network()
        exps = {(s, 0.01): 0.7 for s in net["site_no"]}
        out = qa.dar_loocv(net, exps, [0.01])
        # Sites 1 and 2 are each other's only eligible donor (same HUC8, ratio 0.83/1.2);
        # 3 is alone in its HUC8 and 4 is out of the area band.
        assert sorted(out["site_no"]) == ["1", "2"]
        assert np.allclose(out["log_est"], out["log_obs"])
        assert set(out["donor"]) == {"1", "2"}

    def test_dar_loocv_skips_targets_without_exponent(self):
        out = qa.dar_loocv(self._network(), {}, [0.01])
        assert out.empty and "log_est" in out.columns

    def test_area_exponent(self, libs):
        eq = libs["WA"].equation("1", 0.01)
        assert qa.area_exponent(libs["WA"], "1", 0.01) == eq.coefficients[0]
        assert qa.area_exponent(libs["MT"], "W", 0.01) is not None  # CONTDA
        assert qa.area_exponent(libs["MT"], "W-AC", 0.01) is None  # width only
        assert qa.area_exponent(libs["NV"], "10", 0.002) is None  # not published


class TestCompareBasin:
    """Trimmed inputs for one border basin captured by tools/wave_qa.py."""

    @pytest.fixture(scope="class")
    def case(self):
        return json.loads((FIXTURES / "border_basin.json").read_text(encoding="utf-8"))

    def test_rows_and_flags_are_consistent(self, libs, case):
        b17c = {float(k): tuple(v) for k, v in case["b17c"].items()}
        rows = qa.compare_basin(
            libs,
            case["home"],
            case["neighbor"],
            case["characteristics"],
            tuple(case["outlet"]),
            [tuple(x) for x in case["located_home"]],
            [tuple(x) for x in case["located_neighbor"]],
            b17c,
            case["aeps"],
        )
        assert [r["aep"] for r in rows] == case["aeps"]
        for r, want in zip(rows, case["expected"]):
            assert r["home_q"] == pytest.approx(want["home_q"], rel=1e-9)
            assert r["neighbor_q"] == pytest.approx(want["neighbor_q"], rel=1e-9)
            assert r["flag_hn"] == want["flag_hn"]
            assert r["flag_hn"] == (abs(r["diff_log"]) > r["thr_log"])
            assert r["diff_log"] == pytest.approx(
                math.log10(r["home_q"]) - math.log10(r["neighbor_q"])
            )

    def test_unevaluable_neighbor_is_reported(self, libs, case):
        rows = qa.compare_basin(
            libs,
            case["home"],
            case["neighbor"],
            case["characteristics"],
            tuple(case["outlet"]),
            [tuple(x) for x in case["located_home"]],
            [],
            None,
            [0.01],
        )
        assert rows[0]["neighbor_q"] is None and rows[0]["neighbor_why"].startswith("no equation")
        assert "flag_hn" not in rows[0] and "b17c_q" not in rows[0]


@pytest.mark.requires_network
def test_nss_locates_washington_region_live(libs):
    """NSS bylocation on a 2 km circle in the Yakima basin falls in a WA peak region."""
    from flowfreq.streamstats import locate_regression_regions

    located = locate_regression_regions("WA", qa.circle(46.75, -120.75, 2.0, 16))
    weights, coverage = qa.region_weights(libs["WA"], [(r.code, r.percent_weight) for r in located])
    assert weights and coverage > 0.9
    assert set(weights) <= {"1", "2", "3", "4"}
