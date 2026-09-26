"""Tests for flowfreq.donor_similarity (issue #13)."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import pytest

from flowfreq.donor_similarity import (
    rank_donors_by_similarity,
    screen_donors,
    similarity_donor_chooser,
)
from flowfreq.qppq import loocv_qppq
from flowfreq.streamstats import Characteristic
from flowfreq.transpose import DEFAULT_AREA_RATIO_RANGE


def _char(code: str, value: float, unit: str) -> Characteristic:
    return Characteristic(code=code, name=code, description=code, value=value, unit=unit)


class TestScreenDonors:
    def test_area_only_matches_the_transpose_band(self):
        """Without predictors this is exactly the existing area screen."""
        low, high = DEFAULT_AREA_RATIO_RANGE
        table = screen_donors(100.0, {"in": 100.0 / ((low + high) / 2), "out": 100.0 / (high * 2)})
        by = table.set_index("donor")
        assert by.loc["in", "status"] == "pass"
        assert bool(by.loc["in", "area_ok"])
        assert by.loc["out", "status"] == "fail"
        assert "area ratio" in by.loc["out", "reason"]
        assert list(table["donor"]) == ["in", "out"]
        assert not any(c.startswith("ratio_") for c in table.columns)

    def test_area_ratio_is_target_over_donor(self):
        table = screen_donors(252.0, {"goat": 412.0})
        assert table.loc[0, "area_ratio"] == pytest.approx(252.0 / 412.0)

    def test_precipitation_mismatch_fails_a_donor_the_area_screen_passes(self):
        """Issue #13's point: an area ratio inside the band (252/412, its Goat
        Ck -> Lost R pair) says nothing about precipitation. Predictor values
        are the issue's own suggested-shape example."""
        table = screen_donors(
            252.0,
            {"goat": 412.0},
            target_predictors={"PRECPRIS10": 45.6},
            donor_predictors={"goat": {"PRECPRIS10": 32.2}},
            similarity_band={"PRECPRIS10": (0.85, 1.18)},
        )
        row = table.iloc[0]
        assert row["area_ok"]
        assert not row["similarity_ok"]
        assert row["status"] == "fail"
        assert row["ratio_PRECPRIS10"] == pytest.approx(45.6 / 32.2)
        assert "PRECPRIS10" in row["reason"] and "area" not in row["reason"]

    def test_passing_donor_sorts_first(self):
        table = screen_donors(
            100.0,
            {"wet": 100.0, "dry": 100.0},
            target_predictors={"PRECPRIS10": 40.0},
            donor_predictors={"wet": {"PRECPRIS10": 39.0}, "dry": {"PRECPRIS10": 20.0}},
            similarity_band={"PRECPRIS10": (0.85, 1.18)},
        )
        assert list(table["donor"]) == ["wet", "dry"]
        assert list(table["status"]) == ["pass", "fail"]

    def test_missing_donor_predictor_is_excluded_not_imputed(self):
        table = screen_donors(
            100.0,
            {"a": 100.0, "b": 100.0},
            target_predictors={"PRECPRIS10": 40.0},
            donor_predictors={"a": {"PRECPRIS10": 40.0}},
            similarity_band={"PRECPRIS10": (0.85, 1.18)},
        )
        b = table.set_index("donor").loc["b"]
        assert b["status"] == "excluded"
        assert "PRECPRIS10 missing" in b["reason"]
        assert np.isnan(b["ratio_PRECPRIS10"])
        assert b["similarity_ok"] is None

    def test_unit_mismatch_is_excluded(self):
        table = screen_donors(
            100.0,
            {"a": 100.0},
            target_predictors={"PRECPRIS10": _char("PRECPRIS10", 40.0, "inches")},
            donor_predictors={"a": {"PRECPRIS10": _char("PRECPRIS10", 1000.0, "millimeters")}},
            similarity_band={"PRECPRIS10": (0.85, 1.18)},
        )
        assert table.loc[0, "status"] == "excluded"
        assert "millimeters" in table.loc[0, "reason"]

    def test_target_missing_predictor_raises(self):
        with pytest.raises(ValueError, match="will not be imputed"):
            screen_donors(
                100.0,
                {"a": 100.0},
                target_predictors={},
                donor_predictors={"a": {"PRECPRIS10": 40.0}},
                similarity_band={"PRECPRIS10": (0.85, 1.18)},
            )

    def test_band_without_target_predictors_raises(self):
        with pytest.raises(ValueError, match="target_predictors"):
            screen_donors(100.0, {"a": 100.0}, similarity_band={"PRECPRIS10": (0.85, 1.18)})

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"target_area_sqmi": 0.0, "donor_areas": {"a": 1.0}},
            {"target_area_sqmi": 1.0, "donor_areas": {}},
            {"target_area_sqmi": 1.0, "donor_areas": {"a": -1.0}},
            {"target_area_sqmi": 1.0, "donor_areas": {"a": 1.0}, "area_ratio_range": (2.0, 1.0)},
        ],
    )
    def test_bad_inputs_raise(self, kwargs):
        with pytest.raises(ValueError):
            screen_donors(**kwargs)


BASINS = {
    "twin": {"DRNAREA": 100.0, "PRECPRIS10": 40.0, "ELEV": 5000.0},
    "near": {"DRNAREA": 120.0, "PRECPRIS10": 38.0, "ELEV": 4800.0},
    "far": {"DRNAREA": 900.0, "PRECPRIS10": 20.0, "ELEV": 2500.0},
    "wet_big": {"DRNAREA": 400.0, "PRECPRIS10": 40.0, "ELEV": 5000.0},
}
TARGET = {"DRNAREA": 100.0, "PRECPRIS10": 40.0, "ELEV": 5000.0}
CODES = ["DRNAREA", "PRECPRIS10", "ELEV"]


class TestRankDonorsBySimilarity:
    def test_identical_basin_ranks_first_at_distance_zero(self):
        table = rank_donors_by_similarity(TARGET, BASINS, CODES)
        assert table.loc[0, "donor"] == "twin"
        assert table.loc[0, "distance"] == 0.0
        assert table.loc[0, "rank"] == 1
        assert table.loc[len(table) - 1, "donor"] == "far"
        assert np.all(np.diff(table["distance"].to_numpy()) >= 0)

    def test_distance_is_burn_standardized_euclidean(self):
        scale = {"DRNAREA": 0.5, "PRECPRIS10": 10.0, "ELEV": 1000.0}
        table = rank_donors_by_similarity(TARGET, BASINS, CODES, scale=scale).set_index("donor")
        expected = np.sqrt(
            (np.log10(120.0 / 100.0) / 0.5) ** 2 + (2.0 / 10.0) ** 2 + (200.0 / 1000.0) ** 2
        )
        assert table.loc["near", "distance"] == pytest.approx(expected)
        assert table.loc["near", "z_PRECPRIS10"] == pytest.approx(-0.2)

    def test_default_scale_is_pool_sample_std(self):
        table = rank_donors_by_similarity(
            TARGET, BASINS, ["PRECPRIS10"], log_attributes=()
        ).set_index("donor")
        pool = [40.0, 40.0, 38.0, 20.0, 40.0]
        assert table.loc["far", "z_PRECPRIS10"] == pytest.approx(-20.0 / np.std(pool, ddof=1))

    def test_weights_change_the_ranking(self):
        """With equal weights the area-matched donor wins; weighting only
        precipitation and elevation makes the larger but climatically
        identical basin tie with the twin and beat the near one."""
        donors = {k: BASINS[k] for k in ("near", "wet_big", "far")}
        equal = rank_donors_by_similarity(TARGET, donors, CODES)
        assert equal.loc[0, "donor"] == "near"
        climate = rank_donors_by_similarity(
            TARGET, donors, CODES, weights={"DRNAREA": 0.0, "PRECPRIS10": 1.0, "ELEV": 1.0}
        )
        assert climate.loc[0, "donor"] == "wet_big"
        assert climate.loc[0, "distance"] == 0.0

    def test_missing_donor_attribute_is_excluded_and_reported(self):
        donors = dict(BASINS, gap={"DRNAREA": 100.0, "PRECPRIS10": 40.0})
        table = rank_donors_by_similarity(TARGET, donors, CODES).set_index("donor")
        assert table.loc["gap", "status"] == "excluded"
        assert "ELEV missing" in table.loc["gap", "reason"]
        assert np.isnan(table.loc["gap", "distance"])
        assert np.isnan(table.loc["gap", "rank"])
        assert table.index[-1] == "gap"

    def test_nonpositive_logged_attribute_is_excluded(self):
        donors = dict(BASINS, zero={"DRNAREA": 0.0, "PRECPRIS10": 40.0, "ELEV": 5000.0})
        table = rank_donors_by_similarity(TARGET, donors, CODES).set_index("donor")
        assert table.loc["zero", "status"] == "excluded"
        assert "log10" in table.loc["zero", "reason"]

    def test_unit_mismatch_is_excluded(self):
        target = {"PRECPRIS10": _char("PRECPRIS10", 40.0, "inches")}
        donors = {
            "a": {"PRECPRIS10": _char("PRECPRIS10", 41.0, "inches")},
            "b": {"PRECPRIS10": _char("PRECPRIS10", 30.0, "inches")},
            "mm": {"PRECPRIS10": _char("PRECPRIS10", 1000.0, "millimeters")},
        }
        table = rank_donors_by_similarity(target, donors, ["PRECPRIS10"]).set_index("donor")
        assert table.loc["mm", "status"] == "excluded"
        assert table.loc["a", "rank"] == 1

    def test_target_missing_attribute_raises(self):
        with pytest.raises(ValueError, match="will not be imputed"):
            rank_donors_by_similarity({"DRNAREA": 100.0}, BASINS, CODES)

    def test_too_small_pool_without_scale_raises(self):
        with pytest.raises(ValueError, match="Supply scale"):
            rank_donors_by_similarity(TARGET, {"twin": BASINS["twin"]}, CODES)
        table = rank_donors_by_similarity(
            TARGET, {"near": BASINS["near"]}, CODES, scale={c: 1.0 for c in CODES}
        )
        assert table.loc[0, "rank"] == 1

    @pytest.mark.parametrize(
        "kwargs, match",
        [
            ({"attributes": []}, "empty"),
            ({"attributes": ["ELEV", "ELEV"]}, "duplicates"),
            ({"weights": {"DRNAREA": 1.0}}, "exactly"),
            ({"weights": {"DRNAREA": -1.0, "PRECPRIS10": 1.0, "ELEV": 1.0}}, "non-negative"),
            ({"weights": {c: 0.0 for c in CODES}}, "positive"),
            ({"scale": {"DRNAREA": 1.0}}, "missing"),
            ({"scale": {c: 0.0 for c in CODES}}, "positive"),
            ({"donor_attributes": {}}, "no donors"),
        ],
    )
    def test_bad_inputs_raise(self, kwargs, match):
        args = {"target_attributes": TARGET, "donor_attributes": BASINS, "attributes": CODES}
        args.update(kwargs)
        with pytest.raises(ValueError, match=match):
            rank_donors_by_similarity(**args)

    def test_exclusion_is_logged(self, caplog):
        donors = dict(BASINS, gap={})
        with caplog.at_level(logging.INFO, logger="flowfreq.donor_similarity"):
            rank_donors_by_similarity(TARGET, donors, CODES)
        assert "gap" in caplog.text


def _snowmelt(seed: int, scale: float, phase_days: int, years: int = 6) -> pd.DataFrame:
    index = pd.date_range("2000-10-01", periods=365 * years, freq="D")
    doy = index.dayofyear.to_numpy() + phase_days
    melt = np.exp(-(((doy - 150) % 365) ** 2) / (2 * 45.0**2))
    rng = np.random.default_rng(seed)
    flows = scale * (0.15 + 3.0 * melt) * np.exp(rng.normal(0, 0.25, len(index)))
    return pd.DataFrame({"flow_cfs": flows}, index=index)


class TestSimilarityDonorChooser:
    def test_plugs_into_loocv_and_picks_the_most_similar_site(self):
        """Plumbing check, not evidence of skill: two elevation bands with
        their own melt timing; the chooser must pair each site with its band."""
        sites = {
            "low1": _snowmelt(1, 100.0, 0),
            "low2": _snowmelt(2, 80.0, 0),
            "high1": _snowmelt(3, 60.0, -30),
            "high2": _snowmelt(4, 50.0, -30),
        }
        attrs = {
            "low1": {"DRNAREA": 100.0, "ELEV": 2000.0},
            "low2": {"DRNAREA": 80.0, "ELEV": 2100.0},
            "high1": {"DRNAREA": 90.0, "ELEV": 6000.0},
            "high2": {"DRNAREA": 70.0, "ELEV": 6100.0},
        }
        chooser = similarity_donor_chooser(attrs, ["DRNAREA", "ELEV"])
        table = loocv_qppq(sites, seasonal=False, donor_chooser=chooser).set_index("site")
        assert table.loc["low1", "donor"] == "low2"
        assert table.loc["high2", "donor"] == "high1"

    def test_missing_withheld_site_raises(self):
        chooser = similarity_donor_chooser({}, ["DRNAREA"])
        with pytest.raises(KeyError, match="no basin attributes"):
            chooser("x", pd.DataFrame(), {"y": pd.DataFrame()})

    def test_no_rankable_donor_raises(self):
        chooser = similarity_donor_chooser(
            {"x": {"DRNAREA": 1.0}}, ["DRNAREA"], scale={"DRNAREA": 1.0}
        )
        with pytest.raises(ValueError, match="no donor"):
            chooser("x", pd.DataFrame(), {"y": pd.DataFrame()})
