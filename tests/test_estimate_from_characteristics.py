"""``estimate_from_characteristics``: the regression without StreamStats (FL, MI, NE)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from flowfreq import estimate_from_characteristics
from flowfreq.cli import cli
from flowfreq.regression import evaluate
from flowfreq.regression.equations import OutOfRangeError
from flowfreq.regression.library import EquationsUnavailable, load_state

FIX = Path(__file__).resolve().parent / "fixtures" / "regression"


def _fl_gage():
    g = json.loads((FIX / "florida_gages.json").read_text("utf-8"))["rural"]["02327100"]
    return g["region"], {k: float(v) for k, v in g.items() if k not in ("region", "flows")}


def test_florida_matches_the_equation_and_needs_no_streamstats(monkeypatch):
    def no_network(*a, **k):
        raise AssertionError("StreamStats must not be called")

    monkeypatch.setattr("flowfreq.ungaged.delineate_and_get_characteristics", no_network)
    monkeypatch.setattr("flowfreq.ungaged.locate_regression_regions", no_network)
    region, chars = _fl_gage()
    est = estimate_from_characteristics(
        "FL", chars, region=region, lat=30.4, lon=-84.3, aeps=[0.01]
    )
    expected = evaluate(load_state("FL").equation(region, 0.01), chars).flow_cfs
    assert est.quantiles.loc[0, "flow_cfs"] == pytest.approx(expected)
    assert est.basin_geojson is None
    assert est.provenance["characteristics_source"] == "caller"
    assert est.provenance["regions_given"] == {region: 100.0}
    assert est.library_status == "verified"


def test_michigan_and_nebraska_rows_reproduce():
    mi = json.loads((FIX / "michigan_gages.json").read_text("utf-8"))["04001000"]
    chars = {k: float(v) for k, v in mi.items() if k not in ("region", "flows")}
    est = estimate_from_characteristics("MI", chars, region=mi["region"], lat=46.9, lon=-89.6)
    assert len(est.quantiles) >= 1
    row = json.loads((FIX / "nebraska_gages.json").read_text("utf-8"))["rows"][0]
    chars = {k: float(str(v).replace(",", "")) for k, v in row["characteristics"].items()}
    est = estimate_from_characteristics(
        "NE", chars, region=row["region"], lat=42.0, lon=-102.0, aeps=[0.01]
    )
    expected = evaluate(load_state("NE").equation(row["region"], 0.01), chars).flow_cfs
    assert est.quantiles.loc[0, "flow_cfs"] == pytest.approx(expected)


def test_a_basin_spanning_regions_is_weighted():
    _, chars = _fl_gage()
    one = estimate_from_characteristics("FL", chars, region="1", lat=30.4, lon=-84.3, aeps=[0.01])
    two = estimate_from_characteristics("FL", chars, region="2", lat=30.4, lon=-84.3, aeps=[0.01])
    both = estimate_from_characteristics(
        "FL", chars, region={"1": 60.0, "2": 40.0}, lat=30.4, lon=-84.3, aeps=[0.01]
    )
    q = [e.quantiles.loc[0, "flow_cfs"] for e in (one, both, two)]
    assert min(q[0], q[2]) < q[1] < max(q[0], q[2])
    assert len(both.regions) == 2


def test_refusals_say_why():
    _, chars = _fl_gage()
    with pytest.raises(EquationsUnavailable, match="choose from"):
        estimate_from_characteristics("FL", chars, region="9", lat=30.4, lon=-84.3)
    with pytest.raises(ValueError, match="sum to 100"):
        estimate_from_characteristics("FL", chars, region={"1": 60.0}, lat=30.4, lon=-84.3)
    with pytest.raises(OutOfRangeError):
        estimate_from_characteristics(
            "FL", {**chars, "DRNAREA": 1e7}, region="1", lat=30.4, lon=-84.3
        )
    with pytest.raises(EquationsUnavailable):
        estimate_from_characteristics("FL", {}, region="1", lat=30.4, lon=-84.3)


def test_cli_ungaged_chars():
    region, chars = _fl_gage()
    args = ["ungaged-chars", "--state", "FL", "--region", region, "--lat", "30.4", "--lon", "-84.3"]
    for k, v in chars.items():
        args += ["--char", f"{k}={v}"]
    out = CliRunner().invoke(cli, args)
    assert out.exit_code == 0, out.output
    assert "# Regression estimate at 30.4000, -84.3000 (FL)" in out.output
    bad = CliRunner().invoke(cli, args[:-2] + ["--char", "DRNAREA"])
    assert bad.exit_code != 0 and "CODE=VALUE" in bad.output
