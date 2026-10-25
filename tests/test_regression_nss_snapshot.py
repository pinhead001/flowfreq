"""Equation sets harvested from the USGS NSS web service: the ``nss_snapshot`` tier.

Status ``nss_snapshot`` lives under ``nss_snapshot/`` and, like the other opt-in sets,
is returned only on opt-in (``allow_nss_snapshot=True``), so it never changes the
default ``load_state``. Its advisory says the equations are not checked against the
source report.

``tests/fixtures/nss_snapshot_TN_trimmed.json`` is the 2026-10-25 Tennessee snapshot
trimmed to two regions: GC2015 with three equations, and GC341 marked ``unresolved``.
GC2015's 0.2% AEP equation string is synthetic -- a region-indicator interaction, the
GA/NC/SC form the parser refuses -- so that the fixture carries one refusal; its
``parse_error`` is the parser's message for that string.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import logging
import math
from pathlib import Path
from typing import Any, Dict

import pytest

from flowfreq.regression import (
    EquationsUnavailable,
    asep_to_sep_log,
    available_states,
    awaiting_nss_states,
    load_state,
    nss_snapshot_states,
    provisional_states,
    snapshot_to_library,
)
from flowfreq.regression.equations import RegressionEquation
from flowfreq.regression.library import (
    AWAITING_NSS_SUBDIR,
    NSS_SNAPSHOT_SUBDIR,
    OPT_IN_SUBDIRS,
    PROVISIONAL_SUBDIR,
    STATUSES,
    parse_state,
)

FIXTURE = Path(__file__).parent / "fixtures" / "nss_snapshot_TN_trimmed.json"

EQ = {
    "region_code": "GC2015",
    "aep": 0.01,
    "intercept": 2.0,
    "variables": [{"code": "DRNAREA"}],
    "coefficients": [0.7],
    "citation": {"publication": "Synthetic", "table": "T1"},
}
UNPARSED = [
    {
        "region_code": "GC2015",
        "region_name": "Region 1",
        "statistic_code": "PK0_2AEP",
        "aep": 0.002,
        "equation": "DRNAREA^(0.7+0.001*X)",
        "reason": "refused",
    },
    {
        "region_code": "GC341",
        "region_name": "Area 1",
        "statistic_code": "",
        "aep": None,
        "equation": "",
        "reason": "unresolved: no results",
    },
]


def _write(path: Path, status: str, intercept: float, state: str = "TN", **extra: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "schema_version": 1,
        "state": state,
        "status": status,
        "source_reports": [{"citation": "Synthetic report"}],
        "equations": [{**EQ, "intercept": intercept}],
        **extra,
    }
    path.write_text(json.dumps(data), encoding="utf-8")


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    _write(
        tmp_path / NSS_SNAPSHOT_SUBDIR / "TN.json",
        "nss_snapshot",
        5.0,
        retrieved_utc="2026-10-25T17:20:40+00:00",
    )
    return tmp_path


@pytest.fixture
def trimmed() -> Dict[str, Any]:
    data: Dict[str, Any] = json.loads(FIXTURE.read_text(encoding="utf-8"))
    return data


# --------------------------------------------------------------------------- loader


def test_status_is_defined():
    assert "nss_snapshot" in STATUSES
    assert OPT_IN_SUBDIRS["nss_snapshot"] == NSS_SNAPSHOT_SUBDIR == "nss_snapshot"


def test_default_lookup_ignores_the_snapshot(data_dir):
    with pytest.raises(EquationsUnavailable, match="no regression library file"):
        load_state("TN", data_dir=data_dir)
    # The other opt-ins do not reach it either.
    with pytest.raises(EquationsUnavailable):
        load_state("TN", data_dir=data_dir, allow_provisional=True, allow_awaiting_nss=True)


def test_opt_in_returns_it_with_advisory(data_dir, caplog):
    with caplog.at_level(logging.WARNING, logger="flowfreq.regression.library"):
        lib = load_state("TN", data_dir=data_dir, allow_nss_snapshot=True)
    assert lib.nss_snapshot and not lib.provisional and not lib.awaiting_nss
    assert lib.equation("GC2015", 0.01).intercept == 5.0
    assert lib.retrieved_utc == "2026-10-25T17:20:40+00:00"
    assert lib.advisory == (
        "NSS SNAPSHOT (not report-verified): the TN equations were harvested from the USGS "
        "NSS web service on 2026-10-25 and parsed by flowfreq. They have not been checked "
        "against the source report and carry no covariance or model-error variance; "
        "intervals use NSS's average standard error of prediction where NSS reports one. "
        "0 published equations or regions could not be parsed and are listed in "
        "StateLibrary.unparsed."
    )
    assert "NSS SNAPSHOT (not report-verified)" in caplog.text


def test_published_set_wins_without_the_flag(data_dir):
    _write(data_dir / "TN.json", "verified", 2.0)
    lib = load_state("TN", data_dir=data_dir)
    assert lib.status == "verified" and not lib.nss_snapshot and lib.advisory == ""
    assert lib.equation("GC2015", 0.01).intercept == 2.0
    assert load_state("TN", data_dir=data_dir, allow_nss_snapshot=True).nss_snapshot


def test_opt_in_without_a_file_returns_the_published_set(tmp_path):
    _write(tmp_path / "TN.json", "verified", 2.0)
    lib = load_state("TN", data_dir=tmp_path, allow_nss_snapshot=True)
    assert lib.status == "verified" and lib.advisory == ""


def test_awaiting_nss_wins_over_snapshot_when_both_opted_in(data_dir):
    _write(data_dir / PROVISIONAL_SUBDIR / "TN.json", "provisional", 3.0)
    _write(data_dir / AWAITING_NSS_SUBDIR / "TN.json", "awaiting_nss", 4.0)
    every = dict(allow_provisional=True, allow_awaiting_nss=True, allow_nss_snapshot=True)
    assert load_state("TN", data_dir=data_dir, **every).awaiting_nss
    (data_dir / AWAITING_NSS_SUBDIR / "TN.json").unlink()
    assert load_state("TN", data_dir=data_dir, **every).provisional
    (data_dir / PROVISIONAL_SUBDIR / "TN.json").unlink()
    assert load_state("TN", data_dir=data_dir, **every).nss_snapshot


def test_snapshot_status_cannot_replace_the_published_file(tmp_path):
    _write(tmp_path / "TN.json", "nss_snapshot", 5.0)
    with pytest.raises(ValueError, match="belongs under nss_snapshot/"):
        load_state("TN", data_dir=tmp_path)


@pytest.mark.parametrize("status, state", [("verified", "TN"), ("nss_snapshot", "KY")])
def test_snapshot_file_must_say_so(tmp_path, status, state):
    _write(tmp_path / NSS_SNAPSHOT_SUBDIR / "TN.json", status, 5.0, state=state)
    with pytest.raises(ValueError, match="must have status 'nss_snapshot'"):
        load_state("TN", data_dir=tmp_path, allow_nss_snapshot=True)


def test_unparsed_is_loaded_and_counted_in_advisory(tmp_path):
    _write(tmp_path / NSS_SNAPSHOT_SUBDIR / "TN.json", "nss_snapshot", 5.0, unparsed=UNPARSED)
    lib = load_state("TN", data_dir=tmp_path, allow_nss_snapshot=True)
    assert lib.unparsed == UNPARSED
    assert "2 published equations or regions could not be parsed" in lib.advisory
    # No retrieval time recorded: the advisory says so rather than inventing one.
    assert "on an unrecorded date" in lib.advisory


def test_unparsed_defaults_to_empty_for_other_statuses(tmp_path):
    _write(tmp_path / "TN.json", "verified", 2.0)
    assert load_state("TN", data_dir=tmp_path).unparsed == []


def test_listings_keep_the_sets_apart(data_dir):
    assert nss_snapshot_states(data_dir) == ["TN"]
    assert available_states(data_dir) == {}
    assert provisional_states(data_dir) == [] and awaiting_nss_states(data_dir) == []
    assert nss_snapshot_states(data_dir / "nowhere") == []


# ------------------------------------------------------------------------ generator


def test_asep_to_sep_log():
    expected = math.sqrt(math.log(1 + 0.467**2)) / math.log(10)
    assert asep_to_sep_log(46.7) == pytest.approx(expected)
    # Small errors: SE_P(%) ~ 100 * ln(10) * sep_log.
    assert asep_to_sep_log(1.0) == pytest.approx(0.01 / math.log(10), rel=1e-4)


def test_snapshot_to_library_builds_equations_and_unparsed(trimmed):
    lib = parse_state(snapshot_to_library(trimmed))
    assert lib.status == "nss_snapshot" and lib.state == "TN"
    assert lib.retrieved_utc == "2026-10-25T17:20:40+00:00"
    assert lib.issue is None
    assert lib.regions == ["GC2015"] and lib.aeps("GC2015") == [0.5, 0.01]

    eq = lib.equation("GC2015", 0.5)
    assert eq.intercept == pytest.approx(0.6866362692622934)
    assert eq.coefficients == pytest.approx((0.735, 4.068))
    assert eq.region_name == "TN_Peakflow_Region_1_SIR_2024_5130"
    assert eq.citation.publication == (
        "Ladd, D.E., and Ensminger, P.A. 2025, Estimating the magnitude and frequency of "
        "floods at ungaged locations on streams in Tennessee through the 2013 water year: "
        "U.S. Geological Survey Scientific Investigations Report 2024–5130, 19 p."
    )
    assert eq.citation.table == "NSS region 1062 (GC2015), citation 210"
    assert eq.citation.year == 2025
    assert eq.citation.url.startswith("http")
    assert eq.sep_log == pytest.approx(asep_to_sep_log(46.7))
    area = eq.variables[0]
    assert (area.code, area.units, area.minimum, area.maximum) == ("DRNAREA", "mi^2", 0.2, 21400.0)

    bad = trimmed["regions"][0]["equations"][2]
    assert [u["reason"] for u in lib.unparsed] == [
        bad["parse_error"],
        "unresolved: Estimate returned the region with no results",
    ]
    assert lib.unparsed[0] == {
        "region_code": "GC2015",
        "region_name": "TN_Peakflow_Region_1_SIR_2024_5130",
        "statistic_code": "PK0_2AEP",
        "aep": 0.002,
        "equation": bad["equation"],
        "reason": bad["parse_error"],
    }
    assert lib.unparsed[1]["region_code"] == "GC341" and lib.unparsed[1]["equation"] == ""


def test_source_reports_and_notes(trimmed):
    data = snapshot_to_library(trimmed)
    assert data is not None
    cites = trimmed["citations"]
    # One entry per citation a region uses; GC341, though unresolved, still cites 129.
    assert data["source_reports"] == [
        {
            "citation": f"{cites[cid]['author']} {cites[cid]['title']}",
            "url": cites[cid]["citation_url"],
            "nss": {
                "citation_id": int(cid),
                "last_year_of_data": cites[cid]["last_year_of_data"],
            },
        }
        for cid in ("210", "129")
    ]
    assert data["issue"] is None
    assert "data/nss_snapshots/TN_2026-10-25.json" in data["notes"]
    assert "2026-10-25T17:20:40+00:00" in data["notes"]
    assert "2 equations stored" in data["notes"] and "2 unparsed" in data["notes"]


def test_equations_round_trip_through_json(trimmed):
    data = snapshot_to_library(trimmed)
    assert data is not None
    again = parse_state(json.loads(json.dumps(data)))
    assert again.equations == parse_state(data).equations


def test_parse_uses_the_equation_string_not_the_stored_parse(trimmed):
    stale = copy.deepcopy(trimmed)
    stale["regions"][0]["equations"][0]["parsed"] = None
    stale["regions"][0]["equations"][0]["parse_error"] = "refused by an older parser"
    lib = parse_state(snapshot_to_library(stale))
    assert lib.equation("GC2015", 0.5).intercept == pytest.approx(0.6866362692622934)
    assert len(lib.unparsed) == 2


@pytest.mark.parametrize("asep", [0, 0.0, None])
def test_asep_zero_or_missing_leaves_sep_log_unset(trimmed, asep):
    trimmed["regions"][0]["equations"][0]["asep_pct"] = asep
    lib = parse_state(snapshot_to_library(trimmed))
    assert lib.equation("GC2015", 0.5).sep_log is None
    assert lib.equation("GC2015", 0.01).sep_log is not None


def test_snapshot_without_regions_returns_none():
    snapshot = {"state": "DC", "retrieved_utc": "2026-10-25T17:20:40+00:00", "regions": []}
    assert snapshot_to_library(snapshot) is None


def test_non_aep_statistics_are_skipped_and_counted(trimmed):
    extra = copy.deepcopy(trimmed["regions"][0]["equations"][0])
    extra.update(statistic_code="PKUNIT", aep=None)
    trimmed["regions"][0]["equations"].append(extra)
    data = snapshot_to_library(trimmed)
    assert data is not None
    lib = parse_state(data)
    assert len(lib.equations) == 2 and len(lib.unparsed) == 2
    assert "1 non-AEP statistic skipped" in data["notes"]


def test_regression_equation_to_dict_round_trips(trimmed):
    lib = parse_state(snapshot_to_library(trimmed))
    for eq in lib.equations:
        assert RegressionEquation.from_dict(eq.to_dict()) == eq


# --------------------------------------------------------------------------- generator tool


def _load_build_tool():
    path = Path(__file__).resolve().parent.parent / "tools" / "build_nss_snapshot_library.py"
    spec = importlib.util.spec_from_file_location("build_nss_snapshot_library", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run_build_tool(tmp_path: Path, snapshots: Any, *extra: str) -> int:
    published = tmp_path / "published"
    published.mkdir(exist_ok=True)
    (published / "WY.json").write_text("{}", encoding="utf-8")
    argv = [a for p in snapshots for a in ("--snapshot", str(p))]
    argv += ["--out", str(tmp_path / "out"), "--published-dir", str(published), *extra]
    result: int = _load_build_tool().main(argv)
    return result


def test_build_tool_writes_skips_and_reports(tmp_path, capsys, trimmed):
    snaps = tmp_path / "snaps"
    snaps.mkdir()
    (snaps / "TN.json").write_text(json.dumps(trimmed), encoding="utf-8")
    (snaps / "WY.json").write_text(json.dumps(dict(trimmed, state="WY")), encoding="utf-8")
    dc = {"state": "DC", "retrieved_utc": "2026-10-25T17:20:40+00:00", "regions": []}
    (snaps / "DC.json").write_text(json.dumps(dc), encoding="utf-8")
    out = tmp_path / "out"

    assert _run_build_tool(tmp_path, [snaps / f"{c}.json" for c in ("TN", "WY", "DC")]) == 0

    assert capsys.readouterr().out.splitlines() == [
        "TN: 2 equations, 2 unparsed",
        "skipped WY: published set exists",
        "skipped DC: no NSS peak-flow regions",
    ]
    assert sorted(p.name for p in out.iterdir()) == ["TN.json"]
    text = (out / "TN.json").read_text(encoding="utf-8")
    assert text.endswith("}\n")
    assert json.loads(text) == snapshot_to_library(trimmed)
    assert len(parse_state(json.loads(text)).equations) == 2


def test_build_tool_reports_bad_snapshots_and_continues(tmp_path, capsys, trimmed):
    snaps = tmp_path / "snaps"
    snaps.mkdir()
    bad_json = snaps / "bad.json"
    bad_json.write_text("{not json", encoding="utf-8")
    no_state = snaps / "nostate.json"
    no_state.write_text(json.dumps({"regions": []}), encoding="utf-8")
    first = snaps / "TN_a.json"
    first.write_text(json.dumps(trimmed), encoding="utf-8")
    second = snaps / "TN_b.json"
    second.write_text(json.dumps(dict(trimmed, retrieved_utc="2026-10-26T00:00:00+00:00")))

    assert _run_build_tool(tmp_path, [bad_json, no_state, first, second]) == 1

    lines = capsys.readouterr().out.splitlines()
    assert lines[0].startswith(f"error {bad_json}: cannot read snapshot (JSONDecodeError")
    assert lines[1].startswith(f"error {no_state}: cannot read snapshot (KeyError")
    assert lines[2] == "TN: 2 equations, 2 unparsed"
    assert lines[3] == f"error TN: two snapshots given, {first} and {second}; skipped {second}"
    written = json.loads((tmp_path / "out" / "TN.json").read_text(encoding="utf-8"))
    assert written == snapshot_to_library(trimmed)
