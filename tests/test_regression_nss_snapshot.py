"""Equation sets harvested from the USGS NSS web service: the ``nss_snapshot`` tier.

Status ``nss_snapshot`` lives under ``nss_snapshot/`` and, like the other opt-in sets,
is returned only on opt-in (``allow_nss_snapshot=True``), so it never changes the
default ``load_state``. Its advisory says the equations are not checked against the
source report.

``tests/fixtures/nss_snapshot_TN_trimmed.json`` is the 2026-10-25 Tennessee snapshot
trimmed to two regions: GC2015 with three equations, and GC341. Two things in it are
invented, as its ``note`` says: GC341 was flipped to ``unresolved`` with a made-up
reason, and GC2015's 0.2% AEP equation string is synthetic -- a power whose exponent
is not affine in the variables (a squared characteristic), a form the parser refuses --
so that the fixture carries one refusal; its ``parse_error`` is the parser's message for that string. (The real string
is ``65.8*(DRNAREA^0.701)*(CLIMFAC2YR^2.928)``.) The printed-constant case is built
inline by ``test_constant_q0_is_stored_as_a_constant``.
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
REPO_ROOT = Path(__file__).resolve().parent.parent

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


def test_stored_order_is_stable(trimmed):
    # Reverse the harvest order: the stored order must not follow it.
    trimmed["regions"].reverse()
    for region in trimmed["regions"]:
        region["equations"] = list(reversed(region.get("equations") or []))
    data = snapshot_to_library(trimmed)
    assert data is not None
    assert [(e["region_code"], e["aep"]) for e in data["equations"]] == [
        ("GC2015", 0.5),
        ("GC2015", 0.01),
    ]
    assert [(u["region_code"], u["statistic_code"]) for u in data["unparsed"]] == [
        ("GC2015", "PK0_2AEP"),
        ("GC341", ""),
    ]


def test_constant_q0_is_stored_as_a_constant(trimmed):
    # NSS's form of a printed Q=0 (CA GC30): parse_equation refuses it, and the
    # parse_constant_equation retry is the only way form="constant" enters a library.
    q0 = copy.deepcopy(trimmed["regions"][0]["equations"][0])
    q0.update(
        statistic_code="PK10AEP",
        aep=0.1,
        equation="0*(0*(DRNAREA))^(0)",
        asep_pct=30.0,
        parsed=None,
        parse_error=None,
    )
    trimmed["regions"][0]["equations"].append(q0)
    data = snapshot_to_library(trimmed)
    assert data is not None
    lib = parse_state(data)
    eq = lib.equation("GC2015", 0.1)
    assert eq.form == "constant" and eq.constant_cfs == 0.0
    # A positive ASEp does not give a constant an error.
    assert eq.sep_log is None
    assert eq.citation.table == "NSS region 1062 (GC2015), citation 210"
    assert eq.citation.year == 2025
    assert lib.aeps("GC2015") == [0.5, 0.1, 0.01]
    # Appended last in the harvest, stored between the 50% and 1% equations.
    assert [e["aep"] for e in data["equations"]] == [0.5, 0.1, 0.01]
    assert [u["statistic_code"] for u in lib.unparsed] == ["PK0_2AEP", ""]


@pytest.mark.parametrize(
    "author, title, year",
    [
        # LA citation 62: the title opens with the report number, not the year.
        ("TR No.", "60, (1998) Floods in Louisiana, magnitude and frequency", 1998),
        # PR citation 172: the year is in the author field.
        (
            "Ryan, P.J., and others, 2021,",
            "Magnitude and frequency of floods in Puerto Rico through the 2019 water year",
            2021,
        ),
        ("Anon.", "Flood-frequency equations, no date given", None),
    ],
)
def test_citation_year_is_the_first_year_in_author_and_title(trimmed, author, title, year):
    trimmed["citations"]["210"].update(author=author, title=title)
    lib = parse_state(snapshot_to_library(trimmed))
    assert lib.equation("GC2015", 0.5).citation.year == year


def test_parse_state_reads_nulls_as_empty(tmp_path):
    data: Dict[str, Any] = {
        "schema_version": 1,
        "state": "TN",
        "status": "nss_snapshot",
        "source_reports": [{"citation": "Synthetic report"}],
        "equations": [EQ],
        "retrieved_utc": None,
        "notes": None,
        "unparsed": None,
    }
    lib = parse_state(data)
    assert lib.retrieved_utc == "" and lib.notes == "" and lib.unparsed == []
    path = tmp_path / NSS_SNAPSHOT_SUBDIR / "TN.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(data), encoding="utf-8")
    advisory = load_state("TN", data_dir=tmp_path, allow_nss_snapshot=True).advisory
    assert "on an unrecorded date" in advisory and "None" not in advisory


def test_notes_count_one_unparsed_in_the_singular(trimmed):
    trimmed["regions"] = trimmed["regions"][:1]
    data = snapshot_to_library(trimmed)
    assert data is not None
    assert "1 unparsed equation or region," in data["notes"]


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
    assert "2 equations stored" in data["notes"]
    assert "2 unparsed equations or regions," in data["notes"]


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


def test_asep_below_one_leaves_sep_log_unset_and_notes_the_region(trimmed):
    trimmed["regions"][0]["equations"][0]["asep_pct"] = 0.5
    data = snapshot_to_library(trimmed)
    lib = parse_state(data)
    assert lib.equation("GC2015", 0.5).sep_log is None
    assert lib.equation("GC2015", 0.01).sep_log is not None
    assert "GC2015" in data["notes"] and "not converted" in data["notes"]


def test_asep_as_fraction_region_is_converted_from_percent(trimmed):
    # #213: in the ASEP_AS_FRACTION regions NSS stores the report's SEP percent / 100.
    trimmed["regions"][0]["region_code"] = "GC1889"
    trimmed["regions"][0]["equations"][0]["asep_pct"] = 0.446
    data = snapshot_to_library(trimmed)
    lib = parse_state(data)
    assert lib.equation("GC1889", 0.5).sep_log == pytest.approx(asep_to_sep_log(44.6))
    assert "GC1889" in data["notes"] and "SIR 2024-5031" in data["notes"]
    assert "not converted" not in data["notes"]


@pytest.mark.parametrize(
    "code, region, aep, sep_pct",
    [
        # SIR 2021-5046 Table 2 (SEP, percent) and SIR 2024-5031 Table 3.
        ("LA", "GC1888", 0.5, 45.0),
        ("MS", "GC1912", 0.01, 55.0),
        ("AR", "GC1917", 0.002, 61.0),
        ("LA", "GC1889", 0.5, 44.60),
        ("LA", "GC1889", 0.002, 52.69),
    ],
)
def test_shipped_sep_log_matches_the_report(code, region, aep, sep_pct):
    eq = load_state(code, allow_nss_snapshot=True).equation(region, aep)
    assert eq.sep_log == pytest.approx(asep_to_sep_log(sep_pct), rel=1e-9)


@pytest.mark.parametrize("code", nss_snapshot_states())
def test_every_stored_sep_log_is_plausible(code):
    # Bounds are wider than the 0.05-0.6 first proposed: genuine percentage ASEps in the
    # sweep give AL 0.039-0.048 (ASEp ~9-11%) and CA GC1534 0.61-0.90 (ASEp ~150-856%).
    # The unconverted bad values would have been ~0.002, far below the floor.
    for eq in load_state(code, allow_nss_snapshot=True).equations:
        if eq.sep_log is not None:
            assert 0.03 <= eq.sep_log <= 1.0, (code, eq.region_code, eq.aep, eq.sep_log)


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


def test_build_tool_skips_a_snapshot_with_no_parseable_equations(tmp_path, capsys, trimmed):
    """Every region unresolved (Indiana before its gating fix): no empty library is written."""
    unresolved = copy.deepcopy(trimmed)
    unresolved["state"] = "IN"
    for region in unresolved["regions"]:
        region["status"] = "unresolved"
        region["unresolved_reason"] = "Estimate returned no results"
    snap = tmp_path / "IN.json"
    snap.write_text(json.dumps(unresolved), encoding="utf-8")

    assert _run_build_tool(tmp_path, [snap]) == 0

    assert capsys.readouterr().out.splitlines() == ["skipped IN: no parseable equations"]
    assert not (tmp_path / "out" / "IN.json").exists()


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
    second.write_text(
        json.dumps(dict(trimmed, retrieved_utc="2026-10-26T00:00:00+00:00")), encoding="utf-8"
    )

    assert _run_build_tool(tmp_path, [bad_json, no_state, first, second]) == 1

    lines = capsys.readouterr().out.splitlines()
    assert lines[0].startswith(f"error {bad_json}: cannot read snapshot (JSONDecodeError")
    assert lines[1].startswith(f"error {no_state}: cannot read snapshot (KeyError")
    # Every snapshot is read before any is written (an extra NSS region is merged into
    # its state's first, #246), so the duplicate is reported before the write.
    assert lines[2] == f"error TN: two snapshots given, {first} and {second}; skipped {second}"
    assert lines[3] == "TN: 2 equations, 2 unparsed"
    written = json.loads((tmp_path / "out" / "TN.json").read_text(encoding="utf-8"))
    assert written == snapshot_to_library(trimmed)


# --- the committed library: 41 jurisdictions harvested 2026-10-25 -------------------


def _midpoint(v) -> float:
    lo, hi = v.minimum, v.maximum
    if lo is None or hi is None:
        return 1.0
    return math.sqrt(lo * hi) if lo > 0 else (lo + hi) / 2


@pytest.mark.parametrize("code", nss_snapshot_states())
def test_every_snapshot_library_loads_and_evaluates(code):
    from flowfreq.regression.equations import evaluate

    lib = load_state(code, allow_nss_snapshot=True)
    assert lib.nss_snapshot and lib.equations
    for eq in lib.equations:
        mid = {v.code: _midpoint(v) for v in (*eq.applicability, *eq.variables)}
        flow = evaluate(eq, mid).flow_cfs
        if eq.form == "constant":
            # CA GC30 at AEP 0.5 is a genuine printed Q = 0.
            assert flow >= 0, eq
        else:
            assert flow > 0, eq


def test_snapshot_states_are_the_41():
    assert len(nss_snapshot_states()) == 41


def test_unparsed_entries_are_the_known_2():
    # Observed in the 2026-10-25 sweep: 26 entries, of which 24 (GA, NC and SC, 8 each:
    # the SIR 2023-5006 region-indicator equations) parse since #211. The 2 left are
    # New York's two urban regions (WRI84_4350 Method 1 and Method 2 Urban Peak,
    # GC1077/GC1078), for which NSS Estimate returned no results; the library lists
    # each region, not equations, under unparsed. They are out of scope (#212): Sauer
    # and others' (1983) urban adjustments, each a power of a rural estimate.
    by_state = {
        c: len(load_state(c, allow_nss_snapshot=True).unparsed) for c in nss_snapshot_states()
    }
    assert {c: n for c, n in by_state.items() if n} == {"NY": 2}


@pytest.mark.parametrize("code, region", [("GA", "GC1934"), ("SC", "GC1935"), ("NC", "GC1936")])
def test_region_indicator_equations_match_nss(code, region):
    # #211: every SIR 2023-5006 equation reproduces NSS's own Estimate value at the
    # characteristics the harvest submitted, to NSS's 3-significant-figure rounding.
    from flowfreq.regression.equations import evaluate

    snapshot = json.loads(
        (REPO_ROOT / "data" / "nss_snapshots" / f"{code}_2026-10-25.json").read_text("utf-8")
    )
    (raw,) = [r for r in snapshot["regions"] if r["region_code"] == region]
    submitted = {p["code"]: p["submitted"] for p in raw["parameters"]}
    lib = load_state(code, allow_nss_snapshot=True)
    assert len(raw["equations"]) == 8
    for rec in raw["equations"]:
        eq = lib.equation(region, rec["aep"])
        assert any(v.times == "PCTREG2" for v in eq.variables)
        assert "PCTREG2" in [v.code for v in eq.applicability]
        assert "PCTREG4" not in [v.code for v in eq.variables]
        flow = evaluate(eq, submitted).flow_cfs
        assert flow == pytest.approx(rec["nss_value"], rel=5e-3), (code, rec["aep"])
        # 0 percent in region 4 (PCTREG4^0 in the 50% AEP form) is in range and evaluates.
        assert evaluate(eq, {**submitted, "PCTREG4": 0.0}).flow_cfs > 0


def _rrb_snapshot(trimmed):
    """A one-region snapshot shaped like NSS's ``RRB`` harvest (GC1775, citation 19)."""
    rrb = copy.deepcopy(trimmed)
    region = copy.deepcopy(rrb["regions"][0])
    region.update(region_code="GC1775", region_name="Rainy_River_Basin_Peakflow_SIR_2019_5012",
                  region_id=57, citation_id=19)  # fmt: skip
    rrb.update(state="RRB", regions=[region])
    rrb["citations"] = {"19": {"author": "Sanocki, C.A., and others,", "title": "2019, Rainy "
                               "River Basin peak flows: SIR 2019-5012", "citation_url": ""}}  # fmt: skip
    return rrb


def test_merge_snapshots_adds_the_extra_regions_and_citations(trimmed):
    from flowfreq.regression.nss_library import merge_snapshots

    merged = merge_snapshots(trimmed, [_rrb_snapshot(trimmed)])
    assert merged["state"] == "TN"
    codes = [r["region_code"] for r in merged["regions"]]
    assert codes[-1] == "GC1775" and codes[:-1] == [r["region_code"] for r in trimmed["regions"]]
    assert "19" in merged["citations"]
    lib = parse_state(snapshot_to_library(merged))
    eq = lib.equation("GC1775", 0.5)
    assert eq.citation.table == "NSS region 57 (GC1775), citation 19"
    assert "NSS region RRB's regions are merged in" in snapshot_to_library(merged)["notes"]
    # The inputs are not modified.
    assert all(r["region_code"] != "GC1775" for r in trimmed["regions"])


def test_merge_snapshots_refuses_a_clashing_region(trimmed):
    from flowfreq.regression.nss_library import merge_snapshots

    with pytest.raises(ValueError, match="GC2015"):
        merge_snapshots(trimmed, [copy.deepcopy(trimmed)])


def test_build_tool_merges_an_extra_region_into_its_state(tmp_path, capsys, trimmed):
    mn = copy.deepcopy(trimmed)
    mn["state"] = "MN"
    rrb = _rrb_snapshot(trimmed)
    paths = []
    for snap in (mn, rrb):
        p = tmp_path / f"{snap['state']}.snapshot.json"
        p.write_text(json.dumps(snap), encoding="utf-8")
        paths.append(p)
    out = tmp_path / "out"
    assert _run_build_tool(tmp_path, paths) == 0
    assert sorted(f.name for f in out.iterdir()) == ["MN.json"]
    lib = parse_state(json.loads((out / "MN.json").read_text(encoding="utf-8")))
    assert "GC1775" in lib.regions
    assert "RRB merged into MN" in capsys.readouterr().out


def test_build_tool_refuses_an_extra_region_without_its_state(tmp_path, capsys, trimmed):
    p = tmp_path / "RRB.snapshot.json"
    p.write_text(json.dumps(_rrb_snapshot(trimmed)), encoding="utf-8")
    assert _run_build_tool(tmp_path, [p]) == 1
    assert "RRB" in capsys.readouterr().out and not (tmp_path / "out" / "RRB.json").exists()


def test_minnesota_carries_the_rainy_river_basin_equations():
    # #246: GC1775 (SIR 2019-5012) is NSS region RRB's, merged into MN's library, and
    # reproduces NSS's own Estimate values at the harvest's submitted inputs.
    from flowfreq.regression.equations import evaluate

    snapshot = json.loads(
        (REPO_ROOT / "data" / "nss_snapshots" / "RRB_2026-10-31.json").read_text("utf-8")
    )
    (raw,) = snapshot["regions"]
    submitted = {p["code"]: p["submitted"] for p in raw["parameters"]}
    lib = load_state("MN", allow_nss_snapshot=True)
    assert len(raw["equations"]) == 8
    for rec in raw["equations"]:
        eq = lib.equation("GC1775", rec["aep"])
        assert "2019" in eq.citation.publication and "Rainy" in eq.citation.publication
        flow = evaluate(eq, submitted).flow_cfs
        assert flow == pytest.approx(rec["nss_value"], rel=5e-3), rec["aep"]
