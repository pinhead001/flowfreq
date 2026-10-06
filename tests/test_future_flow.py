"""Tests for flowfreq.future_flow.

Factor values are synthetic except in the NCHRP 15-61 worked-example tests, which
transcribe the Guide's published tables (and, for Table 7.1, the USGS equations it used)
into the test. None of those values ships as data.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import pandas as pd
import pytest

from flowfreq.future_flow import (
    ChangeFactorSet,
    apply_change_factors,
    apply_with_precedence,
    available_factor_sets,
    ensemble_ratio_summary,
    factor_set_from_dict,
    index_flood_projection,
    national_source_review,
    regression_change_factors,
    select_factor_set,
)
from flowfreq.regression.equations import Citation, OutOfRangeError, RegressionEquation, Variable

GUIDANCE_DOC = Path(__file__).resolve().parents[1] / "docs" / "FUTURE_FLOW_GUIDANCE.md"


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


def _round_sig(x: float, sig: int = 3) -> float:
    return round(x, sig - 1 - int(math.floor(math.log10(abs(x)))))


def _half_unit(x: float, sig: int = 3) -> float:
    """Half a unit in the last place of ``x`` printed to ``sig`` significant figures."""
    return 0.5 * 10.0 ** (int(math.floor(math.log10(abs(x)))) - sig + 1)


# --------------------------------------------------------------------------
# Precedence: the national result is reported alongside a state override
# --------------------------------------------------------------------------


def test_apply_with_precedence_reports_national_alongside_override():
    nat = _set()
    wa = _set(
        name="synthetic-wa",
        level="state",
        geography=("WA",),
        factors=(1.10, 1.20, 1.40),
        legal_status="required",
    )
    fq = apply_with_precedence(QUANTILES, nat, [wa], state="WA")
    t = fq.table
    assert fq.factor_set is wa and fq.overrode == "synthetic-national"
    assert fq.national is not None and fq.national.factor_set is nat
    assert list(t["future_flow_cfs"][:3]) == pytest.approx([110.0, 240.0, 560.0])
    assert list(t["national_future_flow_cfs"][:3]) == pytest.approx([105.0, 220.0, 480.0])
    assert list(t["national_change_factor"][:3]) == pytest.approx([1.05, 1.10, 1.20])
    assert list(t["national_extrapolated"]) == [False, False, False, True]
    assert list(t["regulatory_flow_cfs"]) == list(QUANTILES["flow_cfs"])  # never overwritten
    prov = fq.provenance
    assert prov["basis"] == "state override of national default 'synthetic-national'"
    assert prov["national_alongside"] == "synthetic-national"
    assert prov["legal_status"] == "required"
    assert any("overrides" in n for n in fq.notes)


def test_apply_with_precedence_national_only_where_no_state_applies():
    nat = _set()
    wa = _set(name="synthetic-wa", level="state", geography=("WA",))
    fq = apply_with_precedence(QUANTILES, nat, [wa], state="OR")
    assert fq.factor_set is nat and fq.national is None and fq.overrode is None
    assert "national_future_flow_cfs" not in fq.table.columns
    assert fq.provenance["basis"] == "national default"
    assert fq.provenance["national_alongside"] is None


def test_apply_with_precedence_without_national_default():
    wa = _set(name="synthetic-wa", level="state", geography=("WA",))
    fq = apply_with_precedence(QUANTILES, None, [wa], state="WA")
    assert fq.factor_set is wa and fq.national is None
    assert fq.provenance["basis"].startswith("state guidance")
    with pytest.raises(LookupError):
        apply_with_precedence(QUANTILES, None, [wa], state="OR")
    with pytest.raises(ValueError, match="2 state sets"):
        apply_with_precedence(
            QUANTILES, None, [wa, _set(name="b", level="state", geography=("WA",))], state="WA"
        )


# --------------------------------------------------------------------------
# NCHRP 15-61 Ch. 6: ensemble ratio summary, Tables 6.6-6.9 (Design Practices, pp. 54-56)
# --------------------------------------------------------------------------

# Table 6.6 (2000-2049) and Table 6.7 (2050-2099): future/baseline 24-hour precipitation
# quantile ratios for 12 BCCA-downscaled CMIP5 models, RCP6.0, Denver CO (p. 54).
NCHRP_6_6 = {
    0.5: [1.21, 1.00, 0.99, 0.98, 1.10, 0.97, 1.10, 1.10, 1.10, 1.04, 1.02, 0.90],
    0.1: [1.17, 0.91, 0.98, 1.01, 1.23, 1.03, 1.06, 1.02, 1.05, 1.01, 0.87, 1.05],
    0.04: [1.13, 0.86, 1.00, 1.02, 1.30, 1.08, 1.04, 0.97, 1.00, 1.00, 0.78, 1.18],
    0.01: [1.07, 0.79, 1.02, 1.04, 1.41, 1.17, 1.03, 0.89, 0.91, 0.98, 0.67, 1.41],
}
NCHRP_6_7 = {
    0.5: [1.11, 1.04, 1.04, 1.03, 1.19, 0.98, 0.97, 1.24, 1.19, 0.96, 1.12, 0.97],
    0.1: [1.15, 0.97, 1.07, 0.96, 1.26, 1.00, 1.03, 1.15, 1.15, 0.95, 0.95, 1.06],
    0.04: [1.13, 0.90, 1.10, 0.92, 1.32, 1.01, 1.10, 1.09, 1.12, 0.93, 0.87, 1.11],
    0.01: [1.09, 0.81, 1.15, 0.88, 1.41, 1.03, 1.21, 0.99, 1.07, 0.88, 0.76, 1.20],
}
# NOAA Atlas 14 24-hour depths, inches (Table 6.2, p. 51; repeated in Tables 6.8-6.9).
NCHRP_6_2 = {0.5: 1.84, 0.1: 3.00, 0.04: 3.72, 0.01: 4.89}
# Tables 6.8 and 6.9 (p. 56): mean, SD, 5% CL, 95% CL of ratios; projected depth (in).
NCHRP_6_8 = {
    0.5: (1.04, 0.082, 1.01, 1.08, 1.92),
    0.1: (1.03, 0.098, 0.99, 1.08, 3.10),
    0.04: (1.03, 0.137, 0.97, 1.10, 3.84),
    0.01: (1.03, 0.220, 0.93, 1.14, 5.05),
}
NCHRP_6_9 = {
    0.5: (1.07, 0.099, 1.02, 1.12, 1.97),
    0.1: (1.06, 0.101, 1.01, 1.11, 3.18),
    0.04: (1.05, 0.128, 0.99, 1.11, 3.94),
    0.01: (1.04, 0.189, 0.95, 1.13, 5.18),
}


@pytest.mark.parametrize("ratios,published", [(NCHRP_6_6, NCHRP_6_8), (NCHRP_6_7, NCHRP_6_9)])
def test_reproduces_nchrp_15_61_tables_6_8_and_6_9(ratios, published):
    """Tables 6.6/6.7 summarise to Tables 6.8/6.9.

    The inputs are the published ratios, rounded to two decimals, so the
    statistics agree to that rounding rather than exactly.
    """
    s = ensemble_ratio_summary(ratios, historical=NCHRP_6_2).set_index("aep")
    for aep, (mean, sd, lo, hi, projected) in published.items():
        row = s.loc[aep]
        assert row["n_models"] == 12
        assert round(row["mean_ratio"], 2) == pytest.approx(mean)
        assert row["sd_ratio"] == pytest.approx(sd, abs=0.0015)
        assert row["lower_cl"] == pytest.approx(lo, abs=0.008)
        assert row["upper_cl"] == pytest.approx(hi, abs=0.008)
        assert row["projected"] == pytest.approx(projected, abs=0.006)
    # Eq. 6.3: the 0.1-AEP ratio is used for the rarer AEPs (table footnotes).
    assert list(s["substituted"]) == [False, False, True, True]
    assert s.loc[0.01, "applied_ratio"] == s.loc[0.1, "mean_ratio"]


def test_ensemble_text_example_6_9():
    # p. 56: "the projected 0.01 AEP precipitation estimate in Table 6.9 is 1.06 x 4.89 = 5.18".
    s = ensemble_ratio_summary(NCHRP_6_7, historical=NCHRP_6_2).set_index("aep")
    assert round(s.loc[0.01, "applied_ratio"], 2) == 1.06
    assert round(s.loc[0.01, "projected"], 2) == 5.18


def test_ensemble_ratio_summary_errors():
    with pytest.raises(ValueError, match="Eq. 6.3"):
        ensemble_ratio_summary({0.5: [1.0, 1.1], 0.01: [1.0, 1.2]})
    with pytest.raises(ValueError, match="two model ratios"):
        ensemble_ratio_summary({0.5: [1.0]})
    with pytest.raises(ValueError, match="positive"):
        ensemble_ratio_summary({0.5: [1.0, -1.0]})
    with pytest.raises(ValueError, match="missing"):
        ensemble_ratio_summary({0.5: [1.0, 1.1]}, historical={0.1: 2.0})
    with pytest.raises(ValueError, match="confidence"):
        ensemble_ratio_summary({0.5: [1.0, 1.1]}, confidence=1.5)


# --------------------------------------------------------------------------
# NCHRP 15-61 Ch. 7: projected MAP in a USGS regression equation, Table 7.1 (p. 80)
# --------------------------------------------------------------------------

# Lumia, Freehafer and Smith (2006), USGS SIR 2006-5112, Table 1 (p. 32), Region 5
# (southwestern New York): Q = c * A^a * SL^b * P^d. The Guide's Eq. 7.3 (p. 79) is the
# 50-year row. Ranges are the Guide's Step 2 (p. 80).
LUMIA_R5 = {  # AEP: (c, a, b, d)
    0.5: (0.083, 0.965, 0.431, 1.305),
    0.2: (0.322, 0.965, 0.498, 0.995),
    0.1: (0.597, 0.967, 0.538, 0.853),
    0.04: (1.05, 0.972, 0.581, 0.724),
    0.02: (1.46, 0.976, 0.610, 0.651),
    0.01: (1.91, 0.980, 0.636, 0.590),
}
NY_HIST = {"DRNAREA": 100.0, "CSL1085LFP": 30.0, "PRECIP": 36.0}


def _lumia_r5():
    cite = Citation("Lumia, Freehafer and Smith (2006), USGS SIR 2006-5112", "Table 1", 2006)
    vars_ = (
        Variable("DRNAREA", minimum=1.7, maximum=4770),
        Variable("CSL1085LFP", minimum=2.76, maximum=223),
        Variable("PRECIP", minimum=31.6, maximum=49.8),
    )
    return [
        RegressionEquation("NY5", aep, math.log10(c), vars_, (a, b, d), cite)
        for aep, (c, a, b, d) in LUMIA_R5.items()
    ]


def _table_7_1():
    return regression_change_factors(
        _lumia_r5(),
        NY_HIST,
        {"PRECIP": 1.25 * 36.0},  # Eq. 7.2 with an ensemble MAP ratio of 1.25 (p. 80)
        name="nchrp-15-61-sec-7.2",
        scenario="hypothetical MAP ratio 1.25",
        horizon="unspecified",
        citation="NCHRP 15-61 Design Practices, Section 7.2",
    )


def test_reproduces_nchrp_15_61_table_7_1_from_the_equations():
    """The Ch. 7 procedure, run on Lumia's equations, gives the Guide's Table 7.1.

    Five of six columns agree at the table's three significant figures. At the
    0.2 AEP the table's 5,060/6,320 are both about 4 percent below what Region 5's
    printed 5-year equation gives (5,270/6,580), while their ratio (1.249) is the
    equation's 1.25**0.995. That is an inconsistency inside the published table,
    pinned here so a later change to the transcription is noticed.
    """
    proj = _table_7_1()
    t = proj.table.set_index("aep")
    published = dict(zip(NCHRP_7_1_AEP, zip(NCHRP_7_1_HIST, NCHRP_7_1_PROJ)))
    for aep, (hist, fut) in published.items():
        # The factor lies within what the table's three-figure rounding allows.
        lo = (fut - _half_unit(fut)) / (hist + _half_unit(hist))
        hi = (fut + _half_unit(fut)) / (hist - _half_unit(hist))
        assert lo <= t.loc[aep, "change_factor"] <= hi
        if aep == 0.2:
            assert t.loc[aep, "historical_flow_cfs"] == pytest.approx(5272, abs=1)
            assert t.loc[aep, "historical_flow_cfs"] / hist == pytest.approx(1.04, abs=0.005)
            continue
        assert _round_sig(t.loc[aep, "historical_flow_cfs"]) == hist
        assert _round_sig(t.loc[aep, "projected_flow_cfs"]) == fut
    # p. 80: 10,700 -> 12,400 cfs, "an estimated 16-percent increase" at the 0.02 AEP.
    assert round(t.loc[0.02, "change_factor"] - 1, 2) == 0.16
    assert proj.changed == {"PRECIP": (36.0, 45.0)}
    fs = proj.factor_set
    assert fs.level == "site" and "Ch. 7" in fs.citation
    q = pd.DataFrame({"aep": list(t.index), "flow_cfs": list(t["historical_flow_cfs"])})
    out = apply_change_factors(q, fs).table
    assert list(out["future_flow_cfs"]) == pytest.approx(list(t["projected_flow_cfs"]))


def test_regression_change_factors_refuses_out_of_range_projection():
    # Step 3: the projected MAP must stay inside the equation's 31.6-49.8 in range.
    with pytest.raises(OutOfRangeError, match="PRECIP"):
        regression_change_factors(
            _lumia_r5(),
            NY_HIST,
            {"PRECIP": 50.0},
            name="n",
            scenario="s",
            horizon="h",
            citation="c",
        )


def test_regression_change_factors_needs_a_projected_variable():
    eqs = _lumia_r5()
    no_precip = RegressionEquation(
        "X", 0.01, 1.0, (Variable("DRNAREA"),), (0.8,), Citation("synthetic", "T1")
    )
    kw = dict(name="n", scenario="s", horizon="h", citation="c")
    with pytest.raises(ValueError, match="cannot estimate future"):
        regression_change_factors([no_precip], NY_HIST, {"PRECIP": 40.0}, **kw)
    with pytest.raises(ValueError, match="nothing to project"):
        regression_change_factors(eqs, NY_HIST, {"PRECIP": 36.0}, **kw)
    with pytest.raises(ValueError, match="duplicate"):
        regression_change_factors(eqs + eqs[:1], NY_HIST, {"PRECIP": 40.0}, **kw)
    with pytest.raises(KeyError):
        regression_change_factors(eqs, NY_HIST, {"T_MEAN": 40.0}, **kw)
    with pytest.raises(ValueError, match="no equations"):
        regression_change_factors([], NY_HIST, {"PRECIP": 40.0}, **kw)


# --------------------------------------------------------------------------
# NCHRP 15-61 Ch. 8: index-flood method, Table 8.1 (p. 87)
# --------------------------------------------------------------------------

# Pajarito Creek at Newkirk, NM (07225000); A = 55 mi2, MAP = 14.1 in (p. 86).
NCHRP_8_1_AEP = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01)
NCHRP_8_1_HIST = (940, 1840, 2490, 3330, 3990, 4650)
NCHRP_8_1_RATIO = (0.38, 0.74, 1.00, 1.34, 1.60, 1.87)
NCHRP_8_1_PROJ = (1040, 2030, 2750, 3680, 4410, 5140)


def _eq_8_1():
    # Guide Eq. 8.1 (p. 87), from Asquith and Thompson (2008): Q0.1 = 111 A^0.5311 MAP^0.5469.
    return RegressionEquation(
        "NM-AT2008",
        0.1,
        math.log10(111),
        (Variable("DRNAREA"), Variable("PRECIP")),
        (0.5311, 0.5469),
        Citation("NCHRP 15-61 Design Practices, quoting Asquith and Thompson (2008)", "Eq. 8.1"),
    )


def test_reproduces_nchrp_15_61_table_8_1():
    """Ch. 7 gives the projected index flood; Ch. 8 scales the gauged curve by it."""
    hist = {"DRNAREA": 55.0, "PRECIP": 14.1}
    reg = regression_change_factors(
        [_eq_8_1()],
        hist,
        {"PRECIP": 1.2 * 14.1},  # "a projected 20-percent increase in MAP"
        name="pajarito-index",
        scenario="hypothetical MAP +20 percent",
        horizon="unspecified",
        citation="NCHRP 15-61 Design Practices, Section 8.2",
    ).table.iloc[0]
    # p. 87: "an estimated 0.1 AEP discharge of 3,964 ft3/s ... (after the increase) of 4,380".
    assert round(reg["historical_flow_cfs"]) == 3964
    assert _round_sig(reg["projected_flow_cfs"]) == 4380
    curve = pd.DataFrame({"aep": NCHRP_8_1_AEP, "flow_cfs": NCHRP_8_1_HIST})
    fq = index_flood_projection(
        curve,
        historical_index_estimate=reg["historical_flow_cfs"],
        projected_index_estimate=reg["projected_flow_cfs"],
        name="pajarito-ffc",
        scenario="hypothetical MAP +20 percent",
        horizon="unspecified",
        citation="NCHRP 15-61 Design Practices, Table 8.1",
    )
    t = fq.table
    assert [round(r, 2) for r in t["flood_ratio"]] == list(NCHRP_8_1_RATIO)
    assert [_round_sig(q) for q in t["future_flow_cfs"]] == list(NCHRP_8_1_PROJ)
    assert list(t["regulatory_flow_cfs"]) == list(NCHRP_8_1_HIST)
    assert not t["extrapolated"].any()
    assert fq.factor_set.level == "site" and "Ch. 8" in fq.factor_set.citation
    assert fq.provenance["basis"] == "site-specific derivation"


def test_index_flood_projection_errors():
    curve = pd.DataFrame({"aep": [0.5, 0.01], "flow_cfs": [100.0, 300.0]})
    kw = dict(name="n", scenario="s", horizon="h", citation="c")
    with pytest.raises(ValueError, match="index AEP"):
        index_flood_projection(
            curve, historical_index_estimate=1.0, projected_index_estimate=1.1, **kw
        )
    with pytest.raises(ValueError, match="positive"):
        index_flood_projection(
            curve, historical_index_estimate=0.0, projected_index_estimate=1.1, **kw
        )
    with pytest.raises(KeyError):
        index_flood_projection(
            curve.rename(columns={"flow_cfs": "q"}),
            historical_index_estimate=1.0,
            projected_index_estimate=1.1,
            **kw,
        )


# --------------------------------------------------------------------------
# State survey: Wave 1 and Wave 2 states have no tabulated factors
# --------------------------------------------------------------------------

SURVEYED_STATES = ("WA", "OR", "ID", "MT", "CO", "UT", "WY", "NM", "AZ", "NV")


def test_no_state_factor_set_ships():
    assert available_factor_sets(level="state") == ()


@pytest.mark.parametrize("code", SURVEYED_STATES)
def test_guidance_doc_records_each_surveyed_state(code):
    text = GUIDANCE_DOC.read_text(encoding="utf-8")
    m = re.search(rf"^### [^\n]*\({code}\)[^\n]*$(.*?)(?=^##|\Z)", text, re.M | re.S)
    assert m, f"no section for {code}"
    section = m.group(1)
    assert "no tabulated factors" in section.lower()
    assert "Surveyed 2026-10" in section
    # The status table carries the same verdict.
    row = re.search(rf"^\| [^|]*\({code}\)[^\n]*$", text, re.M)
    assert row and "no tabulated factors" in row.group(0).lower()
