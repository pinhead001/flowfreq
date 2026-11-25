"""Bulletin17B's inputs from an NWIS peak record, and the Task 5 review's minor points."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import pytest

from flowfreq.bulletin17b import Bulletin17B, Bulletin17BError, bulletin17b_kwargs, removed_peaks
from tests.fixtures.b17b import record

# 18 ordinary years, used to build records.
BASE = [1000.0, 1300, 800, 1150, 900, 1250, 700, 1050, 1400]
BASE = BASE + [q * 1.03 for q in BASE]
YEARS = list(range(1950, 1968))


def _steps(**kw):
    return Bulletin17B(**kw).run_analysis().steps


def synthetic(extra, **kw):
    flows = BASE + list(extra)
    return Bulletin17B(flows, list(range(1950, 1950 + len(flows))), **kw).run_analysis()


# --- minors from the Task 5 review ---------------------------------------------------------


def test_truncation_level_is_not_called_the_largest_flow_when_flows_are_unknown():
    """Below-base years have no usable flow, so only known flows are "the largest"."""
    res = Bulletin17B(BASE, YEARS, below_base_years=[1968]).run_analysis()
    text = next(s for s in res.steps if s.startswith("Made the conditional"))
    assert "largest" not in text and "smallest peak kept" in text
    res = Bulletin17B(BASE + [0.0], YEARS + [1969], below_base_years=[1968]).run_analysis()
    text = next(s for s in res.steps if s.startswith("Made the conditional"))
    assert "largest known flow truncated" in text


def test_hand_computation_docstring_names_the_conditional_and_k_skews():
    doc = Bulletin17B.__doc__ or ""
    part = doc[doc.index("hand_computation") :]
    assert "conditional" in part and "eq. 5-4" in part


def test_hand_step_says_which_curve_k_is_read_for():
    res = synthetic([0.0], hand_computation=True)
    text = next(s for s in res.steps if s.startswith("The conditional curve"))
    assert "K is read at skew" not in text
    assert "K for the conditional curve is read at skew" in text and "in eq. 5-4" in text


def test_steps_call_the_synthetic_skew_by_its_name():
    res = synthetic([0.0], generalized_skew=0.1)
    text = next(s for s in res.steps if s.startswith("Weighted"))
    assert "synthetic skew" in text and "station skew" not in text
    res = synthetic([0.0])
    text = next(s for s in res.steps if s.startswith("No generalized skew"))
    assert "synthetic skew" in text and "station skew" not in text
    res = synthetic([], generalized_skew=0.1)
    assert any(s.startswith("Weighted the station skew") for s in res.steps)


@pytest.mark.parametrize("name", ["water_years", "below_base_years"])
@pytest.mark.parametrize(
    "bad",
    [["1987"], np.array([True, False]), [np.True_], "1987", b"1987"],
    ids=["str items", "numpy mask", "numpy bool", "str", "bytes"],
)
def test_str_and_boolean_years_refuse(name, bad):
    kw = {"peak_flows": BASE, "water_years": YEARS}
    if name == "water_years":
        kw["water_years"] = bad
        kw["peak_flows"] = BASE[: len(bad)]
    else:
        kw["below_base_years"] = bad
    with pytest.raises(Bulletin17BError):
        Bulletin17B(**kw)


def test_the_year_count_in_the_too_few_message_is_years():
    with pytest.raises(Bulletin17BError, match="of the 11 years of the systematic record"):
        Bulletin17B(
            BASE[:9] + [0.0], list(range(1950, 1960)), below_base_years=[1960]
        ).run_analysis()


def test_a_flow_given_for_a_below_base_year_is_named_as_ignored():
    res = Bulletin17B(BASE, YEARS, below_base_years=[1950]).run_analysis()
    text = next(s for s in res.steps if "not used" in s)
    assert "1950" in text and "below" in text
    res = Bulletin17B(BASE[1:], YEARS[1:], below_base_years=[1950]).run_analysis()
    assert not any("not used" in s for s in res.steps)


# --- bulletin17b_kwargs --------------------------------------------------------------------


def codes_for(**by_year):
    """Peak codes for YEARS, with ``y1955="7"``-style overrides."""
    out = [None] * len(YEARS)
    for key, code in by_year.items():
        out[YEARS.index(int(key[1:]))] = code
    return out


def test_uncoded_record_passes_through():
    kw = bulletin17b_kwargs(BASE, YEARS, None)
    assert kw["peak_flows"] == BASE and kw["water_years"] == YEARS
    assert kw["historic_peaks"] == [] and kw["historic_period"] is None
    assert kw["below_base_years"] == []
    assert set(kw) == {
        "peak_flows",
        "water_years",
        "historic_peaks",
        "historic_period",
        "below_base_years",
        "unobserved_years",
    }
    assert kw["unobserved_years"] == []


def test_code_4_is_below_base():
    kw = bulletin17b_kwargs(BASE, YEARS, codes_for(y1955="4"))
    assert kw["below_base_years"] == [1955]
    assert 1955 not in kw["water_years"] and len(kw["peak_flows"]) == 17


def test_code_4_with_no_printed_flow_is_below_base():
    flows = list(BASE)
    flows[5] = float("nan")
    kw = bulletin17b_kwargs(flows, YEARS, codes_for(y1955="4"))
    assert kw["below_base_years"] == [1955]


@pytest.mark.parametrize("code", ["3", "O", "6", "C", "2,6", "3,4"])
def test_codes_that_remove_the_peak(code):
    kw = bulletin17b_kwargs(BASE, YEARS, codes_for(y1955=code))
    assert 1955 not in kw["water_years"] and kw["below_base_years"] == []
    assert len(kw["peak_flows"]) == 17


def test_a_removed_peak_needs_no_flow():
    flows = list(BASE)
    flows[5] = float("nan")
    kw = bulletin17b_kwargs(flows, YEARS, codes_for(y1955="6"))
    assert len(kw["peak_flows"]) == 17


@pytest.mark.parametrize("code", ["8", "7,8", "7,4", "4,8"])
def test_codes_without_a_b17b_treatment_refuse(code):
    with pytest.raises(Bulletin17BError, match="no Bulletin 17B treatment"):
        bulletin17b_kwargs(BASE, YEARS, codes_for(y1955=code))


def _historic_record():
    """Systematic 1950-1967 and two historic peaks, years 1900 and 1925."""
    years = [1900, 1925] + YEARS
    flows = [3000.0, 2600.0] + BASE
    codes = ["7", "7"] + [None] * 18
    return flows, years, codes


def test_code_7_with_threshold_row_gives_the_period():
    flows, years, codes = _historic_record()
    row = {"start_year": 1890, "end_year": 1949, "threshold_cfs": 2500.0}
    kw = bulletin17b_kwargs(flows, years, codes, [row])
    assert kw["historic_peaks"] == [(1900, 3000.0), (1925, 2600.0)]
    assert kw["historic_period"] == (1890, 1967)
    assert kw["water_years"] == YEARS and kw["peak_flows"] == BASE


def test_code_7_without_threshold_starts_at_first_historic_year():
    flows, years, codes = _historic_record()
    kw = bulletin17b_kwargs(flows, years, codes)
    assert kw["historic_period"] == (1900, 1967)


def test_two_historic_periods_refuse():
    flows, years, codes = _historic_record()
    rows = [
        {"start_year": 1890, "end_year": 1910, "threshold_cfs": 2500.0},
        {"start_year": 1911, "end_year": 1949, "threshold_cfs": 2500.0},
    ]
    with pytest.raises(Bulletin17BError, match="one historic period"):
        bulletin17b_kwargs(flows, years, codes, rows)


def test_a_threshold_that_rises_inside_the_record_refuses():
    flows, years, codes = _historic_record()
    rows = [
        {"start_year": 1890, "end_year": 1949, "threshold_cfs": 2500.0},
        {"start_year": 1960, "end_year": 1967, "threshold_cfs": 1200.0},
    ]
    with pytest.raises(Bulletin17BError, match="inside the systematic record"):
        bulletin17b_kwargs(flows, years, codes, rows)


def test_a_threshold_that_changes_nothing_inside_the_record_is_ignored():
    flows, years, codes = _historic_record()
    rows = [
        {"start_year": 1890, "end_year": 1949, "threshold_cfs": 2500.0},
        {"start_year": 1950, "end_year": 1967, "threshold_cfs": 0.0},
        {"start_year": 1960, "end_year": 1967, "threshold_cfs": 500.0},
    ]
    kw = bulletin17b_kwargs(flows, years, codes, rows)
    assert kw["historic_period"] == (1890, 1967)


def test_a_row_running_into_the_record_at_a_high_threshold_refuses():
    flows, years, codes = _historic_record()
    row = {"start_year": 1890, "end_year": 1967, "threshold_cfs": 2500.0}
    with pytest.raises(Bulletin17BError, match="inside the systematic record"):
        bulletin17b_kwargs(flows, years, codes, [row])


def test_a_gap_between_the_historic_row_and_the_record_refuses():
    flows, years, codes = _historic_record()
    row = {"start_year": 1890, "end_year": 1940, "threshold_cfs": 2500.0}
    with pytest.raises(Bulletin17BError, match="1941 to 1949"):
        bulletin17b_kwargs(flows, years, codes, [row])


def test_a_historic_peak_before_the_row_refuses():
    flows, years, codes = _historic_record()
    row = {"start_year": 1910, "end_year": 1949, "threshold_cfs": 2500.0}
    with pytest.raises(Bulletin17BError, match="1900"):
        bulletin17b_kwargs(flows, years, codes, [row])


def test_a_row_with_no_historic_peaks_still_gives_the_period():
    row = {"start_year": 1930, "end_year": 1949, "threshold_cfs": 5000.0}
    kw = bulletin17b_kwargs(BASE, YEARS, None, [row])
    assert kw["historic_peaks"] == [] and kw["historic_period"] == (1930, 1967)


def test_a_code_7_peak_in_a_gap_of_the_record_is_accepted():
    years = [y for y in YEARS if y != 1955]
    flows = [q for y, q in zip(YEARS, BASE) if y != 1955]
    kw = bulletin17b_kwargs([2900.0] + flows, [1955] + years, ["7"] + [None] * 17)
    assert kw["historic_peaks"] == [(1955, 2900.0)]
    assert kw["historic_period"] == (1950, 1967)
    Bulletin17B(**kw)


def test_a_code_7_peak_in_a_year_that_has_a_systematic_peak_refuses():
    with pytest.raises(Bulletin17BError, match="more than once.*1955"):
        bulletin17b_kwargs(BASE + [2900.0], YEARS + [1955], [None] * 18 + ["7"])


def test_no_peaks_left_refuses():
    with pytest.raises(Bulletin17BError, match="no peaks"):
        bulletin17b_kwargs([1.0, 2.0], [1950, 1951], ["6", "3"])
    with pytest.raises(Bulletin17BError, match="systematic"):
        bulletin17b_kwargs([1.0, 2.0], [1950, 1951], ["7", "7"])


def test_duplicate_years_after_code_handling_refuse_but_a_removed_one_does_not():
    with pytest.raises(Bulletin17BError, match="more than once.*1950"):
        bulletin17b_kwargs([1.0, 2.0], [1950, 1950], None)
    kw = bulletin17b_kwargs(BASE + [9999.0], YEARS + [1950], [None] * 18 + ["6"])
    assert kw["water_years"] == YEARS


@pytest.mark.parametrize(
    "rows",
    [
        [{"start_year": 1890, "end_year": 1949}],
        ["a row"],
        [{"start_year": 1890.5, "end_year": 1949, "threshold_cfs": 1.0}],
        [{"start_year": 1950, "end_year": 1940, "threshold_cfs": 1.0}],
        [{"start_year": 1890, "end_year": 1949, "threshold_cfs": -1.0}],
        [{"start_year": 1890, "end_year": 1949, "threshold_cfs": "x"}],
    ],
)
def test_bad_threshold_rows_refuse(rows):
    with pytest.raises(Bulletin17BError):
        bulletin17b_kwargs(BASE, YEARS, None, rows)


def test_codes_must_match_the_peaks():
    with pytest.raises(Bulletin17BError, match="peak codes"):
        bulletin17b_kwargs(BASE, YEARS, ["7"])
    with pytest.raises(Bulletin17BError, match="peak flows"):
        bulletin17b_kwargs(BASE[:5], YEARS, None)


def test_bad_years_and_negative_flows_refuse():
    with pytest.raises(Bulletin17BError):
        bulletin17b_kwargs(BASE, ["1950"] * 18, None)
    with pytest.raises(Bulletin17BError, match="negative"):
        bulletin17b_kwargs([-1.0] + BASE[1:], YEARS, None)


def test_pandas_style_codes_are_read_as_peak_codes_does():
    codes = [float("nan")] * 18
    codes[5] = "4"
    codes[6] = "2,6"
    kw = bulletin17b_kwargs(BASE, YEARS, codes)
    assert kw["below_base_years"] == [1955] and 1956 not in kw["water_years"]


# --- end to end ----------------------------------------------------------------------------


def test_big_sandy_as_nwis_gives_task_4s_answers():
    """Code 7 on 1897, 1919 and 1927, with one threshold row, runs as Appendix 6 does."""
    bs = record("appendix6_big_sandy")
    years = [y for y, _ in bs.historic_peaks] + list(bs.water_years)
    flows = [q for _, q in bs.historic_peaks] + list(bs.peaks_cfs)
    codes = ["7"] * 3 + [None] * len(bs.water_years)
    row = {"start_year": 1897, "end_year": 1929, "threshold_cfs": 18500.0}
    kw = bulletin17b_kwargs(flows, years, codes, [row])
    assert kw["historic_period"] == bs.historic_period == (1897, 1973)
    res = Bulletin17B(**kw, generalized_skew=-0.2).run_analysis()
    h = res.historic
    assert (h.h, h.z, h.l) == (77, 3, 0) and h.w == pytest.approx(1.68182, abs=1e-5)
    assert res.mean_log == pytest.approx(3.71581, abs=1e-5)
    assert res.std_log == pytest.approx(0.28898, abs=1e-5)
    assert res.skew_weighted == pytest.approx(-0.00409, abs=1.5e-4)
    # The same answers as handing the class the record's own historic information.
    direct = Bulletin17B(
        bs.peaks_cfs,
        bs.water_years,
        historic_peaks=bs.historic_peaks,
        historic_period=bs.historic_period,
        generalized_skew=-0.2,
    ).run_analysis()
    assert res.quantiles["flow_cfs"].tolist() == direct.quantiles["flow_cfs"].tolist()


def test_below_base_code_runs_end_to_end():
    kw = bulletin17b_kwargs(BASE, YEARS, codes_for(y1955="4"))
    res = Bulletin17B(**kw).run_analysis()
    assert res.n_years == 18 and res.n_below_base == 1


# --- fix round 1 ---------------------------------------------------------------------------


def test_a_positive_row_over_unmeasured_years_refuses_when_the_codes_are_off():
    """Big Sandy with its code 7 peaks fitted as systematic: the row would vanish."""
    bs = record("appendix6_big_sandy")
    years = [y for y, _ in bs.historic_peaks] + list(bs.water_years)
    flows = [q for _, q in bs.historic_peaks] + list(bs.peaks_cfs)
    row = {"start_year": 1897, "end_year": 1929, "threshold_cfs": 18500.0}
    with pytest.raises(Bulletin17BError, match="no measured peak"):
        bulletin17b_kwargs(flows, years, None, [row])


def test_a_row_at_its_covered_minimum_over_a_gap_refuses():
    years = list(range(1985, 1988)) + list(range(2003, 2005))
    flows = [900.0, 1100.0, 1000.0, 1200.0, 950.0]
    row = {"start_year": 1985, "end_year": 2004, "threshold_cfs": 900.0}
    with pytest.raises(Bulletin17BError, match="1988 to 2002|no measured peak"):
        bulletin17b_kwargs(flows, years, None, [row])


def test_a_zero_row_is_not_a_historic_row():
    flows, years, codes = _historic_record()
    row = {"start_year": 1800, "end_year": 1949, "threshold_cfs": 0.0}
    kw = bulletin17b_kwargs(flows, years, codes, [row])
    assert kw["historic_period"] == (1900, 1967)


def test_the_peakfq_shape_of_two_rows_gives_the_period():
    flows, years, codes = _historic_record()
    rows = [
        {"start_year": 1890, "end_year": 2025, "threshold_cfs": 0.0},
        {"start_year": 1890, "end_year": 1949, "threshold_cfs": 2500.0},
    ]
    assert bulletin17b_kwargs(flows, years, codes, rows)["historic_period"] == (1890, 1967)


def test_a_duplicated_row_counts_once():
    flows, years, codes = _historic_record()
    row = {"start_year": 1890, "end_year": 1949, "threshold_cfs": 2500.0}
    kw = bulletin17b_kwargs(flows, years, codes, [row, dict(row)])
    assert kw["historic_period"] == (1890, 1967)


def test_removed_peaks_lists_year_flow_and_codes_as_given(caplog):
    codes = codes_for(y1955="2,6", y1956="O", y1957="3")
    flows = list(BASE)
    flows[8] = float("nan")
    with caplog.at_level(logging.INFO, logger="flowfreq.bulletin17b"):
        got = removed_peaks(flows, YEARS, codes)
        bulletin17b_kwargs(BASE, YEARS, codes)
    assert [(y, c) for y, _, c in got] == [(1955, "2,6"), (1956, "O"), (1957, "3")]
    assert got[0][1] == BASE[5]
    removal = [r for r in caplog.records if "removed" in r.getMessage()]
    assert removal and all(str(y) in removal[0].getMessage() for y in (1955, 1956, 1957))
    assert removed_peaks(BASE, YEARS, None) == []


def test_a_historic_peak_after_the_record_extends_the_period():
    flows, years, codes = _historic_record()
    kw = bulletin17b_kwargs(flows + [3100.0], years + [1990], codes + ["7"])
    assert kw["historic_period"] == (1900, 1990)
    assert kw["historic_peaks"][-1] == (1990, 3100.0)
    row = {"start_year": 1890, "end_year": 1949, "threshold_cfs": 2500.0}
    kw = bulletin17b_kwargs(flows + [3100.0], years + [1990], codes + ["7"], [row])
    assert kw["historic_period"] == (1890, 1990)
    Bulletin17B(**kw)


def test_a_kept_peak_with_no_flow_is_dropped_and_logged(caplog):
    flows = list(BASE)
    flows[3] = float("nan")
    with caplog.at_level(logging.WARNING, logger="flowfreq.bulletin17b"):
        kw = bulletin17b_kwargs(flows, YEARS, None)
    assert 1953 not in kw["water_years"] and len(kw["peak_flows"]) == 17
    assert any("1953" in r.getMessage() and "no flow" in r.getMessage() for r in caplog.records)


def test_a_historic_peak_with_no_flow_is_dropped_and_logged(caplog):
    flows, years, codes = _historic_record()
    flows[0] = float("nan")
    with caplog.at_level(logging.WARNING, logger="flowfreq.bulletin17b"):
        kw = bulletin17b_kwargs(flows, years, codes)
    assert kw["historic_peaks"] == [(1925, 2600.0)] and kw["historic_period"] == (1925, 1967)
    assert any("1900" in r.getMessage() for r in caplog.records)


def test_the_gap_check_starts_at_the_record_before_removal():
    codes = codes_for(y1950="6")
    row = {"start_year": 1890, "end_year": 1949, "threshold_cfs": 2500.0}
    kw = bulletin17b_kwargs(BASE, YEARS, codes, [row])
    assert kw["historic_period"] == (1890, 1967) and 1950 not in kw["water_years"]


def test_integer_codes_are_read_like_strings():
    flows, years, _ = _historic_record()
    codes = [np.int64(7), 7] + [None] * 18
    kw = bulletin17b_kwargs(flows, years, codes)
    assert [y for y, _ in kw["historic_peaks"]] == [1900, 1925]
    kw = bulletin17b_kwargs(BASE, YEARS, [4] + [None] * 17)
    assert kw["below_base_years"] == [1950]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"peak_codes": pd.DataFrame({"c": ["7"] * 18})},
        {"peak_codes": "7"},
        {"peak_codes": 7},
        {"peak_codes": [{"a": 1}] + [None] * 17},
        {"peak_codes": None, "perception_thresholds": pd.DataFrame({"start_year": [1890]})},
        {"peak_codes": None, "perception_thresholds": 5},
        {"peak_codes": None, "perception_thresholds": "rows"},
    ],
)
def test_wrong_types_raise_bulletin17b_error(kwargs):
    with pytest.raises(Bulletin17BError):
        bulletin17b_kwargs(BASE, YEARS, **kwargs)


def test_an_infinite_threshold_over_peaks_refuses():
    row = {"start_year": 1950, "end_year": 1960, "threshold_cfs": float("inf")}
    with pytest.raises(Bulletin17BError, match="unobserved"):
        bulletin17b_kwargs(BASE, YEARS, None, [row])


def test_an_infinite_threshold_over_unmeasured_years_outside_the_period_is_ignored(caplog):
    flows, years, codes = _historic_record()
    rows = [
        {"start_year": 1890, "end_year": 1949, "threshold_cfs": 2500.0},
        {"start_year": 1800, "end_year": 1850, "threshold_cfs": float("inf")},
    ]
    with caplog.at_level(logging.INFO, logger="flowfreq.bulletin17b"):
        kw = bulletin17b_kwargs(flows, years, codes, rows)
    assert kw["historic_period"] == (1890, 1967)
    assert any("1800" in r.getMessage() and "no record" in r.getMessage() for r in caplog.records)


def test_a_non_integer_threshold_year_says_whole_numbers():
    row = {"start_year": 1890.5, "end_year": 1949, "threshold_cfs": 1.0}
    with pytest.raises(Bulletin17BError, match="threshold start years must be whole numbers"):
        bulletin17b_kwargs(BASE, YEARS, None, [row])


# --- fix round 2 ---------------------------------------------------------------------------


def test_a_base_row_over_below_base_years_is_accepted():
    codes = codes_for(y1955="4", y1956="4")
    row = {"start_year": 1950, "end_year": 1967, "threshold_cfs": 721.0}
    kw = bulletin17b_kwargs(BASE, YEARS, codes, [row])
    assert kw["below_base_years"] == [1955, 1956] and kw["historic_period"] is None
    high = {"start_year": 1950, "end_year": 1967, "threshold_cfs": 722.0}
    with pytest.raises(Bulletin17BError, match="change of threshold"):
        bulletin17b_kwargs(BASE, YEARS, codes, [high])


def _after_record():
    """Systematic 1950-1967 and a code 7 peak in 1975."""
    return BASE + [3000.0], YEARS + [1975], [None] * 18 + ["7"]


@pytest.mark.parametrize("end", [1974, 1975, 1980])
def test_a_row_after_the_record_to_a_historic_peak_is_the_historic_row(end):
    flows, years, codes = _after_record()
    row = {"start_year": 1968, "end_year": end, "threshold_cfs": 2500.0}
    kw = bulletin17b_kwargs(flows, years, codes, [row])
    assert kw["historic_period"] == (1950, 1975)
    assert kw["historic_peaks"] == [(1975, 3000.0)]
    Bulletin17B(**kw)


def test_a_row_before_and_a_row_after_refuse_as_two_historic_periods():
    flows, years, codes = _after_record()
    rows = [
        {"start_year": 1940, "end_year": 1949, "threshold_cfs": 2500.0},
        {"start_year": 1968, "end_year": 1975, "threshold_cfs": 2500.0},
    ]
    with pytest.raises(Bulletin17BError, match="more than one historic period"):
        bulletin17b_kwargs(flows, years, codes, rows)


def test_a_row_after_the_record_with_no_historic_peak_refuses_with_its_reason():
    row = {"start_year": 1968, "end_year": 1975, "threshold_cfs": 2500.0}
    with pytest.raises(Bulletin17BError, match="reaches no historic peak"):
        bulletin17b_kwargs(BASE, YEARS, None, [row])


def test_a_gap_between_the_record_and_a_row_after_it_refuses():
    flows, years, codes = _after_record()
    row = {"start_year": 1971, "end_year": 1975, "threshold_cfs": 2500.0}
    with pytest.raises(Bulletin17BError, match="1968 to 1970"):
        bulletin17b_kwargs(flows, years, codes, [row])


def _gapped(drop):
    years = [y for y in YEARS if y not in drop]
    return [q for y, q in zip(YEARS, BASE) if y not in drop], years


def test_infinite_rows_over_empty_years_override_a_base_row(caplog):
    flows, years = _gapped({1955, 1960})
    base = {"start_year": 1950, "end_year": 1967, "threshold_cfs": 700.0}
    rows = [
        base,
        {"start_year": 1955, "end_year": 1955, "threshold_cfs": float("inf")},
        {"start_year": 1960, "end_year": 1960, "threshold_cfs": float("inf")},
    ]
    with caplog.at_level(logging.INFO, logger="flowfreq.bulletin17b"):
        kw = bulletin17b_kwargs(flows, years, None, rows)
    assert kw["water_years"] == years and kw["historic_period"] is None
    assert any("1955" in r.getMessage() for r in caplog.records)
    with pytest.raises(Bulletin17BError, match="no measured peak"):
        bulletin17b_kwargs(flows, years, None, [base])


def test_an_infinite_row_inside_the_historic_period_is_accepted():
    flows, years, codes = _historic_record()
    rows = [
        {"start_year": 1890, "end_year": 1949, "threshold_cfs": 2500.0},
        {"start_year": 1930, "end_year": 1940, "threshold_cfs": float("inf")},
    ]
    kw = bulletin17b_kwargs(flows, years, codes, rows)
    assert kw["historic_period"] == (1890, 1967)
    assert kw["unobserved_years"] == list(range(1930, 1941))
    res = Bulletin17B(**kw).run_analysis()
    assert res.historic is not None and res.historic.h == 78 - 11


def test_overlapping_infinite_rows_count_a_year_once():
    flows, years, codes = _historic_record()
    inf = float("inf")
    rows = [
        {"start_year": 1890, "end_year": 1949, "threshold_cfs": 2500.0},
        {"start_year": 1930, "end_year": 1940, "threshold_cfs": inf},
        {"start_year": 1938, "end_year": 1942, "threshold_cfs": inf},
    ]
    kw = bulletin17b_kwargs(flows, years, codes, rows)
    assert kw["unobserved_years"] == list(range(1930, 1943))


def test_an_infinite_row_over_a_kept_peak_still_refuses():
    rows = [{"start_year": 1960, "end_year": 1961, "threshold_cfs": float("inf")}]
    with pytest.raises(Bulletin17BError, match="unobserved|infinite"):
        bulletin17b_kwargs(BASE, YEARS, None, rows)


def test_only_removed_years_next_to_the_record_close_the_gap():
    years = [1906] + list(range(1930, 1948))
    flows = [900.0] + BASE
    codes = ["6"] + [None] * 18
    row = {"start_year": 1897, "end_year": 1905, "threshold_cfs": 2500.0}
    with pytest.raises(Bulletin17BError, match="covered by no period"):
        bulletin17b_kwargs(flows, years, codes, [row])
    years = [1929] + list(range(1930, 1948))
    row = {"start_year": 1897, "end_year": 1928, "threshold_cfs": 2500.0}
    kw = bulletin17b_kwargs(flows, years, codes, [row])
    assert kw["historic_period"] == (1897, 1947)


def test_unmeasured_years_are_listed_as_years_and_runs():
    flows, years = _gapped({1955, 1960, 1961})
    row = {"start_year": 1950, "end_year": 1967, "threshold_cfs": 700.0}
    with pytest.raises(Bulletin17BError) as err:
        bulletin17b_kwargs(flows, years, None, [row])
    assert "1955, 1960 to 1961" in str(err.value) and "and others" not in str(err.value)
