"""The transcribed Bulletin 17B fixtures load and reproduce the bulletin's own sums."""

from __future__ import annotations

import csv

import numpy as np
import pytest

from tests.fixtures.b17b import FIXTURES, record


def test_records_have_the_printed_lengths():
    assert len(record("example1").peaks_cfs) == 24
    assert len(record("example2").peaks_cfs) == 39
    assert len(record("example3").peaks_cfs) == 38
    ex4 = record("example4")
    assert len(ex4.peaks_cfs) == 42 and sum(q == 0 for q in ex4.peaks_cfs) == 6
    bs = record("appendix6_big_sandy")
    assert len(bs.peaks_cfs) == 44 and bs.water_years[0] == 1930 and bs.water_years[-1] == 1973
    assert bs.historic_peaks == [(1897, 25000.0), (1919, 21000.0), (1927, 18500.0)]
    assert bs.historic_period == (1897, 1973)


def test_example1_sums_reproduce():  # Table 12-2: sum X = 80.84043
    x = np.log10(record("example1").peaks_cfs)
    assert x.sum() == pytest.approx(80.84043, abs=5e-5)


def test_kn_table_is_complete():
    with open(FIXTURES / "kn_table.csv", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))[1:]
    assert len(rows) == 140 and rows[0] == ["10", "2.036"] and rows[-1] == ["149", "3.148"]
