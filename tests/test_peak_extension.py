"""Record extension of annual peaks (Bulletin 17C Appendix 8), on the Bulletin's own example.

Suwanee Creek at Suwanee, Georgia (20 years, 1985 to 2004) extended from the Etowah River
at Canton (1892 to 2004). Table 8-1 is the concurrent record, table 8-3 the Etowah River
before it and table 8-2 the 13 values the Bulletin generates, by its equation 8-25.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from flowfreq.peak_extension import PeakExtension, extend_peaks, variance_terms

#: Table 8-1: water year, Etowah River, Suwanee Creek (cfs).
CONCURRENT = [
    (1985, 5030, 1440),
    (1986, 3090, 386),
    (1987, 12200, 2150),
    (1988, 9340, 948),
    (1989, 9080, 1220),
    (1990, 27100, 3760),
    (1991, 5940, 1320),
    (1992, 7660, 696),
    (1993, 10900, 2540),
    (1994, 9420, 1190),
    (1995, 10500, 2650),
    (1996, 19500, 4350),
    (1997, 11300, 2360),
    (1998, 15000, 2900),
    (1999, 5530, 816),
    (2000, 8900, 862),
    (2001, 9270, 2090),
    (2002, 7100, 1260),
    (2003, 13600, 2940),
    (2004, 15300, 3270),
]
#: Table 8-3: the Etowah River, 1892 to 1984, in order (cfs).
ETOWAH_BEFORE = [
    36700, 9380, 4360, 12400, 5300, 10100, 10900, 20800, 12300, 19600, 20800, 20400, 9620,
    8870, 12900, 6950, 10600, 14900, 6880, 11100, 18400, 8350, 7440, 9300, 36100, 14300, 6320,
    29500, 36100, 21500, 10200, 11700, 4960, 8740, 6530, 9460, 10700, 13600, 15300, 5200,
    10400, 28200, 9600, 6530, 22700, 15300, 19700, 6360, 8900, 8820, 13300, 10100, 10600,
    5180, 32300, 14500, 8500, 17200, 8500, 7790, 19500, 8140, 15500, 12600, 7300, 15500, 5440,
    7230, 6320, 19300, 20900, 22600, 25000, 8740, 19000, 15900, 11000, 11900, 6590, 5790,
    14200, 11200, 12300, 11900, 19000, 18500, 17200, 21600, 14900, 5500, 24100, 8970, 12800,
]  # fmt: skip
#: Table 8-2: the extended record, 1972 to 1984 (cfs).
EXTENDED = [2830, 2010, 2300, 2200, 4310, 4150, 3730, 5180, 3040, 722, 6070, 1460, 2440]


@pytest.fixture(scope="module")
def records():
    assert len(ETOWAH_BEFORE) == 93
    etowah = pd.Series(
        ETOWAH_BEFORE + [e for _, e, _ in CONCURRENT],
        index=list(range(1892, 1985)) + [y for y, _, _ in CONCURRENT],
        dtype=float,
    )
    suwanee = pd.Series([s for _, _, s in CONCURRENT], index=[y for y, _, _ in CONCURRENT])
    return suwanee.astype(float), etowah


@pytest.fixture(scope="module")
def extension(records) -> PeakExtension:
    return extend_peaks(*records)


def test_sample_statistics_are_the_bulletins(extension):
    assert (extension.n_concurrent, extension.n_additional) == (20, 93)
    assert extension.correlation == pytest.approx(0.8519, abs=5e-5)
    assert extension.mean_short == pytest.approx(3.215, abs=5e-4)
    assert extension.mean_index_concurrent == pytest.approx(3.984, abs=5e-4)
    assert extension.mean_index_additional == pytest.approx(4.079, abs=5e-4)
    assert extension.sd_short == pytest.approx(0.279, abs=5e-4)
    assert extension.sd_index_concurrent == pytest.approx(0.214, abs=5e-4)
    assert extension.sd_index_additional == pytest.approx(0.219, abs=5e-4)


def test_augmented_moments_and_length_are_the_bulletins(extension):
    assert extension.mean_augmented == pytest.approx(3.3024, abs=5e-5)
    assert extension.variance_augmented == pytest.approx(0.08105, abs=5e-6)
    assert extension.n_extended == 13
    assert round(extension.effective_length_variance - 20) == 13
    assert extension.effective_length_mean > 20


def test_the_move_line_is_equation_8_25(extension):
    assert extension.intercept == pytest.approx(3.4364, abs=5e-5)
    assert extension.slope == pytest.approx(1.4416, abs=5e-5)
    assert extension.mean_index_extended == pytest.approx(4.1414, abs=5e-5)


def test_the_extended_record_is_table_8_2(extension):
    table = extension.extended
    assert list(table["water_year"]) == list(range(1972, 1985))
    for value, printed in zip(table["flow_cfs"], EXTENDED):
        # Three significant figures, as printed.
        assert value == pytest.approx(printed, rel=6e-3)
    assert extension.predict([5500.0])[0] == pytest.approx(722, rel=6e-3)


def test_the_extended_record_has_the_augmented_moments(records, extension):
    """What equations 8-23 and 8-24 are for: n1 observed + ne generated reproduce 8-7 and 8-8."""
    suwanee, _ = records
    logs = np.log10(np.concatenate([suwanee.to_numpy(), extension.extended["flow_cfs"]]))
    assert logs.mean() == pytest.approx(extension.mean_augmented, abs=1e-12)
    assert logs.var(ddof=1) == pytest.approx(extension.variance_augmented, abs=1e-12)


def test_years_and_length_can_be_chosen(records):
    chosen = extend_peaks(*records, years=list(range(1960, 1973)))
    assert list(chosen.extended["water_year"]) == list(range(1960, 1973))
    assert chosen.mean_augmented == pytest.approx(3.3024, abs=5e-5)  # the augmentation is unchanged
    fewer = extend_peaks(*records, n_extend=8)
    assert list(fewer.extended["water_year"]) == list(range(1977, 1985))


def test_short_only_years_are_counted_not_fitted(records):
    suwanee, etowah = records
    longer = pd.concat([suwanee, pd.Series({2005: 1500.0, 2006: 900.0})])
    ext = extend_peaks(longer, etowah)
    assert ext.n_short_only == 2 and ext.n_concurrent == 20
    assert ext.intercept == pytest.approx(3.4364, abs=5e-5)


def test_variance_terms_are_finite_and_symmetric_in_nothing():
    a, b, c = variance_terms(20, 93)
    assert all(np.isfinite(v) for v in (a, b, c))
    assert (a, b, c) != variance_terms(21, 93)


def test_a_weak_relation_is_noted_or_refused(records):
    suwanee, etowah = records
    rng = np.random.default_rng(3)
    noisy = suwanee * 10 ** rng.normal(0.0, 0.35, len(suwanee))
    try:
        ext = extend_peaks(noisy, etowah)
    except ValueError as exc:
        assert "too little to extend" in str(exc) or "not possible" in str(exc)
    else:
        assert any("not above the 0.80" in note for note in ext.notes)


@pytest.mark.parametrize(
    "change,said",
    [
        (lambda s, e: (s.iloc[:8], e), "at least 10"),
        (lambda s, e: (s, e.loc[1985:]), "fewer than 2 years outside"),
        (lambda s, e: (s, e, {"years": [1990, 1991]}), "not index-gage years outside"),
        (lambda s, e: (s, e, {"n_extend": 200}), "more than the 93 index years"),
        (lambda s, e: (s, e, {"years": [1980, 1981], "n_extend": 3}), "n_extend is 3"),
        (lambda s, e: (pd.Series(1000.0, index=s.index), e), "no variance"),
    ],
)
def test_refusals(records, change, said):
    args = change(*records)
    kwargs = args[2] if len(args) == 3 else {}
    with pytest.raises(ValueError, match=said):
        extend_peaks(args[0], args[1], **kwargs)


def test_not_a_series_is_a_type_error(records):
    with pytest.raises(TypeError, match="pandas Series"):
        extend_peaks([1.0, 2.0], records[1])


def test_concurrent_relation_sets_the_two_parts_side_by_side(records):
    from flowfreq.peak_extension import concurrent_relation

    suwanee, etowah = records
    table = concurrent_relation(suwanee, etowah)
    assert list(table["period"]) == ["1985-1994", "1995-2004", "1985-2004"]
    assert list(table["n"]) == [10, 10, 20]
    assert table.loc[2, "correlation"] == pytest.approx(0.8519, abs=5e-5)
    ratios = suwanee / etowah.loc[suwanee.index]
    # The median is taken of the log10 ratio, so an even count gives a geometric middle.
    assert table.loc[2, "median_ratio"] == pytest.approx(
        float(10 ** np.median(np.log10(ratios))), rel=1e-9
    )
    chosen = concurrent_relation(suwanee, etowah, split_year=1990)
    assert list(chosen["period"]) == ["1985-1989", "1990-2004", "1985-2004"]


def test_concurrent_relation_sees_a_relation_that_moved(records):
    from flowfreq.peak_extension import concurrent_relation

    suwanee, etowah = records
    shifted = suwanee.copy()
    shifted.loc[1995:] *= 2.0
    table = concurrent_relation(shifted, etowah, split_year=1995)
    assert table.loc[1, "median_ratio"] > 1.5 * table.loc[0, "median_ratio"]
    with pytest.raises(ValueError, match="leaves one part"):
        concurrent_relation(suwanee, etowah, split_year=1985)
    with pytest.raises(ValueError, match="fewer than 2 concurrent"):
        concurrent_relation(suwanee.iloc[:1], etowah.loc[:1984])
