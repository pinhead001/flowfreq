# Bulletin 17B procedure for comparing with Bulletin 17C — design

Date: 2026-11-21. Status: for the owner's approval (builder check done). TODO.md, review of
2026-11-19, item 7.

## Intent

Many effective flood studies and older reports used Bulletin 17B (1981/1982), usually through
PeakFQ 5.x. Bulletin 17C replaced it in 2018. When flowfreq's B17C answer for a gage differs from
the number in such a study, the owner needs to show how much of the difference comes from the
change in procedure, as opposed to added years of record or a new regional skew. So flowfreq gets
a faithful Bulletin 17B analysis that runs on the same peaks as its Bulletin 17C analysis, and a
report can put the two side by side.

The old `B17CEngine` (removed in #405) was not this. It did only B17B's first step, a moments fit
with station skew. `MethodOfMoments` in `bulletin17c.py` is closer: it has a Grubbs-Beck
low-outlier test and a conditional-probability adjustment. But it uses B17C's skew variance, it
has no synthetic statistics, no historic adjustment and no high-outlier test, and it leaves zero
flows out of the adjustment. `Bulletin17C(method=MOM)` stays as it is, because its numbers are
already in use.

## Goals

- A `Bulletin17B` class in a new module `flowfreq/bulletin17b.py` that runs the Bulletin 17B
  procedure as printed (Revised September 1981, editorial corrections March 1982), in the order
  of the Appendix 12 flow chart (p. 12-3):
  - moments of the base-10 logs (eqs. 2, 3a, 4a);
  - outlier tests with the Appendix 4 K_N table, in B17B's order: high first when the station
    skew is above +0.4, low first below -0.4, and both on the unadjusted systematic statistics
    in between;
  - high outliers: treated as historic peaks when a historic period is given, otherwise kept in
    the systematic record;
  - the historic adjustment (Appendix 6, eqs. 6-1 to 6-4a, with L counting zeros, below-base
    years and low outliers), and the low-outlier test against the adjusted statistics (eq. 8b)
    when it was made;
  - zero flows, below-base years and low outliers handled by the conditional-probability
    adjustment with synthetic statistics (Appendix 5, eqs. 5-1a/5-1b to 5-5), refused when more
    than 25 percent of the record is truncated;
  - the station skew weighted with a generalized skew by eq. 5, with the station skew's MSE
    from eq. 6 (Wallis, Matalas and Slack; the historically adjusted skew and H when the
    historic adjustment was made). The generalized skew's MSE defaults to 0.302 (Plate I);
  - quantiles by `log Q = X̄ + K S` with the exact Pearson III K (`core.kfactor`);
  - Appendix 9 confidence limits (eqs. 9-4/9-5/9-6, systematic N even with historic data),
    90 percent two-sided by default as in the examples (z = 1.64485 each side);
  - the expected-probability adjustment (Appendix 11, eq. 11-1) as an extra column, never
    applied to the quantiles themselves.
- A record of each step the procedure took (which test ran first, which peaks were outliers,
  whether the historic adjustment and the conditional adjustment were made, each set of
  statistics), so a report can say why B17B and B17C differ.
- A `b17b_ffa` report section kind, accepted by the `comparison` section, so a report shows the
  B17B and B17C quantiles of the same gage in one table and on one curve plot.
- Tests against every worked example the bulletin prints: Appendix 12 examples 1-4, the
  Appendix 6 historic example (Big Sandy River at Bruceton, TN, this project's primary site),
  and the Appendix 9 and Appendix 11 hand checks.

## Non-goals

- Appendix 7 (two-station record extension; flowfreq has `extend_peaks` for B17C), Appendix 8
  (comparisons with similar watersheds), precipitation-based estimates, the mixed-population
  and regulated-flow guidance. They are judgement steps outside the computation.
- Digitising Plate I. The caller supplies the generalized skew; the docstring says where it
  came from in B17B practice.
- Reproducing PeakFQ 5.x output digit for digit. No PeakFQ 5 output is available to test
  against, and vendored peakfq 8.1 has no B17B mode. If the owner finds PeakFQ 5 output for a
  gage, it becomes another test, and any difference is investigated then.
- A chat-tool (MCP) entry point. It can follow once the class and the section exist.

## Anti-goals

- Changing any `Bulletin17C` result, including `method=MOM`.
- A B17B result that could be mistaken for current guidance. Every table, figure caption and
  section title that shows one says "Bulletin 17B (superseded)".
- Silently falling back to something else when B17B gives no answer (more than 25 percent
  truncated, a synthetic skew outside the eq. 5-3 range of -2.0 to +2.5, N above the K_N
  table). These raise with the reason, as `MethodOfMoments` does.

## Constraints

- Python 3.11-3.14, mypy-clean, black/isort, NumPy docstrings, no `print`, logging only
  (CLAUDE.md conventions).
- Inputs in the same shapes as `Bulletin17C`: peak flows with water years, the historic period
  and historic peaks, a user low-outlier threshold. A peak record with NWIS codes reaches it
  through the same `peak_code_kwargs` path where that is meaningful (code 7, historic).
- The worked-example data (transcribed from the scanned bulletin, every value as printed) live
  under `tests/fixtures/b17b/` as JSON, with the bulletin page of each number.

## Approaches

1. **A new `Bulletin17B` class (chosen by the owner, 2026-11-21).** New module, built from
   `MethodOfMoments`' fitting code, with everything B17B adds. `method=MOM` is untouched. Costs
   some duplication of the moments and conditional-probability code, which can be shared
   through small helper functions rather than inheritance.
2. A `procedure="17B"` switch on `MethodOfMoments`. Less new code, but one class with two rule
   books is hard to verify, and a forgotten flag gives B17C-flavoured numbers labelled B17B.
3. Turn `MethodOfMoments` into B17B. Simplest, but it silently changes every existing
   `method=MOM` result.

## Decisions the owner made

- Approach 1, a new class (2026-11-21).
- Validate against the Bulletin 17B worked examples (2026-11-21: "show what changed in method";
  PeakFQ output if any turns up).
- **The generalized skew on the B17B side (2026-11-21, option 3).** By default a `b17b_ffa`
  section uses the same regional skew and MSE as the B17C side (`regional_skew`, or the
  `regional_skew_state` / `regional_skew_from_site` lookups `gaged_ffa` has), so the difference
  is procedure only, and the section says so in a line. A `generalized_skew` with
  `generalized_skew_mse` (default 0.302, Plate I) overrides it with the old study's value. A
  `comparison` can then carry three columns, B17C, B17B with today's skew and B17B with the old
  study's skew, splitting the difference into procedure and skew. With no skew given and no
  lookup, B17B uses station skew only and says so, as `gaged_ffa` does.

## Decisions taken as the designer's calls (the owner may override)

- **Exact by default, hand-computation mode for the tests.** The bulletin's examples round the
  weighted skew to the nearest tenth before reading K, read the station skew's MSE off Table 1,
  and round the conditional skew to a tenth. PeakFQ and any program compute these exactly, and
  a comparison with B17C should not carry hand-computation rounding. So the default is exact
  (eq. 6 evaluated, K at the unrounded skew). A `hand_computation=True` option reproduces the
  examples' roundings so the tests can match the printed numbers to the printed precision; the
  exact mode is also tested against the same examples with a stated tolerance. The Appendix 6
  Big Sandy example supports exact as the default: it does not round, evaluating K at the
  weighted skew -0.00409 (eq. 6 at H = 77 gives MSE 0.07074; adjusted M̃ 3.71581, S̃ 0.28898,
  G̃ 0.0418; generalized skew -0.2 at 0.302; 1% flood 24,391 cfs). It is tested in exact mode.
  Table 1's |G| = 0.9 row was computed with the "|G| > 0.90" form of A, against eq. 6 as printed;
  eq. 6 as printed is what the code evaluates.
- Eq. 6 uses N = systematic years for the plain and low-outlier cases, H with the historically
  adjusted skew after the historic adjustment (the bulletin says so, p. 13), and the total
  years n for a conditional record (Examples 3 and 4 do this: eq. 6 at three decimals gives their
  0.183 at n = 38 and 0.163 at n = 42, and Example 2's 0.073 at H = 82).
- Eqs. 5-4/5-5 take K at the synthetic skew G_s, as Appendix 5 says. In hand-computation mode,
  the examples' own choices are followed (Example 3 uses the conditional skew's K).
- Eq. 8b's K_H is read for the systematic N, not H (the table stops at 149 and H is often
  larger). Open point 1 below.
- Confidence limits use the final curve's mean and standard deviation (synthetic ones after a
  conditional adjustment) with the systematic N, per p. 9-3.
- Expected probability uses eq. 11-1 with N (Table 11-1 is the same formula tabulated at N-1;
  the hand checks are tested against the table).
- Plotting positions: Weibull, with Appendix 6's historically weighted positions (eqs. 6-6 to
  6-8) when the historic adjustment was made, for the figure only.

- **Historic peaks (Z, H).** B17B takes one historic period of H years with its historic peaks.
  Systematic peaks above the smallest historic peak move into Z, as B17B practice and PeakFQ 5
  did; high outliers move into Z only when a historic period is given. From B17C inputs: a
  single historic period (one perception threshold before the systematic record, with the
  historic peaks inside it) translates directly; more than one historic period, or thresholds
  that change within the systematic record, are refused with the reason.
- **Censored and B17C-only peak types.** Zero flows and NWIS code 4 (less-than, i.e. below the
  gage base) peaks are truncated (counted in L, and in the 25 percent limit). Code 8 interval
  peaks and historic interval peaks have no B17B treatment and are refused. The
  regulated-record refusal applies as in `gaged_ffa`.
- **Record end year.** No new option: to compare like with like, give the B17B section the
  peaks the old study had (truncated `peaks`, or `site_no` with an `end_year` if `gaged_ffa`
  takes one by then).
- **Refusal inside a report.** The class raises. A `b17b_ffa` section, alone or in a
  `comparison`, catches that refusal and renders "Bulletin 17B gives no answer for this record:
  <reason>" in place of its table, so the B17C side and the rest of the report still build.
- **The comparison plot** is drawn only when a `b17b_ffa` estimate is in the comparison, so
  existing comparison reports are unchanged.
- **Eq. 8b's K_H** is read at the systematic N unless the Appendix 6 example shows otherwise.
  Example 2 cannot settle it: its station skew is within ±0.4, so eq. 8b never runs there.

## Open points for the builder (minor)

1. Eq. 8b's K_H: test the Appendix 6 Big Sandy example and Example 2 both ways and keep the one
   that reproduces them; if neither discriminates, use N and document it.
2. Whether to tabulate Table 1 or always evaluate eq. 6: evaluate, and keep Table 1 only as a
   test fixture for hand-computation mode.
3. The K_N table: the transcribed Appendix 4 values for N = 10-149, refusing outside that range
   (the existing `core.grubbs_beck_critical_value` approximation is not the bulletin's table).
4. Where the shared helpers live (moments, conditional probability): `bulletin17b.py` or a small
   private module both work.
5. The section's option names follow `gaged_ffa` (`site_no`, `site_name`, `peaks`,
   `regional_skew`, `regional_skew_mse`, `regional_skew_state`, `regional_skew_from_site`) plus
   `generalized_skew`, `generalized_skew_mse`, `historic_period`, `low_outlier_threshold`,
   `hand_computation`.

## Source material

- Bulletin 17B PDF: https://pubs.usgs.gov/unnumbered/70275162/report.pdf (scanned).
- Transcription made on 2026-11-21 (to be committed as test fixtures): the procedure with
  equations and page labels, the Appendix 4 K_N table, Appendix 12 examples 1-4, the Appendix 6
  Big Sandy example, Tables 1, 9-1 and 11-1, and a list of 14 ambiguities in the bulletin
  itself (skew rounding, which N goes into eq. 6, which skew's K goes into the synthetic
  statistics, the confidence-level labels, a misprinted 98900 for 9890 in Table 12-10).
