# Bulletin 17B transcription

Source: `b17b.pdf`, "Guidelines for Determining Flood Flow Frequency, Bulletin #17B of the Hydrology
Subcommittee, Revised September 1981, Editorial Corrections March 1982". 194 PDF pages of scanned images.
The PDF has an OCR text layer (Adobe Paper Capture). That layer was used only to locate pages; every
number here was read from the page images, at 200-220 dpi crops for the tables.

## PDF page to bulletin page map

| PDF pages | Bulletin labels | Content |
|---|---|---|
| 1 | (cover) | Cover |
| 2-6 | i-v | Foreword, Hydrology Subcommittee, page-revision notice |
| 7-8 | vi-vii | Contents |
| 9-36 | 1-28 | Main text (**label n = PDF n+8**). V.B.2 eqs. 1-4 on pp. 9-10 (PDF 17-18); V.B.4 weighting on pp. 12-13 (PDF 20-21); Table 1 on p. 14 (PDF 22); V.B.9 outliers on pp. 17-19 (PDF 25-27); VI on pp. 22-25 (PDF 30-33); VII.B plotting position on pp. 26-27 (PDF 34-35) |
| 37-40 | 1-1 to 1-4 | App. 1 References |
| 41-48 | 2-1 to 2-8 | App. 2 Glossary and Notation |
| 49-75 | 3-1 to 3-27 (approx.) | App. 3 K-value tables (positive skew PDF 50-62, negative skew PDF 63-75) |
| 76 | 4-1 | App. 4 Outlier test K_N |
| 77-80 | 5-1 to 5-4 | App. 5 Conditional probability adjustment |
| 81-87 | 6-1 to 6-7 | App. 6 Historic data (equations p. 6-3 = PDF 83; Big Sandy example Figs. 6-1/6-2 on PDF 84-86; Fig. 6-3 on PDF 87) |
| 88-96 | 7-1 to 7-9 | App. 7 Two-station comparison |
| 97-98 | 8-1 to 8-2 | App. 8 Weighting of independent estimates |
| 99-108 | 9-1 to 9-10 | App. 9 Confidence limits (Table 9-1 on PDF 105-108) |
| 109-112 | 10-1 to 10-4 | App. 10 Risk |
| 113-118 | 11-1 to 11-6 | App. 11 Expected probability (Table 11-1 on PDF 118) |
| 119-155 | 12-1 to 12-37 | App. 12 Flow diagrams and examples (**label 12-k = PDF 118+k**) |
| 156 | 13-1 | App. 13 Computer program |
| 157-191 | 14-1 to 14-35 | App. 14 Flood flow frequency techniques report summary |
| 192-194 | (Plate I) | Generalized skew map and its text |

Appendix 12 detail:

| Bulletin | PDF | Content |
|---|---|---|
| 12-2, 12-3 | 120, 121 | Flow diagrams |
| 12-5 | 123 | Table 12-1, annual peaks for the four example stations |
| 12-6 to 12-14 | 124-132 | Example 1 |
| 12-15 to 12-22 | 133-140 | Example 2 |
| 12-23 to 12-31 | 141-149 | Example 3 |
| 12-32 to 12-37 | 150-155 | Example 4 |

## Files

- `procedure.md`: the B17B procedure with equations and page labels. Covers moments; generalized and weighted skew with eq. 6 (A/B, 0.90 and 1.50 breakpoints, Plate I MSE 0.302); the outlier tests and their order; Appendix 5 conditional probability and synthetic statistics; Appendix 6 historic weighting and plotting positions; Appendix 9 confidence limits; Appendix 11 expected probability; the Appendix 12 flow diagrams. Ends with a list of ambiguities and source misprints.
- `kn_table.csv`: Appendix 4 K_N for N = 10 to 149 (140 rows, every printed row). Columns: `sample_size_N,K_N`.
- `examples/example1.json` to `example4.json`: the Appendix 12 examples. All numeric values from the bulletin are **strings exactly as printed**, so leading-dot forms like ".2456" and trailing zeros survive. Integer counts (N, H, Z, years) are JSON integers. A value the bulletin does not print is `null`. Each file has:
  - `station`
  - `input`: the Table 12-1 peaks as printed, plus historic, zero-year and generalized-skew information
  - `computations`: every printed intermediate, with page labels and the printed equation text
  - the final frequency table, plus the confidence-limit and expected-probability tables in Example 1
  - `uncertain_digits`

| File | Example | Station | Record |
|---|---|---|---|
| example1.json | 1. Fitting the Log-Pearson Type III Distribution | Fishkill Creek at Beacon, NY, 01-3735 | 1945-1968, N = 24 |
| example2.json | 2. Adjusting for a High Outlier | Floyd River at James, IA, 06-6005 | 1935-1973, N = 39 (38 after removing 1953); historic period 1892-1973, H = 82, Z = 1 |
| example3.json | 3. Testing and Adjusting for a Low Outlier | Back Creek nr Jones Springs, WV, 01-6140 | 1929-31 and 1939-1973, N = 38 (37 after removing the 1969 low outlier); 1936 estimate 22000* not used |
| example4.json | 4. Adjusting for Zero Flood Years | Orestimba Creek nr Newman, CA, 11-2745 | 1932-1973, n = 42; 6 zero years; 36 non-zero; 35 after removing the 1955 low outlier |

Stations are printed in the old 6-digit form (e.g. "01-3735"). The modern 8-digit numbers (e.g. 01373500) are not printed and have not been added.

## Verification

- Each table was read twice: once at page resolution and once from a 200-220 dpi crop.
- Transcribed data were also checked for internal consistency. These checks were never used to change a value.
  - Table 12-1 peaks reproduce every printed mean, std and skew to 4 decimals: Example 2 before and after removing 1953; Example 3 before and after removing 1969; Example 4 non-zero and after removing 1955.
  - Table 12-2 X, X^2 and X^3 agree with log10 of the peaks. Their sums equal the printed sums.
  - Every final-table log Q equals mean + K·S within 1 in the 4th decimal.
  - Table 12-4 K^U and K^L agree with the Appendix 9 formulas.
  - Table 12-5 agrees with Table 11-1 at N-1 = 23.
  - K_N agrees with the Grubbs-Beck approximation `-0.9043 + 3.345 sqrt(log N) - 0.4046 log N` within 0.004 at all 140 rows, and is strictly increasing.

## Uncertain digits

**No digit in the transcribed tables or computations was illegible.** All `uncertain_digits` lists are empty.
Things worth knowing:

- Appendix 4, N = 98: printed "3,011" with a comma-like decimal mark. Transcribed as 3.011.
- Appendix 4, N = 43: printed "2.7l0" with a slightly malformed "1". Read as 2.710, which is consistent with neighbouring rows.
- **Misprint, not an illegible digit:** Table 12-10 (p. 12-33), Pd = .05, Q printed **98900** against log Q 3.9952 (about 9890). Kept as printed in `example4.json`, with a `Q_cfs_note` and an entry in `printed_inconsistencies`.
- Ex. 4 p. 12-36 text says "Step 4" where the figure legend says "Step 5". This affects labels only.
- Table 1, |G| = 0.0 and N = 10: the last digit is faint in the dot-matrix print and could be 6 or 8. I read 0.468, which agrees with eq. 6 (10^-0.33 = 0.4677). This is the only Table 1 cell with any doubt.
- Table 9-1, level .05, N = 80, P .005: printed "3:016" with a smudged decimal point. Transcribed as 3.016, which agrees with the exact value.
- Appendix 6, Fig. 6-2, eq. 6-3b: the second "1.68182" has an overstruck digit. Read as 1.68182.
- Appendix 6, Fig. 6-1: 1950's log is printed 3.99475, where the correct rounding is 3.99476. Departures for 1897 and 1919 differ by 0.00001 between Figs. 6-1 and 6-2. These are source rounding, kept as printed.
- Appendix 6, eq. 6-5 table, 2 percent row: (S~)(K) is clearly printed 0.59289. But 0.28898 x 2.05159 = 0.59287, and the printed log Q 4.30868 uses 0.59287. This is a source slip, kept as printed.

## Second batch (added on request)

- `examples/appendix6_big_sandy.json`: the Appendix 6 worked historic example (Figures 6-1 and 6-2, pp. 6-4 to 6-6, PDF 84-86). Big Sandy River at Bruceton, TN, 3-6065, DA 205 sq mi. It contains:
  - the full Figure 6-1 table: year, Q, log, departure, weight, E, weighted order m~, Weibull PP
  - historic peaks and systematic peaks by year
  - every Figure 6-1 continuation sum and result: eqs. 6-2a, 6-3a and 6-4a; eq. 6 MSE; weighted skew; the eq. 6-5 frequency table at 12 exceedance percents; eqs. 6-6 to 6-8
  - all of Figure 6-2: standard statistics, the historic-peaks table, and eqs. 6-1, 6-2b, 6-3b and 6-4b
  - Inputs: N = 44, Z = 3, H = 77. L is not printed; W = (77-3)/44 = 1.68182 implies L = 0. Generalized skew is -0.2 with MSE 0.302.
  - **Historic years:** 1897 (25,000), 1919 (21,000), 1927 (18,500).
  - **Systematic years:** 1930-1973, all 44 present, no gaps. "Record: 1897, 1919, 1927, 1930-1973 (47 years)"; historical period 1897-1973.
- `table1_skew_mse.csv`: Table 1 (p. 14, PDF 22), MSE of station skew for |G| = 0.0-3.0 (step 0.1) and N = 10-100 (step 10). 31 x 10 values.
- `table9_1_confidence_deviates.csv`: Table 9-1 (pp. 9-7 to 9-10, PDF 105-108), confidence-limit deviates for the normal distribution (zero skew). Confidence levels .01, .05, .10, .25, .75, .90, .95, .99; N = 10, 15, 20, 25, 30, 40, ..., 100; 12 exceedance probabilities from .002 to .990. 96 rows x 12 values.
- `table11_1_expected_probability.csv`: Table 11-1 (p. 11-6, PDF 118), P_N against P_inf. N-1 = 1-30, 40, 60, 120 and inf; P_inf = .50, .30, .10, .05, .01, .001, .0001.

### Verification of the second batch

**Big Sandy.**
- log10 of every Q matches the printed log, except 1950 (see below).
- The systematic statistics recompute to M = 3.69094, S = 0.26721 and G = -0.1874 against the printed -0.18746.
- The sums of squared and cubed departures recompute to the printed 3.09755 and -0.37648.

**Table 1.** Every cell matches eq. 6 within 0.001 except the whole |G| = 0.9 row. That row was computed with the |G| > 0.90 branch of A (`-0.52 + 0.30|G|`), e.g. 10^-0.25 = 0.562. Eq. 6 as printed says `A = -0.33 + 0.08|G|` for `|G| <= 0.90`, which would give 0.552 at N = 10. This is a real inconsistency between Table 1 and eq. 6 at the breakpoint; I read the values carefully and they are as printed.

**Table 9-1.** All 1152 values agree with the exact non-central t quantile within 0.0025, except one printed value. The quantile is `nct.isf(level, N-1, z_p sqrt(N))/sqrt(N)` (scipy, run from the project venv read-only).
- The exception is at level .01, N = 10, P .200: printed 2.243, exact 2.246.
- Its mirror image, level .99, N = 10, P .800, is printed -2.243 (exact -2.246).
- Both are clearly printed. This is a source rounding or computation difference.
- **The confidence level in Table 9-1 is the exceedance probability of the non-central t deviate.** So level .05 gives the upper limit and .95 the lower, as used on p. 9-5.

**Table 11-1.** Agrees with eq. 11-1 (exact Student t) within 4% everywhere, except four .0001-column cells that are within 5%: N-1 = 11 (.0021), 12 (.0018), 30 (.00046) and 120 (.00017). These are source rounding, and the printed digits are clear.

### Appendix 3 (not transcribed)

Appendix 3 is the Harter table. The footnote on p. 3-1 cites H. Leon Harter, Technometrics 11(1):177-187, 1969 and 13(1):203-204, 1971. Negative-skew values were obtained by inverting the positive table and changing signs (Central Technical Unit, SCS).

The grid:
- **Skew G:** 0.0 to 9.0 and -0.0 to -9.0, in steps of 0.1. That is 181 columns, 7 per page across PDF 50-75 (pp. 3-2 to 3-27).
- **Exceedance P, 31 rows:** .9999, .9995, .9990, .9980, .9950, .9900, .9800, .9750, .9600, .9500, .9000, .8000, .7000, .6000, .5704, .5000, .4296, .4000, .3000, .2000, .1000, .0500, .0400, .0250, .0200, .0100, .0050, .0020, .0010, .0005, .0001.
- **Precision:** 5 decimals.
- About 5,600 values, so it was skipped as instructed.
- Eq. 3-1, the Wilson-Hilferty approximation `K = (2/G){[(K_n - G/6)(G/6) + 1]^3 - 1}`, is given as an approximation for -1 <= G <= 1. The bulletin says "use of the table is preferred".

## Not transcribed (out of scope)

- Appendix 3 K tables (described above).
- The Appendix 12 figures and Figure 6-3 (plots only).
