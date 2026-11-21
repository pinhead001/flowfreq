# Bulletin 17B procedure, as printed (Revised Sept 1981, editorial corrections March 1982)

Each item gives the bulletin page label as `[p. x]`. To convert a label to a PDF page, see README.md: main text page n = PDF page n+8.
Equations use LaTeX-style plain text. Logs are base 10 throughout. Where the bulletin's own wording matters,
it is quoted.

> **Section numbering note.** In this edition the weighting of skew is **Section V.B.4** (V.B.3 is
> "Estimating Generalized Skew"). Outliers are **V.B.9**, historic flood data **V.B.10**, incomplete
> record **V.B.6**, zero flood years **V.B.7**. Confidence limits, risk and expected probability are
> **Section VI.A, VI.B, VI.C**. There is no Section V.C.5: V.C is "Refinements to Frequency Curve"
> (comparisons with similar watersheds, precipitation estimates).

---

## 1. Distribution and fitting (V.B.1-2) [pp. 9-10]

Log-Pearson Type III fitted by moments of the base-10 logarithms.

- Eq. 1 [p. 9]: `\log Q = \bar X + K S`, where K is the Pearson III frequency factor for the skew and exceedance probability P (Appendix 3).
- Eq. 2 [p. 10]: `\bar X = \sum X / N`
- Eq. 3a: `S = [ \sum (X-\bar X)^2 / (N-1) ]^{0.5}`
- Eq. 3b: `S = [ ( \sum X^2 - (\sum X)^2/N ) / (N-1) ]^{0.5}`
- Eq. 4a: `G = N \sum (X-\bar X)^3 / [ (N-1)(N-2) S^3 ]`
- Eq. 4b: `G = [ N^2 \sum X^3 - 3N \sum X \sum X^2 + 2 (\sum X)^3 ] / [ N(N-1)(N-2) S^3 ]`

X is the log of the annual peak, N the number of items. The bulletin warns that 3b/4b lose more precision than 3a/4a with limited significant digits [p. 10].

## 2. Generalized skew (V.B.3) [pp. 10-12]

A regional study needs at least 40 stations (or all within 100 miles) with 25 or more years of record. Three methods are compared: an isoline map, a prediction equation, and the mean of station skews (about 20 stations where hydrology is heterogeneous). The method with the smaller MSE is used if that MSE is "significantly smaller" than the variance of the station skews; otherwise the mean is used, with the variance as `MSE_{\bar G}` [p. 12]. Without a detailed study, `\bar G` is read from Plate I.

## 3. Weighted skew (V.B.4) [pp. 12-15]

- Eq. 5 [p. 12]: `G_W = [ MSE_{\bar G} \cdot G + MSE_G \cdot \bar G ] / [ MSE_{\bar G} + MSE_G ]`
  - G = station skew, `\bar G` = generalized skew, `MSE_{\bar G}` = MSE of the generalized skew, `MSE_G` = MSE of the station skew.
- **Default generalized-skew MSE:** "When generalized skews are read from Plate I, the value of `MSE_{\bar G}` = 0.302 should be used in equation 5." [p. 13]
- Eq. 6 [p. 13], from Wallis, Matalas and Slack:
  `MSE_G \approx 10^{[A - B \log_{10}(N/10)]}`
  - `A = -0.33 + 0.08|G|` if `|G| \le 0.90`; `A = -0.52 + 0.30|G|` if `|G| > 0.90`
  - `B = 0.94 - 0.26|G|` if `|G| \le 1.50`; `B = 0.55` if `|G| > 1.50`
  - |G| is the absolute station skew (used as the population skew estimate), and N is the record length in years.
  - "If the historic adjustment described in Appendix 6 has been applied, the historically adjusted skew, `\tilde G`, and historic period, H, are to be used for G and N, respectively, in equation 6." [p. 13]
- Table 1 [p. 14] tabulates eq. 6 for |G| = 0.0-3.0 and N = 10-100 (`table1_skew_mse.csv`). **Its |G| = 0.9 row uses the `|G| > 0.90` branch of A** (0.562 at N = 10, i.e. 10^-0.25), contrary to the `<= 0.90` in eq. 6, which would give 0.552. All other rows match eq. 6 as printed.
- Caution [p. 15]: with |station skew| > 2 and long records, eq. 6 gives the station value little weight. Eq. 5 may weight the generalized skew improperly if the generalized and station skews differ by more than 0.5. In both cases the analyst should examine the data and possibly give more weight to the station skew.
- **In the Appendix 12 examples the weighted skew is always rounded to the nearest tenth** before K values are read from Appendix 3 (Ex. 1: 0.6678 to 0.7; Ex. 2: 0.0745 to 0.1; Ex. 3: 0.5590 to 0.6; Ex. 4: -0.4487 to -0.4).

## 4. Broken record (V.B.5) and incomplete record (V.B.6) [p. 15]

- A broken record is analyzed as continuous, with length equal to the sum of the segments, unless the watershed changed physically between segments.
- Incomplete record: missing high peaks should be estimated (usually by the data agency). Where a crest gage's bottom is not reached in some years, "use of the conditional probability adjustment is recommended as described in Appendix 5."

## 5. Zero flood years (V.B.7) [p. 15]

The log of zero is minus infinity, so "The conditional probability adjustment is recommended for determining frequency curves for records with zero flood years as described in Appendix 5."

## 6. Outliers (V.B.9, Appendix 4, flow chart p. 12-3) [pp. 17-19]

**Which test first** [p. 17]:
- "If the station skew is greater than +0.4, tests for high outliers are considered first."
- "If the station skew is less than -0.4 tests for low outliers are considered first."
- "Where the station skew is between ± 0.4, tests for both high and low outliers should be applied before eliminating any outliers from the data set." (Ex. 2 [p. 12-15]: "the tests for both high outliers and low outliers are based on the systematic record statistics before any adjustments are made.")

**High outliers**, eq. 7 [p. 17]: `X_H = \bar X + K_N S`
- `\bar X` = "mean logarithm of systematic peaks (X's) excluding zero flood events, peaks below gage base, and outliers previously detected"; S = standard deviation of the X's; `K_N` = Appendix 4 value for sample size N.
- Peaks with log > `X_H` are high outliers. They are compared with historic flood data and with flood information at nearby sites. "If information is available which indicated a high outlier(s) is the maximum in an extended period of time, the outlier(s) is treated as historic flood data as described in Section V.B.10. If useful historic information is not available to adjust for high outliers, then they should be retained as part of the systematic record." [p. 17]

**Low outliers**, eq. 8a [p. 18]: `X_L = \bar X - K_N S` (terms as in eq. 7).
- Eq. 8b, used if a historic adjustment has already been made: `X_L = \tilde M - K_H \tilde S`. Here `K_H` = Appendix 4 value "for period used to compute `\tilde M` and `\tilde S`", and `\tilde M`, `\tilde S` are the historically adjusted mean and standard deviation.
- "Flood peaks considered low outliers are deleted from the record and the conditional probability adjustment described in Appendix 5 is applied." [p. 18]
- The test is applied **once**. All values above the eq. 7 threshold or below the eq. 8a/8b threshold are outliers. It is equivalent to a one-sided 10 percent significance test with normal-distribution K values. If several values sit very close to a threshold, sensitivity testing may be desirable. [p. 18]

**Appendix 4** [p. 4-1]: one-sided 10 percent significance K_N for N = 10 to 149 (see `kn_table.csv`).

## 7. Historic flood data (V.B.10, Appendix 6) [p. 19, pp. 6-1 to 6-3]

Historic data should be used unless a comparison of the systematic and historically adjusted analyses, the magnitude of the observed peaks, or other factors suggest otherwise. The adjustment assumes that the systematic record represents the intervening (H-Z) years.

Symbols [p. 6-2]: X = logs of systematic peaks **excluding zero flood events, peaks below base, high or low outliers**; `X_z` = log of a historic peak, including a high outlier that has historic information; N = number of X's; Z = number of historic peaks including high outliers with historic information; H = years in the historic period; L = "number of low values to be excluded, such as: number of zeros, number of incomplete record years (below measurable base), and low outliers which have been identified"; M, S, G = mean, std, skew of the X's; E = event rank, 1 to (Z+N); a = plotting-position constant (Weibull 0, Beard 0.3, Hazen 0.5).

Equations [p. 6-3]:
- 6-1: `W = (H - Z) / (N + L)` (systematic record weight)
- 6-2a: `\tilde M = [ W \sum X + \sum X_z ] / (H - WL)`
- 6-3a: `\tilde S^2 = [ W \sum (X-\tilde M)^2 + \sum (X_z-\tilde M)^2 ] / (H - WL - 1)`
- 6-4a: `\tilde G = \frac{H - WL}{(H-WL-1)(H-WL-2)} \cdot \frac{ W \sum (X-\tilde M)^3 + \sum (X_z-\tilde M)^3 }{ \tilde S^3 }`
- 6-2b (from previously computed statistics): `\tilde M = [ W N M + \sum X_z ] / (H - WL)`
- 6-3b: `\tilde S^2 = [ W(N-1)S^2 + WN(M-\tilde M)^2 + \sum (X_z-\tilde M)^2 ] / (H - WL - 1)`
- 6-4b: `\tilde G = \frac{H - WL}{(H-WL-1)(H-WL-2)\tilde S^3} [ \frac{W(N-1)(N-2)S^3 G}{N} + 3W(N-1)(M-\tilde M)S^2 + WN(M-\tilde M)^3 + \sum (X_z-\tilde M)^3 ]`
- 6-5: `\log Q = \tilde M + K \tilde S`
- 6-6: `\tilde m = E` for `1 \le E \le Z`
- 6-7: `\tilde m = W E - (W-1)(Z + 0.5)` for `(Z+1) \le E \le (Z+N+L)`
- 6-8: `\tilde{PP} = \frac{\tilde m - a}{H + 1 - 2a} \cdot 100` (percent)

Notes:
- Appendix 6 step 4 [p. 6-1] says its example "does not include the final step of weighting with the generalized skew. The historically adjusted skew developed by this procedure is appropriate to use in developing a generalized skew." Weighting then uses `\tilde G` and H in eq. 6 (Section 3 above; Ex. 2 uses H = 82).
- Plotting positions (VII.B) [pp. 26-27]: no formula is required for the calculations. The Weibull position `P = m/(N+1)` (a = 0) is used in the examples and the computer program.
- Figures 6-1 and 6-2 [pp. 6-4 to 6-6] give a worked historic example, transcribed in `examples/appendix6_big_sandy.json`.
  - Big Sandy River at Bruceton, TN, 3-6065. Record 1897, 1919, 1927, 1930-1973; historical period 1897-1973. N = 44, Z = 3, H = 77, W = 1.68182.
  - M~ = 3.71581, S~ = 0.28898, G~ = 0.0418.
  - Eq. 6 is applied with N = H = 77, giving MSE 0.07074. Weighted with generalized skew -0.2, G_W = -0.00409.
  - Unlike Appendix 12, G_W is **not** rounded to a tenth: the K values are for -0.00409. Figure 6-1 thus does include the weighting step, despite step 4's text.

## 8. Conditional probability adjustment and synthetic statistics (Appendix 5) [pp. 5-1 to 5-3]

Applies when the record is truncated by peaks below a gage base, zero years, and/or low outliers. "These procedures should only be used when not over 25 percent of the total record has been truncated." The truncation level is the minimum discharge that excludes the below-base peaks, zero flows and all low outliers, and no other discharges. Before the adjustment, statistics are computed for the above-truncation peaks, outlier detection and historic information are completed, and all flow-chart steps except the last are done [p. 5-1].

1. Eq. 5-1a: `\tilde P = N / n` (N = peaks above the truncation level, n = total years of record). With historic information, eq. 5-1b instead: `\tilde P = (H - WL) / H`, where L = number of peaks truncated and W is from eq. 6-1 [p. 5-2].
2. Eq. 5-2: `P = \tilde P \times P_d` (P_d = exceedance probabilities on the conditional curve, which is computed from the above-truncation statistics).
3. Interpolate (graphically or mathematically) the adjusted curve at exceedance probabilities 0.01, 0.10 and 0.50.
4. Synthetic statistics:
   - 5-3: `G_s = -2.50 + 3.12 \frac{\log(Q_{.01}/Q_{.10})}{\log(Q_{.10}/Q_{.50})}`
   - 5-4: `S_s = \frac{\log(Q_{.01}/Q_{.50})}{K_{.01} - K_{.50}}`
   - 5-5: `\bar X_s = \log(Q_{.50}) - K_{.50} S_s`
   - K.01 and K.50 are the Pearson III deviates for skew `G_s` [p. 5-3]. Eq. 5-3 is an approximation for skews between +2.5 and -2.0.
5. Compare the synthetic-statistics curve with the observed peaks, using plotting positions on the total years n or H.
6. "The minimum additional requirement to arrive at a final frequency curve is the determination of the weighted skew." In the examples `G_s` is weighted with the generalized skew by eq. 5, and the final curve is `\log Q = \bar X_s + K_{G_W} S_s`.

## 9. Confidence limits (Appendix 9) [pp. 9-1 to 9-5]

- 9-3a/b: `U_{P,c}(X) = \bar X + S K^U_{P,c}`, `L_{P,c}(X) = \bar X + S K^L_{P,c}`. `\bar X` and S are the mean and std "of the final estimated log Pearson Type III frequency curve".
- 9-4a: `K^U_{P,c} = \frac{K_{G_W,P} + \sqrt{K_{G_W,P}^2 - ab}}{a}`
- 9-4b: `K^L_{P,c} = \frac{K_{G_W,P} - \sqrt{K_{G_W,P}^2 - ab}}{a}`
- 9-5: `a = 1 - \frac{z_c^2}{2(N-1)}`
- 9-6: `b = K_{G_W,P}^2 - \frac{z_c^2}{N}`
- `z_c` = standard normal deviate with cumulative probability c (exceedance 1-c). One-sided limits at level c. The union of two one-sided 95 percent limits is a two-sided 90 percent interval (eq. 9-2).
- "The systematic record length N is deemed to control the statistical reliability of the estimated frequency curve and is to be used for calculating confidence limits even when historic information has been used to estimate the frequency curve." [p. 9-3]
- Worked check [pp. 9-3 to 9-4]: mean 3.00, S 0.25, skew 0.20, N 50, P 0.01. z_c = 1.645, K = 2.4723, a = 0.9724, b = 6.058, K^U = 3.026, K^L = 2.059, U = 3.756, L = 3.515. Natural units 3270 and 5700 cfs; Q.01 = 4150 cfs.
- Table 9-1 (`table9_1_confidence_deviates.csv`) gives the exact zero-skew non-central-t deviates. Its "confidence level" is the exceedance probability of the deviate, so .05 gives K^U and .95 gives K^L for a 95 percent one-sided limit.
- Table 9-1 worked check [pp. 9-4 to 9-5]: mean 3.00, S 0.25, skew 0.00, N 50. K^U = 2.862, K^L = 1.936, U = 3.715, L = 3.484. Natural units 3050 and 5190 cfs; Q.01 = 3820 cfs.

## 10. Expected probability (VI.C, Appendix 11) [pp. 24-25, 11-3 to 11-6]

- The expected-probability adjustment, if made, is the final step and must be documented. Confidence limits, risk and weighting use the curve without it [p. 25].
- 11-1: `P_N = \text{Prob}[ t_{N-1} > K_n (N/(N+1))^{1/2} ]` (normal-theory). K_n = standard normal variate, N = sample size.
- Table 11-1 [p. 11-6] gives P_N against P_inf, entered with N-1. For P > 0.50, P_N = 1 minus the tabulated value for 1-P.
- 11-2a-f [p. 11-4]: P_N for P = .0001: `.0001(1.0 + 1600/N^{1.72})`; .001: `.001(1.0 + 280/N^{1.55})`; .01: `.01(1.0 + 26/N^{1.16})`; .05: `.05(1.0 + 6/N^{1.04})`; .10: `.1(1.0 + 3/N^{1.04})`; .30: `.3(1.0 + 0.46/N^{0.925})`.
- Check [p. 11-4]: P = .01 and N = 20 gives `(.01)(1.0 + 26/32.3)` = 0.018; Table 11-1 gives 0.0174.

## 11. Order of steps (Appendix 12 flow diagrams) [pp. 12-1 to 12-3]

Outline [p. 12-1]: A. determine available data (previous studies, gage records, historic data, similar watersheds, watershed model). B. evaluate data (homogeneity, reliability). C. compute the curve per the flow diagrams.

Main flow diagram [p. 12-2]:
- Zero flood or incomplete record: see Appendix 5 (conditional probability adjustment); for outliers see pp. 17-19 and Appendices 5 and 6.
- Complete record: compute station statistics, then decide on a 2-station comparison. Footnote: "IF SYSTEMATIC RECORD LENGTH IS LESS THAN 50 YEARS THE ANALYST SHOULD CONSIDER WHETHER THE USE OF THE PROCEDURES OF APPENDIX 7 IS APPROPRIATE."
  - Yes: compute the extended record (App. 7), then test and adjust for high outliers / historic information (App. 6).
  - No: test and adjust for outliers / historic information (p. 12-3 chart), then determine the weighted G ("SEE PAGES 12 TO 15", i.e. main-text pp. 12-15).
- Then a decision point (note: is further analysis warranted?), then seek more historic flood data.
- If P < 0.01, do the flood estimate from precipitation (pp. 21-22). If N <= 24, also the precipitation estimate. If N <= 50, comparisons with similar watersheds (App. 8).
- Final curve. Optionally confidence limits (App. 9) and the expected-probability adjustment (App. 11).

Historic and outlier flow chart [p. 12-3]. Start by testing the station skew G:
- **G > 0.4:** test for historic peaks/high outliers. If yes, recompute statistics adjusted for historic peaks/high outliers (App. 6). Then, either way, test for low outliers. If yes, recompute statistics without low outliers, then do the conditional probability adjustment (App. 5). If no, return.
- **G < -0.4:** test for low outliers. If yes, recompute statistics without low outliers. Then test for historic peaks/high outliers. If yes, recompute adjusted for historic peaks/high outliers (App. 6). Then "low outliers omitted?" If yes, do the conditional probability adjustment (App. 5). If no, return.
- **-0.4 <= G <= 0.4:** test for historic peaks/high outliers.
  - If no, go to the low-outlier test of the G > 0.4 branch.
  - If yes, test for low outliers. If yes, recompute statistics without low outliers, then recompute adjusted for historic peaks/high outliers (App. 6). If no, go directly to the App. 6 recomputation. Then "low outliers omitted?" If yes, do the conditional probability adjustment. Then return to p. 12-2.

---

## Ambiguities and points to note

1. **Section labels.** As noted at the top, the weighting section is V.B.4 and the outlier section is V.B.9. Treatment of low outliers is in V.B.9 and Appendix 5; there is no V.C.5.
2. **Record length N in eq. 6 is not uniform across the examples.**
   - Historic case: the text says H, and Ex. 2 uses H = 82, with MSE 0.073 obtained by interpolating Table 1.
   - Conditional (low-outlier) case: Ex. 3 uses n = 38 (total years, including the deleted low outlier) and G_s = 0.5948, giving MSE 0.183.
   - Zero-flow case: Ex. 4 gives MSE 0.163 "From Table I" without stating N or the skew used. Table 1 near |G| = 0.5 at N around 40 is consistent, but this is not stated.
   - The examples take MSE_G from Table 1 interpolation (3 decimals), not from eq. 6 directly. Ex. 1 does evaluate eq. 6: A = -.2716, B = .7502, MSE_G = 0.277.
3. **Weighted skew is rounded to the nearest tenth in every example** before K is read. An implementation that interpolates K at the unrounded G_W will not reproduce the printed tables to the printed precision.
4. **Which skew K.01/K.50 use in eqs. 5-4/5-5.** Appendix 5 says "skew coefficient G_s". Ex. 3 uses the K values for the rounded *conditional* skew 0.6 (2.75514, -.09945), and its G_s is 0.5948. Ex. 4 uses K values for G_s rounded to -0.5 (1.95472, .08302). In Ex. 3, 0.5948 rounds to 0.6 anyway, so both readings agree to the tenth.
5. **Unrounded G_s in the weighting.** Ex. 3 weights G_s = 0.5948. Ex. 4 weights -.529 (three decimals of -0.5287), not the rounded -0.5.
6. **Conditional skew is rounded to the nearest tenth** ("for ease in computation") before the conditional curve: Ex. 3 0.6311 to 0.6; Ex. 4 -0.4396 to -0.4.
7. **Ordering when |G| <= 0.4.** The flow chart says that with both high outliers (historic) and low outliers present, low outliers are removed first and then the historic adjustment is made. Eq. 8b is for the case "If an adjustment for historic flood data has previously been made". Its K_H is the Appendix 4 value "for period used to compute M~ and S~", which is ambiguous: H, or N? Appendix 4 stops at N = 149.
8. **High outliers without historic information** stay in the systematic record. The flow-chart boxes say "test for historic peaks/high outliers"; a high outlier triggers the App. 6 adjustment only when historic information exists.
9. **L in Appendix 6** counts zeros, below-base years and low outliers. W = (H-Z)/(N+L), and eq. 6-7 runs E up to Z+N+L. So truncated values take ranks in the plotting positions, while N (the number of X's) excludes them.
10. **Eq. 12-22 (Ex. 2)** prints the factor as `H / [(H-WL-1)(H-WL-2) S~^3]`, whereas eq. 6-4b has `(H-WL)`. They are equal there because L = 0.
11. **Confidence-level labelling is inconsistent.**
    - Ex. 1 text [p. 12-10]: "confidence limits for levels of significance of .05 and .95 percent". Both column heads of Table 12-4 read "0.05 UPPER LIMIT CURVE" and "0.05 LOWER LIMIT CURVE". Figure 12-1's legend says ".05 and .95 Confidence Limit Curves". z_c = 1.64485, i.e. one-sided c = 0.95, a two-sided 90 percent band.
    - App. 9 p. 9-5 also says table 9-1 is entered "with confidence level 0.05" for K^U and "0.95" for K^L.
    - N in eq. 9-5/9-6 is the systematic record length (Ex. 1: N = 24). No example shows which N to use for conditional or zero-flow records.
12. **Expected probability, N vs N-1.** Table 11-1 is entered with N-1 (Ex. 1 uses N-1 = 23 for N = 24). Eq. 11-2 uses N, the sample size: the p. 11-4 check uses N = 20 (20^1.16 = 32.3) and compares with Table 11-1 at N-1 = 19.
13. **Ex. 1 skew denominator** uses S = .24561 (S^3 = .014816), while the reported S is 0.2456.
14. **Eq. 5-1a vs 5-1b.** Ex. 4 uses P~ = 35/42, counting both zero years and the low outlier as truncated.

Misprints in the source (kept as printed in the JSON):
- Table 12-10 [p. 12-33], Pd = .05: Q printed **98900**, with log Q 3.9952. It should be about 9890.
- Ex. 4 text [p. 12-36] says "the adjusted frequency derived in Step 4"; Figure 12-8 says "Adjusted Frequency Curve at Step 5". The adjustment is in Step 5.
- Ex. 3 note [p. 12-30]: "producers" for "procedures".
