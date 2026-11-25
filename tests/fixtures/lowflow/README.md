# Low flows against USGS SIR 2026-5021

`sir2026_5021_table1.csv` is three rows of Table 1 of the data release of USGS SIR 2026-5021 (Feaster and others, 2026;
doi:10.5066/P13BHRVE), and `<site>_daily.csv.gz` the daily mean flows of those gages from the USGS Water Data API,
retrieved 2026-11-24. The report's method is this module's: log-Pearson III, by station skew, of the annual minimum n-day
means of April-March years, with a conditional probability for zeros.

Of 32 unregulated gages compared on 2026-11-24 (20 to 92 years; 9 with zero years), the number of years agrees at all, and the moments
of the 7-day minima and the 1Q10, 7Q2, 7Q10, 30Q2 and 30Q3 at 31, to the figures printed (three values at three of those
are 0.005 off, on a rounding tie). With windows that could cross April 1, as before 0.29.0, it was 29. Notes:

- 03487500 (here): the lowest 7 days of 1956 that end in it end on April 1, with six of them in the year before. The module
  credited them to 1956 until 0.29.0 (7Q10 498, 504 printed); windows are now kept inside each year.
- 02330450 (not here): the 30Q2 and 30Q3 were 0.6 and 0.3 percent low for the same reason.
- 02065500 (not here), the one that does not agree: the printed moments of the 7-day minima need about 20 cfs in 2002, where the approved record has
  0.02 cfs for seven days in August. The printed 1Q10 agrees with the record.

Flow duration agrees within 1 percent at q95, q50 and q10 at every gage, at all but one at q90 (1.1 percent) and q1 (1.4),
and at all but three at q99 (1.5 to 7 percent, flows under 1 cfs at two). Not traced.
