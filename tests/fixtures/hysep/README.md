# French Creek near Phoenixville, Pa. (USGS 01472157)

`01472157_daily.csv`: daily mean flow (cfs), 1968-12-01 to 1970-02-28 and 1992-03-01 to 1992-05-31, from the USGS Water
Data API, retrieved 2026-11-22.

It is the station USGS WRIR 96-4040 (Sloto and Crouse, 1996) works through: drainage area 59.1 mi2, interval 2N* of 5
days. The report's text and figures 2 to 4 (April 1992) and its Appendix 2 (monthly summary, local-minimum method, 1969)
are the expected values in `tests/test_regime.py`.

**The Appendix 2 listing runs 11 days early.** Moving the published daily flows back 11 days reproduces all 12 of its
monthly mean streamflows to the 0.01 cfs printed, and no other offset from -12 to 12 days comes within 1.8 cfs of every
month. The test moves the dates the same way. The report does not mention it.
