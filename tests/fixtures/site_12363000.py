"""USGS 12363000 annual peaks, WY 1922-2023.

A parity fixture with a shape none of the others have: a **gap**. Water years
1924-1927 carry no peak, so the record is 98 observations spanning
102 years. Big Sandy, Powder River and Cains Coulee are all
contiguous, and ``wymt_case`` refuses a site that is not.

That gap is the reason the case sets ``fill_missing_years=False``. The default
fill in ``build_emafit_inputs`` censors every unobserved year in range against
the lower perception threshold, which is right when a threshold has been
declared -- somebody asserting a peak that size would have been recorded. No
such assertion exists here: the years are simply unmeasured. Filling them
anyway is not a rounding difference, it changes the answer, at-site skew
+0.435 -> +0.250 and Q100 120,064 -> 119,473 cfs.
``Bulletin17C`` given peaks and water years fits the 98 observations, so the
Fortran must be handed the same 98 rows for the comparison to mean anything.

The record is also strongly right-skewed: the 176,000 cfs peak in
1964 is 1.7x the next largest, which drives the at-site skew to
+0.435 -- well clear of the 0.04 floor at ``emafit.f:763``, so ``detrat`` is
reached here. MGBT finds no low outliers, so unlike Cains Coulee nothing is
censored by the fit either: this is the clean, uncensored, strongly-skewed,
gapped case.

Source: NWIS peak-flow export, retrieved 2026-09. Flows are cfs, keyed by USGS
water year. The station name is deliberately not recorded -- the export carried
only the numbers, and inventing one would put an unverified claim in a fixture
whose whole purpose is to be trustworthy.
"""

#: Annual peak flows in cfs, keyed by water year. No historic peaks; every
#: value is a systematic observation.
SYSTEMATIC_PEAKS = {
    1922: 82_200,
    1923: 88_000,
    1928: 101_000,
    1929: 69_700,
    1930: 38_800,
    1931: 62_300,
    1932: 89_800,
    1933: 91_200,
    1934: 60_200,
    1935: 71_000,
    1936: 71_800,
    1937: 46_000,
    1938: 70_400,
    1939: 65_600,
    1940: 41_300,
    1941: 26_400,
    1942: 56_300,
    1943: 62_800,
    1944: 34_700,
    1945: 52_700,
    1946: 68_400,
    1947: 83_700,
    1948: 102_000,
    1949: 65_200,
    1950: 74_600,
    1951: 69_000,
    1952: 47_200,
    1953: 48_900,
    1954: 69_600,
    1955: 42_100,
    1956: 66_200,
    1957: 50_500,
    1958: 44_600,
    1959: 58_000,
    1960: 46_400,
    1961: 58_300,
    1962: 34_300,
    1963: 29_600,
    1964: 176_000,
    1965: 45_700,
    1966: 42_300,
    1967: 59_400,
    1968: 38_900,
    1969: 39_900,
    1970: 43_900,
    1971: 51_100,
    1972: 59_400,
    1973: 36_600,
    1974: 65_900,
    1975: 77_600,
    1976: 46_400,
    1977: 19_700,
    1978: 37_400,
    1979: 43_200,
    1980: 37_900,
    1981: 40_100,
    1982: 41_000,
    1983: 36_700,
    1984: 34_400,
    1985: 39_900,
    1986: 41_800,
    1987: 37_700,
    1988: 25_600,
    1989: 37_300,
    1990: 47_700,
    1991: 59_500,
    1992: 27_400,
    1993: 40_400,
    1994: 32_100,
    1995: 66_000,
    1996: 49_900,
    1997: 60_800,
    1998: 27_900,
    1999: 44_000,
    2000: 34_000,
    2001: 24_700,
    2002: 45_400,
    2003: 39_700,
    2004: 26_400,
    2005: 37_600,
    2006: 47_000,
    2007: 37_900,
    2008: 53_500,
    2009: 40_500,
    2010: 31_100,
    2011: 49_200,
    2012: 53_100,
    2013: 54_300,
    2014: 50_400,
    2015: 30_700,
    2016: 29_600,
    2017: 47_000,
    2018: 50_200,
    2019: 34_300,
    2020: 50_200,
    2021: 43_600,
    2022: 53_400,
    2023: 35_500,
}

#: Water years inside the record's span that carry no observation.
MISSING_YEARS = [1924, 1925, 1926, 1927]

#: Nationwide B17C generalized skew and its standard error, the app's defaults
#: and what the parity run used.
REGIONAL_SKEW = -0.302
REGIONAL_SKEW_SD = 0.55
