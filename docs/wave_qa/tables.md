## Sample: evaluated site-pairs by kind

| pair | border | inland | straddles |
|---|---|---|---|
| AZ-NM | 1 | 1 | 4 |
| AZ-NV | 0 | 0 | 1 |
| AZ-UT | 0 | 2 | 2 |
| CO-NM | 5 | 2 | 4 |
| CO-UT | 0 | 3 | 0 |
| CO-WY | 1 | 0 | 1 |
| ID-MT | 6 | 12 | 8 |
| ID-NV | 1 | 1 | 2 |
| ID-OR | 5 | 6 | 0 |
| ID-UT | 2 | 0 | 1 |
| ID-WA | 3 | 8 | 0 |
| ID-WY | 1 | 0 | 6 |
| MT-ID | 9 | 23 | 13 |
| MT-WY | 4 | 4 | 22 |
| NM-AZ | 5 | 5 | 6 |
| NM-CO | 8 | 1 | 8 |
| NV-AZ | 2 | 0 | 0 |
| NV-ID | 2 | 1 | 1 |
| NV-OR | 0 | 1 | 4 |
| NV-UT | 1 | 1 | 0 |
| OR-ID | 1 | 4 | 0 |
| OR-NV | 1 | 1 | 1 |
| OR-WA | 15 | 6 | 1 |
| UT-AZ | 3 | 3 | 0 |
| UT-CO | 1 | 1 | 1 |
| UT-ID | 2 | 0 | 1 |
| UT-NV | 1 | 0 | 0 |
| UT-WY | 4 | 0 | 0 |
| WA-ID | 7 | 0 | 5 |
| WA-OR | 25 | 5 | 2 |
| WY-CO | 3 | 1 | 1 |
| WY-ID | 6 | 7 | 2 |
| WY-MT | 5 | 2 | 1 |
| WY-UT | 2 | 0 | 3 |

## Candidate gages dropped before evaluation

| home | reason | n |
|---|---|---|
| AZ | delineation failed | 5 |
| AZ | drainage-area mismatch | 5 |
| CO | delineation failed | 9 |
| CO | no DRNAREA returned | 6 |
| ID | delineation failed | 15 |
| ID | drainage-area mismatch | 7 |
| MT | delineation failed | 3 |
| MT | drainage-area mismatch | 3 |
| NM | delineation failed | 1 |
| NM | drainage-area mismatch | 3 |
| NV | delineation failed | 1 |
| NV | drainage-area mismatch | 1 |
| OR | delineation failed | 5 |
| OR | drainage-area mismatch | 3 |
| UT | delineation failed | 4 |
| UT | drainage-area mismatch | 1 |
| WA | delineation failed | 9 |
| WA | drainage-area mismatch | 6 |
| WY | delineation failed | 2 |
| WY | drainage-area mismatch | 4 |

## B17C fits on the network

| state | gages | failed | degenerate |
|---|---|---|---|
| AZ | 150 | 4 | 1 |
| CO | 258 | 16 | 3 |
| ID | 162 | 2 | 0 |
| MT | 414 | 12 | 2 |
| NM | 223 | 10 | 5 |
| NV | 133 | 5 | 36 |
| OR | 256 | 2 | 2 |
| UT | 132 | 0 | 0 |
| WA | 381 | 3 | 7 |
| WY | 172 | 0 | 0 |

## Border discontinuity by state pair and AEP

| pair | aep | n | median_diff_log | median_abs | max_abs | share_home_higher | median_thr | flagged | median_ratio | systematic |
|---|---|---|---|---|---|---|---|---|---|---|
| AZ-UT | 0.010 | 2 | 0.381 | 0.424 | 0.805 | 0.500 | 0.730 | 0 | 2.406 | False |
| AZ-UT | 0.100 | 2 | 0.164 | 0.187 | 0.352 | 0.500 | 0.668 | 0 | 1.459 | False |
| AZ-UT | 0.500 | 2 | -0.200 | 0.200 | 0.288 | 0.000 | 0.861 | 0 | 0.631 | False |
| ID-MT | 0.010 | 9 | 0.263 | 0.272 | 0.368 | 0.778 | 0.625 | 0 | 1.833 | True |
| ID-MT | 0.100 | 9 | 0.226 | 0.233 | 0.354 | 0.778 | 0.571 | 0 | 1.682 | True |
| ID-MT | 0.500 | 9 | 0.196 | 0.200 | 0.333 | 0.667 | 0.567 | 0 | 1.569 | False |
| ID-OR | 0.010 | 4 | 0.043 | 0.058 | 0.216 | 0.750 | 0.417 | 0 | 1.104 | False |
| ID-OR | 0.100 | 4 | 0.102 | 0.102 | 0.214 | 0.750 | 0.426 | 0 | 1.265 | True |
| ID-OR | 0.500 | 4 | 0.171 | 0.171 | 0.220 | 0.750 | 0.536 | 0 | 1.482 | True |
| ID-UT | 0.010 | 3 | 0.105 | 0.105 | 0.179 | 0.667 | 0.928 | 0 | 1.274 | False |
| ID-UT | 0.100 | 3 | 0.217 | 0.217 | 0.447 | 1.000 | 0.959 | 0 | 1.649 | True |
| ID-UT | 0.500 | 3 | 0.490 | 0.490 | 0.707 | 1.000 | 1.316 | 0 | 3.092 | True |
| MT-ID | 0.010 | 22 | -0.026 | 0.143 | 0.376 | 0.409 | 0.625 | 0 | 0.941 | False |
| MT-ID | 0.100 | 22 | 0.007 | 0.132 | 0.474 | 0.500 | 0.571 | 0 | 1.016 | False |
| MT-ID | 0.500 | 22 | 0.021 | 0.098 | 0.597 | 0.545 | 0.583 | 0 | 1.050 | False |
| NM-AZ | 0.010 | 11 | -0.048 | 0.108 | 0.306 | 0.273 | 0.507 | 0 | 0.895 | False |
| NM-AZ | 0.100 | 11 | -0.065 | 0.098 | 0.290 | 0.364 | 0.486 | 0 | 0.862 | False |
| NM-AZ | 0.500 | 11 | -0.021 | 0.189 | 0.413 | 0.455 | 0.729 | 0 | 0.952 | False |
| NM-CO | 0.010 | 6 | 0.059 | 0.093 | 0.110 | 0.667 | 0.559 | 0 | 1.144 | False |
| NM-CO | 0.100 | 6 | 0.036 | 0.061 | 0.217 | 0.667 | 0.538 | 0 | 1.085 | False |
| NM-CO | 0.500 | 6 | 0.018 | 0.064 | 0.322 | 0.500 | 0.665 | 0 | 1.043 | False |
| OR-ID | 0.010 | 1 | -0.152 | 0.152 | 0.152 | 0.000 | 0.417 | 0 | 0.704 | False |
| OR-ID | 0.100 | 1 | -0.179 | 0.179 | 0.179 | 0.000 | 0.426 | 0 | 0.662 | False |
| OR-ID | 0.500 | 1 | -0.225 | 0.225 | 0.225 | 0.000 | 0.536 | 0 | 0.595 | False |
| OR-WA | 0.010 | 2 | 0.058 | 0.058 | 0.062 | 1.000 | 0.461 | 0 | 1.142 | False |
| OR-WA | 0.100 | 2 | 0.064 | 0.072 | 0.136 | 0.500 | 0.424 | 0 | 1.158 | False |
| OR-WA | 0.500 | 2 | 0.075 | 0.144 | 0.219 | 0.500 | 0.451 | 0 | 1.189 | False |
| UT-AZ | 0.010 | 3 | -0.002 | 0.017 | 0.224 | 0.333 | 0.593 | 0 | 0.996 | False |
| UT-AZ | 0.100 | 3 | -0.030 | 0.030 | 0.177 | 0.000 | 0.659 | 0 | 0.933 | False |
| UT-AZ | 0.500 | 3 | 0.005 | 0.063 | 0.144 | 0.667 | 0.928 | 0 | 1.011 | False |
| UT-ID | 0.010 | 2 | 0.192 | 0.328 | 0.520 | 0.500 | 0.771 | 0 | 1.556 | False |
| UT-ID | 0.100 | 2 | -0.044 | 0.190 | 0.234 | 0.500 | 0.754 | 0 | 0.905 | False |
| UT-ID | 0.500 | 2 | -0.352 | 0.352 | 0.362 | 0.000 | 0.945 | 0 | 0.444 | False |
| WA-ID | 0.010 | 11 | -0.172 | 0.234 | 0.419 | 0.364 | 0.615 | 0 | 0.673 | False |
| WA-ID | 0.100 | 11 | -0.150 | 0.196 | 0.312 | 0.364 | 0.605 | 0 | 0.709 | False |
| WA-ID | 0.500 | 11 | -0.170 | 0.170 | 0.297 | 0.091 | 0.798 | 0 | 0.675 | True |
| WY-CO | 0.010 | 3 | 0.021 | 0.114 | 0.125 | 0.667 | 0.545 | 0 | 1.049 | False |
| WY-CO | 0.100 | 3 | -0.009 | 0.032 | 0.174 | 0.333 | 0.556 | 0 | 0.980 | False |
| WY-CO | 0.500 | 3 | 0.024 | 0.125 | 0.203 | 0.667 | 0.731 | 0 | 1.057 | False |
| WY-MT | 0.010 | 6 | -0.019 | 0.066 | 0.271 | 0.333 | 0.519 | 0 | 0.957 | False |
| WY-MT | 0.100 | 6 | 0.120 | 0.184 | 0.395 | 0.833 | 0.529 | 0 | 1.318 | True |
| WY-MT | 0.500 | 6 | 0.194 | 0.241 | 1.022 | 0.833 | 0.724 | 1 | 1.563 | True |
| WY-UT | 0.010 | 4 | 0.010 | 0.023 | 0.416 | 0.750 | 0.466 | 0 | 1.023 | False |
| WY-UT | 0.100 | 4 | 0.065 | 0.065 | 0.280 | 1.000 | 0.434 | 0 | 1.160 | False |
| WY-UT | 0.500 | 4 | 0.071 | 0.071 | 0.106 | 1.000 | 0.529 | 0 | 1.177 | False |

## Flagged site-pairs (1)

| site_no | home | neighbor | site_name | kind | km | frac_nb | da | home_regions | neighbor_regions | neighbor_how | proxies | aeps | worst_ratio | worst_thr_ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 06298000 | WY | MT | TONGUE RIVER NEAR DAYTON, WY | border | 16.85 | 0.00 | 206.00 | 1:1.00 | SEP:1.00;UYCM:0.00 | region layer | FOREST=LC16FOREST | 0.5 | 10.52 | 7.29 |

## Border site-pairs not evaluable at AEP 0.01

| side | state | reason | n |
|---|---|---|---|
| home | NV | no equation: NV: no peak-flow region located | 10 |
| home | OR | missing 'BSLOPD', 'ELEV' | 2 |
| home | OR | missing 'ELEV' | 5 |
| home | OR | missing 'I24H2Y', 'JANMAXT2K', 'WATCAPORC', 'SOILPERM' | 1 |
| home | OR | missing 'JANAVPRE2K', 'JULAVPRE2K', 'WATCAPORC' | 5 |
| home | WA | missing 'CANOPY_PCT' | 1 |
| home | WY | missing 'JANAVPRE' | 9 |
| neighbor | CO | missing 'BSLDEM10M', 'STATSCLAY' | 2 |
| neighbor | CO | missing 'EL7500' | 6 |
| neighbor | CO | missing 'I6H100Y', 'STATSCLAY', 'OUTLETELEV' | 1 |
| neighbor | CO | missing 'STATSCLAY' | 4 |
| neighbor | ID | missing 'MINBELEV' | 3 |
| neighbor | MT | missing 'EL6000' | 5 |
| neighbor | NM | no equation: NM: no peak-flow region located | 14 |
| neighbor | NV | no equation: NV: no peak-flow region located | 7 |
| neighbor | OR | missing 'ASPECT' | 5 |
| neighbor | OR | missing 'BSLOPD' | 6 |
| neighbor | OR | missing 'BSLOPD', 'I24H2Y' | 9 |
| neighbor | OR | missing 'I24H2Y', 'JANMAXT2K', 'WATCAPORC', 'SOILPERM' | 3 |
| neighbor | OR | missing 'JANAVPRE2K', 'JULAVPRE2K', 'WATCAPORC' | 9 |
| neighbor | WA | missing 'CANOPY_PCT' | 8 |
| neighbor | WY | no equation: WY: no peak-flow region located | 39 |

## Equations vs. gage B17C

| role | state | aep | n | bias_log | rmse_log | median_abs_log | bias_pct | flagged |
|---|---|---|---|---|---|---|---|---|
| home | AZ | 0.010 | 9 | 0.045 | 0.319 | 0.204 | 11.015 | 1 |
| home | AZ | 0.100 | 9 | 0.152 | 0.297 | 0.149 | 41.956 | 1 |
| home | AZ | 0.500 | 9 | 0.219 | 0.343 | 0.340 | 65.519 | 0 |
| home | CO | 0.010 | 16 | 0.090 | 0.267 | 0.191 | 22.900 | 1 |
| home | CO | 0.100 | 16 | 0.059 | 0.199 | 0.163 | 14.532 | 0 |
| home | CO | 0.500 | 16 | 0.014 | 0.206 | 0.139 | 3.380 | 0 |
| home | ID | 0.010 | 50 | 0.034 | 0.310 | 0.121 | 8.143 | 6 |
| home | ID | 0.100 | 50 | 0.019 | 0.276 | 0.128 | 4.421 | 8 |
| home | ID | 0.500 | 50 | -0.005 | 0.264 | 0.141 | -1.061 | 4 |
| home | MT | 0.010 | 66 | -0.016 | 0.262 | 0.156 | -3.515 | 7 |
| home | MT | 0.100 | 66 | -0.024 | 0.262 | 0.136 | -5.332 | 6 |
| home | MT | 0.500 | 66 | -0.018 | 0.330 | 0.149 | -4.099 | 5 |
| home | NM | 0.010 | 31 | -0.012 | 0.240 | 0.147 | -2.711 | 4 |
| home | NM | 0.100 | 31 | -0.033 | 0.200 | 0.083 | -7.346 | 2 |
| home | NM | 0.500 | 31 | -0.002 | 0.252 | 0.111 | -0.471 | 1 |
| home | OR | 0.010 | 12 | 0.002 | 0.214 | 0.099 | 0.412 | 1 |
| home | OR | 0.100 | 12 | -0.021 | 0.194 | 0.151 | -4.764 | 0 |
| home | OR | 0.500 | 12 | -0.062 | 0.182 | 0.119 | -13.327 | 0 |
| home | UT | 0.010 | 17 | 0.129 | 0.400 | 0.178 | 34.608 | 2 |
| home | UT | 0.100 | 17 | 0.050 | 0.340 | 0.183 | 12.300 | 2 |
| home | UT | 0.500 | 17 | -0.012 | 0.377 | 0.349 | -2.627 | 2 |
| home | WA | 0.010 | 41 | 0.044 | 0.250 | 0.169 | 10.622 | 4 |
| home | WA | 0.100 | 41 | 0.021 | 0.217 | 0.154 | 5.007 | 5 |
| home | WA | 0.500 | 41 | -0.013 | 0.220 | 0.128 | -2.883 | 3 |
| home | WY | 0.010 | 18 | 0.030 | 0.255 | 0.119 | 7.208 | 2 |
| home | WY | 0.100 | 18 | 0.083 | 0.249 | 0.125 | 21.001 | 2 |
| home | WY | 0.500 | 18 | 0.099 | 0.276 | 0.131 | 25.696 | 3 |
| neighbor | AZ | 0.010 | 22 | 0.004 | 0.299 | 0.192 | 0.971 | 2 |
| neighbor | AZ | 0.100 | 22 | -0.041 | 0.284 | 0.189 | -9.101 | 2 |
| neighbor | AZ | 0.500 | 22 | 0.005 | 0.343 | 0.246 | 1.122 | 1 |
| neighbor | CO | 0.010 | 10 | 0.158 | 0.302 | 0.159 | 43.780 | 2 |
| neighbor | CO | 0.100 | 10 | 0.066 | 0.260 | 0.088 | 16.533 | 1 |
| neighbor | CO | 0.500 | 10 | 0.018 | 0.285 | 0.139 | 4.183 | 1 |
| neighbor | ID | 0.010 | 78 | 0.020 | 0.237 | 0.133 | 4.825 | 7 |
| neighbor | ID | 0.100 | 78 | 0.005 | 0.226 | 0.133 | 1.161 | 6 |
| neighbor | ID | 0.500 | 78 | -0.008 | 0.224 | 0.124 | -1.773 | 6 |
| neighbor | MT | 0.010 | 25 | -0.248 | 0.372 | 0.260 | -43.529 | 6 |
| neighbor | MT | 0.100 | 25 | -0.227 | 0.296 | 0.232 | -40.740 | 5 |
| neighbor | MT | 0.500 | 25 | -0.212 | 0.313 | 0.200 | -38.632 | 1 |
| neighbor | OR | 0.010 | 10 | -0.175 | 0.238 | 0.139 | -33.101 | 2 |
| neighbor | OR | 0.100 | 10 | -0.234 | 0.295 | 0.211 | -41.671 | 3 |
| neighbor | OR | 0.500 | 10 | -0.314 | 0.376 | 0.262 | -51.526 | 3 |
| neighbor | UT | 0.010 | 17 | -0.060 | 0.383 | 0.141 | -12.977 | 0 |
| neighbor | UT | 0.100 | 17 | -0.058 | 0.264 | 0.257 | -12.493 | 1 |
| neighbor | UT | 0.500 | 17 | -0.040 | 0.394 | 0.306 | -8.862 | 4 |
| neighbor | WA | 0.010 | 12 | -0.362 | 0.409 | 0.316 | -56.530 | 5 |
| neighbor | WA | 0.100 | 12 | -0.449 | 0.529 | 0.495 | -64.410 | 7 |
| neighbor | WA | 0.500 | 12 | -0.565 | 0.700 | 0.727 | -72.786 | 7 |

## Equation vs. B17C disagreements beyond the threshold (76)

`worst_ratio` is equation / B17C at the AEP with the largest log difference.

| role | state | region | site_no | site_name | kind | aeps | worst_ratio |
|---|---|---|---|---|---|---|---|
| home | AZ | 1 | 09383500 | NUTRIOSO CREEK ABV NELSON RES NR SPRINGERVILLE, AZ | border | 0.2 | 2.95 |
| home | AZ | 4 | 09457000 | SAN SIMON RIVER NEAR SOLOMON, ARIZ. | straddles | 0.2,0.1,0.04,0.02,0.01,0.002 | 4.33 |
| home | CO | RioGrande | 08248500 | SAN ANTONIO RIVER AT MOUTH, NEAR MANASSA, CO. | straddles | 0.002 | 2.79 |
| home | CO | Southwest | 09363100 | SALT CREEK NEAR OXFORD, CO. | border | 0.04,0.02,0.01,0.002 | 5.75 |
| home | ID | 1_2 | 12392100 | TRAPPER CREEK NR CLARK FORK ID | border | 0.1,0.04,0.02,0.01,0.002 | 0.08 |
| home | ID | 1_2 | 12415200 | PLUMMER CREEK TRIB AT PLUMMER ID | border | 0.5,0.2,0.1,0.04,0.02,0.01 | 0.17 |
| home | ID | 7 | 13083000 | TRAPPER CREEK NR OAKLEY ID | border | 0.1,0.04 | 3.75 |
| home | ID | 6_8 | 13116000 | MEDICINE LODGE CREEK AT ELLIS RANCH NR ANGORA ID | straddles | 0.5,0.2,0.1,0.04,0.02 | 3.63 |
| home | ID | 6_8 | 13116500 | MEDICINE LODGE CREEK NR SMALL ID | straddles | 0.5,0.2,0.1,0.04,0.02 | 4.00 |
| home | ID | 6_8 | 13118700 | LITTLE LOST RIVER BL WET CREEK NR HOWE ID | inland | 0.1,0.04,0.02,0.01,0.002 | 4.91 |
| home | ID | 7 | 13153500 | MALAD RIVER NR BLISS ID | inland | 0.04,0.02,0.01,0.002 | 4.32 |
| home | ID | 6_8 | 13302000 | PAHSIMEROI RIVER NR MAY ID | inland | 0.5,0.2,0.1,0.04,0.02,0.01,0.002 | 5.44 |
| home | ID | 6_8 | 13302005 | PAHSIMEROI RIVER AT ELLIS ID | inland | 0.2,0.1,0.04,0.02,0.01,0.002 | 4.46 |
| home | MT | SW | 06014000 | Red Rock River near Dell MT | inland | 0.5,0.2,0.1,0.04,0.02 | 10.76 |
| home | MT | SW | 06018000 | Beaverhead River near Dillon MT | straddles | 0.5,0.2,0.1,0.04,0.02,0.01,0.002 | 8.48 |
| home | MT | SW | 06038550 | Cabin Creek near West Yellowstone MT | border | 0.5,0.2,0.1,0.04,0.02,0.01 | 0.20 |
| home | MT | UYCM | 06187915 | Soda Butte Cr at Park Bndry at Silver Gate | straddles | 0.5,0.2,0.1 | 0.21 |
| home | MT | SEP | 06306300 | Tongue River at State Line nr Decker MT | straddles | 0.5,0.2,0.1,0.04,0.02,0.01 | 0.11 |
| home | MT | SEP | 06325500 | Little Powder River near Broadus MT | straddles | 0.04,0.02,0.01,0.002 | 5.22 |
| home | MT | W | 12303400 | Ross Creek near Troy MT | border | 0.1,0.04,0.02,0.01,0.002 | 0.19 |
| home | MT | W | 12340000 | Blackfoot River near Bonner MT | inland | 0.01,0.002 | 2.50 |
| home | MT | W | 12354000 | St. Regis River near St. Regis, MT | straddles | 0.04,0.02,0.01,0.002 | 0.24 |
| home | NM | 5 | 08253500 | SANTISTEVAN CREEK NEAR COSTILLA, NM | border | 0.01,0.002 | 3.89 |
| home | NM | 2 | 09367860 | CHUSCA WASH NEAR MEXICAN SPRINGS, NM | border | 0.2,0.1,0.04,0.02,0.01 | 0.22 |
| home | NM | 7 | 09431000 | GILA RIVER NEAR CLIFF, NM | inland | 0.01,0.002 | 3.18 |
| home | NM | 8 | 09442740 | TULAROSA RIVER NEAR RESERVE, NM | inland | 0.5,0.2,0.1 | 3.99 |
| home | NM | 7 | 09455800 | STEINS CREEK AT STEINS, NM | border | 0.01,0.002 | 2.83 |
| home | OR | E6 | 10406500 | TROUT CREEK NEAR DENIO,NEV. | border | 0.04,0.02,0.01 | 3.38 |
| home | UT | 6 | 09182000 | CASTLE CREEK ABOVE DIVERSIONS, NEAR MOAB, UTAH | border | 0.5,0.2,0.1,0.04,0.02,0.01,0.002 | 23.48 |
| home | UT | 1 | 10023000 | BIG CREEK NEAR RANDOLPH, UT | border | 0.5,0.2,0.1,0.04,0.02,0.01 | 4.98 |
| home | WA | 4 | 12010600 | LANE CREEK NEAR NASELLE, WA | border | 0.002 | 2.52 |
| home | WA | 4 | 14111400 | KLICKITAT RIVER BL SUMMIT CREEK NEAR GLENWOOD, WA | inland | 0.5,0.2,0.1,0.04,0.02,0.01,0.002 | 5.01 |
| home | WA | 4 | 14112200 | LITTLE KLICKITAT RIVER TRIB NEAR GOLDENDALE, WA | border | 0.1,0.04,0.02,0.01,0.002 | 0.19 |
| home | WA | 4 | 14113000 | KLICKITAT RIVER NEAR PITT, WA | border | 0.5,0.2,0.1,0.04 | 2.44 |
| home | WA | 4 | 14123000 | WHITE SALMON RIVER AT HUSUM, WA | border | 0.5,0.2,0.1,0.04,0.02,0.01,0.002 | 3.28 |
| home | WA | 4 | 14123500 | WHITE SALMON RIVER NEAR UNDERWOOD, WA | border | 0.2,0.1,0.04,0.02,0.01,0.002 | 2.92 |
| home | WY | 1 | 06188000 | Lamar River nr Tower Ranger Station YNP | straddles | 0.5,0.2,0.1,0.04,0.02,0.01,0.002 | 0.31 |
| home | WY | 4 | 06755000 | SOUTH CROW CREEK NEAR HECLA, WY | border | 0.5 | 2.48 |
| home | WY | 6 | 09258980 | MUDDY CREEK BELOW YOUNG DRAW, NEAR BAGGS, WY | border | 0.5,0.2,0.1,0.04,0.02,0.01 | 6.51 |
| neighbor | AZ | 2 | 09367860 | CHUSCA WASH NEAR MEXICAN SPRINGS, NM | border | 0.5,0.2,0.1,0.04,0.02,0.01 | 0.15 |
| neighbor | AZ | 2 | 09367950 | CHACO RIVER NEAR WATERFLOW , NM | straddles | 0.002 | 3.50 |
| neighbor | AZ | 3 | 09419910 | GYPSUM WASH AT NORTHSHORE RD NR LAS VEGAS BAY, NV | border | 0.2,0.1 | 0.19 |
| neighbor | AZ | 4 | 09431000 | GILA RIVER NEAR CLIFF, NM | inland | 0.04,0.02,0.01,0.002 | 3.50 |
| neighbor | CO | RioGrande | 08253500 | SANTISTEVAN CREEK NEAR COSTILLA, NM | border | 0.01,0.002 | 3.88 |
| neighbor | CO | Northwest | 09258980 | MUDDY CREEK BELOW YOUNG DRAW, NEAR BAGGS, WY | border | 0.5,0.2,0.1,0.04,0.02,0.01 | 6.16 |
| neighbor | ID | 6_8 | 06006000 | Red Rock Cr ab Lakes, nr Lakeview, MT | straddles | 0.2,0.1,0.04,0.01 | 3.12 |
| neighbor | ID | 6_8 | 06014000 | Red Rock River near Dell MT | inland | 0.5,0.2,0.1,0.04 | 5.27 |
| neighbor | ID | 6_8 | 06018000 | Beaverhead River near Dillon MT | straddles | 0.5,0.2,0.1,0.04,0.02,0.01,0.002 | 5.01 |
| neighbor | ID | 6_8 | 06191500 | Yellowstone River at Corwin Springs MT | inland | 0.5,0.2 | 0.38 |
| neighbor | ID | 6_8 | 06192500 | Yellowstone River near Livingston, MT | inland | 0.5,0.2 | 0.36 |
| neighbor | ID | 6_8 | 13018350 | FLAT CREEK BELOW CACHE CREEK, NEAR JACKSON, WY | inland | 0.5,0.2,0.1,0.04,0.02,0.01,0.002 | 5.59 |
| neighbor | ID | 6_8 | 13025500 | CROW CREEK NEAR FAIRVIEW, WY | straddles | 0.5,0.2,0.1,0.04,0.02,0.01 | 3.38 |
| neighbor | ID | 3 | 13334450 | ASOTIN CREEK BELOW CONFLUENCE NEAR ASOTIN, WA | border | 0.04,0.02,0.01,0.002 | 2.63 |
| neighbor | ID | 3 | 13334500 | ASOTIN CREEK NEAR ASOTIN, WA | border | 0.2,0.1,0.04,0.02,0.01,0.002 | 3.95 |
| neighbor | ID | 3 | 13348400 | MISSOURI FLAT CREEK TRIB NEAR PULLMAN, WA | border | 0.01,0.002 | 0.29 |
| neighbor | MT | UYCM | 06188000 | Lamar River nr Tower Ranger Station YNP | straddles | 0.2,0.1,0.04,0.02,0.01,0.002 | 0.26 |
| neighbor | MT | SEP | 06298000 | TONGUE RIVER NEAR DAYTON, WY | border | 0.5,0.2 | 0.09 |
| neighbor | MT | W | 12392100 | TRAPPER CREEK NR CLARK FORK ID | border | 0.2,0.1,0.04,0.02,0.01,0.002 | 0.05 |
| neighbor | MT | W | 12392155 | LIGHTNING CREEK AT CLARK FORK ID | straddles | 0.2,0.1,0.04,0.02,0.01,0.002 | 0.21 |
| neighbor | MT | W | 12413000 | NF COEUR D ALENE RIVER AT ENAVILLE ID | inland | 0.2,0.1,0.04,0.02,0.01,0.002 | 0.20 |
| neighbor | MT | W | 12413500 | COEUR D ALENE RIVER NR CATALDO ID | inland | 0.2,0.1,0.04,0.02,0.01,0.002 | 0.22 |
| neighbor | MT | W | 12414500 | ST JOE RIVER AT CALDER, ID | inland | 0.02,0.01,0.002 | 0.36 |
| neighbor | OR | E3 | 13339000 | CLEARWATER RIVER AT KAMIAH ID | inland | 0.5,0.2,0.1,0.04,0.02,0.01 | 0.25 |
| neighbor | OR | E3 | 13340000 | CLEARWATER RIVER AT OROFINO ID | inland | 0.5,0.2,0.1,0.04 | 0.27 |
| neighbor | OR | E3 | 13341000 | NF CLEARWATER RIVER AT AHSAHKA ID | inland | 0.5,0.2,0.1,0.04,0.02,0.01,0.002 | 0.25 |
| neighbor | UT | 6 | 09404110 | HAVASU CREEK AT SUPAI, AZ | inland | 0.5,0.2,0.1,0.04 | 4.47 |
| neighbor | UT | 7 | 09414900 | BEAVER DAM WASH AT BEAVER DAM, AZ | straddles | 0.5 | 4.30 |
| neighbor | UT | 1 | 10027000 | TWIN CREEK AT SAGE, WY | border | 0.5 | 3.11 |
| neighbor | UT | 1 | 10093000 | CUB RIVER NR PRESTON ID | border | 0.5 | 0.30 |
| neighbor | WA | 4 | 14120000 | HOOD RIVER AT TUCKER BRIDGE, NEAR HOOD RIVER, OR | border | 0.5,0.2,0.1 | 0.27 |
| neighbor | WA | 4 | 14138850 | BULL RUN RIVER NEAR MULTNOMAH FALLS, OR | border | 0.5,0.2,0.1,0.04,0.02,0.01,0.002 | 0.09 |
| neighbor | WA | 4 | 14138870 | FIR CREEK NEAR BRIGHTWOOD, OR | border | 0.5,0.2,0.1,0.04,0.02,0.01,0.002 | 0.13 |
| neighbor | WA | 4 | 14139700 | CEDAR CREEK NEAR BRIGHTWOOD, OREG. | border | 0.5,0.2,0.1,0.04,0.02,0.01,0.002 | 0.10 |
| neighbor | WA | 4 | 14139800 | SOUTH FORK BULL RUN RIVER NEAR BULL RUN, OR | border | 0.5,0.2,0.1,0.04,0.02,0.01,0.002 | 0.11 |
| neighbor | WA | 4 | 14141500 | LITTLE SANDY RIVER NEAR BULL RUN, OR | border | 0.5,0.2,0.1,0.04,0.02,0.01,0.002 | 0.13 |
| neighbor | WA | 4 | 14251500 | YOUNGS RIVER NEAR ASTORIA, OREG. | border | 0.5,0.2,0.1,0.04,0.02 | 0.13 |

## Home-state equations vs. B17C by dominant region (AEP 0.01)

| state | region | aep | n | bias_log | rmse_log | median_abs_log | bias_pct |
|---|---|---|---|---|---|---|---|
| AZ | 1 | 0.010 | 1 | 0.330 | 0.330 | 0.330 | 113.636 |
| AZ | 2 | 0.010 | 4 | 0.057 | 0.238 | 0.101 | 14.081 |
| AZ | 3 | 0.010 | 1 | -0.419 | 0.419 | 0.419 | -61.937 |
| AZ | 4 | 0.010 | 3 | 0.090 | 0.367 | 0.204 | 22.972 |
| CO | Foothills | 0.010 | 1 | -0.193 | 0.193 | 0.193 | -35.818 |
| CO | Mountain | 0.010 | 2 | 0.007 | 0.094 | 0.093 | 1.651 |
| CO | Plains | 0.010 | 1 | 0.147 | 0.147 | 0.147 | 40.400 |
| CO | RioGrande | 0.010 | 1 | 0.304 | 0.304 | 0.304 | 101.199 |
| CO | Southwest | 0.010 | 11 | 0.105 | 0.297 | 0.235 | 27.489 |
| ID | 1_2 | 0.010 | 17 | -0.155 | 0.316 | 0.129 | -29.997 |
| ID | 3 | 0.010 | 3 | -0.084 | 0.139 | 0.113 | -17.502 |
| ID | 4 | 0.010 | 6 | -0.044 | 0.060 | 0.035 | -9.579 |
| ID | 5 | 0.010 | 2 | -0.003 | 0.069 | 0.069 | -0.585 |
| ID | 6_8 | 0.010 | 17 | 0.208 | 0.350 | 0.129 | 61.386 |
| ID | 7 | 0.010 | 5 | 0.263 | 0.436 | 0.465 | 83.436 |
| MT | SEP | 0.010 | 8 | 0.030 | 0.342 | 0.252 | 7.113 |
| MT | SW | 0.010 | 18 | 0.047 | 0.307 | 0.173 | 11.319 |
| MT | UYCM | 0.010 | 15 | -0.083 | 0.243 | 0.222 | -17.469 |
| MT | W | 0.010 | 25 | -0.034 | 0.200 | 0.119 | -7.545 |
| NM | 1 | 0.010 | 4 | -0.101 | 0.156 | 0.112 | -20.760 |
| NM | 2 | 0.010 | 6 | -0.155 | 0.325 | 0.195 | -30.084 |
| NM | 5 | 0.010 | 8 | 0.104 | 0.211 | 0.141 | 26.994 |
| NM | 7 | 0.010 | 8 | 0.060 | 0.220 | 0.120 | 14.826 |
| NM | 8 | 0.010 | 4 | -0.104 | 0.280 | 0.295 | -21.342 |
| NM | 9 | 0.010 | 1 | 0.073 | 0.073 | 0.073 | 18.354 |
| OR | 2 | 0.010 | 2 | -0.142 | 0.169 | 0.142 | -27.955 |
| OR | E1 | 0.010 | 1 | -0.066 | 0.066 | 0.066 | -14.087 |
| OR | E3 | 0.010 | 6 | -0.015 | 0.185 | 0.124 | -3.400 |
| OR | E5 | 0.010 | 2 | -0.033 | 0.045 | 0.033 | -7.343 |
| OR | E6 | 0.010 | 1 | 0.529 | 0.529 | 0.529 | 237.706 |
| UT | 1 | 0.010 | 5 | 0.020 | 0.266 | 0.174 | 4.635 |
| UT | 3 | 0.010 | 3 | 0.400 | 0.413 | 0.405 | 151.229 |
| UT | 6 | 0.010 | 5 | 0.309 | 0.571 | 0.126 | 103.634 |
| UT | 7 | 0.010 | 4 | -0.162 | 0.238 | 0.170 | -31.159 |
| WA | 1 | 0.010 | 14 | -0.021 | 0.170 | 0.124 | -4.793 |
| WA | 2 | 0.010 | 2 | 0.090 | 0.327 | 0.314 | 22.963 |
| WA | 4 | 0.010 | 25 | 0.077 | 0.279 | 0.169 | 19.306 |
| WY | 1 | 0.010 | 15 | -0.014 | 0.222 | 0.117 | -3.119 |
| WY | 4 | 0.010 | 1 | 0.121 | 0.121 | 0.121 | 32.064 |
| WY | 6 | 0.010 | 2 | 0.315 | 0.457 | 0.331 | 106.477 |

## DAR transposition LOOCV, all targets by AEP

| aep | n | bias_log | rmse_log | median_abs_log | bias_pct |
|---|---|---|---|---|---|
| 0.002 | 980.000 | -0.000 | 0.425 | 0.196 | -0.115 |
| 0.010 | 980.000 | -0.001 | 0.356 | 0.168 | -0.297 |
| 0.020 | 980.000 | -0.002 | 0.333 | 0.156 | -0.382 |
| 0.040 | 980.000 | -0.002 | 0.315 | 0.149 | -0.492 |
| 0.100 | 980.000 | -0.003 | 0.302 | 0.148 | -0.683 |
| 0.200 | 980.000 | -0.004 | 0.305 | 0.146 | -0.890 |
| 0.500 | 980.000 | -0.006 | 0.341 | 0.157 | -1.381 |

## DAR transposition LOOCV by state

| state | aep | n | bias_log | rmse_log | median_abs_log | bias_pct |
|---|---|---|---|---|---|---|
| AZ | 0.010 | 44 | -0.034 | 0.368 | 0.228 | -7.556 |
| AZ | 0.100 | 44 | -0.027 | 0.258 | 0.161 | -6.106 |
| AZ | 0.500 | 44 | -0.053 | 0.326 | 0.182 | -11.515 |
| CO | 0.010 | 180 | -0.011 | 0.377 | 0.181 | -2.437 |
| CO | 0.100 | 180 | -0.020 | 0.260 | 0.157 | -4.457 |
| CO | 0.500 | 180 | -0.026 | 0.294 | 0.151 | -5.750 |
| ID | 0.010 | 88 | -0.005 | 0.285 | 0.089 | -1.172 |
| ID | 0.100 | 88 | -0.005 | 0.254 | 0.077 | -1.151 |
| ID | 0.500 | 88 | -0.007 | 0.252 | 0.088 | -1.656 |
| MT | 0.010 | 203 | 0.010 | 0.428 | 0.196 | 2.258 |
| MT | 0.100 | 203 | 0.016 | 0.361 | 0.148 | 3.853 |
| MT | 0.500 | 203 | 0.025 | 0.450 | 0.215 | 5.944 |
| OR | 0.010 | 150 | 0.006 | 0.304 | 0.153 | 1.319 |
| OR | 0.100 | 150 | 0.002 | 0.285 | 0.148 | 0.517 |
| OR | 0.500 | 150 | -0.001 | 0.290 | 0.146 | -0.253 |
| UT | 0.010 | 68 | -0.033 | 0.353 | 0.229 | -7.254 |
| UT | 0.100 | 68 | -0.030 | 0.308 | 0.156 | -6.729 |
| UT | 0.500 | 68 | -0.027 | 0.373 | 0.237 | -6.034 |
| WA | 0.010 | 247 | 0.008 | 0.326 | 0.159 | 1.897 |
| WA | 0.100 | 247 | 0.003 | 0.308 | 0.157 | 0.637 |
| WA | 0.500 | 247 | -0.006 | 0.319 | 0.165 | -1.299 |
