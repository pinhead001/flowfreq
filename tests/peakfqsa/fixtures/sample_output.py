"""
Sample PeakfqSA output text for parser tests.

Simulates the output format from a PeakfqSA analysis of the Big Sandy River
station. Values match the expected results from the PeakfqSA User Manual.
"""

# Realistic PeakfqSA .out file for Big Sandy River
BIG_SANDY_OUTPUT = """\
 PeakfqSA - Bulletin 17C Flood-Frequency Analysis
 Version 1.0

 Station: BIG SANDY RIVER AT BRUCETON, TN, 1890-1973
 Analysis Period: 1890 - 1973

 Number of peaks:              47
 Number of systematic peaks:   44
 Number of historical peaks:    3

 --- Fitted Log-Pearson Type III Parameters ---
 Mean            =  3.717272
 Std. Dev        =  0.289200
 Skew (weighted) = -0.118702
 Skew (at-site)  = -0.126000
 Mean (at-site)  =  3.720000
 Std Dev (at-site)=  0.290000
 MSE of Skew     =  0.120000
 Regional Skew   = -0.500000
 Regional MSE    =  0.302500

 Skew Option: Weighted
 Weight Option: HWN

 --- Low Outlier Analysis ---
 Method: MGBT
 Low outlier threshold: 0 PILFs detected, threshold = 1234.56

 --- Frequency Curve ---
 AEP        Estimate     Conf_Low     Conf_Up      Variance
 0.9950       871.25       540.00      1302.00       0.0150
 0.9900      1045.59       672.00      1522.00       0.0130
 0.9500      1706.18      1204.00      2310.00       0.0090
 0.9000      2203.77      1620.00      2890.00       0.0070
 0.8000      2990.15      2310.00      3780.00       0.0055
 0.6667      3957.50      3150.00      4890.00       0.0042
 0.5000      5284.36      4320.00      6450.00       0.0033
 0.2000      9166.15      7560.00     11200.00       0.0030
 0.1000     12134.65      9766.00     15218.32       0.0035
 0.0400     16276.60     12900.00     21200.00       0.0050
 0.0200     19617.73     15154.99     29124.18       0.0070
 0.0100     23158.65     17388.03     37986.08       0.0095
 0.0050     26912.12     19800.00     46200.00       0.0125
 0.0020     32217.14     22800.00     60500.00       0.0170

 --- End of Analysis ---
"""

# Minimal output with just parameters and a few quantiles
MINIMAL_OUTPUT = """\
 PeakfqSA - Bulletin 17C Flood-Frequency Analysis

 Station: TEST_STATION
 Analysis Period: 2000 - 2020

 Number of peaks:              21
 Number of systematic peaks:   21
 Number of historical peaks:    0

 --- Fitted Log-Pearson Type III Parameters ---
 Mean            =  2.500000
 Std. Dev        =  0.300000
 Skew (weighted) =  0.100000

 Skew Option: Station

 --- Low Outlier Analysis ---
 Method: NONE

 --- Frequency Curve ---
 AEP        Estimate     Conf_Low     Conf_Up      Variance
 0.0100      5000.00      3500.00      7200.00       0.0100
 0.0200      4200.00      3000.00      6100.00       0.0085

 --- End of Analysis ---
"""

# Output with station skew option (no regional skew info)
STATION_SKEW_OUTPUT = """\
 PeakfqSA - Bulletin 17C Flood-Frequency Analysis

 Station: STATION_SKEW_TEST
 Analysis Period: 1950 - 2000

 Number of peaks:              51
 Number of systematic peaks:   51
 Number of historical peaks:    0

 --- Fitted Log-Pearson Type III Parameters ---
 Mean            =  3.000000
 Std. Dev        =  0.250000
 Skew (weighted) =  0.300000
 Skew (at-site)  =  0.300000

 Skew Option: Station
 Weight Option: HWN

 --- Frequency Curve ---
 AEP        Estimate     Conf_Low     Conf_Up      Variance
 0.0100      8000.00      6000.00     11000.00       0.0080
 0.0200      6500.00      5000.00      8800.00       0.0065

 --- End of Analysis ---
"""
