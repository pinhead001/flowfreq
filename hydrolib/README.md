# HydroLib - Python Library for Hydrologic Analysis

A comprehensive Python library for downloading USGS streamflow data, performing Bulletin 17C flood frequency analysis, and generating professional technical reports.

## Features

- **USGS Data Retrieval**: Download mean daily flow and annual peak flow data from USGS NWIS
- **Hydrograph Analysis**: Generate summary hydrographs with day of water year on x-axis
- **Bulletin 17C Analysis**: 
  - Log-Pearson Type III distribution fitting
  - Station skew or weighted regional skew support
  - Multiple Grubbs-Beck test for low outlier detection
  - 90% confidence intervals (5% limits)
- **Professional Figures**: Publication-quality plots
- **Technical Reports**: Automated report generation in Markdown format

## Installation

```bash
pip install numpy pandas matplotlib scipy requests
```

## Quick Start

```python
from hydrolib import USGSgage, Bulletin17C, Hydrograph, HydroReport, analyze_gage

# Option 1: Full automated analysis
results = analyze_gage(
    site_no='01638500',          # USGS gage number
    regional_skew=-0.05,         # Optional: regional skew coefficient
    regional_skew_mse=0.12,      # Optional: regional skew MSE
    output_dir='./output'
)

# Option 2: Step-by-step analysis
gage = USGSgage('01638500')
daily_data = gage.download_daily_flow(start_date='2010-01-01')
peak_data = gage.download_peak_flow()

# Run B17C analysis
b17c = Bulletin17C(peak_data['peak_flow_cfs'].values)
results = b17c.run_analysis()

# Generate flood frequency curve
b17c.plot_frequency_curve(gage.site_name, gage.site_no, save_path='ffa_plot.png')
```

## Classes

### USGSgage
Handles USGS gage data retrieval.

```python
gage = USGSgage('01638500')

# Download mean daily flow
daily_data = gage.download_daily_flow(
    start_date='2000-01-01',  # Optional
    end_date='2020-12-31'     # Optional
)

# Download annual peak flows
peak_data = gage.download_peak_flow()

# Access site information
print(gage.site_name)       # Station name
print(gage.drainage_area)   # Drainage area (sq mi)
```

### Hydrograph
Static methods for hydrograph plotting.

```python
# Daily time series plot
Hydrograph.plot_daily_timeseries(
    daily_data,
    site_name='Station Name',
    site_no='01638500',
    save_path='timeseries.png'
)

# Summary hydrograph (day of water year)
Hydrograph.plot_summary_hydrograph(
    daily_data,
    site_name='Station Name',
    site_no='01638500',
    save_path='hydrograph.png',
    percentiles=[10, 25, 50, 75, 90]
)
```

### Bulletin17C
Performs flood frequency analysis.

```python
# Station skew only
b17c = Bulletin17C(peak_flows_array)

# With weighted regional skew
b17c = Bulletin17C(
    peak_flows_array,
    regional_skew=-0.05,      # Regional skew coefficient
    regional_skew_mse=0.12    # Mean squared error of regional skew
)

# Run complete analysis
results = b17c.run_analysis()

# Access results
print(results['mean_log'])           # Mean of log Q
print(results['std_log'])            # Std dev of log Q
print(results['skew_station'])       # Station skew
print(results['skew_weighted'])      # Weighted skew (if regional provided)
print(results['low_outlier_threshold'])  # Grubbs-Beck threshold

# Access quantile tables
quantiles = results['quantiles']           # Flow estimates
conf_limits = results['confidence_limits'] # With 5%/95% limits

# Generate plot
b17c.plot_frequency_curve(
    site_name='Station Name',
    site_no='01638500',
    save_path='frequency_curve.png'
)
```

### HydroReport
Generates technical reports.

```python
report = HydroReport(gage, b17c)

# Generate all figures
figures = report.generate_all_figures(output_dir='./output')

# Generate markdown report
report.save_report('report.md')
```

## Output Tables

### Flood Frequency Quantiles
| AEP (%) | Return Period (yr) | Flow (cfs) | 5% Limit | 95% Limit |
|---------|-------------------|------------|----------|-----------|
| 10.0    | 10                | 77,422     | 62,411   | 96,043    |
| 4.0     | 25                | 111,195    | 85,497   | 144,617   |
| 2.0     | 50                | 141,752    | 105,260  | 190,895   |
| 1.0     | 100               | 177,426    | 127,390  | 247,115   |
| 0.2     | 500               | 284,553    | 189,779  | 426,656   |

## Generated Figures

1. **Daily Flow Time Series** - Complete daily flow record with log scale
2. **Summary Hydrograph** - Day of water year with percentile bands
3. **Annual Peak Flow** - Bar chart of annual peaks
4. **Flood Frequency Curve** - LP3 curve with confidence limits

## Technical Details

### Bulletin 17C Implementation

The library implements key elements of USGS Bulletin 17C:

1. **Log-Pearson Type III Distribution**: Statistics computed in log10 space
2. **Weighted Skew**: Combines station and regional skew using MSE weighting
3. **K-Factor Calculation**: Wilson-Hilferty approximation for LP3 quantiles
4. **Grubbs-Beck Test**: Multiple Grubbs-Beck test for low outlier detection
5. **Confidence Limits**: Approximate variance for quantile estimates

### Weighted Skew Calculation

When regional skew values are provided:

```
MSE_station = (6/n) * (1 + (6/n)*G² + (15/n²)*G⁴)

w_regional = MSE_station / (MSE_station + MSE_regional)
w_station = MSE_regional / (MSE_station + MSE_regional)

G_weighted = w_station * G_station + w_regional * G_regional
```

### Confidence Limits

90% confidence intervals (5% and 95% limits) computed as:

```
Var(X_p) = s² * [1/n + K² * (1 + 0.75*G²) / (2*(n-1))]
SE = s * sqrt(Var_factor)
CI = Q ± z_α/2 * SE
```

## References

England, J.F., Jr., et al., 2019, Guidelines for determining flood flow frequency—Bulletin 17C: U.S. Geological Survey Techniques and Methods, book 4, chap. B5, 148 p., https://doi.org/10.3133/tm4B5

## License

MIT License

## Author

HydroLib - Python Hydrologic Analysis Library
encyResults`: Complete analysis output

### Decorators Used
- `@lru_cache`: Caches `kfactor()` and `grubbs_beck_critical_value()` for performance
- `@cached_property`: Lazy evaluation of `log_flows`, `period_of_record`
- `@property`/`@setter`: Encapsulated access to gage data
- `@staticmethod`: Stateless utility methods in `Hydrograph`
- `@classmethod`: Factory methods in `FlowInterval`
- `@abstractmethod`: Interface definition in `FloodFrequencyAnalysis`

## Utility Functions

```python
from hydrolib import kfactor, grubbs_beck_critical_value

# K-factor for LP3 distribution (cached)
K = kfactor(skew=0.3, aep=0.01)  # 1% AEP

# Grubbs-Beck critical value (cached)
K_n = grubbs_beck_critical_value(n=50)
```

## Results Access

```python
results = b17c.run_analysis(method='ema')

# Statistics
results.mean_log          # Mean of log Q
results.std_log           # Std dev of log Q
results.skew_station      # Station skew
results.skew_weighted     # Weighted skew
results.skew_used         # Skew used in analysis

# Sample counts
results.n_peaks           # Total observations
results.n_systematic      # Systematic record count
results.n_historical      # Historical observations
results.n_censored        # Censored observations
results.n_low_outliers    # Low outliers (MGB)

# Thresholds
results.low_outlier_threshold
results.mgb_critical_value

# EMA-specific
results.ema_iterations    # Number of iterations
results.ema_converged     # Convergence status

# Tables
results.quantiles         # DataFrame of quantiles
results.confidence_limits # DataFrame with 5%/95% limits
```

## References

England, J.F., Jr., et al., 2019, Guidelines for determining flood flow frequency—Bulletin 17C: U.S. Geological Survey Techniques and Methods, book 4, chap. B5, 148 p., https://doi.org/10.3133/tm4B5

## Version History

- **v2.0.0**: Added EMA algorithm, historical flood handling, refactored with decorators
- **v1.0.0**: Initial release with MOM analysis

## License

MIT License
