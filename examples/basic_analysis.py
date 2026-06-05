"""
Basic flood frequency analysis example.

Demonstrates both MOM and EMA methods with synthetic data.
"""

import matplotlib
import numpy as np

from hydrolib import Bulletin17C, ExpectedMomentsAlgorithm, MethodOfMoments

matplotlib.use("Agg")  # Use non-interactive backend

# Generate synthetic peak flow data
np.random.seed(42)
n_years = 50
water_years = np.arange(1971, 2021)

# Log-Pearson Type III parameters
mean_log = 4.5  # ~31,600 cfs median
std_log = 0.25
skew = 0.3

# Generate LP3 random variates
alpha = 4 / skew**2
z = (np.random.gamma(alpha, 1, n_years) - alpha) / np.sqrt(alpha)
peak_flows = 10 ** (mean_log + std_log * z)

print("=" * 60)
print("HYDROLIB FLOOD FREQUENCY ANALYSIS EXAMPLE")
print("=" * 60)

# Example 1: Method of Moments
print("\n1. METHOD OF MOMENTS (MOM)")
print("-" * 40)

mom = MethodOfMoments(peak_flows, regional_skew=0.05, regional_skew_mse=0.12)
mom_results = mom.run_analysis()

print(f"Mean(log Q): {mom_results.mean_log:.4f}")
print(f"Std(log Q): {mom_results.std_log:.4f}")
print(f"Station Skew: {mom_results.skew_station:.4f}")
print(f"Weighted Skew: {mom_results.skew_weighted:.4f}")

# Example 2: EMA with historical floods
print("\n2. EXPECTED MOMENTS ALGORITHM (EMA)")
print("-" * 40)

historical_peaks = [
    (1936, 150000),  # Large historical flood
    (1955, 120000),
]

ema = ExpectedMomentsAlgorithm(
    peak_flows,
    water_years=water_years,
    regional_skew=0.05,
    regional_skew_mse=0.12,
    historical_peaks=historical_peaks,
)
ema_results = ema.run_analysis()

print(f"Mean(log Q): {ema_results.mean_log:.4f}")
print(f"Std(log Q): {ema_results.std_log:.4f}")
print(f"Station Skew: {ema_results.skew_station:.4f}")
print(f"Weighted Skew: {ema_results.skew_weighted:.4f}")
print(f"Historical peaks: {ema_results.n_historical}")
print(f"EMA iterations: {ema_results.ema_iterations}")
print(f"Converged: {ema_results.ema_converged}")

# Example 3: Unified interface
print("\n3. UNIFIED Bulletin17C INTERFACE")
print("-" * 40)

b17c = Bulletin17C(
    peak_flows,
    water_years=water_years,
    regional_skew=0.05,
    regional_skew_mse=0.12,
    historical_peaks=historical_peaks,
)

# Run with EMA (default)
results = b17c.run_analysis(method="ema")

# Print key quantiles
print("\nFlood Frequency Estimates:")
print(f"{'Return Period':>15} {'Flow (cfs)':>15} {'90% CI':>25}")
print("-" * 57)

cl = results.confidence_limits
for aep in [0.10, 0.04, 0.02, 0.01, 0.002]:
    row = cl[cl["aep"] == aep].iloc[0]
    rp = int(1 / aep)
    ci = f"({row['lower_5pct']:,.0f} - {row['upper_5pct']:,.0f})"
    print(f"{rp:>12} yr {row['flow_cfs']:>15,.0f} {ci:>25}")

# Save plot
print("\nSaving flood frequency curve...")
fig = b17c.plot_frequency_curve(
    site_name="Synthetic Example", site_no="00000000", save_path="flood_frequency_curve.png"
)
print("Saved: flood_frequency_curve.png")

print("\n" + "=" * 60)
print("ANALYSIS COMPLETE")
print("=" * 60)
