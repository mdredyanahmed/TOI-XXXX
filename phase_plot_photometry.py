# Import necessary packages
from PyAstronomy.pyasl import foldAt             # For phase-folding time series data
import matplotlib.pyplot as plt                  # For plotting
from matplotlib import gridspec                  # For more complex subplot layouts
import seaborn as sns                            # Optional, good for plot aesthetics
import numpy as np                               # For numerical operations
import pandas as pd                              # For reading CSV files
from scipy.signal import savgol_filter           # Not used here, but could smooth data
from scipy.stats import norm                     # Not used here, but useful for statistics

# Optional: Use dark background style for plots
# plt.style.use('dark_background')

# --- Load photometric time series data ---
data = pd.read_csv("TESS.csv")

# Extract time, flux, and flux error columns
time = data['#time']
flux = data['flux']
flux_err = data['flux_err']

# --- Phase-fold the data ---
# Using known period and T0 (mid-transit time)
phases = foldAt(time, 3.7246957, T0=2459891.63471, centralzero=True)

# Sort data by phase for a smooth plot
sortIndi = np.argsort(phases)
phases = phases[sortIndi]
flux = flux[sortIndi]
flux_err = flux_err[sortIndi]

# --- Load GP model fit ---
modelcsv = pd.read_csv("model_sector_all_sector_TESS.csv")
model = modelcsv['model'][sortIndi]

# --- Create plot layout with two rows (flux + residuals) ---
fig = plt.figure(figsize=(8, 10))
gs = gridspec.GridSpec(2, 1, height_ratios=[5, 2])  # First row taller than second

# === Top panel: Phase-folded light curve ===
ax0 = plt.subplot(gs[0])

# Plot the model (e.g., GP fit)
ax0.plot(phases, model, 'r', alpha=1, lw=1, zorder=2, label="Model (GP matern 32)")

# Plot TESS data with error bars
ax0.errorbar(phases, flux, flux_err, color='orange', alpha=0.8,
             markersize=1, capsize=1, zorder=1, label="Data (TESS-SPOC)")

# X-axis limits around the transit
ax0.set_xlim(-0.04, 0.04)
ax0.set_ylabel("Relative flux", fontsize=12)

# Add legend
ax0.legend(ncol=2, loc='upper left', fontsize=12)

# === Bottom panel: Residuals ===
ax1 = plt.subplot(gs[1])  # Shared subplot axis for residuals

# Compute residuals = data - model
residuals = flux - model

# Plot residuals as scatter
ax1.scatter(phases, residuals, c='orange', s=2, alpha=0.7, label='Residuals')

# Add horizontal line at y=0
ax1.axhline(0, color='grey', linestyle='--')

# Axis labels and limits
ax1.set_ylabel("Residuals", fontsize=12)
ax1.set_xlabel("Phase", fontsize=12)
ax1.set_xlim(-0.04, 0.04)

# Adjust layout and save the figure
plt.tight_layout()
plt.savefig('phaseplot_TESS.png', dpi=400)
plt.show()

