
# === Imports ===
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
from PyAstronomy.pyasl import foldAt
from scipy.interpolate import interp1d

# === Plot style for publication ===
mpl.rcParams.update({
    'font.size': 14,
    'font.family': 'serif',
    'axes.labelsize': 16,
    'axes.titlesize': 18,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'axes.linewidth': 1.2,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'xtick.top': True,
    'ytick.right': True,
    'grid.alpha': 0.3,
    'grid.linestyle': '--'
})

# === Load TESS photometric data ===
data = pd.read_csv("TESS.csv")
time = data['#time'].values
flux = data['flux'].values
flux_err = data['flux_err'].values

# === Load model output from allesfitter ===
modelcsv = pd.read_csv("model_allsector_tess.csv")
model_time = modelcsv['time'].values
model = modelcsv['model'].values           # full model = transit + baseline
baseline = modelcsv['baseline'].values     # GP/systematics-only component

# === Having relative flux by removing baseline ===
# Here we subtract the baseline 

flux_corrected = flux - baseline
flux_err_corrected = flux_err  # assuming uncertainty remains the same

# === Phase-fold observed light curve ===
period = 3.7246960
T0 = 2459891.63476  # BJD_TDB

phases = foldAt(time, period, T0=T0, centralzero=True)
sortIndi = np.argsort(phases)

phases = phases[sortIndi]
flux_corrected = flux_corrected[sortIndi]
flux_err_corrected = flux_err_corrected[sortIndi]

# === Phase-fold model for  plotting ===
phases_model = foldAt(model_time, period, T0=T0, centralzero=True)
sortIndi_model = np.argsort(phases_model)
phases_model = phases_model[sortIndi_model]
model_smooth = model[sortIndi_model]

# === Bin data in phase (10-minute bins) ===
bin_width_minutes = 10
bin_width_phase = (bin_width_minutes / 60.0) / 24.0 / period  # convert to phase units

bins = np.arange(phases.min(), phases.max() + bin_width_phase, bin_width_phase)
bin_centers = 0.5 * (bins[1:] + bins[:-1])
digitized = np.digitize(phases, bins)

binned_flux, binned_flux_err = [], []
for i in range(1, len(bins)):
    in_bin = digitized == i
    if np.any(in_bin):
        weights = 1.0 / flux_err_corrected[in_bin]**2
        binned_flux.append(np.average(flux_corrected[in_bin], weights=weights))
        binned_flux_err.append(np.sqrt(1.0 / np.sum(weights)))
    else:
        binned_flux.append(np.nan)
        binned_flux_err.append(np.nan)

binned_flux = np.array(binned_flux)
binned_flux_err = np.array(binned_flux_err)

# === Create figure with two panels (light curve + residuals) ===
fig = plt.figure(figsize=(8, 10))
gs = gridspec.GridSpec(2, 1, height_ratios=[5, 2], hspace=0.1)

# === Top panel: Phase-folded light curve ===
ax0 = plt.subplot(gs[0])
ax0.scatter(phases, flux_corrected, color='#E69F00', alpha=0.4, s=8, label="Data (TESS - SPOC)")
ax0.errorbar(bin_centers, binned_flux, yerr=binned_flux_err, fmt='o',
             markersize=2, capsize=0.5, color='#0072B2', alpha=1, label='10-minute-binned', zorder=1)
ax0.plot(phases_model, model_smooth, 'r-', lw=1, alpha=1, label="Model (GP matern 3/2)", zorder=2)

ax0.set_xlim(-0.04, 0.04)
ax0.set_ylabel("Relative Flux - Baseline", fontsize=16)
ax0.legend(ncol=2, loc='upper left', frameon=False)
ax0.minorticks_on()

# === Bottom panel: Residuals ===
ax1 = plt.subplot(gs[1], sharex=ax0)

residuals = flux - (model + baseline)  # observed minus full model
residuals = residuals[sortIndi]        # apply same sorting as phase

ax1.scatter(phases, residuals, c='#E69F00', s=8, alpha=0.4)
ax1.axhline(0, color='grey', linestyle='--', lw=1)
ax1.set_ylabel("Residuals", fontsize=16)
ax1.set_xlabel("Phase", fontsize=16)
ax1.minorticks_on()

# === Show/save ===
plt.tight_layout()
plt.savefig("phaseplot_TESS_binned.png", dpi=600, bbox_inches='tight')
plt.show()

# # Import necessary packages
# from PyAstronomy.pyasl import foldAt             # For phase-folding time series data
# import matplotlib.pyplot as plt                  # For plotting
# from matplotlib import gridspec                  # For more complex subplot layouts
# import seaborn as sns                            # Optional, good for plot aesthetics
# import numpy as np                               # For numerical operations
# import pandas as pd                              # For reading CSV files
# from scipy.signal import savgol_filter           # Not used here, but could smooth data
# from scipy.stats import norm                     # Not used here, but useful for statistics

# # Optional: Use dark background style for plots
# # plt.style.use('dark_background')

# # --- Load photometric time series data ---
# data = pd.read_csv("TESS.csv")

# # Extract time, flux, and flux error columns
# time = data['#time']
# flux = data['flux']
# flux_err = data['flux_err']

# # --- Phase-fold the data ---
# # Using known period and T0 (mid-transit time)
# phases = foldAt(time, 3.7246957, T0=2459891.63471, centralzero=True)

# # Sort data by phase for a smooth plot
# sortIndi = np.argsort(phases)
# phases = phases[sortIndi]
# flux = flux[sortIndi]
# flux_err = flux_err[sortIndi]

# # --- Load GP model fit ---
# modelcsv = pd.read_csv("model_sector_all_sector_TESS.csv")
# model = modelcsv['model'][sortIndi]

# # --- Create plot layout with two rows (flux + residuals) ---
# fig = plt.figure(figsize=(8, 10))
# gs = gridspec.GridSpec(2, 1, height_ratios=[5, 2])  # First row taller than second

# # === Top panel: Phase-folded light curve ===
# ax0 = plt.subplot(gs[0])

# # Plot the model (e.g., GP fit)
# ax0.plot(phases, model, 'r', alpha=1, lw=1, zorder=2, label="Model (GP matern 32)")

# # Plot TESS data with error bars
# ax0.errorbar(phases, flux, flux_err, color='orange', alpha=0.8,
#              markersize=1, capsize=1, zorder=1, label="Data (TESS-SPOC)")

# # X-axis limits around the transit
# ax0.set_xlim(-0.04, 0.04)
# ax0.set_ylabel("Relative flux", fontsize=12)

# # Add legend
# ax0.legend(ncol=2, loc='upper left', fontsize=12)

# # === Bottom panel: Residuals ===
# ax1 = plt.subplot(gs[1])  # Shared subplot axis for residuals

# # Compute residuals = data - model
# residuals = flux - model

# # Plot residuals as scatter
# ax1.scatter(phases, residuals, c='orange', s=2, alpha=0.7, label='Residuals')

# # Add horizontal line at y=0
# ax1.axhline(0, color='grey', linestyle='--')

# # Axis labels and limits
# ax1.set_ylabel("Residuals", fontsize=12)
# ax1.set_xlabel("Phase", fontsize=12)
# ax1.set_xlim(-0.04, 0.04)

# # Adjust layout and save the figure
# plt.tight_layout()
# plt.savefig('phaseplot_TESS.png', dpi=400)
# plt.show()

