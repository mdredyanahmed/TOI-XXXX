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
data = pd.read_csv("data_ground_based.csv")
time = data['#time'].values
flux = data['flux'].values
flux_err = data['flux_err'].values

# === Load model output from allesfitter ===
modelcsv = pd.read_csv("model_alldays_groundbased.csv")
model_time = modelcsv['time'].values
model = modelcsv['model'].values           # full model = transit + baseline
baseline = modelcsv['baseline'].values     # GP/systematics-only component

# === Subtract baseline to isolate transit signal ===
flux_corrected = flux - baseline
flux_err_corrected = flux_err  

# === Phase-fold observed light curve ===
period = 3.724771
T0 = 2459485.6502 # BJD_TDB

phases = foldAt(time, period, T0=T0, centralzero=True)
sortIndi = np.argsort(phases)

phases = phases[sortIndi]
flux_corrected = flux_corrected[sortIndi]
flux_err_corrected = flux_err_corrected[sortIndi]

# === Phase-fold model ===
phases_model = foldAt(model_time, period, T0=T0, centralzero=True)
sortIndi_model = np.argsort(phases_model)
phases_model = phases_model[sortIndi_model]
model_smooth = model[sortIndi_model]

# === Bin data in phase (10-minute bins) ===
bin_width_minutes = 10
bin_width_phase = (bin_width_minutes / 60.0) / 24.0 / period  # convert to phase

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

# Define colors
color_data = '#990099'     # Violet
color_binned_face = 'green'  # White fill
color_binned_edge = '#117733'  # Green border for visibility
color_model = '#D62728'     # Red
color_resid = '#990099'     # Violet (same as raw)

# === Top panel: Phase-folded light curve ===
ax0 = plt.subplot(gs[0])

# Raw data points (violet)
ax0.scatter(phases, flux_corrected, color=color_data, alpha=0.5, s=10, label="Data (Ground-based observations)")

# Binned points (white with green edge)
ax0.errorbar(bin_centers, binned_flux, yerr=binned_flux_err, fmt='o',
             markersize=6, capsize=1.5, color=color_binned_edge, 
             mfc=color_binned_face, mec=color_binned_edge, label="10-minute-Binned")

# Model curve (red)
ax0.plot(phases_model, model_smooth, color=color_model, lw=2, alpha=1, label=" Model (GP Matern 3/2)")

ax0.set_ylabel("Relative Flux-Baseline", fontsize=16)
ax0.legend(ncol=1, loc='best', frameon=False)
ax0.set_xlim(-0.03, 0.03)
ax0.minorticks_on()

# === Bottom panel: Residuals (violet) ===
ax1 = plt.subplot(gs[1], sharex=ax0)
residuals = flux - (model+baseline)
residuals = residuals[sortIndi]
ax1.scatter(phases, residuals, c=color_resid, s=10, alpha=0.5, label="Residuals")
ax1.axhline(0, color='black', linestyle='--', lw=1)
ax1.set_ylabel("Residuals", fontsize=16)
ax1.set_xlabel("Phase", fontsize=16)
ax0.set_xlim(-0.03, 0.03)
ax1.minorticks_on()

# === Final layout ===
plt.tight_layout()
plt.savefig("phaseplot_groundbased.png", dpi=600, bbox_inches='tight')
plt.show()

