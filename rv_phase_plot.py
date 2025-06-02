# === Import required libraries ===
from PyAstronomy.pyasl import foldAt         # For folding time-series data on a given period
import numpy as np                           # For numerical operations
import pandas as pd                          # For reading and manipulating CSV data
import matplotlib.pyplot as plt              # For plotting
from matplotlib import gridspec              # For more customizable subplot layout




# === Load and prepare radial velocity data ===
data = pd.read_csv("Raphael.csv")            # Read RV data from CSV
print(data.columns)                          # Print column names to verify structure
data.columns = data.columns.str.strip()      # Remove any leading/trailing spaces from column names

# Extract necessary columns
time = data['#time']
flux = data['RV']
flux_err = data['RV_err']

# === Phase-fold the data using known period and T0 ===
phases = foldAt(time, 3.7246957, T0=2459891.63471)  # Period and T0 from your fit

# Sort data by phase for clean plotting
sortIndi = np.argsort(phases)
phases = phases[sortIndi]
flux = flux[sortIndi]
flux_err = flux_err[sortIndi]

# === Load the model and residuals from CSV ===
modelcsv = pd.read_csv("RV_model.csv")
model = modelcsv['model']
modeltime = modelcsv['model_time']
baseline = modelcsv['baseline']
residuals = modelcsv['residuals']

# === Fold the model data using same period and T0 ===
phase_model = foldAt(modeltime, 3.7246957, T0=2459891.63471)

# Sort model and residuals by phase
sortIndi2 = np.argsort(phase_model)
phase_model = phase_model[sortIndi2]
model = model[sortIndi2] + baseline[sortIndi2]  # Combine model and baseline
residuals = residuals[sortIndi2]

# === Start plotting ===
fig = plt.figure(figsize=(8, 10))
gs = gridspec.GridSpec(2, 1, height_ratios=[5, 2])  # Two vertically stacked subplots

# --- Top plot: Radial velocity data and model ---
ax0 = plt.subplot(gs[0])
ax0.plot(phase_model, model, 'r', alpha=1, label="Model (GP Matérn 3/2)", zorder=2)
ax0.errorbar(phases, flux, flux_err, fmt='b.', alpha=1, label="Data (TRES)",
             markersize=5, capsize=2, zorder=1)

ax0.set_ylabel("Radial Velocity [km/s]", fontsize=12)
ax0.legend(ncol=2, loc='upper left', fontsize=12)

# --- Bottom plot: Residuals ---
ax1 = plt.subplot(gs[1])
residuals = residuals[sortIndi]  # Match residuals to data phase order
ax1.scatter(phases, residuals, c='b', label='Residuals', alpha=0.6)
ax1.axhline(0, color='grey', linestyle='--')  # Horizontal line at 0 for reference

ax1.set_xlabel("Phase", fontsize=12) 
ax1.set_ylabel("Residuals", fontsize=12)

# === Final adjustments ===
plt.tight_layout()
plt.savefig('RVplot_phase.png', dpi=400)  # Save figure as high-res image
plt.show()                                # Display the plot








#for dark background
#from PyAstronomy.pyasl import foldAt
#import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt
#from matplotlib import gridspec
#from scipy.signal import savgol_filter
#from scipy.stats import norm
#import seaborn as sns
#
## Use dark background
#plt.style.use('dark_background')
#
## --- Load data ---
#data = pd.read_csv("Raphael.csv")
#data.columns = data.columns.str.strip()
#
#time = data['#time']
#flux = data['RV']
#flux_err = data['RV_err']
#
## Fold the data on the period
#period = 3.7246871
#t0 = 2459798.51786
#phases = foldAt(time, period, T0=t0)
#
## Sort by phase
#sortIndi = np.argsort(phases)
#phases = phases[sortIndi]
#flux = flux[sortIndi]
#flux_err = flux_err[sortIndi]
#
## --- Load model data ---
#modelcsv = pd.read_csv("RV_model.csv")
#model = modelcsv['model']
#model_time = modelcsv['model_time']
#baseline = modelcsv['baseline']
#residuals = modelcsv['residuals']
#
## Fold and sort model
#phase_model = foldAt(model_time, period, T0=t0)
#sortIndi2 = np.argsort(phase_model)
#phase_model = phase_model[sortIndi2]
#model = model[sortIndi2] + baseline[sortIndi2]
#residuals = residuals[sortIndi2]
#residuals = residuals[sortIndi]  # Match to data phase order
#
## --- Create plot ---
#fig = plt.figure(figsize=(8, 10))
#gs = gridspec.GridSpec(2, 1, height_ratios=[5, 2])
#
## Main panel
#ax0 = plt.subplot(gs[0])
#ax0.plot(phase_model, model, color='orange', lw=2, label="Model (GP Matérn 3/2)", zorder=2)
#ax0.errorbar(phases, flux, flux_err, fmt='o', color='cyan', ecolor='white',
#             elinewidth=1, capsize=2, markersize=4, label="Data (TRES)", zorder=1)
#
#ax0.set_ylabel("Radial Velocity [km/s]", fontsize=12)
#ax0.legend(loc='upper left', fontsize=12)
##ax0.grid(alpha=0.2, color='white', linestyle='--')
#
## Residuals panel
#ax1 = plt.subplot(gs[1])
#ax1.axhline(0, color='white', linestyle='--', lw=1, alpha=0.7)
#ax1.scatter(phases, residuals, color='lime', alpha=0.7, s=15, label='Residuals')
#
#ax1.set_xlabel("Phase", fontsize=12)
#ax1.set_ylabel("Residuals", fontsize=12)
##ax1.grid(alpha=0.2, color='white', linestyle='--')
#
## Axis and tick styling
#for ax in [ax0, ax1]:
#    ax.tick_params(colors='white', labelsize=10)
#    for spine in ax.spines.values():
#        spine.set_color('white')
#
## Final layout and save
#plt.tight_layout()
#plt.savefig('RVplot_phase_dark.png', dpi=400, facecolor='black')
#plt.show()
