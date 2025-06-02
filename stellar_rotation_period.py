import numpy as np
import matplotlib.pyplot as plt
import lightkurve as lk
from PyAstronomy.pyasl import foldAt

# -----------------------------
# STEP 1: Download and Prepare Light Curve Data
# -----------------------------

# Define the TIC ID of the target
tic_id = 461591646  

# Search and download TESS light curves for sectors 59 and 60 from the SPOC pipeline
lc = lk.search_lightcurve(f"TIC {tic_id}", mission="TESS", author='SPOC', sector=[59, 60]).download()
print(lc)  # Print info about the downloaded light curve

# Remove NaN values and normalize the flux
lc = lc.remove_nans().flatten(window_length=4321, polyorder=2).normalize()


# Extract time (in days) and normalized flux arrays
time, flux = lc.time.value, lc.flux.value

# -----------------------------
# STEP 2: Define Transit Parameters and Mask Out Transit
# -----------------------------

# Known transit parameters
t0 = 2459485.6499        # Mid-transit time (BJD)
P = 3.7246957            # Orbital period (days)
T_T_days = 3.078 / 24    # Transit duration (in days)

# Compute orbital phase for each observation
phase = ((time - t0) / P) % 1

# Mask data outside of the transit window to isolate "out-of-transit" data
out_of_transit_mask = (np.abs(phase) > (T_T_days / (2 * P)))

# Extract time and flux values for only out-of-transit points
time_out = time[out_of_transit_mask]
flux_out = flux[out_of_transit_mask]


# STEP 3: Phase-Fold with Secondary Period and Bin the Data
# -----------------------------

# Define a second period for phase folding (e.g., suspected secondary modulation)
second_period = 5.1  # in days

# Fold the out-of-transit time data over the second period
time_folded_second = foldAt(time_out, second_period, t0)

# Choose bin size (2 hours in days)
bin_size = 2 / 24  # 2 hours in days

# Create bin edges across the secondary period
bin_edges = np.arange(0, second_period, bin_size)

# Assign each data point to a bin
digitized = np.digitize(time_folded_second, bin_edges)

# Calculate mean flux within each bin
binned_flux = [np.mean(flux_out[digitized == i]) for i in range(1, len(bin_edges))]

# Compute bin centers for plotting
binned_time = [np.mean(bin_edges[i:i+2]) for i in range(len(bin_edges) - 1)]

# -----------------------------
# STEP 4: Plot the Results
# -----------------------------

plt.figure(figsize=(8, 6))

# Plot original phase-folded out-of-transit data
plt.plot(time_folded_second, flux_out, '.', label="Out-of-transit Data (TESS)", 
         color="gray", lw=0.5, alpha=0.5, zorder=1)

# Plot binned data over phase
plt.scatter(binned_time, binned_flux, label="2-hour Binning", 
            color="red", alpha=1, zorder=2)

# Add reference line at normalized flux = 1
plt.axhline(1, color='black', linestyle='--')

# Label axes and add legend
plt.xlabel("Phase (5.1 day period)", fontsize=12)
plt.ylabel("Normalized Flux", fontsize=12)
plt.legend(loc='best')

# Set plot limits and layout
plt.ylim(0.99290, 1.00750)
plt.tight_layout()

# Save the figure to a file
plt.savefig("outer_transit_data_s59_60.png", dpi=400)

# Display the plot
plt.show()

