import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import allesfitter

# === User settings ===
data_folder = 'allesfit_eccentric_model'   # Folder with allesfitter output
instrument = 'Raphael'                      # Instrument name used in allesfitter data
data_key = 'rv'                            # Key for radial velocity data

# Load allesfitter results
alles = allesfitter.allesclass(data_folder)

# Extract observation times and measured radial velocities (RV)
times = alles.data[instrument]['time']         # Observation time array
rv_measurements = alles.data[instrument][data_key]  # Measured RV values

# Calculate RV errors by combining white noise and jitter terms in quadrature
rv_white_noise = alles.data[instrument]['white_noise_' + data_key]
rv_jitter = alles.posterior_params_median['jitter_' + data_key + '_' + instrument]
rv_errors = np.sqrt(rv_white_noise**2 + rv_jitter**2)

# Get the baseline RV (e.g. systemic velocity)
rv_baseline = alles.get_posterior_median_baseline(instrument, data_key)

# Get the median model RV values at observed times
rv_model_at_data = alles.get_posterior_median_model(instrument, data_key)

# Calculate residuals between observed RV and model + baseline
rv_residuals = rv_measurements - (rv_model_at_data + rv_baseline)

print("Number of residuals:", len(rv_residuals))

# Create a finely spaced time grid to evaluate the model smoothly
time_grid = np.linspace(2459113.834178, 2459861.987281, 1000)

# Evaluate the model and baseline on the new time grid
model_on_grid = alles.get_posterior_median_model(instrument, data_key, time_grid)
baseline_on_grid = alles.get_posterior_median_baseline(instrument, data_key, xx=time_grid)

# Prepare dictionary of results to save in a table
output_data = {
    'time': times,
    'baseline': baseline_on_grid,
    'rv_observed': rv_measurements,
    'model_time': time_grid,
    'model_rv': model_on_grid,
    'residuals': rv_residuals
}

# Convert dictionary to a pandas DataFrame (transpose so columns are keys)
df = pd.DataFrame.from_dict(output_data, orient='index').T

# Save DataFrame to CSV file
df.to_csv('RV_model.csv', index=False)

