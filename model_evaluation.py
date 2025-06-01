import csv


import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
import allesfitter 
datadir = 'allesfit_eccentric_model' #change this to what you need 
inst = 'TESS' #change this to 
#what you need 
key = 'flux' #change this to what you need

#::: initialize the allesclass
# === Initialize the allesfitter class with the specified directory ===
alles = allesfitter.allesclass(datadir)

# === Load the raw TESS time-series data ===
# This assumes 'TESS.csv' contains at least 1 column (time)
data = np.loadtxt(f"{datadir}/TESS.csv", delimiter=",")
time = data[:, 0]  # Extract time from the first column
print("Time shape for the run is", time.shape)

# === Retrieve the best-fit model and baseline from the posterior ===
# `xx=time` means we want the model evaluated at these specific time points
model = alles.get_one_posterior_model(inst, key, xx=time)
baseline = alles.get_one_posterior_baseline(inst, key, xx=time)
print("Model shape for the run is", model.shape)


# Create a dictionary of arrays
data = { 'baseline':baseline,'model': model,'time':time}

# Create a DataFrame from the dictionary
df = pd.DataFrame(data)

# Save the DataFrame to a CSV file
df.to_csv('model_sector_all_sector_TESS.csv', index=False)

