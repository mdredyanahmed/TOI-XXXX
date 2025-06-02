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
alles = allesfitter.allesclass(datadir)


parameters = alles.posterior_params
# Define additional parameters to include in the dataset
additional_parameters = ['b_rr', 'b_rsuma', 'b_cosi', 'b_epoch', 'host_ldc_q1_TESS', 'host_ldc_q2_TESS', 
                         'ln_err_flux_TESS', 'ln_jitter_rv_Raphael', 'baseline_offset_rv_Raphael']

# Update the parameters dictionary by adding the new ones
parameters.update(dict(zip(additional_parameters, list(parameters.values()))))

# Prepare the data to be saved
keys = list(parameters.keys())
values = list(parameters.values())
X = np.column_stack((keys, values))
X=X.T
# Save to CSV file
np.savetxt('allelclass_eccentric_values.csv', X, delimiter=',', fmt='%s')
