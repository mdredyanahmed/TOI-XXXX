import pandas as pd
import allesfitter
#::: your settings
datadir = 'allesfit_eccentric_model' #change this to what you need
inst = 'TESS' #change this to what you need
key = 'flux' #change this to what you need

#::: initialize the allesclass
alles = allesfitter.allesclass(datadir)

parameters = alles.posterior_params
df = pd.DataFrame(parameters)

print (df.head())
# Save to CSV: each parameter is a column, each row is a posterior sample
df.to_csv('posterior_samples_from_posterior_params.csv', index=False)

print(f"Saved {len(df)} samples to 'posterior_samples_from_posterior_params.csv'")


