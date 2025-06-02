import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('seaborn-v0_8-white')  # Clean white background, no gridlines

# File paths for different metallicities
model_files = {
    '[M/H] = 0.0': 'evolution/nc_m0.0_age',
    '[M/H] = -0.5': 'evolution/nc_m-0.5_age',
    '[M/H] = +0.5': 'evolution/nc_m+0.5_age',
}
browndwarf_file = 'browndwarf_final.csv'

# Read and process model data files
model_dfs = {}
for label, filepath in model_files.items():
    df = pd.read_csv(filepath, sep=r'\s+', header=None)
    df.columns = df.iloc[1]
    df = df[2:].copy().reset_index(drop=True)
    df = df.apply(pd.to_numeric, errors='coerce').dropna()
    df.columns = ['age(Gyr)', 'M/MSun', 'logL/Lsun', 'Teff(K)', 'logg', 'R/Rsun']
    model_dfs[label] = df

# Read brown dwarf data
browndwarf_df = pd.read_csv(browndwarf_file)
browndwarf_mass = browndwarf_df['mass_intermsof_jupitermass']
browndwarf_radius = browndwarf_df['radius_intermsof_jupiterradius']
browndwarf_mass_err = browndwarf_df['massuncertainity_intermsof_jupitermass']
browndwarf_radius_err = browndwarf_df['radiusuncertainity_intermsof_jupiterradius']
bd_names = browndwarf_df['Name']

# Target ages
target_ages = [3.4, 0.49, 13]

# Linestyles by metallicity
metallicity_styles = {
    '[M/H] = 0.0': '-',
    '[M/H] = -0.5': '--',
    '[M/H] = +0.5': ':'
}

# Custom color mapping
custom_colors = {
    ('[M/H] = 0.0', 0.49): 'black',
    ('[M/H] = 0.0', 3.4): '#228B22',
    ('[M/H] = 0.0', 13): '#1f77b4',
    ('[M/H] = -0.5', 3.4): '#FF8C00',
    ('[M/H] = +0.5', 3.4): '#800080',
}

# Create plot
fig, ax = plt.subplots(figsize=(8.5, 6.5))

# Plot model curves
for label, df in model_dfs.items():
    for target_age in target_ages:
        if label != '[M/H] = 0.0' and target_age != 3.4:
            continue
        nearest_age = df['age(Gyr)'].sub(target_age).abs().idxmin()
        age_value = df.loc[nearest_age, 'age(Gyr)']
        df_nearest = df[df['age(Gyr)'] == age_value]
        linestyle = metallicity_styles[label]
        color = custom_colors[(label, target_age)]
        ax.plot(df_nearest['M/MSun'] * 1047.56, df_nearest['R/Rsun'],
                linestyle=linestyle, color=color,
                label=f'{label}, Age = {age_value:.2f} Gyr, Sonora Models (Diamondback)', lw=2)

# Brown dwarfs
ax.errorbar(browndwarf_mass, browndwarf_radius,
            xerr=browndwarf_mass_err, yerr=browndwarf_radius_err,
            fmt='o', color='gray', markersize=3,
            ecolor='gray', elinewidth=0.7, capsize=2,
            label='Transiting brown dwarfs')

# TOI 2155 b
# Plot the planet with error bars
planet_mass = 81
planet_radius = 0.94
mass_err = 0.951 
radius_err = 0.021



ax.errorbar(planet_mass, planet_radius,
            xerr=mass_err, yerr=radius_err,
            fmt='*', color='#d62728', markersize=12,
            ecolor='#d62728', elinewidth=1.2, capsize=2,
            label='TOI 2155b')

# Axis labels and formatting
ax.set_xlabel('Mass [$M_\\mathrm{J}$]', fontsize=12)
ax.set_ylabel('Radius [$R_\\mathrm{J}$]', fontsize=12)
ax.set_xlim(8, 90)
ax.set_ylim(0.5, 1.650)
ax.tick_params(axis='both', which='major', labelsize=11, direction='out', length=6, width=1)
ax.tick_params(axis='both', which='minor', direction='out', length=3, width=0.8)
ax.minorticks_on()

# Legend
legend = ax.legend(fontsize=8, frameon=True, loc='best')
legend.get_frame().set_edgecolor('gray')
legend.get_frame().set_alpha(0.6)

plt.tight_layout()
plt.savefig('TOI2155_evolution.png', dpi=400)
plt.show()

 #Dark background for better contrast
#import pandas as pd
#import matplotlib.pyplot as plt
#
#plt.style.use('dark_background')  # Dark background for better contrast
#
## File paths for different metallicities
#model_files = {
#    '[M/H] = 0.0': 'evolution/nc_m0.0_age',
#    '[M/H] = -0.5': 'evolution/nc_m-0.5_age',
#    '[M/H] = +0.5': 'evolution/nc_m+0.5_age',
#}
#browndwarf_file = 'browndwarf_final.csv'
#
## Read and process model data files
#model_dfs = {}
#for label, filepath in model_files.items():
#    df = pd.read_csv(filepath, sep=r'\s+', header=None)
#    df.columns = df.iloc[1]
#    df = df[2:].copy().reset_index(drop=True)
#    df = df.apply(pd.to_numeric, errors='coerce').dropna()
#    df.columns = ['age(Gyr)', 'M/MSun', 'logL/Lsun', 'Teff(K)', 'logg', 'R/Rsun']
#    model_dfs[label] = df
#
## Read brown dwarf data
#browndwarf_df = pd.read_csv(browndwarf_file)
#browndwarf_mass = browndwarf_df['mass_intermsof_jupitermass']
#browndwarf_radius = browndwarf_df['radius_intermsof_jupiterradius']
#browndwarf_mass_err = browndwarf_df['massuncertainity_intermsof_jupitermass']
#browndwarf_radius_err = browndwarf_df['radiusuncertainity_intermsof_jupiterradius']
#bd_names = browndwarf_df['Name']
#
## Target ages
#target_ages = [3.4, 0.49, 13]
#
## Linestyles by metallicity
#metallicity_styles = {
#    '[M/H] = 0.0': '-',
#    '[M/H] = -0.5': '--',
#    '[M/H] = +0.5': ':'
#}
#
## Custom color mapping (vivid, colorblind-friendly)
#custom_colors = {
#    ('[M/H] = 0.0', 0.49): '#FFDD57',  # bright yellow
#    ('[M/H] = 0.0', 3.4): '#00CC96',  # teal green
#    ('[M/H] = 0.0', 13): '#19D3F3',   # bright cyan
#    ('[M/H] = -0.5', 3.4): '#FFA15A', # orange
#    ('[M/H] = +0.5', 3.4): '#AB63FA', # purple
#}
#
## Create plot
#fig, ax = plt.subplots(figsize=(8.5, 6.5))
#
## Plot model curves
#for label, df in model_dfs.items():
#    for target_age in target_ages:
#        if label != '[M/H] = 0.0' and target_age != 3.4:
#            continue
#        nearest_age = df['age(Gyr)'].sub(target_age).abs().idxmin()
#        age_value = df.loc[nearest_age, 'age(Gyr)']
#        df_nearest = df[df['age(Gyr)'] == age_value]
#        linestyle = metallicity_styles[label]
#        color = custom_colors[(label, target_age)]
#        ax.plot(df_nearest['M/MSun'] * 1047.56, df_nearest['R/Rsun'],
#                linestyle=linestyle, color=color,
#                label=f'{label}, Age = {age_value:.2f} Gyr, Sonora Models (Diamondback)', lw=2.2)
#
## Brown dwarfs
#ax.errorbar(browndwarf_mass, browndwarf_radius,
#            xerr=browndwarf_mass_err, yerr=browndwarf_radius_err,
#            fmt='o', color='lightgray', markersize=3.5,
#            ecolor='lightgray', elinewidth=0.7, capsize=2,
#            label='Transiting brown dwarfs')
#
## TOI 2155 b
#planet_mass = 81
#planet_radius = 0.94
#mass_err = 0.951 
#radius_err = 0.021
#
#ax.errorbar(planet_mass, planet_radius,
#            xerr=mass_err, yerr=radius_err,
#            fmt='*', color='#EF553B', markersize=13,
#            ecolor='#EF553B', elinewidth=1.5, capsize=2,
#            label='TOI 2155b')
#
## Axis labels and formatting
#ax.set_xlabel('Mass [$M_\\mathrm{J}$]', fontsize=12)
#ax.set_ylabel('Radius [$R_\\mathrm{J}$]', fontsize=12)
#ax.set_xlim(8, 90)
#ax.set_ylim(0.5, 1.650)
#ax.tick_params(axis='both', which='major', labelsize=11, direction='out', length=6, width=1.1, colors='white')
#ax.tick_params(axis='both', which='minor', direction='out', length=3, width=0.8, colors='white')
#ax.minorticks_on()
#
## Legend
#legend = ax.legend(fontsize=8.5, frameon=True, loc='best')
#legend.get_frame().set_edgecolor('white')
#legend.get_frame().set_alpha(0.5)
#
#plt.tight_layout()
#plt.savefig('TOI2155_evolution_dark.png', dpi=400)
#plt.show()
#
