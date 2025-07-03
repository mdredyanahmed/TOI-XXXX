import pandas as pd
import matplotlib.pyplot as plt

# === Light background plot ===
#plt.style.use('dark_background')  # Clean white background, no gridlines

# File paths for Sonora models
model_files = {
    '[M/H] = 0.0': 'nc+0.0_co1.0_age_bobcat_cleaned.csv',
    '[M/H] = -0.5': 'nc-0.5_co1.0_age_bobcat_cleaned.csv',
    '[M/H] = +0.5': 'nc+0.5_co1.0_age_bobcat_cleaned.csv',
}

# File paths for Baraffe models
baraffe_files = {
    'Age = 1 Gyr': 'baraffe_1gyr.csv',

    'Age =  3 Gyr': 'baraffe_3gyr_interpolated.csv',
        'Age =  5 Gyr': 'baraffe_5gyr.csv'
}

# Read Sonora models
model_dfs = {}
for label, filepath in model_files.items():
    df = pd.read_csv(filepath)
    df = df.apply(pd.to_numeric, errors='coerce').dropna()
    model_dfs[label] = df

# Read Baraffe models
baraffe_dfs = {}
for label, filepath in baraffe_files.items():
    df = pd.read_csv(filepath)
    df = df.apply(pd.to_numeric, errors='coerce').dropna()
    baraffe_dfs[label] = df

# Read brown dwarf data
browndwarf_file = 'transiting_brown_dwarf_list_march2025.txt'
browndwarf_df = pd.read_csv(browndwarf_file, delim_whitespace=True)
browndwarf_mass = browndwarf_df['Mj']
browndwarf_radius = browndwarf_df['Rj']
browndwarf_mass_err = browndwarf_df['Mj_err1']
browndwarf_radius_err = browndwarf_df['Rj_err1']

# Separate brown dwarfs and low mass stars
bd_mask = browndwarf_mass <= 84
star_mask = browndwarf_mass > 84

# Target ages for Sonora models
target_ages = [3, 4]

# Linestyles and colors for Sonora
metallicity_styles = {
    '[M/H] = 0.0': '-',
    '[M/H] = -0.5': '--',
    '[M/H] = +0.5': ':'
}
custom_colors = {
    ('[M/H] = 0.0', 4): '#E69F00',
    ('[M/H] = 0.0', 3): '#56B4E9',
    ('[M/H] = -0.5', 3): '#D55E00',
    ('[M/H] = +0.5', 3): '#CC79A7',
}

# Linestyles and colors for Baraffe models
baraffe_styles = {
    'Age = 1 Gyr': ':',
  
    'Age =  3 Gyr': '-',
        'Age =  5 Gyr': '--'
}
baraffe_colors = {
    'Age = 1 Gyr': 'black',    # Use basic green
   # Use orange
    'Age =  3 Gyr': 'green' ,
        'Age =  5 Gyr': 'orange',     # Use red
}
# Plot
fig, ax = plt.subplots(figsize=(14, 12))

# Plot Sonora models
#for label, df in model_dfs.items():
#    for target_age in target_ages:
#        if label != '[M/H] = 0.0' and target_age != 3:
#            continue
#        nearest_age = df['age(Gyr)'].sub(target_age).abs().idxmin()
#        age_value = df.loc[nearest_age, 'age(Gyr)']
#        df_nearest = df[df['age(Gyr)'] == age_value]
#        linestyle = metallicity_styles[label]
#        color = custom_colors.get((label, target_age), 'black')
#        ax.plot(df_nearest['M/MSun'] * 1047.56, df_nearest['R/Rsun'] * 9.731,
#                linestyle=linestyle, color=color,
#                label=f'{label}, Age = {age_value:.2f} Gyr, Sonora Bobcat Model (2021) ', lw=2)

# Plot Baraffe models
for label, df in baraffe_dfs.items():
    ax.plot(df['Mass'] * 1047.56, df['Radius'] * 9.731,
            linestyle=baraffe_styles[label], color=baraffe_colors[label],
            label=f'{label} (Baraffe+2003)', lw=2)

# Plot brown dwarfs
ax.errorbar(browndwarf_mass[bd_mask], browndwarf_radius[bd_mask],
            xerr=browndwarf_mass_err[bd_mask], yerr=browndwarf_radius_err[bd_mask],
            fmt='o', color='blue', markersize=6,
            ecolor='blue', elinewidth=0.7, capsize=2,
            label='Transiting Brown Dwarfs')

# Plot low mass stars with different marker
ax.errorbar(browndwarf_mass[star_mask], browndwarf_radius[star_mask],
            xerr=browndwarf_mass_err[star_mask], yerr=browndwarf_radius_err[star_mask],
            fmt='*', color='#8B006A', markersize=12,
            ecolor='#8B006A', elinewidth=0.7, capsize=2,
            label='Low-Mass Stars')

# Plot TOI 2155 b
planet_mass = 80.60
planet_radius = 0.970
mass_err = 0.928
radius_err_lower = [0.020]
radius_err_upper = [0.018]
yerr = [radius_err_lower, radius_err_upper]

ax.errorbar(planet_mass, planet_radius,
            xerr=mass_err, yerr=yerr,
            fmt='o', color='#d62728', markersize=6,
            ecolor='#d62728', elinewidth=1.2, capsize=2,
            label='TOI 2155b')

# Region separators
ax.axvline(x=12, color='red', linestyle='--')
#ax.axvline(x=42, color='red', linestyle='--')
ax.axvline(x=84, color='red', linestyle='--')

#ax.text(13, 1.450, 'Low-mass Brown dwarfs', fontsize=13, fontweight='bold')
ax.text(45, 1.450, 'Brown dwarfs', fontsize=20, fontweight='bold')
ax.text(86, 1.450, 'Low Mass Stars', fontsize=20, fontweight='bold')

# Axis labels and formatting
ax.set_xlabel('Mass [$M_\\mathrm{J}$]', fontsize=20)
ax.set_ylabel('Radius [$R_\\mathrm{J}$]', fontsize=20)
ax.set_xlim(8, 105)
ax.set_ylim(0.5, 1.8)
ax.tick_params(axis='both', which='major', labelsize=11, direction='out', length=6, width=1)
ax.tick_params(axis='both', which='minor', direction='out', length=3, width=0.8)
ax.minorticks_on()

# Legend
legend = ax.legend(fontsize=11, frameon=True, loc='best')
legend.get_frame().set_edgecolor('gray')
legend.get_frame().set_alpha(0.6)

plt.tight_layout()
plt.savefig('TOI2155_evolution.png', format='png', dpi=400, bbox_inches='tight')
plt.show()






#import pandas as pd
#import matplotlib.pyplot as plt
#
#plt.style.use('dark_background')
## === Light background plot ===
##plt.style.use('seaborn-v0_8-white')  # Clean white background, no gridlines
#
## File paths for Sonora models
#model_files = {
#    '[M/H] = 0.0': 'nc+0.0_co1.0_age_bobcat_cleaned.csv',
#    '[M/H] = -0.5': 'nc-0.5_co1.0_age_bobcat_cleaned.csv',
#    '[M/H] = +0.5': 'nc+0.5_co1.0_age_bobcat_cleaned.csv',
#}
#
## File paths for Baraffe models
#baraffe_files = {
#    'Age = 1 Gyr': 'baraffe_1gyr.csv',
#    'Age =  5 Gyr': 'baraffe_5gyr.csv',
#    # 'Baraffe 10 Gyr': 'baraffe_10gyr.csv',
#    'Age =  3 Gyr':'baraffe_3gyr_interpolated.csv'
#}
#
## Read Sonora models
#model_dfs = {}
#for label, filepath in model_files.items():
#    df = pd.read_csv(filepath)
#    df = df.apply(pd.to_numeric, errors='coerce').dropna()
#    model_dfs[label] = df
#
## Read Baraffe models
#baraffe_dfs = {}
#for label, filepath in baraffe_files.items():
#    df = pd.read_csv(filepath)
#    df = df.apply(pd.to_numeric, errors='coerce').dropna()
#    baraffe_dfs[label] = df
#
## Read brown dwarf data
#browndwarf_file = 'transiting_brown_dwarf_list_march2025.txt'
#browndwarf_df = pd.read_csv(browndwarf_file, delim_whitespace=True)
#browndwarf_mass = browndwarf_df['Mj']
#print(browndwarf_mass)
#browndwarf_radius = browndwarf_df['Rj']
#browndwarf_mass_err = browndwarf_df['Mj_err1']  # assuming this is the upper error
#browndwarf_radius_err = browndwarf_df['Rj_err1']  # assuming this is the upper error
#
## Target ages for Sonora models
#target_ages = [3, 4, 13]
#
## Linestyles and colors for Sonora
#metallicity_styles = {
#    '[M/H] = 0.0': '-',
#    '[M/H] = -0.5': '--',
#    '[M/H] = +0.5': ':'
#}
#custom_colors = {
#    ('[M/H] = 0.0', 4): '#E69F00',
#    ('[M/H] = 0.0', 3): '#56B4E9',
#    ('[M/H] = 0.0', 13): '#ff0000',
#    ('[M/H] = -0.5', 3): '#D55E00',
#    ('[M/H] = +0.5', 3): '#CC79A7',
#}
#
## Linestyles and colors for Baraffe models
#baraffe_styles = {
#    'Age = 1 Gyr': ':',
#    'Age =  5 Gyr': '--',
#    'Age =  3 Gyr': '-'
#}
#baraffe_colors = {
#    'Age = 1 Gyr': '#9467bd',  # Purple
#    'Age =  5 Gyr': '#8c564b',  # Brown
# 'Age =  3 Gyr': '#009E73'  # Bright Red
#
#}
#
## Plot
#fig, ax = plt.subplots(figsize=(8.5, 6.5))
#
## Plot Sonora models
#for label, df in model_dfs.items():
#    for target_age in target_ages:
#        if label != '[M/H] = 0.0' and target_age != 3:
#            continue
#        nearest_age = df['age(Gyr)'].sub(target_age).abs().idxmin()
#        age_value = df.loc[nearest_age, 'age(Gyr)']
#        df_nearest = df[df['age(Gyr)'] == age_value]
#        linestyle = metallicity_styles[label]
#        color = custom_colors.get((label, target_age), 'black')
#        ax.plot(df_nearest['M/MSun'] * 1047.56, df_nearest['R/Rsun'] * 9.731,
#                linestyle=linestyle, color=color,
#                label=f'{label}, Age = {age_value:.2f} Gyr, Sonora Bobcat Model', lw=2)
#
## Plot Baraffe models
#for label, df in baraffe_dfs.items():
#    ax.plot(df['Mass'] * 1047.56, df['Radius'] * 9.731,  # assuming Mass in Msun, Radius in Rsun
#            linestyle=baraffe_styles[label], color=baraffe_colors[label],
#            label=f'{label} (Baraffe 2003)', lw=2)
#
## Plot transiting brown dwarfs
#ax.errorbar(browndwarf_mass, browndwarf_radius,
#            xerr=browndwarf_mass_err, yerr=browndwarf_radius_err,
#            fmt='o', color='gray', markersize=3,
#            ecolor='gray', elinewidth=0.7, capsize=2,
#            label='Transiting brown dwarfs')
#
## Plot TOI 2155 b
#planet_mass = 80.60
#planet_radius = 0.970
#mass_err = 0.928
## radius_err = 0.018
#radius_err_lower = [0.020]  # lower uncertainty
#radius_err_upper = [0.018]
#yerr = [radius_err_lower, radius_err_upper]
#
#ax.errorbar(planet_mass, planet_radius,
#            xerr=mass_err, yerr=yerr,
#            fmt='*', color='#d62728', markersize=12,
#            ecolor='#d62728', elinewidth=1.2, capsize=2,
#            label='TOI 2155b')
#
## Axis labels and formatting
#ax.set_xlabel('Mass [$M_\\mathrm{J}$]', fontsize=14)
#ax.set_ylabel('Radius [$R_\\mathrm{J}$]', fontsize=14)
#ax.set_xlim(8, 90)
#ax.set_ylim(0.5, 1.650)
#ax.tick_params(axis='both', which='major', labelsize=11, direction='out', length=6, width=1)
#ax.tick_params(axis='both', which='minor', direction='out', length=3, width=0.8)
#ax.minorticks_on()
#
## Legend
#legend = ax.legend(fontsize=7.5, frameon=True, loc='lower left')
#legend.get_frame().set_edgecolor('gray')
#legend.get_frame().set_alpha(0.6)
#
#plt.tight_layout()
#plt.savefig("bd_evolution_dark.png",dpi=400)
#plt.show()
#
