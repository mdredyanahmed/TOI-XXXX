import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ---------------------------
# Plot Aesthetics
# ---------------------------
plt.style.use('seaborn-v0_8-white')
params = {
    'legend.fontsize': 12,
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'grid.color': 'k',
    'grid.linestyle': ':',
    'grid.linewidth': 0.5,
    'mathtext.fontset': 'stix',
    'mathtext.rm': 'DejaVu serif',
    'font.family': 'DejaVu serif',
    'font.serif': 'Times New Roman',
}
rcParams.update(params)

# ---------------------------
# Load Data
# ---------------------------
browndwarf_file = 'browndwarf_list_with_density_asymmetric.csv'
df = pd.read_csv(browndwarf_file)

# Extract columns
mass = df['Mj']
mass_err_upper = df['Mj_err1']
mass_err_lower = df['Mj_err2']
density = df['Density_cgs']
density_err_upper = df['Density_err_upper_cgs']
density_err_lower = df['Density_err_lower_cgs']
radius = df['Rj']
temperature = df['Teff']

# ---------------------------
# Plotting
# ---------------------------
fig, ax = plt.subplots(figsize=(14, 10))

# Main scatter plot (Brown dwarfs)
scatter = ax.scatter(
    mass, density,
    edgecolors='none',
    s=radius * 25,  # size scaled by radius
    c=temperature,
    cmap='jet',
    alpha=0.9
)

# Colorbar
clb = plt.colorbar(scatter, ax=ax, orientation='vertical')
clb.set_label('Temperature (K)', labelpad=15, rotation=270, fontsize=20)

# Asymmetric error bars for brown dwarfs
ax.errorbar(
    mass,
    density,
    xerr=[mass_err_lower, mass_err_upper],
    yerr=[density_err_lower, density_err_upper],
    linestyle='None',
    color='gray',
    alpha=0.7,
    capsize=5
)

# Plot low-mass stars (masses > 84 Mj)
star_mask = mass > 84

# Error bars for stars
ax.errorbar(
    mass[star_mask], density[star_mask],
    xerr=[mass_err_lower[star_mask], mass_err_upper[star_mask]],
    yerr=[density_err_lower[star_mask], density_err_upper[star_mask]],
    linestyle='None',
    ecolor='gray',
    elinewidth=0.7,
    capsize=2,
    alpha=0.7
)

# Scatter stars with color following temperature
scatter_stars = ax.scatter(
    mass[star_mask], density[star_mask],
    marker='*',
    c=temperature[star_mask],
    cmap='jet',
    s=160,  # Adjust star size
    edgecolors='black',
    linewidths=0.7,
    label='Low-Mass Stars'
)

# Labels and limits
ax.set_xlabel('Mass [$M_{J}$]', fontsize=20)
ax.set_ylabel('Density [g cm$^{-3}$]', fontsize=20)
ax.set_xlim(0, 120)
ax.set_ylim(-30, max(density + density_err_upper) * 1.1)

# Reference vertical lines
ax.axvline(x=12, color='red', linestyle='--')
# ax.axvline(x=42, color='red', linestyle='--')
ax.axvline(x=84, color='red', linestyle='--')
# ax.text(13, max(density) + 45, 'Low-mass Brown dwarfs', fontsize=13, fontweight='bold')
ax.text(45, max(density) + 45, ' Brown dwarfs', fontsize=20, fontweight='bold')
ax.text(86, max(density) + 45, 'Low Mass Stars', fontsize=20, fontweight='bold')

# Annotate TOI-2155b if present
if 'TOI-2155b' in df['Name'].values:
    idx = df[df['Name'] == 'TOI-2155b'].index[0]
    toi_color = plt.cm.jet((temperature[idx] - min(temperature)) / (max(temperature) - min(temperature)))

    ax.scatter(
        mass[idx], density[idx],
        s=radius[idx] * 25,
        color=toi_color,
        marker='o',
        linewidths=1.2,
        zorder=5
    )

    ax.scatter(
        mass[idx], density[idx],
        s=radius[idx] * 20,
        color='red',
        marker='x',
        linewidths=1.5,
        zorder=6,
        label='TOI-2155b'
    )

    ax.annotate(
        '',
        xy=(mass[idx], density[idx]),
        xytext=(mass[idx] + 4, density[idx] + (max(density) * 0.15)),
        arrowprops=dict(arrowstyle='->', color='red', linewidth=2)
    )
    ax.text(mass[idx], density[idx] + (max(density) * 0.1), 'TOI-2155b', fontsize=20, fontweight='bold', ha='center')

# Radius size legend
desired_radii = np.array([0.5, 0.7, 0.9, 1.0, 1.5, 2.0, 3.0])
desired_sizes = desired_radii * 40

handles = [
    plt.scatter([], [], s=size, edgecolor='none', c='gray', alpha=0.9)
    for size in desired_sizes
]
labels = [f'{r} $R_J$' for r in desired_radii]










ax.legend(handles, labels, loc="upper left", title="Radius [$R_{J}$]", fontsize=20)



plt.tight_layout()
plt.savefig('mass_density_plot.png', dpi=400, bbox_inches='tight')
plt.show()



#import os
#import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt
#import matplotlib.ticker as ticker
#import seaborn as sns
#from matplotlib import rcParams
#plt.style.use('dark_background')
## ---------------------------
## Set Plotting Aesthetics
## ---------------------------
#params = {
#    'legend.fontsize': 12,
#    'axes.labelsize': 12,
#    'axes.titlesize': 12,
#    'xtick.labelsize': 12,
#    'ytick.labelsize': 12,
#    'grid.color': 'k',
#    'grid.linestyle': ':',
#    'grid.linewidth': 0.5,
#    'mathtext.fontset': 'stix',
#    'mathtext.rm': 'DejaVu serif',
#    'font.family': 'DejaVu serif',
#    'font.serif': 'Times New Roman',
#}
#rcParams.update(params)
#
## ---------------------------
## Load Data
## ---------------------------
#data = pd.read_csv("browndwarf_final_earlier.csv")
#
## ---------------------------
## Calculate Density Uncertainty
## ---------------------------
#mass = data['mass_intermsof_jupitermass']
#radius = data['radius_intermsof_jupiterradius']
#mass_unc = data['massuncertainity_intermsof_jupitermass']
#radius_unc = data['radiusuncertainity_intermsof_jupiterradius']
#density = data['Density']
#temp = data['TEMP']
#
## Propagate error: fractional uncertainties added, radius cubed → multiply radius term by 3
#density_unc = ((mass_unc / mass) + 3 * (radius_unc / radius)) * density
#
## Append density uncertainty to DataFrame
#data['density_uncertainity'] = density_unc
#
## ---------------------------
## Prepare Arrays for Plotting
## ---------------------------
#m = mass.to_numpy()
#d = density.to_numpy()
#r = radius.to_numpy()
#t = temp.to_numpy()
#
## ---------------------------
## Begin Plotting
## ---------------------------
#fig, ax = plt.subplots(figsize=(10, 8))
#
## Main scatter plot with temperature as color and radius as size
#scatter = ax.scatter(
#    m, d,
#    edgecolors='none',
#    s=r * 40,            # Scale point size with radius
#    c=t,                 # Color by temperature
#    cmap='jet',
#    alpha=0.9
#)
#
## Colorbar setup
#clb = plt.colorbar(scatter, ax=ax, orientation='vertical')
#clb.set_label('Temperature (K)', labelpad=15, rotation=270)
#
## Add error bars
#ax.errorbar(
#    x=m,
#    y=d,
#    yerr=data['density_uncertainity'],
#    xerr=mass_unc,
#    linestyle='None',
#    color='gray',
#    alpha=0.4,
#    capsize=5
#)
#
## Axis labels and limits
#ax.set_xlabel('Mass [$M_{J}$]')
#ax.set_ylabel('Density [g cm$^{-3}$]')
#ax.set_xlim(0, 120)
#ax.set_ylim(-30, max(d) + 70)  # Add padding above data
#
## ---------------------------
## Reference Lines and Annotations
## ---------------------------
#ax.axvline(x=94, color='red', linestyle='--')  # Division between BDs and low-mass stars
#ax.text(10, 125, 'Brown dwarfs', fontsize=10, fontweight='bold')
#ax.text(94.5, 75, 'Low mass stars', fontsize=10, fontweight='bold')
#
## Annotate TOI-2155b with arrow and label
#ax.annotate(
#    '',
#    xy=(81, 111),       # Arrowhead
#    xytext=(90, 150),   # Start of arrow
#    arrowprops=dict(arrowstyle='->', color='red', linewidth=2)
#)
#ax.text(81, 127, 'TOI-2155b', fontsize=8, fontweight='bold', ha='center')
#
## ---------------------------
## Custom Radius Size Legend
## ---------------------------
#desired_radii = np.array([0.5, 0.7, 0.8, 0.9, 1, 2, 3])  # Radii in Jupiter units
## Convert to plot size scale (based on original scaling)
#desired_sizes = (desired_radii / 11.21) * 490
#
## Create dummy handles for legend
#handles = [
#    plt.scatter([], [], s=size, edgecolor='none', c='gray', alpha=0.6)
#    for size in desired_sizes
#]
#labels = [f'{radius} $R_J$' for radius in desired_radii]
#
## Draw legend on plot
#legend = ax.legend(handles, labels, loc="upper right", title="Radius [$R_{J}$]", fontsize=12)
#
## ---------------------------
## Save and Show
## ---------------------------
#plt.savefig('mass_density_plot.png', format='png', dpi=400, bbox_inches='tight')
#plt.show()
#
