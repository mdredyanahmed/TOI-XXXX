
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import Normalize

# ---------------------------
# Plot Aesthetics
# ---------------------------
plt.style.use('seaborn-v0_8-white')
params = {
    'legend.fontsize': 12,
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
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
df = pd.read_csv('browndwarf_list_with_density_asymmetric.csv')

mass = df['Mj']
mass_err_upper = df['Mj_err1']
mass_err_lower = df['Mj_err2']
density = df['Density_cgs']
density_err_upper = df['Density_err_upper_cgs']
density_err_lower = df['Density_err_lower_cgs']
radius = df['Rj']
temperature = df['Teff']

# Classification masks
star_mask = mass >= 82
bd_mask = mass <= 81.7

# Select densities for BDs
bd_density = density[bd_mask]

# Mean and median density
bd_density_mean = np.mean(bd_density)
bd_density_median = np.median(bd_density)
print("Mean density for BDs:", bd_density_mean)
print("Median density for BDs:", bd_density_median)

# ---------------------------
# Normalization for Temperature Colors
# ---------------------------
norm = Normalize(vmin=3000, vmax=9000)
cmap = 'jet'

# ---------------------------
# Plotting
# ---------------------------
fig, ax = plt.subplots(figsize=(12, 8))

# Brown dwarfs
scatter_bd = ax.scatter(
    mass[bd_mask], density[bd_mask],
    edgecolors='none',
    s=radius[bd_mask] * 80,
    c=temperature[bd_mask],
    cmap=cmap,
    norm=norm,
    alpha=0.9,
    label='Brown Dwarfs'
)

# Error bars for brown dwarfs
ax.errorbar(
    mass[bd_mask], density[bd_mask],
    xerr=[mass_err_lower[bd_mask], mass_err_upper[bd_mask]],
    yerr=[density_err_lower[bd_mask], density_err_upper[bd_mask]],
    linestyle='None',
    color='gray',
    alpha=0.7,
    capsize=5
)

# Low-mass stars as star markers
ax.scatter(
    mass[star_mask], density[star_mask],
    marker='*',
    c=temperature[star_mask],
    cmap=cmap,
    norm=norm,
    s=160,
    edgecolors='black',
    linewidths=0.7,
    label='Low-Mass Stars'
)

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

# Colorbar
clb = plt.colorbar(scatter_bd, ax=ax, orientation='vertical')
clb.set_label('Temperature (K)', labelpad=15, fontsize=18)

# Labels and limits
ax.set_xlabel('Mass [$M_{J}$]', fontsize=25)
ax.set_ylabel('Density [g cm$^{-3}$]', fontsize=25)
ax.set_xlim(0, 120)
ax.set_ylim(-30, max(density + density_err_upper) * 1.1)

# Vertical classification lines
ax.axvline(x=12, color='red', linestyle='--')
ax.axvline(x=81.7, color='red', linestyle='--')
ax.text(45, max(density) + 45, 'Brown Dwarfs', fontsize=18, fontweight='bold')
ax.text(86, max(density) + 45, 'Low-Mass Stars', fontsize=18, fontweight='bold')

# Annotate TOI-2155b if present
if 'TOI-2155b' in df['Name'].values:
    idx = df[df['Name'] == 'TOI-2155b'].index[0]
    ax.scatter(
        mass[idx], density[idx],
        s=300,
        color='red',
        marker='x',
        linewidths=2,
        zorder=6,
        label='TOI-2155b'
    )
    ax.annotate(
        '',
        xy=(mass[idx], density[idx]),
        xytext=(mass[idx] + 4, density[idx] + (max(density) * 0.15)),
        arrowprops=dict(arrowstyle='->', color='blue', linewidth=2)
    )
    ax.text(mass[idx], density[idx] + (max(density) * 0.1), 'TOI-2155b', fontsize=20, fontweight='bold', ha='center')

# Radius size legend
desired_radii = np.array([0.5, 0.7, 0.9, 1.0, 1.5, 2.0, 3.0])
desired_sizes = desired_radii * 80
handles = [plt.scatter([], [], s=size, edgecolor='none', c='gray', alpha=0.9) for size in desired_sizes]
labels = [f'{r} $R_J$' for r in desired_radii]

ax.tick_params(axis='both', which='major', labelsize=22, direction='out', length=6, width=1)
ax.tick_params(axis='both', which='minor', labelsize=18, direction='out', length=3, width=0.8)
ax.legend(handles, labels, loc="upper left", title="Radius [$R_J$]", fontsize=16, title_fontsize=20)

plt.tight_layout()
plt.savefig('mass_density_plot_clean.png', dpi=600, bbox_inches='tight')
plt.show()

# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import rcParams
# from matplotlib.colors import Normalize

# # ---------------------------
# # Plot Aesthetics
# # ---------------------------
# plt.style.use('seaborn-v0_8-white')
# params = {
#     'legend.fontsize': 12,
#     'axes.labelsize': 12,
#     'axes.titlesize': 12,
#     'xtick.labelsize': 12,
#     'ytick.labelsize': 12,
#     'grid.color': 'k',
#     'grid.linestyle': ':',
#     'grid.linewidth': 0.5,
#     'mathtext.fontset': 'stix',
#     'mathtext.rm': 'DejaVu serif',
#     'font.family': 'DejaVu serif',
#     'font.serif': 'Times New Roman',
# }
# rcParams.update(params)

# # ---------------------------
# # Load Data
# # ---------------------------
# df = pd.read_csv('browndwarf_list_with_density_asymmetric.csv')

# mass = df['Mj']
# mass_err_upper = df['Mj_err1']
# mass_err_lower = df['Mj_err2']
# density = df['Density_cgs']
# density_err_upper = df['Density_err_upper_cgs']
# density_err_lower = df['Density_err_lower_cgs']
# radius = df['Rj']
# temperature = df['Teff']

# # Classification masks
# star_mask = mass > 84
# bd_mask = mass <= 84
# hot_bd_mask = bd_mask & (temperature >= 6000) & (temperature <= 7200)
# num_hot_bds = hot_bd_mask.sum()
# print(f"Number of brown dwarfs with F-type host stars (6000–7200 K): {num_hot_bds}")

# # ---------------------------
# # Normalization for Temperature Colors
# # ---------------------------
# norm = Normalize(vmin=3000, vmax=9000)
# cmap = 'jet'

# # ---------------------------
# # Plotting
# # ---------------------------
# fig, ax = plt.subplots(figsize=(14, 10))

# # Brown dwarfs (excluding hot BDs)
# scatter_bd = ax.scatter(
#     mass[bd_mask & ~hot_bd_mask], density[bd_mask & ~hot_bd_mask],
#     edgecolors='none',
#     s=radius[bd_mask & ~hot_bd_mask] * 60,
#     c=temperature[bd_mask & ~hot_bd_mask],
#     cmap=cmap,
#     norm=norm,
#     alpha=0.9,
#     label='Brown Dwarfs'
# )

# # Error bars for brown dwarfs
# ax.errorbar(
#     mass[bd_mask], density[bd_mask],
#     xerr=[mass_err_lower[bd_mask], mass_err_upper[bd_mask]],
#     yerr=[density_err_lower[bd_mask], density_err_upper[bd_mask]],
#     linestyle='None',
#     color='gray',
#     alpha=0.7,
#     capsize=5
# )

# # Error bars for stars
# ax.errorbar(
#     mass[star_mask], density[star_mask],
#     xerr=[mass_err_lower[star_mask], mass_err_upper[star_mask]],
#     yerr=[density_err_lower[star_mask], density_err_upper[star_mask]],
#     linestyle='None',
#     ecolor='gray',
#     elinewidth=0.7,
#     capsize=2,
#     alpha=0.7
# )

# # Low-mass stars as star markers
# ax.scatter(
#     mass[star_mask], density[star_mask],
#     marker='*',
#     c=temperature[star_mask],
#     cmap=cmap,
#     norm=norm,
#     s=160,
#     edgecolors='black',
#     linewidths=0.7,
#     label='Low-Mass Stars'
# )

# # Hot BDs (6000–7200 K) as triangles
# ax.scatter(
#     mass[hot_bd_mask], density[hot_bd_mask],
#     marker='^',
#     c=temperature[hot_bd_mask],
#     cmap=cmap,
#     norm=norm,
#     s=radius[hot_bd_mask] * 60,
#     edgecolors='black',
#     linewidths=0.8,
#     zorder=6,
#     label='BDs with Teff 6000–7200 K'
# )

# # Colorbar
# clb = plt.colorbar(scatter_bd, ax=ax, orientation='vertical')
# clb.set_label('Temperature (K)', labelpad=15, rotation=270, fontsize=20)

# # Labels and limits
# ax.set_xlabel('Mass [$M_{J}$]', fontsize=20)
# ax.set_ylabel('Density [g cm$^{-3}$]', fontsize=20)
# ax.set_xlim(0, 120)
# ax.set_ylim(-30, max(density + density_err_upper) * 1.1)

# # Vertical classification lines
# ax.axvline(x=12, color='red', linestyle='--')
# ax.axvline(x=84, color='red', linestyle='--')
# ax.text(45, max(density) + 45, 'Brown Dwarfs', fontsize=20, fontweight='bold')
# ax.text(86, max(density) + 45, 'Low-Mass Stars', fontsize=20, fontweight='bold')

# # Annotate TOI-2155b if present
# if 'TOI-2155b' in df['Name'].values:
#     idx = df[df['Name'] == 'TOI-2155b'].index[0]
#     ax.scatter(
#         mass[idx], density[idx],
#         s=radius[idx] * 60,
#         color='red',
#         marker='x',
#         linewidths=2,
#         zorder=6,
#         label='TOI-2155b'
#     )
#     ax.annotate(
#         '',
#         xy=(mass[idx], density[idx]),
#         xytext=(mass[idx] + 4, density[idx] + (max(density) * 0.15)),
#         arrowprops=dict(arrowstyle='->', color='blue', linewidth=2)
#     )
#     ax.text(mass[idx], density[idx] + (max(density) * 0.1), 'TOI-2155b', fontsize=20, fontweight='bold', ha='center')

# # Radius size legend
# desired_radii = np.array([0.5, 0.7, 0.9, 1.0, 1.5, 2.0, 3.0])
# desired_sizes = desired_radii * 70
# handles = [plt.scatter([], [], s=size, edgecolor='none', c='gray', alpha=0.9) for size in desired_sizes]
# labels = [f'{r} $R_J$' for r in desired_radii]

# # Add triangle legend entry for hot BDs
# avg_size = np.mean(radius[hot_bd_mask]) * 25 if np.any(hot_bd_mask) else 150
# triangle_handle = plt.Line2D([], [], marker='^', color='black',
#                              markerfacecolor='gray', markersize=np.sqrt(avg_size),
#                              linestyle='None')
# handles += [triangle_handle]
# labels += ['BDs with F-type host stars']

# ax.legend(handles, labels, loc="upper left", title="Radius [$R_J$]", fontsize=16, title_fontsize=20)

# plt.tight_layout()
# plt.savefig('mass_density_plot_with_hotbd_triangle.png', dpi=400, bbox_inches='tight')
# plt.show()



# #import os
# #import numpy as np
# #import pandas as pd
# #import matplotlib.pyplot as plt
# #import matplotlib.ticker as ticker
# #import seaborn as sns
# #from matplotlib import rcParams
# #plt.style.use('dark_background')
# ## ---------------------------
# ## Set Plotting Aesthetics
# ## ---------------------------
# #params = {
# #    'legend.fontsize': 12,
# #    'axes.labelsize': 12,
# #    'axes.titlesize': 12,
# #    'xtick.labelsize': 12,
# #    'ytick.labelsize': 12,
# #    'grid.color': 'k',
# #    'grid.linestyle': ':',
# #    'grid.linewidth': 0.5,
# #    'mathtext.fontset': 'stix',
# #    'mathtext.rm': 'DejaVu serif',
# #    'font.family': 'DejaVu serif',
# #    'font.serif': 'Times New Roman',
# #}
# #rcParams.update(params)
# #
# ## ---------------------------
# ## Load Data
# ## ---------------------------
# #data = pd.read_csv("browndwarf_final_earlier.csv")
# #
# ## ---------------------------
# ## Calculate Density Uncertainty
# ## ---------------------------
# #mass = data['mass_intermsof_jupitermass']
# #radius = data['radius_intermsof_jupiterradius']
# #mass_unc = data['massuncertainity_intermsof_jupitermass']
# #radius_unc = data['radiusuncertainity_intermsof_jupiterradius']
# #density = data['Density']
# #temp = data['TEMP']
# #
# ## Propagate error: fractional uncertainties added, radius cubed → multiply radius term by 3
# #density_unc = ((mass_unc / mass) + 3 * (radius_unc / radius)) * density
# #
# ## Append density uncertainty to DataFrame
# #data['density_uncertainity'] = density_unc
# #
# ## ---------------------------
# ## Prepare Arrays for Plotting
# ## ---------------------------
# #m = mass.to_numpy()
# #d = density.to_numpy()
# #r = radius.to_numpy()
# #t = temp.to_numpy()
# #
# ## ---------------------------
# ## Begin Plotting
# ## ---------------------------
# #fig, ax = plt.subplots(figsize=(10, 8))
# #
# ## Main scatter plot with temperature as color and radius as size
# #scatter = ax.scatter(
# #    m, d,
# #    edgecolors='none',
# #    s=r * 40,            # Scale point size with radius
# #    c=t,                 # Color by temperature
# #    cmap='jet',
# #    alpha=0.9
# #)
# #
# ## Colorbar setup
# #clb = plt.colorbar(scatter, ax=ax, orientation='vertical')
# #clb.set_label('Temperature (K)', labelpad=15, rotation=270)
# #
# ## Add error bars
# #ax.errorbar(
# #    x=m,
# #    y=d,
# #    yerr=data['density_uncertainity'],
# #    xerr=mass_unc,
# #    linestyle='None',
# #    color='gray',
# #    alpha=0.4,
# #    capsize=5
# #)
# #
# ## Axis labels and limits
# #ax.set_xlabel('Mass [$M_{J}$]')
# #ax.set_ylabel('Density [g cm$^{-3}$]')
# #ax.set_xlim(0, 120)
# #ax.set_ylim(-30, max(d) + 70)  # Add padding above data
# #
# ## ---------------------------
# ## Reference Lines and Annotations
# ## ---------------------------
# #ax.axvline(x=94, color='red', linestyle='--')  # Division between BDs and low-mass stars
# #ax.text(10, 125, 'Brown dwarfs', fontsize=10, fontweight='bold')
# #ax.text(94.5, 75, 'Low mass stars', fontsize=10, fontweight='bold')
# #
# ## Annotate TOI-2155b with arrow and label
# #ax.annotate(
# #    '',
# #    xy=(81, 111),       # Arrowhead
# #    xytext=(90, 150),   # Start of arrow
# #    arrowprops=dict(arrowstyle='->', color='red', linewidth=2)
# #)
# #ax.text(81, 127, 'TOI-2155b', fontsize=8, fontweight='bold', ha='center')
# #
# ## ---------------------------
# ## Custom Radius Size Legend
# ## ---------------------------
# #desired_radii = np.array([0.5, 0.7, 0.8, 0.9, 1, 2, 3])  # Radii in Jupiter units
# ## Convert to plot size scale (based on original scaling)
# #desired_sizes = (desired_radii / 11.21) * 490
# #
# ## Create dummy handles for legend
# #handles = [
# #    plt.scatter([], [], s=size, edgecolor='none', c='gray', alpha=0.6)
# #    for size in desired_sizes
# #]
# #labels = [f'{radius} $R_J$' for radius in desired_radii]
# #
# ## Draw legend on plot
# #legend = ax.legend(handles, labels, loc="upper right", title="Radius [$R_{J}$]", fontsize=12)
# #
# ## ---------------------------
# ## Save and Show
# ## ---------------------------
# #plt.savefig('mass_density_plot.png', format='png', dpi=400, bbox_inches='tight')
# #plt.show()
# #
