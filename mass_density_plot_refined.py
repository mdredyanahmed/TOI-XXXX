import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from matplotlib import rcParams

# ---------------------------
# Set Plotting Aesthetics
# ---------------------------
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
data = pd.read_csv("browndwarf_final_earlier.csv")

# ---------------------------
# Calculate Density Uncertainty
# ---------------------------
mass = data['mass_intermsof_jupitermass']
radius = data['radius_intermsof_jupiterradius']
mass_unc = data['massuncertainity_intermsof_jupitermass']
radius_unc = data['radiusuncertainity_intermsof_jupiterradius']
density = data['Density']
temp = data['TEMP']

# Propagate error: fractional uncertainties added, radius cubed → multiply radius term by 3
density_unc = ((mass_unc / mass) + 3 * (radius_unc / radius)) * density

# Append density uncertainty to DataFrame
data['density_uncertainity'] = density_unc

# ---------------------------
# Prepare Arrays for Plotting
# ---------------------------
m = mass.to_numpy()
d = density.to_numpy()
r = radius.to_numpy()
t = temp.to_numpy()

# ---------------------------
# Begin Plotting
# ---------------------------
fig, ax = plt.subplots(figsize=(10, 8))

# Main scatter plot with temperature as color and radius as size
scatter = ax.scatter(
    m, d,
    edgecolors='none',
    s=r * 40,            # Scale point size with radius
    c=t,                 # Color by temperature
    cmap='jet',
    alpha=0.9
)

# Colorbar setup
clb = plt.colorbar(scatter, ax=ax, orientation='vertical')
clb.set_label('Temperature (K)', labelpad=15, rotation=270)

# Add error bars
ax.errorbar(
    x=m,
    y=d,
    yerr=data['density_uncertainity'],
    xerr=mass_unc,
    linestyle='None',
    color='k',
    alpha=0.4,
    capsize=5
)

# Axis labels and limits
ax.set_xlabel('Mass [$M_{J}$]')
ax.set_ylabel('Density [g cm$^{-3}$]')
ax.set_xlim(0, 120)
ax.set_ylim(-30, max(d) + 70)  # Add padding above data

# ---------------------------
# Reference Lines and Annotations
# ---------------------------
ax.axvline(x=94, color='red', linestyle='--')  # Division between BDs and low-mass stars
ax.text(10, 125, 'Brown dwarfs', fontsize=10, fontweight='bold')
ax.text(94.5, 75, 'Low mass stars', fontsize=10, fontweight='bold')

# Annotate TOI-2155b with arrow and label
ax.annotate(
    '',
    xy=(81, 121),       # Arrowhead
    xytext=(90, 160),   # Start of arrow
    arrowprops=dict(arrowstyle='->', color='red', linewidth=2)
)
ax.text(81, 127, 'TOI-2155b', fontsize=8, fontweight='bold', ha='center')

# ---------------------------
# Custom Radius Size Legend
# ---------------------------
desired_radii = np.array([0.5, 0.7, 0.8, 0.9, 1, 2, 3])  # Radii in Jupiter units
# Convert to plot size scale (based on original scaling)
desired_sizes = (desired_radii / 11.21) * 490

# Create dummy handles for legend
handles = [
    plt.scatter([], [], s=size, edgecolor='none', c='k', alpha=0.6)
    for size in desired_sizes
]
labels = [f'{radius} $R_J$' for radius in desired_radii]

# Draw legend on plot
legend = ax.legend(handles, labels, loc="upper right", title="Radius [$R_{J}$]", fontsize=12)

# ---------------------------
# Save and Show
# ---------------------------
plt.savefig('mass_density_plot.png', format='png', dpi=400, bbox_inches='tight')
plt.show()

