import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({
    'font.size': 14,
    'font.family': 'serif',
    'axes.labelsize': 16,
    'axes.titlesize': 18,
    'xtick.labelsize': 25,
    'ytick.labelsize': 25,
    'legend.fontsize': 14,
    'axes.linewidth': 1.2,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'xtick.top': True,
    'ytick.right': True,
    'grid.alpha': 0.3,
    'grid.linestyle': '--'
})
# === Load & Clean Data ===
df = pd.read_csv('browndwarf_list_with_density_asymmetric.csv')
df.columns = df.columns.str.strip()

# === Extract Variables ===
mass = df['Mj']
mass_err_low = df['Mj_err1']
mass_err_high = df['Mj_err2']
period = df['P']
ecc = df['e']

# === Masks ===
bd_mask = mass <= 84       # Brown dwarfs
star_mask = mass > 84      # Low-mass stars

# === Plot Setup ===
fig, ax = plt.subplots(figsize=(9, 6))
ax.errorbar(mass, period, 
            xerr=[mass_err_low, mass_err_high], 
            fmt='o', 
            ecolor='gray', elinewidth=1, capsize=3,
            markersize=6, alpha=0.8)

# === Desert Boundaries ===
ax.axvline(13, color='red', linestyle='--', linewidth=1)
ax.axvline(84, color='red', linestyle='--', linewidth=1)
ax.axhline(100, color='black', linestyle='dotted', linewidth=1)

# === Plot Brown Dwarfs ===
sc_bd = ax.scatter(mass[bd_mask], period[bd_mask],
                   c=ecc[bd_mask], cmap='plasma', s=60,
                   edgecolor='k', zorder=3, label='Brown Dwarfs')

# === Plot Low-Mass Stars (as Stars) ===
sc_star = ax.scatter(mass[star_mask], period[star_mask],
                     c=ecc[star_mask], cmap='plasma', s=100,
                     edgecolor='k', marker='*', zorder=4,
                     label='Low-Mass Stars')

# === TOI-2155 Marker ===
toi2155_mass = 80.6
toi2155_period = 3.7246960

ax.scatter(toi2155_mass, toi2155_period,
           color='red', edgecolor='k', s=200, marker='X', zorder=5, label='TOI-2155b')

# ax.annotate('TOI-2155',
            # xy=(toi2155_mass, toi2155_period),
            # xytext=(1.2 * toi2155_mass, 1.1 * toi2155_period),
            # arrowprops=dict(arrowstyle='->', color='red'),
            # fontsize=10, color='red')
# 
# === Colorbar ===
cb = plt.colorbar(sc_bd, ax=ax)
cb.set_label('Eccentricity',fontsize=25)

# === Axes and Labels ===
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Mass [$M_\\mathrm{J}$]',fontsize=25)
ax.set_ylabel('Period [days]',fontsize=25)
# ax.set_title('The Brown Dwarf Orbital Desert', fontsize=14)

plt.ylim(1, 200)
plt.xlim(10, 120)

# === Legend & Layout ===
ax.legend(loc='best')
plt.tight_layout()
plt.savefig('brown_dwarf_ecc.png', dpi=600, bbox_inches='tight')

print(f"{len(df[bd_mask])} brown dwarfs and {len(df[star_mask])} low-mass stars in the plot.")
plt.show()

