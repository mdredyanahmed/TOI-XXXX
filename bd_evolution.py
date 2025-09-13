import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# === Load model and data files ===
model_files = {
    '[M/H] = 0.0': 'nc+0.0_co1.0_age_bobcat_cleaned.csv',
    '[M/H] = -0.5': 'nc-0.5_co1.0_age_bobcat_cleaned.csv',
    '[M/H] = +0.5': 'nc+0.5_co1.0_age_bobcat_cleaned.csv',
}

baraffe_files = {
    'Age = 1 Gyr': 'baraffe_1gyr.csv',
    'Age = 2.27 Gyr': 'baraffe_2.27gyr_interpolated.csv',
    'Age = 3.2 Gyr': 'baraffe_3.2gyr_interpolated.csv',
    'Age = 5.1 Gyr': 'baraffe_5.1gyr_interpolated.csv',
    'Age = 10 Gyr': 'baraffe_10gyr.csv'
}

browndwarf_file = 'transiting_brown_dwarf_list_march2025.txt'

# === Read files ===
model_dfs = {label: pd.read_csv(fp).apply(pd.to_numeric, errors='coerce').dropna() 
             for label, fp in model_files.items()}
baraffe_dfs = {label: pd.read_csv(fp).apply(pd.to_numeric, errors='coerce').dropna() 
               for label, fp in baraffe_files.items()}
browndwarf_df = pd.read_csv(browndwarf_file, delim_whitespace=True)

# --- Extract brown dwarf data ---
browndwarf_mass = browndwarf_df['Mj']
browndwarf_radius = browndwarf_df['Rj']
browndwarf_mass_err = browndwarf_df['Mj_err1']
browndwarf_radius_err = browndwarf_df['Rj_err1']
bd_mask = browndwarf_mass <= 81.7
star_mask = browndwarf_mass >= 81.8

# --- Target parameters for TOI-2155 ---
target_age = 3.2      # Gyr
target_mh = 0.13

# --- Helper to access model dfs properly ---
def get_model_df(mh):
    if mh == 0.5:
        key = '[M/H] = +0.5'
    elif mh == 0.0:
        key = '[M/H] = 0.0'
    elif mh == -0.5:
        key = '[M/H] = -0.5'
    else:
        raise ValueError(f"No model available for [M/H]={mh}")
    return model_dfs[key]

# --- Bilinear interpolation function ---
def bilinear_interpolate(target_age, target_mh):
    mh_values = np.array([-0.5, 0.0, 0.5])
    # Find bounding metallicities
    mh_low = mh_values[mh_values <= target_mh].max()
    mh_high = mh_values[mh_values >= target_mh].min()
    
    df_low_mh = get_model_df(mh_low)
    df_high_mh = get_model_df(mh_high)

    def get_age_bracket(df, age_target):
        ages = df['age(Gyr)'].unique()
        age1 = ages[ages <= age_target].max()
        age2 = ages[ages >= age_target].min()
        return age1, age2

    age1_low, age2_low = get_age_bracket(df_low_mh, target_age)
    age1_high, age2_high = get_age_bracket(df_high_mh, target_age)

    # Interpolate radius along mass for each corner
    def radius_at_age(df, age):
        df_age = df[df['age(Gyr)'] == age]
        return interp1d(df_age['M/MSun'], df_age['R/Rsun'], bounds_error=False, fill_value=np.nan)

    f11 = radius_at_age(df_low_mh, age1_low)
    f12 = radius_at_age(df_low_mh, age2_low)
    f21 = radius_at_age(df_high_mh, age1_high)
    f22 = radius_at_age(df_high_mh, age2_high)

    # Common mass grid
    mass_min = max(df_low_mh['M/MSun'].min(), df_high_mh['M/MSun'].min())
    mass_max = min(df_low_mh['M/MSun'].max(), df_high_mh['M/MSun'].max())
    mass_common = np.linspace(mass_min, mass_max, 300)

    # Interpolate in age
    R_low = f11(mass_common) + (f12(mass_common) - f11(mass_common)) * (target_age - age1_low) / (age2_low - age1_low)
    R_high = f21(mass_common) + (f22(mass_common) - f21(mass_common)) * (target_age - age1_high) / (age2_high - age1_high)

    # Interpolate in metallicity
    R_target = R_low + (R_high - R_low) * (target_mh - mh_low) / (mh_high - mh_low)

    return mass_common, R_target

# --- Compute bilinear interpolation ---
mass_common_sonora, radius_interp = bilinear_interpolate(target_age, target_mh)
mass_common_jup = mass_common_sonora * 1047.56
radius_interp_jup = radius_interp * 9.731

# === Plotting ===
plt.style.use('seaborn-v0_8-white')
fig, ax = plt.subplots(figsize=(10, 8))

# --- Plot Baraffe 3.2 Gyr line with uncertainty band ---
df_32 = baraffe_dfs['Age = 3.2 Gyr']
df_227 = baraffe_dfs['Age = 2.27 Gyr']
df_51 = baraffe_dfs['Age = 5.1 Gyr']

mass_32 = df_32['Mass'] * 1047.56
radius_32 = df_32['Radius'] * 9.731
mass_227 = df_227['Mass'] * 1047.56
radius_227 = df_227['Radius'] * 9.731
mass_51 = df_51['Mass'] * 1047.56
radius_51 = df_51['Radius'] * 9.731

# Plot 3.2 Gyr line
ax.plot(mass_32, radius_32, linestyle='-', color='#2ca02c',
        label=r'Baraffe+2003, Age = $3.2^{+1.9}_{-0.9}$ Gyr', lw=3, zorder=1)

# Fill between 2.27 and 5.1 Gyr
common_mass = np.linspace(max(min(mass_227), min(mass_51)), min(max(mass_227), max(mass_51)), 300)
interp_227 = interp1d(mass_227, radius_227, bounds_error=False, fill_value=np.nan)
interp_51 = interp1d(mass_51, radius_51, bounds_error=False, fill_value=np.nan)
ax.fill_between(common_mass, interp_227(common_mass), interp_51(common_mass),
                color='#2ca02c', alpha=0.25, zorder=0)

# --- Plot Sonora bilinear interpolation ---
ax.plot(mass_common_jup, radius_interp_jup, color='purple', lw=2.5, linestyle='--',
        label=f'Sonora Bobcat+2021, [M/H]={target_mh}, Age={target_age} Gyr', zorder=4)

# --- Brown Dwarfs ---
ax.errorbar(browndwarf_mass[bd_mask], browndwarf_radius[bd_mask],
            xerr=browndwarf_mass_err[bd_mask], yerr=browndwarf_radius_err[bd_mask],
            fmt='o', color='#1f77b4', markersize=1.5, ecolor='#1f77b4', elinewidth=0.8, capsize=1,
            label='Brown Dwarfs', alpha=0.5, zorder=2)

# --- Low-Mass Stars ---
ax.errorbar(browndwarf_mass[star_mask], browndwarf_radius[star_mask],
            xerr=browndwarf_mass_err[star_mask], yerr=browndwarf_radius_err[star_mask],
            fmt='*', color='#ff7f0e', markersize=4, ecolor='#ff7f0e', elinewidth=0.8, capsize=2,
            label='Low-Mass Stars', alpha=0.9, zorder=2)

# --- Highlight TOI-2155b ---
planet_mass, planet_radius = 81.1, 0.975
mass_err = 1.1
radius_err_lower = 0.008
radius_err_upper = 0.008

ax.errorbar(
    planet_mass, planet_radius,
    xerr=mass_err, yerr=[[radius_err_lower], [radius_err_upper]],
    fmt='.', color='crimson', markersize=4, ecolor='crimson', elinewidth=1.5, capsize=2,
    label='TOI-2155b', zorder=10,
    markeredgecolor='black', markeredgewidth=1
)
ax.annotate('TOI-2155b', xy=(planet_mass, planet_radius),
            xytext=(planet_mass + 5, planet_radius + 0.3),
            arrowprops=dict(facecolor='crimson', arrowstyle='->', lw=2),
            fontsize=16, fontweight='bold', color='crimson', zorder=11)

# --- Vertical lines for BD regime ---
ax.axvline(x=13, color='red', linestyle='--', lw=1)
ax.axvline(x=81.7, color='red', linestyle='--', lw=1)
ax.text(40, 1.55, 'Brown Dwarfs', fontsize=16, fontweight='bold')
ax.text(86, 1.55, 'Low-Mass Stars', fontsize=15, fontweight='bold')

# --- Axis labels and limits ---
ax.set_xlabel('Mass [$M_\\mathrm{J}$]', fontsize=24)
ax.set_ylabel('Radius [$R_\\mathrm{J}$]', fontsize=24)
ax.set_xlim(8, 105)
ax.set_ylim(0.5, 1.8)
ax.tick_params(axis='both', which='major', labelsize=20, direction='out', length=6, width=1)
ax.tick_params(axis='both', which='minor', labelsize=20, direction='out', length=3, width=0.8)
ax.minorticks_on()

# --- Legend ---
legend = ax.legend(fontsize=13, frameon=True, loc='lower center', ncol=2, handlelength=3)
legend.get_frame().set_edgecolor('gray')
legend.get_frame().set_alpha(0.9)

plt.tight_layout()
#plt.savefig("bd_evolution_plot_with_sonora_bilinear_interpolation.png", dpi=600, bbox_inches='tight')
plt.show()


# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np
# from scipy.interpolate import interp1d
# import matplotlib.patches as patches

# # === Load model and data files ===
# model_files = {
#     '[M/H] = 0.0': 'nc+0.0_co1.0_age_bobcat_cleaned.csv',
#     '[M/H] = -0.5': 'nc-0.5_co1.0_age_bobcat_cleaned.csv',
#     '[M/H] = +0.5': 'nc+0.5_co1.0_age_bobcat_cleaned.csv',
# }
# baraffe_files = {
#     'Age = 1 Gyr': 'baraffe_1gyr.csv',
#     'Age = 2.27 Gyr': 'baraffe_2.27gyr_interpolated.csv',
#     'Age = 3.2 Gyr': 'baraffe_3.2gyr_interpolated.csv',
#     'Age = 5.1 Gyr': 'baraffe_5.1gyr_interpolated.csv',
#     'Age = 10 Gyr': 'baraffe_10gyr.csv'
# }
# browndwarf_file = 'transiting_brown_dwarf_list_march2025.txt'

# # === Read files ===
# model_dfs = {label: pd.read_csv(fp).apply(pd.to_numeric, errors='coerce').dropna() for label, fp in model_files.items()}
# baraffe_dfs = {label: pd.read_csv(fp).apply(pd.to_numeric, errors='coerce').dropna() for label, fp in baraffe_files.items()}
# browndwarf_df = pd.read_csv(browndwarf_file, delim_whitespace=True)

# # === Extract brown dwarf data ===
# browndwarf_mass = browndwarf_df['Mj']
# browndwarf_radius = browndwarf_df['Rj']
# browndwarf_mass_err = browndwarf_df['Mj_err1']
# browndwarf_radius_err = browndwarf_df['Rj_err1']
# bd_mask = browndwarf_mass <= 84
# star_mask = browndwarf_mass > 84

# # === Plot settings ===
# fig, ax = plt.subplots(figsize=(10, 8))

# baraffe_colors = {
#     'Age = 3.2 Gyr': '#2ca02c',  # Bright green
# }

# # === Plot only the 3.2 Gyr Baraffe line ===
# df_32 = baraffe_dfs['Age = 3.2 Gyr']
# mass_32 = df_32['Mass'] * 1047.56
# radius_32 = df_32['Radius'] * 9.731
# ax.plot(mass_32, radius_32,
#         linestyle='-', color=baraffe_colors['Age = 3.2 Gyr'],
#         label='Baraffe+2003, Age = 3.2 Gyr', lw=2)

# # === Fill between 2.27 and 5.1 Gyr models as uncertainty region ===
# df_227 = baraffe_dfs['Age = 2.27 Gyr']
# df_51 = baraffe_dfs['Age = 5.1 Gyr']
# mass_227 = df_227['Mass'] * 1047.56
# radius_227 = df_227['Radius'] * 9.731
# mass_51 = df_51['Mass'] * 1047.56
# radius_51 = df_51['Radius'] * 9.731

# # Interpolate to a common mass grid
# common_mass = np.linspace(
#     max(min(mass_227), min(mass_51)),
#     min(max(mass_227), max(mass_51)),
#     300
# )
# interp_227 = interp1d(mass_227, radius_227, bounds_error=False, fill_value=np.nan)
# interp_51 = interp1d(mass_51, radius_51, bounds_error=False, fill_value=np.nan)

# ax.fill_between(
#     common_mass, interp_227(common_mass), interp_51(common_mass),
#     color=baraffe_colors['Age = 3.2 Gyr'], alpha=0.15
# )

# # === Brown dwarfs and stars ===
# ax.errorbar(browndwarf_mass[bd_mask], browndwarf_radius[bd_mask],
#             xerr=browndwarf_mass_err[bd_mask], yerr=browndwarf_radius_err[bd_mask],
#             fmt='o', color='gray', markersize=3, ecolor='gray', elinewidth=0.5, capsize=2,
#             label='Transiting Brown Dwarfs', alpha=0.6)
# ax.errorbar(browndwarf_mass[star_mask], browndwarf_radius[star_mask],
#             xerr=browndwarf_mass_err[star_mask], yerr=browndwarf_radius_err[star_mask],
#             fmt='*', color='gray', markersize=7, ecolor='gray', elinewidth=0.5, capsize=2,
#             label='Low-Mass Stars', alpha=0.6)

# # === TOI-2155 b ===
# planet_mass, planet_radius = 80.5, 0.975
# mass_err, radius_err_lower, radius_err_upper = 1.1, [0.008], [0.008]

# # Error bars only
# ax.errorbar(
#     planet_mass, planet_radius,
#     xerr=mass_err, yerr=[radius_err_lower, radius_err_upper],
#     fmt='none', ecolor='red', elinewidth=1, capsize=2, zorder=2,label='TOI-2155b'
# )

# # Annotate TOI-2155b
# ax.annotate('TOI-2155b', xy=(planet_mass, planet_radius),
#             xytext=(planet_mass + 5, planet_radius + 0.2),
#             arrowprops=dict(facecolor='black', arrowstyle='->', lw=0.5),
#             fontsize=12, fontweight='bold',label='TOI-2155b')

# # === Add red rectangle around TOI-2155b ===
# rect = patches.Rectangle(
#     (78.5, 0.93),  # lower left corner (mass, radius)
#     width=4.5,     # mass range
#     height=0.1,    # radius range
#     linewidth=1.5,
#     edgecolor='black',
#     linestyle='--',
#     facecolor='none',
#     zorder=4
# )
# ax.add_patch(rect)

# # === Decorations ===
# ax.axvline(x=13, color='red', linestyle='--')
# ax.axvline(x=81.6, color='red', linestyle='--')
# ax.text(45, 1.450, 'Brown Dwarfs', fontsize=14, fontweight='bold')
# ax.text(85, 1.450, 'Low-Mass Stars', fontsize=14, fontweight='bold')
# ax.set_xlabel('Mass [$M_\\mathrm{J}$]', fontsize=25)
# ax.set_ylabel('Radius [$R_\\mathrm{J}$]', fontsize=25)
# ax.set_xlim(8, 105)
# ax.set_ylim(0.5, 1.8)
# ax.tick_params(axis='both', which='major', labelsize=22, direction='out', length=6, width=1)
# ax.tick_params(axis='both', which='minor', labelsize=18, direction='out', length=3, width=0.8)
# ax.minorticks_on()
# legend = ax.legend(fontsize=11, frameon=True, loc='best', ncol=2)
# legend.get_frame().set_edgecolor('gray')
# legend.get_frame().set_alpha(0.6)

# plt.tight_layout()
# plt.savefig("bd_evolution_plot.png", dpi=600)
# plt.show()


# # import pandas as pd
# # import matplotlib.pyplot as plt

# # # === Light background plot ===
# # #plt.style.use('dark_background')  # Clean white background, no gridlines

# # # File paths for Sonora models
# # model_files = {
# #     '[M/H] = 0.0': 'nc+0.0_co1.0_age_bobcat_cleaned.csv',
# #     '[M/H] = -0.5': 'nc-0.5_co1.0_age_bobcat_cleaned.csv',
# #     '[M/H] = +0.5': 'nc+0.5_co1.0_age_bobcat_cleaned.csv',
# # }

# # # File paths for Baraffe models
# # baraffe_files = {
# #     'Age = 1 Gyr': 'baraffe_1gyr.csv',

# #     'Age =  3 Gyr': 'baraffe_3gyr_interpolated.csv',
# #         'Age =  5 Gyr': 'baraffe_5gyr.csv'
# # }

# # # Read Sonora models
# # model_dfs = {}
# # for label, filepath in model_files.items():
# #     df = pd.read_csv(filepath)
# #     df = df.apply(pd.to_numeric, errors='coerce').dropna()
# #     model_dfs[label] = df

# # # Read Baraffe models
# # baraffe_dfs = {}
# # for label, filepath in baraffe_files.items():
# #     df = pd.read_csv(filepath)
# #     df = df.apply(pd.to_numeric, errors='coerce').dropna()
# #     baraffe_dfs[label] = df

# # # Read brown dwarf data
# # browndwarf_file = 'transiting_brown_dwarf_list_march2025.txt'
# # browndwarf_df = pd.read_csv(browndwarf_file, delim_whitespace=True)
# # browndwarf_mass = browndwarf_df['Mj']
# # browndwarf_radius = browndwarf_df['Rj']
# # browndwarf_mass_err = browndwarf_df['Mj_err1']
# # browndwarf_radius_err = browndwarf_df['Rj_err1']

# # # Separate brown dwarfs and low mass stars
# # bd_mask = browndwarf_mass <= 84
# # star_mask = browndwarf_mass > 84

# # # Target ages for Sonora models
# # target_ages = [3, 4]

# # # Linestyles and colors for Sonora
# # metallicity_styles = {
# #     '[M/H] = 0.0': '-',
# #     '[M/H] = -0.5': '--',
# #     '[M/H] = +0.5': ':'
# # }
# # custom_colors = {
# #     ('[M/H] = 0.0', 4): '#E69F00',
# #     ('[M/H] = 0.0', 3): '#56B4E9',
# #     ('[M/H] = -0.5', 3): '#D55E00',
# #     ('[M/H] = +0.5', 3): '#CC79A7',
# # }

# # # Linestyles and colors for Baraffe models
# # baraffe_styles = {
# #     'Age = 1 Gyr': ':',
  
# #     'Age =  3 Gyr': '-',
# #         'Age =  5 Gyr': '--'
# # }
# # baraffe_colors = {
# #     'Age = 1 Gyr': 'black',    # Use basic green
# #    # Use orange
# #     'Age =  3 Gyr': 'green' ,
# #         'Age =  5 Gyr': 'orange',     # Use red
# # }
# # # Plot
# # fig, ax = plt.subplots(figsize=(14, 12))

# # # Plot Sonora models
# # #for label, df in model_dfs.items():
# # #    for target_age in target_ages:
# # #        if label != '[M/H] = 0.0' and target_age != 3:
# # #            continue
# # #        nearest_age = df['age(Gyr)'].sub(target_age).abs().idxmin()
# # #        age_value = df.loc[nearest_age, 'age(Gyr)']
# # #        df_nearest = df[df['age(Gyr)'] == age_value]
# # #        linestyle = metallicity_styles[label]
# # #        color = custom_colors.get((label, target_age), 'black')
# # #        ax.plot(df_nearest['M/MSun'] * 1047.56, df_nearest['R/Rsun'] * 9.731,
# # #                linestyle=linestyle, color=color,
# # #                label=f'{label}, Age = {age_value:.2f} Gyr, Sonora Bobcat Model (2021) ', lw=2)

# # # Plot Baraffe models
# # for label, df in baraffe_dfs.items():
# #     ax.plot(df['Mass'] * 1047.56, df['Radius'] * 9.731,
# #             linestyle=baraffe_styles[label], color=baraffe_colors[label],
# #             label=f'{label} (Baraffe+2003)', lw=2)

# # # Plot brown dwarfs
# # ax.errorbar(browndwarf_mass[bd_mask], browndwarf_radius[bd_mask],
# #             xerr=browndwarf_mass_err[bd_mask], yerr=browndwarf_radius_err[bd_mask],
# #             fmt='o', color='blue', markersize=6,
# #             ecolor='blue', elinewidth=0.7, capsize=2,
# #             label='Transiting Brown Dwarfs')

# # # Plot low mass stars with different marker
# # ax.errorbar(browndwarf_mass[star_mask], browndwarf_radius[star_mask],
# #             xerr=browndwarf_mass_err[star_mask], yerr=browndwarf_radius_err[star_mask],
# #             fmt='*', color='#8B006A', markersize=12,
# #             ecolor='#8B006A', elinewidth=0.7, capsize=2,
# #             label='Low-Mass Stars')

# # # Plot TOI 2155 b
# # planet_mass = 80.60
# # planet_radius = 0.970
# # mass_err = 0.928
# # radius_err_lower = [0.020]
# # radius_err_upper = [0.018]
# # yerr = [radius_err_lower, radius_err_upper]

# # ax.errorbar(planet_mass, planet_radius,
# #             xerr=mass_err, yerr=yerr,
# #             fmt='o', color='#d62728', markersize=6,
# #             ecolor='#d62728', elinewidth=1.2, capsize=2,
# #             label='TOI 2155b')

# # # Region separators
# # ax.axvline(x=12, color='red', linestyle='--')
# # #ax.axvline(x=42, color='red', linestyle='--')
# # ax.axvline(x=84, color='red', linestyle='--')

# # #ax.text(13, 1.450, 'Low-mass Brown dwarfs', fontsize=13, fontweight='bold')
# # ax.text(45, 1.450, 'Brown dwarfs', fontsize=20, fontweight='bold')
# # ax.text(86, 1.450, 'Low Mass Stars', fontsize=20, fontweight='bold')

# # # Axis labels and formatting
# # ax.set_xlabel('Mass [$M_\\mathrm{J}$]', fontsize=20)
# # ax.set_ylabel('Radius [$R_\\mathrm{J}$]', fontsize=20)
# # ax.set_xlim(8, 105)
# # ax.set_ylim(0.5, 1.8)
# # ax.tick_params(axis='both', which='major', labelsize=11, direction='out', length=6, width=1)
# # ax.tick_params(axis='both', which='minor', direction='out', length=3, width=0.8)
# # ax.minorticks_on()

# # # Legend
# # legend = ax.legend(fontsize=11, frameon=True, loc='best')
# # legend.get_frame().set_edgecolor('gray')
# # legend.get_frame().set_alpha(0.6)

# # plt.tight_layout()
# # plt.savefig('TOI2155_evolution.png', format='png', dpi=400, bbox_inches='tight')
# # plt.show()






# # #import pandas as pd
# # #import matplotlib.pyplot as plt
# # #
# # #plt.style.use('dark_background')
# # ## === Light background plot ===
# # ##plt.style.use('seaborn-v0_8-white')  # Clean white background, no gridlines
# # #
# # ## File paths for Sonora models
# # #model_files = {
# # #    '[M/H] = 0.0': 'nc+0.0_co1.0_age_bobcat_cleaned.csv',
# # #    '[M/H] = -0.5': 'nc-0.5_co1.0_age_bobcat_cleaned.csv',
# # #    '[M/H] = +0.5': 'nc+0.5_co1.0_age_bobcat_cleaned.csv',
# # #}
# # #
# # ## File paths for Baraffe models
# # #baraffe_files = {
# # #    'Age = 1 Gyr': 'baraffe_1gyr.csv',
# # #    'Age =  5 Gyr': 'baraffe_5gyr.csv',
# # #    # 'Baraffe 10 Gyr': 'baraffe_10gyr.csv',
# # #    'Age =  3 Gyr':'baraffe_3gyr_interpolated.csv'
# # #}
# # #
# # ## Read Sonora models
# # #model_dfs = {}
# # #for label, filepath in model_files.items():
# # #    df = pd.read_csv(filepath)
# # #    df = df.apply(pd.to_numeric, errors='coerce').dropna()
# # #    model_dfs[label] = df
# # #
# # ## Read Baraffe models
# # #baraffe_dfs = {}
# # #for label, filepath in baraffe_files.items():
# # #    df = pd.read_csv(filepath)
# # #    df = df.apply(pd.to_numeric, errors='coerce').dropna()
# # #    baraffe_dfs[label] = df
# # #
# # ## Read brown dwarf data
# # #browndwarf_file = 'transiting_brown_dwarf_list_march2025.txt'
# # #browndwarf_df = pd.read_csv(browndwarf_file, delim_whitespace=True)
# # #browndwarf_mass = browndwarf_df['Mj']
# # #print(browndwarf_mass)
# # #browndwarf_radius = browndwarf_df['Rj']
# # #browndwarf_mass_err = browndwarf_df['Mj_err1']  # assuming this is the upper error
# # #browndwarf_radius_err = browndwarf_df['Rj_err1']  # assuming this is the upper error
# # #
# # ## Target ages for Sonora models
# # #target_ages = [3, 4, 13]
# # #
# # ## Linestyles and colors for Sonora
# # #metallicity_styles = {
# # #    '[M/H] = 0.0': '-',
# # #    '[M/H] = -0.5': '--',
# # #    '[M/H] = +0.5': ':'
# # #}
# # #custom_colors = {
# # #    ('[M/H] = 0.0', 4): '#E69F00',
# # #    ('[M/H] = 0.0', 3): '#56B4E9',
# # #    ('[M/H] = 0.0', 13): '#ff0000',
# # #    ('[M/H] = -0.5', 3): '#D55E00',
# # #    ('[M/H] = +0.5', 3): '#CC79A7',
# # #}
# # #
# # ## Linestyles and colors for Baraffe models
# # #baraffe_styles = {
# # #    'Age = 1 Gyr': ':',
# # #    'Age =  5 Gyr': '--',
# # #    'Age =  3 Gyr': '-'
# # #}
# # #baraffe_colors = {
# # #    'Age = 1 Gyr': '#9467bd',  # Purple
# # #    'Age =  5 Gyr': '#8c564b',  # Brown
# # # 'Age =  3 Gyr': '#009E73'  # Bright Red
# # #
# # #}
# # #
# # ## Plot
# # #fig, ax = plt.subplots(figsize=(8.5, 6.5))
# # #
# # ## Plot Sonora models
# # #for label, df in model_dfs.items():
# # #    for target_age in target_ages:
# # #        if label != '[M/H] = 0.0' and target_age != 3:
# # #            continue
# # #        nearest_age = df['age(Gyr)'].sub(target_age).abs().idxmin()
# # #        age_value = df.loc[nearest_age, 'age(Gyr)']
# # #        df_nearest = df[df['age(Gyr)'] == age_value]
# # #        linestyle = metallicity_styles[label]
# # #        color = custom_colors.get((label, target_age), 'black')
# # #        ax.plot(df_nearest['M/MSun'] * 1047.56, df_nearest['R/Rsun'] * 9.731,
# # #                linestyle=linestyle, color=color,
# # #                label=f'{label}, Age = {age_value:.2f} Gyr, Sonora Bobcat Model', lw=2)
# # #
# # ## Plot Baraffe models
# # #for label, df in baraffe_dfs.items():
# # #    ax.plot(df['Mass'] * 1047.56, df['Radius'] * 9.731,  # assuming Mass in Msun, Radius in Rsun
# # #            linestyle=baraffe_styles[label], color=baraffe_colors[label],
# # #            label=f'{label} (Baraffe 2003)', lw=2)
# # #
# # ## Plot transiting brown dwarfs
# # #ax.errorbar(browndwarf_mass, browndwarf_radius,
# # #            xerr=browndwarf_mass_err, yerr=browndwarf_radius_err,
# # #            fmt='o', color='gray', markersize=3,
# # #            ecolor='gray', elinewidth=0.7, capsize=2,
# # #            label='Transiting brown dwarfs')
# # #
# # ## Plot TOI 2155 b
# # #planet_mass = 80.60
# # #planet_radius = 0.970
# # #mass_err = 0.928
# # ## radius_err = 0.018
# # #radius_err_lower = [0.020]  # lower uncertainty
# # #radius_err_upper = [0.018]
# # #yerr = [radius_err_lower, radius_err_upper]
# # #
# # #ax.errorbar(planet_mass, planet_radius,
# # #            xerr=mass_err, yerr=yerr,
# # #            fmt='*', color='#d62728', markersize=12,
# # #            ecolor='#d62728', elinewidth=1.2, capsize=2,
# # #            label='TOI 2155b')
# # #
# # ## Axis labels and formatting
# # #ax.set_xlabel('Mass [$M_\\mathrm{J}$]', fontsize=14)
# # #ax.set_ylabel('Radius [$R_\\mathrm{J}$]', fontsize=14)
# # #ax.set_xlim(8, 90)
# # #ax.set_ylim(0.5, 1.650)
# # #ax.tick_params(axis='both', which='major', labelsize=11, direction='out', length=6, width=1)
# # #ax.tick_params(axis='both', which='minor', direction='out', length=3, width=0.8)
# # #ax.minorticks_on()
# # #
# # ## Legend
# # #legend = ax.legend(fontsize=7.5, frameon=True, loc='lower left')
# # #legend.get_frame().set_edgecolor('gray')
# # #legend.get_frame().set_alpha(0.6)
# # #
# # #plt.tight_layout()
# # #plt.savefig("bd_evolution_dark.png",dpi=400)
# # #plt.show()
# # #
