import pandas as pd
import matplotlib.pyplot as plt

# Load raw stitched normalized data
raw_df = pd.read_csv('TOI2155_spoc_stitched_normalized.csv')

# Load detrended flux data
detrended_df = pd.read_csv('TESS.csv')
detrended_df['#time'] = detrended_df['#time'] - 2457000  # Adjust time as you did

# Define the time ranges (sectors) you want to plot
time_ranges = [
    (2715, 2770),
    (2910, 2963),
    (3285, 3315),
    (3430, 3480),
    (3635, 3665)
]

# Colors for plotting
raw_color = '#0077b6'        # strong cyan-blue
detrended_color = '#ff6f3c'  # muted orange

# Create subplots, one for each time range
fig, axes = plt.subplots(len(time_ranges), 1, figsize=(12, 4 * len(time_ranges)), sharex=False)

for ax, (t_start, t_end) in zip(axes, time_ranges):
    # Filter raw data for current time range
    raw_mask = (raw_df['#time'] >= t_start) & (raw_df['#time'] <= t_end)
    raw_slice = raw_df[raw_mask]

    # Filter detrended data for current time range
    detrended_mask = (detrended_df['#time'] >= t_start) & (detrended_df['#time'] <= t_end)
    detrended_slice = detrended_df[detrended_mask]

    # Plot raw data (cyan-blue dots)
    ax.plot(raw_slice['#time'], raw_slice['pdcsap_flux'], 'o', color=raw_color, markersize=1, alpha=1, label='Raw Data (TESS-SPOC)')

    # Plot detrended data (muted orange dots)
    ax.plot(detrended_slice['#time'], detrended_slice['flux'], 'o', color=detrended_color, markersize=0.5, alpha=0.85, label='Detrended Data (TESS-SPOC)')

    ax.set_xlim(t_start, t_end)
    
    
    # ax.set_title(f'Time range: {t_start} - {t_end}', fontsize=14)
    # ax.grid(True, linestyle='--', alpha=0.3)

axes[-1].set_xlabel('Time (BJD - 2457000)', fontsize=12)
axes[2].set_ylabel('Relative Flux', fontsize=12)
axes[0].legend(loc='lower center', fontsize=12)
plt.tight_layout()
plt.savefig("allsec_detrended_plot.png", dpi=400)
plt.show()

#dark background
#import pandas as pd
#import matplotlib.pyplot as plt
#
## Load raw stitched normalized data
#raw_df = pd.read_csv('TOI2155_spoc_stitched_normalized.csv')
#
## Load detrended flux data
#detrended_df = pd.read_csv('TESS.csv')
#detrended_df['#time'] = detrended_df['#time'] - 2457000  # Adjust time as you did
#
## Define the time ranges (sectors) you want to plot
#time_ranges = [
#    (2715, 2770),
#    (2910, 2963),
#    (3285, 3315),
#    (3430, 3480),
#    (3635, 3665)
#]
#
## Colors for plotting (bright colors for dark bg)
#raw_color = '#00d6ff'        # bright cyan
#detrended_color = '#ffa94d'  # bright orange
#
## Dark background style
#plt.style.use('dark_background')
#
## Create subplots, one for each time range
#fig, axes = plt.subplots(len(time_ranges), 1, figsize=(12, 4 * len(time_ranges)), sharex=False)
#fig.patch.set_facecolor('#121212')  # very dark background
#
#for ax, (t_start, t_end) in zip(axes, time_ranges):
#    # Set axis background color
#    ax.set_facecolor('#121212')
#    
#    # Filter raw data for current time range
#    raw_mask = (raw_df['#time'] >= t_start) & (raw_df['#time'] <= t_end)
#    raw_slice = raw_df[raw_mask]
#
#    # Filter detrended data for current time range
#    detrended_mask = (detrended_df['#time'] >= t_start) & (detrended_df['#time'] <= t_end)
#    detrended_slice = detrended_df[detrended_mask]
#
#    # Plot raw data (bright cyan dots)
#    ax.plot(raw_slice['#time'], raw_slice['pdcsap_flux'], 'o', color=raw_color, markersize=2, alpha=1, label='Raw Data (TESS-SPOC)')
#
#    # Plot detrended data (bright orange dots)
#    ax.plot(detrended_slice['#time'], detrended_slice['flux'], 'o', color=detrended_color, markersize=1, alpha=0.9, label='Detrended Data (TESS-SPOC)')
#
#    ax.set_xlim(t_start, t_end)
#    
#    # Grid lines with subtle color
#    # ax.grid(True, linestyle='--', linewidth=0.5, color='#444444', alpha=0.7)
#    #
#    # Set axis spine and tick colors to light gray
#    for spine in ax.spines.values():
#        spine.set_color('#cccccc')
#    ax.tick_params(colors='#cccccc', labelsize=11)
#
## Labels with light gray color
#axes[-1].set_xlabel('Time (BJD - 2457000)', fontsize=13, color='#cccccc')
#axes[2].set_ylabel('Relative Flux', fontsize=13, color='#cccccc')
#
## Legend on top subplot with custom colors and frame
#axes[0].legend(loc='lower center', fontsize=12, frameon=True, facecolor='#222222', edgecolor='#444444', labelcolor='#cccccc')
#
#plt.tight_layout()
#plt.savefig("allsec_detrended_dark_bg.png", dpi=400)
#plt.show()
#
