#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase-folded photometry + residuals per instrument, styled data points.
"""

import matplotlib.pyplot as plt
from allesfitter import allesclass

# Load your fit results
alles = allesclass('allesfit_eccentric_all_instrument_ground_based')

instruments = ['WST', 'WBRO', 'LCOGT']
companion = 'b'
colors = {
    'WST': 'tab:cyan',
    'WBRO': 'tab:red',
    'LCOGT': 'tab:green',
}

for inst in instruments:
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(6, 6),
        gridspec_kw={'height_ratios': [3, 1]},
        sharex=True
    )
    fig.subplots_adjust(hspace=0)

    # Plot data + model (top) and residuals (bottom)
    #data_label = f'Data ({inst})'
    alles.plot(inst, companion, 'phase', ax=ax_top)
    alles.plot(inst, companion, 'phase_residuals', ax=ax_bot)
    # After plotting
    ax_top.set_title(f'{inst}', fontsize=14, loc='center')  # Dynamic title per instrument

#    ax_top.title.set_visible(False)
    # Style and label the data line in the top panel
    if len(ax_top.lines) >= 2:
        data_line = ax_top.lines[1]
 #	 data_line.set_label(f'Data ({inst})')
        data_line.set_color(colors[inst])
        data_line.set_marker('.')
        data_line.set_linestyle('none')
        data_line.set_markersize(6)
#        ax_top.legend(fontsize='small')

   # if len(ax_top.lines) >= 1:
    #    #model_line = ax_top.lines[0]
     #   model_line.set_label('Model (Allesfitter)')
      #  model_line.set_color('red')
       # model_line.set_linestyle('-')
       # model_line.set_linewidth(1.5)
       # ax_top.legend(fontsize='small')
#yle residuals panel
    ax_bot.axhline(0, color='gray', linestyle='--', linewidth=1)
    for line in ax_bot.lines:
        line.set_color(colors[inst])
        line.set_marker('.')
        line.set_linestyle('none')
        line.set_markersize(6)

    # Labels and formatting
    ax_bot.set_xlabel('Phase')
    ax_top.set_ylabel('Relative flux-Baseline')
    ax_bot.set_ylabel('Residuals')
    ax_top.set_xlim(-0.05, 0.05)
    ax_bot.title.set_visible(False)
#ax_top.set_title('')
#    ax_top.grid(False)
 #   ax_bot.grid(False)
    plt.tight_layout()
    # Save figure
    fig.savefig(f'{inst}_phase_residuals_styled.png', bbox_inches='tight',dpi=600)
