import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import corner
from math import log10, floor

# === Load data ===
filename = 'without_gp_posterior_samples_RM_with_model_TESS.csv'  # Replace with your file path
data = pd.read_csv(filename)

param_names = [
    'b_rr',
    'b_rsuma',
    'b_cosi',
    'b_epoch',
    'b_period',
    'b_K',
    'b_f_c',
    'b_f_s',
    'host_ldc_q1_TESS',
    'host_ldc_q2_TESS',
    'ln_err_flux_TESS',
    'ln_jitter_rv_Raphael',
     'baseline_offset_flux_TESS',
    'baseline_offset_rv_Raphael',
    # 'baseline_gp_offset_flux_TESS',
    # 'baseline_gp_matern32_lnsigma_flux_TESS',
    # 'baseline_gp_matern32_lnrho_flux_TESS'
]

samples = data[param_names].copy()

# === LaTeX labels dictionary ===
latex_labels = {
    'b_rr': r"$R_b / R_\star$",
    'b_rsuma': r"$(R_\star + R_b) / a_b$",
    'b_cosi': r"$\cos i_b$",
    'b_epoch': r"$T_{0;b}$",
    'b_period': r"$P_b$",
    'b_K': r"$K_b$",
    'b_f_c': r"$\sqrt{e}\cos \omega_b$",
    'b_f_s': r"$\sqrt{e}\sin \omega_b$",
    'host_ldc_q1_TESS': r"$q_{1;\mathrm{TESS}}$",
    'host_ldc_q2_TESS': r"$q_{2;\mathrm{TESS}}$",
    'ln_err_flux_TESS': r"$\log \sigma_\mathrm{TESS}$",
    'ln_jitter_rv_Raphael': r"$\ln \sigma_\mathrm{jitter} (RV_\mathrm{TRES})$",
    'baseline_offset_rv_Raphael': r"$\Delta_\mathrm{TRES}$",
    'baseline_offset_flux_TESS': r"$\Delta_\mathrm{TESS}$",
    'baseline_offset_flux_Raphael': r"$\Delta_\mathrm{TRES}$",
    
    'baseline_gp_matern32_lnsigma_flux_TESS': r"$\mathrm{gp \ln \sigma (TESS)}$",
    'baseline_gp_matern32_lnrho_flux_TESS': r"$\mathrm{gp \ln \rho (TESS)}$"
}

# === Units dictionary ===
units = {
    'b_rr': '',
    'b_rsuma': '',
    'b_cosi': '',
    'b_epoch': 'BJD',
    'b_period': 'days',
    'b_K': 'km/s',
    'b_f_c': '',
    'b_f_s': '',
    'host_ldc_q1_TESS': '',
    'host_ldc_q2_TESS': '',
    'ln_err_flux_TESS': r'$\log \mathrm{rel.\ flux.}$',
    'ln_jitter_rv_Raphael': 'km/s',
     'baseline_offset_flux_TESS': '',
    'baseline_offset_rv_Raphael': '',
    # 'baseline_gp_offset_flux_TESS': '',
    # 'baseline_gp_matern32_lnsigma_flux_TESS': '',
    # 'baseline_gp_matern32_lnrho_flux_TESS': ''
}

# === Function to format median and uncertainties with 2 significant figures ===
def format_sigfigs(p16, p50, p84, sigfigs=2):
    err_plus = p84 - p50
    err_minus = p50 - p16
    err = max(err_plus, err_minus)

    if err == 0:
        return f"{p50:.3f} +0.00 -0.00"

    order = floor(log10(err))
    digits = sigfigs - 1 - order
    fmt = f".{max(digits, 0)}f"

    return f"{format(p50, fmt)}^{{+{format(err_plus, fmt)}}}_{{-{format(err_minus, fmt)}}}"

# === Prepare labels for corner axes ===
labels = [latex_labels.get(p, p) for p in param_names]

# === Make corner plot without default titles ===
fig = corner.corner(
    samples,
    labels=labels,
    show_titles=False,
    label_kwargs={"fontsize": 16},
    quantiles=[0.16, 0.5, 0.84],
    smooth=0.9,
    color="C0"
)

# === Add custom LaTeX formatted titles with uncertainties above histograms ===
axes = np.array(fig.axes).reshape(len(param_names), len(param_names))
for i, col in enumerate(param_names):
    arr = samples[col].values
    p16, p50, p84 = np.percentile(arr, [16, 50, 84])
    formatted_val = format_sigfigs(p16, p50, p84, sigfigs=2)
    axes[i, i].set_title(f"{latex_labels[col]} = ${formatted_val}$", fontsize=16)

plt.tight_layout()

# === Save the corner plot figure ===
plot_filename = "corner_plot_with_uncertainties.pdf"
fig.savefig(plot_filename)
print(f"Corner plot saved as '{plot_filename}'")

plt.show()

# === Prepare and save LaTeX table with percentiles ===
lines = []
lines.append(r"\begin{tabular}{l l l l}")
lines.append(r"\hline")
lines.append(r"\multicolumn{4}{c}{\textit{Fitted parameters}} \\")
lines.append(r"\hline")

for param in param_names:
    arr = samples[param].values
    p16, p50, p84 = np.percentile(arr, [16, 50, 84])
    val_str = format_sigfigs(p16, p50, p84, sigfigs=2)
    unit_str = units.get(param, '')
    label = latex_labels.get(param, param)
    lines.append(f"{label} & ${val_str}$ & {unit_str} & fit \\\\")

lines.append(r"\hline")
lines.append(r"\end{tabular}")

latex_filename = "fitted_parameters_table.tex"
with open(latex_filename, "w") as f:
    f.write("\n".join(lines))

print(f"LaTeX table saved as '{latex_filename}'")

