import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Constants for all siRNA versions:
k_prod_values = [0.00756, 0.00756, 0.00189, 0.002835, 0.00756, 0.004725, 0.00756, 0.00756, 0.004725, 0.00756]
k_deg_values = [0.084, 0.084, 0.021, 0.0315, 0.084, 0.0525, 0.084, 0.084, 0.0525, 0.084]
k_fold_values = [64.65405566, 37.69417324, 149.5186658, 78.90681233, 42.69218927, 47.16369143, 46.77382026, 78.25454168, 19.24303282, 37.38258033]

# Yeast growth constants:
T_d = 90
k_growth = np.log(2) / T_d
N_max = 100
N0 = 1

# Dicer cleavage constants:
k_cat = 7e-5  # s^-1
Dicer_molecules = 749032
cell_volume = 3.5e-14  # L
Avogadro_num = 6.022e23
Dicer_conc = Dicer_molecules / (Avogadro_num * cell_volume)

# Convert Dicer concentration to ng/ml
MW_DICER1 = 219e3  # g/mol
Dicer_conc = Dicer_conc * MW_DICER1 * 1e6  # ng/ml

# Initial shRNA conditions:
R0 = 0
F0 = 0
C0 = 0  # Initial condition for cleaved shRNA

# Time array:
T_max = 1440
times = np.linspace(0, T_max, 500)

# ODE for yeast growth:
def model_yeast(N, t):
    return k_growth * N * (1 - N / N_max)

def N_to_OD_yeast(N):
    return 10 * (N/N_max)

# Adjusted ODE system (in the order [R, F, N, C]):
def model(y, t, k_prod, k_deg, k_fold):
    R, F, N, C = y
    p = 2
    dR = k_prod * N - k_deg * R - k_fold * R
    dF = k_fold * R - k_cat * Dicer_conc * F
    dC = k_cat * Dicer_conc * F
    dN = k_growth * N * (1 - (N / N_max)**p)
    return [dR, dF, dN, dC]

# Adjusting font sizes for improved readability
plt.rcParams.update({'font.size': 18, 'axes.labelsize': 18, 'axes.titlesize': 18, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 12})
fig, ax4 = plt.subplots(1, 1, figsize=(12, 10))

fig.patch.set_facecolor('none')
ax4.patch.set_facecolor('none')

color_value = '#ffe8c0'

# Plotting yeast growth
N_values_yeast = odeint(model_yeast, N0, times)
ax4.plot(times, N_to_OD_yeast(N_values_yeast), color='tab:red', label='Yeast Population (OD)', linestyle='--')
ax4.set_xlabel('Time, min', color=color_value)
ax4.set_ylabel('Yeast culture OD', color=color_value)
ax4.tick_params(axis='y', labelcolor=color_value)

ax5 = ax4.twinx()
color = 'tab:green'
ax5.set_ylabel('Folded shRNA, ng/ml', color=color_value)
ax5.tick_params(axis='y', labelcolor=color_value)

# Plotting cleaved shRNA
ax6 = ax4.twinx()
color = 'tab:blue'
ax6.spines['right'].set_position(('outward', 100))
ax6.set_ylabel('siRNA, ng/ml', color=color_value, labelpad=20)
ax6.tick_params(axis='y', labelcolor=color_value)

groups = {
    "siRNA-3": [2],
    "siRNA-4": [3],
    "siRNA-6, siRNA-9": [5, 8],
    "siRNA-1, siRNA-2, siRNA-5, siRNA-7, siRNA-8, siRNA-10": [0, 1, 4, 6, 7, 9]
}
group_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

folded_labels_plotted = []
cleaved_labels_plotted = []

for (group_name, indices), group_color in zip(groups.items(), group_colors):
    for idx in indices:
        result = odeint(model, [R0, F0, N0, C0], times, args=(k_prod_values[idx], k_deg_values[idx], k_fold_values[idx]))
        R_values, F_values, N_values, C_values = result.T
        if group_name not in folded_labels_plotted:
            ax5.plot(times, F_values, color=group_color, label=f"{group_name} Folded")
            folded_labels_plotted.append(group_name)
        else:
            ax5.plot(times, F_values, color=group_color)
        
        if group_name not in cleaved_labels_plotted:
            ax6.plot(times, C_values, color=group_color, linestyle='-.', label=f"{group_name} Cleaved")
            cleaved_labels_plotted.append(group_name)
        else:
            ax6.plot(times, C_values, color=group_color, linestyle='-.')

ax4.legend(loc='upper left')
ax5.legend(loc='lower left', bbox_to_anchor=(0, -0.4))
ax6.legend(loc='lower right', bbox_to_anchor=(1, -0.3))

def set_axis_colors(ax):
    ax.set_xlabel(ax.get_xlabel(), color=color_value)
    ax.set_ylabel(ax.get_ylabel(), color=color_value)
    ax.tick_params(axis='y', labelcolor=color_value, colors=color_value)
    ax.tick_params(axis='x', colors=color_value)
    for spine in ax.spines.values():
        spine.set_edgecolor(color_value)

set_axis_colors(ax4)
set_axis_colors(ax5)
set_axis_colors(ax6)

plt.tight_layout()
plt.savefig('step2.png', transparent=True)
plt.show()

cleaved_values_at_24h = []

for idx in range(10):  # Since there are 10 siRNA versions
    result = odeint(model, [R0, F0, N0, C0], times, args=(k_prod_values[idx], k_deg_values[idx], k_fold_values[idx]))
    R_values, F_values, N_values, C_values = result.T
    idx_24h = np.argmin(np.abs(times - 1440))
    cleaved_values_at_24h.append(C_values[idx_24h])

siRNA_labels = ["siRNA-1", "siRNA-2", "siRNA-3", "siRNA-4", "siRNA-5", "siRNA-6", "siRNA-7", "siRNA-8", "siRNA-9", "siRNA-10"]

print("Cleaved shRNA values at 24h (ng/ml) for each siRNA version:")
for i, value in enumerate(cleaved_values_at_24h):
    print(f"{siRNA_labels[i]}: {value:.4f}")
