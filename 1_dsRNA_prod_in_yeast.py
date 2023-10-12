import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.lines import Line2D

# Constants for all siRNA versions:
k_prod_values = [0.00756, 0.00756, 0.00189, 0.002835, 0.00756, 0.004725, 0.00756, 0.00756, 0.004725, 0.00756]
k_deg_values = [0.084, 0.084, 0.021, 0.0315, 0.084, 0.0525, 0.084, 0.084, 0.0525, 0.084]
k_fold_values = [64.65405566, 37.69417324, 149.5186658, 78.90681233, 42.69218927, 47.16369143, 46.77382026, 78.25454168, 19.24303282, 37.38258033]

# Yeast growth constants:
T_d = 90  # doubling time in minutes
k_growth = np.log(2) / 90
# k_growth = 0.0077
N_max = 100
N_intermediate = 60
N0 = 1

# Initial shRNA conditions:
R0 = 0
F0 = 0

# Time array:
T_max = 1440
times = np.linspace(0, T_max, 500)

# ODE for yeast growth:
def model_yeast(N, t):
    dN = k_growth * N * (1 - N / N_max)
    return dN

def N_to_OD_yeast(N):
    return 10 * (N/N_max)

# ODE system (in the order [R, F, N]):
def model(y, t, k_prod, k_deg, k_fold):
    R, F, N = y

    p = 2

    dR = k_prod * N - k_deg * R - k_fold * R
    dF = k_fold * R
    # dN = k_growth * N * (1 - N / N_max)
    #dN = k_growth * N * np.log(N_max/N)
    dN = k_growth * N * (1 - (N / N_max)**p)
    
    return [dR, dF, dN]

def N_to_OD(N):
    return 10 * (np.log10(N/N0) / np.log10(N_max/N0))

# Adjusting font sizes for improved readability
plt.rcParams.update({'font.size': 18, 'axes.labelsize': 18, 'axes.titlesize': 18, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 12})

# Combined figure with subplots
fig, (ax1, ax4) = plt.subplots(2, 1, figsize=(12, 12))

fig.patch.set_facecolor('none')
ax1.patch.set_facecolor('none')
ax4.patch.set_facecolor('none')

# Solve the ODE for the yeast growth
N_values_yeast = odeint(model_yeast, N0, times)
ax1.plot(times, N_to_OD_yeast(N_values_yeast), color='tab:red', label='Yeast Population (OD)', linestyle='--')

color = 'tab:red'
ax1.set_xlabel('Time, min')
ax1.set_ylabel('Yeast Population log (OD)', color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_yscale('log')

ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel('Unfolded shRNA, ng/ml', color=color)
ax2.tick_params(axis='y', labelcolor=color)

for i in range(10):
    result = odeint(model, [R0, F0, N0], times, args=(k_prod_values[i], k_deg_values[i], k_fold_values[i]))
    R_values = result[:, 0]
    ax2.plot(times, R_values, label=f'shRNA-{i+1}')

ax1.legend(loc='upper left')
ax2.legend(loc='upper right')

ax1.set_yscale('log', basey=10)
ax1.yaxis.set_major_formatter(plt.ScalarFormatter())
ax1.get_yaxis().set_minor_formatter(plt.NullFormatter())

color_value = '#ffe8c0'

# Set color for ax1
ax1.set_xlabel('Time, min', color=color_value)
ax1.set_ylabel('Yeast culture OD', color=color_value)
ax1.tick_params(axis='y', labelcolor=color_value, colors=color_value)
ax1.tick_params(axis='x', colors=color_value)

# Set color for ax2
ax2.set_ylabel('Unfolded shRNA, ng/ml', color=color_value)
ax2.tick_params(axis='y', labelcolor=color_value, colors=color_value)

# Solve the ODE for the yeast growth
N_values_yeast = odeint(model_yeast, N0, times)
ax4.plot(times, N_to_OD_yeast(N_values_yeast), color='tab:red', label='Yeast Population (OD)', linestyle='--')

color = 'tab:red'
ax4.set_xlabel('Time, min')
ax4.set_ylabel('Yeast culture OD', color=color)
ax4.tick_params(axis='y', labelcolor=color)
ax4.set_yscale('log')

ax5 = ax4.twinx()
color = 'tab:green'
ax5.set_ylabel('Folded shRNA, ng/ml', color=color)
ax5.tick_params(axis='y', labelcolor=color)

# Define the groups
groups = {
    "shRNA-3": [2],  # siRNA-3
    "shRNA-4": [3],    # siRNA-4
    "shRNA-6, shRNA-9": [5, 8], # siRNA-6, siRNA-9
    "shRNA-1, shRNA-2, shRNA-5, shRNA-7, shRNA-8, shRNA-10": [0, 1, 4, 6, 7, 9]  # siRNA-1, siRNA-2, siRNA-5, siRNA-7, siRNA-8, siRNA-10
}

group_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
# group_colors = ['brown', '#ae5a21', '#ca6a28', 'yellow']

legend_elements = []

for (group_name, indices), group_color in zip(groups.items(), group_colors):
    for idx in indices:
        result = odeint(model, [R0, F0, N0], times, args=(k_prod_values[idx], k_deg_values[idx], k_fold_values[idx]))
        R_values, F_values, N_values = result.T
        ax5.plot(times, F_values, color=group_color, label=f'Folded shRNA-{idx+1}' if group_name not in [elem.get_label() for elem in legend_elements] else "")
    legend_elements.append(plt.Line2D([0], [0], color=group_color, label=group_name))

ax4.legend(loc='upper left')
ax5.legend(handles=legend_elements, loc='upper right')

ax4.set_yscale('log', basey=10)
ax4.yaxis.set_major_formatter(plt.ScalarFormatter())
ax4.get_yaxis().set_minor_formatter(plt.NullFormatter())

# Set color for ax4
ax4.set_xlabel('Time, min', color=color_value)
ax4.set_ylabel('Yeast culture OD', color=color_value)
ax4.tick_params(axis='y', labelcolor=color_value, colors=color_value)
ax4.tick_params(axis='x', colors=color_value)

# Set color for ax5
ax5.set_ylabel('Folded shRNA, ng/ml', color=color_value)
ax5.tick_params(axis='y', labelcolor=color_value, colors=color_value)

def set_axis_colors(ax):
    ax.set_xlabel(ax.get_xlabel(), color=color_value)
    ax.set_ylabel(ax.get_ylabel(), color=color_value)
    ax.tick_params(axis='y', labelcolor=color_value, colors=color_value)
    ax.tick_params(axis='x', colors=color_value)
    for spine in ax.spines.values():
        spine.set_edgecolor(color_value)

set_axis_colors(ax1)
set_axis_colors(ax2)
set_axis_colors(ax4)
set_axis_colors(ax5)

ax1.text(-0.09, 0.95, 'A', transform=ax1.transAxes, fontsize=24, va='top', color = color_value)
ax4.text(-0.09, 0.95, 'B', transform=ax4.transAxes, fontsize=24, va='top', color = color_value)

plt.tight_layout()
plt.savefig('step1.png', transparent=True)
plt.show()

# Extract the folded siRNA values at the 1440 min timepoint for each siRNA.
final_folded_values = []
for i in range(10):
    result = odeint(model, [R0, F0, N0], times, args=(k_prod_values[i], k_deg_values[i], k_fold_values[i]))
    final_folded_value = result[-1, 1]  # -1 refers to the last timepoint, 1 refers to the Folded shRNA value
    final_folded_values.append(final_folded_value)

print("Folded shRNA amount for each siRNA version:")
for i, value in enumerate(final_folded_values):
    print(f"v{[i]}: {value:.4f}")