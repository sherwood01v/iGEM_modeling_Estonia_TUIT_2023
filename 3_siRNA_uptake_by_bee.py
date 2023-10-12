import numpy as np
import matplotlib.pyplot as plt

# Constants
S_0_values = [
    43.54770403, 42.60853786, 10.65261762, 16.00765907, 42.66378902,
    26.68768909, 42.67119476, 42.65832787, 26.64382571, 42.56382243
]  # Updated list of S_0 values
k_int_const = 0.045
k_deg_const = 0.084

# Function to compute M(t)
def compute_Mt(k_int, k_deg, S_0, t):
    return k_int * S_0 * t * np.exp(-k_deg * t)

# Time range (e.g., from 0 to 500, with 1000 points)
t_range = np.linspace(0, 200, 1000)

# Define font properties
font = {
    'family': 'serif',
    'color': 'black',
    'weight': 'normal',
    'size': 14,
}

# Plot
plt.rcParams.update({'font.size': 18, 'axes.labelsize': 18, 'axes.titlesize': 18, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 12})
plt.figure(figsize=(12, 8))

colors = plt.cm.cividis(np.linspace(0, 1, len(S_0_values)))
color_value = '#ffe8c0'

peak_values = []

for index, S_0 in enumerate(S_0_values):
    M_t_values = compute_Mt(k_int_const, k_deg_const, S_0, t_range)
    max_value = np.max(M_t_values)
    peak_values.append(max_value)
    plt.plot(t_range, M_t_values, color=colors[index], label=f'siRNA-{index+1}')

plt.xlabel('Time (min)', fontdict=font, color=color_value)
plt.ylabel('siRNA concentration (nM)', fontdict=font, color=color_value)
plt.legend()
plt.grid(False)

ax = plt.gca()
for spine in ax.spines.values():
    spine.set_edgecolor(color_value)

ax.tick_params(axis='both', colors=color_value)

plt.tight_layout()
plt.savefig('step3.png', transparent=True)
plt.show()

# Print peak values
for idx, peak in enumerate(peak_values, 1):
    print(f'Peak value for siRNA-{idx}')