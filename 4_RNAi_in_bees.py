from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def model(y, t, k1, k2, kv):
    S, A, R, V = y
    dR_dt = k1 * S * A - k2 * R * V
    dV_dt = kv * V - k2 * R * V
    dS_dt = -k1 * S * A
    dA_dt = -k1 * S * A
    return [dS_dt, dA_dt, dR_dt, dV_dt]

# Initial concentrations
# A0 = 33696  # nM in BEES
A0 = 6335  #nM in YEAST
R0 = 0.0  # initial RISC is zero
conc_microgram_per_m3 = 2.044e-20
avogadro_number = 6.022e23
molecular_weight_g_per_mol = 2.3e3
V0 = (conc_microgram_per_m3 * avogadro_number) / (molecular_weight_g_per_mol * 1e9)

# Rate constants
k1 = 10 * 1000 * 60
k2 = 0.426
kv = 0.0125 #IN BEES
kv = 0.0280 #IN YEAST

# Time array
t = np.linspace(0, 150, 1000)

# List of siRNA concentrations
# siRNA_concentrations = [8.58, 8.40, 2.10, 3.15, 8.41, 5.26, 8.41, 8.41, 5.25, 8.39]
# siRNA_concentrations = [43.54770403, 42.60853786, 10.65261762, 16.00765907, 42.66378902, 26.68768909, 42.67119476, 42.65832787, 26.64382571, 42.56382243]

plt.rcParams.update({'font.size': 18, 'axes.labelsize': 18, 'axes.titlesize': 18, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 12})
plt.figure(figsize=(10,6))
color_value = '#ffe8c0'

# Loop through each siRNA concentration, solve the model, and plot the result
# for S0 in siRNA_concentrations:
#     result = odeint(model, [S0, A0, R0, V0], t, args=(k1, k2, kv))
#     S, A, R, V = result.T
#     plt.plot(t, V, label=f'siRNA = {S0:.2f} nM')
for i in np.arange(0.1, 20, 0.1):
    if i==0.1 or i==0.2 or i==0.6 or i==2.0 or i==12.0 or i==16.8:
        S0 = i
        result = odeint(model, [S0, A0, R0, V0], t, args=(k1, k2, kv))
        S, A, R, V = result.T
        plt.plot(t, V, label=f'siRNA = {S0:.2f} nM')

ax = plt.gca()
for spine in ax.spines.values():
    spine.set_edgecolor(color_value)
ax.tick_params(axis='both', colors=color_value)

plt.xlabel("Time (min)", color=color_value)
plt.ylabel("Concentration of viral RNA (nM)", color=color_value)
plt.grid(False)
plt.legend()
plt.tight_layout()
plt.savefig('step4a.png', transparent=True)
plt.show()