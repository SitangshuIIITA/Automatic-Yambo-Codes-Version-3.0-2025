#!/usr/bin/env python3
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# ============================================================
# Base directory structure (same as the spectral map script)
# ============================================================
base_dir = os.getcwd()
hself_dir = os.path.join(base_dir, "Hself_LW")

# ============================================================
# Constants
# ============================================================
meV_to_THz = 0.15192674  # convert meV → THz
kB = 8.617333262e-5      # eV/K
temps_to_plot = [77, 300, 500]  # Temperatures to compare

# ============================================================
# Helper: Compute Fermi-weighted scattering per phonon mode
# ============================================================
def compute_weighted_scattering(T):
    subdir = os.path.join(hself_dir, f"Hself_LW-{T}")
    file_pattern = os.path.join(subdir, f"quadrupole.linewidth.hlself.{T:.3f}K")
    matches = glob.glob(file_pattern)

    if not matches:
        print(f"❌ File not found for {T}K: {file_pattern}")
        return None, None

    data = np.loadtxt(matches[0])
    E = data[:, 2]        # eV
    imode = data[:, 3].astype(int)
    imsigma = data[:, 4]  # meV
    rate = 2 * imsigma * meV_to_THz  # scattering rate in THz

    # Define energy window near VBM
    VBM = np.max(E)
    E_window = 3 * kB * T
    filtered = (E >= VBM - E_window) & (E <= VBM)
    if np.sum(filtered) == 0:
        print(f"⚠️ No states in 3kBT window for {T}K.")
        return None, None

    # Fermi-Dirac weighting
    E_rel = E[filtered] - VBM
    f_occ = 1 / (np.exp(E_rel / (kB * T)) + 1)

    # Weighted *sum* of scattering per mode (matches cmap)
    unique_modes = np.unique(imode[filtered])
    weighted_rate = []
    for m in unique_modes:
        mask = imode[filtered] == m
        wsum = np.sum(rate[filtered][mask] * f_occ[mask])
        weighted_rate.append(wsum)

    return unique_modes, np.array(weighted_rate)


# ============================================================
# Compute for selected temperatures
# ============================================================
all_data = {}
for T in temps_to_plot:
    modes, rate_T = compute_weighted_scattering(T)
    if modes is not None:
        all_data[T] = (modes, rate_T)

if not all_data:
    raise RuntimeError("No valid scattering data found for requested temperatures.")

# ============================================================
# Plot: Bar chart comparison (center middle temperature)
# ============================================================
plt.rcParams.update({
    'font.size': 16,
    'axes.linewidth': 1.5,
    'font.family': 'serif'
})

fig, ax = plt.subplots(figsize=(12, 5))
colors = cm.cividis(np.linspace(0.2, 0.9, len(temps_to_plot)))

bar_width = 0.25
n_temps = len(temps_to_plot)

# Shift bars so that the *middle temperature* is centered on the tick
offsets = np.arange(n_temps) - (n_temps - 1) / 2
offsets = offsets * bar_width  # e.g. [-0.15, 0.0, +0.15] for 3 temps

for i, (T, offset) in enumerate(zip(temps_to_plot, offsets)):
    if T in all_data:
        modes, sum_rates = all_data[T]
        ax.bar(
            modes + offset,
            sum_rates,
            width=bar_width,
            label=f"{T} K",
            color=colors[i],
            alpha=0.9
        )

ax.set_xlabel("Phonon Mode Index", fontsize=18)
ax.set_ylabel(r"$\sum_i \tau_i^{-1} f_i$ (THz)", fontsize=18)
#ax.set_title("Mode-Resolved Scattering Rate Summation (3$k_BT$ below VBM)", fontsize=18)
# ============================================================
# Legend (with extra text note)
# ============================================================
legend_labels = [f"{T} K" for T in temps_to_plot]
bars = [plt.Rectangle((0, 0), 1, 1, color=colors[i], alpha=0.9) for i in range(len(temps_to_plot))]
bars.append(plt.Rectangle((0, 0), 0, 0, alpha=0))  # dummy entry for note
legend_labels.append(r"All energies within" "\n" r"$3k_BT$ window below VBM")

legend = ax.legend(
    bars, legend_labels,
    loc='upper left',
    fontsize=16,
    frameon=True
)
legend.get_frame().set_alpha(0.4)



plt.tight_layout()
out_path = os.path.join(base_dir, "VBM_scattering_quadrupole_mode_bar_77K_300K_500K.png")
plt.savefig(out_path, dpi=400)
#plt.show()

print(f"\n✅ Bar chart saved as: {out_path}")

