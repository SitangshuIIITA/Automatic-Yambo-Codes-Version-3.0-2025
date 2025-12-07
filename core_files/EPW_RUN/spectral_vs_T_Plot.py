#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os

# ---------------------------------------------------
# File path (relative to EPW_plots)
# ---------------------------------------------------
folder = "e-Mobility_vs_T/e-Mobility_quadrupole"
fname = os.path.join(folder, "inv_taucb_freq_0.fmt")

if not os.path.exists(fname):
    raise FileNotFoundError(f"File '{fname}' not found in {folder}!")

# ---------------------------------------------------
# Read formatted data
# ---------------------------------------------------
data = []
with open(fname, 'r') as f:
    for line in f:
        if line.strip().startswith('#') or not line.strip():
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            kpt = int(parts[0])
            ibnd = int(parts[1])
            energy_Ry = float(parts[2])
            freq_meV = float(parts[3])
            tau_Ry = float(parts[4])
            data.append([kpt, ibnd, energy_Ry, freq_meV, tau_Ry])
        except ValueError:
            continue

data = np.array(data)
if data.size == 0:
    raise ValueError("No valid numerical data found in file!")

# ---------------------------------------------------
# Unit conversion
# ---------------------------------------------------
factor = 20670.6944033  # Ry time to ps
tau_ps = data[:, 4] * factor
inv_tau_ps = 1.0 / tau_ps
inv_tau_ps[np.isinf(inv_tau_ps)] = np.nan

# ---------------------------------------------------
# Unique k-points and bands
# ---------------------------------------------------
unique_k = np.unique(data[:, 0])
unique_b = np.unique(data[:, 1])

# Choose k-point index (example)
k_select = 1
mask_k = data[:, 0] == k_select

# ---------------------------------------------------
# Plot 1: Inverse lifetime vs. phonon frequency
# ---------------------------------------------------
plt.figure(figsize=(8, 6))
colors = plt.cm.plasma(np.linspace(0, 1, len(unique_b)))

for i, ib in enumerate(unique_b):
    mask = mask_k & (data[:, 1] == ib)
    if np.any(mask):
        freq = data[mask, 3]
        invtau = inv_tau_ps[mask]
        idx = np.argsort(freq)
        freq = freq[idx]
        invtau = invtau[idx]
        plt.plot(freq, invtau, lw=2.0, color=colors[i], label=f'Band {int(ib)}')

plt.xlabel("Phonon frequency (meV)", fontsize=14)
plt.ylabel(r"$1/\tau$ (1/ps)", fontsize=14)
plt.title(f"Hole inverse lifetime vs phonon energy (k = {k_select})", fontsize=15)
plt.legend(frameon=False)
plt.grid(True, ls='--', alpha=0.4)
plt.tight_layout()
plt.savefig(f"hole_inv_tau_k{k_select}.png", dpi=300)

# ---------------------------------------------------
# Plot 2: Derivative d(1/tau)/dw
# ---------------------------------------------------
plt.figure(figsize=(8, 6))
for i, ib in enumerate(unique_b):
    mask = mask_k & (data[:, 1] == ib)
    if np.any(mask):
        freq = data[mask, 3]
        invtau = inv_tau_ps[mask]
        idx = np.argsort(freq)
        freq = freq[idx]
        invtau = invtau[idx]
        derivative = np.gradient(invtau, freq)
        plt.plot(freq, derivative, lw=2.0, color=colors[i], label=f'Band {int(ib)}')

plt.xlabel("Phonon frequency (meV)", fontsize=14)
plt.ylabel(r"$\partial(1/\tau) / \partial\omega$ (1/ps·meV$^{-1}$)", fontsize=14)
plt.title(f"Spectral derivative at k = {k_select}", fontsize=15)
plt.legend(frameon=False)
plt.grid(True, ls='--', alpha=0.4)
plt.tight_layout()
plt.savefig(f"hole_dinv_tau_dw_k{k_select}.png", dpi=300)
plt.show()

print(f"✅ Saved: hole_inv_tau_k{k_select}.png and hole_dinv_tau_dw_k{k_select}.png in current folder.")

