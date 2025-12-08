#!/usr/bin/env python3
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm

base_dir = os.getcwd()
hself_dir = os.path.join(base_dir, "Hself_LW")

# ============================================================
# Constants
# ============================================================
meV_to_THz = 0.15192674  # convert meV → THz
kB = 8.617333262e-5      # eV/K
bar_width = 0.10
Nx = 1500

# ============================================================
# Locate all temperature subfolders
# ============================================================
subdirs = sorted(glob.glob(os.path.join(hself_dir, "Hself_LW-*")))
if not subdirs:
    raise FileNotFoundError(f"No Hself_LW-* folders found in {hself_dir}")

all_T, all_modes, all_data = [], set(), []

# ============================================================
# Loop over each temperature folder
# ============================================================
for sub in subdirs:
    try:
        T_str = sub.split("-")[-1]
        T = float(T_str)
    except ValueError:
        print(f"⚠️ Skipping folder {sub} (invalid T format)")
        continue

    file_pattern = os.path.join(sub, f"quadrupole.linewidth.hlself.{T:.3f}K")
    matches = glob.glob(file_pattern)

    if not matches:
        print(f"❌ File not found for {T}K: {file_pattern}")
        continue

    fname = matches[0]
    data = np.loadtxt(fname)
    E = data[:, 2]        # eV
    imode = data[:, 3].astype(int)
    imsigma = data[:, 4]  # meV
    rate = 2 * imsigma * meV_to_THz

    VBM = np.max(E)
    E_window = 3 * kB * T
    filtered = (E >= VBM - E_window) & (E <= VBM)
    if np.sum(filtered) == 0:
        print(f"⚠️ No states in 3kBT window for {T}K.")
        continue

    E_rel = E[filtered] - VBM
    f_occ = 1 / (np.exp(E_rel / (kB * T)) + 1)
    all_modes.update(np.unique(imode[filtered]))
    all_data.append((T, imode[filtered], rate[filtered], f_occ))

all_T = np.array(sorted([d[0] for d in all_data]))
all_modes = np.array(sorted(all_modes))
x_dense = np.linspace(all_modes[0] - 0.5, all_modes[-1] + 0.5, Nx)

# ============================================================
# Build Z matrix (Temperature × Phonon mode)
# ============================================================
Z = np.zeros((len(all_T), Nx))

for iT, (T, modes_T, rate_T, f_occ_T) in enumerate(all_data):
    for m, r, w in zip(modes_T, rate_T, f_occ_T):
        profile = np.exp(-((x_dense - m) ** 2) / (2 * bar_width ** 2))
        Z[iT, :] += r * w * profile

# ============================================================
# Plot spectral map
# ============================================================
plt.rcParams.update({'font.size': 16, 'axes.linewidth': 1.2, 'font.family': 'serif'})
fig, ax = plt.subplots(figsize=(10, 6))

norm = Normalize(vmin=0, vmax=np.max(Z))
im = ax.imshow(
    Z, origin='lower', aspect='auto',
    extent=[x_dense[0], x_dense[-1], all_T[0], all_T[-1]],
    cmap=cm.cividis, norm=norm
)

ax.set_xlabel("Phonon Mode Index")
ax.set_ylabel("Temperature (K)")
ax.set_xticks(all_modes)
#plt.title("Fermi-weighted Scattering Rate per Phonon Mode (3 kB*T below VBM)")

cbar = plt.colorbar(im, ax=ax, pad=0.02)
cbar.set_label(r"Fermi-weighted Scattering Rate $\tau^{-1}$ (THz)")

plt.tight_layout()

# ============================================================
# Save plot and summary
# ============================================================
out_path = os.path.join(base_dir, "VBM_quadrupole_scattering_phonon_mode_vs_T.png")
plt.savefig(out_path, dpi=400)

summary_text = """
This Python script generates temperature-dependent, phonon-mode-resolved scattering maps
for hole self-energies obtained from EPW calculations. It automatically searches
subdirectories "Hself_LW-{T}" for files named "quadrupole.linewidth.hlself.{T}K", extracts
the energy (E_i), phonon mode index (m_i), and linewidth (ImΣ_i), and converts them into
scattering rates in terahertz (THz):

    τ_i^{-1} = 2 * ImΣ_i

Only states within an energy window of 3k_B*T below the valence-band maximum (VBM)
are considered:

    E_VBM - 3k_B*T ≤ E_i ≤ E_VBM

These states are weighted using the Fermi–Dirac distribution:

    f_i = 1 / (exp((E_i - E_VBM)/(k_B*T)) + 1)

Each phonon mode’s contribution is smoothed using a Gaussian smearing of width σ,
producing a continuous spectral distribution over the mode index. The final weighted
intensity map is computed as:

    Z(T, m) = Σ_i [ τ_i^{-1} * f_i * exp(-(m - m_i)^2 / (2σ^2)) ]

The resulting 2D colormap Z(T, m) reveals how Fermi-weighted scattering rates vary with
both temperature and phonon mode index. The final figure is saved as
"quadrupole_phonon_mode_vs_T_spectral_FD_smooth.png" in the working directory.
"""

# Save the summary
summary_file = "vbm_scattering_summary.txt"
with open(summary_file, "w") as f:
    f.write(summary_text)

print("\n✅ Scattering summary written successfully.")
print(f"📝 Summary saved as: {summary_file}")

