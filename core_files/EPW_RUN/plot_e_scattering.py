#!/usr/bin/env python3
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm

# ============================================================
# Base directory
# ============================================================
base_dir = os.getcwd()
eself_dir = os.path.join(base_dir, "Elself_LW")  # electron self-energy directory

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
subdirs = sorted(glob.glob(os.path.join(eself_dir, "Elself_LW-*")))
if not subdirs:
    raise FileNotFoundError(f"No Elself_LW-* folders found in {eself_dir}")

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

    file_pattern = os.path.join(sub, f"quadrupole.linewidth.elself.{T:.3f}K")
    matches = glob.glob(file_pattern)

    if not matches:
        print(f"❌ File not found for {T}K: {file_pattern}")
        continue

    fname = matches[0]
    data = np.loadtxt(fname)
    E = data[:, 2]        # eV
    imode = data[:, 3].astype(int)
    imsigma = data[:, 4]  # meV
    rate = 2 * imsigma * meV_to_THz  # scattering rate in THz

    # Find CBM and energy window (3kBT above CBM)
    CBM = np.min(E)
    E_window = 3 * kB * T
    filtered = (E >= CBM) & (E <= CBM + E_window)
    if np.sum(filtered) == 0:
        print(f"⚠️ No states in 3kBT window for {T}K.")
        continue

    # Fermi-Dirac weighting (for electrons)
    E_rel = E[filtered] - CBM
    f_occ = 1 / (np.exp(-E_rel / (kB * T)) + 1)  # electrons occupy above CBM

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
#ax.set_title("Electron Scattering Rates (3 $k_B T$ above CBM)")

cbar = plt.colorbar(im, ax=ax, pad=0.02)
cbar.set_label(r"Fermi-weighted Scattering Rate $\tau^{-1}$ (THz)")

plt.tight_layout()

# ============================================================
# Save plot and summary
# ============================================================
out_path = os.path.join(base_dir, "CBM_quadrupole_scattering_phonon_mode_vs_T.png")
plt.savefig(out_path, dpi=400)

summary_text = r"""
This Python script generates temperature-dependent, phonon-mode-resolved scattering maps
for electron self-energies obtained from EPW calculations. It automatically searches
subdirectories "Eself_LW-{T}" for files named "quadrupole.linewidth.elself.{T}K", extracts
the energy ($E_i$), phonon mode index ($m_i$), and linewidth (ImΣ$_i$), and converts them
into scattering rates in terahertz (THz):

    $\tau_i^{-1} = 2 \, \mathrm{Im}\Sigma_i$

Only states within an energy window of $3 k_B T$ above the conduction-band minimum (CBM)
are considered:

    $E_\mathrm{CBM} \le E_i \le E_\mathrm{CBM} + 3 k_B T$

These states are weighted using the Fermi–Dirac distribution:

    $f_i = \frac{1}{1 + \exp(-(E_i - E_\mathrm{CBM})/(k_B T))}$

Each phonon mode’s contribution is smoothed using a Gaussian of width σ, producing a
continuous spectral distribution over the mode index. The final weighted intensity map is:

    $Z(T, m) = \sum_i [ \tau_i^{-1} f_i \exp(-(m - m_i)^2 / (2\sigma^2)) ]$

The resulting colormap $Z(T, m)$ reveals how Fermi-weighted electron–phonon scattering
rates vary with temperature and phonon mode index.
"""

summary_file = "cbm_scattering_summary.txt"
with open(summary_file, "w") as f:
    f.write(summary_text)

print("\n✅ Electron scattering summary written successfully.")
print(f"📝 Summary saved as: {summary_file}")
print(f"📊 Plot saved as: {out_path}")

