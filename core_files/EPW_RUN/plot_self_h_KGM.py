#!/usr/bin/env python3
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# ------------------------------
# Current folder
# ------------------------------
path = "./"  # EPW_plots folder

# ------------------------------
# Read global bash variable "temps"
# ------------------------------
temps_str = os.environ.get("temps")
if temps_str is None:
    raise RuntimeError("Environment variable 'temps' not set. Please source the bash script first.")
temps = [int(t) for t in temps_str.split()]
print("Temperatures found:", temps)

# ------------------------------
# Read VBM, CBM from wannier_dis_details.txt
# ------------------------------
dis_file = os.path.join(path, "wannier_dis_details.txt")
VBM = None
CBM = None
with open(dis_file, "r") as f:
    for line in f:
        if "Detected VBM, CBM from nscf" in line:
            parts = line.strip().split(":")[1].split(",")
            VBM = float(parts[0])
            CBM = float(parts[1])
            break

if VBM is None or CBM is None:
    raise RuntimeError("Could not find VBM/CBM in wannier_dis_details.txt")
print(f"Detected VBM = {VBM:.4f}, CBM = {CBM:.4f} from nscf")

# ------------------------------
# Read klist.dat to get k-points
# ------------------------------
klist_file = os.path.join(path, "klist.dat")
kcoords = []
with open(klist_file, "r") as f:
    for line in f:
        if "coord.:" in line:
            parts = line.split()
            ik = int(parts[2])
            x = float(parts[4])
            y = float(parts[5])
            z = float(parts[6])
            kcoords.append((ik, x, y, z))

kcoords = np.array(kcoords)
unique_k = kcoords[:, 0].astype(int)
n_k = len(unique_k)

# ------------------------------
# Compute cumulative k-path distance
# ------------------------------
dist = np.zeros(n_k)
for i in range(1, n_k):
    dk = np.linalg.norm(kcoords[i, 1:] - kcoords[i - 1, 1:])
    dist[i] = dist[i - 1] + dk

x_path = dist  # cumulative path distance

# ------------------------------
# Identify high-symmetry points using exact k-coordinates
# ------------------------------
K_ik = unique_k[0]           # first k-point
M_ik = unique_k[-1]          # last k-point
Gamma_idx = np.argmin(np.linalg.norm(kcoords[:, 1:], axis=1))  # closest to origin
Gamma_ik = int(kcoords[Gamma_idx, 0])

# Coordinates
K_coord = kcoords[unique_k == K_ik, 1:][0]
Gamma_coord = kcoords[unique_k == Gamma_ik, 1:][0]
M_coord = kcoords[unique_k == M_ik, 1:][0]

# Euclidean distances along BZ path
K_Gamma_dist = np.linalg.norm(Gamma_coord - K_coord)
Gamma_M_dist = np.linalg.norm(M_coord - Gamma_coord)

print(f"High-symmetry points: K={K_ik}, Γ={Gamma_ik}, M={M_ik}")
print(f"K → Γ distance = {K_Gamma_dist:.6f}")
print(f"Γ → M distance = {Gamma_M_dist:.6f}")

# Positions along cumulative path
highsym_positions = [x_path[unique_k == K_ik][0],
                     x_path[unique_k == Gamma_ik][0],
                     x_path[unique_k == M_ik][0]]
highsym_labels = ['K', r'$\Gamma$', 'M']

# ------------------------------
# Loop over temperatures
# ------------------------------
for T in temps:
    folder = os.path.join("Hself_LW", f"Hself_LW-{T}")
    pattern = os.path.join(folder, f"*.linewidth.hlself.{T}.000K")
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"No linewidth files found for T={T} in {folder}. Skipping.")
        continue

    for filename in files:
        data = np.loadtxt(filename, comments="#")
        ik_data = data[:, 0].astype(int)
        ibnd_data = data[:, 1].astype(int)
        E_data = data[:, 2]  # Ry → eV (convert if needed)
        ImSigma_data = data[:, 4]  # meV

        # Shift energies: VBM → 0, then add CBM-VBM offset
        E_data = E_data #+  ( -CBM + VBM)

        # Filter along full path
        mask = (ik_data >= unique_k[0]) & (ik_data <= unique_k[-1])
        ik_data = ik_data[mask]
        ibnd_data = ibnd_data[mask]
        E_data = E_data[mask]
        ImSigma_data = ImSigma_data[mask]

        unique_bands = np.unique(ibnd_data)
        n_b = len(unique_bands)

        # Build E_map and Sigma_map
        E_map = np.zeros((n_b, n_k))
        Sigma_map = np.zeros_like(E_map)
        k_to_idx = {k: i for i, k in enumerate(unique_k)}

        for i, b in enumerate(unique_bands):
            for k in np.unique(ik_data):
                sel = (ibnd_data == b) & (ik_data == k)
                if np.any(sel):
                    idx = k_to_idx[k]
                    E_map[i, idx] = np.mean(E_data[sel])
                    Sigma_map[i, idx] = np.sum(ImSigma_data[sel])

        # ------------------------------
        # Plot bands with Σ
        # ------------------------------
        plt.rcParams.update({
            'font.size': 22,
            'axes.linewidth': 0.8,
            'font.family': 'serif',
            'xtick.direction': 'in',
            'ytick.direction': 'in',
            'xtick.top': True,
            'ytick.right': True,
        })
        fig, ax = plt.subplots(figsize=(12, 5))
        cmap = plt.cm.Oranges

        for i in range(n_b):
            E_fine = np.interp(np.linspace(0, n_k - 1, 800),
                               np.arange(n_k), E_map[i, :])
            Sigma_fine = np.interp(np.linspace(0, n_k - 1, 800),
                                   np.arange(n_k), Sigma_map[i, :])
            x_fine = np.interp(np.linspace(0, n_k - 1, 800),
                               np.arange(n_k), x_path)

            Sigma_norm = (Sigma_fine / np.max(Sigma_fine)
                          if np.max(Sigma_fine) > 0 else np.zeros_like(Sigma_fine))
            lw = 3.0 + 5.0 * (Sigma_norm ** 1e-3)
            color = cmap(Sigma_norm)

            points = np.array([x_fine, E_fine]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc = LineCollection(segments, colors=color, linewidths=lw)
            ax.add_collection(lc)

            # Dashed band line
            ax.plot(x_path, E_map[i, :], color='k', linestyle='--',
                    linewidth=0.6, alpha=0.4)

        # Colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap)
        sm.set_array(Sigma_map)
        cbar = plt.colorbar(sm, ax=ax, label=r'Im$\Sigma_{n\mathbf{k}}\left(\omega,T\right)$ (meV)')

        # Axis
        ax.set_ylabel('Energy (eV)', fontsize=22)
        ax.set_xticks(highsym_positions)
        ax.set_xticklabels(highsym_labels)
        ax.set_xlim(x_path[0], x_path[-1])

        # Vertical line at Γ
        ax.axvline(x=highsym_positions[1], color='k', linestyle='--', linewidth=1.0, alpha=0.8)

        plt.tight_layout()

        outname = os.path.join(folder, f"{os.path.basename(filename)}.png")
        plt.savefig(outname, dpi=300)
        plt.close()
        print("Saved:", outname)

