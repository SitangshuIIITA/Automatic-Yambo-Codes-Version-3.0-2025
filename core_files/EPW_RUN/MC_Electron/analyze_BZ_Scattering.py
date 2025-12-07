#!/usr/bin/env python3
import os, glob, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

plt.rcParams.update({
    'font.size': 28,
    'axes.linewidth': 1.2,
    'font.family': 'serif',
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'figure.dpi': 300,
    'lines.markersize': 7,
    'legend.handlelength': 2.6
})

# ------------------------------------------------------------
# Read Fermi level (robust)
# ------------------------------------------------------------
fermi_file = "report_parsed.txt"
E_fermi = 0.0

if os.path.isfile(fermi_file):
    with open(fermi_file, "r") as fh:
        for line in fh:
            if "Fermi energy" in line:
                m = re.search(r"=\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)", line)
                if m:
                    try:
                        E_fermi = float(m.group(1))
                        print(f"[INFO] Fermi level loaded: E_F = {E_fermi:.6f} eV")
                    except:
                        print("[WARN] Could not parse Fermi level; using 0.0 eV")
else:
    print("[WARN] report_parsed.txt not found; using E_F = 0.0 eV")


# ------------------------------------------------------------
# Extract field from folder name
# ------------------------------------------------------------
def extract_field(folder):
    m = re.search(r"E_([0-9.+\-eE]+)", folder)
    if m:
        try:
            return float(m.group(1))
        except:
            return None
    return None


# ------------------------------------------------------------
# Scan kcloud folder
# ------------------------------------------------------------
base = "kcloud"
if not os.path.isdir(base):
    print("ERROR: kcloud/ missing.")
    exit(1)

field_dirs = sorted([d for d in glob.glob(os.path.join(base, "E_*")) if os.path.isdir(d)])

if len(field_dirs) == 0:
    print("ERROR: No E_* folders under kcloud/")
    exit(1)

fields = []
for d in field_dirs:
    fv = extract_field(os.path.basename(d))
    if fv is not None:
        fields.append((fv, d))

fields_sorted = sorted(fields, key=lambda x: x[0])

print("[INFO] Fields detected:")
for fv, d in fields_sorted:
    print(f"  {fv:.3e}  ->  {d}")


# ------------------------------------------------------------
# Load CB1.txt distributions (with smart reference detection)
# ------------------------------------------------------------
all_distributions = []
field_values = []
ref_mode = None  # will be "vbm" or "fermi" or "none"

for fval, folder in fields_sorted:

    cb_files = sorted(glob.glob(os.path.join(folder, "CB*.txt")))
    if len(cb_files) == 0:
        print(f"[WARN] No CB*.txt in {folder}, skipping.")
        continue

    cbfile = os.path.join(folder, "CB1.txt")
    if not os.path.isfile(cbfile):
        cbfile = cb_files[0]     # fallback
        print(f"[WARN] CB1.txt missing in {folder}, using {cbfile}")

    ik_list = []
    E_list  = []
    N_list  = []

    with open(cbfile, "r") as fh:
        for line in fh:
            if line.strip().startswith("#") or len(line.strip()) == 0:
                continue
            p = line.split()
            try:
                ik = int(float(p[0]))
                en = float(p[1])
                cnt = int(float(p[2]))
            except:
                continue

            ik_list.append(ik)
            E_list.append(en)
            N_list.append(cnt)

    if len(ik_list) == 0:
        print(f"[WARN] No valid data in {cbfile}, skipping.")
        continue

    iks = np.array(ik_list)
    Es_raw = np.array(E_list)
    Ns = np.array(N_list)

    # ---------- Auto-detection of reference system ----------
    if len(all_distributions) == 0:
        # first field → detect reference
        max_en = np.max(Es_raw)
        vbm_tol = 1e-4

        print(f"[DEBUG] Raw energy range in first field: {Es_raw.min():.6f} to {Es_raw.max():.6f} eV")

        if abs(max_en) <= vbm_tol:
            ref_mode = "vbm"
            print(f"[INFO] Detected VBM-referenced energies (max ≈ 0). Using VBM reference (no shift).")
        else:
            ref_mode = "fermi"
            print(f"[INFO] Detected absolute energies; using Fermi reference: subtract E_F = {E_fermi:.6f} eV")

    # ---------- Apply chosen reference ----------
    if ref_mode == "vbm":
        Es = Es_raw.copy()
    elif ref_mode == "fermi":
        Es = Es_raw - E_fermi
    else:
        Es = Es_raw.copy()

    print(f"[DEBUG] After referencing: min={Es.min():.6f}, max={Es.max():.6f} eV")

    all_distributions.append((iks, Es, Ns))
    field_values.append(fval)


# ------------------------------------------------------------
# Band reference (use already-shifted first field)
# ------------------------------------------------------------
_, Eband_ref, _ = all_distributions[0]

# ------------------------------------------------------------
# Auto-detect K, Γ, M k-point indices from report_parsed.txt
# ------------------------------------------------------------
def detect_K_G_M_indices(parsed_file="report_parsed.txt"):
    if not os.path.isfile(parsed_file):
        print("[WARN] Cannot auto-detect K,Γ,M (report_parsed.txt missing). Using defaults 0,53,100.")
        return 0, 53, 100

    kpoints = []
    with open(parsed_file, "r") as fh:
        read_flag = False
        for line in fh:
            if "K-POINT INDEX MAP" in line:
                read_flag = True
                continue
            if read_flag:
                if line.strip().startswith("---"):
                    continue
                if len(line.strip()) == 0:
                    break
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        ik_pos = int(parts[1])
                        kx = float(parts[2])
                        ky = float(parts[3])
                        kpoints.append((ik_pos, kx, ky))
                    except:
                        pass

    if len(kpoints) == 0:
        print("[WARN] Could not parse kpoint map. Using defaults 0,53,100.")
        return 0, 53, 100

    # Convert list to numpy arrays
    ik_arr = np.array([p[0] for p in kpoints])
    kx_arr = np.array([p[1] for p in kpoints])
    ky_arr = np.array([p[2] for p in kpoints])

    # ---- detect Γ (closest to 0,0)
    G_idx = ik_arr[np.argmin(kx_arr**2 + ky_arr**2)]

    # ---- detect K (closest to 1/3, 1/3)
    K_idx = ik_arr[np.argmin((kx_arr - 1/3)**2 + (ky_arr - 1/3)**2)]

    # ---- detect M (closest to 0.5, 0) OR (0, 0.5)
    M_candidates = np.vstack([
        (kx_arr - 0.5)**2 + (ky_arr - 0.0)**2,
        (kx_arr - 0.0)**2 + (ky_arr - 0.5)**2
    ])
    M_idx = ik_arr[np.argmin(M_candidates.min(axis=0))]

    print(f"[INFO] Auto-detected symmetry points: K={K_idx}, Γ={G_idx}, M={M_idx}")

    return K_idx, G_idx, M_idx

# ------------------------------------------------------------
# Figure setup
# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 5))
fig.tight_layout()

sc = ax.scatter([], [], s=[], edgecolors='k')
band_line, = ax.plot([], [], 'k:', linewidth=1.5)

text_field = ax.text(0.02, 0.89, "", transform=ax.transAxes, fontsize=28)

ax.set_ylabel("Energy (eV)")

all_ik = all_distributions[0][0]
ax.set_xlim(np.min(all_ik), np.max(all_ik))
ax.set_ylim(Eband_ref.min() - 0.1, Eband_ref.max() + 0.1)


# ------------------------------------------------------------
# Replace x-axis ticks with K, Γ, M
# (Indices from your parsed k-mesh)
# ------------------------------------------------------------
K_idx, G_idx, M_idx = detect_K_G_M_indices("report_parsed.txt")


ax.set_xticks([K_idx, G_idx, M_idx])
ax.set_xticklabels(["K", r"$\Gamma$", "M"])

ax.axvline(K_idx, color="gray", linestyle="--", linewidth=1.2, alpha=0.8)
ax.axvline(G_idx, color="gray", linestyle="--", linewidth=1.2, alpha=0.8)
ax.axvline(M_idx, color="gray", linestyle="--", linewidth=1.2, alpha=0.8)


# ------------------------------------------------------------
# Colors
# ------------------------------------------------------------
cmap = plt.get_cmap("viridis")
colors = cmap(np.linspace(0, 1, len(all_distributions)))


# ------------------------------------------------------------
# Frame update
# ------------------------------------------------------------
def update(i):
    iks, Es, Ns = all_distributions[i]

    if np.max(Ns) > 0:
        sizes = (Ns / np.max(Ns)) * 200
    else:
        sizes = np.full_like(Ns, 10)

    sc.set_offsets(np.column_stack((iks, Es)))
    sc.set_sizes(sizes)
    sc.set_color(colors[i])

    band_line.set_data(iks, Eband_ref)

    # ---- Scientific notation with SAME significant digits in kV/cm ----
 # ---- Scientific notation with FIXED 3 decimal places in mantissa (e.g. 1.000×10^{0}) ----
    fv = field_values[i] * 1e-5   # convert V/m → kV/cm

# Convert to scientific notation with controlled digits
    s = f"{fv:.6e}"               # e.g. "1.000000e+00", "3.728000e+02"
    mant_str, exp_str = s.split("e")

# Convert mantissa to float to avoid issues, then format with EXACT decimals:
    mant = float(mant_str)
    mant_formatted = f"{mant:.3f}"   # <-- FIXED DECIMAL PLACES (3 decimals)

    exp = int(exp_str)

# Final LaTeX-style engineering notation
    sci_str = fr"{mant_formatted}×10$^{{{exp}}}$"

    text_field.set_text(f"Field = {sci_str} kV/cm")


    return sc, band_line, text_field



# ------------------------------------------------------------
# Create GIF
# ------------------------------------------------------------
anim = FuncAnimation(fig, update, frames=len(all_distributions),
                     interval=900, blit=False)

writer = PillowWriter(fps=2)
anim.save("kcloud_animation.gif", writer=writer)

print("[OK] Saved animation: kcloud_animation.gif")

