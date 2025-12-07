#!/usr/bin/env python3
"""
Plotting script for population_out/E_*.txt files.

Creates:
  1. ik vs Enk bubble plot (weights = normalized Abs/Em)
  2. Histogram: normalized counts vs ik

Now includes:
 • Energy shifted by Fermi level.
 • AUTO high-symmetry x-axis (K → Γ → M path).
"""
# NOTE:
# All phonon mode numbers shown in legends and labels are 1-based.
# (The internal computation still uses 0-based indices, but displayed
# mode numbers are shifted by +1 for human readability.)


import os
import re
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------- STYLE ----------
plt.rcParams.update({
    'font.size': 28,
    'axes.linewidth': 1.2,
    'font.family': 'serif',
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'xtick.bottom': True,
    'ytick.right': True,
    'figure.dpi': 300,
    'lines.markersize': 7,
    'legend.handlelength': 2.6
})

# ------------------------------------------------------------
# Read Fermi level
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
# Read parsed k-point map from report_parsed.txt (robust)
# returns list of tuples (ik_pos, kx, ky)
# ------------------------------------------------------------
def read_kpoint_map(parsed_file="report_parsed.txt"):
    if not os.path.isfile(parsed_file):
        return []

    kpoints = []
    with open(parsed_file, "r") as fh:
        read_flag = False
        for line in fh:
            if "K-POINT INDEX MAP" in line:
                read_flag = True
                continue
            if not read_flag:
                continue
            if line.strip().startswith("---") or len(line.strip()) == 0:
                # continue reading until blank section ends
                # keep reading; blank lines inside map are ignored
                continue
            parts = line.split()
            # Expecting lines like:  <label>  ik  kx  ky  ...
            if len(parts) >= 4:
                try:
                    ik_pos = int(parts[1])
                    kx = float(parts[2])
                    ky = float(parts[3])
                    kpoints.append((ik_pos, kx, ky))
                except:
                    continue
    return kpoints


# ------------------------------------------------------------
# Auto-detect K, Γ, M indices (existing behaviour)
# ------------------------------------------------------------
def detect_K_G_M_indices(parsed_file="report_parsed.txt"):
    kpoints = read_kpoint_map(parsed_file)
    if len(kpoints) == 0:
        print("[WARN] Cannot auto-detect K,Γ,M. Using defaults.")
        return 0, 53, 100

    ik_arr = np.array([p[0] for p in kpoints])
    kx_arr = np.array([p[1] for p in kpoints])
    ky_arr = np.array([p[2] for p in kpoints])

    G_idx = ik_arr[np.argmin(kx_arr**2 + ky_arr**2)]
    K_idx = ik_arr[np.argmin((kx_arr - 1/3)**2 + (ky_arr - 1/3)**2)]

    M_candidates = np.vstack([
        (kx_arr - 0.5)**2 + (ky_arr - 0.0)**2,
        (kx_arr - 0.0)**2 + (ky_arr - 0.5)**2
    ])
    M_idx = ik_arr[np.argmin(M_candidates.min(axis=0))]

    print(f"[INFO] Auto-detected K={K_idx}, Γ={G_idx}, M={M_idx}")
    return K_idx, G_idx, M_idx


# ------------------------------------------------------------
# Convert raw ik → symmetry axis coordinate (robust)
# Uses actual kx,ky coordinates from report_parsed.txt and
# computes cumulative distance along the k-point path, then
# normalizes to [0,1]. Returns (s_array, pts) where pts is
# list of (label, normalized_position).
# ------------------------------------------------------------
def build_symmetry_axis(df, parsed_file="report_parsed.txt"):
    iks = df["ik"].values.astype(int)

    kpoints = read_kpoint_map(parsed_file)
    if len(kpoints) == 0:
        # fallback: linear mapping based on ik min/max
        pts = [("K", int(np.min(iks))), ("$\Gamma$", int((np.min(iks)+np.max(iks))//2)), ("M", int(np.max(iks)))]
        s = (iks - np.min(iks)).astype(float)
        if s.ptp() == 0:
            return np.zeros_like(s), pts
        s = (s - s.min()) / s.ptp()
        return s, pts

    # build dict of ik->(kx,ky)
    kp_dict = {int(p[0]): (float(p[1]), float(p[2])) for p in kpoints}

    # we will consider path from first_k to last_k present in map
    sorted_map = sorted(kpoints, key=lambda x: x[0])  # sort by ik index
    idxs = [int(p[0]) for p in sorted_map]
    coords = [(float(p[1]), float(p[2])) for p in sorted_map]

    # compute cumulative distances along the sorted map
    cum = [0.0]
    for i in range(1, len(coords)):
        dx = coords[i][0] - coords[i-1][0]
        dy = coords[i][1] - coords[i-1][1]
        dist = math.hypot(dx, dy)
        cum.append(cum[-1] + dist)
    cum = np.array(cum)

    total = cum[-1] - cum[0]
    if total == 0.0:
        # degenerate -> fallback linear
        pts = [("K", idxs[0]), ("$\Gamma$", idxs[len(idxs)//2]), ("M", idxs[-1])]
        s = (iks - idxs[0]).astype(float)
        if s.ptp() == 0:
            return np.zeros_like(s), pts
        s = (s - s.min()) / s.ptp()
        return s, pts

    # map every ik in the map to its cumulative distance
    ik_to_cum = {idxs[i]: cum[i] for i in range(len(idxs))}

    # For ik values that are *between* indices in the map, we interpolate
    # We'll create arrays of map indices and cum values for interpolation
    map_idx_arr = np.array(idxs)
    map_cum_arr = np.array(cum)

    # compute s for each requested ik in df:
    s_vals = []
    for ik in iks:
        if ik in ik_to_cum:
            s_val = ik_to_cum[ik]
        else:
            # interpolate: find where ik would fit in map_idx_arr
            # use nearest two map indices around ik
            pos = np.searchsorted(map_idx_arr, ik)
            if pos == 0:
                s_val = map_cum_arr[0]
            elif pos >= len(map_idx_arr):
                s_val = map_cum_arr[-1]
            else:
                i0 = pos - 1
                i1 = pos
                x0 = map_idx_arr[i0]; x1 = map_idx_arr[i1]
                y0 = map_cum_arr[i0]; y1 = map_cum_arr[i1]
                if x1 == x0:
                    s_val = y0
                else:
                    frac = (ik - x0) / (x1 - x0)
                    s_val = y0 + frac * (y1 - y0)
        s_vals.append(s_val)

    s_arr = np.array(s_vals)
    # normalize to [0,1] using min/max of the mapped path
    s_norm = (s_arr - map_cum_arr[0]) / (map_cum_arr[-1] - map_cum_arr[0])

    # build pts list with normalized positions for K, Γ, M
    # detect K,G,M indices using detect function
    K_idx, G_idx, M_idx = detect_K_G_M_indices(parsed_file)
    tick_positions = []
    for lab, idx in [("K", K_idx), ("$\Gamma$", G_idx), ("M", M_idx)]:
        # get cumulative for this idx (interpolate if needed)
        if idx in ik_to_cum:
            cval = ik_to_cum[idx]
        else:
            pos = np.searchsorted(map_idx_arr, idx)
            if pos == 0:
                cval = map_cum_arr[0]
            elif pos >= len(map_idx_arr):
                cval = map_cum_arr[-1]
            else:
                i0 = pos - 1; i1 = pos
                x0 = map_idx_arr[i0]; x1 = map_idx_arr[i1]
                y0 = map_cum_arr[i0]; y1 = map_cum_arr[i1]
                if x1 == x0:
                    cval = y0
                else:
                    frac = (idx - x0) / (x1 - x0)
                    cval = y0 + frac * (y1 - y0)
        norm_pos = (cval - map_cum_arr[0]) / (map_cum_arr[-1] - map_cum_arr[0])
        tick_positions.append((lab, norm_pos))

    return s_norm, tick_positions


# ------------------------------------------------------------
# Read mapping output (SHIFT energies by E_F)
# ------------------------------------------------------------
def read_population_file(path):
    df = []
    modes = []

    with open(path, "r") as fh:
        for line in fh:
            if line.startswith("# Dominant modes:"):
                modes = eval(line.split(":")[1].strip())

            if line.lstrip().startswith("#") or len(line.strip()) == 0:
                continue

            parts = line.split()
            ik = int(parts[0]) - 1
            Enk = float(parts[1]) - E_fermi  # SHIFT BY FERMI

            row = {"ik": ik, "energy": Enk}

            offset = 2
            for m in modes:
                row[f"abs_m{m}"] = float(parts[offset])
                row[f"em_m{m}"] = float(parts[offset + 1])
                row[f"iqA_m{m}"] = int(parts[offset + 2])
                row[f"iqE_m{m}"] = int(parts[offset + 3])
                offset += 4

            df.append(row)

    return pd.DataFrame(df), modes


# Bubble plot (all modes, clear styling)
# ------------------------------------------------------------
def plot_bubble(df, modes, outpath, field_value):

    # ---- symmetry transformation ----
    s, pts = build_symmetry_axis(df)
    Enk = df["energy"].values

    # ------------------------------------------------------------
    # Gradient colormap for modes (NOW modes is defined!)
    # ------------------------------------------------------------
    cmap = plt.get_cmap("PiYG")
    m_min, m_max = min(modes), max(modes)

    mode_colors = {}
    for m in modes:
        frac = (m - m_min) / (m_max - m_min + 1e-12)
        mode_colors[m] = cmap(frac)
    # ------------------------------------------------------------

    plt.figure(figsize=(13, 5))

    # ---- dotted band ----
    order = np.argsort(s)
    plt.plot(s[order], Enk[order], "k--", alpha=0.65, zorder=0)

    # ---- symmetry ticks ----
    xticks = [p[1] for p in pts]
    pad = 0.0
    plt.xlim(-pad, 1 + pad)

    xlabels = [p[0] for p in pts]
    plt.xticks(xticks, xlabels)
    for x in xticks:
        plt.axvline(x, color="gray", linestyle="--", linewidth=1.0, zorder=0)

    #plt.xlim(-0.03, 1.03)

    # ---- offsets to avoid overlap ----
    offset_step = 1e-2

    for i, m in enumerate(modes):

        dx = (i - (len(modes)-1)/2.0) * offset_step
        x_shifted = np.clip(s + dx, -0.03, 1.03)

        absw = df[f"abs_m{m}"].values
        emw  = df[f"em_m{m}"].values

        abs_size = (absw / max(absw.max(), 1e-12)) * 400
        em_size  = (emw  / max(emw.max(), 1e-12)) * 400

        # Abs bubbles
        plt.scatter(
            x_shifted, Enk, s=abs_size, color=mode_colors[m],
            alpha=0.5, edgecolor="none", linewidth=0.8,
            marker="o", label=f"Abs (mode {m+1})"
        )

        # Em bubbles (same colormap but triangle)
        plt.scatter(
            x_shifted, Enk, s=em_size, color=mode_colors[m],
            alpha=0.5, edgecolor="black", linewidth=0.2,
            marker="^", label=f"Em (mode {m+1})"
        )

    # ---- Deduplicate legend ----
    handles, labels = plt.gca().get_legend_handles_labels()
    uniq = {}
    h2, l2 = [], []
    for h, l in zip(handles, labels):
        if l not in uniq:
            uniq[l] = True
            h2.append(h)
            l2.append(l)

    plt.legend(
    h2, l2,
    frameon=True,
    ncol=1,
    fontsize=22,
    framealpha=0.30,          # <-- transparency of legend box
    facecolor="white",        # <-- background color
    edgecolor="black",        # <-- border
    loc="upper right"
    )


    plt.ylabel(r"E - E$_{\mathrm{F}}$ (eV)")
        # ---- Field label (converted to kV/cm, formatted) ----
    fv = field_value * 1e-5  # V/m → kV/cm
    s = f"{fv:.6e}"          # scientific notation
    mant_str, exp_str = s.split("e")
    mant = float(mant_str)
    exp = int(exp_str)
    sci_str = rf"{mant:.3f}×10$^{{{exp}}}$ kV/cm"

    plt.text(
        0.02, 0.92,
        f"Field = {sci_str}",
        transform=plt.gca().transAxes,
        fontsize=26,
        verticalalignment='top',
        bbox=dict(
            facecolor="white",
            edgecolor="black",
            alpha=0.35,      # transparency of the box
            pad=6            # padding inside box
        )
    )

    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()

    print("[OK] Bubble plot →", outpath)



# ------------------------------------------------------------
# Histogram plot (now using secondary y-axis)
# ------------------------------------------------------------
def plot_hist(df, modes, outpath, field_value, global_max):

    # apply symmetry mapping
    s, pts = build_symmetry_axis(df)
    ik = s

    order = np.argsort(ik)

    fig, ax1 = plt.subplots(figsize=(13, 5))

    width = 0.02  # small width in normalized axis units

    # ---- Primary axis: energy band ----
    ax1.plot(
        ik[order],
        df["energy"].values[order],
        "k--",
        alpha=0.7
    )
    ax1.set_ylabel(r"E - E$_{\mathrm{F}}$ (eV)")

    # ---- symmetry ticks ----
    xticks = [p[1] for p in pts]
    xlabels = [p[0] for p in pts]
    pad = 0.03
    ax1.set_xlim(-pad, 1 + pad)
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xlabels)

    for x in xticks:
        ax1.axvline(x, color="gray", linestyle="--", linewidth=1.0)

    # ------------------------------------------------------------
    # SECONDARY AXIS for histogram data
    # ------------------------------------------------------------
    ax2 = ax1.twinx()
    ax2.set_ylabel("Normalized Count")

    # Draw histogram bars on secondary axis
    for m in modes:
        ax2.bar(
            ik - width/2,
            df[f"abs_m{m}"].values / global_max,
            width=width,
            alpha=0.75,
            label=f"Abs (mode {m+1})"
        )
        ax2.bar(
            ik + width/2,
            df[f"em_m{m}"].values / global_max,
            width=width,
            alpha=0.75,
            label=f"Em (mode {m+1})"
        )

    # ------------------------------------------------------------
    # Deduplicate legend entries (for both axes)
    # ------------------------------------------------------------
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()

    handles = handles1 + handles2
    labels = labels1 + labels2

    uniq = {}
    h2, l2 = [], []
    for h, l in zip(handles, labels):
        if l not in uniq:
            uniq[l] = True
            h2.append(h)
            l2.append(l)
    ax2.set_ylim(0, 1.05)
    ax2.legend(
        h2, l2,
        frameon=True,
        ncol=1,
        fontsize=22,
        framealpha=0.30,
        facecolor="white",
        edgecolor="black",
        loc="upper right"
    )
        # ---- Field label (converted to kV/cm, formatted) ----
    fv = field_value * 1e-5  # V/m → kV/cm
    s = f"{fv:.6e}"          # scientific notation
    mant_str, exp_str = s.split("e")
    mant = float(mant_str)
    exp = int(exp_str)
    sci_str = rf"{mant:.3f}×10$^{{{exp}}}$ kV/cm"
    plt.text(
        0.02, 0.92,
        f"Field = {sci_str}",
        transform=plt.gca().transAxes,
        fontsize=26,
        verticalalignment='top',
        bbox=dict(
            facecolor="white",
            edgecolor="black",
            alpha=0.35,      # transparency of the box
            pad=6            # padding inside box
        )
    )


    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)

    print("[OK] Histogram →", outpath)


import imageio
from matplotlib.gridspec import GridSpec


# ============================================================
# Combined TOP (bubble) + BOTTOM (histogram) figure per field
# ============================================================
def plot_combined(df, modes, outpath, field_value):

    # ---- Build symmetry axis once ----
    s, pts = build_symmetry_axis(df)
    Enk = df["energy"].values

    # ---- Colormap for modes (same as bubble) ----
    cmap = plt.get_cmap("viridis")
    m_min, m_max = min(modes), max(modes)
    mode_colors = {m: cmap((m - m_min) / (m_max - m_min + 1e-12)) for m in modes}

    # ============================================================
    # Create figure with two rows (top: bubble, bottom: histogram)
    # ============================================================
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 1, height_ratios=[3, 2], hspace=0.25)

    # ============================================================
    # ---------------------- TOP PANEL ---------------------------
    # ============================================================
    ax1 = fig.add_subplot(gs[0])

    # band curve
    order = np.argsort(s)
    ax1.plot(s[order], Enk[order], "k--", alpha=0.65, zorder=0)

    # symmetry ticks
    xticks = [p[1] for p in pts]
    ax1.set_xlim(-0.03, 1.03)
    ax1.set_xticks(xticks)
    ax1.set_xticklabels([p[0] for p in pts])
    for x in xticks:
        ax1.axvline(x, color="gray", linestyle="--", linewidth=1.0, zorder=0)

    # plot bubbles
    for i, m in enumerate(modes):
        absw = df[f"abs_m{m}"].values
        emw = df[f"em_m{m}"].values

        abs_size = (absw / max(absw.max(), 1e-12)) * 450
        em_size  = (emw  / max(emw.max(), 1e-12)) * 450

        ax1.scatter(s, Enk, s=abs_size, color=mode_colors[m],
                    alpha=0.55, edgecolor="black", marker="o",
                    label=f"Abs (mode {m+1})")

        ax1.scatter(s, Enk, s=em_size, color=mode_colors[m],
                    alpha=0.55, edgecolor="black", marker="^",
                    label=f"Em (mode {m+1})")

    # deduplicate legend
    handles, labels = ax1.get_legend_handles_labels()
    uniq = {}
    h2, l2 = [], []
    for h, l in zip(handles, labels):
        if l not in uniq:
            uniq[l] = True
            h2.append(h)
            l2.append(l)

    ax1.legend(h2, l2, frameon=True, framealpha=0.30,
               facecolor="white", edgecolor="black",
               fontsize=20, ncol=2, loc="upper right")

    # field label inside box
    sci = f"{field_value/1e5:.3g}"  # convert to kV/cm
    ax1.text(
        0.02, 0.93,
        f"Field = {sci} kV/cm",
        transform=ax1.transAxes,
        fontsize=24,
        bbox=dict(facecolor="white", edgecolor="black", alpha=0.3),
    )

    ax1.set_ylabel(r"$E - E_F$ (eV)")

    # ============================================================
    # --------------------- BOTTOM PANEL -------------------------
    # ============================================================
    ax2 = fig.add_subplot(gs[1])
    ax3 = ax2.twinx()  # secondary axis for histogram

    # energy line again (on left axis)
    ax2.plot(s[order], Enk[order], "k--", alpha=0.6)

    # symmetry ticks
    ax2.set_xlim(-0.03, 1.03)
    ax2.set_xticks(xticks)
    ax2.set_xticklabels([p[0] for p in pts])
    for x in xticks:
        ax2.axvline(x, color="gray", linestyle="--", linewidth=1.0)

    # histograms on ax3 (right axis)
    width = 0.02
    for m in modes:
        ax3.bar(s - width/2, df[f"abs_m{m}"].values / global_max,
                width=width, alpha=0.6, label=f"Abs m{m+1}")
        ax3.bar(s - width/2, df[f"em_m{m}"].values / global_max,
                width=width, alpha=0.6, label=f"Em m{m+1}")

    ax2.set_ylabel(r"$E - E_F$ (eV)")
    ax3.set_ylabel("Normalized Count")
    ax3.set_ylim(0, 1.05)

    # legend (right axis)
    handles, labels = ax3.get_legend_handles_labels()
    ax3.legend(handles, labels, frameon=True, framealpha=0.30,
               facecolor="white", edgecolor="black",
               fontsize=18, ncol=1, loc="upper right")

    # SAVE FIGURE
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)

    print("[OK] Combined plot →", outpath)



# ============================================================
# Make GIF from all combined figures
# ============================================================
def make_gif(png_folder="plots_population", outname="combined.gif", duration=2000):
    frames = []
    for f in sorted(os.listdir(png_folder)):
        if f.endswith("_combined.png"):
            frames.append(imageio.imread(os.path.join(png_folder, f)))

    if len(frames) == 0:
        print("[WARN] No combined PNG files found.")
        return

    imageio.mimsave(outname, frames, duration=duration)
    print("[OK] GIF saved →", outname)

from PIL import Image

def combine_existing_pngs(bubble_path, hist_path, outpath):

    if not os.path.isfile(bubble_path):
        print(f"[ERR] Missing bubble PNG: {bubble_path}")
        return

    if not os.path.isfile(hist_path):
        print(f"[ERR] Missing hist PNG: {hist_path}")
        return

    img1 = Image.open(bubble_path)
    img2 = Image.open(hist_path)

    # Create stacked canvas
    w = max(img1.width, img2.width)
    h = img1.height + img2.height

    combined = Image.new("RGB", (w, h), (255, 255, 255))
    combined.paste(img1, (0, 0))
    combined.paste(img2, (0, img1.height))

    combined.save(outpath)
    print("[OK] Combined PNG →", outpath)

import re
from PIL import Image

# Extract numeric field value from filename like "E_1.000e+03_combined.png"
def extract_field(fname):
    m = re.search(r"E_([0-9.eE+-]+)", fname)
    if not m:
        return None
    try:
        return float(m.group(1))
    except:
        return None


def make_gif_from_combined(outdir="plots_population", outname="combined_fields.gif", duration=3000):

    print("\n[INFO] Creating GIF...")

    # collect all combined PNGs
    pngs = []
    for f in os.listdir(outdir):
        if f.endswith("_combined.png"):
            fval = extract_field(f)
            if fval is not None:
                pngs.append((fval, os.path.join(outdir, f)))
            else:
                print(f"[WARN] Cannot extract field from: {f}")

    if not pngs:
        print("[ERROR] No combined PNGs found!")
        return

    # sort by field
    pngs_sorted = sorted(pngs, key=lambda x: x[0])
    files_sorted = [p for (_, p) in pngs_sorted]

    frames = []
    for p in files_sorted:
        try:
            im = Image.open(p).convert("RGB")
            frames.append(im)
        except Exception as e:
            print(f"[ERR] Could not read {p}: {e}")

    if not frames:
        print("[ERROR] No frames loaded!")
        return

    outpath = os.path.join(outdir, outname)

    frames[0].save(
        outpath,
        save_all=True,
        append_images=frames[1:],
        duration=duration,
        loop=0
    )

    print(f"[OK] GIF saved → {outpath}")


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------
def main():

    input_dir = "population_out"
    output_dir = "plots_population"
    os.makedirs(output_dir, exist_ok=True)

    files = sorted([f for f in os.listdir(input_dir)
                    if f.startswith("E_") and f.endswith(".txt")])

    if not files:
        print("[ERROR] No mapping files found.")
        return

    # -------- GLOBAL NORMALIZATION --------
    global_max = 0.0
    for fname in files:
        df_tmp, modes_tmp = read_population_file(os.path.join(input_dir, fname))
        for m in modes_tmp:
            global_max = max(global_max,
                             df_tmp[f"abs_m{m}"].max(),
                             df_tmp[f"em_m{m}"].max())
    if global_max == 0:
        global_max = 1.0


    combined_pngs = []

    for fname in files:
        path = os.path.join(input_dir, fname)

        m = re.search(r"E_([0-9.eE+-]+)\.txt", fname)
        field_value = float(m.group(1)) if m else 0.0

        df, modes = read_population_file(path)
        base = fname.replace(".txt", "")

        # Output file names
        bubble_out   = os.path.join(output_dir, base + "_bubble.png")
        hist_out     = os.path.join(output_dir, base + "_hist.png")
        combined_out = os.path.join(output_dir, base + "_combined.png")

        # Generate two individual plots
        plot_bubble(df, modes, bubble_out, field_value)
        plot_hist(df, modes, hist_out, field_value, global_max)


        # Combine them (stack bubble above histogram)
        combine_existing_pngs(bubble_out, hist_out, combined_out)
        combined_pngs.append(combined_out)

    # ------------------------------------------------------------
    # MAKE GIF FROM ALL COMBINED PNGs (SORTED BY FIELD)
    # ------------------------------------------------------------
    if combined_pngs:
        print("[INFO] Creating GIF from combined PNGs...")

        # ---- SORT BY NUMERIC ELECTRIC FIELD ----
        combined_pngs = sorted(
            combined_pngs,
            key=lambda p: extract_field(os.path.basename(p))
        )

        frames = [imageio.imread(p) for p in combined_pngs]
        #imageio.mimsave("combined.gif", frames, duration=0.8)
        imageio.mimsave(
            "combined.gif",
            frames,
            duration=600,
            loop=0   # 0 = infinite repeat
        )


        print("[OK] GIF saved → combined.gif")
    else:
        print("[WARN] No combined PNGs were created.")


# ------------------------------------------------------------
# Run main
# ------------------------------------------------------------
if __name__ == "__main__":
    main()

