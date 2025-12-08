#!/usr/bin/env python3
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image

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

    print(f"[INFO] Auto-detected symmetry points: K={K_idx}, Γ={G_idx}, M={M_idx}")
    return K_idx, G_idx, M_idx


# Global max Y (will be computed in main)
global_maxY = 1

# --------------------------------------------------------
# Read kcloud CBn_abs/em files
# --------------------------------------------------------
def read_kcloud_file(path):
    data = []
    with open(path, "r") as f:
        for line in f:
            if line.strip().startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                ik = int(parts[0])
            except:
                continue
            try:
                energy = float(parts[1])
            except:
                energy = 0.0
            try:
                cnt = int(parts[2])
            except:
                cnt = 0
            data.append([ik, energy, cnt])
    if len(data) == 0:
        return pd.DataFrame(columns=["ik", "energy_eV", "count"])
    return pd.DataFrame(data, columns=["ik", "energy_eV", "count"])


# --------------------------------------------------------
# Classification rules
# --------------------------------------------------------
def classify(abs_count, em_count, z):
    if abs_count == 0 and em_count == 0:
        return "none"
    if z > 3:
        return "abs_dominated"
    if z < -3:
        return "em_dominated"
    return "mixed"

# --------------------------------------------------------
# Process single field for one CBn
# --------------------------------------------------------
def process_pair(abs_file, em_file, outdir_global, field_name):
    global global_maxY

    df_abs = read_kcloud_file(abs_file).rename(columns={"count": "abs_count"})
    df_em = read_kcloud_file(em_file).rename(columns={"count": "em_count"})

    df = pd.merge(df_abs, df_em, on="ik", how="outer", suffixes=("_abs", "_em")).fillna(0)

    df["energy_eV"] = df.apply(
        lambda r: r.energy_eV_abs if (r.energy_eV_abs != 0 and not pd.isna(r.energy_eV_abs)) else r.energy_eV_em,
        axis=1,
    )

    df["abs_count"] = df["abs_count"].fillna(0).astype(int)
    df["em_count"] = df["em_count"].fillna(0).astype(int)

    df["total"] = df["abs_count"] + df["em_count"]
    df["net"] = df["abs_count"] - df["em_count"]
    df["diff_z"] = df["net"] / np.sqrt(df["total"] + 1)

    def safe_ratio(a, b):
        if b == 0:
            return np.inf if a > 0 else np.nan
        return a / b

    df["ratio_abs_em"] = df.apply(lambda r: safe_ratio(r.abs_count, r.em_count), axis=1)
    df["class"] = df.apply(lambda r: classify(r.abs_count, r.em_count, r.diff_z), axis=1)

    cbname = os.path.basename(abs_file).split("_")[0]
    prefix = f"{field_name}_{cbname}"

    K_idx, G_idx, M_idx = detect_K_G_M_indices("report_parsed.txt")

    # -------------------------
    # Plotting
    # -------------------------
    plt.figure(figsize=(12, 5))
    width = 0.4
    ik = df["ik"].values

    order = np.argsort(ik)
    ik_plot = ik[order]

    abs_raw = df["abs_count"].values[order]
    em_raw  = df["em_count"].values[order]

    # ----------------------------------------
    # NORMALIZATION (per field)
    # ----------------------------------------
    max_val = max(abs_raw.max(), em_raw.max(), 1)
    abs_plot = abs_raw / max_val
    em_plot  = em_raw  / max_val

    # ----------------------------------------
    # BAR PLOT (linear normalized)
    # ----------------------------------------
    plt.bar(ik_plot - width / 2, abs_plot, width=width, color="orange", alpha=0.75, label="Absorption")
    plt.bar(ik_plot + width / 2, em_plot, width=width, color="green",  alpha=0.75, label="Emission")

    plt.xticks([K_idx, G_idx, M_idx], ["K", r"$\Gamma$", "M"])
    plt.axvline(G_idx, color="gray", linestyle="--", linewidth=1.2, alpha=0.8)

    plt.ylabel("Scattering Count")
    plt.margins(x=0.005)

    # ---------------------------
    # LEFT-TOP ELECTRIC FIELD TEXT
    # ---------------------------
    m = re.search(r"E_([0-9.eE+-]+)", field_name)
    if m:
        fraw = float(m.group(1))             # V/m
        fv = fraw * 1e-5                    # kV/cm
        s = f"{fv:.6e}"
        mant, exp = s.split("e")
        mant = float(mant)
        mant_fmt = f"{mant:.3f}"
        exp = int(exp)
        sci_str = fr"Field = {mant_fmt}×10$^{{{exp}}}$ kV/cm"
    else:
        sci_str = "E = ?"

    plt.text(
        0.02, 0.95, sci_str,
        transform=plt.gca().transAxes,
        fontsize=28,
        verticalalignment='top',
        #bbox=dict(facecolor='white', edgecolor='black', alpha=0.80)
    )

    # ---------------------------
    # RIGHT-TOP LEGEND (Abs/Em)
    # ---------------------------
    plt.legend(
        loc="upper right",
        frameon=False,
        facecolor="white",
        edgecolor="black",
        framealpha=0.85,
        fontsize=28
    )

    plt.tight_layout()

    plot_path = os.path.join(outdir_global, f"{prefix}_abs_em_combined.png")
    plt.savefig(plot_path)
    plt.close()

    print(f"[OK] Processed {prefix}")
    return plot_path

    # ---------------------------
    # Structured summary
    # ---------------------------
    txt_path = os.path.join(outdir_global, f"{prefix}_summary.txt")

    with open(txt_path, "w") as f:
        f.write(f"# Summary for {cbname} at field {field_name}\n")
        f.write("# ------------------------------------------------------------------------------------------------------------------------------\n")
        f.write("# ik     energy_eV    abs_cnt  em_cnt  total   net   diff_z   ratio_abs_em   class\n")
        f.write("# ------------------------------------------------------------------------------------------------------------------------------\n")

        for _, r in df.iterrows():
            ratio = r["ratio_abs_em"]
            if not np.isfinite(ratio):
                ratio = float("inf")

            f.write(
                f"{int(r['ik']):4d}   "
                f"{r['energy_eV']:12.6f}   "
                f"{int(r['abs_count']):8d}   "
                f"{int(r['em_count']):8d}   "
                f"{int(r['total']):8d}   "
                f"{int(r['net']):8d}   "
                f"{r['diff_z']:10.4f}   "
                f"{ratio:14.4f}   "
                f"{r['class']}\n"
            )

    print(f"[OK] Processed {prefix}")
    return plot_path

# --------------------------------------------------------
# Extract numeric field value from filename
# --------------------------------------------------------
def extract_field(filename):
    # matches 'E_1.000e+05_CB...' and returns 1.000e+05 as float
    m = re.search(r"E_([0-9.+-eE]+)_CB", filename)
    if m:
        try:
            return float(m.group(1))
        except:
            return None
    return None

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


# --------------------------------------------------------
# MAIN
# --------------------------------------------------------
def main():
    global global_maxY

    base = "./kcloud"
    outdir_global = "./plots"
    os.makedirs(outdir_global, exist_ok=True)

    if not os.path.isdir(base):
        print("[ERROR] kcloud/ not found in current directory.")
        return

    print("[INFO] Scanning kcloud/...")

    # -------------------------
    # First pass: determine global max count across all abs/em files
    # -------------------------
    gm = 0
    for fld in sorted(os.listdir(base)):
        fld_path = os.path.join(base, fld)
        if not os.path.isdir(fld_path):
            continue
        for fname in os.listdir(fld_path):
            if fname.endswith(".txt") and ("_abs" in fname or "_em" in fname):
                p = os.path.join(fld_path, fname)
                try:
                    df_temp = read_kcloud_file(p)
                    if not df_temp.empty:
                        cmax = int(df_temp["count"].max())
                        if cmax > gm:
                            gm = cmax
                except Exception:
                    pass

    if gm <= 0:
        gm = 1
    global_maxY = gm

    print(f"[INFO] global_maxY = {global_maxY} (used for y-axis limits)")

    all_generated_pngs = []

    # Loop and process
    for fld in sorted(os.listdir(base)):
        fld_path = os.path.join(base, fld)
        if not os.path.isdir(fld_path):
            continue

        field_name = fld
        cb_files = {}

        for fname in os.listdir(fld_path):
            if fname.endswith(".txt"):
                if "_abs" in fname:
                    cb = fname.split("_abs")[0]
                    cb_files.setdefault(cb, {})["abs"] = os.path.join(fld_path, fname)
                if "_em" in fname:
                    cb = fname.split("_em")[0]
                    cb_files.setdefault(cb, {})["em"] = os.path.join(fld_path, fname)

        for cb, files in cb_files.items():
            if "abs" in files and "em" in files:
                try:
                    png_path = process_pair(files["abs"], files["em"], outdir_global, field_name)
                    all_generated_pngs.append(png_path)
                except Exception as e:
                    print(f"[ERROR] processing {cb} in {fld_path}: {e}")
            else:
                print(f"[WARN] Missing pair for {cb} inside {fld_path}")

    # ----------------------------------------------------
    # MAKE GIF
    # ----------------------------------------------------
    print("\n[INFO] Creating GIF...")

    valid_pngs = []
    for p in all_generated_pngs:
        fname = os.path.basename(p)
        fval = extract_field(fname)
        if fval is not None:
            valid_pngs.append((fval, p))
        else:
            print(f"[WARN] Could not extract field value from: {fname}")

    valid_pngs_sorted = sorted(valid_pngs, key=lambda x: x[0])
    all_generated_pngs_sorted = [p for _, p in valid_pngs_sorted]

    if len(all_generated_pngs_sorted) == 0:
        print("[ERROR] No valid images available for GIF.")
        return

    # Load frames; ensure consistent image mode
    frames = []
    for p in all_generated_pngs_sorted:
        try:
            im = Image.open(p).convert("RGBA")
            frames.append(im)
        except Exception as e:
            print(f"[WARN] Could not open image {p}: {e}")

    if len(frames) == 0:
        print("[ERROR] No images could be opened for GIF.")
        return

    gif_path = os.path.join(outdir_global, "combined_fields.gif")
    frames[0].save(
        gif_path,
        save_all=True,
        append_images=frames[1:],
        duration=600,
        loop=0,
        optimize=False,
    )

    print(f"[OK] GIF saved to {gif_path}")


if __name__ == "__main__":
    main()

