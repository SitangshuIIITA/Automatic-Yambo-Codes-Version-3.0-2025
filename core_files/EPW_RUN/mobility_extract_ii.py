#!/usr/bin/env python3
import os
import re
import numpy as np

# ------------------------------------------------------------
# Parse a single mobility output file (quadrupole type)
# ------------------------------------------------------------
def parse_mobility_file(path):
    """Return (mu_SERTA, mu_BTE, rH_SERTA, rH_BTE) from one *_quad.out file"""
    if not os.path.exists(path):
        print(f"⚠️ File not found: {path}")
        return None

    with open(path, "r") as f:
        data = f.read()

    # Split into SERTA and BTE sections
    sections = re.split(r"=+\s*BTE\s*=+", data, flags=re.IGNORECASE)
    serta_block = sections[0]
    bte_block = sections[1] if len(sections) > 1 else ""

    def extract(block):
        mu_vals, mu_hall_vals = [], []
        temp_blocks = re.split(r"Temperature:\s+([\d.]+)\s*K", block)
        for i in range(1, len(temp_blocks), 2):
            content = temp_blocks[i + 1]
            m = re.search(r"Mobility tensor without magnetic field.*?\n(.*)", content)
            if not m:
                continue
            line = m.group(1).strip()
            if "|" not in line:
                continue

            left, right = line.split("|", 1)
            try:
                mu = float(left.strip().split()[0])
            except Exception:
                mu = np.nan
            try:
                hall_vals = re.findall(r"[-+]?\d*\.\d+E[+-]?\d+", right)
                mu_hall = max([abs(float(v)) for v in hall_vals]) if hall_vals else np.nan
            except Exception:
                mu_hall = np.nan

            mu_vals.append(mu)
            mu_hall_vals.append(mu_hall)

        mu_vals, mu_hall_vals = np.array(mu_vals), np.array(mu_hall_vals)
        with np.errstate(divide="ignore", invalid="ignore"):
            rH = np.where(mu_vals != 0, mu_hall_vals / mu_vals, np.nan)
        return np.nanmean(mu_vals), np.nanmean(mu_hall_vals), np.nanmean(rH)

    mu_s, muh_s, rH_s = extract(serta_block)
    mu_b, muh_b, rH_b = extract(bte_block)
    return mu_s, mu_b, rH_s, rH_b


# ------------------------------------------------------------
# Collect all mobility data and write one summary file
# ------------------------------------------------------------
def collect_mobility(base_dir, carrier_type):
    summary_file = f"{carrier_type}_mobility_vs_ii.txt"
    entries = []

    for subdir in sorted(os.listdir(base_dir)):
        subpath = os.path.join(base_dir, subdir)
        if not os.path.isdir(subpath):
            continue

        # Extract concentration from folder name (e.g., 1.0d10)
        conc_match = re.search(r"(\d+\.\d+d\d+)", subdir)
        if not conc_match:
            continue
        conc_str = conc_match.group(1)
        conc_val = float(conc_str.replace("d", "E"))

        # Find file inside folder (e.g., 1.0d10_quad.out)
        quad_file = os.path.join(subpath, f"{conc_str}_quad.out")
        parsed = parse_mobility_file(quad_file)
        if parsed is None:
            continue

        mu_s, mu_b, rHs, rHb = parsed
        entries.append((conc_val, mu_s, mu_b, rHs, rHb))

    # Sort by concentration
    entries.sort(key=lambda x: x[0])

    # Write summary
    with open(summary_file, "w") as f:
        f.write("# Concentration(cm^-2)   mu_SERTA(Quad)   mu_BTE(Quad)   r_H_SERTA(Quad)   r_H_BTE(Quad)\n")
        f.write("=" * 100 + "\n")
        for conc, mu_s, mu_b, rHs, rHb in entries:
            f.write(f"{conc:16.3E}{mu_s:18.4f}{mu_b:18.4f}{rHs:18.6f}{rHb:18.6f}\n")

    print(f"✅ {summary_file} written with {len(entries)} entries.")


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
if __name__ == "__main__":
    ROOT = "./"

    hole_dir = os.path.join(ROOT, "II_H-mobility")
    if os.path.exists(hole_dir):
        collect_mobility(hole_dir, "hole")

    electron_dir = os.path.join(ROOT, "II_E-mobility")
    if os.path.exists(electron_dir):
        collect_mobility(electron_dir, "electron")

    print("✅ All mobility summaries generated successfully.")

