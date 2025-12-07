import re
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--nscf", required=True, help="NSCF output file")
parser.add_argument("--num_wann", type=int, required=True, help="Target number of Wannier functions")
args = parser.parse_args()

with open(args.nscf) as f:
    lines = f.readlines()

bands = []
count, in_occ = 0, False
eigs, all_eigs = [], []

for line in lines:
    if re.match(r"^\s*k\s*=", line):
        if count:
            bands.append(count)
            all_eigs.extend(eigs)
        count, in_occ = 0, False
        eigs = []
    elif "occupation numbers" in line:
        in_occ = True
    elif line.strip() == "":
        in_occ = False
    elif not in_occ:
        try:
            vals = [float(x) for x in line.split()]
            eigs.extend(vals)
            count += len(vals)
        except ValueError:
            pass

if count:
    bands.append(count)
    all_eigs.extend(eigs)

# Per-k stats
print(f"Parsed {len(bands)} k-points; bands per k point: min={min(bands)}, max={max(bands)}, mean={np.mean(bands):.2f}")

# Global band energies
emin = min(all_eigs)
emax = max(all_eigs)
print(f"Global energy range: {emin:.3f} to {emax:.3f} eV")

# Find VBM and CBM
vbm, cbm = None, None
for line in lines:
    if "highest occupied" in line and "lowest unoccupied" in line:
        parts = line.split()
        vbm, cbm = float(parts[-2]), float(parts[-1])
        break

if vbm is None or cbm is None:
    print("⚠️ Could not detect VBM/CBM automatically")
else:
    print(f"Detected VBM, CBM from nscf: {vbm:.4f}, {cbm:.4f}")

    # Suggest disentanglement windows
    dis_win_min = emin
    dis_win_max = emax
    dis_froz_min = vbm - 0.5   # 0.5 eV below VBM
    dis_froz_max = cbm + 0.5   # 0.5 eV above CBM

    print("\nSuggested windows for wannierization:")
    print(f"  dis_win_min  = {dis_win_min:.2f}")
    print(f"  dis_win_max  = {dis_win_max:.2f}")
    print(f"  dis_froz_min = {dis_froz_min:.2f}")
    print(f"  dis_froz_max = {dis_froz_max:.2f}")

    print(f"\n⚡ These are starting points. Adjust if 'More states than target WFs' appears.")

