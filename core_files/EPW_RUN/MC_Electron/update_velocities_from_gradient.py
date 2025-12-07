#!/usr/bin/env python3
"""
update_velocities_from_gradient.py (fixed)

Parses celldm(1) and the 'reciprocal axes' block from scf.out (QE-style),
computes the reciprocal lattice in SI units (1/m) using the correct unit:
    printed_reciprocal * (2*pi / alat_m)

Then computes event-wise velocities and writes them in-place to scattering_events.h5:
 - v_group_cart, v_group_mag, v_scattered_cart, v_scattered_mag

Also stores B_reciprocal and B_inverse in the HDF5 inside group 'scattering_events'.
Writes an aligned text table event_velocities.txt.

Author: ChatGPT (fixed)
"""

import numpy as np
import h5py
import os
import re
import math

hbar = 1.054571817e-34        # J*s
meV_to_J = 1.602176634e-22    # meV -> J
a0 = 0.529177210903e-10       # Bohr in meters

h5path = "scattering_events.h5"
scf_files = ["scf.out", "scf.ot", "pw.out"]

def find_file(paths):
    for p in paths:
        if os.path.exists(p):
            return p
    return None

def parse_alat_from_scf(txt):
    """Try to parse celldm(1) or lattice parameter (alat) from scf text.
       Return alat in meters.
    """
    # Try celldm(1)
    m = re.search(r"celldm\s*\(\s*1\s*\)\s*=\s*([0-9Ee+.-]+)", txt)
    if m:
        val = float(m.group(1))
        # celldm(1) printed by QE is in atomic units (Bohr)
        alat_bohr = val
        alat_m = alat_bohr * a0
        return alat_m, "bohr"

    # Try lattice parameter (alat) = X a.u.
    m2 = re.search(r"lattice parameter\s*\(alat\)\s*=\s*([0-9Ee+.-]+)\s*(\w*)", txt)
    if m2:
        val = float(m2.group(1))
        unit = m2.group(2).strip().lower()
        if unit.startswith("a.u") or unit.startswith("bohr") or unit == "":
            alat_bohr = val
            alat_m = alat_bohr * a0
            return alat_m, "bohr"
        # if it's in Angstrom (rare), convert directly
        if unit.startswith("ang") or unit.startswith("a"):
            alat_m = val * 1e-10
            return alat_m, "angstrom"

    # attempt to infer from CELL_PARAMETERS printed in units 'alat'
    m3 = re.search(r"CELL_PARAMETERS\s*\(.*\)\s*\n((?:\s*[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s*\n){3})", txt)
    if m3:
        block = m3.group(1)
        rows = []
        for L in block.strip().splitlines():
            nums = re.findall(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?", L)
            if len(nums) >= 3:
                rows.append([float(nums[0]), float(nums[1]), float(nums[2])])
        if len(rows) >= 1:
            # Norm of first vector times alat (unknown) -> cannot infer alat without further info.
            pass

    return None, None

def parse_reciprocal_block(txt):
    """Extract printed reciprocal vectors (the numeric rows) near 'reciprocal axes'.
       Returns numpy array shape (3,3) or None.
    """
    # find 'reciprocal axes' block
    idx = txt.lower().find("reciprocal axes")
    if idx < 0:
        idx = txt.lower().find("reciprocal")
    if idx < 0:
        return None, None

    window = txt[idx: idx + 1000]  # take a chunk
    # Try to capture explicit lines b(1) = ( ... )
    rows = []
    for line in window.splitlines():
        m = re.search(r"b\(\s*\d+\s*\)\s*=\s*\(\s*([^\)]+)\)", line, flags=re.IGNORECASE)
        if m:
            nums = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?", m.group(1))]
            if len(nums) >= 3:
                rows.append(nums[:3])
    # fallback: capture numeric lines in the window
    if len(rows) < 3:
        rows = []
        for line in window.splitlines():
            vals = re.findall(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?", line)
            if len(vals) >= 3:
                try:
                    rows.append([float(vals[0]), float(vals[1]), float(vals[2])])
                except:
                    pass
            if len(rows) >= 3:
                break

    if len(rows) >= 3:
        return np.array(rows[:3], dtype=float), window.lower()
    return None, window.lower()

# -------------------------
# Main: parse scf.out and compute B
# -------------------------
scf_file = find_file(scf_files)
if scf_file is None:
    raise SystemExit("No SCF file found (scf.out / scf.ot / pw.out). Place file in cwd.")

txt = open(scf_file, "r", errors="ignore").read()
alat_m, alat_unit = parse_alat_from_scf(txt)
b_print, window_lower = parse_reciprocal_block(txt)

if b_print is None:
    raise SystemExit("Could not find reciprocal axes in SCF file.")

# Decide how to interpret b_print
# If the window says 'in units 2 pi/alat' we must use scale = 2*pi / alat_m
if "2 pi/alat" in window_lower or "2pi/alat" in window_lower or "2 pi / alat" in window_lower:
    if alat_m is None:
        raise SystemExit("Reciprocal axes are in units 2 pi/alat but 'alat' not found in scf.out.")
    scale = 2.0 * math.pi / alat_m
    B = b_print * scale
    units_expl = "2*pi/alat -> converted using celldm(1)"
else:
    # If the magnitudes look like ~1 (order 1) and we have celldm(1), prefer 2pi/alat as above.
    # Else, if numbers look like 1/Angstrom, detect and convert.
    meanabs = np.mean(np.abs(b_print))
    if alat_m is not None and meanabs < 100.0:
        # prefer 2pi/alat since QE often prints that way even without explicit string
        scale = 2.0 * math.pi / alat_m
        B = b_print * scale
        units_expl = "assumed 2*pi/alat (based on celldm(1))"
    else:
        # fallback: assume numbers are printed in 1/Angstrom
        B = b_print * 1e10
        units_expl = "assumed 1/Angstrom"

# Compute inverse
B_inv = np.linalg.inv(B)

print("[INFO] scf file:", scf_file)
print("[INFO] alat (m):", alat_m, "unit:", alat_unit)
print("[INFO] Interpretation:", units_expl)
print("[INFO] B_reciprocal (1/m):\n", B)
print("[INFO] B_inverse (m):\n", B_inv)

# ----------------------------------------
# Load scattering events HDF5 and update
# ----------------------------------------
with h5py.File(h5path, "r+") as h5:
    if "scattering_events" not in h5:
        raise SystemExit("scattering_events group not found in HDF5.")
    grp = h5["scattering_events"]

    # store or replace reciprocal matrices
    if "B_reciprocal" in grp:
        del grp["B_reciprocal"]
    if "B_inverse" in grp:
        del grp["B_inverse"]

    grp.create_dataset("B_reciprocal", data=B)
    grp.create_dataset("B_inverse", data=B_inv)
    print("[OK] Stored B_reciprocal and B_inverse in HDF5.")

    # remove obsolete datasets if present (preserve others)
    for old in ["full_k_frac", "full_velocities_m_s"]:
        if old in grp:
            del grp[old]

    # load required arrays
    k   = grp["k_frac"][...]       # (N,3)
    kq  = grp["kq_frac"][...]      # (N,3)
    E   = grp["enk_meV"][...]      # (N,)
    Eq  = grp["enkq_meV"][...]     # (N,)

    N = k.shape[0]

    v_init   = np.zeros((N,3))
    v_init_M = np.zeros(N)
    v_scat   = np.zeros((N,3))
    v_scat_M = np.zeros(N)

    tiny = 1e-50

    # compute event-wise velocities using correct B (fractional -> cart: k_cart = B.T @ k_frac)
    for i in range(N):
        dk_frac = kq[i] - k[i]
        dk_cart = B.T @ dk_frac     # convert fractional Δk -> cartesian (1/m)

        norm2 = np.dot(dk_cart, dk_cart)

        if norm2 <= tiny:
            v = np.zeros(3)
        else:
            dE = (Eq[i] - E[i]) * meV_to_J   # Joules
            gradk = (dE / norm2) * dk_cart   # J/m
            v = gradk / hbar                 # m/s

        v_init[i]   = v
        v_init_M[i] = np.linalg.norm(v)
        v_scat[i]   = v.copy()
        v_scat_M[i] = v_init_M[i]

    # write velocity datasets (replace if exist)
    for name in ["v_group_cart", "v_group_mag", "v_scattered_cart", "v_scattered_mag"]:
        if name in grp:
            del grp[name]
    grp.create_dataset("v_group_cart",      data=v_init)
    grp.create_dataset("v_group_mag",       data=v_init_M)
    grp.create_dataset("v_scattered_cart",  data=v_scat)
    grp.create_dataset("v_scattered_mag",   data=v_scat_M)

    print(f"[OK] Updated velocities for {N} events in-place.")

# Write TXT output
txt_path = "event_velocities.txt"
with open(txt_path, "w") as ftxt:
    ftxt.write("# Event-wise velocities (initial and scattered)\n")
    ftxt.write("# Units: velocities in m/s\n")
    header = (
        f"{'#event':>6}  "
        f"{'vx_init':>15} {'vy_init':>15} {'vz_init':>15} {'vmag_init':>15}  "
        f"{'vx_scatt':>15} {'vy_scatt':>15} {'vz_scatt':>15} {'vmag_scatt':>15}\n"
        "#--------------------------------------------------------------------------------\n"
    )
    ftxt.write(header)

    with h5py.File(h5path, "r") as h5:
        grp = h5["scattering_events"]
        v_init = grp["v_group_cart"][...]
        v_init_M = grp["v_group_mag"][...]
        v_scat = grp["v_scattered_cart"][...]
        v_scat_M = grp["v_scattered_mag"][...]

    for i in range(N):
        ftxt.write(
            f"{i:6d}  "
            f"{v_init[i,0]:15.6e} {v_init[i,1]:15.6e} {v_init[i,2]:15.6e} {v_init_M[i]:15.6e}  "
            f"{v_scat[i,0]:15.6e} {v_scat[i,1]:15.6e} {v_scat[i,2]:15.6e} {v_scat_M[i]:15.6e}\n"
        )

print(f"[OK] TXT velocity table written to {txt_path}")

