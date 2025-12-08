#!/usr/bin/env python3
"""
epw_parser_and_prep_fixed.py

EPW parser + preprocessor — FIXED VERSION (NO ROTATION)

This script:
 - parses hole_inv_taucb_0.fmt (tau & energies)
 - parses hole_IBTEvel_sup_0.fmt (velocities in a.u., eig, weight)
 - converts velocities: a.u. -> m/s (NO rotation applied — EPW already outputs Cartesian)
 - writes epw_preproc.h5 with datasets expected by the MC drivers

Minimal addition: parses hole_m_effective.fmt (optional) and writes T_Mstar & Mstar_T datasets.
"""

import numpy as np
import h5py
import re
import argparse
import sys
from pathlib import Path

# Constants
RY_TO_EV = 13.605693
AU_TO_MS = 2.18769126364e6     # 1 a.u. velocity -> m/s

# ---------- parsers ----------
def parse_tau_file(fname):
    temps, kpts, bands, energy, invtime = [], [], [], [], []
    with open(fname, "r") as f:
        for line in f:
            if line.strip().startswith("#") or not line.strip():
                continue
            vals = line.split()
            if len(vals) < 5:
                continue
            try:
                itemp = int(vals[0])
                ik = int(vals[1])
                ib = int(vals[2])
                Ery = float(vals[3])
                tRy = float(vals[4])   # this is inverse-time in Ry (EPW format)
            except ValueError:
                continue
            temps.append(itemp)
            kpts.append(ik)
            bands.append(ib)
            energy.append(Ery * RY_TO_EV)   # Ry -> eV
            invtime.append(tRy)

    if len(temps) == 0:
        raise RuntimeError(f"No valid data read from {fname}")

    temps = np.array(temps, dtype=int)
    kpts = np.array(kpts, dtype=int)
    bands = np.array(bands, dtype=int)
    energy = np.array(energy, dtype=float)
    invtime = np.array(invtime, dtype=float)

    ntemp = int(temps.max())
    nband = int(bands.max())
    Nk = int(kpts.max())

    # allocate
    Gamma3d = np.zeros((ntemp, nband, Nk), dtype=float)   # rates (1/ps)
    Tau3d = np.zeros((ntemp, nband, Nk), dtype=float)     # lifetimes (ps)
    E3d = np.zeros((ntemp, nband, Nk), dtype=float)

    # constant: multiply tRy by this to get 1/ps (EPW header)
    CONV_INV = 20670.6944033

    # fill arrays (file is 1-based)
    for it, ik, ib, e, tRy in zip(temps, kpts, bands, energy, invtime):
        i_t = int(it) - 1
        i_b = int(ib) - 1
        i_k = int(ik) - 1
        rate_1_per_ps = tRy * CONV_INV
        Gamma3d[i_t, i_b, i_k] = rate_1_per_ps
        if rate_1_per_ps > 0.0:
            Tau3d[i_t, i_b, i_k] = 1.0 / rate_1_per_ps   # ps
        else:
            Tau3d[i_t, i_b, i_k] = 0.0
        E3d[i_t, i_b, i_k] = e

    return {"tau3d": Tau3d, "Gamma3d": Gamma3d, "E3d": E3d, "dims": (ntemp, nband, Nk)}



def parse_velocity_file(fname):
    """
    Parse IBTEvel file. Velocity components are in a.u. (already Cartesian in EPW).
    We convert them directly to m/s. NO ROTATION APPLIED.
    """
    kpts, bands = [], []
    vx_list, vy_list, vz_list = [], [], []

    with open(fname, "r") as f:
        for line in f:
            if line.strip().startswith("#") or not line.strip():
                continue
            vals = line.split()
            if len(vals) < 7:
                continue
            try:
                ik = int(vals[0])
                ib = int(vals[1])
                v1 = float(vals[2])
                v2 = float(vals[3])
                v3 = float(vals[4])
            except ValueError:
                continue

            kpts.append(ik)
            bands.append(ib)
            vx_list.append(v1)
            vy_list.append(v2)
            vz_list.append(v3)

    if len(kpts) == 0:
        raise RuntimeError(f"No valid velocity data read from {fname}")

    kpts = np.array(kpts, dtype=int)
    bands = np.array(bands, dtype=int)
    vx = np.array(vx_list, dtype=float)
    vy = np.array(vy_list, dtype=float)
    vz = np.array(vz_list, dtype=float)

    nband = int(bands.max())
    Nk = int(kpts.max())

    Vx3d = np.zeros((1, nband, Nk), dtype=float)
    Vy3d = np.zeros((1, nband, Nk), dtype=float)
    Vz3d = np.zeros((1, nband, Nk), dtype=float)

    # convert a.u. -> m/s (NO ROTATE)
    vx_ms = vx * AU_TO_MS
    vy_ms = vy * AU_TO_MS
    vz_ms = vz * AU_TO_MS

    for ik, ib, vxv, vyv, vzv in zip(kpts, bands, vx_ms, vy_ms, vz_ms):
        i_k = int(ik) - 1
        i_b = int(ib) - 1
        Vx3d[0, i_b, i_k] = vxv
        Vy3d[0, i_b, i_k] = vyv
        Vz3d[0, i_b, i_k] = vzv

    # Debug print
    sample_idx = np.linspace(0, len(vx_ms) - 1, min(4, len(vx_ms)), dtype=int)
    print(f"[INFO] Parsed velocities (a.u. → m/s). rotation=False")
    for idx in sample_idx:
        print(f"  idx {idx}: V_raw(a.u.)→m/s = [{vx_ms[idx]:.3e}, {vy_ms[idx]:.3e}, {vz_ms[idx]:.3e}]")

    return {"Vx3d": Vx3d, "Vy3d": Vy3d, "Vz3d": Vz3d, "dims": (1, nband, Nk)}


def parse_effective_mass_file(fname="hole_m_effective.fmt"):
    """
    Minimal parser for hole_m_effective.fmt that extracts the temperature blocks
    and computes an in-plane averaged effective mass M*(T). Returns (temps, mstars).
    """
    if not Path(fname).exists():
        raise FileNotFoundError(fname)

    lines = [ln.rstrip("\n") for ln in open(fname, "r") if ln.strip() and not ln.strip().startswith("#")]
    temps = []
    mstars = []
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) < 4:
            i += 1
            continue
        try:
            T = float(parts[0])
            row0 = [float(x.replace("D", "E")) for x in parts[1:4]]
            row1 = [float(x.replace("D", "E")) for x in lines[i+1].split()]
            row2 = [float(x.replace("D", "E")) for x in lines[i+2].split()]
            tensor = np.array([row0, row1, row2], dtype=float)
            inv_mass_avg = 0.5 * (tensor[0, 0] + tensor[1, 1])
            m_eff = (1.0 / inv_mass_avg) if (inv_mass_avg != 0.0) else np.inf
            temps.append(T)
            mstars.append(m_eff)
            i += 3
        except Exception:
            i += 1
            continue

    temps = np.array(temps, dtype=float)
    mstars = np.array(mstars, dtype=float)

    if temps.size == 0:
        raise RuntimeError(f"No effective mass blocks parsed from {fname}")

    print("   Parsed effective mass tensor blocks:")
    for T, m in zip(temps, mstars):
        print(f"      {T:6.1f} K  →  m* = {m:.3f} m_e")

    return temps, mstars


def write_h5(outname, tau_data, vel_data, mstar_T=None, T_Mstar=None):
    with h5py.File(outname, "w") as h5:
        h5.create_dataset("Gamma_3d", data=tau_data["Gamma3d"])
        h5.create_dataset("Tau_3d", data=tau_data["tau3d"])
        h5.create_dataset("Energies_3d", data=tau_data["E3d"])

        grp = h5.require_group("Velocities")
        grp.create_dataset("Vx3d", data=vel_data["Vx3d"])
        grp.create_dataset("Vy3d", data=vel_data["Vy3d"])
        grp.create_dataset("Vz3d", data=vel_data["Vz3d"])

        # optional effective mass arrays: T_Mstar and Mstar_T (keeps earlier naming)
        if (mstar_T is not None) and (T_Mstar is not None):
            h5.create_dataset("Mstar_T", data=mstar_T)
            h5.create_dataset("T_Mstar", data=T_Mstar)
            h5["Mstar_T"].attrs["unit"] = "m_e"
            h5["Mstar_T"].attrs["description"] = "Conductivity effective mass (in-plane averaged)"

        ntemp, nband, Nk = tau_data["dims"]
        h5.attrs.update({
            "ntemp": ntemp,
            "nband": nband,
            "Nk": Nk,
            "units": "Gamma:1/ps, Tau:ps, Velocity:m/s, Energy:eV"
        })


# ---------- CLI ----------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tau", default="hole_inv_tau_0.fmt")
    p.add_argument("--vel", default="hole_IBTEvel_sup_0.fmt")
    p.add_argument("--mstar", default="hole_m_effective.fmt")
    p.add_argument("--out", default="epw_preproc.h5")
    args = p.parse_args()

    print("=== EPW parser & preprocessor (NO ROTATION) ===")

    print("Parsing tau/energy file...")
    tau_data = parse_tau_file(args.tau)
    print(f"  dims: ntemp={tau_data['dims'][0]}, nband={tau_data['dims'][1]}, Nk={tau_data['dims'][2]}")

    print("Parsing velocity file...")
    vel_data = parse_velocity_file(args.vel)

    # try to parse effective mass file if present, add T arrays to HDF5
    mstar_T = None
    T_Mstar = None
    if Path(args.mstar).exists():
        try:
            T_Mstar, mstar_T = parse_effective_mass_file(args.mstar)
            print(f"   Parsed effective mass file and will add T_Mstar / Mstar_T to HDF5.")
        except Exception as e:
            print(f"[WARN] Could not parse {args.mstar}: {e}")


    print("Writing HDF5 ...")
    write_h5(args.out, tau_data, vel_data, mstar_T, T_Mstar)
    print(f"Saved to {args.out}")

if __name__ == "__main__":
    main()

