#!/usr/bin/env python3
"""
EPW Monte Carlo Parser and Preprocessor
---------------------------------------
Reads EPW output files:
    - electron_inv_taucb_0.fmt   (relaxation times)
    - electron_IBTEvel_sup_0.fmt (band velocities and energies)
Builds combined arrays of:
    Energies [eV], Velocities [m/s], and Scattering rates [1/ps]
and saves them to epw_preproc.h5.

This file is meant to be run once before launching the MC driver.
"""

import numpy as np
import h5py
import re
from pathlib import Path

# --- Physical constants ---
RY_TO_EV = 13.605693
BOHR_VEL = 2.1876912633e6     # a.u. velocity -> m/s
RY_TAU_TO_PSEC = 20670.6944033  # from EPW header

def parse_tau_file(fname="electron_inv_taucb_0.fmt"):
    """Parse relaxation times vs (temp, kpt, band) from EPW format file."""
    temps, kpts, bands, energy, tau = [], [], [], [], []

    with open(fname, "r") as f:
        for line in f:
            # Skip headers and empty lines
            if line.strip().startswith("#") or not line.strip():
                continue

            vals = line.split()
            # Your file has: itemp, kpt, ibnd, energy[Ry], tau[Ry]
            if len(vals) < 5:
                continue

            try:
                itemp = int(vals[0])
                ik = int(vals[1])
                ib = int(vals[2])
                Ery = float(vals[3])
                tRy = float(vals[4])
            except ValueError:
                # skip lines with any parsing issue
                continue

            temps.append(itemp)
            kpts.append(ik)
            bands.append(ib)
            energy.append(Ery * RY_TO_EV)  # convert Ry → eV
            tau.append(tRy * RY_TAU_TO_PSEC)  # convert Ry → ps

    # Convert to numpy arrays
    temps = np.array(temps)
    kpts = np.array(kpts)
    bands = np.array(bands)
    energy = np.array(energy)
    tau = np.array(tau)

    # Safety check: ensure file had valid data
    if len(temps) == 0:
        raise RuntimeError(
            f"[ERROR] No relaxation time data found in {fname}. "
            "Check file structure or header format."
        )

    # Determine dimensions
    ntemp = temps.max()
    nband = bands.max()
    Nk = kpts.max()

    # Initialize arrays
    tau3d = np.zeros((ntemp, nband, Nk))
    E3d = np.zeros((ntemp, nband, Nk))

    # Fill arrays
    for t, k, b, e, val in zip(temps, kpts, bands, energy, tau):
        tau3d[t - 1, b - 1, k - 1] = val
        E3d[t - 1, b - 1, k - 1] = e

    # Compute scattering rate Γ = 1/τ
    with np.errstate(divide="ignore", invalid="ignore"):
        Gamma3d = np.where(tau3d > 0, 1.0 / tau3d, 0.0)

    return {
        "tau3d": tau3d,
        "Gamma3d": Gamma3d,
        "E3d": E3d,
        "dims": (ntemp, nband, Nk),
    }


def parse_velocity_file(fname="electron_IBTEvel_sup_0.fmt"):
    """Parse velocity file (ik, ibnd, vx, vy, vz, eig, weight) and convert from a.u. → m/s."""
    import numpy as np

    AU_TO_MS =  2.18769126364e6  # 1 a.u. velocity = 2.18769126364e6 m/s

    kpts, bands = [], []
    vx, vy, vz = [], [], []

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
                v1, v2, v3 = float(vals[2]), float(vals[3]), float(vals[4])
            except ValueError:
                continue

            # Convert each component from a.u. → m/s
            v1 *= AU_TO_MS
            v2 *= AU_TO_MS
            v3 *= AU_TO_MS

            kpts.append(ik)
            bands.append(ib)
            vx.append(v1)
            vy.append(v2)
            vz.append(v3)

    if not kpts:
        raise RuntimeError(f"[ERROR] No valid data read from {fname}. Check file structure.")

    kpts = np.array(kpts)
    bands = np.array(bands)
    vx = np.array(vx)
    vy = np.array(vy)
    vz = np.array(vz)

    nband = int(bands.max())
    Nk = int(kpts.max())

    Vx3d = np.zeros((1, nband, Nk))
    Vy3d = np.zeros((1, nband, Nk))
    Vz3d = np.zeros((1, nband, Nk))

    for k, b, vx_, vy_, vz_ in zip(kpts, bands, vx, vy, vz):
        # k and b are 1-based in the input files; convert to 0-based indices
        Vx3d[0, int(b) - 1, int(k) - 1] = vx_
        Vy3d[0, int(b) - 1, int(k) - 1] = vy_
        Vz3d[0, int(b) - 1, int(k) - 1] = vz_

    # Debug: print sample values to confirm units, but print as Band, kpt
    n_entries = len(vx)
    n_samples = min(3, n_entries)
    sample_idx = np.random.choice(n_entries, size=n_samples, replace=False)

    print(f"[INFO] Velocity data successfully parsed and converted from a.u. → m/s")
    for idx in sample_idx:
        ib = int(bands[idx])
        ik = int(kpts[idx])
        print(f"    Band {ib}, kpt {ik}: Vx = {vx[idx]:.3e} m/s, Vy = {vy[idx]:.3e} m/s, Vz = {vz[idx]:.3e} m/s")

    print(f"[DEBUG] Parsed {Nk} k-points × {nband} bands from {fname}\n")

    return {
        "Vx3d": Vx3d,
        "Vy3d": Vy3d,
        "Vz3d": Vz3d,
        "dims": (1, nband, Nk)
    }



def parse_effective_mass_file(fname="electron_m_effective.fmt"):
    temps, mstars = [], []
    lines = [l for l in open(fname) if l.strip() and not l.startswith("#")]
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        try:
            # First value in line is temperature
            T = float(parts[0])
            # Read 3x3 tensor: first line may have extra cols
            tensor = []
            tensor.append([float(x.replace("E","e")) for x in parts[1:4]])
            for j in range(1, 3):  # read next 2 lines
                tensor.append([float(x.replace("E","e")) for x in lines[i+j].split()])
            tensor = np.array(tensor)
            # In-plane average
            inv_mass_avg = 0.5 * (tensor[0,0] + tensor[1,1])
            m_eff = 1.0 / inv_mass_avg if inv_mass_avg != 0 else np.inf
            temps.append(T)
            mstars.append(m_eff)
            i += 3  # move to next block
        except:
            i += 1

    temps = np.array(temps)
    mstars = np.array(mstars)
    print("   Parsed effective mass tensor blocks:")
    for T, m in zip(temps, mstars):
        print(f"      {T:6.1f} K  →  m* = {m:.3f} m_e")
    return temps, mstars




def quick_build_example():
    """Main entry point."""
    tau_data = parse_tau_file()
    vel_data = parse_velocity_file()

    outname = Path("epw_preproc.h5")
    with h5py.File(outname, "w") as h5:
        # Core quantities
        h5.create_dataset("Gamma_3d", data=tau_data["Gamma3d"])
        h5.create_dataset("Tau_3d", data=tau_data["tau3d"])
        h5.create_dataset("Energies_3d", data=tau_data["E3d"])

        # Velocity components
        h5.create_dataset("Velocities/Vx3d", data=vel_data["Vx3d"])
        h5.create_dataset("Velocities/Vy3d", data=vel_data["Vy3d"])
        h5.create_dataset("Velocities/Vz3d", data=vel_data["Vz3d"])

        # Flattened arrays
        h5.create_dataset("Gamma_flat", data=tau_data["Gamma3d"].ravel())
        h5.create_dataset("Tau_flat", data=tau_data["tau3d"].ravel())

        # Parse and add effective mass vs temperature
        try:
            Tvals, Mstar_T = parse_effective_mass_file()
            h5.create_dataset("Mstar_T", data=Mstar_T)
            h5.create_dataset("T_Mstar", data=Tvals)
            h5["Mstar_T"].attrs["unit"] = "m_e"
            h5["Mstar_T"].attrs["description"] = "Conductivity effective mass (in-plane averaged)"
            print(f"   Added M*(T) with {len(Tvals)} temperature points.")
        except Exception as e:
            print(f"[WARN] Skipped effective mass file: {e}")

        # Metadata
        ntemp, nband, Nk = tau_data["dims"]
        h5.attrs.update({
            "ntemp": ntemp,
            "nband": nband,
            "Nk": Nk,
            "units": "Gamma: 1/ps, Tau: ps, Velocity: m/s, Energy: eV, M*: m_e",
            "description": "EPW-preprocessed data for Monte Carlo transport"
        })

    print(f"✅ Saved combined data to {outname} with shape (ntemp={ntemp}, nband={nband}, Nk={Nk})")
    print("   Gamma[1/ps], Tau[ps], Velocity[Vx,Vy,Vz in m/s], Energies[eV], and M*(T) written.")


    # ---- DEBUG SUMMARY ----
    Tvals = None
    try:
        with h5py.File(outname, "r") as h5:
            if "T_Mstar" in h5:
                Tvals = h5["T_Mstar"][:]
    except Exception:
        pass

    # Always use the first T index (itemp = 0)
    T_display = Tvals[0] if Tvals is not None and len(Tvals) > 0 else 0.0
    print(f"\n[DEBUG] Data summary for T = {T_display:.2f} K")
    print("Format: Band idx | Energy [eV] | Gamma [1/ps] | Velocities [m/s]")

    # Select first few bands and k-points for readable output
    n_show_bands = min(3, nband)
    n_show_kpts = min(5, Nk)

    for ib in range(n_show_bands):
        for ik in range(n_show_kpts):
            E = tau_data["E3d"][0, ib, ik]
            Gamma = tau_data["Gamma3d"][0, ib, ik]
            Vx = vel_data["Vx3d"][0, ib, ik]
            Vy = vel_data["Vy3d"][0, ib, ik]
            Vz = vel_data["Vz3d"][0, ib, ik]
            print(
                f"Band {ib+1}, kpt {ik+1}: "
                f"E = {E:8.4f} | Gamma = {Gamma:9.4e} | "
                f"V = [{Vx:8.2e}, {Vy:8.2e}, {Vz:8.2e}]"
            )

if __name__ == "__main__":
    print("This is a parser + MC-skeleton file. Import functions and call quick_build_example().")

