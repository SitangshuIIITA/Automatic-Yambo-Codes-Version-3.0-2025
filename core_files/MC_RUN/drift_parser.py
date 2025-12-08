#!/usr/bin/env python3
"""
parse_epw_quadrupole.py

Parse EPW quadrupole output (epw_quadrupole.out) together with q.dat (k-path
fractional coordinates) and write a compact HDF5 containing:
 - k-path coords (in order from q.dat)
 - discovered iq coords
 - band list and energies E(band, ik)
 - scattering rows with fields:
    (ik_pos, iq_pos, ibnd, jbnd, imode, omega_meV, g_meV, Re_meV, Im_meV, enk_eV, enkq_eV)
 - meta: Fermi energy, ibndmin/ebndmax (if present), band energy stats

New: computes group velocities from energies and reciprocal lattice vectors
and stores them at /bands/velocities_m_s (nb x Nk x 3). Also writes parsed_log.txt.

Usage:
    python parse_epw_quadrupole.py --quad epw_quadrupole.out --qdat q.dat --out epw_quadrupole_parsed.h5
"""
import argparse
import re
import numpy as np
import h5py
from collections import defaultdict, OrderedDict
import sys
import os
import math

# Physical constants
eV_to_J = 1.602176634e-19
hbar = 1.054571817e-34  # Js

# ---------------------------
# Regex patterns (robust)
# ---------------------------
IK_RE = re.compile(r"\bik\s*=\s*(\d+)\s+coord\.*:\s*([-\d\.\sEe+]+)")
IQ_RE = re.compile(r"\biq\s*=\s*(\d+)\s+coord\.*:\s*([-\d\.\sEe+]+)")
FERMI_RE = re.compile(r"Fermi energy .*?=\s*([-\d\.\sEe+]+)\s*eV", re.IGNORECASE)
IBNDMIN_RE = re.compile(r"ibndmin\s*=\s*(\d+)\s+ebndmin\s*=\s*([^\n]+)", re.IGNORECASE)
IBNDMAX_RE = re.compile(r"ibndmax\s*=\s*(\d+)\s+ebndmax\s*=\s*([^\n]+)", re.IGNORECASE)

ROW_RE = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+([-\d\.Ee+]+)\s+([-\d\.Ee+]+)\s+([-\d\.Ee+]+)\s+([-\d\.Ee+]+)\s+([-\d\.Ee+]+)\s+([-\d\.Ee+]+)\s+([-\d\.Ee+]+)"
)

# ---------------------------
# Helpers
# ---------------------------
def read_recip_from_scf(scf_path, verbose=True):
    """
    Reads QE scf.out and returns reciprocal vectors (3x3) in 1/m.
    """
    alat_bohr = None
    bvec_bohr = []

    with open(scf_path, "r") as f:
        for line in f:
            if "celldm(1)" in line:
                parts = line.replace("=", " ").split()
                for i,p in enumerate(parts):
                    if p.strip().startswith("celldm(1)"):
                        try:
                            alat_bohr = float(parts[i+1])
                        except:
                            pass
            if "reciprocal axes:" in line:
                for _ in range(3):
                    L = next(f)
                    nums = re.findall(r"[-+]?\d*\.\d+|\d+", L)
                    if len(nums) >= 3:
                        bvec_bohr.append([float(nums[0]), float(nums[1]), float(nums[2])])

    if alat_bohr is None or len(bvec_bohr) != 3:
        raise RuntimeError("Could not parse reciprocal lattice vectors from scf.out")

    BOHR_TO_M = 0.52917721092e-10
    alat_m = alat_bohr * BOHR_TO_M
    factor = (2.0 * np.pi) / alat_m
    bvec = np.array(bvec_bohr, dtype=float) * factor  # now in 1/m

    if verbose:
        print("[info] Parsed reciprocal lattice vectors (1/m):")
        print(bvec)

    return bvec

def read_qdat(qdat_path):
    """Read q.dat (k-path file) and return ordered fractional coords (N,3)."""
    coords = []
    with open(qdat_path, "r") as fh:
        first = fh.readline()
        def parse_line(line):
            parts = line.strip().split()
            if len(parts) < 3:
                return None
            try:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                return (x, y, z)
            except:
                return None

        maybe = parse_line(first)
        if maybe is not None:
            coords.append(maybe)
            for line in fh:
                p = parse_line(line)
                if p:
                    coords.append(p)
        else:
            for line in fh:
                p = parse_line(line)
                if p:
                    coords.append(p)
    coords = np.array(coords, dtype=float)
    if coords.ndim == 1 and coords.size == 3:
        coords = coords.reshape((1,3))
    return coords

def parse_quadrupole_file(fname, verbose=True):
    """
    Parse epw.out (streaming).
    """
    entries = []
    ik_ordered = []
    iq_ordered = []
    ik_coords = {}
    iq_coords = {}
    band_set = set()
    fermi = None
    ibndmin = None
    ebndmin = None
    ibndmax = None
    ebndmax = None

    cur_ik = None
    cur_iq = None

    with open(fname, "r") as fh:
        for line in fh:
            if fermi is None:
                m = FERMI_RE.search(line)
                if m:
                    try:
                        fermi = float(m.group(1))
                        if verbose:
                            print(f"[meta] Found QE/EPW's Fermi energy: {fermi} eV")
                    except:
                        pass
            if ibndmin is None:
               m = IBNDMIN_RE.search(line)
               if m:
                   ibndmin = int(m.group(1))
                   numtxt = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", m.group(2))
                   ebndmin = float(numtxt[0]) if numtxt else None
                   if verbose:
                       print(f"[meta] Found EPW's ibndmin={ibndmin}, ebndmin={ebndmin} eV")
            if ibndmax is None:
               m = IBNDMAX_RE.search(line)
               if m:
                   ibndmax = int(m.group(1))
                   numtxt = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", m.group(2))
                   ebndmax = float(numtxt[0]) if numtxt else None
                   if verbose:
                       print(f"[meta] Found EPW's ibndmax={ibndmax}, ebndmax={ebndmax} eV")

            m_ik = IK_RE.search(line)
            if m_ik:
                cur_ik = int(m_ik.group(1))
                coords = [float(x) for x in m_ik.group(2).strip().split()[:3]]
                if cur_ik not in ik_coords:
                    ik_ordered.append(cur_ik)
                    ik_coords[cur_ik] = np.array(coords, dtype=float)
                continue
            m_iq = IQ_RE.search(line)
            if m_iq:
                cur_iq = int(m_iq.group(1))
                coords = [float(x) for x in m_iq.group(2).strip().split()[:3]]
                if cur_iq not in iq_coords:
                    iq_ordered.append(cur_iq)
                    iq_coords[cur_iq] = np.array(coords, dtype=float)
                continue

            m_row = ROW_RE.match(line)
            if m_row and (cur_ik is not None) and (cur_iq is not None):
                ibnd = int(m_row.group(1))
                jbnd = int(m_row.group(2))
                imode = int(m_row.group(3))
                enk = float(m_row.group(4))
                enkq = float(m_row.group(5))
                omega = float(m_row.group(6))
                g_sym = float(m_row.group(7))
                g_abs = float(m_row.group(8))
                re_g = float(m_row.group(9))
                im_g = float(m_row.group(10))

                band_set.add(ibnd)
                band_set.add(jbnd)

                entries.append({
                    "ik": int(cur_ik),
                    "iq": int(cur_iq),
                    "ibnd": ibnd,
                    "jbnd": jbnd,
                    "imode": imode,
                    "enk": enk,
                    "enkq": enkq,
                    "omega_meV": omega,
                    "g_sym_meV": g_sym,
                    "g_abs_meV": g_abs,
                    "re_meV": re_g,
                    "im_meV": im_g
                })
                continue
    band_list = sorted(band_set)
    return {
        "rows": entries,
        "ik_ordered": ik_ordered,
        "ik_coords": ik_coords,
        "iq_ordered": iq_ordered,
        "iq_coords": iq_coords,
        "band_list": band_list,
        "fermi": fermi,
        "ibndmin": ibndmin,
        "ebndmin": ebndmin,
        "ibndmax": ibndmax,
        "ebndmax": ebndmax
    }

# ---------------------------
# VELOCITY computation
# ---------------------------
def frac_to_cart(frac_coords, bvecs):
    """frac_coords (N,3) * bvecs (3,3) -> cart k (N,3)"""
    return np.dot(frac_coords, bvecs)

def compute_velocities(energies_eV, kpath_frac, bvecs):
    """
    Compute v = (1/ħ) ∇_k E for each band and k-point.
    energies_eV: (nb, Nk) array with NaNs for missing entries
    kpath_frac: (Nk,3) fractional coords
    bvecs: (3,3) reciprocal vectors (1/m)

    Returns velocities_m_s: (nb, Nk, 3)
    """
    nb, Nk = energies_eV.shape
    # convert kpath to cartesian (1/m)
    kcart = frac_to_cart(kpath_frac, bvecs)  # Nk x 3
    vel = np.zeros((nb, Nk, 3), dtype=float)

    for ib in range(nb):
        E = energies_eV[ib, :]  # eV
        for ik in range(Nk):
            # find left and right indices with valid energies
            # prefer central difference il < ik < ir
            il = ik - 1
            ir = ik + 1
            # move left until valid or out
            while il >= 0 and np.isnan(E[il]):
                il -= 1
            while ir < Nk and np.isnan(E[ir]):
                ir += 1
            # if both valid and distinct, do central-ish diff using those
            if il >= 0 and ir < Nk and il != ir:
                # ensure k coords differ
                dkx = kcart[ir,0] - kcart[il,0]
                dky = kcart[ir,1] - kcart[il,1]
                # if both dk components are zero (degenerate coords), fallback
                if abs(dkx) < 1e-40 and abs(dky) < 1e-40:
                    vel[ib, ik, :] = 0.0
                    continue
                # compute derivative components carefully:
                # dE/dkx approx (E_ir - E_il) / (k_ir_x - k_il_x)
                # handle division by zero component-wise
                v_x = 0.0
                v_y = 0.0
                # convert energy difference to J
                dE_J = (E[ir] - E[il]) * eV_to_J
                if abs(dkx) > 1e-40:
                    dEdkx = (E[ir] - E[il]) * eV_to_J / dkx
                    v_x = dEdkx / hbar
                if abs(dky) > 1e-40:
                    dEdky = (E[ir] - E[il]) * eV_to_J / dky
                    v_y = dEdky / hbar
                vel[ib, ik, 0] = v_x
                vel[ib, ik, 1] = v_y
                vel[ib, ik, 2] = 0.0
            else:
                # fallback: one-sided difference
                if il >= 0:
                    dkx = kcart[ik,0] - kcart[il,0]
                    dky = kcart[ik,1] - kcart[il,1]
                    if (abs(dkx) < 1e-40 and abs(dky) < 1e-40) or np.isnan(E[il]) or np.isnan(E[ik]):
                        vel[ib, ik, :] = 0.0
                        continue
                    v_x = 0.0
                    v_y = 0.0
                    if abs(dkx) > 1e-40:
                        dEdkx = (E[ik] - E[il]) * eV_to_J / dkx
                        v_x = dEdkx / hbar
                    if abs(dky) > 1e-40:
                        dEdky = (E[ik] - E[il]) * eV_to_J / dky
                        v_y = dEdky / hbar
                    vel[ib, ik, 0] = v_x
                    vel[ib, ik, 1] = v_y
                    vel[ib, ik, 2] = 0.0
                elif ir < Nk:
                    dkx = kcart[ir,0] - kcart[ik,0]
                    dky = kcart[ir,1] - kcart[ik,1]
                    if (abs(dkx) < 1e-40 and abs(dky) < 1e-40) or np.isnan(E[ir]) or np.isnan(E[ik]):
                        vel[ib, ik, :] = 0.0
                        continue
                    v_x = 0.0
                    v_y = 0.0
                    if abs(dkx) > 1e-40:
                        dEdkx = (E[ir] - E[ik]) * eV_to_J / dkx
                        v_x = dEdkx / hbar
                    if abs(dky) > 1e-40:
                        dEdky = (E[ir] - E[ik]) * eV_to_J / dky
                        v_y = dEdky / hbar
                    vel[ib, ik, 0] = v_x
                    vel[ib, ik, 1] = v_y
                    vel[ib, ik, 2] = 0.0
                else:
                    vel[ib, ik, :] = 0.0
    return vel

# ---------------------------
# Write human-readable parsed_log.txt
# ---------------------------
def write_parsed_log(out_h5, parsed, qdat_coords, recip_vectors, kpath_coords,
                     q_coords, energies, velocities, rows_array, band_list,
                     summary, ik_ordered, iq_ordered, ik_pos_map, iq_pos_map,
                     n_k=10, n_q=10, n_energies=8, n_rows=20, n_vel=8):

    """
    Write a human-readable parsed_log.txt next to out_h5 summarizing key contents.
    """
    out_dir = os.path.dirname(os.path.abspath(out_h5))
    log_path = os.path.join(out_dir, "report_parsed.txt")
    with open(log_path, "w") as fh:
        fh.write("Parsed EPW summary log\n")
        fh.write("Generated by parse_epw.py\n")
        fh.write("===================================\n\n")

        # reciprocal lattice
        fh.write("===== RECIPROCAL LATTICE (b_vectors) =====\n")
        if recip_vectors is not None:
            fh.write("Units: 1/m\n")
            for i, vec in enumerate(recip_vectors):
                fh.write(f"b[{i}] = {vec[0]: .6e}  {vec[1]: .6e}  {vec[2]: .6e}\n")
        else:
            fh.write("reciprocal_lattice/b_vectors: NOT STORED (scf.out not found)\n")
        fh.write("\n")

        # k-points
        fh.write("===== K-POINTS (kpath/coords) (fractional) =====\n")
        Nk = kpath_coords.shape[0]
        fh.write(f"Total Nk = {Nk}\n")
        fh.write("index |   kx       |   ky       |   kz\n")
        fh.write("---------------------------------------\n")
        for i in range(min(n_k, Nk)):
            kx, ky, kz = kpath_coords[i]
            fh.write(f"{i:5d} | {kx: .6f} | {ky: .6f} | {kz: .6f}\n")
        if Nk > n_k:
            fh.write(f"... ({Nk - n_k} more k-points)\n")
        fh.write("\n")

        # q-points
        fh.write("===== Q-POINTS (q/coords) (fractional) =====\n")
        Nq = q_coords.shape[0]
        fh.write(f"Total Nq = {Nq}\n")
        fh.write("index |   qx       |   qy       |   qz\n")
        fh.write("---------------------------------------\n")
        for i in range(min(n_q, Nq)):
            qx, qy, qz = q_coords[i]
            fh.write(f"{i:5d} | {qx: .6f} | {qy: .6f} | {qz: .6f}\n")
        if Nq > n_q:
            fh.write(f"... ({Nq - n_q} more q-points)\n")
        fh.write("\n")

	# ===== IK / IQ INDEX MAPS =====
        fh.write("===== K-POINT INDEX MAP =====\n")
        fh.write("ik_raw  ik_pos   kx         ky         kz\n")
        fh.write("----------------------------------------------\n")
        for ik_raw in ik_ordered:
            pos = ik_pos_map[ik_raw]
            kx, ky, kz = kpath_coords[pos]
            fh.write(f"{ik_raw:6d}  {pos:6d}   {kx:9.6f}  {ky:9.6f}  {kz:9.6f}\n")
        fh.write("\n")

        fh.write("===== Q-POINT INDEX MAP =====\n")
        fh.write("iq_raw  iq_pos   qx         qy         qz\n")
        fh.write("----------------------------------------------\n")
        for iq_raw in iq_ordered:
            pos = iq_pos_map[iq_raw]
            qx, qy, qz = q_coords[pos]
            fh.write(f"{iq_raw:6d}  {pos:6d}   {qx:9.6f}  {qy:9.6f}  {qz:9.6f}\n")
        fh.write("\n")




        # band list
        fh.write("===== BAND LIST =====\n")
        fh.write("band_list (nbands) = " + ", ".join([str(int(b)) for b in band_list]) + "\n")
        fh.write("\n")

        # sample energies per band (first n_energies kpoints)
        fh.write(f"===== SAMPLE ENERGIES (first {n_energies} k-points) per band (energies_eV) =====\n")
        nb = energies.shape[0]
        Nk = energies.shape[1]
        fh.write("band | " + " | ".join([f"k{j}" for j in range(min(n_energies, Nk))]) + "\n")
        fh.write("-" * (10 + 14*min(n_energies, Nk)) + "\n")
        for ib in range(nb):
            row_vals = []
            for ik in range(min(n_energies, Nk)):
                val = energies[ib, ik]
                if np.isnan(val):
                    row_vals.append("   NaN   ")
                else:
                    row_vals.append(f"{val:8.5f}")
            fh.write(f"{int(band_list[ib]):4d} | " + " | ".join(row_vals) + "\n")
        fh.write("\n")


	# sample velocities (swap axes: k-points as rows, bands as columns)
        fh.write(f"===== SAMPLE VELOCITIES (m/s) (first {n_vel} k-points) =====\n")

        # Header
        fh.write("k \\ band |")
        for ib in range(nb):
            fh.write(f"{'b'+str(int(band_list[ib])):^24s}")
        fh.write("\n")
        fh.write("-" * (10 + 24*nb) + "\n")

        # Rows: loop over k points
        for ik in range(min(n_vel, Nk)):
            fh.write(f"k{ik:<8d}|")
            for ib in range(nb):
                vx, vy, vz = velocities[ib, ik, :]
                fh.write(f"({vx:>10.2e},{vy:>10.2e}) ")
            fh.write("\n")

        fh.write("\n")

      
        # meta: fermi, ibndmin/ebndmax
        fh.write("===== META =====\n")
        fermi = parsed.get("fermi", None)
        if fermi is not None:
            fh.write(f"Fermi energy (E_fermi_eV) = {fermi:.6f} eV\n")
        else:
            fh.write("Fermi energy: NOT FOUND\n")
        ibndmin = parsed.get("ibndmin", None)
        ebndmin = parsed.get("ebndmin", None)
        ibndmax = parsed.get("ibndmax", None)
        ebndmax = parsed.get("ebndmax", None)
        fh.write(f"ibndmin = {ibndmin}, ebndmin = {ebndmin}\n")
        fh.write(f"ibndmax = {ibndmax}, ebndmax = {ebndmax}\n")
        fh.write("\n")

        # scattering rows (first n_rows)
        fh.write(f"===== SCATTERING ROWS (first {n_rows}) =====\n")
        fh.write("row | ik_pos | iq_pos | ibnd | jbnd | imode | omega(meV) | g_abs(meV) | re_meV | im_meV | enk_eV | enkq_eV\n")
        fh.write("-" * 120 + "\n")
        Nr = rows_array.shape[0]
        for i in range(min(n_rows, Nr)):
            r = rows_array[i]
            fh.write(f"{i:4d} | {int(r['ik_pos']):6d} | {int(r['iq_pos']):6d} | {int(r['ibnd']):4d} | {int(r['jbnd']):4d} | {int(r['imode']):5d} | "
                     f"{r['omega_meV']:10.5f} | {r['g_abs_meV']:10.5f} | {r['re_meV']:8.5f} | {r['im_meV']:8.5f} | {r['enk_eV']:8.5f} | {r['enkq_eV']:8.5f}\n")
        if Nr > n_rows:
            fh.write(f"... ({Nr - n_rows} more rows)\n")
        fh.write("\n")

        # summary
        fh.write("===== SUMMARY =====\n")
        fh.write(f"n_rows = {summary.get('n_rows', Nr)}\n")
        fh.write(f"nbands = {summary.get('nbands', nb)}\n")
        fh.write(f"Nk = {summary.get('Nk', Nk)}\n")
        fh.write(f"Nq = {summary.get('Nq', Nq)}\n")
        fh.write("\n")
    print(f"[info] Wrote human-readable log to report_parsed.txt")

# ---------------------------
# Build arrays and HDF5
# ---------------------------
def build_and_write_h5(parsed, qdat_coords, out_h5, recip_vectors=None, tol=1e-6, compress=True, verbose=True):
    rows = parsed["rows"]
    ik_ordered = parsed["ik_ordered"]
    iq_ordered = parsed["iq_ordered"]
    ik_coords_map = parsed["ik_coords"]
    iq_coords_map = parsed["iq_coords"]
    band_list = parsed["band_list"]
    fermi = parsed["fermi"]
    ibndmin = parsed["ibndmin"]
    ebndmin = parsed["ebndmin"]
    ibndmax = parsed["ibndmax"]
    ebndmax = parsed["ebndmax"]

    Nk_qdat = qdat_coords.shape[0]
    Nk_file = len(ik_ordered)
    if Nk_qdat == Nk_file:
        ik_pos_map = {ik_ordered[i]: i for i in range(Nk_file)}
        kpath_coords = qdat_coords.copy()
    else:
        if verbose:
            print(f"[warn] q.dat length ({Nk_qdat}) != ik count ({Nk_file}). Using ik encountered order.")
        ik_pos_map = {ik_ordered[i]: i for i in range(Nk_file)}
        kpath_coords = np.array([ik_coords_map[ik] for ik in ik_ordered], dtype=float)

    Nq = len(iq_ordered)
    iq_pos_map = {iq_ordered[i]: i for i in range(Nq)}
    q_coords = np.array([iq_coords_map[iq] for iq in iq_ordered], dtype=float) if Nq > 0 else np.zeros((0,3),dtype=float)

    nb = len(band_list)
    band_to_idx = {b: i for i, b in enumerate(band_list)}

    Nk = len(ik_pos_map)
    energies = np.full((nb, Nk), np.nan, dtype=float)

    for r in rows:
        ibnd = r["ibnd"]
        ik = r["ik"]
        if ik not in ik_pos_map:
            continue
        ib_idx = band_to_idx[ibnd]
        ikp = ik_pos_map[ik]
        if np.isnan(energies[ib_idx, ikp]):
            energies[ib_idx, ikp] = r["enk"]

    nmissing = np.isnan(energies).sum()
    if nmissing > 0 and verbose:
        print(f"[warn] {nmissing} missing band-energy entries (nb={nb}, Nk={Nk}).")

    band_stats = {}
    E_min = np.nanmin(energies, axis=1)
    E_max = np.nanmax(energies, axis=1)
    E_mean = np.nanmean(energies, axis=1)
    for i, b in enumerate(band_list):
        band_stats[b] = {"E_min": float(E_min[i]), "E_max": float(E_max[i]), "E_mean": float(E_mean[i])}

    valence_bands = []
    conduction_bands = []
    ambiguous_bands = []
    if fermi is not None:
        for i, b in enumerate(band_list):
            if E_max[i] <= (fermi + tol):
                valence_bands.append(int(b))
            elif E_min[i] >= (fermi - tol):
                conduction_bands.append(int(b))
            else:
                ambiguous_bands.append(int(b))
    else:
        ambiguous_bands = [int(b) for b in band_list]

    dtype_row = np.dtype([
        ("ik_pos", np.int32),
        ("iq_pos", np.int32),
        ("ibnd", np.int32),
        ("jbnd", np.int32),
        ("imode", np.int32),
        ("omega_meV", np.float64),
        ("g_abs_meV", np.float64),
        ("re_meV", np.float64),
        ("im_meV", np.float64),
        ("enk_eV", np.float64),
        ("enkq_eV", np.float64),
        ("k_frac", np.float64, (3,)),      
        ("kq_frac", np.float64, (3,)),     
    ])
    rows_array = np.zeros((len(rows),), dtype=dtype_row)
    for i, r in enumerate(rows):
        ik = r["ik"]
        iq = r["iq"]
        rows_array[i]["ik_pos"] = ik_pos_map.get(ik, -1)
        rows_array[i]["iq_pos"] = iq_pos_map.get(iq, -1)
        rows_array[i]["ibnd"] = int(r["ibnd"])
        rows_array[i]["jbnd"] = int(r["jbnd"])
        rows_array[i]["imode"] = int(r["imode"])
        rows_array[i]["omega_meV"] = float(r["omega_meV"])
        rows_array[i]["g_abs_meV"] = float(r["g_abs_meV"])
        rows_array[i]["re_meV"] = float(r["re_meV"])
        rows_array[i]["im_meV"] = float(r["im_meV"])
        rows_array[i]["enk_eV"] = float(r["enk"])
        rows_array[i]["enkq_eV"] = float(r["enkq"])

            # --- ADD K AND K+Q FRACTIONAL COORDINATES ---
        k_frac = kpath_coords[ rows_array[i]["ik_pos"] ]     # (3,)
        q_frac = q_coords[ rows_array[i]["iq_pos"] ]         # (3,)
        kq_frac = (k_frac + q_frac) % 1.0                    # wrap in BZ

        rows_array[i]["k_frac"] = k_frac
        rows_array[i]["kq_frac"] = kq_frac



    # compute velocities if reciprocal lattice present
    velocities = np.zeros((nb, Nk, 3), dtype=float)
    if recip_vectors is not None:
        try:
            velocities = compute_velocities(energies, kpath_coords, recip_vectors)
        except Exception as e:
            print(f"[warn] velocity computation failed: {e}. velocities will be zeros.")

    if verbose:
        print(f"[info] Writing HDF5 to {out_h5} ... (rows={len(rows_array)}, nb={nb}, Nk={Nk}, Nq={Nq})")
    with h5py.File(out_h5, "w") as fh:
        if recip_vectors is not None:
           rl = fh.create_group("reciprocal_lattice")
           rl.create_dataset("b_vectors", data=recip_vectors)
           rl.attrs["units"] = "1/m"

        meta = fh.create_group("meta")
        meta.attrs["source_file"] = str(out_h5)
        meta.attrs["parsed_from"] = "epw.out"
        if fermi is not None:
            meta.attrs["E_fermi_eV"] = float(fermi)
        if ibndmin is not None:
            meta.attrs["ibndmin"] = int(ibndmin)
            meta.attrs["ebndmin_eV"] = float(ebndmin)
        if ibndmax is not None:
            meta.attrs["ibndmax"] = int(ibndmax)
            meta.attrs["ebndmax_eV"] = float(ebndmax)

        grp_k = fh.create_group("kpath")
        grp_k.create_dataset("coords", data=kpath_coords, compression="gzip" if compress else None)
        grp_k.create_dataset("ik_list", data=np.array(ik_ordered, dtype=np.int32), compression="gzip" if compress else None)

        grp_q = fh.create_group("q")
        grp_q.create_dataset("coords", data=q_coords, compression="gzip" if compress else None)
        grp_q.create_dataset("iq_list", data=np.array(iq_ordered, dtype=np.int32), compression="gzip" if compress else None)

        grp_b = fh.create_group("bands")
        grp_b.create_dataset("band_list", data=np.array(band_list, dtype=np.int32))
        grp_b.create_dataset("energies_eV", data=energies, compression="gzip" if compress else None)
        grp_b.create_dataset("velocities_m_s", data=velocities, compression="gzip" if compress else None)

        bs = grp_b.create_group("band_stats")
        for b in band_list:
            g = bs.create_group(f"band_{int(b)}")
            st = band_stats[b]
            g.attrs["E_min_eV"] = st["E_min"]
            g.attrs["E_max_eV"] = st["E_max"]
            g.attrs["E_mean_eV"] = st["E_mean"]

        cls = meta.create_group("classification")
        cls.create_dataset("valence_bands", data=np.array(valence_bands, dtype=np.int32))
        cls.create_dataset("conduction_bands", data=np.array(conduction_bands, dtype=np.int32))
        cls.create_dataset("ambiguous_bands", data=np.array(ambiguous_bands, dtype=np.int32))

        scat = fh.create_group("scattering")
        scat.create_dataset("rows", data=rows_array, compression="gzip" if compress else None)
           # --- Store raw ik/iq lists and index maps for MC k+q reconstruction ---
        scat.create_dataset("ik_list_raw", data=np.array(ik_ordered, dtype=np.int32))
        scat.create_dataset("iq_list_raw", data=np.array(iq_ordered, dtype=np.int32))

        # ik_pos_map → array of size Nk giving EPW → parser index mapping
        ik_pos_arr = np.full((Nk,), -1, dtype=np.int32)
        for raw, pos in ik_pos_map.items():
            ik_pos_arr[pos] = raw
        scat.create_dataset("ik_pos_map", data=ik_pos_arr)

        # iq_pos_map → array of size Nq
        iq_pos_arr = np.full((Nq,), -1, dtype=np.int32)
        for raw, pos in iq_pos_map.items():
            iq_pos_arr[pos] = raw
        scat.create_dataset("iq_pos_map", data=iq_pos_arr)

        # --- Cartesian k and q for MC driver ---
        if recip_vectors is not None:
           kpts_cart = frac_to_cart(kpath_coords, recip_vectors)   # Nk x 3
           scat.create_dataset("kpts_cart", data=kpts_cart)

           if q_coords.shape[0] > 0:
               q_cart = frac_to_cart(q_coords, recip_vectors)      # Nq x 3
               scat.create_dataset("q_cart", data=q_cart)
        else:
            
            # fallback: store zeros so MC driver gracefully disables k+q
           scat.create_dataset("kpts_cart", data=np.zeros((Nk,3), dtype=float))
           scat.create_dataset("q_cart", data=np.zeros((Nq,3), dtype=float))


        summary_grp = fh.create_group("summary")
        summary_grp.attrs["n_rows"] = len(rows_array)
        summary_grp.attrs["nbands"] = nb
        summary_grp.attrs["Nk"] = Nk
        summary_grp.attrs["Nq"] = Nq

    if verbose:
        print(f"[ok] Wrote {out_h5}.")
        print(f"      valence_bands = {valence_bands}")
        print(f"      conduction_bands = {conduction_bands}")
        if ambiguous_bands:
            print(f"      ambiguous_bands = {ambiguous_bands} (bands crossing Fermi or uncertain)")

    # Write the human-readable parsed_log.txt next to out_h5
    summary = {"n_rows": len(rows_array), "nbands": nb, "Nk": Nk, "Nq": Nq}
    write_parsed_log(out_h5, parsed, qdat_coords, recip_vectors,
                 kpath_coords, q_coords, energies, velocities,
                 rows_array, band_list, summary,
                 ik_ordered, iq_ordered, ik_pos_map, iq_pos_map)

# ---------------------------
# CLI
# ---------------------------
def main():
    p = argparse.ArgumentParser(description="Parse EPW parsed into HDF5")
    p.add_argument("--quad", required=True, help="epw.out file")
    p.add_argument("--qdat", required=True, help="q.dat (k-path fractional coords)")
    p.add_argument("--out", default="epw_parsed.h5", help="output HDF5 filename")
    p.add_argument("--tol", type=float, default=1e-4, help="tolerance (eV) for Fermi compare")
    p.add_argument("--no-compress", action="store_true", help="disable gzip compression")
    args = p.parse_args()

    qcoords = read_qdat(args.qdat)
    if qcoords is None or qcoords.size == 0:
        print("[error] q.dat appears empty or unreadable.", file=sys.stderr)
        sys.exit(1)

    quad_dir = os.path.dirname(os.path.abspath(args.quad))
    scf_path = os.path.join(quad_dir, "scf.out")

    if os.path.exists(scf_path):
        print(f"[info] Reading reciprocal lattice from QE-scf")
        recip = read_recip_from_scf(scf_path, verbose=True)
    else:
        print("[warn] scf.out not found — reciprocal lattice vectors will not be stored.")
        recip = None

    parsed = parse_quadrupole_file(args.quad, verbose=True)

    build_and_write_h5(
        parsed,
        qcoords,
        args.out,
        recip_vectors=recip,
        tol=args.tol,
        compress=(not args.no_compress),
        verbose=True
    )

if __name__ == "__main__":
    main()

