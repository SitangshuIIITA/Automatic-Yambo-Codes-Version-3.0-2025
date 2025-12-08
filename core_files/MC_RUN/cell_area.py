#!/usr/bin/env python3
"""
compute_scattering_stats.py — cleaned & consolidated version

- Computes per-event Golden-rule scattering rates from parsed EPW HDF5,
  aggregates them, enforces detailed balance (per-key) and then redistributes
  corrected aggregated Abs/Em rates back to events proportionally.
- Writes:
   - scattering_events.h5   (HDF5 with per-event Gamma_abs_w, Gamma_em_w, Gamma_w, velocities when available)
   - scattering_events.txt  (plain text one-to-one)
   - scattering_rates.txt   (aggregated per-(ik,ib,jb,type))
   - report_scatterings.txt (summary)
   - Also writes A_cell and A_BZ into the EPW-parsed HDF5 file
"""

from __future__ import annotations
import argparse, os, sys, time, math
from collections import defaultdict
import h5py, re
import numpy as np

BOHR_TO_M = 0.529177210903e-10  # meters


def parse_scf_out(path):
    """
    Parse QE scf.out file in units:
      - a(i) in units of alat
      - b(i) in units of (2π/alat)
    """

    a_vecs = []
    b_vecs = []
    alat_bohr = None

    num_regex = re.compile(r"[-+]?\d*\.\d+(?:[Ee][-+]?\d+)?|[-+]?\d+")

    with open(path) as f:
        lines = f.readlines()

    for i, line in enumerate(lines):

        # Extract celldm(1) => alat (Bohr)
        if "celldm(1)" in line:
            nums = num_regex.findall(line)
            alat_bohr = float(nums[0])   # celldm(1)
            # alat in meters
            alat_m = alat_bohr * BOHR_TO_M

        # Crystal axes in units of alat
        if "crystal axes:" in line.lower():
            for j in range(1, 4):
                nums = num_regex.findall(lines[i + j])
                a_vecs.append([float(n) for n in nums[-3:]])

        # Reciprocal axes in units 2π/alat
        if "reciprocal axes:" in line.lower():
            for j in range(1, 4):
                nums = num_regex.findall(lines[i + j])
                b_vecs.append([float(n) for n in nums[-3:]])

    return np.array(a_vecs), np.array(b_vecs), alat_m


def compute_2D_areas(a_alat, b_2pi_alat, alat_m):
    """
    Inputs:
      a_alat      : real-space lattice vectors in units of 'alat'
      b_2pi_alat  : reciprocal lattice vectors in units 2π/alat
      alat_m      : alat in meters

    Returns:
      A_real  : 2D real-space area (m^2)
      A_bz    : 2D BZ area (1/m^2)
    """

    # Convert real-space a-vectors → meters
    a = a_alat * alat_m

    # Convert reciprocal b-vectors → 1/m
    factor = (2*np.pi) / alat_m
    b = b_2pi_alat * factor

    # Only a1, a2 for 2D layer
    a1 = a[0]
    a2 = a[1]
    A_real = np.linalg.norm(np.cross(a1, a2))

    b1 = b[0]
    b2 = b[1]
    A_bz = np.linalg.norm(np.cross(b1, b2))

    return A_real, A_bz


def write_areas_to_h5(h5_path, A_cell, A_BZ):
    """
    Writes A_cell and A_BZ into:
       /meta/areas/A_cell
       /meta/areas/A_BZ
    """

    with h5py.File(h5_path, "a") as f:
        meta = f.require_group("meta")
        areas = meta.require_group("areas")

        # Delete existing datasets if present
        for key in ["A_cell", "A_BZ"]:
            if key in areas:
                del areas[key]

        areas.create_dataset("A_cell", data=A_cell)
        areas.create_dataset("A_BZ", data=A_BZ)

    print(f"✓ Wrote A_cell and A_BZ into {h5_path} under /meta/areas/")


if __name__ == "__main__":

    scf_path = "scf.out"
    h5_path  = "epw_quadrupole_parsed.h5"   # <-- MODIFY IF NEEDED

    # Parse SCF
    a_alat, b_2pi_alat, alat_m = parse_scf_out(scf_path)

    # Compute areas
    A_real, A_bz = compute_2D_areas(a_alat, b_2pi_alat, alat_m)

    # Display
    print("Real-space a-vectors (units of alat):")
    print(a_alat)
    print("\nReciprocal b-vectors (units of 2π/alat):")
    print(b_2pi_alat)
    print("===================================================")
    print(f"alat = {alat_m:.6e} m")
    print(f"A_cell (real space) = {A_real:.6e} m^2")
    print(f"A_BZ   (reciprocal) = {A_bz:.6e} 1/m^2")
    print("Check: A_cell * A_BZ =", A_real * A_bz)
    print("Should equal (2π)^2 =", (2*np.pi)**2)

    # Write into HDF5
    write_areas_to_h5(h5_path, A_real, A_bz)

