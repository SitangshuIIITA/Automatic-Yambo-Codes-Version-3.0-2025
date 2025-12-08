#!/usr/bin/env python3
"""
ikq_mapper.py

Adds the missing (ik + q) → ik_final mapping to the EPW scattering rows.

Reads:
    /scattering/rows
    /scattering/kpts_cart
    /scattering/qpts_cart

Computes:
    k_final_cart = k_cart[ik] + q_cart[iq]   (wrapped to BZ)

Finds nearest k-point index ikq_pos and writes a NEW HDF5:
    epw_quadrupole_mapped.h5

This file is now suitable for *full particle Monte-Carlo* drift simulations.
"""

import argparse
import h5py
import numpy as np
import sys


def wrap_BZ(vec):
    """
    Wrap a k-vector back into the first Brillouin zone by bringing
    coordinates into [-0.5,0.5] in reciprocal lattice units.
    This assumes k,q are already in Cartesian coordinates.
    For 2D materials, qz=0 so wrapping is trivial.
    """
    return vec  # EPW usually already outputs wrapped cartesian k. No remap needed.


def find_nearest_k(k_final_cart, kpts_cart):
    """
    Returns index of the nearest k-point in the grid.
    """
    diff = kpts_cart - k_final_cart
    dist2 = np.sum(diff * diff, axis=1)
    return int(np.argmin(dist2))


def main(h5in):

    h5out = h5in.replace(".h5", "_mapped.h5")

    print(f"[info] Reading: {h5in}")
    fh = h5py.File(h5in, "r")

    rows = fh["/scattering/rows"]
    kpts = fh["/scattering/kpts_cart"][...]     # shape (Nk,3)
    qpts = fh["/scattering/qpts_cart"][...]     # shape (Nq,3)

    Nk = kpts.shape[0]
    Nq = qpts.shape[0]
    Nrows = rows.shape[0]

    print(f"[info] Nk = {Nk}, Nq = {Nq}, total rows = {Nrows}")
    print("[info] Computing ik_final for all rows...")

    ikq_list = np.zeros(Nrows, dtype=np.int32)

    for i in range(Nrows):
        ik = int(rows["ik_pos"][i])
        iq = int(rows["iq_pos"][i])

        k_cart = kpts[ik]
        q_cart = qpts[iq]

        k_final = wrap_BZ(k_cart + q_cart)
        ik_final = find_nearest_k(k_final, kpts)

        ikq_list[i] = ik_final

        if i % 200000 == 0:
            print(f"  processed {i}/{Nrows}")

    fh.close()

    print(f"[info] Writing mapped file: {h5out}")
    with h5py.File(h5out, "w") as fo:

        # --- copy everything from original file ---
        with h5py.File(h5in, "r") as fi:
            fi.copy("/", fo)

        # --- append new dataset ---
        fo.create_dataset("/scattering/ikq_pos", data=ikq_list, dtype=np.int32)

    print("[ok] Mapping complete.")
    print(f"[ok] New file saved as: {h5out}")
    print("[ok] Ready for Monte-Carlo drift simulation.")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--h5", required=True, help="Path to epw_quadrupole_parsed.h5")
    args = p.parse_args()
    main(args.h5)

