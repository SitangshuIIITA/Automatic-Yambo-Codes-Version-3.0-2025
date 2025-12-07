#!/usr/bin/env python3
"""
Filtered EMC → EPW mapping script:
 - Supports ab_iq_mode, em_iq_mode with shape (N_iq, N_mode)
 - Computes GLOBAL dominant iq for each dominant mode
 - Prints normalized Abs_m, Em_m (0–1), and iq indices (0 instead of -1)
"""

import numpy as np
import h5py
import os, sys

NPZ_FILE  = "fb_emc_results_fast_mpi.npz"
EPW_FILE  = "epw_quadrupole_parsed.h5"
OUTDIR    = "population_out"
THRESHOLD = 0.05

# -------- column widths for aligned output --------
COL_W_IK  = 8     # width for ik_raw
COL_W_FLT = 14    # width for energy and normalized floats
COL_W_INT = 15    # width for integer columns (iq indices)

# ---------------- file checks ----------------
if not os.path.exists(NPZ_FILE):
    print("[ERROR] NPZ not found:", NPZ_FILE); sys.exit(1)
if not os.path.exists(EPW_FILE):
    print("[ERROR] EPW file not found:", EPW_FILE); sys.exit(1)

# ---------------- load NPZ --------------------
d = np.load(NPZ_FILE, allow_pickle=True)

E_fields        = d["E_fields"]
ik_list_emc     = np.array(d["ik_list"], dtype=int)
final_pops      = d["final_populations"]

em_per_mode_all = d["emission_per_mode"]
ab_per_mode_all = d["absorption_per_mode"]

em_iq_mode_all  = d["emission_per_iq_mode"]   # shape = (N_iq, N_mode)
ab_iq_mode_all  = d["absorption_per_iq_mode"]

used_bands = np.array(d["used_bands"], dtype=int)

# -------- detect active compact band -----------
total_by_band = None
for fp in final_pops:
    arr = np.array(fp, dtype=int)
    total_by_band = arr.sum(axis=0) if total_by_band is None else total_by_band + arr.sum(axis=0)

ib_active = int(np.argmax(total_by_band))
physical_ib = int(used_bands[ib_active])

# ---------------- load EPW ---------------------
with h5py.File(EPW_FILE, "r") as f:
    ik_list_raw = np.array(f["scattering"]["ik_list_raw"], dtype=int)
    energies    = np.array(f["bands"]["energies_eV"])

n_full = len(ik_list_raw)

def map_ik_to_epw_pos(ik_emc):
    m = np.where(ik_list_raw == ik_emc)[0]
    return int(m[0]) if m.size else None

os.makedirs(OUTDIR, exist_ok=True)

# ================= MAIN LOOP ==================
for i_field, E in enumerate(E_fields):

    print(f"\n[info] Field {i_field+1}/{len(E_fields)}  E={E:.3e}")

    # per-field arrays
    em_mode = np.array(em_per_mode_all[i_field], dtype=float)
    ab_mode = np.array(ab_per_mode_all[i_field], dtype=float)

    ab_iq_mode = np.array(ab_iq_mode_all[i_field], dtype=float)  # shape (N_iq, N_mode)
    em_iq_mode = np.array(em_iq_mode_all[i_field], dtype=float)

    # handle shapes robustly
    if ab_iq_mode.ndim != 2 or em_iq_mode.ndim != 2:
        raise RuntimeError("ab_iq_mode / em_iq_mode must be 2D arrays (N_iq, N_mode)")

    N_iq, N_mode = ab_iq_mode.shape

    # -------- find dominant modes (> THRESHOLD %) ----------
    abs_tot = ab_mode.sum(axis=0)
    em_tot  = em_mode.sum(axis=0)

    dom_modes = list(set(np.where(abs_tot >= THRESHOLD * abs_tot.sum())[0]) |
                     set(np.where(em_tot  >= THRESHOLD * em_tot.sum())[0]))

    # sort for consistent column order
    dom_modes = sorted(dom_modes)

    print("[info] Dominant modes:", dom_modes)

    # -------- compute GLOBAL iq for each mode -----
    iq_abs_global = {}
    iq_em_global  = {}

    for m in dom_modes:
        abs_vec = ab_iq_mode[:, m]
        em_vec  = em_iq_mode[:, m]

        # choose global iq index (0 if no events for that mode)
        iq_abs_global[m] = int(np.argmax(abs_vec)) if abs_vec.sum() > 0 else 0
        iq_em_global[m]  = int(np.argmax(em_vec))  if em_vec.sum() > 0 else 0

    print("[info] Global iq_abs:", iq_abs_global)
    print("[info] Global iq_em:",  iq_em_global)

    # -------- NORMALIZATION CONSTANTS ------------
    max_abs = ab_mode.max() if ab_mode.max() > 0 else 1.0
    max_em  = em_mode.max() if em_mode.max() > 0 else 1.0

    # -------- output file ------------------------
    outname = os.path.join(OUTDIR, f"E_{E:.3e}.txt")
    with open(outname, "w") as fh:

        fh.write(f"# This is phonon mode (iq) resolved final steady state scattering counts for E={E:.3e} V/m\n")
        fh.write(f"# Dominant modes: {dom_modes}\n")
        fh.write(f"# Careful: Modes are in 0 indexed numbers.\n")
        fh.write("# Columns: Normalized counts for Absorbtion & Emission, iq (Abs), iq (Em) (iq are global)\n")
        fh.write(f"# compact_ib={ib_active}, physical_ib={physical_ib}\n")

        # ---- WRITE ALIGNED HEADER ----
        header_fields = []

        # ik_raw and Enk_eV
        header_fields.append(f"{'#ik':>{COL_W_IK}}")
        header_fields.append(f"{'Enk (eV)':>{COL_W_FLT}}")

        # For each dominant mode, add four aligned columns
        for m in dom_modes:
            header_fields.append(f"{('Abs_mmode:'+str(m)):>{COL_W_INT}}")
            header_fields.append(f"{('Em_mode:'+str(m)):>{COL_W_INT}}")
            header_fields.append(f"{('iq_abs_mode:'+str(m)):>{COL_W_INT}}")
            header_fields.append(f"{('iq_em_mode:'+str(m)):>{COL_W_INT}}")

        fh.write(" ".join(header_fields) + "\n\n")

        # ---- loop over EPW grid ----
        for idx_epw in range(n_full):

            ik_raw = int(ik_list_raw[idx_epw])
            Enk = float(energies[physical_ib - 1, idx_epw])

            pos = np.where(ik_list_emc == ik_raw)[0]

            row_values = []

            if pos.size:
                ikpos = int(pos[0])
                for m in dom_modes:

                    # normalized floats (0..1)
                    abs_m = ab_mode[ikpos, m] / max_abs
                    em_m  = em_mode[ikpos, m] / max_em

                    # global iq indices (0 if none)
                    iqA = iq_abs_global[m]
                    iqE = iq_em_global[m]

                    # ensure -1 replaced with 0 (already ensured above)
                    if iqA < 0: iqA = 0
                    if iqE < 0: iqE = 0

                    row_values += [
                        abs_m,
                        em_m,
                        iqA,
                        iqE
                    ]

            else:
                # replace -1 by 0 and zero-normalized counts
                for m in dom_modes:
                    row_values += [0.0, 0.0, 0, 0]

            # ---- FORMAT A DATA ROW ----
            row_items = []

            row_items.append(f"{ik_raw:{COL_W_IK}d}")
            row_items.append(f"{Enk:{COL_W_FLT}.6f}")

            for v in row_values:
                if isinstance(v, float):
                    row_items.append(f"{v:{COL_W_FLT}.6f}")
                else:
                    row_items.append(f"{v:{COL_W_INT}d}")

            fh.write(" ".join(row_items) + "\n")

    print("[write]", outname)

print("\n[done] Mapping complete.\n")

