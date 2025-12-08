#!/usr/bin/env python3
"""
fb_emc_fast.py

Full-band EMC engine (speed-optimized defaults)

Speed fixes included:
 - MAX_STEPS cap per particle to prevent enormous step counts
 - Interpolation disabled by default
 - Fast nearest-ik lookup (vectorized) for small ik grids
 - Default T_sim reduced for quick tests

Behavior:
 - Uses scattering_events.h5 table (expects same structure you provided).
 - Uses correct dk_cart -> fractional conversion via B_reciprocal.
 - Shifts DFT energies so CBM -> 0 (energies in meV).

This version adds v_x(t) snapshot recording and writes a text file
vx_time_E_<field>.txt for each field. Otherwise behavior is unchanged.
"""

from __future__ import annotations
import h5py
import numpy as np
import math
from dataclasses import dataclass
from typing import Tuple, List
from scipy.spatial import cKDTree
import time
import argparse
import matplotlib.pyplot as plt
import os

# Physical constants
q = 1.602176634e-19
hbar = 1.054571817e-34
kB = 1.380649e-23
eV_to_J = 1.602176634e-19

def detect_steady_state(snapshot_times, v_snapshot, rel_tol=0.01, bins_required=10):
    """
    Detect the time at which steady state is reached.

    Steady state = when the relative change in v_x between consecutive bins
    stays below rel_tol for bins_required consecutive bins.

    Returns:
        t_ss (float or None): time where steady state starts.
    """
    count = 0
    for i in range(1, len(v_snapshot)):
        v0 = v_snapshot[i-1]
        v1 = v_snapshot[i]
        if v0 == 0:
            continue
        rel_change = abs(v1 - v0) / max(abs(v0), 1e-30)
        if rel_change < rel_tol:
            count += 1
            if count >= bins_required:
                return snapshot_times[i - bins_required]
        else:
            count = 0
    return None


# ---------- Data container ----------
@dataclass
class EventsTable:
    k_frac: np.ndarray
    kq_frac: np.ndarray
    ik: np.ndarray
    ib: np.ndarray
    jb: np.ndarray
    imode: np.ndarray
    Gamma_w: np.ndarray
    Gamma_em: np.ndarray
    Gamma_abs: np.ndarray
    enk_meV: np.ndarray
    enkq_meV: np.ndarray
    v_group: np.ndarray
    v_scattered: np.ndarray
    B_reciprocal: np.ndarray
    enk_above_CBM_meV: np.ndarray = None
    enkq_above_CBM_meV: np.ndarray = None

# ---------- Load HDF5 ----------
def load_h5(path: str) -> EventsTable:
    with h5py.File(path, "r") as f:
        g = f['scattering_events']
        k_frac = np.array(g['k_frac'])
        kq_frac = np.array(g['kq_frac'])
        ik = np.array(g['ik']).astype(int)
        ib = np.array(g['ib']).astype(int)
        jb = np.array(g['jb']).astype(int)
        imode = np.array(g['imode']).astype(int) if 'imode' in g else np.zeros_like(ik)
        Gamma_w = np.array(g['Gamma_w']).astype(float)
        Gamma_em = np.array(g['Gamma_em_w']).astype(float) if 'Gamma_em_w' in g else np.zeros_like(Gamma_w)
        Gamma_abs = np.array(g['Gamma_abs_w']).astype(float) if 'Gamma_abs_w' in g else np.zeros_like(Gamma_w)
        enk_meV = np.array(g['enk_meV']).astype(float)
        enkq_meV = np.array(g['enkq_meV']).astype(float)
        v_group = np.array(g['v_group_cart']).astype(float)
        v_scattered = np.array(g['v_scattered_cart']).astype(float)
        B_reciprocal = np.array(g['B_reciprocal'])
    # shift energies so CBM -> 0
    E_CBM_meV = np.min(enk_meV)
    enk_above = enk_meV - E_CBM_meV
    enkq_above = enkq_meV - E_CBM_meV
    return EventsTable(k_frac, kq_frac, ik, ib, jb, imode, Gamma_w, Gamma_em, Gamma_abs,
                       enk_meV, enkq_meV, v_group, v_scattered, B_reciprocal,
                       enk_above, enkq_above)

# ---------- Build state index ----------
def build_state_index(events: EventsTable):
    state_to_evs = {}
    for j in range(len(events.ik)):
        key = (int(events.ik[j]), int(events.ib[j]))
        state_to_evs.setdefault(key, []).append(j)
    state_Gamma_tot = {k: float(np.sum(events.Gamma_w[v])) for k,v in state_to_evs.items()}
    state_cum = {}
    for k, evs in state_to_evs.items():
        rates = events.Gamma_w[evs]
        tot = rates.sum()
        if tot > 0.0:
            probs = rates / tot
            cum = np.concatenate(([0.0], np.cumsum(probs)))
            state_cum[k] = (np.array(evs, dtype=int), cum)
        else:
            state_cum[k] = (np.array(evs, dtype=int), np.array([0.0]))
    # convert lists to arrays in state_to_evs
    for k in list(state_to_evs.keys()):
        state_to_evs[k] = np.array(state_to_evs[k], dtype=int)
    return state_to_evs, state_Gamma_tot, state_cum

# ---------- Build ik grid and KD-tree (cartesian) ----------
def build_kgrid_cart(events: EventsTable):
    unique_ik = np.unique(events.ik)
    ik_to_rep = {}
    for ik in unique_ik:
        idx = np.where(events.ik == ik)[0][0]
        ik_to_rep[int(ik)] = idx
    ik_list = np.array(list(ik_to_rep.keys()), dtype=int)
    k_frac_grid = np.array([events.k_frac[ik_to_rep[ik]] for ik in ik_list])  # fractional
    # prepare Brec_cart (fractional -> cart). Heuristic: if entries small assume 1/Ang -> convert to 1/m
    Braw = np.array(events.B_reciprocal)
    B = Braw.reshape(3,3) if np.asarray(Braw).size == 9 else np.array(Braw)
    if np.max(np.abs(B)) < 100:
        Brec_cart = B * 1e10
    else:
        Brec_cart = B.copy()
    # cart k points
    k_cart = (Brec_cart @ k_frac_grid.T).T
    # KD-tree
    tree = cKDTree(k_cart)
    return ik_list, ik_to_rep, k_frac_grid, k_cart, tree, Brec_cart

# ---------- Particle dataclass ----------
@dataclass
class Particle:
    ik: int
    ib: int
    k_frac: np.ndarray
    v: np.ndarray
    E_meV: float
    t: float = 0.0
    time_vx: float = 0.0
    time_E: float = 0.0
    scat_count: int = 0

# ---------- Init particles ----------
def init_particles(N:int, events:EventsTable, ik_list, ik_to_rep, T=300.0):
    reps = np.array([ik_to_rep[ik] for ik in ik_list])
    E_J = events.enk_meV[reps] * 1e-3 * eV_to_J
    w = np.exp(-E_J / (kB * T))
    w /= w.sum()
    rng = np.random.default_rng()
    parts = []
    for _ in range(N):
        pick = rng.choice(len(ik_list), p=w)
        ik0 = int(ik_list[pick]); rep = ik_to_rep[ik0]
        parts.append(Particle(ik0, int(events.ib[rep]), events.k_frac[rep].copy(),
                              events.v_group[rep].copy(), float(events.enk_above_CBM_meV[rep])))
    return parts

# ---------- Fast nearest-ik (vectorized) ----------
def nearest_ik_from_frac(k_frac_query: np.ndarray, ik_list: np.ndarray, k_frac_grid: np.ndarray) -> int:
    # works well when ik_count is small-to-moderate
    diff = (k_frac_grid - k_frac_query + 0.5) % 1.0 - 0.5
    d2 = np.sum(diff*diff, axis=1)
    arg = int(np.argmin(d2))
    return int(ik_list[arg])

# ---------- Core EMC runner (optimized) ----------
def run_field(events:EventsTable, state_Gamma_tot, state_cum, ik_list, ik_to_rep,
              k_frac_grid, k_cart_grid, tree, Brec_cart,
              particles:List[Particle], E_field:float, T_sim:float,
              interpolate:bool=False, MAX_STEPS:int=5000, use_kdtree_threshold:int=2000):
    """
    Run EMC for one field.
      - interpolate: whether to interpolate v(k). Default False (fast).
      - MAX_STEPS: hard cap on per-particle iteration count.
      - use_kdtree_threshold: if number of ik > threshold, use KD-tree; otherwise vectorized lookup.
    """

    invB = np.linalg.inv(Brec_cart)
    rng = np.random.default_rng()
    ik_count = len(ik_list)
    use_kdtree = ik_count > use_kdtree_threshold   # KD-tree only if many ik

    # --------- v(t) snapshot arrays (added, non-invasive) -----------
    snapshots = 200
    snapshot_times = np.linspace(0.0, T_sim, snapshots)
    bin_time  = np.zeros(snapshots)
    bin_vtime = np.zeros(snapshots)

            # ---- Population accumulation: (time, ik, band) ----
    N_snap = snapshots
    N_ik = len(ik_list)
    N_band = np.max(events.ib) + 1
    pop_bins = np.zeros((N_snap, N_ik, N_band), dtype=int)

    # map ik → position (0..N_ik-1)
    ik_to_pos = {ik: idx for idx, ik in enumerate(ik_list)}

    def _accumulate_population(t0, dt, ik_pos, ib):
        
        """Record population in all snapshot bins overlapping [t0, t0+dt]."""
        if dt <= 0:
            return
        t1 = t0 + dt
        f0 = max(0.0, min(1.0, t0 / T_sim))
        f1 = max(0.0, min(1.0, t1 / T_sim))
        i0 = int(f0 * (N_snap - 1))
        i1 = int(f1 * (N_snap - 1))

        if i0 == i1:
           pop_bins[i0, ik_pos, ib] += 1
        else:
            pop_bins[i0:i1+1, ik_pos, ib] += 1

    def _accumulate_snapshot(t0, dt, vx):
        if dt <= 0:
            return
        t1 = t0 + dt
        f0 = max(0.0, min(1.0, t0 / T_sim))
        f1 = max(0.0, min(1.0, t1 / T_sim))
        i0 = int(f0 * (snapshots - 1))
        i1 = int(f1 * (snapshots - 1))
        if i0 == i1:
            bin_time[i0]  += dt
            bin_vtime[i0] += vx * dt
        else:
            t_edge = snapshot_times[i0+1] if i0+1 < snapshots else T_sim
            dt0 = min(dt, t_edge - t0)
            bin_time[i0]  += dt0
            bin_vtime[i0] += vx * dt0
            for ib in range(i0+1, i1):
                tb0 = snapshot_times[ib]
                tb1 = snapshot_times[ib+1] if ib+1 < snapshots else T_sim
                dtb = tb1 - tb0
                bin_time[ib]  += dtb
                bin_vtime[ib] += vx * dtb
            tb0 = snapshot_times[i1]
            dt_last = max(0.0, t1 - tb0)
            bin_time[i1]  += dt_last
            bin_vtime[i1] += vx * dt_last

    # Reset particle statistics
    for p in particles:
        p.t = 0.0
        p.time_vx = 0.0
        p.time_E = 0.0
        p.scat_count = 0

    # ---------------- Main particle loop ----------------
    for p in particles:
        steps = 0

        while p.t < T_sim and steps < MAX_STEPS:
            steps += 1

            key = (int(p.ik), int(p.ib))
            Gtot = state_Gamma_tot.get(key, 0.0)

            # No scattering available → free fly to end
            if Gtot <= 0.0:
                dt = T_sim - p.t
                # accumulate snapshot
                _accumulate_snapshot(p.t, dt, p.v[0])
                p.time_vx += p.v[0] * dt
                p.time_E  += p.E_meV * dt
                # population logging
                ik_pos = ik_to_pos[p.ik]
                _accumulate_population(p.t, dt, ik_pos, p.ib)
                
                p.t += dt
                
                break

            # Sample free-flight
            dt = -math.log(rng.random()) / Gtot
            if p.t + dt > T_sim:
                dt = T_sim - p.t

            # ---------- FREE FLIGHT k-update ----------
            dk_cart = np.array([-(q * E_field / hbar) * dt, 0.0, 0.0])   # 1/m
            dk_frac = invB.dot(dk_cart)
            p.k_frac = (p.k_frac + dk_frac) % 1.0

            # ---------- nearest ik selection ----------
            if use_kdtree:
                k_cart_point = Brec_cart @ p.k_frac
                _, idx = tree.query(k_cart_point)
                p.ik = int(ik_list[idx])
                rep = ik_to_rep[p.ik]
            else:
                p.ik = nearest_ik_from_frac(p.k_frac, ik_list, k_frac_grid)
                rep  = ik_to_rep[p.ik]

            # ---------- update velocity ----------
            if not interpolate:
                # FAST MODE
                p.v = events.v_group[rep].copy()

            else:
                # ------- INTERPOLATION MODE (corrected) -------
                k_cart_point = Brec_cart @ p.k_frac

                # query nearest 6 points (may return scalar)
                dists, idxs = tree.query(k_cart_point, k=6)

                # ensure arrays
                if np.isscalar(idxs):
                    idxs = np.array([idxs])
                    dists = np.array([dists])
                else:
                    idxs = np.array(idxs)
                    dists = np.array(dists)

                # corresponding ik values
                neighbors_ik = ik_list[idxs]

                # representative indices
                rep_indices = [ik_to_rep[int(k)] for k in neighbors_ik]
                v_neighbors = events.v_group[rep_indices]

                # weights
                d2 = dists ** 2
                if len(d2) > 1 and np.any(d2[1:] > 0):
                    denom = np.mean(d2[1:]) + 1e-30
                else:
                    denom = d2[0] + 1e-30

                weights = np.exp(-d2 / denom)
                weights /= weights.sum()

                # interpolated velocity
                p.v = np.average(v_neighbors, axis=0, weights=weights)

            # update energy after free flight
            p.E_meV = float(events.enk_above_CBM_meV[rep])

            # integrate observables (and accumulate vx(t))
            _accumulate_snapshot(p.t, dt, p.v[0])
            p.time_vx += p.v[0] * dt
            p.time_E  += p.E_meV * dt
            p.t       += dt

            if p.t >= T_sim:
                break

            # ---------- SCATTERING ----------
            evs, cum = state_cum[key]
            if len(evs) == 0:
                continue

            u = rng.random()
            idx_in = max(0, np.searchsorted(cum, u) - 1)
            ev = int(evs[idx_in])

            # new k_frac from scattering
            p.k_frac = events.kq_frac[ev].copy()

            # nearest ik for final state
            if use_kdtree:
                k_cart_point = Brec_cart @ p.k_frac
                _, idx2 = tree.query(k_cart_point)
                p.ik = int(ik_list[idx2])
                rep2 = ik_to_rep[p.ik]
            else:
                p.ik = nearest_ik_from_frac(p.k_frac, ik_list, k_frac_grid)
                rep2 = ik_to_rep[p.ik]

            # band + velocity + energy update
            p.ib     = int(events.jb[ev])
            p.v      = events.v_scattered[ev].copy()
            p.E_meV  = float(events.enkq_above_CBM_meV[ev])
            p.scat_count += 1

        # END particle

    # -------------- observables ---------------
    Tref = 300.0  # the simulation temperature
    kBT_meV = (kB * Tref) * 1000 / (1.602176634e-19)   # convert J to meV

    total_time = len(particles) * T_sim
    drift_raw  = sum(pp.time_vx for pp in particles) / total_time
    meanE_raw  = sum(pp.time_E  for pp in particles) / total_time

    # total kinetic energy = DFT band energy + thermal kinetic energy
    meanE_total = meanE_raw #+ kBT_meV

    mean_s = np.mean([pp.scat_count for pp in particles])

    # finalize v(t) snapshot
    v_snapshot = np.zeros_like(bin_time)
    for i in range(len(bin_time)):
        if bin_time[i] > 0:
            v_snapshot[i] = bin_vtime[i] / bin_time[i]
        else:
            v_snapshot[i] = 0.0

    #return drift_raw, meanE_total, mean_s, snapshot_times, v_snapshot
    return drift_raw, meanE_total, mean_s, snapshot_times, v_snapshot, pop_bins

# ---------- Minimal main ----------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fast Full-band EMC runner')
    parser.add_argument('--h5', type=str, default='scattering_events.h5', help='HDF5 events file')
    parser.add_argument('--npar', type=int, default=1000, help='Number of particles (ensemble)')
    parser.add_argument('--tsim', type=float, default=5e-12, help='Simulation time per particle (s) (reduced default for speed)')
    parser.add_argument('--Emin', type=float, default=1e3, help='Min field V/m')
    parser.add_argument('--Emax', type=float, default=1e8, help='Max field V/m')
    parser.add_argument('--nE', type=int, default=6, help='Number of field points')
    parser.add_argument('--interp', action='store_true', help='Enable v(k) interpolation (slow)')
    parser.add_argument('--maxsteps', type=int, default=20000, help='Max steps per particle (cap)')
    args = parser.parse_args()

    t0 = time.time()
    print("Loading HDF5 ...")
    events = load_h5(args.h5)
    print("Building state index ...")
    state_to_evs, state_Gamma_tot, state_cum = build_state_index(events)
    print("Building k-grid and KD-tree ...")
    ik_list, ik_to_rep, k_frac_grid, k_cart_grid, tree, Brec_cart = build_kgrid_cart(events)
    print(f"Unique ik count = {len(ik_list)}")
    print("Initializing particles ...")
    particles = init_particles(args.npar, events, ik_list, ik_to_rep, T=300.0)

    E_fields = np.logspace(math.log10(args.Emin), math.log10(args.Emax), args.nE)
    drifts = []
    meanEs = []
    meanScats = []

    for E in E_fields:
        t1 = time.time()

        base_out = "Drift_out"
        os.makedirs(base_out, exist_ok=True)

        field_folder = os.path.join(base_out, f"E_{E:.3e}")
        os.makedirs(field_folder, exist_ok=True)
    	# ----------------------------------------

        dv, me, ms, snap_t, snap_v, pop_bins = run_field(
            events, state_Gamma_tot, state_cum,
            ik_list, ik_to_rep,
            k_frac_grid, k_cart_grid,
            tree, Brec_cart,
            particles, E, args.tsim,
            interpolate=args.interp, MAX_STEPS=args.maxsteps
        )

        # store snapshot_times once (important)
        if "snapshot_times" not in locals():
            snapshot_times = snap_t

        dt = time.time() - t1
        print(f"E={E:.3e} V/m -> drift {dv:.3e} m/s, meanE {me:.3f} meV, avg_scats {ms:.2f} (time {dt:.2f}s)")

        drifts.append(dv)
        meanEs.append(me)
        meanScats.append(ms)

        # store population grids
        if "all_pops" not in locals():
            all_pops = []
        all_pops.append(pop_bins)

        # write vx(t)
        fname = os.path.join(field_folder, "vx_time.txt")
        try:
            with open(fname, "w") as f:
                f.write("# time(s)    vx(m/s)\n")
                for tval, vxval in zip(snap_t, snap_v):
                    f.write(f"{tval: .6e}   {vxval: .6e}\n")
            print(f"Saved vx(t) -> {fname}")

            # steady state
            t_ss = detect_steady_state(snap_t, snap_v)
            # Output path
            ssfile = os.path.join(field_folder, "steady_state.txt")
            try:
               with open(ssfile, "w") as fss:
                   fss.write(f"# Steady state info for E={E:.3e} V/m\n")

                   if t_ss is not None:
            	       # ------------------- steady state reached -------------------
                       vx_ss = snap_v[np.argmin(abs(snap_t - t_ss))]
                       fss.write(f"steady state time in s    : {t_ss:.6e}\n")
                       fss.write(f"steady state vx in ms^-1  : {vx_ss:.6e}\n")
                       fss.write("steady state status:       : reached\n")
                       print(f"    Steady-state reached at t = {t_ss:.3e} s")

                   else:
                       # ----------------- steady state not reached -----------------
                       fss.write("steady state_time (s)      :  not_reached\n")
                       fss.write("steady state vx ms^-1      :  not_reached\n")
                       fss.write("steady state status        :  not_reached\n")
                       print("    WARNING: steady state not detected.")

                   # always include these fields
                       fss.write(f"drift velocity in ms^-1   :    {dv:.6e}\n")
                       fss.write(f"mean energy in meV        :    {me:.6f}\n")
                       fss.write(f"avgerage scatter count    :    {ms:.2f}\n")

               print(f"    Saved steady-state info -> {ssfile}")

            except Exception as exc:
                print("Failed to write steady-state file:", exc)

        except Exception as exc:
            print("Failed to write vx(t) file:", exc)        
                
    # Save master NPZ including new population data
    np.savez(
        'fb_emc_results_fast.npz',
        E_fields=E_fields,
        drifts=drifts,
        meanEs=meanEs,
        meanScats=meanScats,
        populations=np.array(all_pops, dtype=object),
        snapshot_times=snapshot_times
    )

    # simple plots
    plt.figure(figsize=(10,4))
    plt.subplot(1,2,1)
    plt.semilogx(E_fields, drifts, marker='o')
    plt.xlabel('Electric field |E| (V/m)')
    plt.ylabel('Drift velocity <v_x> (m/s)')
    plt.grid(True)

    plt.subplot(1,2,2)
    plt.semilogx(E_fields, meanEs, marker='o', color='C1')
    plt.xlabel('Electric field |E| (V/m)')
    plt.ylabel('Mean energy above CBM (meV)')
    plt.grid(True)

    plt.tight_layout()
    plt.savefig('fb_emc_results_fast.png')
    print("Saved fb_emc_results_fast.npz and fb_emc_results_fast.png")
    print("Total runtime: %.2f s" % (time.time() - t0))

