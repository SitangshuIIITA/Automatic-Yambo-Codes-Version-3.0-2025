#!/usr/bin/env python3
"""
fb_emc_fast_mpi.py

MPI-enabled Full-band EMC runner (parallel over particles)

This file is a minimally invasive MPI port of your original fb_emc_fast.py:
 - Uses mpi4py to distribute particles across ranks
 - Rank 0 loads the HDF5 file and broadcasts the EventsTable
 - Each rank runs run_field() on its subset of particles
 - Reduction gathers numeric sums/arrays and rank 0 assembles final outputs

Physics and numerical logic are unchanged.
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
import sys

# MPI
try:
    from mpi4py import MPI
except Exception as exc:
    print("ERROR: mpi4py is required for this script. Install with 'pip install mpi4py'.")
    raise

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Physical constants
q = 1.602176634e-19
hbar = 1.054571817e-34
kB = 1.380649e-23
eV_to_J = 1.602176634e-19

def detect_steady_state(snapshot_times, v_snapshot, rel_tol=0.03, bins_required=25):
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
    Braw = np.array(events.B_reciprocal)
    B = Braw.reshape(3,3) if np.asarray(Braw).size == 9 else np.array(Braw)
    if np.max(np.abs(B)) < 100:
        Brec_cart = B * 1e10
    else:
        Brec_cart = B.copy()
    k_cart = (Brec_cart @ k_frac_grid.T).T
    tree = cKDTree(k_cart)
    # ensure contiguous arrays and float64 for speed
    k_frac_grid = np.ascontiguousarray(k_frac_grid, dtype=np.float64)
    ik_list = np.ascontiguousarray(ik_list, dtype=np.int64)
    return ik_list, ik_to_rep, k_frac_grid, k_cart, tree, Brec_cart

# ---------- Fast nearest-ik (vectorized, replacement) ----------
def nearest_ik_from_frac(k_frac_query: np.ndarray, ik_list: np.ndarray, k_frac_grid: np.ndarray) -> int:
    kq = np.asarray(k_frac_query, dtype=np.float64).reshape(1, 3)
    kgrid = k_frac_grid  # assumed contiguous float64
    diff = kgrid - kq
    diff -= np.round(diff)
    d2 = np.einsum('ij,ij->i', diff, diff)
    arg = int(np.argmin(d2))
    return int(ik_list[arg])

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

# ---------- Core EMC runner (optimized) ----------
def run_field(events:EventsTable, state_Gamma_tot, state_cum, ik_list, ik_to_rep,
              k_frac_grid, k_cart_grid, tree, Brec_cart,
              particles:List[Particle], E_field:float, T_sim:float,
              interpolate:bool=False, MAX_STEPS:int=5000, use_kdtree_threshold:int=2000):
    """
    Runs EMC for the local subset of particles.
    Returns:
      - sum_time_vx: sum over particles of (integral v_x dt)
      - sum_time_E:  sum over particles of (integral E_meV dt)
      - sum_scats:   sum of scatter counts over particles
      - local_total_time: len(particles) * T_sim
      - mean_s_local: average scatter per particle (for convenience)
      - snapshot_times (shared)
      - bin_time: per-bin accumulated time (for reductions)
      - bin_vtime: per-bin accumulated v*dt (for reductions)
      - pop_bins: local population bin counts (for reductions)
    """
    
    # --- NEW: emission / absorption counters ---
    emission_count = 0
    absorption_count = 0

     
    invB = np.linalg.inv(Brec_cart)
    rng = np.random.default_rng()
    ik_count = len(ik_list)
    use_kdtree = ik_count > use_kdtree_threshold

    # snapshot config
    snapshots = 200
    snapshot_times = np.linspace(0.0, T_sim, snapshots)
    bin_time  = np.zeros(snapshots, dtype=np.float64)
    bin_vtime = np.zeros(snapshots, dtype=np.float64)

    N_snap = snapshots
    N_ik = len(ik_list)
    N_band = int(np.max(events.ib) + 1)

    pop_bins = np.zeros((N_snap, N_ik, N_band), dtype=np.int64)

    ik_to_pos = {ik: idx for idx, ik in enumerate(ik_list)}

    def _accumulate_population(t0, dt, ik_pos, ib):
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
        i0 = int(f0 * (N_snap - 1))
        i1 = int(f1 * (N_snap - 1))
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

    for p in particles:
        steps = 0
        while p.t < T_sim and steps < MAX_STEPS:
            steps += 1
            key = (int(p.ik), int(p.ib))
            Gtot = state_Gamma_tot.get(key, 0.0)
            if Gtot <= 0.0:
                dt = T_sim - p.t
                _accumulate_snapshot(p.t, dt, p.v[0])
                p.time_vx += p.v[0] * dt
                p.time_E  += p.E_meV * dt
                ik_pos = ik_to_pos[p.ik]
                _accumulate_population(p.t, dt, ik_pos, p.ib)
                p.t += dt
                break
            dt = -math.log(rng.random()) / Gtot
            if p.t + dt > T_sim:
                dt = T_sim - p.t
            # free flight k-update: precompute can be used by caller per-field if desired
            dk_cart = np.array([-(q * E_field / hbar) * dt, 0.0, 0.0])
            dk_frac = invB.dot(dk_cart)
            p.k_frac = (p.k_frac + dk_frac) % 1.0

            if use_kdtree:
                k_cart_point = Brec_cart @ p.k_frac
                _, idx = tree.query(k_cart_point)
                p.ik = int(ik_list[idx])
                rep = ik_to_rep[p.ik]
            else:
                p.ik = nearest_ik_from_frac(p.k_frac, ik_list, k_frac_grid)
                rep  = ik_to_rep[p.ik]

            if not interpolate:
                p.v = events.v_group[rep].copy()
            else:
                k_cart_point = Brec_cart @ p.k_frac
                dists, idxs = tree.query(k_cart_point, k=6)
                if np.isscalar(idxs):
                    idxs = np.array([idxs])
                    dists = np.array([dists])
                else:
                    idxs = np.array(idxs)
                    dists = np.array(dists)
                neighbors_ik = ik_list[idxs]
                rep_indices = [ik_to_rep[int(k)] for k in neighbors_ik]
                v_neighbors = events.v_group[rep_indices]
                d2 = dists ** 2
                if len(d2) > 1 and np.any(d2[1:] > 0):
                    denom = np.mean(d2[1:]) + 1e-30
                else:
                    denom = d2[0] + 1e-30
                weights = np.exp(-d2 / denom)
                weights /= weights.sum()
                p.v = np.average(v_neighbors, axis=0, weights=weights)

            p.E_meV = float(events.enk_above_CBM_meV[rep])
            _accumulate_snapshot(p.t, dt, p.v[0])
            p.time_vx += p.v[0] * dt
            p.time_E  += p.E_meV * dt
            p.t       += dt
            if p.t >= T_sim:
                break

            evs, cum = state_cum[key]
            if len(evs) == 0:
                continue
            u = rng.random()
            idx_in = max(0, np.searchsorted(cum, u) - 1)
            ev = int(evs[idx_in])
            
            # NEW: classify scattering event
            # Proper emission/absorption classification
            if events.Gamma_em[ev] > 0:
                emission_count += 1
            elif events.Gamma_abs[ev] > 0:
                absorption_count += 1
            else:
                # Should never happen, but for safety
                emission_count += 1
            
            p.k_frac = events.kq_frac[ev].copy()
            if use_kdtree:
                k_cart_point = Brec_cart @ p.k_frac
                _, idx2 = tree.query(k_cart_point)
                p.ik = int(ik_list[idx2])
                rep2 = ik_to_rep[p.ik]
            else:
                p.ik = nearest_ik_from_frac(p.k_frac, ik_list, k_frac_grid)
                rep2 = ik_to_rep[p.ik]
            p.ib     = int(events.jb[ev])
            p.v      = events.v_scattered[ev].copy()
            p.E_meV  = float(events.enkq_above_CBM_meV[ev])
            p.scat_count += 1

    # local totals and summaries
    local_total_time = len(particles) * T_sim
    sum_time_vx = sum(pp.time_vx for pp in particles)
    sum_time_E  = sum(pp.time_E  for pp in particles)
    sum_scats   = sum(pp.scat_count for pp in particles)
    mean_s = np.mean([pp.scat_count for pp in particles]) if len(particles) > 0 else 0.0

    # local v_snapshot assembly (not normalized)
    v_snapshot_local = np.zeros_like(bin_time)
    for i in range(len(bin_time)):
        if bin_time[i] > 0:
            v_snapshot_local[i] = bin_vtime[i] / bin_time[i]
        else:
            v_snapshot_local[i] = 0.0
 
    # --- add: collect final states (ik, ib) for each local particle ---
    final_states = [(int(pp.ik), int(pp.ib)) for pp in particles]

    # Return local totals and arrays for reduction
    return {
        'sum_time_vx': sum_time_vx,
        'sum_time_E': sum_time_E,
        'sum_scats': sum_scats,
        'local_total_time': local_total_time,
        'mean_s': mean_s,
        'snapshot_times': snapshot_times,
        'bin_time': bin_time,
        'bin_vtime': bin_vtime,
        'pop_bins': pop_bins,
        'final_states': final_states,
         # NEW ENTRIES
        'emission_count': emission_count,
        'absorption_count': absorption_count
        
        
    }

# ---------- Main with MPI ----------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='MPI Full-band EMC runner')
    parser.add_argument('--h5', type=str, default='scattering_events.h5', help='HDF5 events file')
    parser.add_argument('--npar', type=int, default=2000, help='Number of particles (ensemble)')
    parser.add_argument('--tsim', type=float, default=5e-12, help='Simulation time per particle (s)')
    parser.add_argument('--Emin', type=float, default=1e3, help='Min field V/m')
    parser.add_argument('--Emax', type=float, default=1e8, help='Max field V/m')
    parser.add_argument('--nE', type=int, default=6, help='Number of field points')
    parser.add_argument('--interp', action='store_true', help='Enable v(k) interpolation (slow)')
    parser.add_argument('--maxsteps', type=int, default=1000, help='Max steps per particle (cap)')
    args = parser.parse_args()

    if rank == 0:
        t0 = time.time()
        print(f"[rank 0] Loading HDF5 ...")
        events = load_h5(args.h5)
        
        # ---------------------------------------------------
        # AUTO-DETECT USED BANDS (NEW CODE)
        # ---------------------------------------------------
        used_bands = np.unique(events.ib)
        band_map = {old: new for new, old in enumerate(used_bands)}

        events.ib = np.array([band_map[b] for b in events.ib])
        events.jb = np.array([band_map[b] for b in events.jb])

        events.used_bands = used_bands
        events.band_map = band_map

        print("Auto-band detection:")
        print("  Physical bands detected:", used_bands)
        print("  Remap old → new:", band_map)
        
        
        print(f"[rank 0] Building state index ...")
        state_to_evs, state_Gamma_tot, state_cum = build_state_index(events)
        print(f"[rank 0] Building k-grid and KD-tree ...")
        ik_list, ik_to_rep, k_frac_grid, k_cart_grid, tree, Brec_cart = build_kgrid_cart(events)
        print(f"[rank 0] Unique ik count = {len(ik_list)}")
    else:
        events = None
        state_to_evs = None
        state_Gamma_tot = None
        state_cum = None
        ik_list = None
        ik_to_rep = None
        k_frac_grid = None
        k_cart_grid = None
        tree = None
        Brec_cart = None

    # Broadcast events and indices to all ranks (EventsTable is picklable)
    events = comm.bcast(events, root=0)
    state_to_evs = comm.bcast(state_to_evs, root=0)
    state_Gamma_tot = comm.bcast(state_Gamma_tot, root=0)
    state_cum = comm.bcast(state_cum, root=0)
    ik_list = comm.bcast(ik_list, root=0)
    ik_to_rep = comm.bcast(ik_to_rep, root=0)
    k_frac_grid = comm.bcast(k_frac_grid, root=0)
    k_cart_grid = comm.bcast(k_cart_grid, root=0)
    tree = comm.bcast(tree, root=0)
    Brec_cart = comm.bcast(Brec_cart, root=0)

    # Initialize particles (each rank will create the same total list then slice)
    if rank == 0:
        print(f"[rank 0] Initializing total particles (npar={args.npar}) ...")
    particles_all = init_particles(args.npar, events, ik_list, ik_to_rep, T=300.0)
    # slice across ranks
    particles_local = particles_all[rank::size]
    n_local = len(particles_local)
    if rank == 0:
        print(f"[rank 0] Particles per rank (approx): { [len(particles_all[i::size]) for i in range(size)] }")

    # Prepare fields
    E_fields = np.logspace(math.log10(args.Emin), math.log10(args.Emax), args.nE)

    # per-field storage on rank 0
    drifts = []
    meanEs = []
    meanScats = []
    all_pops = []
    snapshot_times = None
    all_final_pops = []   # <--- ADD THIS LINE
    em_totals = []
    abs_totals = []

    # Loop over fields (each rank computes locally and reduces)
    for E in E_fields:
        if rank == 0:
            print(f"[rank 0] Starting field E = {E:.3e} V/m")
            field_folder = os.path.join("Drift_out", f"E_{E:.3e}")
            os.makedirs(field_folder, exist_ok=True)
        # Barrier to synchronise output directories (optional)
        comm.Barrier()

        # Each rank runs on its local subset
        local_res = run_field(
            events, state_Gamma_tot, state_cum,
            ik_list, ik_to_rep,
            k_frac_grid, k_cart_grid,
            tree, Brec_cart,
            particles_local, E, args.tsim,
            interpolate=args.interp, MAX_STEPS=args.maxsteps
        )

        # Reduce scalar totals
        sum_time_vx_total = comm.reduce(local_res['sum_time_vx'], op=MPI.SUM, root=0)
        sum_time_E_total  = comm.reduce(local_res['sum_time_E'], op=MPI.SUM, root=0)
        sum_scats_total   = comm.reduce(local_res['sum_scats'], op=MPI.SUM, root=0)
        local_total_time  = comm.reduce(local_res['local_total_time'], op=MPI.SUM, root=0)
        mean_s_global     = comm.reduce(local_res['mean_s'], op=MPI.SUM, root=0)  # will divide by size below to get average of means (not needed)
        em_total          = comm.reduce(local_res['emission_count'],  op=MPI.SUM, root=0)
        abs_total         = comm.reduce(local_res['absorption_count'], op=MPI.SUM, root=0)
        
        if rank == 0:
            em_totals.append(em_total)
            abs_totals.append(abs_total)

        # Reduce bin_time and bin_vtime arrays (element-wise sum)
        bin_time_total = None
        bin_vtime_total = None
        pop_bins_total = None

        # create buffers on root for Reduce ops for numpy arrays
        if rank == 0:
            bin_time_total = np.zeros_like(local_res['bin_time'])
            bin_vtime_total = np.zeros_like(local_res['bin_vtime'])
            pop_bins_total = np.zeros_like(local_res['pop_bins'])
        # Use MPI.Reduce with numpy arrays (works with mpi4py)
        comm.Reduce([local_res['bin_time'], MPI.DOUBLE], [bin_time_total, MPI.DOUBLE], op=MPI.SUM, root=0)
        comm.Reduce([local_res['bin_vtime'], MPI.DOUBLE], [bin_vtime_total, MPI.DOUBLE], op=MPI.SUM, root=0)
        comm.Reduce([local_res['pop_bins'], MPI.LONG],    [pop_bins_total,    MPI.LONG],    op=MPI.SUM, root=0)

        # snapshot_times identical across ranks: pick from local_res
        if snapshot_times is None:
            snapshot_times = local_res['snapshot_times']


        # --- gather final states from all ranks ---
        all_final_states = comm.gather(local_res['final_states'], root=0)

        # --- root: build histogram ---
        if rank == 0:
            # flatten list
            flat = [item for sub in all_final_states for item in sub]

            # ik_list is already broadcasted, so create index map
            ik_pos_map = {int(k): idx for idx, k in enumerate(ik_list)}

            N_ik = len(ik_list)
            N_band = int(np.max(events.ib) + 1)
            final_pop = np.zeros((N_ik, N_band), dtype=np.int64)

            for ik_val, ib_val in flat:
                if ik_val in ik_pos_map:  # safety check
                    ikpos = ik_pos_map[ik_val]
                    final_pop[ikpos, ib_val] += 1

            # store per-field final population
            all_final_pops.append(final_pop)


        # Rank 0 finalizes observables and writes outputs
        if rank == 0:
            total_time = local_total_time
            drift = sum_time_vx_total / total_time
            meanE = sum_time_E_total / total_time
            mean_s = sum_scats_total / args.npar  # average scat per particle across all
            drifts.append(drift)
            meanEs.append(meanE)
            meanScats.append(mean_s)
            all_pops.append(pop_bins_total)

            # finalize v_snapshot from totals
            v_snapshot_tot = np.zeros_like(bin_time_total)
            for i in range(len(bin_time_total)):
                if bin_time_total[i] > 0:
                    v_snapshot_tot[i] = bin_vtime_total[i] / bin_time_total[i]
                else:
                    v_snapshot_tot[i] = 0.0

            # write vx(t)
            fname = os.path.join(field_folder, "vx_time.txt")
            try:
                with open(fname, "w") as f:
                    f.write("# time(s)    vx(m/s)\n")
                    for tval, vxval in zip(snapshot_times, v_snapshot_tot):
                        f.write(f"{tval: .6e}   {vxval: .6e}\n")
                print(f"[rank 0] Saved vx(t) -> {fname}")
            except Exception as exc:
                print("[rank 0] Failed to write vx(t) file:", exc)

            # steady state detection and writing
            t_ss = detect_steady_state(snapshot_times, v_snapshot_tot)
            ssfile = os.path.join(field_folder, "steady_state.txt")
            try:
                with open(ssfile, "w") as fss:
                    fss.write(f"# Steady state info for E={E:.3e} V/m\n")
                    if t_ss is not None:
                        vx_ss = v_snapshot_tot[np.argmin(abs(snapshot_times - t_ss))]
                        fss.write(f"steady state time in s    : {t_ss:.6e}\n")
                        fss.write(f"steady state vx in ms^-1  : {vx_ss:.6e}\n")
                        fss.write("steady state status:       : reached\n")
                        print(f"[rank 0]    Steady-state reached at t = {t_ss:.3e} s")
                    else:
                        fss.write("steady state_time (s)      :  not_reached\n")
                        fss.write("steady state vx ms^-1      :  not_reached\n")
                        fss.write("steady state status        :  not_reached\n")
                        print("[rank 0]    WARNING: steady state not detected.")
                    fss.write(f"drift velocity in ms^-1   :    {drift:.6e}\n")
                    fss.write(f"mean energy in meV        :    {meanE:.6f}\n")
                    fss.write(f"avgerage scatter count    :    {mean_s:.2f}\n")
                print(f"[rank 0]    Saved steady-state info -> {ssfile}")
            except Exception as exc:
                print("[rank 0] Failed to write steady-state file:", exc)

        # sync ranks before next field
        comm.Barrier()

    # end field loop

    # Final save (rank 0)
    if rank == 0:
        np.savez(
            'fb_emc_results_fast_mpi.npz',
            E_fields=E_fields,
            drifts=np.array(drifts),
            meanEs=np.array(meanEs),
            meanScats=np.array(meanScats),
            populations=np.array(all_pops, dtype=object),
            snapshot_times=snapshot_times,
            final_populations=np.array(all_final_pops, dtype=object),   # <-- added
            ik_list=np.array(ik_list, dtype=int),
            used_bands=events.used_bands,
            band_map=events.band_map,
            emission_counts = np.array(em_totals),
            absorption_counts = np.array(abs_totals)
        )

        # simple plots
        plt.figure(figsize=(10,4))
        plt.subplot(1,2,1)
        plt.semilogx(E_fields, drifts, marker='o')
        plt.xlabel('Electric field |E| (V/m)')
        plt.ylabel('Drift velocity <v_x> (m/s)')
        plt.grid(True)

        plt.subplot(1,2,2)
        plt.semilogx(E_fields, meanEs, marker='o')
        plt.xlabel('Electric field |E| (V/m)')
        plt.ylabel('Mean energy above CBM (meV)')
        plt.grid(True)

        plt.tight_layout()
        plt.savefig('fb_emc_results_fast_mpi.png')
        print("[rank 0] Saved fb_emc_results_fast_mpi.npz and fb_emc_results_fast_mpi.png")
        print("[rank 0] Total runtime (rank 0): %.2f s" % (time.time() - t0))

    # ensure all ranks finish cleanly
    comm.Barrier()
    if rank == 0:
        print("[rank 0] MPI run complete.")

