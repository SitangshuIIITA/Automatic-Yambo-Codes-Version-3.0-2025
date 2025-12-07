#!/usr/bin/env python3
"""
test_clean.py — SINGLE-FILE MONTE-CARLO ENGINE (Complete)

This is a complete, self-contained version that uses ONLY scattering_events.h5.
It reconstructs the unique k-grid, builds channel lists from event rows,
uses event-level post-scatter velocities, and runs an ensemble Monte Carlo
with an adaptive drift-velocity stability check and optional diagnostics.

Usage examples:
  python3 test_clean.py --events scattering_events.h5 --field 1e5 --npart 2000 --steps 20000
  python3 test_clean.py --events scattering_events.h5 --fields 1e5 1e6 1e7 --diagnostics

Outputs:
  - Console logs with stability diagnostics
  - scatter_events.txt (logged scatter events up to --max-log-events)
  - (optionally) vd_results.txt when --vd-output is set
"""

import argparse
import os
import sys
import math
import numpy as np
import h5py

# optional accelerators
try:
    from numba import njit
    have_numba = True
except Exception:
    have_numba = False

# Physical constants
hbar = 1.054571817e-34
kB_eV_per_K = 8.617333262145e-5
e_charge = 1.602176634e-19

# ---------------- utilities ----------------

def reconstruct_kgrid_from_events(g):
    """Return Nk, Nq, idx_k, idx_q, ikq_map, k_cart
    """
    k_frac = np.asarray(g['k_frac'][...])
    kq_frac = np.asarray(g['kq_frac'][...])
    q_frac  = np.asarray(g['q_frac'][...])

    # Load B_reciprocal (the TRUE reciprocal matrix)
    B = np.asarray(g['B_reciprocal'][...])  # <-- correct

    # round & unique K
    k_round = np.round(k_frac, 9)
    uniq_k, idx_k = np.unique(k_round, axis=0, return_inverse=True)
    Nk = int(uniq_k.shape[0])

    # Compute Cartesian k correctly
    k_cart = uniq_k.dot(B.T)               # <-- correct formula

    # round & unique Q
    q_round = np.round(q_frac, 9)
    uniq_q, idx_q = np.unique(q_round, axis=0, return_inverse=True)
    Nq = int(uniq_q.shape[0])

    # Build lookup for k+q → unique k index
    kq_round = np.round(kq_frac, 9)
    lookup = {tuple(row): i for i, row in enumerate(uniq_k)}
    ikq = np.empty(kq_round.shape[0], dtype=np.int32)
    for i, r in enumerate(kq_round):
        ikq[i] = lookup.get(tuple(r), -1)
        if ikq[i] < 0:
            d = np.sum((uniq_k - r)**2, axis=1)
            ikq[i] = int(np.argmin(d))

    return Nk, Nq, idx_k.astype(np.int32), idx_q.astype(np.int32), ikq.astype(np.int32), k_cart


# ---------------- channels ----------------

def build_channels_from_events(g, idx_k, idx_q, ikq, lcb, mode_flag='both', verbose=False):
    # read arrays (expect shapes matching number of events)
    ib = np.asarray(g['ib'][...]).astype(np.int32)
    jb = np.asarray(g['jb'][...]).astype(np.int32)
    imode = np.asarray(g['imode'][...]).astype(np.int32)
    Gamma = np.asarray(g['Gamma_w'][...]).astype(np.float64)
    Gabs = np.asarray(g['Gamma_abs_w'][...]).astype(np.float64) if 'Gamma_abs_w' in g else np.zeros_like(Gamma)
    Gem = np.asarray(g['Gamma_em_w'][...]).astype(np.float64) if 'Gamma_em_w' in g else np.zeros_like(Gamma)
    vsc = np.asarray(g['v_scattered_cart'][...]).astype(np.float64)
    enkq = np.asarray(g['enkq_meV'][...]).astype(np.float64)

    # event->unique mappings from reconstruct
    ik_map = idx_k
    iq_map = idx_q

    # select conduction-related events
    mask = (ib >= lcb) & (jb >= lcb)
    if mode_flag == 'intra':
        mask &= (ib == jb)
    elif mode_flag == 'inter_up':
        mask &= (jb > ib)
    elif mode_flag == 'inter_down':
        mask &= (jb < ib)

    ik_sel = ik_map[mask]
    iq_sel = iq_map[mask]
    imode_sel = imode[mask]
    ib_sel = ib[mask]; jb_sel = jb[mask]
    Gamma_sel = Gamma[mask]; Gabs_sel = Gabs[mask]; Gem_sel = Gem[mask]
    vsc_sel = vsc[mask]; enkq_sel = enkq[mask]
    ikq_sel = ikq[mask]

    # sort by initial ik to build contiguous blocks
    order = np.argsort(ik_sel)
    ik_sel = ik_sel[order]; iq_sel = iq_sel[order]; imode_sel = imode_sel[order]
    ib_sel = ib_sel[order]; jb_sel = jb_sel[order]
    Gamma_sel = Gamma_sel[order]; Gabs_sel = Gabs_sel[order]; Gem_sel = Gem_sel[order]
    vsc_sel = vsc_sel[order]; enkq_sel = enkq_sel[order]; ikq_sel = ikq_sel[order]

    if ik_sel.size == 0:
        raise RuntimeError('No scattering channels after filtering.')

    Nk_used = int(np.max(ik_sel)) + 1
    ch_start = np.full(Nk_used, -1, dtype=np.int32)
    ch_len = np.zeros(Nk_used, dtype=np.int32)

    uniq, counts = np.unique(ik_sel, return_counts=True)
    pos = 0
    for u, c in zip(uniq, counts):
        ch_start[int(u)] = pos
        ch_len[int(u)] = int(c)
        pos += int(c)

    ch_cumsum = np.zeros_like(Gamma_sel)
    for u in uniq:
        s = ch_start[int(u)]; L = ch_len[int(u)]
        ch_cumsum[s:s+L] = np.cumsum(Gamma_sel[s:s+L])

    channels = {
        'ch_ik': ik_sel,
        'ch_iq': iq_sel,
        'ch_imode': imode_sel,
        'ch_ib': ib_sel,
        'ch_jb': jb_sel,
        'ch_Gamma_w': Gamma_sel,
        'ch_Gamma_abs_w': Gabs_sel,
        'ch_Gamma_em_w': Gem_sel,
        'ch_vscat': vsc_sel,
        'ch_enkq_meV': enkq_sel,
        'ch_ikq': ikq_sel,
        'ch_start': ch_start,
        'ch_len': ch_len,
        'ch_cumsum': ch_cumsum
    }
    if verbose:
        print(f"[ok] Built {ik_sel.size} channels over {Nk_used} k-points (mode={mode_flag}).")
    return channels

# ---------------- per-kpoint properties ----------------

def build_kpoint_properties(g, idx_k, Nk):
    v_group = np.asarray(g['v_group_cart'][...]).astype(np.float64)
    enk = np.asarray(g['enk_meV'][...]).astype(np.float64)

    vel = np.zeros((Nk, 3), dtype=np.float64)
    energy_meV = np.zeros((Nk,), dtype=np.float64)
    count = np.zeros((Nk,), dtype=np.int32)

    for ev_idx, kidx in enumerate(idx_k):
        kid = int(kidx)
        if kid < 0 or kid >= Nk:
            continue
        vel[kid] += v_group[ev_idx]
        energy_meV[kid] += enk[ev_idx]
        count[kid] += 1

    mask = count > 0
    vel[mask] /= count[mask, None]
    energy_meV[mask] /= count[mask]
    energy_eV = energy_meV / 1000.0
    return vel, energy_eV

# ---------------- selection kernels ----------------

def select_channels_python(particles_ik, to_scatter_idx, rand_r, ch_start, ch_len, ch_cumsum):
    m = to_scatter_idx.shape[0]
    chosen_idx = np.empty((m,), dtype=np.int32)
    for ii in range(m):
        p = int(to_scatter_idx[ii])
        ik = int(particles_ik[p])
        s = int(ch_start[ik]); L = int(ch_len[ik])
        if L <= 0:
            chosen_idx[ii] = -1
            continue
        target = rand_r[ii]
        chosen = -1
        for j in range(L):
            idx = s + j
            if target <= ch_cumsum[idx]:
                chosen = idx
                break
        if chosen == -1:
            chosen = s + L - 1
        chosen_idx[ii] = int(chosen)
    return chosen_idx

if have_numba:
    @njit
    def select_channels_numba(particles_ik, to_scatter_idx, rand_r, ch_start, ch_len, ch_cumsum):
        m = to_scatter_idx.shape[0]
        chosen_idx = np.empty((m,), dtype=np.int32)
        for ii in range(m):
            p = int(to_scatter_idx[ii])
            ik = int(particles_ik[p])
            s = int(ch_start[ik]); L = int(ch_len[ik])
            if L <= 0:
                chosen_idx[ii] = -1
                continue
            target = rand_r[ii]
            chosen = -1
            for j in range(L):
                idx = s + j
                if target <= ch_cumsum[idx]:
                    chosen = idx
                    break
            if chosen == -1:
                chosen = s + L - 1
            chosen_idx[ii] = int(chosen)
        return chosen_idx

# ---------------- Monte Carlo (adaptive stability) ----------------

def run_mc(ch, vel, energy, k_cart, npart=2000, steps=20000, dt=1e-15,
           E_field=1e5, T=300, seed=None, scatter_log_fh=None, max_log_events=5000,
           verbose=False, steps_min=5000, steps_extend=None, steps_max=100000,
           stability_threshold=0.10):
    """Run ensemble MC. Returns (mean_vd, logged, energy_out, stability_report).
    """
    if steps_extend is None:
        steps_extend = steps_min

    if seed is not None:
        np.random.seed(seed)

    Nk = int(k_cart.shape[0])
    # thermal init probabilities
    E_CBM = float(np.nanmin(energy))
    kBT = kB_eV_per_K * float(T)
    probs = np.exp(-(energy - E_CBM) / kBT)
    probs /= np.sum(probs)

    # initialize particles
    part_iks = np.random.choice(Nk, size=npart, p=probs).astype(np.int32)
    part_kcart = k_cart[part_iks].astype(np.float64)
    part_vel = vel[part_iks].astype(np.float64)
    part_bands = np.zeros((npart,), dtype=np.int32)
    part_last_event = np.zeros((npart,), dtype=np.int8)

    # channel arrays
    ch_start = ch['ch_start'].astype(np.int32)
    ch_len = ch['ch_len'].astype(np.int32)
    ch_cumsum = ch['ch_cumsum'].astype(np.float64)
    ch_ikq = ch['ch_ikq'].astype(np.int32)
    ch_vscat = ch['ch_vscat'].astype(np.float64)
    ch_enkq_meV = ch['ch_enkq_meV'].astype(np.float64)
    ch_jb = ch['ch_jb'].astype(np.int32) if 'ch_jb' in ch else np.zeros((ch['ch_ik'].size,), dtype=np.int32)
    ch_Gamma_w = ch['ch_Gamma_w'].astype(np.float64)
    ch_Gamma_abs_w = ch.get('ch_Gamma_abs_w', None)
    ch_Gamma_em_w = ch.get('ch_Gamma_em_w', None)

    # total gamma per k
    total_Gamma = np.zeros((Nk,), dtype=np.float64)
    for ik_i in range(Nk):
        s = ch_start[ik_i]; L = ch_len[ik_i]
        if s >= 0 and L > 0:
            total_Gamma[ik_i] = float(np.sum(ch_Gamma_w[s:s+L]))

    # histograms
    max_band_idx = int(np.max(ch_jb)) if ch_jb.size > 0 else 0
    nbands_total = max_band_idx + 1
    abs_hist = np.zeros((nbands_total, Nk), dtype=np.int64)
    em_hist = np.zeros((nbands_total, Nk), dtype=np.int64)

    # k-drift per dt
    delta_k_per_dt = -e_charge * np.array([E_field, 0.0, 0.0]) / hbar
    delta_k_step = delta_k_per_dt * dt

    v_drift_time = np.zeros((0,), dtype=np.float64)
    logged = 0
    total_steps_done = 0

    def run_block(n_steps):
        nonlocal part_kcart, part_iks, part_vel, part_bands, part_last_event, logged
        vbuf = np.zeros((n_steps,), dtype=np.float64)
        for t_local in range(n_steps):
            part_kcart += delta_k_step
            # map to nearest k (brute-force; small Nk expected)
            d = ((part_kcart[:, None, :] - k_cart[None, :, :])**2).sum(axis=2)
            part_iks = np.argmin(d, axis=1).astype(np.int32)
            part_vel = vel[part_iks]

            total_G = total_Gamma[part_iks]
            p_scatter = 1.0 - np.exp(-total_G * dt)

            rand_unif = np.random.random(npart)
            to_scatter = rand_unif < p_scatter

            if np.any(to_scatter):
                scatter_idx = np.nonzero(to_scatter)[0].astype(np.int32)
                rvals = np.random.random(scatter_idx.size) * total_G[scatter_idx]
                if have_numba:
                    chosen = select_channels_numba(part_iks, scatter_idx, rvals, ch_start, ch_len, ch_cumsum)
                else:
                    chosen = select_channels_python(part_iks, scatter_idx, rvals, ch_start, ch_len, ch_cumsum)

                for ii, p in enumerate(scatter_idx):
                    ch_idx = int(chosen[ii])
                    if ch_idx < 0:
                        continue
                    old_ik = int(part_iks[p])
                    ik_prime = int(ch_ikq[ch_idx])

                    Ga = float(ch_Gamma_abs_w[ch_idx]) if ch_Gamma_abs_w is not None else 0.0
                    Ge = float(ch_Gamma_em_w[ch_idx]) if ch_Gamma_em_w is not None else 0.0
                    Gsum = Ga + Ge
                    if Gsum <= 0:
                        ev_type = 0
                    else:
                        if Ga == 0.0:
                            ev_type = 2
                        elif Ge == 0.0:
                            ev_type = 1
                        else:
                            ev_type = 1 if (np.random.random() < (Ga / Gsum)) else 2

                    new_band = int(ch_jb[ch_idx]) - 1 if ch_jb.size>0 else 0
                    if 0 <= new_band < abs_hist.shape[0] and 0 <= ik_prime < Nk:
                        if ev_type == 1:
                            abs_hist[new_band, ik_prime] += 1
                        elif ev_type == 2:
                            em_hist[new_band, ik_prime] += 1

                    # update particle state
                    new_vel = ch_vscat[ch_idx]
                    part_iks[p] = ik_prime
                    part_kcart[p] = k_cart[ik_prime].copy()
                    part_vel[p] = new_vel.copy()
                    part_bands[p] = new_band
                    part_last_event[p] = ev_type

                    if (scatter_log_fh is not None) and ((max_log_events is None) or (logged < max_log_events)):
                        k_old = k_cart[old_ik]; k_new = k_cart[ik_prime]
                        scatter_log_fh.write(
                            f"{E_field: .6e} {total_steps_done + t_local:6d} {int(p):6d} {old_ik:6d} {int(ch['ch_iq'][ch_idx]):6d} {ik_prime:6d} "
                            f"{k_old[0]: .6e} {k_old[1]: .6e} {k_old[2]: .6e} "
                            f"{k_new[0]: .6e} {k_new[1]: .6e} {k_new[2]: .6e}\n"
                        )
                        logged += 1

            mean_v = np.mean(part_vel, axis=0)
            vbuf[t_local] = float(mean_v[0])
        return vbuf

    # adaptive loop
    while total_steps_done < steps:
        if total_steps_done == 0:
            run_now = min(steps_min, steps - total_steps_done)
        else:
            run_now = min(steps_extend, steps - total_steps_done)
        vblock = run_block(run_now)
        v_drift_time = np.concatenate((v_drift_time, vblock))
        total_steps_done += run_now

        n_check = max(1, int(0.30 * v_drift_time.size))
        last30 = v_drift_time[-n_check:]
        mean_last30 = float(np.mean(last30))
        std_last30 = float(np.std(last30))
        rel_std = std_last30 / (abs(mean_last30) + 1e-30)

        n20 = max(1, int(0.20 * v_drift_time.size))
        last20 = v_drift_time[-n20:]
        mean_last20 = float(np.mean(last20))
        rel_diff = abs(mean_last20 - mean_last30) / (abs(mean_last30) + 1e-30)

        if verbose:
            print(f"[stability] steps={total_steps_done} mean_last30={mean_last30:.6e} rel_std={rel_std:.3e} rel_diff20_30={rel_diff:.3e}")

        if (rel_std < stability_threshold) and (rel_diff < stability_threshold):
            if verbose:
                print(f"[stability] Converged after {total_steps_done} steps (threshold={stability_threshold:.3f})")
            break

        if (total_steps_done >= steps_max) or (total_steps_done >= steps):
            if verbose:
                print(f"[stability] Failed to converge within max steps={steps_max} or requested steps={steps}.")
            break

        if verbose:
            print(f"[stability] Not converged, extending by {steps_extend} steps...")

    mean_vd = float(np.mean(v_drift_time[int(0.7 * v_drift_time.size):])) if v_drift_time.size>0 else 0.0

    energy_out = {
        'final_part_iks': part_iks.copy(),
        'final_part_bands': part_bands.copy(),
        'final_part_last_event': part_last_event.copy(),
        'abs_hist': abs_hist,
        'em_hist': em_hist,
        'v_drift_time': v_drift_time
    }

    stability_report = {
        'steps_done': int(total_steps_done),
        'mean_last30': mean_last30 if 'mean_last30' in locals() else float('nan'),
        'std_last30': std_last30 if 'std_last30' in locals() else float('nan'),
        'rel_std': rel_std if 'rel_std' in locals() else float('nan'),
        'rel_diff20_30': rel_diff if 'rel_diff' in locals() else float('nan')
    }

    return mean_vd, logged, energy_out, stability_report

# ---------------- main ----------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--events', required=True)
    p.add_argument('--steps', type=int, default=10000)
    p.add_argument('--npart', type=int, default=20000)
    p.add_argument('--T', type=float, default=100.0)
    p.add_argument('--field', type=float, default=None)
    p.add_argument('--fields', nargs='*', type=float, default=None)
    p.add_argument('--mode', default='both')
    p.add_argument('--seed', type=int, default=None)
    p.add_argument('--max-log-events', type=int, default=5000)
    p.add_argument('--steps-min', type=int, default=5000)
    p.add_argument('--steps-extend', type=int, default=20000)
    p.add_argument('--steps-max', type=int, default=100000)
    p.add_argument('--stability-threshold', type=float, default=0.1)
    p.add_argument('--diagnostics', action='store_true')
    p.add_argument('--vd-output', type=str, default=None, help='Write vd(E) results to this file')
    args = p.parse_args()

    events_path = args.events
    if not os.path.exists(events_path):
        print('[error] events HDF5 not found:', events_path)
        sys.exit(1)

    field_list = None
    if args.fields is not None and len(args.fields) > 0:
        field_list = args.fields
    elif args.field is not None:
        field_list = [args.field]
    else:
        # sensible default if no field provided
        field_list = [1e5]

    results = []

    with h5py.File(events_path, 'r') as f:
        if '/scattering_events' not in f:
            print('[error] /scattering_events not found in HDF5')
            sys.exit(1)
        g = f['/scattering_events']

        Nk, Nq, idx_k, idx_q, ikq, k_cart = reconstruct_kgrid_from_events(g)
        ib = np.asarray(g['ib'][...])
        lcb = int(np.min(ib))

        ch = build_channels_from_events(g, idx_k, idx_q, ikq, lcb, mode_flag=args.mode, verbose=True)
        vel, energy = build_kpoint_properties(g, idx_k, Nk)

        # diagnostics
        if args.diagnostics:
            print('================ Diagnostics ================')
            print('Nk (unique k):', Nk)
            print('Number of channels:', ch['ch_ik'].size)
            print('-- v_group x statistics --')
            vx = vel[:,0]
            print('mean vx =', np.mean(vx), ' min =', np.min(vx), ' max =', np.max(vx))
            print('-- k_cart x-range --')
            try:
                print('min kx =', np.min(k_cart[:,0]), ' max kx =', np.max(k_cart[:,0]))
            except Exception:
                pass
            # total Gamma per k
            total_Gamma_diag = np.zeros((Nk,))
            for ik_i in range(Nk):
                s = ch['ch_start'][ik_i]; L = ch['ch_len'][ik_i]
                if s>=0 and L>0:
                    total_Gamma_diag[ik_i] = float(np.sum(ch['ch_Gamma_w'][s:s+L]))
            print('-- total Gamma per k --')
            print('mean =', np.mean(total_Gamma_diag), ' min =', np.min(total_Gamma_diag), ' max =', np.max(total_Gamma_diag))
            print('===========================================\n')

        scatter_log = 'scatter_events.txt'
        scatter_fh = open(scatter_log, 'w')
        scatter_fh.write('# field[V/m] timestep particle old_ik iq ik_prime k_old_x k_old_y k_old_z k_new_x k_new_y k_new_z\n')

        print('Running fields:', field_list)
        for Ef in field_list:
            print(f'[run] Field = {Ef} V/m')
            vd, logged, energy_out, stability_report = run_mc(
                ch, vel, energy, k_cart,
                npart=args.npart, steps=args.steps, dt=1e-16,
                E_field=Ef, T=args.T, seed=args.seed,
                scatter_log_fh=scatter_fh, max_log_events=args.max_log_events,
                verbose=True, steps_min=args.steps_min, steps_extend=args.steps_extend,
                steps_max=args.steps_max, stability_threshold=args.stability_threshold
            )
            print(f"[ok] Ef={Ef:.3e} -> vd={vd:.6e} m/s logged={logged}")
            print(f"[stability-final] steps={stability_report['steps_done']} rel_std={stability_report['rel_std']:.3e} rel_diff={stability_report['rel_diff20_30']:.3e}")
            # append energy_out so diagnostics can use final particle states
            results.append((Ef, vd, logged, stability_report, energy_out))

        # close scatter log
        scatter_fh.close()

        # ---------------- Diagnostics per-field (initial vs final) ----------------
        print('\n===== POST-RUN DIAGNOSTICS =====')
        # prepare Boltzmann sampling for initial distribution
        E_CBM = float(np.nanmin(energy))
        probs_init = np.exp(-(energy - E_CBM) / (kB_eV_per_K * float(args.T)))
        probs_init /= probs_init.sum()
        Nsample = int(min(10000, max(1000, args.npart)))

        for (Ef, vd, logged, st, energy_out) in results:
            print(f"\n--- Field {Ef:.3e} diagnostics ---")
            # initial sample (Boltzmann)
            init_k = np.random.choice(int(Nk), size=Nsample, p=probs_init)
            init_vx = vel[init_k, 0]
            init_Ee = energy[init_k]  # eV
            mean_init_v = float(np.mean(init_vx))
            mean_init_E = float(np.mean(init_Ee))

            # final MC sample
            final_iks = energy_out.get('final_part_iks')
            if final_iks is None or final_iks.size == 0:
                print('[diagnostics] No final particle info available.')
                continue

            final_vx = vel[final_iks, 0]
            final_Ee = energy[final_iks]      # energy in eV

            # reference energy: conduction-band minimum
            E_CBM = np.min(energy)

            # shift energies to CBM reference
            mean_init_E_rel  = float(np.mean(mean_init_E - E_CBM))   # fix init
            mean_final_E_rel = float(np.mean(final_Ee - E_CBM))      # fix final

            mean_final_v = float(np.mean(final_vx))

            # corrected effective temperature
            T_eff = mean_final_E_rel / (1.5 * kB_eV_per_K)

            print(f"mean_init_vx={mean_init_v:.6e} m/s mean_final_vx={mean_final_v:.6e} m/s")
            print(f"mean_init_E_rel={mean_init_E_rel:.6e} eV mean_final_E_rel={mean_final_E_rel:.6e} eV")
            print(f"T_lattice={args.T:.2f} K T_eff={T_eff:.2f} K (rel diff={(T_eff-args.T)/args.T:.3e})")


            # energy distributions (initial vs final)
            try:
                plt.figure(figsize=(5,4))
                plt.hist(init_Ee, bins=60, density=True, alpha=0.6, label='initial')
                plt.hist(final_Ee, bins=60, density=True, alpha=0.6, label='final')
                plt.xlabel('Energy (eV)')
                plt.ylabel('Probability density')
                plt.title(f'Energy distribution E={Ef:.3e} V/m')
                plt.legend()
                fname = f'energy_dist_E_{int(Ef):d}.png'
                plt.tight_layout()
                plt.savefig(fname, dpi=150)
                plt.close()
                print('[diagnostics] Saved', fname)
            except Exception as e:
                print('[diagnostics] Failed to plot energy distribution:', e)

            # valley occupancy initial vs final (counts per k-index)
            try:
                init_counts = np.bincount(init_k, minlength=Nk)
                final_counts = np.bincount(final_iks, minlength=Nk)
                plt.figure(figsize=(6,3))
                idx = np.arange(Nk)
                width = 0.4
                plt.bar(idx - width/2, init_counts, width=width, label='initial')
                plt.bar(idx + width/2, final_counts, width=width, label='final')
                plt.xlabel('k-index (valley)')
                plt.ylabel('counts')
                plt.title(f'Valley occupancy E={Ef:.3e} V/m')
                plt.legend()
                fname = f'valley_occ_E_{int(Ef):d}.png'
                plt.tight_layout()
                plt.savefig(fname, dpi=150)
                plt.close()
                print('[diagnostics] Saved', fname)
            except Exception as e:
                print('[diagnostics] Failed valley occupancy plot:', e)

            # thermalization check
            rel_T = abs(T_eff - args.T) / args.T
            if rel_T < 0.10:
                print('[thermal-check] PASSED: T_eff within 10% of lattice T')
            else:
                print('[thermal-check] WARNING: T_eff differs from lattice T by >10%')

        print('===== END DIAGNOSTICS =====\\n')

    if args.vd_output is not None:
        with open(args.vd_output, 'w') as fh:
            fh.write('# Field[V/m]  vd[m/s]  logged  steps_done  rel_std  rel_diff20_30\n')
            for Ef, vd, logged, st, energy_out in results:
                fh.write(f"{Ef:.6e} {vd:.6e} {logged:d} {st['steps_done']:d} {st['rel_std']:.6e} {st['rel_diff20_30']:.6e}\n")

if __name__ == '__main__':
    main()

