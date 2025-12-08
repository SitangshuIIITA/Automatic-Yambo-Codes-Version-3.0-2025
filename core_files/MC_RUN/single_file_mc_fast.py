#!/usr/bin/env python3
"""
single_file_mc_fast_fixed_with_checks.py

Unified script (fixed) + stability checks (both global and per-field)
Inserted optional runtime checks guarded by --run-checks.
"""
import argparse, os, sys, time, math
import h5py, numpy as np

# optional accel
try:
    from scipy.spatial import cKDTree as KDTree
    have_kdtree = True
except Exception:
    KDTree = None; have_kdtree = False
try:
    from numba import njit
    have_numba = True
except Exception:
    have_numba = False

# ---- Physical constants (same as your reference) ----
hbar = 1.054571817e-34  # J s
eV_to_J = 1.602176634e-19
meV_to_J = eV_to_J * 1e-3
kB_eV_per_K = 8.617333262145e-5
kB_meV_per_K = kB_eV_per_K * 1000.0
pi = math.pi
# electron charge magnitude
e = 1.602176634e-19

def lorentzian_meV(delta_meV, gamma_meV):
    return (gamma_meV / math.pi) / (delta_meV * delta_meV + gamma_meV * gamma_meV)

# (compute_rates, build_channels_from_rows, select_channels_numba/select_channels_python,
#  run_mc) remain unchanged in core logic; see original script for details.
# For brevity, we include the full implementations (identical to your previous file) and
# then add the checks as requested.

# ---------- compute_rates: copied (structure) from mc_rate_intraband.py ----------
def detect_lcb_from_reference(energies, cbm_ref_eV):
    mins = np.nanmin(energies, axis=1)
    idx = int(np.argmin(np.abs(mins - cbm_ref_eV)))
    return idx + 1, mins

def compute_rates(h5path, gamma_meV=5.0, T=300.0, kbt_multiplier=3.0,
                  cbm_ref_eV=None, nq=101, mode_flag='both',
                  energy_tol_meV=None, verbose=True):
    if energy_tol_meV is None:
        energy_tol_meV = float(gamma_meV)

    with h5py.File(h5path, 'r') as fh:
        energies = fh['/bands/energies_eV'][...]
        nb, Nk = energies.shape

        # determine LCB
        if cbm_ref_eV is not None:
            lcb, mins = detect_lcb_from_reference(energies, cbm_ref_eV)
            if verbose:
                print(f"[info] CBM reference provided: {cbm_ref_eV:.6f} eV")
                for b_idx, m in enumerate(mins, start=1):
                    print(f"  band {b_idx:2d}: min = {m: .6f}")
                print(f"[info] Selected LCB = band {lcb}")
        else:
            try:
                conduction = fh['/meta/classification/conduction_bands'][...]
                if conduction.size == 0:
                    raise RuntimeError("No conduction bands in metadata")
                lcb = int(conduction[0])
                if verbose:
                    print(f"[info] Using EPW's CBM = band {lcb}")
            except Exception:
                raise RuntimeError("No conduction_bands metadata and no --cbm-ref provided. Please supply --cbm-ref")

        energy_idx = lcb - 1
        Eband = energies[energy_idx, :]
        E_CBM_eV = float(np.nanmin(Eband))
        E_CBM_meV = 1000.0 * E_CBM_eV

        if verbose:
            print(f"[info] Using EPW band {lcb} as CBM")
            print(f"[info] CBM (band {lcb}) = {E_CBM_eV:.6f} eV ({E_CBM_meV:.3f} meV)")

        rows = fh['/scattering/rows']
        # determine Nq from rows if user passed nq <= 0
        if nq is None or nq <= 0:
            try:
                max_iq = int(np.max(rows['iq_pos']))
                N_q = float(max_iq + 1)
            except Exception:
                N_q = float(101)
        else:
            N_q = float(nq)
        w_q = 1.0 / N_q
        if verbose:
            print(f"[info] Using q-weight w_q = 1/{int(N_q)} = {w_q:.9f}")

        # Accumulators
        rates_acc = {}
        counts_acc = {}
        Gamma_abs_acc = {}
        Gamma_em_acc = {}
        enk_sum = {}
        enkq_sum = {}
        imode_gamma_acc = {}

        total_rows = int(rows.shape[0])
        used_rows = 0
        rejected_rows = 0
        max_g_unfiltered = 0.0
        max_g_unfiltered_info = None
        max_g_used = 0.0
        max_g_used_info = None

        kBT_meV = kB_meV_per_K * T
        if verbose:
            print(f"[info] kBT = {kBT_meV:.3f} meV, window = {kbt_multiplier * kBT_meV:.3f} meV")
            print(f"[info] energy_tol (meV) = {energy_tol_meV:.6f}")

        gamma = float(gamma_meV)
        two_pi_over_hbar = 2.0 * math.pi / hbar

        allow_intra = mode_flag in ('intra', 'both')
        allow_inter_up = mode_flag in ('inter_up', 'both')
        allow_inter_down = mode_flag in ('inter_down', 'both')

        # iterate rows (now applying strict energy tolerance for ALL transitions)
        for r in rows:
            g_abs_meV = float(r['g_abs_meV'])
            if g_abs_meV > max_g_unfiltered:
                max_g_unfiltered = g_abs_meV
                max_g_unfiltered_info = (
                    int(r['ik_pos']), int(r['iq_pos']), int(r['imode']),
                    float(r['enk_eV']), float(r['enkq_eV']), float(r['omega_meV']),
                    int(r['ibnd']), int(r['jbnd'])
                )
            ib = int(r['ibnd']); jb = int(r['jbnd'])
            # Only conduction manifold
            if ib < lcb or jb < lcb:
                rejected_rows += 1; continue
            # identify transitions
            is_intra = (ib == jb)
            is_inter_up = (ib == lcb and jb > ib)
            is_inter_down = (ib > jb and jb == lcb)
            if is_intra and not allow_intra:
                rejected_rows += 1; continue
            if is_inter_up and not allow_inter_up:
                rejected_rows += 1; continue
            if is_inter_down and not allow_inter_down:
                rejected_rows += 1; continue

            ik = int(r['ik_pos']); iq = int(r['iq_pos']); im = int(r['imode'])
            enk_eV = float(r['enk_eV']); enkq_eV = float(r['enkq_eV']); omega_meV = float(r['omega_meV'])
            if g_abs_meV > max_g_used:
                max_g_used = g_abs_meV
                max_g_used_info = (ik, iq, im, enk_eV, enkq_eV, omega_meV, ib, jb)

            # Energy window (in meV from CBM)
            enk_meV = 1000.0 * enk_eV
            if abs(enk_meV - E_CBM_meV) > (kbt_multiplier * kBT_meV):
                rejected_rows += 1; continue

            # Bose factor
            if omega_meV <= 0:
                n_q = 0.0
            else:
                x = omega_meV / kBT_meV
                n_q = 1.0 / (math.exp(x) - 1.0) if x < 700 else 0.0

            # Lorentz deltas
            enkq_meV = 1000.0 * enkq_eV
            delta_abs = enkq_meV - enk_meV - omega_meV
            delta_em = enkq_meV - enk_meV + omega_meV

            # *** STRICT ENERGY TOLERANCE FOR ALL TRANSITIONS (Option A) ***
            if (abs(delta_abs) > energy_tol_meV) and (abs(delta_em) > energy_tol_meV):
                rejected_rows += 1; continue

            # Golden rule (units EXACTLY as in mc_rate_intraband.py)
            g_J = g_abs_meV * meV_to_J
            prefactor = two_pi_over_hbar * (g_J * g_J)

            L_abs = lorentzian_meV(delta_abs, gamma) / meV_to_J
            L_em = lorentzian_meV(delta_em, gamma) / meV_to_J

            Gamma_abs = prefactor * (n_q + 1.0) * L_abs
            Gamma_em = prefactor * n_q * L_em
            Gamma = Gamma_abs + Gamma_em

            # Weight
            Gamma_w = Gamma * w_q
            Gamma_abs_w = Gamma_abs * w_q
            Gamma_em_w = Gamma_em * w_q

            if is_intra:
                typ = "intra"
            elif is_inter_up:
                typ = "inter_up"
            elif is_inter_down:
                typ = "inter_down"
            else:
                typ = "other"

            key = (ik, ib, jb, typ)
            rates_acc[key] = rates_acc.get(key, 0.0) + Gamma_w
            Gamma_abs_acc[key] = Gamma_abs_acc.get(key, 0.0) + Gamma_abs_w
            Gamma_em_acc[key] = Gamma_em_acc.get(key, 0.0) + Gamma_em_w
            counts_acc[key] = counts_acc.get(key, 0) + 1
            enk_sum[key] = enk_sum.get(key, 0.0) + enk_eV
            enkq_sum[key] = enkq_sum.get(key, 0.0) + enkq_eV

            if key not in imode_gamma_acc:
                imode_gamma_acc[key] = {}
            imode_gamma_acc[key][im] = imode_gamma_acc[key].get(im, 0.0) + Gamma_w

            used_rows += 1

        # build results
        results = []
        for key in sorted(rates_acc.keys(), key=lambda x: (x[0], x[3], x[1], x[2])):
            ik, ib, jb, typ = key
            Gamma = rates_acc[key]
            tau = 1.0 / Gamma if Gamma > 0 else float("inf")
            E_k_meV = energies[energy_idx, ik] * 1000.0
            avg_enk = enk_sum[key] / counts_acc[key]
            avg_enkq = enkq_sum[key] / counts_acc[key]
            Gabs = Gamma_abs_acc[key]
            Gem = Gamma_em_acc[key]
            # event classification
            if (Gabs + Gem) == 0:
                event = "both"
            else:
                diff = abs(Gabs - Gem)
                if diff < 0.05 * (Gabs + Gem):
                    event = "both"
                elif Gabs > Gem:
                    event = "abs"
                else:
                    event = "em"
            im_dict = imode_gamma_acc.get(key, {})
            dominant_imode = max(im_dict, key=lambda m: im_dict[m]) if len(im_dict) > 0 else -1
            results.append((
                ik, E_k_meV, Gamma, tau,
                Gabs, Gem, counts_acc[key],
                avg_enk, avg_enkq, ib, jb, typ, event, dominant_imode
            ))

        if verbose:
            print()
            print("[info] ----- SCATTERING ROW STATISTICS -----")
            print(f"[info] Total rows:     {total_rows}")
            print(f"[info] Used rows:      {used_rows}")
            print(f"[info] Rejected rows:  {total_rows - used_rows}")
            print(f"[info] Reject %:       {100*(total_rows-used_rows)/total_rows:.2f}%")
            print()

        # write max_g_report exactly same format
        with open("max_g_report.txt", "w") as fh:
            fh.write(f"max_g_unfiltered_meV = {max_g_unfiltered:.6e}\n")
            if max_g_unfiltered_info:
                fh.write(str(max_g_unfiltered_info) + "\n")
            fh.write(f"\nmax_g_used_meV = {max_g_used:.6e}\n")
            if max_g_used_info:
                fh.write(str(max_g_used_info) + "\n")

        meta = {'lcb': lcb, 'band_idx': energy_idx, 'Eband': Eband, 'Nq': int(N_q)}
        return results, meta

# ---------- Build per-row channels with *the same* filters/formulas ----------
def build_channels_from_rows(h5path, lcb, gamma_meV, T, kbt_multiplier, energy_tol_meV, nq, mode_flag, verbose=True):
    with h5py.File(h5path,'r') as fh:
        rows = fh['/scattering/rows'][...]
        energies = fh['/bands/energies_eV'][...]
    if nq is None or nq <= 0:
        Nq = int(np.max(rows['iq_pos'])) + 1
    else:
        Nq = int(nq)
    w_q = 1.0 / float(Nq)
    kBT_meV = kB_meV_per_K * T
    gamma = float(gamma_meV)
    two_pi_over_hbar = 2.0 * math.pi / hbar

    ch_ik=[]; ch_iq=[]; ch_imode=[]; ch_omega=[]; ch_gabs=[]; ch_Gamma_w=[]
    for r in rows:
        g_abs_meV = float(r['g_abs_meV'])
        ib=int(r['ibnd']); jb=int(r['jbnd'])
        if ib < lcb or jb < lcb:
            continue
        is_intra = (ib==jb)
        is_inter_up = (ib==lcb and jb>ib)
        is_inter_down = (ib>jb and jb==lcb)
        if is_intra and mode_flag not in ('intra','both'): continue
        if is_inter_up and mode_flag not in ('inter_up','both'): continue
        if is_inter_down and mode_flag not in ('inter_down','both'): continue

        ik=int(r['ik_pos']); iq=int(r['iq_pos']); im=int(r['imode'])
        enk_eV=float(r['enk_eV']); enkq_eV=float(r['enkq_eV']); omega_meV=float(r['omega_meV'])
        enk_meV = 1000.0*enk_eV
        if abs(enk_meV - (1000.0*np.nanmin(energies[lcb-1,:]))) > (kbt_multiplier * kBT_meV):
            continue
        if omega_meV <= 0:
            n_q = 0.0
        else:
            x = omega_meV / kBT_meV
            n_q = 1.0 / (math.exp(x) - 1.0) if x < 700 else 0.0
        enkq_meV = 1000.0*enkq_eV
        delta_abs = enkq_meV - enk_meV - omega_meV
        delta_em  = enkq_meV - enk_meV + omega_meV
        if (abs(delta_abs) > energy_tol_meV) and (abs(delta_em) > energy_tol_meV):
            continue
        gJ = g_abs_meV * meV_to_J
        pref = two_pi_over_hbar * (gJ*gJ)
        L_abs = lorentzian_meV(delta_abs, gamma) / meV_to_J
        L_em  = lorentzian_meV(delta_em,  gamma) / meV_to_J
        Gamma_abs = pref * (n_q + 1.0) * L_abs
        Gamma_em  = pref * n_q * L_em
        Gamma = Gamma_abs + Gamma_em
        Gamma_w = Gamma * w_q
        ch_ik.append(ik); ch_iq.append(iq); ch_imode.append(im)
        ch_omega.append(omega_meV); ch_gabs.append(g_abs_meV); ch_Gamma_w.append(Gamma_w)

    if len(ch_ik)==0:
        raise RuntimeError("No channels found after filtering.")
    ch_ik = np.array(ch_ik, dtype=np.int32)
    ch_iq = np.array(ch_iq, dtype=np.int32)
    ch_imode = np.array(ch_imode, dtype=np.int32)
    ch_omega = np.array(ch_omega, dtype=np.float64)
    ch_gabs = np.array(ch_gabs, dtype=np.float64)
    ch_Gamma_w = np.array(ch_Gamma_w, dtype=np.float64)
    order = np.argsort(ch_ik, kind='stable')
    ch_ik = ch_ik[order]; ch_iq = ch_iq[order]; ch_imode = ch_imode[order]
    ch_omega = ch_omega[order]; ch_gabs = ch_gabs[order]; ch_Gamma_w = ch_Gamma_w[order]
    unique_iks, counts = np.unique(ch_ik, return_counts=True)
    Nk = int(energies.shape[1])
    ch_start = np.full((Nk,), -1, dtype=np.int32); ch_len = np.zeros((Nk,), dtype=np.int32)
    pos = 0
    for ik, c in zip(unique_iks, counts):
        ch_start[ik] = pos; ch_len[ik] = c; pos += c
    ch_cumsum = np.zeros_like(ch_Gamma_w)
    for ik in unique_iks:
        s = ch_start[ik]; L = ch_len[ik]
        if L>0:
            ch_cumsum[s:s+L] = np.cumsum(ch_Gamma_w[s:s+L])
    return {
        'ch_ik': ch_ik, 'ch_iq': ch_iq, 'ch_imode': ch_imode,
        'ch_omega': ch_omega, 'ch_gabs': ch_gabs, 'ch_Gamma_w': ch_Gamma_w,
        'ch_start': ch_start, 'ch_len': ch_len, 'ch_cumsum': ch_cumsum
    }

# ---------- simple selection kernel (identical logic) ----------
if have_numba:
    @njit
    def select_channels_numba(particles_ik, to_scatter_idx, rand_r, ch_start, ch_len, ch_cumsum, ch_iq):
        m = to_scatter_idx.shape[0]
        chosen_iq = np.empty((m,), dtype=np.int32)
        for ii in range(m):
            p = to_scatter_idx[ii]
            ik = particles_ik[p]
            s = ch_start[ik]; L = ch_len[ik]
            if L <= 0:
                chosen_iq[ii] = -1; continue
            target = rand_r[ii]
            chosen = -1
            for j in range(L):
                idx = s + j
                if target <= ch_cumsum[idx]:
                    chosen = idx; break
            if chosen == -1:
                chosen = s + L - 1
            chosen_iq[ii] = ch_iq[chosen]
        return chosen_iq
else:
    def select_channels_python(particles_ik, to_scatter_idx, rand_r, ch_start, ch_len, ch_cumsum, ch_iq):
        m = to_scatter_idx.shape[0]
        chosen_iq = np.empty((m,), dtype=np.int32)
        for ii in range(m):
            p = int(to_scatter_idx[ii])
            ik = int(particles_ik[p])
            s = int(ch_start[ik]); L = int(ch_len[ik])
            if L <= 0:
                chosen_iq[ii] = -1; continue
            target = rand_r[ii]
            chosen = -1
            for j in range(L):
                idx = s + j
                if target <= ch_cumsum[idx]:
                    chosen = idx; break
            if chosen == -1:
                chosen = s + L - 1
            chosen_iq[ii] = int(ch_iq[chosen])
        return chosen_iq

# ---------- MC runner (vectorized, KDTree) ----------
def run_mc(channels, energies_band_eV, velocities_band, kpts_cart, ikq_map,
           npart=200, steps=2000, dt=1e-15, E_field=1e4, T=300.0, seed=None,
           scatter_log_fh=None, max_log_events=None, verbose=False,
           energy_out=None):

    # allow energy diagnostics
    if energy_out is None:
        energy_out = {}

    if seed is not None:
        np.random.seed(seed)

    Nk = kpts_cart.shape[0]

    # --- initial distribution (Boltzmann) ---
    E_CBM_eV = float(np.nanmin(energies_band_eV))
    energies_rel = energies_band_eV - E_CBM_eV
    beta = 1.0 / (kB_eV_per_K * T)
    probs = np.exp(-energies_rel * beta)
    probs /= np.sum(probs)

    init_iks = np.random.choice(Nk, size=npart, p=probs)
    part_iks = init_iks.astype(np.int32)
    part_kcart = kpts_cart[part_iks].astype(np.float64)
    part_vel = velocities_band[part_iks].astype(np.float64)

    # scattering channel arrays
    ch_start  = channels['ch_start'].astype(np.int32)
    ch_len    = channels['ch_len'].astype(np.int32)
    ch_cumsum = channels['ch_cumsum'].astype(np.float64)
    ch_iq     = channels['ch_iq'].astype(np.int32)

    # total Gamma per ik
    total_Gamma = np.zeros((Nk,), dtype=np.float64)
    for ik in range(Nk):
        s = ch_start[ik]
        L = ch_len[ik]
        if s >= 0 and L > 0:
            total_Gamma[ik] = float(np.sum(channels['ch_Gamma_w'][s:s+L]))
        else:
            total_Gamma[ik] = 0.0

    # k-step
    if have_kdtree:
        kd = KDTree(kpts_cart)
    else:
        kd = None

    delta_k_per_dt = -e * np.array([E_field,0.0,0.0]) / hbar
    delta_k_step   = delta_k_per_dt * dt

    # storage
    v_drift_time = np.zeros((steps,), dtype=np.float64)
    logged = 0

    # energy accumulation (steady state)
    E_sum = 0.0
    E_count = 0
    final_E = None

    # --- time loop ---
    for t in range(steps):

        # drift
        part_kcart += delta_k_step

        # map to nearest k
        if kd is not None:
            _, new_iks = kd.query(part_kcart, k=1)
            part_iks = new_iks.astype(np.int32)
        else:
            d = ((part_kcart[:,None,:] - kpts_cart[None,:,:])**2).sum(axis=2)
            part_iks = np.argmin(d, axis=1).astype(np.int32)

        part_vel = velocities_band[part_iks]

        # scattering probability
        total_G = total_Gamma[part_iks]
        p_scatter = 1.0 - np.exp(-total_G * dt)

        rand_unif = np.random.random(npart)
        to_scatter = rand_unif < p_scatter

        if np.any(to_scatter):
            scatter_idx = np.nonzero(to_scatter)[0].astype(np.int32)
            rvals = np.random.random(scatter_idx.size) * total_G[scatter_idx]

            if have_numba:
                chosen_iqs = select_channels_numba(part_iks, scatter_idx, rvals,
                                                   ch_start, ch_len, ch_cumsum, ch_iq)
            else:
                chosen_iqs = select_channels_python(part_iks, scatter_idx, rvals,
                                                    ch_start, ch_len, ch_cumsum, ch_iq)

            for idx_local, p in enumerate(scatter_idx):
                iq = int(chosen_iqs[idx_local])
                if iq < 0:
                    continue
                old_ik = int(part_iks[p])
                if ikq_map is not None:
                    ik_prime = int(ikq_map[old_ik, iq])
                else:
                    ik_prime = old_ik

                # log
                if scatter_log_fh is not None and ((max_log_events is None) or logged < max_log_events):
                    k_old = kpts_cart[old_ik]; k_new = kpts_cart[ik_prime]
                    scatter_log_fh.write(
                        f"{E_field: .6e} {t:6d} {int(p):6d} {old_ik:6d} {iq:6d} {ik_prime:6d} "
                        f"{k_old[0]: .6e} {k_old[1]: .6e} {k_old[2]: .6e} "
                        f"{k_new[0]: .6e} {k_new[1]: .6e} {k_new[2]: .6e}\n"
                    )
                    logged += 1

                # update particle
                part_iks[p] = ik_prime
                part_kcart[p] = kpts_cart[ik_prime].copy()
                part_vel[p] = velocities_band[ik_prime].copy()

        # drift velocity
        mean_v = np.mean(part_vel, axis=0)
        v_drift_time[t] = float(mean_v[0])

        # -------------- ENERGY DIAGNOSTICS ----------------
        inst_energy_eV = energies_band_eV[part_iks] - E_CBM_eV

        # accumulate steady-state (last 30%)
        if t >= int(0.7 * steps):
            E_sum += float(np.sum(inst_energy_eV))
            E_count += inst_energy_eV.size

        # store final-step distribution
        if t == steps - 1:
            final_E = inst_energy_eV.copy()
            # --- Patch A: capture final particle k-indices for k-cloud output ---
            final_part_iks = part_iks.copy()

    # drift velocity result
    start_avg = int(0.7 * steps)
    mean_vd = float(np.mean(v_drift_time[start_avg:]))

    # populate energy_out
    if energy_out is not None:
        energy_out['E_sum'] = E_sum
        energy_out['E_count'] = E_count
        energy_out['final_E'] = final_E
        # --- Patch A (continued): attach final_part_iks for main() to process ---
        # If final_part_iks isn't defined (e.g., zero steps), set None
        if 'final_part_iks' in locals():
            energy_out['final_part_iks'] = final_part_iks
        else:
            energy_out['final_part_iks'] = None

    return mean_vd, logged


# ---------- main ----------
def main():
    p = argparse.ArgumentParser()

    p.add_argument('--h5', default='epw_quadrupole_parsed.h5')                  # Input parsed EPW HDF5
    p.add_argument('--gamma', type=float, default=10.0)                         # Lorentzian broadening (meV)
    p.add_argument('--T', type=float, default=77.0)                            # Temperature (K)
    p.add_argument('--kbt-mult', type=float, default=10.0)                      # Thermal window multiplier
    p.add_argument('--energy-tol', type=float, default=50)                      # Energy tolerance (meV)
    p.add_argument('--nq', type=int, default=None)                              # Use first N q-points
    p.add_argument('--npart', type=int, default=500)                          # Number of MC particles
    p.add_argument('--steps', type=int, default=10000)                          # Number of MC time-steps
    p.add_argument('--dt', type=float, default=0.01e-16)                        # Time-step (s)
    p.add_argument('--field-start', type=float, default=1e5)                    # Start electric field (V/m)
    p.add_argument('--field-stop', type=float, default=1e8)                     # Stop electric field (V/m)
    p.add_argument('--nfields', type=int, default=8)                            # Number of fields (log-spaced)
    p.add_argument('--seed', type=int, default=None)                            # Random seed
    p.add_argument('--scatter-log', default='scatter_events.txt')               # Scatter-event logfile
    p.add_argument('--max-log-events', type=int, default=10000)                  # Max number of logged events
    p.add_argument('--run-checks', action='store_true')                         # Run stability diagnostics

    args = p.parse_args()

    if args.seed is not None: np.random.seed(args.seed)
    if not os.path.exists(args.h5):
        print("[error] HDF5 not found:", args.h5); sys.exit(1)

    print("[info] Computing aggregated rates...")
    results, meta = compute_rates(args.h5, gamma_meV=args.gamma, T=args.T,
                                  kbt_multiplier=args.kbt_mult, cbm_ref_eV=None,
                                  nq=(args.nq if args.nq is not None else 0),
                                  mode_flag='both', energy_tol_meV=args.energy_tol, verbose=True)
    out_rates = 'scattering_rates.txt'
    with open(out_rates,'w') as fo:
        fo.write("#Intraband/interband scattering results (conduction manifold only). Event = abs/em/both.\n")
        fo.write("#---------------------------------------------------------------------------------------------------------------------------------------------------\n")
        fo.write(f"{'#'}{'ik_pos':>6}  {'E_k_meV':>12}  {'Gamma_s':>14}  {'tau_s':>14}  {'Gamma_abs_s':>14}  {'Gamma_em_s':>14}  {'n_rows':>7}  {'avg_enk_eV':>12}  {'avg_enkq_eV':>12}  {'ibnd':>5}  {'jbnd':>5}  {'type':>10}  {'event':>7}  {'imode':>7}\n")
        fo.write("#---------------------------------------------------------------------------------------------------------------------------------------------------\n")
        for (ik, Ek, G, tau, Gabs, Gem, nrows, avg_enk, avg_enkq, ib, jb, typ, event, imode) in results:
            fo.write(f"{ik:6d}  {Ek:12.6f}  {G:14.6e}  {tau:14.6e}  {Gabs:14.6e}  {Gem:14.6e}  {nrows:7d}  {avg_enk:12.6f}  {avg_enkq:12.6f}  {ib:5d}  {jb:5d}  {typ:>10s}  {event:>7s}  {imode:7d}\n")
    print(f"[ok] Wrote {out_rates}")

    print("[info] Building per-row channels (same filters/formulas)...")
    channels = build_channels_from_rows(args.h5, meta['lcb'], args.gamma, args.T, args.kbt_mult,
                                        (args.energy_tol if args.energy_tol is not None else args.gamma),
                                        args.nq, 'both', verbose=True)
    print("[ok] Built channels.")

    with h5py.File(args.h5,'r') as fh:
        kpts_cart = fh['/scattering/kpts_cart'][...]
        try:
            velocities_all = fh['/bands/velocities_m_s'][...]
        except Exception:
            raise RuntimeError("Velocities missing in HDF5")
        band_list = fh['/bands/band_list'][...]
        band_idx = meta['band_idx'] if 'band_idx' in meta else 0
        energies = fh['/bands/energies_eV'][...]
        ikq_map = fh['/scattering/ikq_map'][...] if '/scattering/ikq_map' in fh else None

    if velocities_all.ndim == 3 and velocities_all.shape[0] > band_idx:
        velocities_band = velocities_all[band_idx, :, :]
    elif velocities_all.ndim == 2:
        velocities_band = velocities_all
    else:
        raise RuntimeError("Velocities shape not recognised.")
    energies_band = energies[band_idx, :]

    # ---------- RUN-TIME CHECKS (global) ----------
    if args.run_checks:
        print("[check] Running global stability checks...")
        ch_start = channels['ch_start']; ch_len = channels['ch_len']; ch_Gamma = channels['ch_Gamma_w']
        Nk = kpts_cart.shape[0]
        total_Gamma = np.zeros((Nk,), dtype=np.float64)
        for ik in range(Nk):
            s = int(ch_start[ik])
            L = int(ch_len[ik])
            if s >= 0 and L > 0:
                total_Gamma[ik] = float(np.sum(ch_Gamma[s:s+L]))
        p = 1.0 - np.exp(-total_Gamma * args.dt)
        print(f"[check] p_scatter mean = {np.mean(p):.3e}, median = {np.median(p):.3e}, max = {np.max(p):.3e}")
        if np.max(p) > 0.3:
            print("[warn] max p_scatter > 0.3. Reduce dt for stability.")
        nonzero = np.sum(ch_len > 0)
        print(f"[check] k-points with channels: {nonzero}/{Nk} ({100.0*nonzero/Nk:.1f}%)")
        if nonzero < 0.7 * Nk:
            print("[warn] Less than 70% k-points have channels. Results may be partially ballistic.")
        tau = np.where(total_Gamma > 0, 1.0/total_Gamma, np.inf)
        finite_tau = tau[np.isfinite(tau)]
        if finite_tau.size > 0:
            pcts = np.percentile(finite_tau, [1,10,50,90,99])
            print("[check] tau percentiles (s): 1% {0:.3e}, 10% {1:.3e}, 50% {2:.3e}, 90% {3:.3e}, 99% {4:.3e}".format(*pcts))
            if pcts[0] < 1e-16:
                print("[warn] Some τ < 1e-16 s (too fast). Check gamma/energy_tol/weights.")
        else:
            print("[warn] No finite lifetimes found (total_Gamma == 0 everywhere).")
        try:
            if have_kdtree:
                kd_tmp = KDTree(kpts_cart)
                dists, idxs = kd_tmp.query(kpts_cart, k=2)
                min_spacing = np.min(dists[:,1])
            else:
                diffs = np.sqrt(np.sum((kpts_cart - kpts_cart.mean(axis=0))**2, axis=1))
                min_spacing = np.min(np.abs(np.diff(np.sort(diffs))))
        except Exception:
            min_spacing = None
        delta_k = np.linalg.norm((-e * np.array([args.field_start,0.0,0.0]) / hbar) * args.dt)
        if min_spacing is not None:
            print(f"[check] Δk per dt = {delta_k:.3e}, min grid spacing ≈ {min_spacing:.3e}")
            if delta_k > 0.2 * min_spacing:
                print("[warn] Δk per dt > 0.2 * min_spacing; consider reducing dt or using finer k-grid.")
        else:
            print(f"[check] Δk per dt = {delta_k:.3e} (min spacing not computed)")
        try:
            small_npart = min(200, max(50, int(0.001*args.npart)))
            small_steps = min(500, max(100, int(0.005*args.steps)))
            print(f"[check] quick repeatability test: npart={small_npart}, steps={small_steps}")
            vd1, _ = run_mc(channels, energies_band, velocities_band, kpts_cart, ikq_map,
                            npart=small_npart, steps=small_steps, dt=args.dt, E_field=args.field_start, T=args.T,
                            seed=(args.seed if args.seed is not None else 12345), scatter_log_fh=None, max_log_events=0)
            vd2, _ = run_mc(channels, energies_band, velocities_band, kpts_cart, ikq_map,
                            npart=small_npart, steps=small_steps, dt=args.dt, E_field=args.field_start, T=args.T,
                            seed=(None if args.seed is not None else 54321), scatter_log_fh=None, max_log_events=0)
            print(f"[check] quick vd repeat: vd1={vd1:.3e}, vd2={vd2:.3e}, diff={abs(vd1-vd2):.3e}")
            if abs(vd1 - vd2) > 0.1 * max(abs(vd1), 1e-12):
                print("[warn] Repeatability difference >10%: consider increasing npart or steps.")
        except Exception as ex:
            print("[warn] quick repeatability test failed:", ex)

    scatter_fh = open(args.scatter_log, 'w')
    scatter_fh.write("# field[V/m] timestep particle old_ik iq ik_prime k_old_x k_old_y k_old_z k_new_x k_new_y k_new_z\n")

    fields = np.logspace(np.log10(args.field_start), np.log10(args.field_stop), args.nfields)
    vd_list = []
    energy_list = []          # <-- added
    energy_std_list = []      # <-- added
    total_logged = 0
    t0_all = time.time()

    for E_field in fields:
        print(f"[info] Running MC for E = {E_field:.3e} V/m")

        if args.run_checks:
            delta_k = np.linalg.norm((-e * np.array([E_field,0.0,0.0]) / hbar) * args.dt)
            if 'min_spacing' in locals() and min_spacing is not None:
                if delta_k > 0.2*min_spacing:
                    print(f"[warn] Field {E_field:.3e}: Δk per dt ({delta_k:.3e}) > 0.2*min_spacing ({0.2*min_spacing:.3e})")
            if 'total_Gamma' in locals():
                p = 1.0 - np.exp(-total_Gamma * args.dt)
                print(f"[check] Field {E_field:.3e}: p_scatter mean={np.mean(p):.3e}, max={np.max(p):.3e}")

        t0 = time.time()

        # ----- new dict to collect energy from run_mc -----
        energy_out = {}

        vd, logged = run_mc(channels, energies_band, velocities_band, kpts_cart, ikq_map,
                            npart=args.npart, steps=args.steps, dt=args.dt, E_field=E_field, T=args.T,
                            seed=args.seed, scatter_log_fh=scatter_fh,
                            max_log_events=args.max_log_events, verbose=True,
                            energy_out=energy_out)

        t1 = time.time()
        print(f"  -> vd = {vd:.6e} m/s  (time {t1-t0:.2f}s)  logged={logged}")

        vd_list.append(vd)
        total_logged += logged

        # ----- finalize energy -----
        if ('E_sum' in energy_out) and (energy_out['E_count'] > 0):
            avg_E = energy_out['E_sum'] / energy_out['E_count']
        else:
            avg_E = float('nan')

        if 'final_E' in energy_out:
            std_E = float(np.std(energy_out['final_E']))
        else:
            std_E = float('nan')

        energy_list.append(avg_E)
        energy_std_list.append(std_E)

        # --- Patch B: write k-cloud per-field file (ik_index, energy_eV, N_particles)
        if ('final_part_iks' in energy_out) and (energy_out['final_part_iks'] is not None):
            final_iks = np.array(energy_out['final_part_iks'], dtype=np.int64)
            Nk_here = energies_band.shape[0]
            counts = np.bincount(final_iks, minlength=Nk_here)
            fname = f"kcloud_E_{E_field:.3e}.txt"
            with open(fname, 'w') as fhk:
                fhk.write("# ik_index   energy_eV   N_particles\n")
                for ik in range(Nk_here):
                    e_val = float(energies_band[ik]) if not np.isnan(energies_band[ik]) else float('nan')
                    cnt = int(counts[ik])
                    fhk.write(f"{ik:6d}  {e_val:12.6f}  {cnt:8d}\n")
            print(f"[ok] Wrote k-cloud file: {fname}")

        scatter_fh.flush()

    scatter_fh.close()
    t1_all = time.time()
    print(f"[info] Total MC time: {t1_all-t0_all:.2f}s")

    if len(fields) >= 2:
        mu_low = (vd_list[1] - vd_list[0]) / (fields[1] - fields[0])
        mu_high = (vd_list[-1] - vd_list[-2]) / (fields[-1] - fields[-2])
    else:
        mu_low = 0.0; mu_high = 0.0

    print(f"[info] Low-field mobility (Δvd/ΔE) = {mu_low:.6e} m^2/Vs = {mu_low*1e4:.6e} cm^2/Vs")
    print(f"[info] High-field mobility (Δvd/ΔE) = {mu_high:.6e} m^2/Vs = {mu_high*1e4:.6e} cm^2/Vs")

    # ---------- rewritten drift_results.txt with energy ----------
    with open('drift_results.txt','w') as fh:
        fh.write('# Field(V/m)        vd(m/s)        AvgEnergy(eV_rel_CBM)        StdEnergy(eV_rel_CBM)\n')
        for F, vd, Eavg, Estd in zip(fields, vd_list, energy_list, energy_std_list):
            fh.write(f"{F:15.6e} {vd:15.6e} {Eavg:15.6e} {Estd:15.6e}\n")
        fh.write(f"\n# Low-field mobility (Δvd/ΔE): {mu_low:.6e} m^2/Vs = {mu_low*1e4:.6e} cm^2/Vs\n")
        fh.write(f"# High-field mobility (Δvd/ΔE): {mu_high:.6e} m^2/Vs = {mu_high*1e4:.6e} cm^2/Vs\n")

    print("[ok] Wrote drift_results.txt and scatter log.")

if __name__ == "__main__":
    main()

