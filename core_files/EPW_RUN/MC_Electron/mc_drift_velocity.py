#!/usr/bin/env python3
"""
mc_drift_epw.py

Clean Monte Carlo drift driver matching EPW SERTA scattering form (no 1/|∇E|,
no BZ area prefactor). Use --broadening to set Lorentzian HWHM (meV).

Usage example:
  python mc_drift_epw.py --h5 epw_quadrupole_parsed.h5 --Ex 1e4 --dt 1e-15 \
      --nsteps 5000 --npart 2000 --broadening 20.0 --rate-scale 1.0
"""
import argparse
import math
import os
import sys
import numpy as np
import h5py
import time
from collections import defaultdict

# Physical constants
eV_to_J = 1.602176634e-19
hbar = 1.054571817e-34
q_e = 1.602176634e-19
kb = 1.380649e-23
pi = math.pi

# -------------------------
# Helpers
# -------------------------
def lorentzian_J(delta_J, gamma_J):
    """Normalized Lorentzian in energy (J^-1) with HWHM = gamma_J."""
    return (gamma_J / pi) / (delta_J*delta_J + gamma_J*gamma_J)

def bose_factor(omega_J, T):
    if omega_J <= 0.0 or T <= 0.0:
        return 0.0
    x = omega_J / (kb * T)
    if x > 700:
        return 0.0
    return 1.0 / (math.exp(x) - 1.0)

def frac_wrap(frac):
    """Wrap fractional coordinates into [0,1)"""
    return frac - np.floor(frac)

# -------------------------
# Scattering (EPW SERTA form)
# -------------------------
def calc_row_rate_SERTA(row, T, broadening_meV, rate_scale=1.0, max_rate=1e15):
    """
    EPW SERTA form: Gamma = (2pi/hbar) * g^2 * [ (Nq+1)L(Δ-ω) + Nq L(Δ+ω) ]
    Inputs:
      - row: structured row with fields (enk_eV, enkq_eV, omega_meV, g_abs_meV)
      - T: temperature (K)
      - broadening_meV: HWHM in meV
      - rate_scale: multiplicative scale applied to final rate (for tuning)
    Returns: rate in s^-1
    """
    # get g (meV -> eV -> J)
    g_meV = float(row['g_abs_meV'])
    g_J = (g_meV * 1e-3) * eV_to_J
    g2 = g_J * g_J
    if g2 <= 0.0:
        return 0.0

    # phonon energy
    omega_meV = float(row['omega_meV'])
    omega_J = (omega_meV * 1e-3) * eV_to_J

    # electron energies
    Ei_J = float(row['enk_eV']) * eV_to_J
    Ef_J = float(row['enkq_eV']) * eV_to_J

    # energy differences (J)
    Delta_em = Ef_J - Ei_J - omega_J  # for emission resonance: Δ ≈ 0 when Ef = Ei - ω
    Delta_ab = Ef_J - Ei_J + omega_J  # for absorption resonance

    gamma_J = (broadening_meV * 1e-3) * eV_to_J

    Nq = bose_factor(omega_J, T)

    pref = (2.0 * pi / hbar) * g2

    L_em = lorentzian_J(Delta_em, gamma_J)
    L_ab = lorentzian_J(Delta_ab, gamma_J)

    rate = pref * ((Nq + 1.0) * L_em + (Nq) * L_ab)
    rate *= float(rate_scale)

    if not np.isfinite(rate) or rate < 0.0:
        return 0.0
    return float(min(rate, max_rate))


# -------------------------
# k/q mapping helpers
# -------------------------
def find_nearest_k_index(kcoords_frac, target_frac):
    """
    Find index of kcoords_frac (Nk x 3) nearest to target_frac (3,) in cartesian
    distance on the unit cell (wrap fractional coordinates).
    This uses minimum image in fractional space approximated in Cartesian by mapping
    virtual distances on the torus.
    """
    # compute wrapped delta: choose minimal periodic image per component
    # for robust distance we convert to simple euclidean with wrap in fractional coords
    diffs = np.abs(kcoords_frac - target_frac)
    diffs = np.minimum(diffs, 1.0 - diffs)
    d2 = np.sum(diffs * diffs, axis=1)
    return int(np.nanargmin(d2))


# -------------------------
# Interpolation & masses (kept simple)
# -------------------------
def interp_state_quantities(energies, velocities, kcoords_frac, kpos_float):
    nb, Nk = energies.shape
    if kpos_float <= 0:
        il = 0
        w = 0.0
    elif kpos_float >= Nk - 1:
        il = Nk - 2
        w = 1.0
    else:
        il = int(math.floor(kpos_float))
        w = kpos_float - il
    ir = il + 1
    E_il = energies[:, il]
    E_ir = energies[:, ir]
    V_il = velocities[:, il, :]
    V_ir = velocities[:, ir, :]
    E_interp = np.where(np.isnan(E_il), E_ir, np.where(np.isnan(E_ir), E_il, (1.0 - w) * E_il + w * E_ir))
    V_interp = np.where(np.isnan(V_il), V_ir, np.where(np.isnan(V_ir), V_il, (1.0 - w) * V_il + w * V_ir))
    return E_interp, V_interp


# -------------------------
# Build scattering table
# -------------------------
def build_scattering_table(rows_array, velocities, band_list, kcoords_frac, qcoords_frac,
                           T=300.0, broadening_meV=20.0, scattering_mode='both', rate_scale=1.0):
    """
    Build scatter_table: dict ( (ibnd, ik_pos) -> list of (row_index, rate) )
    and totals dict with total out rates per (ibnd, ik_pos).
    rows_array expected fields: ik_pos, iq_pos, ibnd, jbnd, omega_meV, g_abs_meV, enk_eV, enkq_eV
    """
    Nrows = rows_array.shape[0]
    table = {}
    totals = defaultdict(float)

    # Precompute rates per row
    rates = np.empty(Nrows, dtype=float)
    for i in range(Nrows):
        r = rows_array[i]
        # energy-conservation gating (optional): we still compute full SERTA but gating small tails could be used
        rates[i] = calc_row_rate_SERTA(r, T, broadening_meV, rate_scale=rate_scale)

    # group by (ibnd, ik_pos)
    # use lexsort for grouping
    ib_arr = rows_array['ibnd'].astype(int)
    ik_arr = rows_array['ik_pos'].astype(int)
    order = np.lexsort((ik_arr, ib_arr))
    ibs = ib_arr[order]
    iks = ik_arr[order]
    rates_s = rates[order]
    rows_s = rows_array[order]

    change = np.empty_like(order, dtype=bool)
    change[0] = True
    change[1:] = (ibs[1:] != ibs[:-1]) | (iks[1:] != iks[:-1])
    starts = np.nonzero(change)[0]
    ends = np.append(starts[1:], len(order))

    for s, e in zip(starts, ends):
        ib = int(ibs[s])
        ik = int(iks[s])
        key = (ib, ik)
        seg_rows = rows_s[s:e]
        seg_rates = rates_s[s:e]
        table[key] = [(int(seg_rows[j].tolist().__len__()), float(seg_rates[j])) for j in range(len(seg_rates))]  # placeholder rows not used directly
        totals[key] = float(seg_rates.sum())

    # Note: we keep the original rows_array available; table stores indices not direct row objects (memory friendly)
    # For simplicity in the driver, we'll keep rows_array and use row indices returned from selection to fetch details.

    return table, totals, rows_array, order  # return order for consistent indexing


# -------------------------
# MC Engine
# -------------------------
class MCDriftEPW:
    def __init__(self, parsed, charge_type='electron', mode='B', scattering_choice='both',
                 T=300.0, dt=1e-15, grid_xy=(12,12), broadening_meV=20.0, rate_scale=1.0):
        self.parsed = parsed
        self.energies = parsed['energies']     # (nb, Nk)
        self.velocities = parsed['velocities'] # (nb, Nk, 3)
        self.band_list = parsed['band_list']   # array of band numbers
        self.kcoords = parsed['kcoords']       # fractional (Nk,3)
        self.qcoords = parsed['qcoords']       # fractional (Nq,3)
        self.rows = parsed['rows']             # structured array with fields used
        self.ik_list = parsed['ik_list']
        self.iq_list = parsed['iq_list']

        self.charge = -q_e if str(charge_type).lower().startswith('e') else +q_e
        self.charge_type = charge_type
        self.mode = mode
        self.scattering_choice = scattering_choice
        self.T = T
        self.dt = dt
        self.nb, self.Nk = self.energies.shape

        # broadening and scale
        self.broadening_meV = float(broadening_meV)
        self.rate_scale = float(rate_scale)

        # read meta: CBM/VBM/Fermi from parsed meta if available
        meta = parsed.get('meta', {})
        self.fermi = float(meta.get('E_fermi_eV', np.nan)) if meta else np.nan
        # parser stores ebndmin_eV and ebndmax_eV and ibndmin/ibndmax
        self.vbm_e = float(meta.get('ebndmin_eV', np.nan))
        self.cbm_e = float(meta.get('ebndmax_eV', np.nan))
        self.vbm_band = int(meta.get('ibndmin', -1)) if 'ibndmin' in meta else -1
        self.cbm_band = int(meta.get('ibndmax', -1)) if 'ibndmax' in meta else -1

        # build scattering table
        self.scatter_table, self.total_rate, self.rows_array, self.row_order = build_scattering_table(
            self.rows, self.velocities, self.band_list, self.kcoords, self.qcoords,
            T=self.T, broadening_meV=self.broadening_meV, scattering_mode=self.scattering_choice,
            rate_scale=self.rate_scale
        )

    def sample_initial_state(self, band='auto'):
        # pick band index (index into band_list array)
        if band == 'auto':
            if self.charge < 0:
                # electrons -> CBM band
                if self.cbm_band in self.band_list:
                    ib_idx = int(np.where(self.band_list == self.cbm_band)[0][0])
                else:
                    ib_idx = int(np.nanargmin(np.nanmin(self.energies, axis=1)))  # fallback
            else:
                # holes -> VBM band
                if self.vbm_band in self.band_list:
                    ib_idx = int(np.where(self.band_list == self.vbm_band)[0][0])
                else:
                    ib_idx = int(np.nanargmax(np.nanmax(self.energies, axis=1)))
        else:
            try:
                ib_idx = int(np.where(self.band_list == int(band))[0][0])
            except Exception:
                ib_idx = 0
        # choose kpos nearest to band minimum/maximum as appropriate
        if self.charge < 0:
            # electrons -> pick k with minimum energy
            try:
                kpos = int(np.nanargmin(self.energies[ib_idx, :]))
            except Exception:
                kpos = 0
        else:
            try:
                kpos = int(np.nanargmax(self.energies[ib_idx, :]))
            except Exception:
                kpos = 0
        return (int(self.band_list[ib_idx]), int(kpos))

    def accelerate_kpos(self, kpos_float, Ex):
        # semiclassical dk/dt = qE / ħ, but we operate on fractional k-path index float space.
        # We will map displacement in cartesian k to nearest point on the path (same logic as original code).
        dk_cart = (self.charge * Ex / hbar) * self.dt * np.array([1.0, 0.0, 0.0])  # m^-1
        # convert fractional path to cart via reciprocal vectors if present in parsed (parser stored reciprocal in /reciprocal_lattice/b_vectors)
        recip = self.parsed.get('recip', None)
        if recip is not None:
            kcoords_cart = np.dot(self.kcoords, recip)
        else:
            # fallback: use fractional coordinates as if cartesian (less accurate)
            kcoords_cart = self.kcoords.copy()

        Nk = self.Nk
        if kpos_float <= 0:
            il = 0
            w = 0.0
        elif kpos_float >= Nk - 1:
            il = Nk - 2
            w = 1.0
        else:
            il = int(math.floor(kpos_float))
            w = kpos_float - il
        ir = il + 1
        kcart_il = kcoords_cart[il]
        kcart_ir = kcoords_cart[ir]
        kcart_cur = (1.0 - w) * kcart_il + w * kcart_ir
        kcart_new = kcart_cur + dk_cart

        # find nearest segment projection same as original algorithm
        best_i = 0
        best_t = 0.0
        best_dist = 1e99
        for i in range(Nk - 1):
            a = kcoords_cart[i]
            b = kcoords_cart[i + 1]
            ab = b - a
            denom = np.dot(ab, ab)
            if denom == 0:
                t = 0.0
            else:
                t = np.dot(kcart_new - a, ab) / denom
                t = max(0.0, min(1.0, t))
            proj = a + t * ab
            d2 = np.sum((proj - kcart_new) ** 2)
            if d2 < best_dist:
                best_dist = d2
                best_i = i
                best_t = t
        new_kpos_float = best_i + best_t
        return max(0.0, min(self.Nk - 1.0, new_kpos_float))

    def select_scattering_event(self, ibnd, ik_pos):
        key = (ibnd, ik_pos)
        lst = self.scatter_table.get(key, [])
        if not lst:
            return None
        # lst currently is list of (row_index_placeholder, rate) but we will instead use direct rates array
        # For simplicity, reconstruct list of indices where rows_array matches key
        # This is slower but safe and clear.
        candidate_idxs = []
        candidate_rates = []
        # iterate through all rows (vectorization possible; keep simple)
        for idx in range(self.rows.shape[0]):
            r = self.rows[idx]
            if int(r['ibnd']) == ibnd and int(r['ik_pos']) == ik_pos:
                candidate_idxs.append(idx)
                candidate_rates.append(calc_row_rate_SERTA(r, self.T, self.broadening_meV, rate_scale=self.rate_scale))
        if len(candidate_rates) == 0:
            return None
        rates = np.array(candidate_rates, dtype=float)
        total = float(rates.sum())
        if total <= 0.0:
            return None
        xi = np.random.random() * total
        cum = np.cumsum(rates)
        idx_sel = int(np.searchsorted(cum, xi))
        if idx_sel >= len(candidate_idxs):
            idx_sel = len(candidate_idxs) - 1
        row_idx = candidate_idxs[idx_sel]
        return int(row_idx)

    def map_kplusq(self, ik_pos, iq_pos):
        """
        Map k' = k(ik_pos) + q(iq_pos) (fractional) and find nearest k index on kcoords.
        iq_pos may be -1 (invalid): return original ik_pos.
        """
        if iq_pos < 0 or iq_pos >= self.qcoords.shape[0]:
            return int(ik_pos)
        k_frac = self.kcoords[int(ik_pos)]
        q_frac = self.qcoords[int(iq_pos)]
        target = frac_wrap(k_frac + q_frac)
        new_idx = find_nearest_k_index(self.kcoords, target)
        return int(new_idx)

    def run_particle_ensemble(self, Ex, nsteps=10000, npart=1000, conv_tol=1e-6,
                              out_interval=100, log_path=None, log_freq=50, debug=False):
        # initialize particles
        particles = []
        for _ in range(npart):
            b, kp = self.sample_initial_state('auto')
            particles.append({'band': b, 'kpos': float(kp)})

        drift_hist = []
        t = 0.0
        v_d = 0.0

        if log_path is not None:
            with open(log_path, 'w') as lf:
                lf.write("===== MC DRIFT SIMULATION LOG =====\n")
                lf.write(f"Mode: {self.mode}\n")
                lf.write(f"Carrier: {'electron' if self.charge < 0 else 'hole'}\n")
                lf.write(f"Electric field Ex = {Ex:.6e} V/m\n")
                lf.write(f"Temperature = {self.T} K\n")
                lf.write(f"Time step dt = {self.dt} s\n")
                lf.write(f"Particles = {npart}\n")
                lf.write(f"Steps = {nsteps}\n")
                lf.write(f"broadening_meV = {self.broadening_meV}\n")
                lf.write(f"rate_scale = {self.rate_scale}\n")
                lf.write('-----------------------------------\n')
                lf.write('\n')

        total_scatter = 0
        emission_count = 0
        absorption_count = 0

        for step in range(nsteps):
            vx_acc = 0.0
            energies_acc = 0.0
            scat_rates_acc = 0.0

            for p in particles:
                kpos_new = self.accelerate_kpos(p['kpos'], Ex)
                p['kpos'] = kpos_new
                E_p, V_p = interp_state_quantities(self.energies, self.velocities, self.kcoords, p['kpos'])
                # find band index
                try:
                    ib_i = int(np.where(self.band_list == p['band'])[0][0])
                except Exception:
                    ib_i = 0
                vvec = V_p[ib_i]
                vx_acc += vvec[0]

                # kinetic energy referenced to CBM/VBM
                E_raw = E_p[ib_i]
                if self.charge < 0:
                    if not np.isnan(self.cbm_e):
                        Ekin = E_raw - self.cbm_e
                    elif not np.isnan(self.fermi):
                        Ekin = E_raw - self.fermi
                    else:
                        Ekin = E_raw
                else:
                    if not np.isnan(self.vbm_e):
                        Ekin = self.vbm_e - E_raw
                    elif not np.isnan(self.fermi):
                        Ekin = self.fermi - E_raw
                    else:
                        Ekin = -E_raw

                energies_acc += Ekin

                ik_nearest = int(round(p['kpos']))
                key = (p['band'], ik_nearest)
                Gamma = float(self.total_rate.get(key, 0.0))
                scat_rates_acc += Gamma
                p_scatter = 1.0 - math.exp(-Gamma * self.dt) if Gamma > 0 else 0.0
                if np.random.random() < p_scatter:
                    row_idx = self.select_scattering_event(p['band'], ik_nearest)
                    if row_idx is not None:
                        rrow = self.rows[row_idx]
                        new_band = int(rrow['jbnd'])
                        # map new kpos from iq -> k'
                        new_ikpos = self.map_kplusq(ik_nearest, int(rrow['iq_pos']))
                        # emission vs absorption classification
                        try:
                            dE = float(rrow['enk_eV']) - float(rrow['enkq_eV'])
                        except Exception:
                            dE = 0.0
                        if dE > 1e-12:
                            emission_count += 1
                        elif dE < -1e-12:
                            absorption_count += 1
                        total_scatter += 1
                        p['band'] = new_band
                        p['kpos'] = float(new_ikpos)

            v_d_instant = vx_acc / float(npart) if npart > 0 else 0.0
            avg_energy = (energies_acc / float(npart)) if npart > 0 else 0.0
            avg_scat_rate = (scat_rates_acc / float(npart)) if npart > 0 else 0.0
            drift_hist.append((t, v_d_instant))
            t += self.dt

            # logging
            if log_path is not None and log_freq > 0 and (step % log_freq == 0):
                converged_flag = False
                if len(drift_hist) >= log_freq + 1:
                    recent = np.mean([v for (_, v) in drift_hist[-log_freq:]])
                    if abs(recent - v_d) < conv_tol:
                        converged_flag = True
                with open(log_path, 'a') as lf:
                    lf.write(f"[step {step}]\n")
                    lf.write(f"   avg_vx      = {v_d_instant:.5e} m/s\n")
                    lf.write(f"   avg_energy  = {avg_energy * 1e3:.5e} meV\n")
                    lf.write(f"   scat_rate   = {avg_scat_rate:.5e} 1/s\n")
                    lf.write(f"   emission    = {emission_count}\n")
                    lf.write(f"   absorption  = {absorption_count}\n")
                    lf.write(f"   converged   = {converged_flag}\n")
                    lf.write('-----------------------------------\n')

            if step % out_interval == 0 and step > 0:
                recent = np.mean([v for (_, v) in drift_hist[-out_interval:]])
                if abs(recent - v_d) < conv_tol:
                    if log_path is not None:
                        with open(log_path, 'a') as lf:
                            lf.write(f"[convergence] drift velocity stable within {conv_tol} at step {step}\n\n")
                    break
                v_d = recent

        # final reporting
        if log_path is not None:
            with open(log_path, 'a') as lf:
                lf.write('===== FINAL RESULTS =====\n')
                last_v = drift_hist[-1][1] if len(drift_hist) > 0 else 0.0
                lf.write(f"Steady-state drift velocity = {last_v:.5e} m/s\n")
                mu = 0.0
                if Ex != 0:
                    mu = (last_v / Ex) * 1e4
                lf.write(f"Mobility = {mu:.5e} cm^2/V/s\n")
                lf.write(f"Total scattering events = {total_scatter}\n")
                lf.write(f"Total emission = {emission_count}\n")
                lf.write(f"Total absorption = {absorption_count}\n")
                lf.write('Simulation completed.\n')
                lf.write('====================================\n')

        return np.array(drift_hist)


# -------------------------
# CLI
# -------------------------
def parse_args():
    p = argparse.ArgumentParser(description='Run MC drift simulation (EPW SERTA scattering)')
    p.add_argument('--h5', required=True, help='parsed HDF5 file')
    p.add_argument('--Ex', type=float, required=True, help='Electric field along x (V/m)')
    p.add_argument('--T', type=float, default=300.0, help='Temperature (K)')
    p.add_argument('--dt', type=float, default=1e-15, help='Time step (s)')
    p.add_argument('--nsteps', type=int, default=20000, help='Number of time steps')
    p.add_argument('--npart', type=int, default=1000, help='Number of particles (mode B)')
    p.add_argument('--conv_tol', type=float, default=1e-6, help='Convergence tol on v_d')
    p.add_argument('--out', default='mc_out.npz', help='Output file prefix')
    p.add_argument('--scat', default='both', choices=['both','emission','absorption'], help='Scattering choice')
    p.add_argument('--log', default='mc_runtime_log.txt', help='Runtime log filename (written next to HDF5)')
    p.add_argument('--log_freq', type=int, default=50, help='Log update frequency (steps)')
    p.add_argument('--grid', nargs=2, type=int, default=[12,12], help='EPW sampling grid (Nkx Nky) used (not required)')
    p.add_argument('--broadening', type=float, default=20.0, help='Lorentzian HWHM in meV (custom)')
    p.add_argument('--rate-scale', type=float, default=1.0, help='Multiplicative scale to apply to computed rates')
    p.add_argument('--debug', action='store_true', help='Enable debug prints (minimal)')
    return p.parse_args()

def main():
    args = parse_args()
    with h5py.File(args.h5, 'r') as fh:
        parsed = {
            'energies': fh['/bands/energies_eV'][()],
            'velocities': fh['/bands/velocities_m_s'][()],
            'band_list': fh['/bands/band_list'][()].astype(int),
            'kcoords': fh['/kpath/coords'][()],
            'ik_list': fh['/kpath/ik_list'][()].astype(int),
            'qcoords': fh['/q/coords'][()] if '/q/coords' in fh else np.zeros((0,3)),
            'iq_list': fh['/q/iq_list'][()].astype(int) if '/q/iq_list' in fh else np.zeros((0,),dtype=int),
            'rows': fh['/scattering/rows'][()],
            'meta': {k: fh['/meta'].attrs[k] for k in fh['/meta'].attrs} if '/meta' in fh else {}
        }

    mc = MCDriftEPW(parsed, charge_type='electron', mode='B', scattering_choice=args.scat,
                    T=args.T, dt=args.dt, grid_xy=(int(args.grid[0]), int(args.grid[1])),
                    broadening_meV=args.broadening, rate_scale=args.rate_scale)

    out_dir = os.path.dirname(os.path.abspath(args.h5)) or '.'
    log_path = os.path.join(out_dir, args.log)

    hist = mc.run_particle_ensemble(args.Ex, nsteps=args.nsteps, npart=args.npart,
                                   conv_tol=args.conv_tol, out_interval=100,
                                   log_path=log_path, log_freq=args.log_freq,
                                   debug=args.debug)

    np.savez(args.out.replace('.npz','_drift.npz'), drift=np.array(hist))
    print(f"[ok] Completed. Results written: {args.out} and log: {log_path}")

if __name__ == '__main__':
    main()

