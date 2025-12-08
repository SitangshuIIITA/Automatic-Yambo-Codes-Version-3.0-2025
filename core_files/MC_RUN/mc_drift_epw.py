#!/usr/bin/env python3
"""
mc_drift_epw.py

Monte-Carlo drift driver for EPW SERTA-style scattering, using EPW-provided |g|
directly (no extra division by ħω and no ad-hoc scaling).

Features:
 - Uses g_abs_meV from EPW directly (converted to J).
 - No area normalization, no 1/|∇E| DOS factor.
 - Energies referenced to CBM (electrons) or VBM (holes) when computing rates.
 - k+q -> nearest-k mapping precomputed.
 - Vectorized grouped scattering table for speed.
 - Debug mode prints top contributing rows.
"""
from __future__ import annotations
import argparse
import math
import os
import sys
from collections import defaultdict
import numpy as np
import h5py

# Physical constants
eV_to_J = 1.602176634e-19
hbar = 1.054571817e-34
q_e = 1.602176634e-19
kb = 1.380649e-23
pi = math.pi

# -------------------------
# Helpers
# -------------------------
def bose_factor(omega_J: float, T: float) -> float:
    if omega_J <= 0.0 or T <= 0.0:
        return 0.0
    x = omega_J / (kb * T)
    if x > 700.0:
        return 0.0
    return 1.0 / (math.exp(x) - 1.0)

def frac_wrap(frac: np.ndarray) -> np.ndarray:
    return frac - np.floor(frac)

# -------------------------
# Scattering using EPW-provided |g| directly
# -------------------------
def calc_row_rate_using_g(row, T: float, broadening_meV: float,
                          cbm_e: float, vbm_e: float, carrier_is_electron: bool,
                          rate_scale: float = 1.0, max_rate: float = 1e15) -> float:
    """
    Compute single-row scattering rate using EPW-provided |g| (g_abs_meV).
    No division by ħ*omega — exactly using |g|^2.
    """
    try:
        g_mev = float(row['g_abs_meV'])
    except Exception:
        return 0.0
    if g_mev == 0.0 or not np.isfinite(g_mev):
        return 0.0

    # convert to Joules
    g_J = (g_mev * 1e-3) * eV_to_J   # meV -> eV -> J
    g2 = g_J * g_J

    omega_meV = float(row['omega_meV'])
    omega_J = (omega_meV * 1e-3) * eV_to_J

    # reference energies to CBM / VBM (eV)
    enk = float(row['enk_eV'])
    enkq = float(row['enkq_eV'])
    if carrier_is_electron:
        Ei_rel = enk - cbm_e
        Ef_rel = enkq - cbm_e
    else:
        Ei_rel = vbm_e - enk
        Ef_rel = vbm_e - enkq

    if not np.isfinite(Ei_rel) or Ei_rel < 0.0:
        Ei_rel = 0.0
    if not np.isfinite(Ef_rel) or Ef_rel < 0.0:
        Ef_rel = 0.0

    Ei_J = Ei_rel * eV_to_J
    Ef_J = Ef_rel * eV_to_J

    # Lorentzian deltas
    Delta_em = Ef_J - Ei_J + omega_J
    Delta_ab = Ef_J - Ei_J - omega_J
    gamma_J = (broadening_meV * 1e-3) * eV_to_J

    Nq = bose_factor(omega_J, T)

    pref = (2.0 * pi / hbar) * g2
    L_em = (gamma_J / pi) / (Delta_em * Delta_em + gamma_J * gamma_J)
    L_ab = (gamma_J / pi) / (Delta_ab * Delta_ab + gamma_J * gamma_J)

    rate = pref * ((Nq + 1.0) * L_em + Nq * L_ab)
    rate *= float(rate_scale)

    if not np.isfinite(rate) or rate <= 0.0:
        return 0.0
    return float(min(rate, max_rate))

# -------------------------
# Vectorized grouped-scattering table builder (uses g^2 directly)
# -------------------------
def build_scattering_table_fast(rows_array, T: float, broadening_meV: float,
                                rate_scale: float = 1.0, max_rate: float = 1e15,
                                cbm_e: float = 0.0, vbm_e: float = 0.0, carrier_is_electron: bool = True):
    """
    Compute per-row rates referencing energies to CBM/VBM, using EPW-provided |g| directly.
    Group rows by (ibnd, ik_pos) and return grouped arrays for fast selection.
    """
    Nrows = rows_array.shape[0]
    if Nrows == 0:
        return {
            'group_keys': np.zeros((0,2), dtype=np.int32),
            'group_ptr': np.array([0], dtype=np.int32),
            'rows_idx_sorted': np.array([], dtype=np.int32),
            'rates_sorted': np.array([], dtype=np.float64),
            'rates_cumsum': np.array([], dtype=np.float64),
            'rows_array_sorted_view': rows_array
        }

    ib = rows_array['ibnd'].astype(np.int32)
    ik = rows_array['ik_pos'].astype(np.int32)
    omega_meV = rows_array['omega_meV'].astype(np.float64)
    g_abs_meV = rows_array['g_abs_meV'].astype(np.float64)
    enk_eV = rows_array['enk_eV'].astype(np.float64)
    enkq_eV = rows_array['enkq_eV'].astype(np.float64)

    # Reference energies (eV)
    if carrier_is_electron:
        Ei_rel = enk_eV - cbm_e
        Ef_rel = enkq_eV - cbm_e
    else:
        Ei_rel = vbm_e - enk_eV
        Ef_rel = vbm_e - enkq_eV

    Ei_rel = np.where(np.isfinite(Ei_rel) & (Ei_rel > 0.0), Ei_rel, 0.0)
    Ef_rel = np.where(np.isfinite(Ef_rel) & (Ef_rel > 0.0), Ef_rel, 0.0)

    Ei_J = Ei_rel * eV_to_J
    Ef_J = Ef_rel * eV_to_J
    g_J = (g_abs_meV * 1e-3) * eV_to_J
    g2 = g_J * g_J
    omega_J = (omega_meV * 1e-3) * eV_to_J
    gamma_J = (broadening_meV * 1e-3) * eV_to_J

    # Nq vector
    Nq = np.zeros_like(omega_J)
    for i, o in enumerate(omega_J):
        Nq[i] = bose_factor(o, T)

    Delta_em = Ef_J - Ei_J + omega_J
    Delta_ab = Ef_J - Ei_J - omega_J

    pref = (2.0 * pi / hbar) * g2
    L_em = (gamma_J / pi) / (Delta_em * Delta_em + gamma_J * gamma_J)
    L_ab = (gamma_J / pi) / (Delta_ab * Delta_ab + gamma_J * gamma_J)

    rates = pref * ((Nq + 1.0) * L_em + Nq * L_ab)
    rates = np.where(np.isfinite(rates) & (rates > 0.0), np.minimum(rates * rate_scale, max_rate), 0.0)

    # sort & group
    order = np.lexsort((ik, ib))
    ibs = ib[order]
    iks = ik[order]
    rates_s = rates[order]
    rows_idx_sorted = order.astype(np.int32)

    if rates_s.size == 0:
        return {
            'group_keys': np.zeros((0,2), dtype=np.int32),
            'group_ptr': np.array([0], dtype=np.int32),
            'rows_idx_sorted': np.array([], dtype=np.int32),
            'rates_sorted': np.array([], dtype=np.float64),
            'rates_cumsum': np.array([], dtype=np.float64),
            'rows_array_sorted_view': rows_array[rows_idx_sorted]
        }

    rates_cumsum = np.cumsum(rates_s)

    change = np.empty(len(order), dtype=bool)
    change[0] = True
    change[1:] = (ibs[1:] != ibs[:-1]) | (iks[1:] != iks[:-1])
    starts = np.nonzero(change)[0]
    ends = np.append(starts[1:], len(order))
    Ng = len(starts)

    group_keys = np.empty((Ng,2), dtype=np.int32)
    group_ptr = np.empty(Ng + 1, dtype=np.int32)
    for i, (s, e) in enumerate(zip(starts, ends)):
        group_keys[i,0] = int(ibs[s]); group_keys[i,1] = int(iks[s])
        group_ptr[i] = int(s)
    group_ptr[Ng] = len(order)

    rows_array_sorted_view = rows_array[rows_idx_sorted]

    return {
        'group_keys': group_keys,
        'group_ptr': group_ptr,
        'rows_idx_sorted': rows_idx_sorted,
        'rates_sorted': rates_s,
        'rates_cumsum': rates_cumsum,
        'rows_array_sorted_view': rows_array_sorted_view
    }

# -------------------------
# Interp / masses
# -------------------------
def interp_state_quantities(energies, velocities, kcoords, kpos_float):
    nb, Nk = energies.shape
    if kpos_float <= 0:
        il = 0; w = 0.0
    elif kpos_float >= Nk - 1:
        il = Nk - 2; w = 1.0
    else:
        il = int(math.floor(kpos_float)); w = kpos_float - il
    ir = il + 1
    E_il = energies[:, il]; E_ir = energies[:, ir]
    V_il = velocities[:, il, :]; V_ir = velocities[:, ir, :]
    E_interp = np.where(np.isnan(E_il), E_ir, np.where(np.isnan(E_ir), E_il, (1.0 - w) * E_il + w * E_ir))
    V_interp = np.where(np.isnan(V_il), V_ir, np.where(np.isnan(V_ir), V_il, (1.0 - w) * V_il + w * V_ir))
    return E_interp, V_interp

# -------------------------
# MC Engine
# -------------------------
class MCDriftEPW:
    def __init__(self, parsed, charge_type='electron', mode='B', scattering_choice='both',
                 T=300.0, dt=1e-15, broadening_meV=20.0, rate_scale=1.0):
        self.parsed = parsed
        self.energies = parsed['energies']
        self.velocities = parsed['velocities']
        self.band_list = parsed['band_list']
        self.kcoords = parsed['kcoords']
        self.qcoords = parsed.get('qcoords', np.zeros((0,3)))
        self.rows = parsed['rows']
        self.ik_list = parsed.get('ik_list', np.arange(self.kcoords.shape[0]))
        self.iq_list = parsed.get('iq_list', np.arange(self.qcoords.shape[0]))

        self.charge = -q_e if str(charge_type).lower().startswith('e') else +q_e
        self.charge_type = charge_type
        self.mode = mode
        self.scattering_choice = scattering_choice
        self.T = float(T)
        self.dt = float(dt)
        self.nb, self.Nk = self.energies.shape

        self.broadening_meV = float(broadening_meV)
        self.rate_scale = float(rate_scale)

        # meta
        meta = parsed.get('meta', {}) or {}
        self.fermi = float(meta.get('E_fermi_eV')) if ('E_fermi_eV' in meta and meta.get('E_fermi_eV') is not None) else np.nan

        # robust CBM/VBM from arrays
        band_E_min = np.nanmin(self.energies, axis=1)
        band_E_max = np.nanmax(self.energies, axis=1)
        tol = 1e-6
        if not np.isnan(self.fermi):
            cb_candidates = [i for i in range(len(self.band_list)) if band_E_min[i] >= (self.fermi - tol)]
            if cb_candidates:
                ib_cb = cb_candidates[np.argmin([band_E_min[i] for i in cb_candidates])]
            else:
                ib_cb = int(np.nanargmin(band_E_min))
        else:
            ib_cb = int(np.nanargmin(band_E_min))
        if not np.isnan(self.fermi):
            vb_candidates = [i for i in range(len(self.band_list)) if band_E_max[i] <= (self.fermi + tol)]
            if vb_candidates:
                ib_vb = vb_candidates[np.argmax([band_E_max[i] for i in vb_candidates])]
            else:
                ib_vb = int(np.nanargmax(band_E_max))
        else:
            ib_vb = int(np.nanargmax(band_E_max))

        self.cbm_band_idx = int(ib_cb); self.cbm_band = int(self.band_list[self.cbm_band_idx]); self.cbm_e = float(band_E_min[self.cbm_band_idx])
        self.vbm_band_idx = int(ib_vb); self.vbm_band = int(self.band_list[self.vbm_band_idx]); self.vbm_e = float(band_E_max[self.vbm_band_idx])

        # build grouped scattering table
        built = build_scattering_table_fast(self.rows, T=self.T, broadening_meV=self.broadening_meV,
                                           rate_scale=self.rate_scale, max_rate=1e15,
                                           cbm_e=self.cbm_e, vbm_e=self.vbm_e, carrier_is_electron=(self.charge < 0))
        self.group_keys = built['group_keys']
        self.group_ptr = built['group_ptr']
        self.rows_idx_sorted = built['rows_idx_sorted']
        self.rates_sorted = built['rates_sorted']
        self.rates_cumsum = built['rates_cumsum']
        self.rows_array_sorted_view = built['rows_array_sorted_view']

        # scatter_table simple mapping and totals
        self.scatter_table = {}
        Ng = self.group_keys.shape[0]
        for gi in range(Ng):
            ib = int(self.group_keys[gi,0]); ik = int(self.group_keys[gi,1])
            s = int(self.group_ptr[gi]); e = int(self.group_ptr[gi+1])
            self.scatter_table[(ib, ik)] = (s, e)
        self.total_rate = {}
        for gi in range(Ng):
            ib = int(self.group_keys[gi,0]); ik = int(self.group_keys[gi,1])
            s = int(self.group_ptr[gi]); e = int(self.group_ptr[gi+1])
            tot = float(self.rates_cumsum[e - 1] - (self.rates_cumsum[s - 1] if s > 0 else 0.0))
            self.total_rate[(ib, ik)] = tot

        # precompute k+q -> nearest-k mapping
        self.Nq = int(self.qcoords.shape[0]) if self.qcoords is not None else 0
        if self.Nq > 0:
            Nk = self.Nk
            kcoords = self.kcoords
            qcoords = self.qcoords
            kplusq_map = np.empty((Nk, self.Nq), dtype=np.int32)
            for ik in range(Nk):
                t = frac_wrap(kcoords[ik] + qcoords)  # (Nq,3)
                dif = np.abs(kcoords[None, :, :] - t[:, None, :])  # (Nq, Nk,3)
                dif = np.minimum(dif, 1.0 - dif)
                dist2 = np.sum(dif * dif, axis=2)  # (Nq, Nk)
                nearest = np.argmin(dist2, axis=1)  # (Nq,)
                kplusq_map[ik, :] = nearest
            self.kplusq_map = kplusq_map
        else:
            self.kplusq_map = None

    def sample_initial_state(self, band='auto'):
        if band == 'auto':
            if self.charge < 0:
                ib_idx = self.cbm_band_idx
            else:
                ib_idx = self.vbm_band_idx
        else:
            try:
                ib_idx = int(np.where(self.band_list == int(band))[0][0])
            except Exception:
                ib_idx = 0
        if self.charge < 0:
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
        dk_cart = (self.charge * Ex / hbar) * self.dt * np.array([1.0, 0.0, 0.0])
        recip = self.parsed.get('recip', None)
        if recip is not None:
            kcoords_cart = np.dot(self.kcoords, recip)
        else:
            kcoords_cart = self.kcoords.copy()
        Nk = self.Nk
        if kpos_float <= 0:
            il = 0; w = 0.0
        elif kpos_float >= Nk - 1:
            il = Nk - 2; w = 1.0
        else:
            il = int(math.floor(kpos_float)); w = kpos_float - il
        ir = il + 1
        kcart_il = kcoords_cart[il]; kcart_ir = kcoords_cart[ir]
        kcart_cur = (1.0 - w) * kcart_il + w * kcart_ir
        kcart_new = kcart_cur + dk_cart
        best_i = 0; best_t = 0.0; best_dist = 1e99
        for i in range(Nk - 1):
            a = kcoords_cart[i]; b = kcoords_cart[i + 1]; ab = b - a
            denom = np.dot(ab, ab)
            if denom == 0:
                t = 0.0
            else:
                t = np.dot(kcart_new - a, ab) / denom
                t = max(0.0, min(1.0, t))
            proj = a + t * ab
            d2 = np.sum((proj - kcart_new) ** 2)
            if d2 < best_dist:
                best_dist = d2; best_i = i; best_t = t
        new_kpos_float = best_i + best_t
        return max(0.0, min(self.Nk - 1.0, new_kpos_float))

    def map_kplusq(self, ik_pos, iq_pos):
        if self.kplusq_map is None or iq_pos < 0 or iq_pos >= self.Nq:
            return int(ik_pos)
        return int(self.kplusq_map[int(ik_pos), int(iq_pos)])

    def select_scattering_event(self, ibnd, ik_pos):
        if self.group_keys.size == 0:
            return None
        mask = (self.group_keys[:, 0] == int(ibnd)) & (self.group_keys[:, 1] == int(ik_pos))
        idx = np.nonzero(mask)[0]
        if idx.size == 0:
            return None
        gi = int(idx[0])
        s = int(self.group_ptr[gi]); e = int(self.group_ptr[gi + 1])
        if e <= s:
            return None
        cumsum = self.rates_cumsum
        base = float(cumsum[s - 1]) if s > 0 else 0.0
        total = float(cumsum[e - 1] - base)
        if total <= 0.0:
            return None
        xi = np.random.random() * total
        rel_idx = int(np.searchsorted(cumsum[s:e] - base, xi))
        if rel_idx >= (e - s):
            rel_idx = (e - s) - 1
        row_sorted_index = int(self.rows_idx_sorted[s + rel_idx])
        return row_sorted_index

    def run_particle_ensemble(self, Ex, nsteps=10000, npart=1000, conv_tol=1e-6,
                              out_interval=100, log_path=None, log_freq=50, debug=False):
        particles = []
        for _ in range(npart):
            b, kp = self.sample_initial_state('auto')
            particles.append({'band': b, 'kpos': float(kp)})
        drift_hist = []
        t = 0.0
        v_d = 0.0

        logf = None
        if log_path is not None:
            out_dir = os.path.dirname(os.path.abspath(log_path)) or '.'
            os.makedirs(out_dir, exist_ok=True)
            logf = open(log_path, 'w')
            logf.write("===== MC DRIFT SIMULATION LOG =====\n")
            logf.write(f"Mode: {self.mode}\n")
            logf.write(f"Carrier: {'electron' if self.charge < 0 else 'hole'}\n")
            logf.write(f"Electric field Ex = {Ex:.6e} V/m\n")
            logf.write(f"Temperature = {self.T} K\n")
            logf.write(f"Time step dt = {self.dt} s\n")
            logf.write(f"Particles = {npart}\n")
            logf.write(f"Steps = {nsteps}\n")
            logf.write(f"broadening_meV = {self.broadening_meV}\n")
            logf.write(f"rate_scale = {self.rate_scale}\n")
            logf.write('-----------------------------------\n\n')
            logf.flush()

        total_scatter = 0
        emission_count = 0
        absorption_count = 0

        # debug diagnostics: top sample rates
        if debug and logf is not None:
            sample_rates = []
            for i in range(min(2000, len(self.rows))):
                try:
                    r = self.rows[i]
                    rr = calc_row_rate_using_g(r, self.T, self.broadening_meV,
                                               self.cbm_e, self.vbm_e, (self.charge < 0),
                                               self.rate_scale)
                    sample_rates.append((i, rr, float(r['g_abs_meV']), float(r['omega_meV']),
                                         float(r['enk_eV']), float(r['enkq_eV'])))
                except Exception:
                    continue
            sample_rates.sort(key=lambda x: x[1], reverse=True)
            logf.write("=== DEBUG: top sample rates (idx, rate, g_meV, omega_meV, enk, enkq) ===\n")
            for t in sample_rates[:80]:
                logf.write(f"{t[0]:6d}  {t[1]:12.6e}  {t[2]:8.3f}  {t[3]:8.3f}  {t[4]:11.6f}  {t[5]:11.6f}\n")
            logf.write("\n")
            logf.flush()

        for step in range(nsteps):
            vx_acc = 0.0
            energies_acc = 0.0
            scat_rates_acc = 0.0

            for p in particles:
                kpos_new = self.accelerate_kpos(p['kpos'], Ex)
                p['kpos'] = kpos_new

                E_p, V_p = interp_state_quantities(self.energies, self.velocities, self.kcoords, p['kpos'])
                try:
                    ib_i = int(np.where(self.band_list == p['band'])[0][0])
                except Exception:
                    ib_i = 0
                vvec = V_p[ib_i]
                vx_acc += vvec[0]

                E_raw = E_p[ib_i]
                if self.charge < 0:
                    Ekin = E_raw - self.cbm_e
                else:
                    Ekin = self.vbm_e - E_raw
                if not np.isfinite(Ekin) or Ekin < 0.0:
                    Ekin = 0.0
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
                        new_ikpos = self.map_kplusq(ik_nearest, int(rrow['iq_pos']))
                        try:
                            dE = float(rrow['enkq_eV']) - float(rrow['enk_eV'])
                        except Exception:
                            dE = 0.0
                        if dE < 0:
                            emission_count += 1
                        elif dE > 0:
                            absorption_count += 1
                        total_scatter += 1
                        p['band'] = new_band
                        p['kpos'] = float(new_ikpos)

            v_d_instant = vx_acc / float(npart) if npart > 0 else 0.0
            avg_energy = (energies_acc / float(npart)) if npart > 0 else 0.0
            avg_scat_rate = (scat_rates_acc / float(npart)) if npart > 0 else 0.0
            drift_hist.append((t, v_d_instant))
            t += self.dt

            if logf is not None and log_freq > 0 and (step % log_freq == 0):
                converged_flag = False
                if len(drift_hist) >= log_freq + 1:
                    recent = np.mean([v for (_, v) in drift_hist[-log_freq:]])
                    if abs(recent - v_d) < conv_tol:
                        converged_flag = True
                logf.write(f"[step {step}]\n")
                logf.write(f"   avg_vx      = {v_d_instant:.5e} m/s\n")
                logf.write(f"   avg_energy  = {avg_energy * 1e3:.5e} meV\n")
                logf.write(f"   scat_rate   = {avg_scat_rate:.5e} 1/s\n")
                logf.write(f"   emission    = {emission_count}\n")
                logf.write(f"   absorption  = {absorption_count}\n")
                logf.write(f"   converged   = {converged_flag}\n")
                logf.write('-----------------------------------\n')
                if step % (5 * log_freq) == 0:
                    logf.flush()

            if step % out_interval == 0 and step > 0:
                recent = np.mean([v for (_, v) in drift_hist[-out_interval:]])
                if abs(recent - v_d) < conv_tol:
                    if logf is not None:
                        logf.write(f"[convergence] drift velocity stable within {conv_tol} at step {step}\n\n")
                        logf.flush()
                    break
                v_d = recent

        if logf is not None:
            logf.write('===== FINAL RESULTS =====\n')
            last_v = drift_hist[-1][1] if len(drift_hist) > 0 else 0.0
            logf.write(f"Steady-state drift velocity = {last_v:.5e} m/s\n")
            mu = 0.0
            if Ex != 0:
                mu = (last_v / Ex) * 1e4
            logf.write(f"Mobility = {mu:.5e} cm^2/V/s\n")
            logf.write(f"Total scattering events = {total_scatter}\n")
            logf.write(f"Total emission = {emission_count}\n")
            logf.write(f"Total absorption = {absorption_count}\n")
            logf.write('Simulation completed.\n')
            logf.write('====================================\n')
            logf.flush()
            logf.close()

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
    p.add_argument('--broadening', type=float, default=20.0, help='Lorentzian HWHM in meV')
    p.add_argument('--rate-scale', type=float, default=1.0, help='Multiplicative scale applied to computed rates (default 1)')
    p.add_argument('--debug', action='store_true', help='Enable debug diagnostics')
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
            'iq_list': fh['/q/iq_list'][()].astype(int) if '/q/iq_list' in fh else np.zeros((0,), dtype=int),
            'rows': fh['/scattering/rows'][()],
            'meta': {k: fh['/meta'].attrs[k] for k in fh['/meta'].attrs} if '/meta' in fh else {},
            'recip': fh['/reciprocal_lattice/b_vectors'][()] if '/reciprocal_lattice/b_vectors' in fh else None
        }

    mc = MCDriftEPW(parsed, charge_type='electron', mode='B', scattering_choice=args.scat,
                    T=args.T, dt=args.dt, broadening_meV=args.broadening, rate_scale=args.rate_scale)

    out_dir = os.path.dirname(os.path.abspath(args.h5)) or '.'
    log_path = os.path.join(out_dir, args.log)

    hist = mc.run_particle_ensemble(args.Ex, nsteps=args.nsteps, npart=args.npart,
                                   conv_tol=args.conv_tol, out_interval=100,
                                   log_path=log_path, log_freq=args.log_freq,
                                   debug=args.debug)

    np.savez(args.out.replace('.npz', '_drift.npz'), drift=np.array(hist))
    print(f"[ok] Completed. Results written: {args.out} and log: {log_path}")

if __name__ == '__main__':
    main()

