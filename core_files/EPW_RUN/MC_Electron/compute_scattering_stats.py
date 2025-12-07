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
"""
from __future__ import annotations
import argparse, os, sys, time, math
from collections import defaultdict
import h5py
import numpy as np

# ---------------- constants ----------------
hbar = 1.054571817e-34        # J s
eV_to_J = 1.602176634e-19
meV_to_J = eV_to_J * 1e-3
kB_eV_per_K = 8.617333262145e-5
kB_meV_per_K = kB_eV_per_K * 1000.0

# ---------------- utilities ----------------
def lorentzian_meV(delta_meV: float, gamma_meV: float) -> float:
    """Lorentzian (area-normalized) in 1/meV units: (gamma/pi) / (delta^2 + gamma^2)"""
    return (gamma_meV / math.pi) / (delta_meV*delta_meV + gamma_meV*gamma_meV)

def detect_lcb_from_reference(energies: np.ndarray, cbm_ref_eV: float):
    mins = np.nanmin(energies, axis=1)
    idx = int(np.argmin(np.abs(mins - cbm_ref_eV)))
    return idx + 1, mins

def compute_rates(h5path,
                  gamma_meV: float = 20.0,
                  T: float = 300.0,
                  kbt_multiplier: float = 20.0,
                  cbm_ref_eV: float | None = None,
                  nq: int | None = None,
                  energy_tol_meV: float | None = None,
                  mode_flag: str = 'both',
                  verbose: bool = True):

    print(">>> USING NEWEST compute_rates() <<<")

    if energy_tol_meV is None:
        energy_tol_meV = float(gamma_meV)

    # === Load file ===
    with h5py.File(h5path, 'r') as fh:
        rows = fh['/scattering/rows']

        # read band classification (for LCB)
        if cbm_ref_eV is not None:
            lcb, _ = detect_lcb_from_reference(fh['/bands/energies_eV'][...], cbm_ref_eV)
            if verbose:
                print(f"[info] CBM reference provided: {cbm_ref_eV:.6f} → LCB = {lcb}")
        else:
            try:
                conduction = fh['/meta/classification/conduction_bands'][...]
                lcb = int(conduction[0])
            except:
                lcb = 1
            if verbose:
                print(f"[info] Using EPW's CBM band index guess = {lcb}")

        # === EPW-referenced CBM ===
        epw_cbm_list = []
        for r in rows:
            if int(r['ibnd']) == lcb:
                epw_cbm_list.append(float(r['enk_eV']))

        if len(epw_cbm_list) == 0:
            # fallback (rare)
            all_enk = [float(r['enk_eV']) for r in rows]
            CBM_epw_eV = min(all_enk)
        else:
            CBM_epw_eV = min(epw_cbm_list)

        CBM_epw_meV = CBM_epw_eV * 1000.0

        if verbose:
            print(f"[info] EPW CBM (from rows, band {lcb}) = {CBM_epw_eV:.6f} eV")

        # === q-weight ===
        if nq is None or nq <= 0:
            try:
                max_iq = int(np.max(rows['iq_pos']))
                Nq = max_iq + 1
            except:
                Nq = 101
        else:
            Nq = nq

        w_q = 1.0 / float(Nq)
        if verbose:
            print(f"[info] q-weight = 1/{Nq} = {w_q}")

        # === constants ===
        kBT_meV = kB_meV_per_K * float(T)
        kBT_eV  = kBT_meV / 1000.0

        gamma_J = gamma_meV * meV_to_J
        two_pi_over_hbar = 2.0 * math.pi / hbar

        allow_intra = mode_flag in ('intra', 'both')
        allow_inter_up = mode_flag in ('inter_up', 'both')
        allow_inter_down = mode_flag in ('inter_down', 'both')

        # === accumulation ===
        used_row_list = []
        total_rows = rows.shape[0]
        used_rows = 0
        rejected_rows = 0

        rates_acc = {}
        Gamma_abs_acc = {}
        Gamma_em_acc = {}
        counts_acc = {}
        enk_sum = {}
        enkq_sum = {}
        imode_gamma_acc = {}

        max_g_used = 0.0
        max_g_unfiltered = 0.0

        # === MAIN LOOP ===
        for r in rows:
            ib = int(r['ibnd'])
            jb = int(r['jbnd'])
            g_meV = float(r['g_abs_meV'])

            if g_meV > max_g_unfiltered:
                max_g_unfiltered = g_meV

            # conduction filter
            if ib < lcb or jb < lcb:
                rejected_rows += 1
                continue

            is_intra = (ib == jb)
            is_inter_up   = (ib == lcb and jb > ib)
            is_inter_down = (ib > jb and jb == lcb)

            if is_intra and not allow_intra:
                rejected_rows += 1; continue
            if is_inter_up and not allow_inter_up:
                rejected_rows += 1; continue
            if is_inter_down and not allow_inter_down:
                rejected_rows += 1; continue

            ik = int(r['ik_pos'])
            iq = int(r['iq_pos'])
            im = int(r['imode'])

            enk_eV  = float(r['enk_eV'])
            enkq_eV = float(r['enkq_eV'])
            omega_meV = float(r['omega_meV'])

            # relative to EPW CBM
            enk_rel_meV = 1000.0 * (enk_eV - CBM_epw_eV)

            if abs(enk_rel_meV) > (kbt_multiplier * kBT_meV):
                rejected_rows += 1
                continue

            # Bose factor
            if omega_meV <= 0:
                n_q = 0.0
            else:
                x = omega_meV / (kBT_meV + 1e-300)
                n_q = 1/(math.exp(x) - 1.0) if x < 700 else 0.0

            # ΔE (EPW-referenced)
            omega_eV = omega_meV / 1000.0
            delta_abs_eV = enkq_eV - enk_eV - omega_eV
            delta_em_eV  = enkq_eV - enk_eV + omega_eV

            delta_abs_meV = delta_abs_eV * 1000.0
            delta_em_meV  = delta_em_eV * 1000.0

            if abs(delta_abs_meV) > energy_tol_meV and abs(delta_em_meV) > energy_tol_meV:
                rejected_rows += 1
                continue

            # --- Golden rule (SI units) ---
            g_J = g_meV * meV_to_J
            prefactor = two_pi_over_hbar * (g_J * g_J)

            dA_J = delta_abs_eV * eV_to_J
            dE_J = delta_em_eV  * eV_to_J

            denomA = dA_J*dA_J + gamma_J*gamma_J
            denomE = dE_J*dE_J + gamma_J*gamma_J

            L_abs = (gamma_J/math.pi)/denomA if denomA > 0 else 0.0
            L_em  = (gamma_J/math.pi)/denomE if denomE > 0 else 0.0

            Gamma_abs = prefactor * n_q       * L_abs
            Gamma_em  = prefactor * (n_q+1.0) * L_em
            Gamma     = Gamma_abs + Gamma_em

            # q-weight
            Gamma_abs_w = Gamma_abs * w_q
            Gamma_em_w  = Gamma_em  * w_q
            Gamma_w     = Gamma     * w_q

            if g_meV > max_g_used:
                max_g_used = g_meV

            # event type
            if abs(delta_abs_meV) <= energy_tol_meV and abs(delta_em_meV) <= energy_tol_meV:
                event_type = "both"
            elif abs(delta_abs_meV) <= energy_tol_meV:
                event_type = "abs"
            else:
                event_type = "em"

            if is_intra:
                intra_inter = "intra"
            elif is_inter_up:
                intra_inter = "inter_up"
            elif is_inter_down:
                intra_inter = "inter_down"
            else:
                intra_inter = "other"

            used_row_list.append({
                'ik': ik, 'iq': iq, 'imode': im, 'ib': ib, 'jb': jb,
                'enk_eV': enk_eV, 'enkq_eV': enkq_eV, 'omega_meV': omega_meV,
                'g_abs_meV': g_meV,
                'Gamma_w': Gamma_w,
                'Gamma_abs_w': Gamma_abs_w,
                'Gamma_em_w': Gamma_em_w,
                'delta_abs_meV': delta_abs_meV,
                'delta_em_meV': delta_em_meV,
                'k_frac': r['k_frac'], 'kq_frac': r['kq_frac'],
                'q_frac': r['q_frac'] if 'q_frac' in r.dtype.names else (r['kq_frac'] - r['k_frac']),
                'intra_inter': intra_inter,
                'event': event_type
            })

            # accumulate
            key = (ik, ib, jb, intra_inter)
            rates_acc[key]      = rates_acc.get(key, 0.0) + Gamma_w
            Gamma_abs_acc[key]  = Gamma_abs_acc.get(key, 0.0) + Gamma_abs_w
            Gamma_em_acc[key]   = Gamma_em_acc.get(key, 0.0) + Gamma_em_w
            counts_acc[key]     = counts_acc.get(key, 0) + 1
            enk_sum[key]        = enk_sum.get(key, 0.0) + enk_eV
            enkq_sum[key]       = enkq_sum.get(key, 0.0) + enkq_eV

            if key not in imode_gamma_acc:
                imode_gamma_acc[key] = {}
            imode_gamma_acc[key][im] = imode_gamma_acc[key].get(im, 0.0) + Gamma_w

            used_rows += 1

        # === Prepare channel-level results ===
        results = []
        for key in rates_acc:
            ik, ib, jb, typ = key
            Gtot = rates_acc[key]
            tau = 1.0/Gtot if Gtot > 0 else float('inf')
            avg_enk = enk_sum[key] / counts_acc[key]
            avg_enkq = enkq_sum[key] / counts_acc[key]
            Gabs = Gamma_abs_acc.get(key, 0.0)
            Gem  = Gamma_em_acc.get(key, 0.0)
            dom = max(imode_gamma_acc[key], key=lambda x: imode_gamma_acc[key][x])

            results.append((ik, avg_enk*1000.0, Gtot, tau,
                            Gabs, Gem, counts_acc[key],
                            avg_enk, avg_enkq, ib, jb,
                            typ, "both", dom))

        # === Detailed balance: corrected totals ===
        corrected_Gabs = {}
        corrected_Gem  = {}

        for key in rates_acc:
            Gtot = rates_acc[key]
            rows_k = [r for r in used_row_list if (r['ik'], r['ib'], r['jb'], r['intra_inter']) == key]
            ome = [r['omega_meV'] for r in rows_k if r['omega_meV'] > 0]

            if len(ome) == 0:
                corrected_Gabs[key] = Gtot/2
                corrected_Gem[key]  = Gtot/2
                continue

            om_avg_meV = float(np.mean(ome))
            om_avg_eV = om_avg_meV/1000.0
            x = om_avg_eV/(kBT_eV + 1e-300)

            if x >= 700:
                corrected_Gabs[key] = Gtot/2
                corrected_Gem[key]  = Gtot/2
                continue

            nbar = 1.0/(math.exp(x) - 1.0)
            if nbar <= 0:
                corrected_Gabs[key] = 0
                corrected_Gem[key]  = Gtot
            else:
                R = (nbar+1)/nbar
                Gabs = Gtot/(1+R)
                corrected_Gabs[key] = Gabs
                corrected_Gem[key]  = Gtot - Gabs

        # === nᵢ-aware redistribution ===
        per_key_idx = {}
        for idx, row in enumerate(used_row_list):
            key = (row['ik'], row['ib'], row['jb'], row['intra_inter'])
            per_key_idx.setdefault(key, []).append(idx)

        changed = 0
        large = 0

        for key, idx_list in per_key_idx.items():
            orig_G = np.array([used_row_list[i]['Gamma_w'] for i in idx_list])
            tgtGa = corrected_Gabs[key]
            tgtGe = corrected_Gem[key]

            # compute nᵢ per row
            n_list = []
            for i in idx_list:
                om = used_row_list[i]['omega_meV']
                if om <= 0:
                    n_list.append(0.0)
                else:
                    xx = om/(kBT_meV+1e-300)
                    n_list.append(1.0/(math.exp(xx)-1.0))
            n_list = np.array(n_list)

            denom = 2*n_list + 1
            denom = np.where(denom == 0, 1e-300, denom)

            Ga_raw = orig_G * (n_list/denom)
            Ge_raw = orig_G * ((n_list+1)/denom)

            sGa = Ga_raw.sum()
            sGe = Ge_raw.sum()

            if sGa > 0:
                Ga_new = Ga_raw * (tgtGa/sGa)
            else:
                Ga_new = np.full(len(idx_list), tgtGa/len(idx_list))

            if sGe > 0:
                Ge_new = Ge_raw * (tgtGe/sGe)
            else:
                Ge_new = np.full(len(idx_list), tgtGe/len(idx_list))

            for j, ii in enumerate(idx_list):
                old = used_row_list[ii]['Gamma_w']
                new = Ga_new[j] + Ge_new[j]
                used_row_list[ii]['Gamma_w'] = new
                used_row_list[ii]['Gamma_abs_w'] = Ga_new[j]
                used_row_list[ii]['Gamma_em_w']  = Ge_new[j]

                if abs(new - old) > 1e-15:
                    changed += 1
                    if abs((new-old)/(old+1e-300)) > 0.05:
                        large += 1

        if verbose:
            print(f"[info] nᵢ-aware detailed-balance enforced. Changed={changed}, large={large}")
            print(f"[info] Total rows={total_rows}, used={used_rows}, rejected={rejected_rows}")

        meta = {"lcb": lcb, "CBM_epw_eV": CBM_epw_eV, "Nq": Nq}
        bookkeeping = {
            "total_rows": total_rows,
            "used_rows": used_rows,
            "rejected_rows": rejected_rows,
            "max_g_unfiltered": max_g_unfiltered,
            "max_g_used": max_g_used,
            "used_row_list": used_row_list,
        }

        return results, meta, bookkeeping


# ---------------- writers ----------------
def write_scattering_events_hdf5(raw_rows, outname="scattering_events.h5"):
    """Write per-event HDF5 (plus attempt to attach full EPW velocities if available)."""
    N = len(raw_rows)
    if N == 0:
        print("[WARN] No events to write to HDF5.")
        return

    # arrays
    ik = np.zeros(N, dtype=np.int32); iq = np.zeros(N, dtype=np.int32)
    ib = np.zeros(N, dtype=np.int32); jb = np.zeros(N, dtype=np.int32); imode = np.zeros(N, dtype=np.int32)
    enk_meV = np.zeros(N, dtype=np.float64); enkq_meV = np.zeros(N, dtype=np.float64)
    omega_meV = np.zeros(N, dtype=np.float64); g_meV = np.zeros(N, dtype=np.float64)
    Gabs_w = np.zeros(N, dtype=np.float64); Gem_w = np.zeros(N, dtype=np.float64); Gw = np.zeros(N, dtype=np.float64)
    k_frac = np.zeros((N,3), dtype=np.float64); kq_frac = np.zeros((N,3), dtype=np.float64); q_frac = np.zeros((N,3), dtype=np.float64)
    intra_inter = np.empty(N, dtype='S12'); event = np.empty(N, dtype='S8')

    # try to load full EPW velocities/kpts (optional)
    vel_full = None; k_full = None
    try:
        with h5py.File("epw_quadrupole_parsed.h5", "r") as fepw:
            if "bands/velocities_m_s" in fepw:
                vel_full = fepw["bands/velocities_m_s"][...]
            if "kpath/kpts_frac" in fepw:
                k_full = fepw["kpath/kpts_frac"][...]
            elif "kpath/coords" in fepw:
                k_full = fepw["kpath/coords"][...]
    except Exception:
        pass

    v_cart = np.zeros((N,3), dtype=np.float64); v_mag = np.zeros(N, dtype=np.float64)

    # fill arrays from raw_rows; recompute Gamma_w = Gamma_abs_w + Gamma_em_w for safety
    for i, r in enumerate(raw_rows):
        ik[i] = int(r['ik']); iq[i] = int(r['iq'])
        ib[i] = int(r['ib']); jb[i] = int(r['jb']); imode[i] = int(r['imode'])
        enk_meV[i] = float(r['enk_eV']) * 1000.0
        enkq_meV[i] = float(r['enkq_eV']) * 1000.0
        omega_meV[i] = float(r['omega_meV'])
        g_meV[i] = float(r.get('g_meV', r.get('g_abs_meV', 0.0)))

        # Force read as floats and recompute Gw
        ga = float(r.get('Gamma_abs_w', 0.0))
        ge = float(r.get('Gamma_em_w', 0.0))
        Gabs_w[i] = ga
        Gem_w[i]  = ge
        Gw[i]     = ga + ge          # <-- always recompute here

        k_frac[i] = np.asarray(r['k_frac'], dtype=float)
        kq_frac[i] = np.asarray(r['kq_frac'], dtype=float)
        q_frac[i] = np.asarray(r['q_frac'], dtype=float)

        intra_inter[i] = str(r.get('intra_inter', 'other')).encode("ascii")
        event[i] = str(r.get('event', 'both')).encode("ascii")

        if vel_full is not None:
            try:
                ib0 = ib[i] - 1
                ik0 = ik[i]
                v_cart[i] = vel_full[ib0, ik0, :]
                v_mag[i] = np.linalg.norm(v_cart[i])
            except Exception:
                pass

    # quick debug summary so you can verify what's being written
    try:
        mean_ga = float(np.mean(Gabs_w))
        mean_ge = float(np.mean(Gem_w)) + 1e-300
        print(f"[DEBUG write_hdf5] events={N}, mean(Ga)={mean_ga:.6e}, mean(Ge)={mean_ge:.6e}, mean Ga/Ge={mean_ga/mean_ge:.6e}")
    except Exception:
        pass

    with h5py.File(outname, "w") as h5:
        g = h5.create_group("scattering_events")
        g.create_dataset("ik", data=ik)
        g.create_dataset("iq", data=iq)
        g.create_dataset("ib", data=ib)
        g.create_dataset("jb", data=jb)
        g.create_dataset("imode", data=imode)
        g.create_dataset("enk_meV", data=enk_meV)
        g.create_dataset("enkq_meV", data=enkq_meV)
        g.create_dataset("omega_meV", data=omega_meV)
        g.create_dataset("g_meV", data=g_meV)
        g.create_dataset("Gamma_abs_w", data=Gabs_w)
        g.create_dataset("Gamma_em_w", data=Gem_w)
        g.create_dataset("Gamma_w", data=Gw)
        g.create_dataset("k_frac", data=k_frac)
        g.create_dataset("kq_frac", data=kq_frac)
        g.create_dataset("q_frac", data=q_frac)
        g.create_dataset("intra_inter", data=intra_inter)
        g.create_dataset("event", data=event)
        g.create_dataset("v_group_cart", data=v_cart)
        g.create_dataset("v_group_mag", data=v_mag)
        if vel_full is not None:
            g.create_dataset("full_velocities_m_s", data=vel_full, compression="gzip")
        if k_full is not None:
            g.create_dataset("full_k_frac", data=k_full, compression="gzip")
    print(f"[OK] Wrote HDF5 events: {outname}")


def write_scattering_events_txt(raw_rows, outname="scattering_events.txt"):
    header = (
    "#Raw electron-phonon scattering events (one-to-one per EPW interaction).\n"
    "# ik iq ib jb im  enk_meV  enkq_meV  omega_meV  g_meV  Gamma_abs_w  Gamma_em_w  Gamma_w  "
    "kx ky kz  qx qy qz  kqx kqy kqz  process  event\n"
    )

    with open(outname, "w") as fo:
        fo.write(header)
        for r in raw_rows:
            ik = int(r['ik']); iq = int(r['iq']); ib = int(r['ib']); jb = int(r['jb']); im = int(r['imode'])
            enk_meV = float(r['enk_eV']) * 1000.0
            enkq_meV = float(r['enkq_eV']) * 1000.0
            omega_meV = float(r['omega_meV'])
            g_meV = float(r.get('g_meV', r.get('g_abs_meV', 0.0)))

            # Force float reads and recompute Gw
            Gabs = float(r.get('Gamma_abs_w', 0.0))
            Gem  = float(r.get('Gamma_em_w', 0.0))
            Gw   = Gabs + Gem   # <-- recompute to be safe

            kx, ky, kz = tuple(np.array(r['k_frac'], dtype=float))
            kqx, kqy, kqz = tuple(np.array(r['kq_frac'], dtype=float))
            qx, qy, qz = tuple(np.array(r['q_frac'], dtype=float))
            intra_inter = str(r.get('intra_inter', 'other'))
            evt = str(r.get('event', 'both'))
            fo.write(f"{ik:4d} {iq:4d} {ib:4d} {jb:4d} {im:4d} "
                     f"{enk_meV:10.4f} {enkq_meV:10.4f} {omega_meV:8.3f} {g_meV:10.4f} "
                     f"{Gabs:12.6e} {Gem:12.6e} {Gw:12.6e} "
                     f"{kx:10.6f} {ky:10.6f} {kz:10.6f} "
                     f"{qx:10.6f} {qy:10.6f} {qz:10.6f} "
                     f"{kqx:10.6f} {kqy:10.6f} {kqz:10.6f} "
                     f"{intra_inter:10s} {evt:6s}\n")
    print(f"[OK] Wrote text events: {outname}")


def write_scattering_rates_txt(results, outname="scattering_rates.txt"):
    header1 = "#Intraband/interband scattering results (conduction manifold only). Event = abs/em/both.\n"
    sepline = "#--------------------------------------------------------------------------\n"
    colhdr = "#ik_pos E_k_meV Gamma_s tau_s Gamma_abs_s Gamma_em_s n_rows avg_enk_eV avg_enkq_eV ibnd jbnd type event imode\n"
    with open(outname, "w") as fo:
        fo.write(header1); fo.write(sepline); fo.write(colhdr); fo.write(sepline)
        for (ik, Ek_meV, G, tau, Gabs, Gem, nrows, avg_enk, avg_enkq, ib, jb, typ, event, imode) in results:
            fo.write(f"{int(ik):7d} {Ek_meV:12.6f} {G:14.6e} {tau:12.6e} "
                     f"{Gabs:14.6e} {Gem:14.6e} {int(nrows):6d} "
                     f"{avg_enk:12.6f} {avg_enkq:12.6f} {int(ib):6d} {int(jb):6d} {str(typ):8s} {str(event):6s} {int(imode):6d}\n")
    print(f"[OK] Wrote aggregated rates: {outname}")

def write_report_scatterings(report_path, args, meta, bookkeeping):
    used_rows = bookkeeping.get('used_rows', 0)
    total_rows = bookkeeping.get('total_rows', 0)
    rejected = bookkeeping.get('rejected_rows', total_rows - used_rows)
    with open(report_path, "w") as f:
        f.write("=== Scattering Summary ===\n")
        f.write(f"Run time (UTC): {time.asctime(time.gmtime())}\n\n")
        f.write("[params]\n")
        f.write(f"  T = {args.T} K\n")
        f.write(f"  gamma = {args.gamma} meV\n")
        f.write(f"  energy_tol = {args.energy_tol} meV\n\n")
        f.write("[counts]\n")
        f.write(f"  raw rows: {total_rows}\n")
        f.write(f"  used rows: {used_rows}\n")
        f.write(f"  rejected rows: {rejected}\n\n")
        f.write("[max_g]\n")
        f.write(f"  max_g_unfiltered = {bookkeeping.get('max_g_unfiltered',0.0):.6e}\n")
        f.write(f"  max_g_used       = {bookkeeping.get('max_g_used',0.0):.6e}\n\n")
        f.write("[done]\n")
    print(f"[OK] Wrote report: {report_path}")

# ---------------- CLI ----------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument('--h5', default='epw_quadrupole_parsed.h5')
    p.add_argument('--gamma', type=float, default=10.0)
    p.add_argument('--T', type=float, default=300.0)
    p.add_argument('--kbt-mult', type=float, default=5.0)
    p.add_argument('--energy-tol', type=float, default=10.0)
    p.add_argument('--nq', type=int, default=None)
    p.add_argument('--cbm-ref', type=float, default=None)
    p.add_argument('--mode-flag', choices=['intra','inter_up','inter_down','both'], default='both')
    p.add_argument('--no-verbose', dest='verbose', action='store_false', default=True)
    args = p.parse_args()

    if not os.path.exists(args.h5):
        print("[error] HDF5 not found:", args.h5); sys.exit(1)

    t0 = time.time()
    results, meta, bookkeeping = compute_rates(args.h5,
                                              gamma_meV=args.gamma,
                                              T=args.T,
                                              kbt_multiplier=args.kbt_mult,
                                              cbm_ref_eV=args.cbm_ref,
                                              nq=args.nq,
                                              energy_tol_meV=args.energy_tol,
                                              mode_flag=args.mode_flag,
                                              verbose=args.verbose)

    # write outputs
    write_scattering_events_txt(bookkeeping['used_row_list'], outname="scattering_events.txt")
    write_scattering_events_hdf5(bookkeeping['used_row_list'], outname="scattering_events.h5")
    write_scattering_rates_txt(results, outname="scattering_rates.txt")
    write_report_scatterings("report_scatterings.txt", args, meta, bookkeeping)

    if args.verbose:
        print(f"[info] Done. Time elapsed: {time.time()-t0:.2f}s")
        print("[info] Wrote: scattering_events.(txt|h5), scattering_rates.txt, report_scatterings.txt")

if __name__ == "__main__":
    main()

