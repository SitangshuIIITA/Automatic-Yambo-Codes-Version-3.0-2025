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
    p.add_argument('--gamma', type=float, default=20.0)
    p.add_argument('--T', type=float, default=300.0)
    p.add_argument('--kbt-mult', type=float, default=300.0)
    p.add_argument('--energy-tol', type=float, default=300.0)
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

