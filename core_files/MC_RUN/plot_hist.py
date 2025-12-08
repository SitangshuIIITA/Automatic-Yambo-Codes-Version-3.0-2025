#!/usr/bin/env python3
"""
plot_epw_analysis.py

Read EPW-parsed HDF5 and produce plots:
 - histogram of |g|
 - histogram (log) of per-row weighted Gamma
 - tau(k) vs E_k
 - top phonon modes bar plot
 - log-log g vs omega scatter (with fit)

Produces PNG files under the current working directory (or /mnt/data if run on the same machine).

Usage:
  python plot_epw_analysis.py --h5 /mnt/data/epw_quadrupole_parsed.h5
"""

import argparse
import h5py
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict

# Physical constants
hbar = 1.054571817e-34  # J s
eV_to_J = 1.602176634e-19
meV_to_J = eV_to_J * 1e-3
kB_eV_per_K = 8.617333262145e-5
kB_meV_per_K = kB_eV_per_K * 1000.0

def lorentzian_meV(delta_meV, gamma_meV):
    """Normalized Lorentzian in meV (HWHM = gamma_meV). Units: 1/meV"""
    return (gamma_meV / math.pi) / (delta_meV*delta_meV + gamma_meV*gamma_meV)

def load_lcb_rows(h5path, cbm_ref=None):
    with h5py.File(h5path, 'r') as fh:
        energies = fh['/bands/energies_eV'][...]
        mins = np.nanmin(energies, axis=1)
        # detect EPW band closest to cbm_ref if provided, else use metadata first conduction band
        if cbm_ref is not None:
            lcb_idx = int(np.argmin(np.abs(mins - cbm_ref)))  # 0-based
            lcb = lcb_idx + 1
        else:
            try:
                conduction = fh['/meta/classification/conduction_bands'][...]
                lcb = int(conduction[0])
                lcb_idx = lcb - 1
            except Exception:
                # fallback to band with largest minimum
                lcb_idx = int(np.argmax(mins))
                lcb = lcb_idx + 1
        rows = fh['/scattering/rows'][...]
    # filter rows for LCB->LCB
    mask = (rows['ibnd'].astype(int) == lcb) & (rows['jbnd'].astype(int) == lcb)
    rows_lcb = rows[mask]
    return energies, lcb, lcb_idx, rows_lcb

def compute_per_row_gamma(rows_lcb, gamma_meV=5.0, T=300.0, qmesh=(12,12,1)):
    """
    Compute per-row Gamma (weighted by w_q=1/Nq) and other arrays.
    Returns dictionary of arrays and constants.
    """
    # arrays
    re = rows_lcb['re_meV'].astype(float)
    im = rows_lcb['im_meV'].astype(float)
    gabs = rows_lcb['g_abs_meV'].astype(float)
    enk_eV = rows_lcb['enk_eV'].astype(float)
    enkq_eV = rows_lcb['enkq_eV'].astype(float)
    omega_meV = rows_lcb['omega_meV'].astype(float)
    ik_pos = rows_lcb['ik_pos'].astype(int)
    iq_pos = rows_lcb['iq_pos'].astype(int)
    imode = rows_lcb['imode'].astype(int)

    # constants
    two_pi_over_hbar = 2.0 * math.pi / hbar
    kBT_meV = kB_meV_per_K * T
    N_q = qmesh[0] * qmesh[1] * qmesh[2]
    w_q = 1.0 / float(N_q)

    # prefactor (convert g to Joules)
    g_J = gabs * meV_to_J
    prefactor = two_pi_over_hbar * (g_J * g_J)  # units: 1/(J s) * J^2 = J/s

    enk_meV = enk_eV * 1000.0
    enkq_meV = enkq_eV * 1000.0

    delta_abs = enkq_meV - enk_meV - omega_meV
    delta_em  = enkq_meV - enk_meV + omega_meV

    L_abs = lorentzian_meV(delta_abs, gamma_meV)
    L_em  = lorentzian_meV(delta_em,  gamma_meV)

    # convert Lorentzian 1/meV -> 1/J
    L_abs_SI = L_abs / meV_to_J
    L_em_SI  = L_em  / meV_to_J

    # Bose occupations
    n_q = np.zeros_like(omega_meV)
    pos_mask = omega_meV > 0
    x = np.zeros_like(omega_meV)
    x[pos_mask] = omega_meV[pos_mask] / kBT_meV
    # avoid overflow for enormous x
    n_q[pos_mask] = np.where(x[pos_mask] < 700.0, 1.0 / (np.exp(x[pos_mask]) - 1.0), 0.0)

    Gamma_abs = prefactor * (n_q + 1.0) * L_abs_SI
    Gamma_em  = prefactor * n_q        * L_em_SI
    Gamma_tot_unw = Gamma_abs + Gamma_em
    Gamma_tot_w = Gamma_tot_unw * w_q

    out = {
        're': re, 'im': im, 'gabs': gabs,
        'enk_eV': enk_eV, 'enkq_eV': enkq_eV, 'omega_meV': omega_meV,
        'ik_pos': ik_pos, 'iq_pos': iq_pos, 'imode': imode,
        'prefactor': prefactor,
        'Gamma_abs': Gamma_abs, 'Gamma_em': Gamma_em,
        'Gamma_tot_unw': Gamma_tot_unw, 'Gamma_tot_w': Gamma_tot_w,
        'Gamma_abs_w': Gamma_abs * w_q, 'Gamma_em_w': Gamma_em * w_q,
        'w_q': w_q, 'N_q': N_q
    }
    return out

def per_k_summary(energies, lcb_idx, per_row):
    """Aggregate per-row weighted Gamma into per-k totals and return DataFrame + results list"""
    iks = per_row['ik_pos']
    Gamma_w = per_row['Gamma_tot_w']
    Gamma_abs_w = per_row['Gamma_abs_w']
    Gamma_em_w = per_row['Gamma_em_w']
    enk = per_row['enk_eV']
    enkq = per_row['enkq_eV']

    rates_per_ik = {}
    counts_per_ik = {}
    Gabs_per_ik = {}
    Gem_per_ik = {}
    enk_sum = {}
    enkq_sum = {}

    for i in range(len(iks)):
        ik = int(iks[i])
        rates_per_ik[ik] = rates_per_ik.get(ik, 0.0) + float(Gamma_w[i])
        counts_per_ik[ik] = counts_per_ik.get(ik, 0) + 1
        Gabs_per_ik[ik] = Gabs_per_ik.get(ik, 0.0) + float(Gamma_abs_w[i])
        Gem_per_ik[ik]  = Gem_per_ik.get(ik, 0.0) + float(Gamma_em_w[i])
        enk_sum[ik] = enk_sum.get(ik, 0.0) + float(enk[i])
        enkq_sum[ik] = enkq_sum.get(ik, 0.0) + float(enkq[i])

    results = []
    for ik in sorted(rates_per_ik.keys()):
        Gamma = rates_per_ik[ik]
        tau = 1.0 / Gamma if Gamma > 0.0 else float('inf')
        E_k_meV = energies[lcb_idx, ik] * 1000.0
        avg_enk = enk_sum[ik] / counts_per_ik[ik]
        avg_enkq = enkq_sum[ik] / counts_per_ik[ik]
        results.append((ik, E_k_meV, Gamma, tau, Gabs_per_ik[ik], Gem_per_ik[ik], counts_per_ik[ik], avg_enk, avg_enkq))
    df = pd.DataFrame(results, columns=['ik_pos','E_k_meV','Gamma_s','tau_s','Gamma_abs_s','Gamma_em_s','n_rows','avg_enk_eV','avg_enkq_eV'])
    return df

def plot_and_save(per_row, df_summary, out_prefix=""):
    # hist of |g|
    plt.figure()
    plt.hist(per_row['gabs'], bins=100)
    plt.xlabel("|g| (meV)")
    plt.ylabel("Count")
    plt.title("Histogram of |g| (LCB->LCB)")
    plt.tight_layout()
    fn_g = out_prefix + "hist_g_lcb.png"
    plt.savefig(fn_g)
    plt.close()

    # histogram of per-row weighted Gamma (log10)
    gamma_w = per_row['Gamma_tot_w']
    safe = np.clip(gamma_w, 1e-50, None)
    plt.figure()
    plt.hist(np.log10(safe), bins=120)
    plt.xlabel("log10(Gamma_row_weighted [s^-1])")
    plt.ylabel("Count")
    plt.title("Histogram of per-row weighted Gamma (log10)")
    plt.tight_layout()
    fn_hist_gamma = out_prefix + "hist_gamma_row_weighted_log.png"
    plt.savefig(fn_hist_gamma)
    plt.close()

    # tau vs E_k (scatter)
    plt.figure()
    plt.plot(df_summary['E_k_meV'], df_summary['tau_s'], marker='o', linestyle='none')
    plt.xlabel("E_k (meV)")
    plt.ylabel("tau (s)")
    plt.title("Lifetime tau(k) for LCB (weighted)")
    plt.yscale('log')
    plt.tight_layout()
    fn_tau = out_prefix + "tau_vs_Ek.png"
    plt.savefig(fn_tau)
    plt.close()

    # top modes bar
    # compute sum per mode
    mode_sum = df.groupby('imode')['Gamma_w'].sum() if False else None  # placeholder removed
    # we'll compute from per_row arrays:
    imode = per_row['imode']
    mode_map = defaultdict(float)
    for i, mid in enumerate(imode):
        mode_map[int(mid)] += float(per_row['Gamma_tot_w'][i])
    mode_items = sorted(mode_map.items(), key=lambda x: -x[1])
    top_modes = mode_items[:10]
    modes = [m for m,_ in top_modes]
    gammas = [g for _,g in top_modes]
    plt.figure(figsize=(8,4))
    plt.bar(range(len(modes)), gammas)
    plt.xticks(range(len(modes)), modes)
    plt.xlabel("phonon mode index (imode)")
    plt.ylabel("Summed Gamma (s^-1)")
    plt.title("Top-10 phonon modes by summed weighted Gamma")
    plt.tight_layout()
    fn_modes = out_prefix + "top_modes_bar.png"
    plt.savefig(fn_modes)
    plt.close()

    # log-log g vs omega with fit (only omega>0)
    mask_pos = per_row['omega_meV'] > 0
    omega_pos = per_row['omega_meV'][mask_pos]
    g_pos = per_row['gabs'][mask_pos]
    log_omega = np.log(omega_pos.clip(min=1e-12))
    log_g = np.log(g_pos.clip(min=1e-12))
    if len(log_omega) > 2:
        p, c = np.polyfit(log_omega, log_g, 1)
    else:
        p, c = (np.nan, np.nan)
    plt.figure(figsize=(6,5))
    plt.scatter(log_omega, log_g, s=4)
    plt.xlabel("log(omega_meV)")
    plt.ylabel("log(g_abs_meV)")
    plt.title(f"log-log g vs omega (slope={p:.3f})")
    plt.tight_layout()
    fn_loglog = out_prefix + "loglog_g_vs_omega.png"
    plt.savefig(fn_loglog)
    plt.close()

    return fn_g, fn_hist_gamma, fn_tau, fn_modes, fn_loglog

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('--h5', required=True, help="Path to epw_quadrupole_parsed.h5")
    p.add_argument('--cbm-ref', type=float, default=-1.108100, help="NSCF CBM reference (eV)")
    p.add_argument('--gamma', type=float, default=5.0, help="Lorentzian HWHM (meV)")
    p.add_argument('--qmesh', type=str, default='12x12x1', help="qmesh as NxMxL (default 12x12x1)")
    p.add_argument('--out-prefix', type=str, default='/mnt/data/', help="output folder/prefix")
    args = p.parse_args()

    qmesh = tuple(int(x) for x in args.qmesh.split('x'))
    energies, lcb, lcb_idx, rows_lcb = load_lcb_rows(args.h5, cbm_ref=args.cbm_ref)
    print(f"Detected EPW band for CBM ref {args.cbm_ref}: EPW band {lcb}")

    per_row = compute_per_row_gamma(rows_lcb, gamma_meV=args.gamma, T=300.0, qmesh=qmesh)

    # attach Gamma_w inside per_row DataFrame for convenience
    df_rows = pd.DataFrame({
        'ik_pos': per_row['ik_pos'],
        'iq_pos': per_row['iq_pos'],
        'imode': per_row['imode'],
        'gabs': per_row['gabs'],
        'omega_meV': per_row['omega_meV'],
        'enk_eV': per_row['enk_eV'],
        'enkq_eV': per_row['enkq_eV'],
        'Gamma_w': per_row['Gamma_tot_w'],
        'Gamma_unw': per_row['Gamma_tot_unw']
    })

    # per-k summary
    df_summary = per_k_summary(energies, lcb_idx, per_row)
    csv_out = args.out_prefix.rstrip('/') + "/lcb_rates_summary.csv"
    df_summary.to_csv(csv_out, index=False)
    print("Wrote per-k summary to:", csv_out)

    # save top-20 rows table
    df_top = df_rows.sort_values('Gamma_w', ascending=False).head(20)
    top_out = args.out_prefix.rstrip('/') + "/top20_weighted_gamma_rows.csv"
    df_top.to_csv(top_out, index=False)
    print("Wrote top-20 weighted rows to:", top_out)

    # create and save plots
    fn_g, fn_hist_gamma, fn_tau, fn_modes, fn_loglog = plot_and_save(per_row, df_summary, out_prefix=args.out_prefix)
    print("Saved plots:")
    print(" -", fn_g)
    print(" -", fn_hist_gamma)
    print(" -", fn_tau)
    print(" -", fn_modes)
    print(" -", fn_loglog)

    print("Done.")

