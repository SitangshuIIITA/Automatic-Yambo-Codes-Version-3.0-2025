#!/usr/bin/env python3
"""
mc_rate_intraband.py

Intraband/interband scattering rates for the lowest conduction band (LCB),
with strict energy-conservation filtering for interband transitions.

Adds classification of scattering event:
    event = "abs"   if Gamma_abs > Gamma_em
    event = "em"    if Gamma_em  > Gamma_abs
    event = "both"  if |Gamma_abs - Gamma_em| < 0.05*(Gamma_abs + Gamma_em)

Now an additional column "event" appears after "type".
Also prints dominant imode (phonon branch contributing the largest Gamma_w).
"""

import argparse
import h5py
import numpy as np
import math
import sys

# Physical constants
hbar = 1.054571817e-34  # J s
eV_to_J = 1.602176634e-19
meV_to_J = eV_to_J * 1e-3
kB_eV_per_K = 8.617333262145e-5
kB_meV_per_K = kB_eV_per_K * 1000.0


def lorentzian_meV(delta_meV, gamma_meV):
    return (gamma_meV / math.pi) / (delta_meV * delta_meV + gamma_meV * gamma_meV)


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
                    print(f"[info] Using metadata LCB = band {lcb}")
            except Exception:
                raise RuntimeError("No conduction_bands metadata and no --cbm-ref provided. Please supply --cbm-ref")

        energy_idx = lcb - 1
        Eband = energies[energy_idx, :]
        E_CBM_eV = float(np.nanmin(Eband))
        E_CBM_meV = 1000.0 * E_CBM_eV

        if verbose:
            print(f"[info] Using EPW band {lcb} as LCB")
            print(f"[info] CBM (band {lcb}) = {E_CBM_eV:.6f} eV ({E_CBM_meV:.3f} meV)")

        rows = fh['/scattering/rows']

        # q-weight
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

        # NEW: store Gamma_w per imode
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

        for r in rows:
            g_abs_meV = float(r['g_abs_meV'])

            # Track unfiltered maxima
            if g_abs_meV > max_g_unfiltered:
                max_g_unfiltered = g_abs_meV
                max_g_unfiltered_info = (
                    int(r['ik_pos']), int(r['iq_pos']), int(r['imode']),
                    float(r['enk_eV']), float(r['enkq_eV']), float(r['omega_meV']),
                    int(r['ibnd']), int(r['jbnd'])
                )

            ib = int(r['ibnd'])
            jb = int(r['jbnd'])

            # Only conduction manifold
            if ib < lcb or jb < lcb:
                rejected_rows += 1
                continue

            # identify transitions
            is_intra = (ib == jb)
            is_inter_up = (ib == lcb and jb > ib)
            is_inter_down = (ib > jb and jb == lcb)

            if is_intra and not allow_intra:
                rejected_rows += 1
                continue
            if is_inter_up and not allow_inter_up:
                rejected_rows += 1
                continue
            if is_inter_down and not allow_inter_down:
                rejected_rows += 1
                continue

            ik = int(r['ik_pos'])
            iq = int(r['iq_pos'])
            im = int(r['imode'])
            enk_eV = float(r['enk_eV'])
            enkq_eV = float(r['enkq_eV'])
            omega_meV = float(r['omega_meV'])

            # Track filtered maxima
            if g_abs_meV > max_g_used:
                max_g_used = g_abs_meV
                max_g_used_info = (ik, iq, im, enk_eV, enkq_eV, omega_meV, ib, jb)

            # Energy window (in meV from CBM)
            enk_meV = 1000.0 * enk_eV
            if abs(enk_meV - E_CBM_meV) > (kbt_multiplier * kBT_meV):
                rejected_rows += 1
                continue

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

            # strict energy tolerance for interband
            if ib != jb:
                if (abs(delta_abs) > energy_tol_meV) and (abs(delta_em) > energy_tol_meV):
                    rejected_rows += 1
                    continue

            # Golden rule
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

            # NEW: accumulate Gamma_w by imode
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

            # ---- Event classification ----
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

            # ---- Dominant imode (based on largest Gamma_w contribution) ----
            im_dict = imode_gamma_acc.get(key, {})
            if len(im_dict) > 0:
                dominant_imode = max(im_dict, key=lambda m: im_dict[m])
            else:
                dominant_imode = -1

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

        # write max_g report
        with open("max_g_report.txt", "w") as fh:
            fh.write(f"max_g_unfiltered_meV = {max_g_unfiltered:.6e}\n")
            if max_g_unfiltered_info:
                fh.write(str(max_g_unfiltered_info)+"\n")
            fh.write("\nmax_g_used_meV = {max_g_used:.6e}\n")
            if max_g_used_info:
                fh.write(str(max_g_used_info)+"\n")

        return results


def write_txt(results, out_txt):
    with open(out_txt, 'w') as fh:
        fh.write("#Intraband/interband scattering results (conduction manifold only). Event = abs/em/both.\n")
        fh.write("#---------------------------------------------------------------------------------------------------------------------------------------------------\n")
        fh.write(
            f"{'#'}"
            f"{'ik_pos':>6}  {'E_k_meV':>12}  {'Gamma_s':>14}  {'tau_s':>14}  "
            f"{'Gamma_abs_s':>14}  {'Gamma_em_s':>14}  {'n_rows':>7}  "
            f"{'avg_enk_eV':>12}  {'avg_enkq_eV':>12}  "
            f"{'ibnd':>5}  {'jbnd':>5}  {'type':>10}  {'event':>7}  {'imode':>7}\n"
        )
        fh.write("#---------------------------------------------------------------------------------------------------------------------------------------------------\n")

        for (ik, Ek, G, tau, Gabs, Gem, nrows, avg_enk, avg_enkq,
             ib, jb, typ, event, imode) in results:

            fh.write(
                f"{ik:6d}  {Ek:12.6f}  {G:14.6e}  {tau:14.6e}  "
                f"{Gabs:14.6e}  {Gem:14.6e}  {nrows:7d}  "
                f"{avg_enk:12.6f}  {avg_enkq:12.6f}  "
                f"{ib:5d}  {jb:5d}  {typ:>10s}  {event:>7s}  {imode:7d}\n"
            )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('--h5', required=True)
    p.add_argument('--out', default='lcb_rates.txt')
    p.add_argument('--gamma', type=float, default=5.0)
    p.add_argument('--T', type=float, default=300.0)
    p.add_argument('--kbt-mult', type=float, default=3.0)
    p.add_argument('--cbm-ref', type=float, default=None)
    p.add_argument('--nq', type=int, default=None)
    p.add_argument('--mode', type=str, default='both',
                   choices=['intra','inter_up','inter_down','both'])
    p.add_argument('--energy-tol', type=float, default=None)
    args = p.parse_args()

    results = compute_rates(
        args.h5, gamma_meV=args.gamma, T=args.T,
        kbt_multiplier=args.kbt_mult, cbm_ref_eV=args.cbm_ref,
        nq=args.nq, mode_flag=args.mode, energy_tol_meV=args.energy_tol,
        verbose=True
    )

    write_txt(results, args.out)
    print(f"[ok] Wrote {len(results)} entries to {args.out}")

