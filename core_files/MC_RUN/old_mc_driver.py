#!/usr/bin/env python3
"""
fbmc_rta_driver_live_fixed.py

Live RTA full-band Monte Carlo driver with corrected Einstein mobility.
Uses epw_preproc.h5 produced by your parser.

Author: Adapted for your pipeline
"""

import numpy as np
import h5py
import time
from math import isfinite
import pickle
import matplotlib.pyplot as plt

# --- Physical constants ---
e_charge = 1.602176634e-19       # C
kB = 1.380649e-23                # J/K
m_e = 9.10938356e-31             # kg
eV_to_J = 1.602176634e-19
fs_to_s = 1e-15
fs_to_ps = 1e-3

# -------------------------
# Helper: load preprocessed HDF5
# -------------------------
def load_preproc(fname="epw_preproc.h5"):
    """Load EPW preprocessed data and print sample values for inspection.

    NOTE: velocities are expected already in m/s in the HDF5 file (parser should convert).
    """
    with h5py.File(fname, "r") as f:
        E3d = f["Energies_3d"][:]    # (ntemp, nband, Nk)
        Gamma3d = f["Gamma_3d"][:]   # (ntemp, nband, Nk)
        Tau3d = f["Tau_3d"][:]       # (ntemp, nband, Nk)
        Vx = f["Velocities/Vx3d"][:] # (1, nband, Nk) typically
        Vy = f["Velocities/Vy3d"][:]
        Vz = f["Velocities/Vz3d"][:]
        T_Mstar = f["T_Mstar"][:] if "T_Mstar" in f else None
        Mstar_T = f["Mstar_T"][:] if "Mstar_T" in f else None
    # --- Compute and print thermal velocity using first T index ---
    if T_Mstar is not None and Mstar_T is not None and len(T_Mstar) > 0 and len(Mstar_T) > 0:
        T0 = float(T_Mstar[0])
        m_eff = float(Mstar_T[0])  # in units of m_e
        v_th = np.sqrt(3 * kB * T0 / (m_eff * m_e))
        print(f"\n[DEBUG] Thermal velocity at {T0:.2f} K (m* = {m_eff:.3f} m_e): v_th = {v_th:.3e} m/s\n")
    else:
        print("\n[DEBUG] Thermal velocity not available — missing T_Mstar or Mstar_T\n")



    # --- Handle shape ---
    Vx_s = Vx.squeeze(0) if Vx.ndim == 3 and Vx.shape[0] == 1 else Vx
    Vy_s = Vy.squeeze(0) if Vy.ndim == 3 and Vy.shape[0] == 1 else Vy
    Vz_s = Vz.squeeze(0) if Vz.ndim == 3 and Vz.shape[0] == 1 else Vz
    V = np.stack([Vx_s, Vy_s, Vz_s], axis=-1)  # (nband, Nk, 3)

    # --- IMPORTANT: DO NOT re-convert velocities here ---
    # Parser should already have converted velocities to m/s.
    # (If you ever change parser, remove this comment and ensure units match.)

    # --- Basic sanity check on velocity magnitudes ---
    mean_speed = np.mean(np.linalg.norm(V.reshape(-1,3), axis=1))
    print(f"[INFO] Mean |V| across dataset = {mean_speed:.3e} m/s")
    if mean_speed < 1e-2:
        print("[WARN] Mean velocity magnitude is extremely small (<1e-2). "
              "This might indicate velocities are still in a.u. or parsing failed.")
    elif mean_speed < 1e3:
        print("[WARN] Mean |V| is below 1e3 m/s — check units and parser conversion.")
    # Typical expected mean: ~1e5 - 1e6 m/s for many 2D materials.

    T_list = np.array(T_Mstar) if T_Mstar is not None else np.arange(E3d.shape[0])

    # --- Debug print (use first temperature index) ---
    t_index = 0
    print(f"\n[DEBUG] Data summary for T = {float(T_list[t_index]):.2f} K")
    print("Format: Band idx | Energy [eV] | Gamma [1/ps] | Velocities [m/s]")
    nb_sample = min(5, E3d.shape[1])
    Nk = E3d.shape[2]
    for b in range(nb_sample):
        for k in range(min(5, Nk)):
            energy = E3d[t_index, b, k]
            gamma = Gamma3d[t_index, b, k]
            vel = V[b, k]
            print(f"Band {b+1}, kpt {k+1}: E = {energy:.4f} | Gamma = {gamma:.4e} | "
                  f"V = [{vel[0]:.2e}, {vel[1]:.2e}, {vel[2]:.2e}]")

    # --- Optional: Plot mean Vx per band ---
    try:
        mean_vx = [np.mean(V[b, :, 0]) for b in range(V.shape[0])]
        plt.figure(figsize=(6,3.5))
        plt.plot(np.arange(1, len(mean_vx)+1), mean_vx, "o-", lw=1.5)
        plt.xlabel("Band index")
        plt.ylabel("⟨Vx⟩ (m/s)")
        plt.title(f"Mean Vx per band at T = {float(T_list[t_index]):.1f} K")
        plt.grid(True)
        plt.tight_layout()
        plt.show()
    except Exception:
        pass

    return {
        "E3d": E3d,
        "Gamma3d": Gamma3d,
        "Tau3d": Tau3d,
        "V": V,
        "T_list": T_list,
        "Mstar_T": Mstar_T
    }

# -------------------------
# Boltzmann weights utility
# -------------------------
def boltzmann_weights(E_eV, T_K):
    E_j = E_eV * eV_to_J
    beta = 1.0 / (kB * T_K)
    E_min = np.min(E_j)
    vals = np.exp(-beta * (E_j - E_min))
    flat = vals.ravel()
    s = flat.sum()
    return flat / s if s != 0 else np.ones_like(flat)/flat.size

def sample_states_from_weights(weights, n):
    return np.random.choice(len(weights), size=n, p=weights)

# -------------------------
# RTA Monte Carlo core (live)
# -------------------------
def run_rta_mc_live(preproc_h5="epw_preproc.h5",
                    itemp_index=0,
                    T_K=None,
                    E_fields_kVcm=[0.1, 0.5, 1.0, 2.0, 5.0],
                    Nparticles=800,
                    dt_fs=0.2,
                    t_tot_fs=1000000.0,   # Use long simulation time for proper mu_Einstein
                    n_cut=40,       # energy cutoff above CBM
                    seed=1,
                    live_update_interval_ps=5.0,
                    steady_state_rel_tol=1e-2,
                    log_file=None):

    def log(msg):
        print(msg)
        if log_file is not None:
            with open(log_file, "a") as f:
                f.write(msg + "\n")

    np.random.seed(seed)
    db = load_preproc(preproc_h5)
    E3d, Gamma3d, Tau3d, V, T_list, Mstar_T = db["E3d"], db["Gamma3d"], db["Tau3d"], db["V"], db["T_list"], db["Mstar_T"]

    ntemp, nband, Nk = E3d.shape
    if itemp_index < 0 or itemp_index >= ntemp:
        raise ValueError(f"itemp_index {itemp_index} out of range (0..{ntemp-1})")

    T_run = float(T_K) if T_K is not None else float(T_list[itemp_index])

    # --- flatten arrays for this temperature ---
    E = E3d[itemp_index,:,:]
    Gamma = Gamma3d[itemp_index,:,:]
    flatE = E.ravel()
    flatV = V.reshape((nband*Nk, 3))
    flatGamma = Gamma.ravel()

    # --- energy cutoff filter ---
    CBM = np.min(flatE)		# assume conduction band minimum
    E_cutoff_eV = n_cut * kB * T_run / e_charge
    valid_idx = np.where(flatE <= CBM + E_cutoff_eV)[0]

    
    if len(valid_idx) == 0:
        raise RuntimeError("No states left after applying energy cutoff. Increase E_cutoff_eV.")

    # Filter arrays
    flatE = flatE[valid_idx]
    flatV = flatV[valid_idx]
    flatGamma = flatGamma[valid_idx]
    weights = boltzmann_weights(flatE.reshape((1, -1)), T_run).ravel()
    weights /= np.sum(weights)  # normalize

    log(f"Filtered states by energy cutoff: {len(valid_idx)} states kept (E <= {CBM + E_cutoff_eV:.3f} eV)")

    # --- Print sample velocities, energies, and scattering rates ---
    print("\nSample data at T = {:.2f} K (after energy filter):".format(T_run))
    print("Format: Band idx | Energy [eV] | Gamma [1/ps] | Velocities [m/s]")
    nb_sample = min(5, len(flatE))
    for i in range(nb_sample):
        vel = flatV[i]
        print(f"State {i}: E = {flatE[i]:.4f} | Gamma = {flatGamma[i]:.4e} | V = [{vel[0]:.4e}, {vel[1]:.4e}, {vel[2]:.4e}]")

    # --- Monte Carlo simulation setup ---
    dt_s = dt_fs * fs_to_s
    dt_ps = dt_fs * 1e-3
    P_state = 1.0 - np.exp(-flatGamma * dt_ps)
    nsteps = int(max(1, int(np.round(t_tot_fs / dt_fs))))
    live_interval_steps = max(1, int(round((live_update_interval_ps * 1e3) / dt_fs)))
    results = {}

    log(f"\nStarting RTA-MC (live) at T = {T_run:.2f} K (itemp_index={itemp_index})")
    if Mstar_T is not None and len(Mstar_T) > itemp_index:
        log(f"  M*(T) from preproc = {float(Mstar_T[itemp_index]):.3f} m_e")
    else:
        log("  No M*(T) available in HDF5 (Drude checks limited).")

    # --- Initialize live plot ---
    plt.ion()
    fig, ax_plot = plt.subplots()
    ax_plot.set_xlabel("Time (ps)")
    ax_plot.set_ylabel("Drift velocity <v_x> (m/s)")
    ax_plot.set_title("Drift velocity vs time")
    line, = ax_plot.plot([], [], lw=2, label="v_drift(t)")
    steady_line = ax_plot.axhline(0, color='r', linestyle='--', label="Steady-state")
    ax_plot.legend()
    plt.show()

    mu_Einstein_cm2Vs = None

    for iE, E_kVcm in enumerate(E_fields_kVcm):
        E_field_Vpm = float(E_kVcm) * 1e5
        log("\n" + "-"*70)
        log(f"Field = {E_kVcm:.3f} kV/cm. Ensemble N = {Nparticles}, dt = {dt_fs} fs, total t = {t_tot_fs} fs")

        p_idx = sample_states_from_weights(weights, Nparticles)
        pos = np.zeros((Nparticles, 2))
        vel = flatV[p_idx].copy()

        mean_vx_history = np.zeros(nsteps)
        mean_x_history = np.zeros(nsteps)
        mean_x2_history = np.zeros(nsteps)

        start_time_wall = time.time()
        for step in range(nsteps):
            # --- scattering ---
            r = np.random.random(Nparticles)
            scattered = r < P_state[p_idx]
            if scattered.any():
                n_sc = np.count_nonzero(scattered)
                new_states = sample_states_from_weights(weights, n_sc)
                p_idx[scattered] = new_states
                vel[scattered] = flatV[p_idx[scattered]]

            # --- propagate ---
            pos += vel[:, :2] * dt_s

            # --- acceleration due to field ---
            if E_field_Vpm != 0.0 and Mstar_T is not None and len(Mstar_T) > itemp_index and isfinite(float(Mstar_T[itemp_index])):
                m_eff_kg = float(Mstar_T[itemp_index]) * m_e
                vel[:, 0] += (-e_charge * E_field_Vpm) / m_eff_kg * dt_s

            # --- statistics ---
            mean_vx_history[step] = np.mean(vel[:,0])
            mean_x_history[step] = np.mean(pos[:,0])
            mean_x2_history[step] = np.mean(pos[:,0]**2)

            # --- live plotting ---
            if step % live_interval_steps == 0 or step == nsteps-1:
                times_ps = np.arange(step+1) * dt_fs * 1e-3
                line.set_data(times_ps, mean_vx_history[:step+1])
                if step > int(0.5 * nsteps):
                    steady_v = np.mean(mean_vx_history[int(0.9*nsteps):step+1])
                    steady_line.set_ydata([steady_v, steady_v])
                ax_plot.relim()
                ax_plot.autoscale_view()
                plt.pause(0.001)

        total_wall = time.time() - start_time_wall
        log(f"Finished field {E_kVcm:.3f} kV/cm. Wall time {total_wall:.1f}s")

        # --- diffusion & mobilities ---
        times_s = np.arange(nsteps) * dt_fs * fs_to_s
        var_x = mean_x2_history - mean_x_history**2
        half = int(0.5 * nsteps)
        if times_s[half:].size >= 2:
            slope, _ = np.polyfit(times_s[half:], var_x[half:], 1)
            D_x = slope / 2.0
        else:
            D_x = 0.0

        # Compute mu_Einstein only once
        if mu_Einstein_cm2Vs is None and D_x > 0:
            mu_Einstein = (e_charge * D_x) / (kB * T_run)
            mu_Einstein_cm2Vs = mu_Einstein * 1e4
            log(f"mu_Einstein (from {E_kVcm} kV/cm) = {mu_Einstein_cm2Vs:.4e} cm²/Vs")

        # Drift mobility
        start_idx = max(int(0.9 * nsteps), 0)
        v_slice = mean_vx_history[start_idx:nsteps]
        v_d = float(np.mean(v_slice)) if len(v_slice) > 0 else 0.0
        mu_drift_cm2Vs = (v_d / E_field_Vpm * 1e4) if (E_field_Vpm != 0.0) else None

        # Drude mobility
        flatTau_ps = 1.0 / (flatGamma + 1e-50)
        tau_eff_ps = float(np.sum(flatTau_ps * weights))
        mu_drude_cm2Vs = None
        if Mstar_T is not None and len(Mstar_T) > itemp_index:
            m_eff_kg = float(Mstar_T[itemp_index]) * m_e
            mu_drude_cm2Vs = (e_charge * tau_eff_ps*1e-12 / m_eff_kg) * 1e4

        results[E_kVcm] = {
            "D_x": D_x,
            "mu_Einstein": mu_Einstein_cm2Vs,
            "v_d": v_d,
            "mu_drift": mu_drift_cm2Vs,
            "mu_drude": mu_drude_cm2Vs,
            "tau_eff_ps": tau_eff_ps,
            "mean_vx_history": mean_vx_history,
            "mean_x_history": mean_x_history,
            "nsteps": nsteps,
            "dt_fs": dt_fs
        }

        mu_drift_str = f"{mu_drift_cm2Vs:.4e}" if mu_drift_cm2Vs is not None else "None"
        mu_ein_str = f"{mu_Einstein_cm2Vs:.4e}" if mu_Einstein_cm2Vs is not None else "None"
        mu_drude_str = f"{mu_drude_cm2Vs:.4e}" if mu_drude_cm2Vs is not None else "None"
        log(f"Summary (field {E_kVcm} kV/cm): v_d = {v_d:.4e} m/s | "
            f"mu_drift = {mu_drift_str} cm²/Vs | "
            f"mu_Einstein = {mu_ein_str} cm²/Vs | mu_drude = {mu_drude_str}")
        log(f" tau_eff (ps) = {tau_eff_ps:.4e} ps")

    plt.ioff()
    plt.show()
    return {"T_run": T_run, "itemp_index": itemp_index, "results": results}


# -------------------------
# Main block
# -------------------------
if __name__ == "__main__":
    itemp_index = 0
    explicit_TK = None
    E_fields_kVcm = [0.0, 0.1]
    Nparticles = 1000
    dt_fs = 0.2
    t_tot_fs = 1000000.0   # <- Long enough for proper mu_Einstein
    seed = 1234
    log_file = "log.txt"

    with open(log_file, "w") as f:
        f.write(f"RTA-MC run started at {time.ctime()}\n\n")

    out = run_rta_mc_live(preproc_h5="epw_preproc.h5",
                         itemp_index=itemp_index,
                         T_K=explicit_TK,
                         E_fields_kVcm=E_fields_kVcm,
                         Nparticles=Nparticles,
                         dt_fs=dt_fs,
                         t_tot_fs=t_tot_fs,
                         seed=seed,
                         live_update_interval_ps=20.0,
                         steady_state_rel_tol=1e-2,
                         log_file=log_file)

    print("\nAll done. Field sweep results:")
    for E_k, res in out["results"].items():
        mu_ein = res['mu_Einstein']
        mu_drift = res['mu_drift']
        mu_drude = res['mu_drude']
        mu_ein_s = f"{mu_ein:.4e}" if mu_ein is not None else "None"
        mu_drift_s = f"{mu_drift:.4e}" if mu_drift is not None else "None"
        mu_drude_s = f"{mu_drude:.4e}" if mu_drude is not None else "None"
        print(f" Field {E_k} kV/cm : mu_Einstein = {mu_ein_s} cm²/Vs ; mu_drift = {mu_drift_s} cm²/Vs ; mu_drude = {mu_drude_s} cm²/Vs")

    results_file = "rta_mc_results.pkl"
    with open(results_file, "wb") as f:
        pickle.dump(out, f)
    print(f"\nResults saved to {results_file}")

