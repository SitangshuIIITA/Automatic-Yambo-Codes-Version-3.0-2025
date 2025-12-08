#!/usr/bin/env python3
"""
new_mc_driver.py

Full-band RTA Monte Carlo driver (weighted sampling) extended to ALL bands,
with selectable sampling modes and additional diagnostics to detect and
fix angular sampling bias.

Added Option D convergence:
 - requires μ stability, diffusion linearity (R^2), and absolute μ slope tolerance
 - prints μ every print_interval

Author: adapted for your pipeline
Modifications: logging to files and mu/vel time outputs (by Sitangshu Bhattacharya)
"""
import h5py
import numpy as np
from numpy.polynomial.polynomial import Polynomial
import time
import os
import datetime

# -------------------------
# User parameters
# -------------------------
H5FILE = "epw_preproc.h5"         # change path if needed
bands_to_use = [1]               # None -> use ALL bands; or list like [0,1]
t_index = 2                      # temperature index to use
Nparticles = 20000                # number of ensemble particles
print_interval = 20000            # steps between diagnostic prints
isotropize_after_scatter = False  # set True to enforce in-plane isotropy after scattering
use_kweights_if_present = False   # use parser k-weights if available

# Convergence & stopping parameters (Option D)
tolerance = 0.01                  # relative change tolerance for stability (1%)
stability_fraction = 0.2          # last X% of mu-history used for stability check
min_dt_fs_floor = 0.5            # minimum dt in fs
seed = 1234                       # RNG seed for reproducibility
max_time_factor = 1e10            # used to bound max steps
debug_scatter = 0                 # 0=off, 1=summary per print_interval, 2=full per-scatter output

# Additional convergence knobs for Option D
conv_min_points = 20               # minimum number of μ diagnostic points before attempting conv check
conv_r2_threshold = 0.98          # require diffusion fit R^2 >= this to accept linear diffusive regime
conv_abs_mu_tol_cm2 = 0.5         # absolute μ slope tolerance (cm^2/Vs) over recent window
conv_mu_slope_window = 30          # number of recent μ points to check slope
conv_print_verbose = True         # print detailed convergence diagnostics

# Choose sampling mode: "v2_tau_boltz", "energy_only", or "binned_v2"
sampling_mode = "energy_only"     # <-- recommended to remove angular bias

# Binning parameter for "binned_v2" mode
n_energy_bins = 40

# -------------------------
# Physical constants
# -------------------------
e_charge = 1.602176634e-19  # C
kB_J = 1.380649e-23         # J/K
kB_eV = 8.617333262145e-5   # eV/K
fs_to_s = 1e-15

np.random.seed(seed)

# -------------------------
# Setup logging directories and files
# -------------------------
LOG_DIR = "log"
OUT_DIR = "MC_OUTPUT"
os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

log_path = os.path.join(LOG_DIR, "log.txt")
mu_time_path = os.path.join(OUT_DIR, "mu_time.txt")
vel_time_path = os.path.join(OUT_DIR, "vel_time.txt")

# Option B ASCII header (user selected)
ASCII_HEADER = r"""
╔══════════════════════════════════════════════╗
║      Molecular Dynamics -Mobility Engine     ║
║      Full-Band Monte Carlo Mobility Code     ║
╚══════════════════════════════════════════════╝
**************************************************************************
* Developed by: Dr. Rekha Verma and Dr. Sitangshu Bhattacharya           *
* Department of Electronics and Communication Engineering		 *
* Indian Institute of Information Technology-Allahabad                   *
* Prayagraj, Uttar Pradesh 211105, India                                 * 
**************************************************************************
"""

# open log file (append mode) and write header with timestamp + run config
log_fh = open(log_path, "w", buffering=1)  # line-buffered
def logger(message: str, to_terminal: bool = False):
    """
    Log message to file (always) and optionally print to terminal.
    Produces cleaner, spaced, short-timestamp logging.
    """
    ts = time.strftime("[%H:%M:%S]", time.gmtime())  # short timestamp

    # ensure newline
    if not message.endswith("\n"):
        message = message #+ "\n"

    # write with blank line after each entry
    log_fh.write(f"{ts} {message}\n")
    log_fh.flush()

    if to_terminal:
        print(message, end="")

def logger_raw(message: str):
    """Write a raw line to log WITHOUT timestamp."""
    if not message.endswith("\n"):
        message += "\n"
    log_fh.write(message)
    log_fh.flush()


# Start log with header and run metadata
logger_raw("="*82)
logger("MD-Mobility run log\n")

# --- ASCII HEADER WITHOUT TIMESTAMPS ---
for line in ASCII_HEADER.strip("\n").split("\n"):
    logger_raw(line)
logger_raw("")   # one blank line after header

# --- Metadata WITH timestamps ---
logger(f"Run timestamp (UTC): {datetime.datetime.utcnow().isoformat()}Z")
logger(f"HDF5 file: {H5FILE}")
logger(f"bands_to_use: {bands_to_use}, t_index: {t_index}, Nparticles: {Nparticles}")
logger(f"print_interval: {print_interval}, sampling_mode: {sampling_mode}")
logger(f"seed: {seed}")
logger_raw("="*82 + "\n")
logger_raw("\n\n")


# Prepare mu_time and vel_time files with header (overwrite previous)
with open(mu_time_path, "w") as f:
    f.write("#     time (s)           μ-xx (cm2/Vs)        μ-yy (cm2/Vs)         μ-zz (cm2/Vs)\n")
with open(vel_time_path, "w") as f:
    f.write("#     time (s)          |v| m/s\n")

# -------------------------
# Read HDF5 and build arrays
# -------------------------
with h5py.File(H5FILE, "r") as f:
    Vx_raw = f["Velocities/Vx3d"][:]
    Vy_raw = f["Velocities/Vy3d"][:]
    Vz_raw = f["Velocities/Vz3d"][:]

    if "Energies_3d" not in f:
        raise RuntimeError("Energies_3d missing from HDF5 (required).")
    E3d = f["Energies_3d"][:]   # (ntemp, nband, Nk) eV
    Tau3d = f["Tau_3d"][:]      # (ntemp, nband, Nk) ps

    # optional k-weights
    kweights = None
    if use_kweights_if_present:
        for name in ("k_weights", "kpoints/weights", "kweights"):
            if name in f:
                kweights = f[name][:]
                break

    # optional T list dataset
    T_list = f["T_Mstar"][:] if "T_Mstar" in f else None

# Squeeze dims if necessary
def _squeeze_first_dim(arr):
    arr = np.array(arr)
    if arr.ndim == 3 and arr.shape[0] == 1:
        return arr.squeeze(0)
    return arr

Vx = _squeeze_first_dim(Vx_raw)
Vy = _squeeze_first_dim(Vy_raw)
Vz = _squeeze_first_dim(Vz_raw)

if Vx.ndim != 2:
    raise RuntimeError(f"Unexpected velocity array shape: {Vx.shape}. Expect (nband, Nk).")

nband, Nk = Vx.shape

# Decide which bands to use
if bands_to_use is None:
    band_indices = list(range(nband))
else:
    band_indices = [int(b) for b in bands_to_use]
    for b in band_indices:
        if b < 0 or b >= nband:
            raise IndexError(f"band index {b} out of range [0,{nband-1}]")

n_b_use = len(band_indices)

# Build state arrays by concatenating across bands:
V_states_list = []
E_states_list = []
Tau_states_list = []
state_band = []
state_kpt = []

for b in band_indices:
    vx_b = Vx[b, :].astype(np.float64)
    vy_b = Vy[b, :].astype(np.float64)
    vz_b = Vz[b, :].astype(np.float64)
    v_b = np.stack([vx_b, vy_b, vz_b], axis=-1)   # (Nk,3)
    E_b = E3d[t_index, b, :].astype(np.float64)   # (Nk,)
    Tau_b = Tau3d[t_index, b, :].astype(np.float64)  # (Nk,)

    for ik in range(Nk):
        V_states_list.append(v_b[ik])
        E_states_list.append(E_b[ik])
        Tau_states_list.append(Tau_b[ik])
        state_band.append(b)
        state_kpt.append(ik)

V_states = np.array(V_states_list)
E_states = np.array(E_states_list)
Tau_states = np.array(Tau_states_list)
state_band = np.array(state_band, dtype=int)
state_kpt = np.array(state_kpt, dtype=int)
Nstates = V_states.shape[0]

# Build weights per state from kweights if present (k-weight is per-k typically)
if kweights is not None:
    kweights = np.array(kweights)
    if len(kweights) == Nk:
        weights_states = np.tile(kweights, n_b_use)
    elif len(kweights) == Nstates:
        weights_states = kweights.copy()
    else:
        weights_states = np.ones(Nstates)
else:
    weights_states = np.ones(Nstates)

# Filtering invalid states (e.g., zero velocity or zero tau)
vel_norm_states = np.linalg.norm(V_states, axis=1)
valid_mask = (Tau_states > 1e-12) & (vel_norm_states > 1e-9)
if not np.any(valid_mask):
    raise RuntimeError("No valid states after filtering (Tau and velocity).")

V_states = V_states[valid_mask, :]
E_states = E_states[valid_mask]
Tau_states = Tau_states[valid_mask]
weights_states = weights_states[valid_mask]
state_band = state_band[valid_mask]
state_kpt = state_kpt[valid_mask]
Nstates_valid = V_states.shape[0]

# Determine temperature
if T_list is not None:
    T_K = float(T_list[t_index])
else:
    T_K = 300.0

# -------------------------
# Sampling construction
# -------------------------
v2 = np.sum(V_states**2, axis=1)        # |v|^2 (m^2/s^2)
tau_ps = Tau_states.copy()              # ps
v2 += 1e-30
tau_ps = np.maximum(tau_ps, 1e-12)

E_ref = np.max(E_states)
boltz = np.exp(-(E_ref - E_states) / (kB_eV * T_K))

if sampling_mode == "v2_tau_boltz":
    p_state_raw = weights_states * v2 * tau_ps * boltz

elif sampling_mode == "energy_only":
    # energy-only sampling: rotationally invariant
    p_state_raw = weights_states * boltz

elif sampling_mode == "binned_v2":
    # compute bin-averaged v2 within each energy bin to enforce isotropy per shell
    bins = np.linspace(E_states.min(), E_states.max(), n_energy_bins + 1)
    v2_avg = np.zeros_like(v2)
    for i in range(n_energy_bins):
        mask = (E_states >= bins[i]) & (E_states < bins[i + 1])
        if mask.any():
            avg = np.mean(v2[mask])
            v2_avg[mask] = avg
    # for any leftover (max edge), set to mean
    zero_mask = (v2_avg == 0.0)
    if zero_mask.any():
        v2_avg[zero_mask] = np.mean(v2) if np.mean(v2) > 0 else 1.0
    p_state_raw = weights_states * v2_avg * tau_ps * boltz

else:
    raise ValueError(f"Unknown sampling_mode: {sampling_mode}")

# normalize, remove nans
p_state_raw = np.nan_to_num(p_state_raw, nan=0.0, posinf=0.0, neginf=0.0)
p_state_raw = np.clip(p_state_raw, 0.0, None)
if p_state_raw.sum() <= 0.0:
    p_state = np.ones_like(p_state_raw) / len(p_state_raw)
else:
    p_state = p_state_raw / p_state_raw.sum()

# --------------
# Diagnostics: angular histogram + per-energy-bin vx2,vy2
# --------------
angles = np.arctan2(V_states[:, 1], V_states[:, 0])  # -pi..pi
hist, bin_edges = np.histogram(angles, bins=12, density=True)
logger_raw("="*82 + "\n")
logger(f"[INFO] Total bands used = {n_b_use}, Nk = {Nk}, Nstates_valid = {Nstates_valid}, T_index = {t_index}, T_K = {T_K:.2f} K")
logger(f"[INFO] Sampling_mode = '{sampling_mode}'")
logger_raw("="*82 + "\n")
logger("[INFO] Angular histogram (12 bins) of available states (density):")
logger(str(np.round(hist, 4).tolist()))

# Top states by sampling prob (log)
logger(f"[INFO] Top 12 states by sampling probability (state_idx, band, k, p, E[eV], |v|, Tau[ps], w):")
top_idx = np.argsort(p_state)[-12:][::-1]
for s in top_idx:
    logger(f"  {s:4d}  band={state_band[s]:2d} k={state_kpt[s]:4d} p={p_state[s]:.3e} "
           f"E={E_states[s]:7.4f} |v|={np.linalg.norm(V_states[s]):7.3e} Tau={Tau_states[s]:7.3e} w={weights_states[s]:.3e}")

# Per-energy-bin vx2 vs vy2 (to detect anisotropy) — logged
nbin = min(30, max(6, int(np.sqrt(Nstates_valid)//2)))
bins = np.linspace(E_states.min(), E_states.max(), nbin + 1)
logger_raw("\n")
logger("[INFO] energy-bin statistics (bin_idx, count, <v_x^2>, <v_y^2>)")
for i in range(nbin):
    mask = (E_states >= bins[i]) & (E_states < bins[i + 1])
    if mask.sum() == 0:
        continue
    vx2 = np.mean(V_states[mask, 0] ** 2)
    vy2 = np.mean(V_states[mask, 1] ** 2)
    logger(f"  {i:2d}  {mask.sum():4d}  {vx2:9.3e}  {vy2:9.3e}")

# Adaptive dt based on min Tau among valid states (ps -> fs)
min_tau_ps = np.min(Tau_states)

# fixed dt chosen by user (floor), ignore large τ
dt_fs = min_dt_fs_floor
logger_raw("\n")
logger(f"[INFO] min Tau (ps) = {min_tau_ps:.3e} -> adaptive dt_fs = {dt_fs:.3f} fs")

# -------------------------
# Monte Carlo setup (states sampling)
# -------------------------
r = np.zeros((Nparticles, 3))                   # positions
state_indices = np.random.choice(np.arange(Nstates_valid), size=Nparticles, p=p_state)
v = V_states[state_indices, :].copy()
tau_s = Tau_states[state_indices] * 1e-12       # convert ps -> s

logger_raw("\n")
# Debug initial states logged
logger("[DEBUG] Initial state indices for first 20 particles (state_idx, band, k, p_state, E, |v|, Tau[ps]):")
for i in range(min(20, Nparticles)):
    s = state_indices[i]
    logger(f"  p_{i:2d}: state={s:4d} band={state_band[s]:2d} k={state_kpt[s]:4d} p={p_state[s]:.3e} "
           f"E={E_states[s]:7.4f} |v|={np.linalg.norm(V_states[s]):7.3e} Tau={Tau_states[s]:7.3e}")

variance = []
mu_history_x = []
mu_history_y = []
mu_history_z = []
mu_time_points = []   # time (s) corresponding to mu history points

max_steps = int(max_time_factor / dt_fs)
converged = False
logger_raw("\n")
logger(f"[INFO] Nparticles = {Nparticles}, print_interval = {print_interval}, max_steps = {max_steps}")
start_time = time.time()

# helper: compute R^2 for linear fit
def linear_fit_with_r2(x, y):
    # Fit degree-1 polynomial using numpy.polyfit for stability, return slope, intercept, r2
    if len(x) < 2:
        return 0.0, 0.0, 0.0
    A = np.vstack([x, np.ones_like(x)]).T
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]
    y_pred = m * x + c
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return m, c, r2

# -------------------------
# Monte Carlo loop
# -------------------------
for step in range(max_steps):
    # Free flight
    r += v * dt_fs * fs_to_s

    # Scattering: per-particle probabilities
    P_scatter = dt_fs * fs_to_s / tau_s
    P_scatter = np.clip(P_scatter, 0.0, 1.0)

    rand = np.random.rand(Nparticles)
    scatter_flags = rand < P_scatter
    n_scatter = int(np.count_nonzero(scatter_flags))

    if n_scatter > 0:
        # choose new state indices according to p_state (weighted)
        new_states = np.random.choice(np.arange(Nstates_valid), size=n_scatter, p=p_state)

        idx_scattered_positions = np.nonzero(scatter_flags)[0]
        v[idx_scattered_positions, :] = V_states[new_states, :]
        tau_s[idx_scattered_positions] = Tau_states[new_states] * 1e-12
        state_indices[idx_scattered_positions] = new_states

        # optional isotropize in-plane
        if isotropize_after_scatter:
            angles = 2.0 * np.pi * np.random.rand(n_scatter)
            vx = v[idx_scattered_positions, 0].copy()
            vy = v[idx_scattered_positions, 1].copy()
            speed_xy = np.sqrt(vx**2 + vy**2) + 1e-20
            v[idx_scattered_positions, 0] = speed_xy * np.cos(angles)
            v[idx_scattered_positions, 1] = speed_xy * np.sin(angles)

        # DEBUG reporting (logged, not printed)
        if debug_scatter != 0 and (step % print_interval == 0):
            E_sc = E_states[new_states]
            v_sc = np.linalg.norm(V_states[new_states], axis=1)
            tau_sc = Tau_states[new_states]
            p_sc_local = p_state[new_states]
            if debug_scatter == 1:
                logger(f"\n[SCAT_SUM] step={step}, n_scatter={n_scatter}")
                logger(f"  sample new_states[:10] = {new_states[:10].tolist()}")
                logger(f"  example (band,k) for first 10: {[ (int(state_band[s]), int(state_kpt[s])) for s in new_states[:10] ]}")
                logger(f"  E: mean={E_sc.mean():.6f} eV min={E_sc.min():.6f} max={E_sc.max():.6f}")
                logger(f"  |v| mean={v_sc.mean():.3e} m/s, Tau mean={tau_sc.mean():.3e} ps")
                local_top = np.argsort(p_sc_local)[-5:][::-1]
                logger("  Top new states (local) by p_state:")
                for jj in local_top:
                    s_idx = new_states[jj]
                    logger(f"    state={s_idx:6d} band={state_band[s_idx]:2d} k={state_kpt[s_idx]:4d} p={p_state[s_idx]:.3e} E={E_states[s_idx]:.6f} |v|={np.linalg.norm(V_states[s_idx]):.3e} Tau={Tau_states[s_idx]:.3e}")
            elif debug_scatter == 2:
                logger(f"\n[SCAT_FULL] step={step}, listing {n_scatter} scattered new_states:")
                for ii_s, s_idx in enumerate(new_states):
                    logger(f"  scatter#{ii_s:4d}: state={s_idx:6d}, band={int(state_band[s_idx])}, k={int(state_kpt[s_idx])}, p={p_state[s_idx]:.3e}, E={E_states[s_idx]:.6f}, |v|={np.linalg.norm(V_states[s_idx]):.3e}, Tau={Tau_states[s_idx]:.3e} ps")

    # Track variance & diagnostics at intervals
    
    if step % print_interval == 0:
        var_xyz = np.var(r, axis=0)
        variance.append(var_xyz)
        # time for each variance snapshot (seconds)
        t_now_s = (len(variance) - 1) * (dt_fs * fs_to_s * print_interval)
        mu_time_points.append(t_now_s)

        # compute diffusion coefficients from fit over ALL recorded variance snapshots
        if len(variance) > 2:
            var_arr = np.array(variance)
            time_s_arr = np.array(mu_time_points)

            # Fit linear slope for each direction: var = 2 D t  => D = slope/2
            sx, cx, r2x = linear_fit_with_r2(time_s_arr, var_arr[:, 0])
            sy, cy, r2y = linear_fit_with_r2(time_s_arr, var_arr[:, 1])
            sz, cz, r2z = linear_fit_with_r2(time_s_arr, var_arr[:, 2])

            D_x = 0.5 * sx
            D_y = 0.5 * sy
            D_z = 0.5 * sz

            mu_x = e_charge * D_x / (kB_J * T_K)
            mu_y = e_charge * D_y / (kB_J * T_K)
            mu_z = e_charge * D_z / (kB_J * T_K)

        else:
            # not enough points -> fallback zero
            D_x = D_y = D_z = 0.0
            mu_x = mu_y = mu_z = 0.0
            r2x = r2y = r2z = 0.0

        # append history
        mu_history_x.append(mu_x)
        mu_history_y.append(mu_y)
        mu_history_z.append(mu_z)

        # Prepare friendly units
        mu_x_cm2 = mu_x * 1e4
        mu_y_cm2 = mu_y * 1e4
        mu_z_cm2 = mu_z * 1e4

        mean_v = np.mean(np.linalg.norm(v, axis=1))

        # Log the full verbose diagnostics into log file (not spam terminal)
        logger_raw("\n")
        logger_raw("="*82)
        logger(f"Step {step:8d}")
        logger("    Variance:")
        logger(f"        var_x = {var_xyz[0]:.3e}      var_y = {var_xyz[1]:.3e}      var_z = {var_xyz[2]:.3e}")
        logger("    Averages:")
        logger(f"        mean|v| = {mean_v:.3e} m/s")
        logger(f"        scatters = {n_scatter}")
        logger(f"        avg P_scatter = {np.mean(P_scatter):.3e}")
        logger("    Mobility estimate:")
        logger(f"        μ_x = {mu_x_cm2:8.3f} cm²/V·s      μ_y = {mu_y_cm2:8.3f} cm²/V·s")
        logger(f"        R2_x = {r2x:.4f}             R2_y = {r2y:.4f}")

        # write mu_time and vel_time snapshots (append)
        with open(mu_time_path, "a") as f_mu:
            f_mu.write(f"{t_now_s:.12e}   {mu_x_cm2:.12e}   {mu_y_cm2:.12e}   {mu_z_cm2:.12e}\n")
            f_mu.flush()
        with open(vel_time_path, "a") as f_v:
            f_v.write(f"{t_now_s:.12e}   {mean_v:.12e}\n")
            f_v.flush()

        # --- Option D convergence checks ---
        conv_ok = False
        reason_msgs = []

        # need at least some history points
        if len(mu_history_x) >= conv_min_points:
            # 1) recent stability check (relative spread) on recent window
            last_n = max(3, int(len(mu_history_x) * stability_fraction))
            recent_mu_x = np.array(mu_history_x[-last_n:])
            recent_mu_y = np.array(mu_history_y[-last_n:])
            rel_change_x = (recent_mu_x.max() - recent_mu_x.min()) / (recent_mu_x.mean() + 1e-30)
            rel_change_y = (recent_mu_y.max() - recent_mu_y.min()) / (recent_mu_y.mean() + 1e-30)

            if rel_change_x < tolerance and rel_change_y < tolerance:
                reason_msgs.append(f"stability_ok (rel_change_x={rel_change_x:.4f}, rel_change_y={rel_change_y:.4f})")
                stability_ok = True
            else:
                stability_ok = False
                reason_msgs.append(f"stability_fail (rel_change_x={rel_change_x:.4f}, rel_change_y={rel_change_y:.4f})")
        else:
            stability_ok = False
            reason_msgs.append(f"stability_insufficient_points ({len(mu_history_x)}<{conv_min_points})")

        # 2) diffusion linearity: require high R^2 for x and y over the recent window of variance snapshots
        # choose last M variance snapshots (but at least 3)
        mvar = max(3, min(len(variance), int(conv_mu_slope_window)))
        if len(variance) >= mvar:
            tv = np.array(mu_time_points[-mvar:])
            var_arr_recent = np.array(variance[-mvar:])[:, :2]  # x,y
            # fit
            sx_r, cx_r, r2x_r = linear_fit_with_r2(tv, var_arr_recent[:, 0])
            sy_r, cy_r, r2y_r = linear_fit_with_r2(tv, var_arr_recent[:, 1])
            if (r2x_r >= conv_r2_threshold) and (r2y_r >= conv_r2_threshold):
                reason_msgs.append(f"diffusive_ok (r2x={r2x_r:.4f}, r2y={r2y_r:.4f})")
                diffusive_ok = True
            else:
                diffusive_ok = False
                reason_msgs.append(f"diffusive_fail (r2x={r2x_r:.4f}, r2y={r2y_r:.4f})")
        else:
            diffusive_ok = False
            reason_msgs.append(f"diffusive_insufficient_variance_points ({len(variance)}<{mvar})")

        # 3) absolute μ slope near zero over recent mu points (in cm^2/Vs)
        if len(mu_history_x) >= conv_mu_slope_window:
            last_window = conv_mu_slope_window
            t_for_mu = np.array(mu_time_points[-last_window:])
            # use mu history for x and y
            mu_x_win = np.array(mu_history_x[-last_window:]) * 1e4
            mu_y_win = np.array(mu_history_y[-last_window:]) * 1e4
            # linear slope (cm2/Vs per second)
            mx, cx_mu, r2mx = linear_fit_with_r2(t_for_mu, mu_x_win)
            my, cy_mu, r2my = linear_fit_with_r2(t_for_mu, mu_y_win)
            # convert slope from per second to per second; we compare absolute change in recent window
            # compute total change over window (slope * dt_window)
            dt_window = t_for_mu[-1] - t_for_mu[0] if t_for_mu[-1] > t_for_mu[0] else 1.0
            change_x = abs(mx * dt_window)
            change_y = abs(my * dt_window)
            if (change_x < conv_abs_mu_tol_cm2) and (change_y < conv_abs_mu_tol_cm2):
                reason_msgs.append(f"mu_slope_ok (Δμ_x={change_x:.4f} cm2/Vs, Δμ_y={change_y:.4f} cm2/Vs)")
                slope_ok = True
            else:
                slope_ok = False
                reason_msgs.append(f"mu_slope_fail (Δμ_x={change_x:.4f} cm2/Vs, Δμ_y={change_y:.4f} cm2/Vs)")
        else:
            slope_ok = False
            reason_msgs.append(f"mu_slope_insufficient_points ({len(mu_history_x)}<{conv_mu_slope_window})")

        # final combine conditions
        if stability_ok and diffusive_ok and slope_ok:
            conv_ok = True
        else:
            conv_ok = False

        if conv_print_verbose:
            logger("  [CONV] checks: " + "; ".join(reason_msgs))

        if conv_ok:
            mu_2D_avg = (mu_x + mu_y) / 2.0
            mu_x_cm2_final = mu_x * 1e4
            mu_y_cm2_final = mu_y * 1e4
            logger_raw("\n")
            logger_raw("="*82 + "\n")
            logger(f"[CONV] Converged at step {step}. μ_x = {mu_x_cm2_final:7.2f} cm²/V·s, μ_y = {mu_y_cm2_final:7.2f} cm²/V·s")
            logger(f"[CONV] In-plane average μ_2D = {mu_2D_avg*1e4:7.2f} cm²/V·s")
            logger_raw("\n")
            logger_raw("="*82 + "\n")
            converged = True
            # Print brief terminal message
            print(f"[CONV] Converged at step {step}. μ_2D = {mu_2D_avg*1e4:.2f} cm^2/Vs")
            break

# -------------------------
# Final results (if not converged earlier)
# -------------------------
if not converged:
    if len(variance) >= 2:
        var_arr = np.array(variance)
        time_s = np.array(mu_time_points)
        sx, cx, r2x = linear_fit_with_r2(time_s, var_arr[:, 0])
        sy, cy, r2y = linear_fit_with_r2(time_s, var_arr[:, 1])
        sz, cz, r2z = linear_fit_with_r2(time_s, var_arr[:, 2])
        D_x = 0.5 * sx
        D_y = 0.5 * sy
        D_z = 0.5 * sz
        mu_x = e_charge * D_x / (kB_J * T_K)
        mu_y = e_charge * D_y / (kB_J * T_K)
        mu_z = e_charge * D_z / (kB_J * T_K)
    else:
        D_x = D_y = D_z = 0.0
        mu_x = mu_y = mu_z = 0.0

    mu_2D_avg = (mu_x + mu_y) / 2.0

elapsed = time.time() - start_time
# Convert diffusion tensor to cgs (cm^2/s)
D_x_cm2 = D_x * 1e4
D_y_cm2 = D_y * 1e4
D_z_cm2 = D_z * 1e4

# Convert mobility tensor to cgs (cm^2/Vs)
mu_x_cm2 = mu_x * 1e4
mu_y_cm2 = mu_y * 1e4
mu_z_cm2 = mu_z * 1e4

mu_2D_avg_cm2 = mu_2D_avg * 1e4


# ==== Drude mobility computation using parsed effective mass ====
m_e = 9.10938356e-31  # electron mass (kg)

# Read the effective mass for the selected temperature
with h5py.File(H5FILE, "r") as f:
    if "Mstar_T" in f:
        mstar_me = float(f["Mstar_T"][t_index])   # unit: m_e
    else:
        mstar_me = None

if mstar_me is not None:
    mstar = mstar_me * m_e

    # average scattering time τ (ps → s)
    tau_avg_s = float(np.mean(Tau_states) * 1e-12)

    # Drude mobility μ = e τ / m*
    mu_drude_si = e_charge * tau_avg_s / mstar
    mu_drude_cm2 = mu_drude_si * 1e4
else:
    mu_drude_si = None
    mu_drude_cm2 = None

# =============================
# Prepare CGS unit conversions (already computed above)
# =============================

# Print and log final results (concise terminal + verbose log)
logger_raw("\n")
logger_raw("="*82 + "\n")
logger_raw("\n")
print("--- Final Results (brief) ---")
print(f"Elapsed time: {elapsed:.1f} s")
print(f"In-plane average μ_2D (MC) = {mu_2D_avg_cm2:.2f} cm^2/(V·s)")
if mu_drude_cm2 is not None:
    print(f"Drude μ (arithmetic mean) = {mu_drude_cm2:.3f} cm^2/(V·s)")
else:
    print("Drude μ not available (Mstar_T missing)")

# Full verbose logging of final results
logger_raw("\n")
logger("--- Final Results ---")
logger(f"Elapsed time: {elapsed:.1f} s")
logger(f" Diffusion tensor (SI) D_xx, D_yy, D_zz = {D_x:.3e}, {D_y:.3e}, {D_z:.3e}  m^2/s")
logger(f" Diffusion tensor (CGS) D_xx, D_yy, D_zz = {D_x_cm2:.3e}, {D_y_cm2:.3e}, {D_z_cm2:.3e}  cm^2/s")
logger(f" Mobility tensor (SI) μ_x, μ_y, μ_z = {mu_x:.3e}, {mu_y:.3e}, {mu_z:.3e}  m^2/(V·s)")
logger(f" Mobility tensor (CGS) μ_x, μ_y, μ_z = {mu_x_cm2:.3f}, {mu_y_cm2:.3f}, {mu_z_cm2:.3f}  cm^2/(V·s)")
logger(f"In-plane average μ_2D (MC) = {mu_2D_avg_cm2:.2f} cm^2/(V·s)")
logger_raw("\n")
logger_raw("="*82 + "\n")
logger_raw("\n")

if mu_drude_si is not None:
    logger("\n--- Drude mobility (arithmetic mean tau) ---")
    logger(f"Drude μ (SI) = {mu_drude_si:.3e} m^2/(V·s)")
    logger(f"Drude μ (CGS)= {mu_drude_cm2:.3f} cm^2/(V·s)")
else:
    logger("Drude mobility not available (Mstar_T missing in HDF5).")

# ---- Drude: weighted variants (if mstar exists) ----
if 'mstar' in locals():
    SI_to_CGS = 1e4
    weights_for_avg = p_state if 'p_state' in globals() else np.ones_like(Tau_states)
    v2_states = np.sum(V_states**2, axis=1)     # |v|^2 (m^2/s^2)
    vx2_states = V_states[:,0]**2               # v_x^2
    den_v2 = np.sum(weights_for_avg * v2_states) + 1e-30
    den_vx2 = np.sum(weights_for_avg * vx2_states) + 1e-30
    tau_v2_weighted = np.sum(weights_for_avg * Tau_states * v2_states) / den_v2
    tau_vx2_weighted = np.sum(weights_for_avg * Tau_states * vx2_states) / den_vx2
    tau_v2_weighted_s = tau_v2_weighted * 1e-12
    tau_vx2_weighted_s = tau_vx2_weighted * 1e-12
    mu_drude_v2_si = e_charge * tau_v2_weighted_s / mstar
    mu_drude_v2_cgs = mu_drude_v2_si * SI_to_CGS
    mu_drude_vx2_si = e_charge * tau_vx2_weighted_s / mstar
    mu_drude_vx2_cgs = mu_drude_vx2_si * SI_to_CGS

    logger_raw("\n")
    logger_raw("="*82 + "\n")
    logger_raw("\n")
    logger("--- Weighted Drude estimates ---")
    logger(f"tau (arithmetic mean)    = {tau_avg_s:.3e} s (used earlier)")
    logger(f"tau (v^2-weighted)       = {tau_v2_weighted_s:.3e} s")
    logger(f"tau (v_x^2-weighted)     = {tau_vx2_weighted_s:.3e} s")
    logger(f"Drude μ (v^2-weighted)   = {mu_drude_v2_si:.3e} m^2/(V·s) = {mu_drude_v2_cgs:.3f} cm^2/(V·s)")
    logger(f"Drude μ (v_x^2-weighted) = {mu_drude_vx2_si:.3e} m^2/(V·s) = {mu_drude_vx2_cgs:.3f} cm^2/(V·s)")
    logger_raw("\n")
    logger_raw("="*82 + "\n")
    logger_raw("\n")
    
else:
    logger("Mstar_T not found in HDF5; cannot compute Drude weighted variants (missing m*).")

# Close log file handle
log_fh.close()

