#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def plot_emc_results(npzfile='fb_emc_results.npz', outfile='drift_field_electrons.png'):

    data = np.load(npzfile, allow_pickle=True)

    E_fields = data["E_fields"]
    drifts   = data["drifts"]
    meanEs   = data["meanEs"]

    # ---------- STYLE SETTINGS ----------
    plt.rcParams.update({
        'font.size': 22,
        'axes.linewidth': 1.4,
        'font.family': 'serif',
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.top': True,
        'xtick.bottom': True,
        'ytick.right': True,
        'figure.dpi': 300,
        'lines.linewidth': 2.5,
        'lines.markersize': 10
    })

    # ---------- PLOTTING ----------
    plt.figure(figsize=(14, 6))

    # --- LEFT PANEL: Drift velocity ---
    ax1 = plt.subplot(1, 2, 1)
    ax1.semilogx(E_fields, drifts,
                 marker='o',
                 color='tab:blue',
                 markerfacecolor='white',
                 markeredgewidth=2.5)

    ax1.set_ylim(top = drifts.max() * 1.05)
    ax1.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
    ax1.tick_params(which='major', width=1.2, length=6)
    ax1.tick_params(which='minor', width=1.0, length=4)

    ax1.set_xlabel(r'Electric field |E| (Vm$^{-1}$)')
    ax1.set_ylabel(r'Drift velocity |v$\mathrm{_d}$| (ms$^{-1}$)')

    # --- RIGHT PANEL: Mean energy ---
    ax2 = plt.subplot(1, 2, 2)
    ax2.semilogx(E_fields, meanEs,
                 marker='o',
                 color='tab:red',
                 markerfacecolor='white',
                 markeredgewidth=2.5)

    ax2.set_ylim(top = meanEs.max() * 1.05)   
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
    ax2.tick_params(which='major', width=1.2, length=6)
    ax2.tick_params(which='minor', width=1.0, length=4)

    ax2.set_xlabel(r'Electric field |E| (Vm$^{-1}$)')
    ax2.set_ylabel(r'Mean energy above CBM (meV)')

    plt.tight_layout()
    plt.savefig(outfile)

    print(f"[PLOT] Saved {outfile}")

def plot_vx_time(vx_file, steady_file, out_png):
    """
    Plot v_x(t) and mark the steady-state time t_ss, the steady-state region,
    and a steady-state ± std band (with std printed in the legend).
    Time is plotted in ps (log scale). Velocities in km/s.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # ---------- STYLE ----------
    plt.rcParams.update({
        'font.size': 20,
        'axes.linewidth': 1.4,
        'font.family': 'serif',
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.top': True,
        'xtick.bottom': True,
        'ytick.right': True,
        'figure.dpi': 300,
        'lines.linewidth': 2.5,
        'lines.markersize': 10
    })

    # ---------- Load vx(t) ----------
    data = np.loadtxt(vx_file)
    t = data[:, 0]
    v = data[:, 1]

    # convert units
    t_ps = t * 1e12      # ps
    v_km = v / 1000.0    # km/s

    # ---------- Load steady-state info ----------
    t_ss = v_ss = v_d = None
    with open(steady_file, "r") as f:
        for line in f:
            if "steady state time in s" in line:
                try:
                    t_ss = float(line.split(":")[1])
                except:
                    pass
            if "steady state vx in ms^-1" in line:
                try:
                    v_ss = float(line.split(":")[1])
                except:
                    pass
            if "drift velocity in ms^-1" in line:
                try:
                    v_d = float(line.split(":")[1])
                except:
                    pass

    # safe conversions for plotting
    t_ss_ps = t_ss * 1e12 if t_ss is not None else None
    v_ss_km = v_ss / 1000.0 if v_ss is not None else None
    v_d_km  = v_d  / 1000.0 if v_d is not None else None

    # ============================================================
    # Compute steady-state std: use only plateau points after t_ss
    # Exclude abrupt collapse by thresholding near v_ss
    # ============================================================
    v_std = None
    if t_ss_ps is not None and v_ss_km is not None:
        # basic mask for times after t_ss
        mask_after = t_ps >= t_ss_ps

        # if many points after, refine to exclude the collapse: keep points within
        # a small window around v_ss (e.g., above v_ss - 0.02 km/s)
        candidate = v_km[mask_after]
        times_candidate = t_ps[mask_after]

        if len(candidate) > 0:
            # threshold value: 0.02 km/s below v_ss (tunable)
            thresh = 0.02
            if v_ss_km is not None:
                refined_mask = candidate > (v_ss_km - thresh)
                v_ss_region = candidate[refined_mask]
                t_region = times_candidate[refined_mask]
            else:
                # fallback: use candidate directly
                v_ss_region = candidate
                t_region = times_candidate

            # final std if we have enough samples
            if len(v_ss_region) > 5:
                v_std = float(v_ss_region.std())
            elif len(v_ss_region) > 0:
                v_std = float(v_ss_region.std())  # small sample: still use it
            else:
                v_std = 0.0

            # store masks for plotting the band restricted to the plateau times
            ss_mask_for_fill = (t_ps >= t_ss_ps) & (v_km > (v_ss_km - thresh))
        else:
            v_std = 0.0
            ss_mask_for_fill = (t_ps >= t_ss_ps)
    else:
        ss_mask_for_fill = (t_ps >= (t_ss_ps if t_ss_ps is not None else t_ps.max()+1))

    # ---------- Begin plot ----------
    plt.figure(figsize=(12, 6))
    plt.plot(t_ps, v_km, lw=1.8, label=r"$v_x(t)$")

    # ---------- t_ss vertical marker + red dot ----------
    if t_ss_ps is not None:
        plt.axvline(t_ss_ps, color="red", ls="--", lw=1.4, label=fr"t$_\mathrm{{ss}}$ ≈ {t_ss_ps:.2f} ps")
        # red dot at (t_ss, v_ss) — only if v_ss exists
        if v_ss_km is not None:
            plt.plot(t_ss_ps, v_ss_km, "o", color="red", markersize=7)

        # labelled text, offset a bit to the right
        plt.text(t_ss_ps * 1.12, v_km.max() * 0.94, r"$t_{\mathrm{ss}}$",
                 color="red", ha="left", va="center", fontsize=20)

    # ---------- Steady-state region shading (time) ----------
    if t_ss_ps is not None:
        plt.axvspan(t_ss_ps, t_ps.max(), color="orange", alpha=0.12,
                    label="Steady-state region")

    # ---------- ± std band (only in plateau times) + dashed mean line ----------
    if (v_std is not None) and (v_ss_km is not None) and (v_std > 0.0):
        # fill_between over the restricted time slice so band is horizontal
        if np.any(ss_mask_for_fill):
            plt.fill_between(
                t_ps[ss_mask_for_fill],
                (v_ss_km - v_std),
                (v_ss_km + v_std),
                color="purple",
                alpha=0.9,
                label = r"$\sigma_{\rm ss}$ = $\pm$ " + f"{v_std:.3f} " r"km s$^{-1}$"

            )

        # dashed mean line (across full x for clarity)
        plt.axhline(v_ss_km, color="green", ls="--", lw=1.8)

    # ---------- v_ss & v_d labeled horizontal lines (legend entries) ----------
    if v_ss_km is not None:
        plt.axhline(v_ss_km, color="green", ls="--", lw=1.2,
                    label=fr"v$_\mathrm{{ss}}$ ≈ {v_ss_km:.2f} km s$^{{-1}}$")

    if v_d_km is not None:
        plt.axhline(v_d_km, color="blue", ls=":", lw=1.6, alpha=0.7,
                    label=fr"v$_\mathrm{{d}}$ ≈ {v_d_km:.2f} km s$^{{-1}}$")

    # ---------- Labels and scale ----------
    plt.xlabel("Time (ps)")
    plt.ylabel(r"Velocity $\mathrm{v}_\mathrm{x}$ (km s$^{-1}$)")
    plt.xscale("log")

    # ---------- Legend (inside plot, lower-left) ----------
    plt.legend(loc="lower left", framealpha=0.85, fontsize=16)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f"[PLOT] Saved {out_png}")


if __name__ == "__main__":
    # === 1. Plot drift velocity vs field ===
    plot_emc_results()

    # === 2. Plot vx(t) for each E-field directory ===
    import os

    base = "Drift_out"

    if not os.path.isdir(base):
        print("[WARN] Drift_out folder not found.")
    else:
        for folder in os.listdir(base):
            sub = os.path.join(base, folder)
            if not os.path.isdir(sub):
                continue

            vx_file = os.path.join(sub, "vx_time.txt")
            steady_file = os.path.join(sub, "steady_state.txt")
            out_png = os.path.join(sub, "vx_time.png")

            if os.path.exists(vx_file) and os.path.exists(steady_file):
                print(f"[PLOT] Generating vx_time for {folder}")
                plot_vx_time(vx_file, steady_file, out_png)


