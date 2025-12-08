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


if __name__ == "__main__":
    plot_emc_results()

