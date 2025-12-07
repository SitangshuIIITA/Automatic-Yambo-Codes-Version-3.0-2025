#!/usr/bin/env python3
###################################################################################
# Combined spectral decomposition and integrated scattering rate (Electron + Hole)
# Top: Electron, Bottom: Hole
# Sitangshu — Oct 2025 Originally developed by:
###################################################################################
# Script to produce the spectral decomposition and integrated scattering rate.
# The spectral decomposition is computed at an energy 3/2 k_B T away from the band edge. 
# Samuel Ponc\'e 
# Version1 - Fall 2018
# Version2 - Fall 2020
# Version3 - November 2022
# If you are using this script, please consider citing S. Ponc\'e et al, Phys. Rev. Research 3, 043022 (2021)
###################################################################################
###################################################################################

import os
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 20,
    'axes.linewidth': 1.2,
    'font.family': 'serif',
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'figure.dpi': 300,
    'lines.markersize': 7,
    'legend.handlelength': 2.6
})

kb = 8.617333262145E-5  # eV/K
Ry2meV = 13605.6980659
hbar = 6.582119514E-16  # eV·s
fact = 1.0 / (1000 * hbar * 1E12)  # → THz

###################################################################################
# This needs to be manually updated and corresponds to the EPW calculation
def process_data(file1, label):
    mob_maxfreq = 160	# meV
    mob_nfreq = 640
    T = 300	# Temperature in K
    c = 0.80	# Gaussian broadening in meV. Should be smaller in practice. 
    nkpt = 425	# Number of k-points - same as in file1
    nbnd = 2	# Number of bands - same as in file1
    tol = 10.0	# meV. This means we are taking states around 3/2kbT +- tol in meV. In practice it should be quite small.

    # placeholder (unfilled) value in Rydberg units (will be converted)
    placeholder = -1000.0

    inv_tau_freq_tmp = np.zeros((nkpt, nbnd, mob_nfreq, 2))
    inv_tau_freq_tmp[:, :, :, 0] = placeholder
    step = mob_maxfreq / mob_nfreq

    with open(file1, 'r') as W:
        for line in W:
            tmp = line.split()
            if not tmp or tmp[0].startswith('#'):
                continue
            ikpt = int(tmp[0]) - 1
            ibnd = int(tmp[1]) - 1
            # column 3 in the file is the phonon frequency (meV in your format)
            freq = float(tmp[3])
            ifreq = int(freq / step) - 1
            # guard ifreq bounds
            if ifreq < 0 or ifreq >= mob_nfreq:
                continue
            # store energy and scattering rate (both in Ry in input assumed)
            inv_tau_freq_tmp[ikpt, ibnd, ifreq, 0] = float(tmp[2])
            inv_tau_freq_tmp[ikpt, ibnd, ifreq, 1] = float(tmp[4])

    # convert to meV
    inv_tau_freq_tmp[:, :, :, 0] *= Ry2meV
    inv_tau_freq_tmp[:, :, :, 1] *= Ry2meV

    # compute placeholder value in meV and build valid mask
    placeholder_meV = placeholder * Ry2meV
    # any entry greater than placeholder_meV/2 is considered valid (placeholder_meV is very negative)
    valid_mask = inv_tau_freq_tmp[:, :, :, 0] > (placeholder_meV / 2.0)

    # ensure we have at least one valid data point before taking min/max
    valid_energies = inv_tau_freq_tmp[:, :, :, 0][valid_mask]
    if valid_energies.size == 0:
        print(f"{label}: ERROR - no valid energy entries found in {file1}. Exiting process_data.")
        # return empty/zero arrays consistent with the rest of the pipeline
        freqs = np.arange(mob_nfreq) * step
        return freqs, np.zeros(mob_nfreq), np.zeros(mob_nfreq), 0.0, 0.0, 0.0, mob_maxfreq

    if 'Electron' in label:
        cbm = np.min(valid_energies)
        print(f'{label}: CBM located at {cbm:.3f} meV')
        # shift energies so CBM -> 0
        inv_tau_freq_tmp[:, :, :, 0] -= cbm
    else:
        vbm = np.max(valid_energies)
        print(f'{label}: VBM located at {vbm:.3f} meV')
        # shift energies so VBM -> 0
        inv_tau_freq_tmp[:, :, :, 0] -= vbm

    energy = 1000 * kb * T * 3.0 / 2
    print(f'{label}: Energy at 3/2 k_BT = {energy:.3f} meV')

    inv_tau_freq = np.zeros(mob_nfreq)
    ndegen = np.zeros(mob_nfreq)

    # choose condition depending on carrier type
    if 'Electron' in label:
        # energies are measured upwards from CBM; match E ~= +energy
        def match_E(E):
            return abs(E - energy) < tol
    else:
        # energies are measured downwards from VBM (after shift they are negative); match E ~= -energy
        def match_E(E):
            return abs(E + energy) < tol

    for ikpt in range(nkpt):
        for ibnd in range(nbnd):
            for ifreq in range(mob_nfreq):
                # skip invalid entries quickly
                if not valid_mask[ikpt, ibnd, ifreq]:
                    continue
                E = inv_tau_freq_tmp[ikpt, ibnd, ifreq, 0]
                if match_E(E):
                    inv_tau_freq[ifreq] += inv_tau_freq_tmp[ikpt, ibnd, ifreq, 1]
                    ndegen[ifreq] += 1

    mask = ndegen > 0
    inv_tau_freq[mask] /= ndegen[mask]

    integral = np.sum(inv_tau_freq) * step
    integral_THz = integral * fact
    print(f'{label}: Integrated τ⁻¹ = {integral_THz:.3f} THz')

    def gauss(x, a, b, c):
        prefactor = a * step / (c * np.sqrt(2 * np.pi))
        return prefactor * np.exp(-0.5 * ((x - b) / c)**2)

    tau_gauss = np.zeros(mob_nfreq)
    for ii in range(mob_nfreq):
        for jj in range(mob_nfreq):
            tau_gauss[ii] += gauss(ii * step, inv_tau_freq[jj] * fact, jj * step, c)

    cum = np.cumsum(tau_gauss)
    if np.max(cum) == 0:
        cum_percent = np.zeros_like(cum)
    else:
        cum_percent = 100 * cum / np.max(cum)

    idx_peak = np.argmax(tau_gauss)
    peak_freq = np.arange(mob_nfreq)[idx_peak] * step
    peak_cum_percent = float(cum_percent[idx_peak]) if cum_percent.size > 0 else 0.0

    print(f'{label}: Peak at {peak_freq:.2f} meV → cumulative {peak_cum_percent:.2f}%')

    freqs = np.arange(mob_nfreq) * step
    return freqs, tau_gauss, cum, integral_THz, peak_freq, peak_cum_percent, mob_maxfreq

###################################################################################
# --- Process both datasets ---
###################################################################################
e_data = process_data(
    'e-Mobility_vs_T/e-Mobility_quadrupole/inv_taucb_freq_0.fmt', 'Electron')
h_data = process_data(
    'h-Mobility_vs_T/h-Mobility_quadrupole/inv_tau_freq_0.fmt', 'Hole')

###################################################################################
# --- Combined Plot (Electron top, Hole bottom) ---
###################################################################################
fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

# ---- Electron ----
freqs, tau_gauss, cum, integral_THz, peak_freq, peak_cum_percent, mob_maxfreq = e_data
ax1 = axes[0]
ax1.fill_between(freqs, tau_gauss, color='hotpink', alpha=0.25)
ax1.plot(freqs, tau_gauss, color='hotpink', lw=2.3, label=r'$\frac{d\tau^{-1}}{d\omega}$')
ax1.set_ylabel(r'$\frac{d\tau^{-1}_{3/2 k_{\mathrm{B}} T}}{d\omega}$', color='hotpink', fontsize=28)
ax1.tick_params(axis='y', labelcolor='hotpink', length=6, width=1)
ax1.tick_params(axis='x', length=6, width=1)
ax1.set_xlim(0, mob_maxfreq)
ax1.set_ylim(-0.03 * np.max(tau_gauss), 1.05 * np.max(tau_gauss))

# Add "CB" at bottom-right corner (inside frame)
ax1.text(mob_maxfreq - 2, 0.025 * np.max(tau_gauss), 'CB',
         color='hotpink', fontsize=22, ha='right', va='bottom', fontweight='bold')

ax2 = ax1.twinx()
ax2.plot(freqs, cum, color='navy', lw=2.0, linestyle='--',
         label=r'$\int \tau^{-1}(\omega)\,d\omega$')
ax2.set_ylabel(r'Cumulative $\tau^{-1}_{3/2 k_{\mathrm{B}} T}$ [THz]',
               color='navy', fontsize=24)
ax2.tick_params(axis='y', labelcolor='navy', length=6, width=1)
ax2.set_ylim(-0.02 * np.max(cum), 1.02 * np.max(cum))

ax1.axvline(peak_freq, color='gray', linestyle=':', lw=1.3)
ax1.text(peak_freq + 2, np.max(tau_gauss) * 0.82,
         f'{peak_cum_percent:.1f}%', color='gray', fontsize=16)

legend_labels = [
    r'$\frac{d\tau^{-1}}{d\omega}$ [THz·meV$^{-1}$]',
    rf'$\int \tau^{{-1}}(\omega)\,d\omega = {integral_THz:.3f}$ THz',
    rf'Cumulative = {peak_cum_percent:.1f}% at {peak_freq:.2f} meV'
]
ax1.legend(legend_labels, loc='upper right', fontsize=20, frameon=False)

# ---- Hole ----
freqs, tau_gauss, cum, integral_THz, peak_freq, peak_cum_percent, mob_maxfreq = h_data
ax3 = axes[1]
ax3.fill_between(freqs, tau_gauss, color='mediumseagreen', alpha=0.25)
ax3.plot(freqs, tau_gauss, color='mediumseagreen', lw=2.3,
         label=r'$\frac{d\tau^{-1}}{d\omega}$')
ax3.set_xlabel(r'Phonon frequency (meV)', fontsize=22)
ax3.set_ylabel(r'$\frac{d\tau^{-1}_{3/2 k_{\mathrm{B}} T}}{d\omega}$',
               color='mediumseagreen', fontsize=28)
ax3.tick_params(axis='y', labelcolor='mediumseagreen', length=6, width=1)
ax3.set_xlim(0, mob_maxfreq)
ax3.set_ylim(-0.03 * np.max(tau_gauss), 1.05 * np.max(tau_gauss))

# Add "VB" at bottom-right corner (inside frame)
ax3.text(mob_maxfreq - 2, 0.025 * np.max(tau_gauss), 'VB',
         color='mediumseagreen', fontsize=22, ha='right', va='bottom', fontweight='bold')

ax4 = ax3.twinx()
ax4.plot(freqs, cum, color='navy', lw=2.0, linestyle='--',
         label=r'$\int \tau^{-1}(\omega)\,d\omega$')
ax4.set_ylabel(r'Cumulative $\tau^{-1}_{3/2 k_{\mathrm{B}} T}$ [THz]',
               color='navy', fontsize=24)
ax4.tick_params(axis='y', labelcolor='navy', length=6, width=1)
ax4.set_ylim(-0.02 * np.max(cum), 1.02 * np.max(cum))

ax3.axvline(peak_freq, color='gray', linestyle=':', lw=1.3)
ax3.text(peak_freq + 2, np.max(tau_gauss) * 0.82,
         f'{peak_cum_percent:.1f}%', color='gray', fontsize=16)

legend_labels = [
    r'$\frac{d\tau^{-1}}{d\omega}$ [THz·meV$^{-1}$]',
    rf'$\int \tau^{{-1}}(\omega)\,d\omega = {integral_THz:.3f}$ THz',
    rf'Cumulative = {peak_cum_percent:.1f}% at {peak_freq:.2f} meV'
]
ax3.legend(legend_labels, loc='upper right', fontsize=20, frameon=False)

plt.tight_layout(pad=0.7)
save_dir = os.path.abspath('.')
outpath = os.path.join(save_dir, 'inv_tau_spectral_combined.png')
plt.savefig(outpath, dpi=300, bbox_inches='tight')
plt.close(fig)

print(f"\n✅ Combined plot saved at:\n   {outpath}\n")
