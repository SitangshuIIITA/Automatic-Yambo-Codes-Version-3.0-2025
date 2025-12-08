#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D

# ============================================================
# Plot style
# ============================================================
plt.rcParams.update({
    'font.size': 26,
    'axes.linewidth': 1.5,
    'font.family': 'serif',
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'lines.markersize': 10,
    'legend.handlelength': 3.6,
    'figure.dpi': 300
})

# ============================================================
# Input files
# ============================================================
elec_file = "mobility_summary_electron.txt"
hole_file = "mobility_summary_hole.txt"

# ============================================================
# Read data
# ============================================================
def read_data(filename):
    data = np.loadtxt(filename, comments="#", skiprows=2)
    T = data[:, 0]
    muS_d, muB_d = data[:, 1], data[:, 2]
    rH_S_d, rH_B_d = data[:, 3], data[:, 4]
    muS_q, muB_q = data[:, 5], data[:, 6]
    rH_S_q, rH_B_q = data[:, 7], data[:, 8]
    return T, muS_d, muB_d, muS_q, muB_q, rH_S_d, rH_B_d, rH_S_q, rH_B_q

T_e, muS_d_e, muB_d_e, muS_q_e, muB_q_e, rH_S_d_e, rH_B_d_e, rH_S_q_e, rH_B_q_e = read_data(elec_file)
T_h, muS_d_h, muB_d_h, muS_q_h, muB_q_h, rH_S_d_h, rH_B_d_h, rH_S_q_h, rH_B_q_h = read_data(hole_file)

# ============================================================
# Power-law fitting
# ============================================================
def power_law(T, A, n):
    return A * T**n

def fit_exponent(T, mu):
    popt, _ = curve_fit(power_law, T, mu, p0=[1e3, -1])
    return popt  # A, n

fits = {}
for label, (T, mu) in {
    "e_S_d": (T_e, muS_d_e), "e_B_d": (T_e, muB_d_e),
    "e_S_q": (T_e, muS_q_e), "e_B_q": (T_e, muB_q_e),
    "h_S_d": (T_h, muS_d_h), "h_B_d": (T_h, muB_d_h),
    "h_S_q": (T_h, muS_q_h), "h_B_q": (T_h, muB_q_h)
}.items():
    fits[label] = fit_exponent(T, mu)

# ============================================================
# Color scheme, markers, and line styles
# ============================================================
colors = {
    "electron_dipole": "orange",
    "electron_quad": "green",
    "hole_dipole": "hotpink",
    "hole_quad": "slategray"
}
markers = {"dipole": "s", "quad": "o"}
linestyles = {"SERTA": "-", "BTE": "--"}

# ============================================================
# Legend setup
# ============================================================
legend_lines = [
    Line2D([0], [0], color='k', lw=2, ls='-', label='SERTA'),
    Line2D([0], [0], color='k', lw=2, ls='--', label='BTE'),
    Line2D([0], [0], color='k', marker='s', lw=0, label='Dipole'),
    Line2D([0], [0], color='k', marker='o', lw=0, label='Dipole + Quadrupole')
]

# ============================================================
# Annotate Tⁿ fits automatically (non-overlapping)
# ============================================================
def annotate_fits(ax, T, mu_dict, fits_dict, colors, prefix):
    sort_idx = np.argsort(T)
    T_sorted = T[sort_idx]

    T_first = T_sorted[0]
    x_base_shift = 0.17 * T_first  # baseline 5% to the right
    T_range = T_sorted[-1] - T_sorted[0] if len(T_sorted) > 1 else T_first
    close_threshold = 0.001 * T_range  # 10% of total T range

    y_offset_up = 0.95
    y_offset_down = 0.81
    x_extra_shift = 0.09 * T_first   # ±2% horizontal nudge

    suffixes = [
        ("S_d", "dipole", "SERTA"),
        ("B_d", "dipole", "BTE"),
        ("S_q", "quad", "SERTA"),
        ("B_q", "quad", "BTE")
    ]

    prev_positions = []  # to store prior annotation positions

    for i, (suffix, label, markertype) in enumerate(suffixes):
        key = f"{prefix}_{suffix}"
        mu = mu_dict[suffix][sort_idx]
        n_val = fits_dict[key][1]
        color = colors[f"{'electron' if prefix == 'e' else 'hole'}_{label}"]

        # Base annotation location
        T_annot = T_first + x_base_shift
        mu_base = mu[0]

        # --- Avoid overlaps ---
        # Default: place above
        y_factor = y_offset_up
        # Add small horizontal staggering pattern (+/-/+/−)
        T_annot += ((-1)**i) * x_extra_shift

        # Check if any previous annotation is too close
        for prev_T, prev_mu in prev_positions:
            if abs(T_annot - prev_T) < close_threshold:
                # If too close horizontally → flip vertical placement
                y_factor = y_offset_down
                T_annot += x_extra_shift * 1.2  # nudge further right

        mu_annot = mu_base * y_factor

        # Keep within y-limits
        ymax = ax.get_ylim()[1]
        if mu_annot > ymax * 0.97:
            mu_annot = ymax * 0.95

        ax.text(
            T_annot,
            mu_annot,
            f"T$^{{{n_val:.2f}}}$",
            color=color,
            fontsize=26,
            ha='left',
            va='bottom'
        )

        prev_positions.append((T_annot, mu_annot))



# ============================================================
# 1️⃣ Mobility (top: electron, bottom: hole)
# ============================================================
fig, axes = plt.subplots(2, 1, figsize=(15, 12), sharex=True)

# --- Electrons
ax = axes[0]
ax.plot(T_e, muS_d_e, marker=markers["dipole"], ls=linestyles["SERTA"], color=colors["electron_dipole"], lw=1.8)
ax.plot(T_e, muB_d_e, marker=markers["dipole"], ls=linestyles["BTE"], color=colors["electron_dipole"], lw=1.8)
ax.plot(T_e, muS_q_e, marker=markers["quad"], ls=linestyles["SERTA"], color=colors["electron_quad"], lw=2.2)
ax.plot(T_e, muB_q_e, marker=markers["quad"], ls=linestyles["BTE"], color=colors["electron_quad"], lw=2.2)
ax.set_ylabel(r"Electron Mobility (cm$^2$/V·s)")

annotate_fits(ax, T_e, {
    "S_d": muS_d_e, "B_d": muB_d_e, "S_q": muS_q_e, "B_q": muB_q_e
}, fits, colors, "e")

# --- Holes
ax = axes[1]
ax.plot(T_h, muS_d_h, marker=markers["dipole"], ls=linestyles["SERTA"], color=colors["hole_dipole"], lw=1.8)
ax.plot(T_h, muB_d_h, marker=markers["dipole"], ls=linestyles["BTE"], color=colors["hole_dipole"], lw=1.8)
ax.plot(T_h, muS_q_h, marker=markers["quad"], ls=linestyles["SERTA"], color=colors["hole_quad"], lw=2.2)
ax.plot(T_h, muB_q_h, marker=markers["quad"], ls=linestyles["BTE"], color=colors["hole_quad"], lw=2.2)
ax.set_xlabel(r"$\mathrm{Temperature\ (K)}$")
ax.set_ylabel(r"Hole Mobility (cm$^2$/V·s)")

annotate_fits(ax, T_h, {
    "S_d": muS_d_h, "B_d": muB_d_h, "S_q": muS_q_h, "B_q": muB_q_h
}, fits, colors, "h")

axes[0].legend(handles=legend_lines, fontsize=22, frameon=False,
               loc='upper right', bbox_to_anchor=(1.0, 1.0), ncol=1)

plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig("mobility_dipole_quadrupole.png", dpi=300)
plt.close(fig)

# ============================================================
# 2️⃣ Hall Factor (top: electron, bottom: hole)
# ============================================================
fig, axes = plt.subplots(2, 1, figsize=(15, 12), sharex=True)

# --- Electrons
ax = axes[0]
ax.plot(T_e, rH_S_d_e, marker=markers["dipole"], ls=linestyles["SERTA"], color=colors["electron_dipole"], lw=1.8)
ax.plot(T_e, rH_B_d_e, marker=markers["dipole"], ls=linestyles["BTE"], color=colors["electron_dipole"], lw=1.8)
ax.plot(T_e, rH_S_q_e, marker=markers["quad"], ls=linestyles["SERTA"], color=colors["electron_quad"], lw=2.2)
ax.plot(T_e, rH_B_q_e, marker=markers["quad"], ls=linestyles["BTE"], color=colors["electron_quad"], lw=2.2)
ax.set_ylabel(r"Hall Factor, $r_{xy}^e$")

# --- Holes
ax = axes[1]
ax.plot(T_h, rH_S_d_h, marker=markers["dipole"], ls=linestyles["SERTA"], color=colors["hole_dipole"], lw=1.8)
ax.plot(T_h, rH_B_d_h, marker=markers["dipole"], ls=linestyles["BTE"], color=colors["hole_dipole"], lw=1.8)
ax.plot(T_h, rH_S_q_h, marker=markers["quad"], ls=linestyles["SERTA"], color=colors["hole_quad"], lw=2.2)
ax.plot(T_h, rH_B_q_h, marker=markers["quad"], ls=linestyles["BTE"], color=colors["hole_quad"], lw=2.2)
ax.set_xlabel(r"$\mathrm{Temperature\ (K)}$")
ax.set_ylabel(r"Hall Factor, $r_{xy}^h$")

axes[0].legend(handles=legend_lines, fontsize=20, frameon=False, loc='lower right', ncol=1)

plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig("hallfactor_dipole_quadrupole.png", dpi=300)
plt.close(fig)

print("✅ Saved: mobility_dipole_quadrupole.png and hallfactor_dipole_quadrupole.png")

