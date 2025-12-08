#!/usr/bin/env python3
# Plot |g_{mn,k}| from EPW along high-symmetry path with NSCF → EPW band mapping

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import elphmod

# -----------------------
# Read environment variables
# -----------------------
try:
    prefix = os.environ['prefix']
    n_band = int(os.environ['n_band'])  # NSCF band number (initial)
    m_band = int(os.environ['m_band'])  # NSCF band number (final)
except KeyError as e:
    raise KeyError(f"Environment variable {e} not found. Please export prefix, n_band, and m_band.")

# Optional: provide NSCF → EPW mapping via bash variable
exclude_bands_str = os.environ.get('exclude_bands', '')
exclude_set = set()
for part in exclude_bands_str.split(','):
    if ':' in part:
        start, end = map(int, part.split(':'))
        exclude_set.update(range(start, end+1))
    elif part.strip() != '':
        exclude_set.add(int(part))

# Define total NSCF bands (adjust as needed)
all_nscf_bands = range(1, 101)
epw_map = {}
epw_idx = 0
for b in all_nscf_bands:
    if b not in exclude_set:
        epw_map[b] = epw_idx
        epw_idx += 1

if n_band not in epw_map or m_band not in epw_map:
    raise ValueError(f"NSCF band {n_band} or {m_band} not in EPW calculation (excluded?)")

n_idx = epw_map[n_band]
m_idx = epw_map[m_band]

# -----------------------
# User settings
# -----------------------
material = prefix
path = 'KGM'
Npath = 80

# filenames exactly as in File 1
epw_L_files = {
    'no_lr': 'epw_no_lr.out',
    'quadrupole': 'epw_quadrupole.out'
}
epmatwp_files = {
    'no_lr': 'no_lr.epmatwp',
    'quadrupole': 'quadrupole.epmatwp'
}
wigner_fmt = 'wigner.fmt'

plt.rcParams.update({
    'font.size': 22,
    'axes.linewidth': 0.8,
    'font.family': 'serif',
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
})

comm = elphmod.MPI.comm

# -----------------------
# Helper checks
# -----------------------
for fname in list(epw_L_files.values()) + list(epmatwp_files.values()) + [wigner_fmt]:
    if not os.path.exists(fname):
        print(f"Warning: file '{fname}' not found. Make sure it exists in the working directory.", file=sys.stderr)

# -----------------------
# High-symmetry path
# -----------------------
q, x, corners = elphmod.bravais.path(path, ibrav=4, N=Npath, moveG=0.005)
el = elphmod.el.Model(material)

plot_models = ['no_lr', 'quadrupole']
g_results = {}

fig, ax = plt.subplots(figsize=(14,8))
EPS_W = 1e-8

for lr in plot_models:
    Lfile = epw_L_files.get(lr)
    if Lfile is None or not os.path.exists(Lfile):
        print(f"epw L file for model '{lr}' not found: {Lfile}. Skipping this model.", file=sys.stderr)
        continue

    quadrupole_fmt = '_quadrupole.fmt' if lr == 'quadrupole' else None
    ph = elphmod.ph.Model('dyn',
                          apply_asr_simple=True,
                          apply_zasr=True,
                          lr=(lr != 'no_lr'),
                          lr2d=True,
                          L=elphmod.elph.read_L(Lfile),
                          quadrupole_fmt=quadrupole_fmt)

    efile = epmatwp_files.get(lr)
    try:
        elph = elphmod.elph.Model(efile, wigner_fmt, el, ph)
    except Exception as e:
        print(f"Failed to load elph model for '{lr}': {e}", file=sys.stderr)
        continue

    g_list = []
    for q1,q2,q3 in q:
        try:
            gq = elph.g(q1,q2,q3,elbnd=True,phbnd=True)
        except Exception as e:
            print(f"elph.g failed at q=({q1},{q2},{q3}): {e}", file=sys.stderr)
            gq = np.zeros((ph.size, 1, 1))
        g_list.append(np.abs(gq))

    g_arr = np.array(g_list)
    n_q = g_arr.shape[0]
    n_ph = ph.size
    w = elphmod.ph.sgnsqrt(elphmod.dispersion.dispersion(ph.D,q))

    g_meV = np.zeros((n_q, n_ph))
    for nu in range(n_ph):
        try:
            g_nu = g_arr[:, nu, n_idx, m_idx]
        except Exception:
            g_nu = np.reshape(g_arr, (n_q, -1))[:, min(nu, g_arr.size // n_q - 1)]
        w_nu = w[:, nu]
        w_nu_safe = np.where(w_nu <= 0, EPS_W, w_nu)
        g_meV[:, nu] = g_nu * np.sqrt(elphmod.misc.Ry * 1e3 / w_nu_safe)

    g_results[lr] = {'g_meV': g_meV, 'n_ph': n_ph, 'x': x}

# -----------------------
# Final plotting
# -----------------------
if 'no_lr' in g_results and 'quadrupole' in g_results:
    n_ph = g_results['quadrupole']['n_ph']
    if n_ph <= 10:
        cmap = plt.cm.get_cmap('tab10', n_ph)
    elif n_ph <= 20:
        cmap = plt.cm.get_cmap('tab20', n_ph)
    else:
        cmap = plt.cm.get_cmap('tab20', 20)
    mode_colors = [cmap(i % cmap.N) for i in range(n_ph)]

    x = g_results['quadrupole']['x']
    g_no_lr = g_results['no_lr']['g_meV']
    g_quad = g_results['quadrupole']['g_meV']

    for nu in range(n_ph):
        color = mode_colors[nu]

        # no_lr = triangle, dashed
        ax.plot(x, g_no_lr[:, nu],
                color=color, linestyle='--',
                marker='^', markersize=6, fillstyle='full',
                linewidth=1.8)

        # quadrupole = open circle, solid
        ax.plot(x, g_quad[:, nu],
                color=color, linestyle='-',
                marker='o', markersize=6, fillstyle='none',
                linewidth=2.0)

    # Two legends: models + modes
    model_handles = [
        plt.Line2D([], [], color='k', marker='^', linestyle='--', markersize=6, label='no long-range'),
        plt.Line2D([], [], color='k', marker='o', linestyle='-', markerfacecolor='none', markersize=6,
                   label='Dipole+quadrupole')
    ]
    mode_handles = [plt.Line2D([], [], color=mode_colors[i], marker='s', linestyle='None', label=f'Mode {i+1}')
                    for i in range(n_ph)]
    legend1 = ax.legend(handles=model_handles, loc='upper right', title='Models', fontsize=18)
    legend2 = ax.legend(handles=mode_handles, loc='upper left', title='Phonon Modes ($\lambda$)', fontsize=18, ncol=2)
    ax.add_artist(legend1)

else:
    print("Missing one of the model datasets; skipping combined plot.", file=sys.stderr)

# -----------------------
# Figure labels, axes, save
# -----------------------
highsym_labels = [r'K', r'$\Gamma$', 'M']
highsym_positions = x[corners]
ax.set_xticks(highsym_positions)
ax.set_xticklabels(highsym_labels)
ax.set_xlim(x[0], x[-1])
#ax.set_yscale('log')
ax.set_xlabel('Phonon q along BZ')
ax.set_ylabel(rf'$|g_{{{n_band},{m_band}}}^{{\lambda}}(\mathbf{{k}}$=VBM,$\mathbf{{q}})|$ (meV)')
#ax.set_title(f'Electron–Phonon Coupling in {material}')

for pos in highsym_positions:
    ax.axvline(x=pos, color='k', linestyle='-', linewidth=0.6)

plt.tight_layout()
outfile = f'{material}_g_{n_band}_{m_band}_noLR_vs_quadrupole_meV.png'
plt.savefig(outfile, dpi=300)
print("Saved:", outfile)

