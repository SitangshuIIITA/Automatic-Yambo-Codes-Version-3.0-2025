#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import elphmod
from PIL import Image

# -----------------------------
# Plot settings
# -----------------------------
csfont = {'fontname':'Arial Narrow'}
plt.rcParams['font.weight'] = 'ultralight'
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['ytick.major.size'] = 8
plt.rcParams.update({'font.size': 50})

plot_coupling_in_mode_basis = True  # rotate g to phonon mode basis

# -----------------------------
# Load saved data
# -----------------------------
x = np.load('x.npy')
Dd = np.load('Dd.npy')       # Z*+Q fit
gd = np.load('gd.npy')       # Z*+Q fit
Dq = np.load('Dq.npy')       # Z* only
gq = np.load('gq.npy')       # Z* only
D0 = np.load('D0.npy')       # DFPT
g0 = np.load('g0.npy')       # DFPT

# Load DFPT reference q-path
q0, x0, w0_ref = elphmod.el.read_bands('dynref.freq')
x0 += x[-1] - x0[-1]  # scale to match x-path
w0 = elphmod.ph.sgnsqrt(w0_ref)

# Remove x-axis offset
x -= x[0]

# -----------------------------
# Diagonalize dynamical matrices
# -----------------------------
wd2, ud = np.linalg.eigh(Dd)
wq2, uq = np.linalg.eigh(Dq)
w02, u0 = np.linalg.eigh(D0)
wd = elphmod.ph.sgnsqrt(wd2)
wq = elphmod.ph.sgnsqrt(wq2)
w0 = elphmod.ph.sgnsqrt(w02)

if plot_coupling_in_mode_basis:
    gd = np.einsum('qx,qxv->qv', gd, ud)
    gq = np.einsum('qx,qxv->qv', gq, uq)
    g0 = np.einsum('qx,qxv->qv', g0, u0)

# -----------------------------
# Number of phonon modes and colors
# -----------------------------
num_modes = Dd.shape[1]
colors = plt.cm.nipy_spectral(np.linspace(0,1,num_modes))

# -----------------------------
# Plot phonon energies
# -----------------------------
plt.figure(figsize=(30,10))

# Plot first mode with legend
plt.plot(x, wd[:,0]*1e3*elphmod.misc.Ry, color=colors[0], lw=3, label='Z*+Q')
plt.plot(x, wq[:,0]*1e3*elphmod.misc.Ry, color=colors[0], lw=3, linestyle='--', label='Z* only')
plt.scatter(x0, w0[:,0]*1e3*elphmod.misc.Ry, color=colors[0], marker='o', s=150, edgecolors='k', zorder=3, label='DFPT Mode 1')

# Remaining modes (no extra legend for Fit/Z*+Q or Fit/Z* only)
for nu in range(1, num_modes):
    plt.plot(x, wd[:,nu]*1e3*elphmod.misc.Ry, color=colors[nu], lw=3)       # Fit Z*+Q
    plt.plot(x, wq[:,nu]*1e3*elphmod.misc.Ry, color=colors[nu], lw=3, linestyle='--')  # Fit Z* only
    plt.scatter(x0, w0[:,nu]*1e3*elphmod.misc.Ry, color=colors[nu], marker='o', s=150, edgecolors='k', zorder=3, label=f'DFPT Mode {nu+1}')

# Axes labels and title
plt.xlabel('Phonon q along BZ', fontsize=50)
plt.ylabel('Phonon energy (meV)', fontsize=50)
#plt.title('Phonon energies: DFPT vs Fit', fontsize=35)

# Axis formatting
ax = plt.gca()
x_offset = 0.005 * (x[-1] - x[0])  # 5% of total x-range
ax.set_xlim(x[0] - x_offset, x[-1] + x_offset)
ax.set_xticks([x[0], x[-1]])
ax.set_xticklabels([r'$\Gamma$', 'M'])
for spine in ax.spines.values():
    spine.set_linewidth(2)

# Legend
plt.legend(fontsize=25, loc='upper right', frameon=False)

# Ticks
plt.tick_params(axis="y", direction='in', length=12, width=4)
plt.tick_params(axis="x", direction='in', length=12, width=4)

# Save
plt.tight_layout()
plt.savefig('fitQa.png', dpi=300, bbox_inches='tight')
plt.close()

# -----------------------------
# Plot electron-phonon couplings
# -----------------------------
plt.figure(figsize=(30,10))

# Plot first mode with legend for Z*+Q and Z* only
plt.plot(x, np.abs(gd[:,0])*elphmod.misc.Ry*1e3, color=colors[0], lw=3, label='Z*+Q')
plt.plot(x, np.abs(gq[:,0])*elphmod.misc.Ry*1e3, color=colors[0], lw=3, linestyle='--', label='Z* only')
plt.scatter(x0, np.abs(g0[:,0])*elphmod.misc.Ry*1e3, color=colors[0], s=150, edgecolors='k', zorder=3, label='DFPT Mode 1')

# Remaining modes (no extra legend for Z*+Q / Z* only)
for nu in range(1, num_modes):
    plt.plot(x, np.abs(gd[:,nu])*elphmod.misc.Ry*1e3, color=colors[nu], lw=3)       # Z*+Q
    plt.plot(x, np.abs(gq[:,nu])*elphmod.misc.Ry*1e3, color=colors[nu], lw=3, linestyle='--')  # Z* only
    plt.scatter(x0, np.abs(g0[:,nu])*elphmod.misc.Ry*1e3, color=colors[nu], s=150, edgecolors='k', zorder=3, label=f'DFPT Mode {nu+1}')

# Axes labels and title
plt.xlabel('Phonon q along BZ', fontsize=50)
plt.ylabel(r'$|g_{n,m}^{\lambda}(\mathbf{k}=\Gamma,\mathbf{q})|$ (meV)', fontsize=50)

#ax.set_ylabel(rf'$|g_{{{n_band},{m_band},\mathrm{{k}}=\Gamma}}^{{\lambda}}|$ (meV)')
#plt.title('Electron-Phonon Coupling: DFPT vs Quadrupole Fit', fontsize=35)

# Axis formatting
ax = plt.gca()
x_offset = 0.005 * (x[-1] - x[0])  # 5% of total x-range
ax.set_xlim(x[0] - x_offset, x[-1] + x_offset)
#
ax.set_xticks([x[0], x[-1]])
ax.set_xticklabels([r'$\Gamma$', 'M'])
for spine in ax.spines.values():
    spine.set_linewidth(2)

# Legend
plt.legend(fontsize=25, loc='upper right', frameon=False)

# Ticks
plt.tick_params(axis="y", direction='in', length=12, width=4)
plt.tick_params(axis="x", direction='in', length=12, width=4)

# Save
plt.tight_layout()
plt.savefig('fitQb.png', dpi=300, bbox_inches='tight')
plt.close()


# -----------------------------
# Merge figures top and bottom
# -----------------------------
img1 = Image.open('fitQa.png')
img2 = Image.open('fitQb.png')

# Width is the maximum of the two images, height is sum
width = max(img1.width, img2.width)
height = img1.height + img2.height

# Create a new blank image
new_img = Image.new('RGB', (width, height), (255, 255, 255))

# Paste images
new_img.paste(img1, (0, 0))                    # Top
new_img.paste(img2, (0, img1.height))          # Bottom

# Save the merged image
new_img.save('fitQ.png')












