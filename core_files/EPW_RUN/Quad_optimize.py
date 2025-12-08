#!/usr/bin/env python3
import elphmod
import os
import numpy as np
import scipy.optimize
from mpi4py import MPI

comm = elphmod.MPI.comm
info = elphmod.MPI.info

# -----------------------------
# Settings
# -----------------------------
plot_coupling_in_mode_basis = True  # rotate g to phonon mode basis
info('Load tight-binding, mass-spring, and coupling models')

# --- Models ---

prefix = os.environ['prefix']  # will raise KeyError if not set

# --- Models ---
el = elphmod.el.Model(prefix)
ph = elphmod.ph.Model('dyn', lr=False)
elph = elphmod.elph.Model(f'work/{prefix}.epmatwp', f'work/{prefix}.wigner', el, ph)
elph.sample_orig()

# -----------------------------
# q-path (Γ–M)
# -----------------------------
path_points = [(0.0, 0.0, 0.0), (0.0, 0.5, 0.0)]
q, x, special = elphmod.bravais.path(path_points, ibrav=4, N=500)
#q, x, special = elphmod.bravais.path([(0, 0, 0), (0.5, 0, 0)], ibrav=4, N=500)

# -----------------------------
# DFPT reference data
# -----------------------------
q0, x0, w0 = elphmod.el.read_bands('dynref.freq')
q0 = 2 * np.pi * np.dot(q0, ph.a.T) / np.linalg.norm(ph.a[0])
x0 += x[-1] - x0[-1]
sqrtM = np.sqrt(np.repeat(ph.M, 3))

# Dynamical matrices
D0 = np.empty((len(q0), ph.size, ph.size), dtype=complex)
for iq in range(len(q0)):
    D0[iq] = elphmod.ph.read_flfrc(f'dynref{iq+1}')[0][1][0]
D0 /= sqrtM[np.newaxis, np.newaxis, :]
D0 /= sqrtM[np.newaxis, :, np.newaxis]

# Electron-phonon couplings
g0 = np.zeros((len(q0), ph.size), dtype=complex)
iq = 0
with open('phref.out') as fh:
    lines = iter(fh)
    for line in lines:
        if 'Printing the electron-phonon matrix elements' in line:
            read_count = 0
            while read_count < ph.size:
                line_data = next(lines).strip()
                if not line_data:
                    continue
                parts = line_data.split()
                if len(parts) < 2:
                    continue
                re, im = parts[:2]
                g0[iq, read_count] = float(re) + 1j * float(im)
                read_count += 1
            iq += 1
g0 /= sqrtM[np.newaxis, :]

# -----------------------------
# Optimize long-range separation
# -----------------------------
ph.lr = True
def objective_L(L):
    ph.L, = L
    ph.update_short_range()
    return ph.sum_force_constants()

scipy.optimize.minimize(objective_L, [1.0], tol=0.01)
elph.update_short_range()

# -----------------------------
# Sampling helper
# -----------------------------
def sample(qpath):
    D = elphmod.dispersion.sample(ph.D, qpath)
    g = elphmod.dispersion.sample(elph.g, qpath, elbnd=True, comm=elphmod.MPI.I)[:, :, 0, 0]
    return D, g

Dd, gd = sample(q)

# -----------------------------
# Quadrupole fitting
# -----------------------------
def error():
    D, g = sample(q0)
    g *= elphmod.misc.Ry ** 1.5
    dD = (abs(D - D0) ** 2).sum()
    dg = (abs(abs(g) - abs(g0)) ** 2).sum()
    return dD, dg

dD0, dg0 = error()

def objective_Q(Q, dD0=dD0, dg0=dg0):
    # --- General quadrupole assignment for ph.nat <= 4 ---
    ph.Q = np.zeros((ph.nat, 3, 3, 3))

    # Mapping of independent Q components
    independent_Q_map = {}
    if ph.nat == 2:
        independent_Q_map = {
            0: [(1,1,1,0)],
            1: [(1,1,1,1), (2,1,1,2)]
        }

    elif ph.nat == 3:
        independent_Q_map = {
            1: [(1,1,1,0)],
            2: [(1,1,1,1), (2,1,1,2)]
        }
    elif ph.nat == 4:
        independent_Q_map = {
            0: [(1,1,1,0)],
            1: [(1,1,1,1)],
            2: [(2,1,1,2)],
            3: [(2,2,1,3)]
        }
    else:
        raise NotImplementedError(f"ph.nat={ph.nat} not yet mapped for Q assignment.")

    # Assign independent components from Q vector
    for atom_idx, comps in independent_Q_map.items():
        for (i,j,k,q_idx) in comps:
            ph.Q[atom_idx,i,j,k] = Q[q_idx]

    # Symmetrization
    ph.Q[:,0,0,1] = ph.Q[:,0,1,0] = ph.Q[:,1,0,0] = -ph.Q[:,1,1,1]
    ph.Q[:,2,0,0] = ph.Q[:,2,1,1]

    if ph.nat == 3:
        ph.Q[0,:2] = ph.Q[2,:2]
        ph.Q[0,2] = -ph.Q[2,2]
    elif ph.nat == 4:
        ph.Q[2,:2] = ph.Q[0,:2]
        ph.Q[3,:2] = ph.Q[1,:2]
        ph.Q[2,2] = -ph.Q[0,2]
        ph.Q[3,2] = -ph.Q[1,2]

    ph.Q[np.abs(ph.Q) < 1e-12] = 0.0

    # Update models
    ph.update_short_range()
    elph.update_short_range()

    # Compute error
    dD, dg = error()
    dD /= dD0
    dg /= dg0
    dg_scaled = dg * 0.5
    info(f'error(D) = {100*dD:.3f}% | error(g) = {100*dg:.3f}%')
    return np.sqrt(dD**2 + dg_scaled**2)


res = scipy.optimize.minimize(objective_Q, np.ones(3), method='Nelder-Mead', tol=1e-5)
info(f'Quadrupole fit done. Optimal Q = {res.x}')

# -----------------------------
# Save fitted results
# -----------------------------
Dq, gq = sample(q)

if comm.rank != 0:
    raise SystemExit

elphmod.ph.write_quadrupole_fmt('quadrupole.fmt', ph.Q)

# --- Diagonalize dynamical matrices for phonon energies ---
wd2, ud = np.linalg.eigh(Dd)   # sample Q=0 (Z* only)
wq2, uq = np.linalg.eigh(Dq)   # optimized quadrupole
w02, u0 = np.linalg.eigh(D0)   # DFPT reference 

# Save data
np.save('x.npy', x)
np.save('Dd.npy', Dd)
np.save('gd.npy', gd)
np.save('Dq.npy', Dq)
np.save('gq.npy', gq)
np.save('D0.npy', D0)
np.save('g0.npy', g0)

info('Optimization complete. Data saved for plotting in plot_fitQ.py')

