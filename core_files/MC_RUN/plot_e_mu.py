#!/usr/bin/env python3
"""
MC Diagnostic Suite (Option B) with multi-page PDF output.

Reads MC_OUTPUT/mu_time.txt, vel_time.txt and log/log.txt,
generates all scientific diagnostic plots and exports a
single clean multi-page PDF: MC_OUTPUT/MC_Diagnostics.pdf

Author: Sitangshu Bhattacharya
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

OUT = "MC_OUTPUT"
LOG = "log"

# --- Load mobility and velocity time-series ---
mu = np.loadtxt(f"{OUT}/mu_time.txt")
vel = np.loadtxt(f"{OUT}/vel_time.txt")

t_s = mu[:,0]
t_ps = t_s * 1e12

mu_x = mu[:,1]
mu_y = mu[:,2]
mu_z = mu[:,3]

v_mean = vel[:,1]

# --- Load variance info from log ---
var_t = []
var_x = []
var_y = []
var_z = []

with open(f"{LOG}/log.txt") as f:
    for line in f:
        if "var_x=" in line:
            parts = line.strip().split()
            for p in parts:
                if "var_x" in p:
                    vx = float(p.split("=")[1].replace(",", ""))
                if "var_y" in p:
                    vy = float(p.split("=")[1].replace(",", ""))
                if "var_z" in p:
                    vz = float(p.split("=")[1].replace(",", ""))
            var_x.append(vx)
            var_y.append(vy)
            var_z.append(vz)
            var_t.append(len(var_t)*(t_s[1]-t_s[0]))

var_t = np.array(var_t)*1e12  # ps

# --- Helper ---
def sliding_R2(t, var, window=20):
    R2 = []
    for i in range(len(var)):
        if i < window:
            R2.append(np.nan)
            continue
        tt = t[i-window:i]
        vv = var[i-window:i]
        A = np.vstack([tt, np.ones_like(tt)]).T
        m, c = np.linalg.lstsq(A, vv, rcond=None)[0]
        pred = m * tt + c
        ss_res = np.sum((vv - pred)**2)
        ss_tot = np.sum((vv - vv.mean())**2)
        R2.append(1 - ss_res/(ss_tot+1e-30))
    return np.array(R2)

# --- PDF output ---
pdf_path = f"{OUT}/MC_Diagnostics.pdf"
pdf = PdfPages(pdf_path)

window = 20

# ==========================================================
# 1. Mobility vs Time
# ==========================================================
plt.figure(figsize=(7,5))
plt.plot(t_ps, mu_x, label="μₓ", lw=2)
plt.plot(t_ps, mu_y, label="μᵧ", lw=2)
plt.xlabel("Time (ps)")
plt.ylabel("μ (cm²/V·s)")
plt.title("Mobility vs Time")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT}/1_mu_vs_time.png", dpi=300)
pdf.savefig()
plt.close()

# ==========================================================
# 2. Mean Velocity vs Time
# ==========================================================
plt.figure(figsize=(7,5))
plt.plot(t_ps, v_mean, lw=2)
plt.xlabel("Time (ps)")
plt.ylabel("Mean Velocity |v| (m/s)")
plt.title("Mean Velocity vs Time")
plt.grid(True)
plt.tight_layout()
plt.savefig(f"{OUT}/2_velocity_vs_time.png", dpi=300)
pdf.savefig()
plt.close()

# ==========================================================
# 3. Stability Error vs Time
# ==========================================================
stab_x = []
stab_y = []

for i in range(len(mu_x)):
    if i < window:
        stab_x.append(np.nan)
        stab_y.append(np.nan)
        continue
    wx = mu_x[i-window:i]
    wy = mu_y[i-window:i]
    stab_x.append((wx.max()-wx.min())/(wx.mean()+1e-30))
    stab_y.append((wy.max()-wy.min())/(wy.mean()+1e-30))

plt.figure(figsize=(7,5))
plt.plot(t_ps, stab_x, label="rel_change μₓ")
plt.plot(t_ps, stab_y, label="rel_change μᵧ")
plt.axhline(0.01, color='r', ls='--', label="1% threshold")
plt.xlabel("Time (ps)")
plt.ylabel("Relative Stability Error")
plt.title("Stability Criterion vs Time")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT}/3_stability_error.png", dpi=300)
pdf.savefig()
plt.close()

# ==========================================================
# 4. Diffusive R² Quality
# ==========================================================
R2_x = sliding_R2(var_t, var_x)
R2_y = sliding_R2(var_t, var_y)

plt.figure(figsize=(7,5))
plt.plot(var_t, R2_x, label="R²ₓ")
plt.plot(var_t, R2_y, label="R²ᵧ")
plt.axhline(0.98, color='r', ls='--', label="0.98 threshold")
plt.xlabel("Time (ps)")
plt.ylabel("R² of Linear Diffusion")
plt.title("Diffusive R² vs Time")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT}/4_diffusive_R2.png", dpi=300)
pdf.savefig()
plt.close()

# ==========================================================
# 5. μ Slope Error vs Time
# ==========================================================
slope_x = []
slope_y = []

for i in range(len(mu_x)):
    if i < window:
        slope_x.append(np.nan)
        slope_y.append(np.nan)
        continue

    tt = t_ps[i-window:i]
    mx = mu_x[i-window:i]
    my = mu_y[i-window:i]

    A = np.vstack([tt, np.ones_like(tt)]).T
    m_x, c_x = np.linalg.lstsq(A, mx, rcond=None)[0]
    m_y, c_y = np.linalg.lstsq(A, my, rcond=None)[0]

    slope_x.append(abs(m_x*(tt[-1]-tt[0])))
    slope_y.append(abs(m_y*(tt[-1]-tt[0])))

plt.figure(figsize=(7,5))
plt.plot(t_ps, slope_x, label="Δμₓ")
plt.plot(t_ps, slope_y, label="Δμᵧ")
plt.axhline(0.5, color='r', ls='--', label="0.5 cm² threshold")
plt.xlabel("Time (ps)")
plt.ylabel("Δμ (cm²/V·s)")
plt.title("Slope Error (Δμ) vs Time")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT}/5_mu_slope_error.png", dpi=300)
pdf.savefig()
plt.close()

# ==========================================================
# 6. Variance Growth
# ==========================================================
plt.figure(figsize=(7,5))
plt.plot(var_t, var_x, label="σₓ²(t)")
plt.plot(var_t, var_y, label="σᵧ²(t)")
plt.xlabel("Time (ps)")
plt.ylabel("Variance σ²")
plt.title("Variance Growth (Diffusion)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT}/6_variance_growth.png", dpi=300)
pdf.savefig()
plt.close()

# Close PDF
pdf.close()

print(f"All plots generated and saved to:\n  {pdf_path}")

