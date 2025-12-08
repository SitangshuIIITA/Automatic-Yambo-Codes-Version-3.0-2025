#!/usr/bin/env python3
import numpy as np
import sys
import os

# --- 1. Get prefix from environment variable ---
prefix = os.environ.get("prefix")
if prefix is None:
    print("❌ Environment variable 'prefix' not set!")
    sys.exit(1)

# --- 2. Read reciprocal lattice vectors from $prefix.mat.out ---
mat_file = f"{prefix}.elec_self.out"
b = []
with open(mat_file, "r") as f:
    for line in f:
        if line.strip().startswith("b("):
            parts = line.split("=")[1].replace("(", "").replace(")", "").split()
            b.append([float(parts[0]), float(parts[1]), float(parts[2])])
b = np.array(b)  # shape (3,3)
b1, b2, b3 = b[0], b[1], b[2]

# --- 3. Read high-symmetry points from $prefix.band_route ---
route_file = f"{prefix}.band_route"
hs_points = []
with open(route_file, "r") as f:
    for line in f:
        if "|" in line:
            parts = line.split("|")
            hs_points.append([float(p.strip()) for p in parts[:3]])
hs_points = np.array(hs_points)

# --- 4. Read qlist.txt ---
q_red = []
with open("klist_BZ.dat", "r") as f:
    for line in f:
        if "coord." in line:
            parts = line.split()
            qx, qy, qz = float(parts[-3]), float(parts[-2]), float(parts[-1])
            q_red.append([qx, qy, qz])
q_red = np.array(q_red)

# --- 5. Convert to Cartesian coordinates ---
q_cart = np.array([q[0]*b1 + q[1]*b2 + q[2]*b3 for q in q_red])

# --- 6. Compute cumulative |q| distances ---
distances = [0.0]
for i in range(1, len(q_cart)):
    dq = np.linalg.norm(q_cart[i] - q_cart[i-1])
    distances.append(distances[-1] + dq)
distances = np.array(distances)

# --- 7. Identify indices of closest q-points to high-symmetry points ---
hs_q_indices = []
for pt in hs_points:
    # Convert hs_point to Cartesian
    pt_cart = pt[0]*b1 + pt[1]*b2 + pt[2]*b3
    # Find nearest q in q_cart
    idx = np.argmin(np.linalg.norm(q_cart - pt_cart, axis=1))
    hs_q_indices.append(idx)

# Optional: Labels for high-symmetry points
hs_labels = ["#G", "#M", "#K", "#G"]  # adjust according to your band_route
hs_dict = {idx: label for idx, label in zip(hs_q_indices, hs_labels)}

# --- 8. Save output with |q| axis and HS labels ---
out_file = "klist_monotonic.txt"
with open(out_file, "w") as f:
    for i, (d, q) in enumerate(zip(distances, q_red)):
        label = hs_dict.get(i, "")
        f.write(f"{d:.6f}  {q[0]:.6f}  {q[1]:.6f}  {q[2]:.6f}  {label}\n")


