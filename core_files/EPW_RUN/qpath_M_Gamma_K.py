
#!/usr/bin/env python3
import numpy as np


# --- Default parameters ---
nq = 20          # number of q-points
offset = 0.001   # offset applied ONLY to Γ
output_file = "qpath_GM.txt"

# --- High-symmetry points ---
G_y = 0.0
M_y = 0.5773502692  # as in your earlier example

# --- Interpolate qy ---
qy = np.linspace(G_y + offset, M_y, nq)

# --- Write DFPT-style file ---
with open(output_file, "w") as f:
    f.write(f"{nq}\n")
    for q in qy:
        f.write(f"0.0 {q:.10f} 0.0 1\n")

print(f"Generated {nq} q-points along Γ→M and saved to {output_file}")


