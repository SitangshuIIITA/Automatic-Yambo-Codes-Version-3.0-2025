import numpy as np

nk1 = 12
nk2 = 12
nq = nk1 * nk2
wq = 1.0 / nq

with open("q_BZ.dat", "w") as f:
    # Header
    f.write(f"{nq} crystal\n")
    
    # Grid points
    for i in range(nk1):
        for j in range(nk2):
            qx = i / nk1
            qy = j / nk2
            qz = 0.0
            f.write(f"{qx:.12f}  {qy:.12f}  {qz:.12f}  {wq:.12f}\n")

print("Generated q.dat with header + 144 fractional q-points.")


