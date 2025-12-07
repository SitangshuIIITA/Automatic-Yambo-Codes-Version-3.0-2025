import numpy as np
import os
import sys

# --- 1. Get prefix from environment variable ---
prefix = os.environ.get("prefix")
if prefix is None:
    print("❌ Environment variable 'prefix' not set!")
    sys.exit(1)

# Input parent folders created by extraction script
main_in_folder = "electron_gkkp"
main_out_folder = f"{prefix}_elec_gkkp"
os.makedirs(main_out_folder, exist_ok=True)

# --- 2. Read qlist_monotonic.txt (with optional HS labels) ---
qlist_file = "qlist_monotonic.txt"
qlist = []
qmag_list = []
hs_labels_list = []

with open(qlist_file, "r") as f:
    for line in f:
        if line.strip() == "":
            continue
        parts = line.split()
        qmag = float(parts[0])
        qx, qy, qz = float(parts[1]), float(parts[2]), float(parts[3])
        label = parts[4] if len(parts) > 4 else ""
        qmag_list.append(qmag)
        qlist.append([qx, qy, qz])
        hs_labels_list.append(label)
qlist = np.array(qlist)
qmag_list = np.array(qmag_list)

# --- 3. Process both intraband and interband ---
for calc_type in ["intraband", "interband"]:
    in_folder = os.path.join(main_in_folder, f"{calc_type}_data")
    if not os.path.isdir(in_folder):
        print(f"⚠️ Skipping {calc_type}: folder '{in_folder}' not found")
        continue

    out_folder = os.path.join(main_out_folder, f"{calc_type}_data")
    os.makedirs(out_folder, exist_ok=True)

    # --- 4. Loop over all phonon mode files ---
    for file_name in os.listdir(in_folder):
        if not file_name.endswith(".dat"):
            continue
        mode_file_path = os.path.join(in_folder, file_name)

        gdata = []
        with open(mode_file_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split()
                iq = int(parts[0])
                qx, qy, qz = float(parts[1]), float(parts[2]), float(parts[3])
                ibnd, jbnd = int(parts[4]), int(parts[5])
                imode = int(parts[6])
                g_abs = float(parts[7])
                g_re  = float(parts[8])
                g_im  = float(parts[9])
                gdata.append({
                    "qx": qx, "qy": qy, "qz": qz,
                    "ibnd": ibnd, "jbnd": jbnd,
                    "imode": imode,
                    "g_abs": g_abs, "g_re": g_re, "g_im": g_im
                })
        gdata = np.array(gdata, dtype=object)

        # --- 5. Reorder according to qlist_monotonic.txt ---
        out_file_path = os.path.join(out_folder, file_name)
        with open(out_file_path, "w") as fout:
            fout.write("# |q|   qx   qy   qz   ibnd   jbnd   imode   |g|   Re|g|   Im|g|   HS_label\n")
            for qmag, q, label in zip(qmag_list, qlist, hs_labels_list):
                dists = np.linalg.norm(
                    np.array([[g["qx"], g["qy"], g["qz"]] for g in gdata]) - q, axis=1
                )
                idx = np.argmin(dists)
                matched = gdata[idx]
                fout.write(f"{qmag:.6f}  {q[0]:.6f}  {q[1]:.6f}  {q[2]:.6f}  "
                           f"{matched['ibnd']:3d}  {matched['jbnd']:3d}  {matched['imode']:3d}  "
                           f"{matched['g_abs']:.6e}  {matched['g_re']:.6e}  {matched['g_im']:.6e}  {label}\n")

    print(f"✅ Reordered {calc_type} g data written in folder '{out_folder}'")

