#!/usr/bin/env python3
import os
import re
import numpy as np

# ------------------------------------------------------------
# Generic parser for a given base directory (electron or hole)
# ------------------------------------------------------------
def parse_emobility_file(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

    with open(path, "r") as f:
        data = f.read()

    # Split into SERTA and BTE parts (EPW labels like ==== BTE ===)
    sections = re.split(r"=+\s*BTE\s*=+", data, flags=re.IGNORECASE)
    serta_block = sections[0]
    bte_block = sections[1] if len(sections) > 1 else ""

    def extract(block):
        temps, mu, mu_hall = [], [], []
        temp_blocks = re.split(r"Temperature:\s+([\d.]+)\s*K", block)
        for i in range(1, len(temp_blocks), 2):
            T = float(temp_blocks[i])
            content = temp_blocks[i+1]

            # Get line after "Mobility tensor without magnetic field"
            m = re.search(r"Mobility tensor without magnetic field.*?\n(.*)", content)
            if not m:
                continue
            line = m.group(1).strip()

            if "|" not in line:
                continue
            left, right = line.split("|", 1)

            try:
                mu_val = float(left.strip().split()[0])
            except:
                mu_val = np.nan

            try:
                hall_vals = re.findall(r"[-+]?\d*\.\d+E[+-]?\d+", right)
                if hall_vals:
                    mu_hall_val = max([abs(float(v)) for v in hall_vals])
                else:
                    mu_hall_val = np.nan
            except:
                mu_hall_val = np.nan

            temps.append(T)
            mu.append(mu_val)
            mu_hall.append(mu_hall_val)

        temps, mu, mu_hall = np.array(temps), np.array(mu), np.array(mu_hall)
        with np.errstate(divide='ignore', invalid='ignore'):
            rH = np.where(mu != 0, mu_hall / mu, np.nan)
        return temps, mu, mu_hall, rH

    t_s, mu_s, muh_s, rH_s = extract(serta_block)
    t_b, mu_b, muh_b, rH_b = extract(bte_block)
    return np.array(t_s), mu_s, mu_b, rH_s, rH_b


def generate_summary(base_dir, output_file):
    DIPOLE_PATH = os.path.join(base_dir, "e-Mobility_dipole_sp" if "e-" in base_dir else "h-Mobility_dipole_sp", 
                               "emob_dipole_sp.out" if "e-" in base_dir else "hmob_dipole_sp.out")
    QUAD_PATH = os.path.join(base_dir, "e-Mobility_quadrupole" if "e-" in base_dir else "h-Mobility_quadrupole", 
                             "emob_quadrupole.out" if "e-" in base_dir else "hmob_quadrupole.out")

    t_d, mu_serta_d, mu_bte_d, rH_serta_d, rH_bte_d = parse_emobility_file(DIPOLE_PATH)
    t_q, mu_serta_q, mu_bte_q, rH_serta_q, rH_bte_q = parse_emobility_file(QUAD_PATH)

    if not np.allclose(t_d, t_q):
        print(f"⚠️ Warning: Temperature mismatch between dipole and quadrupole in {base_dir}")

    temps = t_d

    with open(output_file, "w") as f:
        header = (
            "#  Temperature(K)    mu_SERTA(Dipole)      mu_BTE(Dipole)     r_H_SERTA(Dipole)       r_H_BTE(Dipole)      "
            "mu_SERTA(Quad)        mu_BTE(Quad)       r_H_SERTA(Quad)         r_H_BTE(Quad)\n"
            + "=" * 140 + "\n"
        )
        f.write(header)
        for i, T in enumerate(temps):
            f.write(
                f"{T:12.1f}"
                f"{mu_serta_d[i]:20.4f}{mu_bte_d[i]:20.4f}"
                f"{rH_serta_d[i]:24.6f}{rH_bte_d[i]:24.6f}"
                f"{mu_serta_q[i]:20.4f}{mu_bte_q[i]:20.4f}"
                f"{rH_serta_q[i]:24.6f}{rH_bte_q[i]:24.6f}\n"
            )

    print(f"✅ {output_file} written with {len(temps)} points.")


# ------------------------------------------------------------
# Main Execution
# ------------------------------------------------------------
if __name__ == "__main__":
    BASE_PATH_E = "./e-Mobility_vs_T"
    BASE_PATH_H = "./h-Mobility_vs_T"

    if os.path.exists(BASE_PATH_E):
        generate_summary(BASE_PATH_E, "mobility_summary_electron.txt")

    if os.path.exists(BASE_PATH_H):
        generate_summary(BASE_PATH_H, "mobility_summary_hole.txt")

    print("✅ All mobility summaries generated successfully.")

