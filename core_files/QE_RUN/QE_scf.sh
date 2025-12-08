#!/bin/bash





cd $DIR

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] QE SCF report:"
                                } >> "$DIR/work_report.txt"

    					echo "============================================================="
				    	echo ">>> Running QE SCF [Non-Spin & With Spin] Computation."

# Ensure required variables are set
if [[ -z "$DIR" || -z "$MPIRUN_PW" || -z "$prefix" || -z "$kgrid_nonsoc" || -z "$shifted_grid" ]]; then
    echo "One or more required environment variables are not set."
    echo "Please ensure DIR, MPIRUN_PW, prefix, kgrid_nonsoc, and shifted_grid are all defined."
    exit 1
fi



# Create nscf directory
mkdir -p nscf

PATH2="$DIR/scf_soc"
PATH3="$DIR/nscf"
PATH6="$DIR/scf_nonsoc"

#--------------------------- DFT SCF-SOC Calculations ---------------------------#

cd "$PATH2" || exit 1

# Run SCF with spin-orbit coupling
$MPIRUN_PW "$prefix.scf.in" > scf.out

# Copy results to nscf directory
cp -r ./* "$PATH3"

# Also copy the input file to scf_nonsoc for modification
cp "$prefix.scf.in" "$PATH6"

#--------------------------- DFT SCF-non-SOC Calculations -----------------------#

cd "$PATH6" || exit 1

# Remove spin-related flags
sed '/lspinorb/d' "$prefix.scf.in" > tmp && mv tmp "$prefix.scf.in"
sed '/noncolin/d' "$prefix.scf.in" > tmp && mv tmp "$prefix.scf.in"

# Truncate at K_POINTS line
sed '/K_POINTS/q' "$prefix.scf.in" > tmp && mv tmp "$prefix.scf.in"

# Append new grid
echo "$kgrid_nonsoc $shifted_grid" >> "$prefix.scf.in"

# Run non-SOC SCF
$MPIRUN_PW "$prefix.scf.in" > scf.out

				    	echo ">>> Ending QE SCF [Non-Spin & With Spin] Computation."
					echo "============================================================="

                            {
                            echo "    → SCF done."
                            echo "    → Grid Used               : = $kgrid"
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of SCF report."
                            echo ""
                           } >> "$DIR/work_report.txt"


exit 0
