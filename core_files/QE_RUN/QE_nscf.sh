#!/bin/bash





cd $DIR


                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] QE NSCF report:"
                                } >> "$DIR/work_report.txt"


					echo "============================================================="
				    	echo ">>> Running QE NSCF [With Spin/Non-Spin] Computation."
				    	


# Create bands directory if it doesn't exist
mkdir -p bands

# Define working paths
PATH3="$DIR/nscf"
PATH4="$DIR/bands"
PATH_SC="$DIR/scf_soc"

# Create nscf directory and copy from scf if nscf doesn't exist
if [ ! -d "$PATH3" ]; then
    echo "nscf directory not found, creating it and copying from scf..."
    mkdir -p "$PATH3"
    cp -r "$PATH_SC"/* "$PATH3"/
fi

# Go to nscf directory
cd "$PATH3" || { echo "Error: Could not enter $PATH3."; exit 1; }

# Clean up previous outputs
rm -f scf.out nscf.out

# Check if SCF input exists, or copy from scf directory
if [[ ! -f "$prefix.scf.in" ]]; then
    echo "Input file $prefix.scf.in not found in $PATH3."
    if [[ -f "$PATH_SC/$prefix.scf.in" ]]; then
        echo "Copying $prefix.scf.in from $PATH_SC to $PATH3..."
        cp "$PATH_SC/$prefix.scf.in" "$PATH3/"
    else
        echo "Error: $prefix.scf.in not found in $PATH_SC either. Exiting."
        exit 1
    fi
fi

# Rename SCF input to NSCF input
mv "$prefix.scf.in" "$prefix.nscf.in"

# Replace calculation = 'scf' with 'nscf'
awk '{ if ($0 ~ /calculation *= *'\''scf'\''/) gsub(/scf/, "nscf"); print }' "$prefix.nscf.in" > tmp && mv tmp "$prefix.nscf.in"

# Insert number of bands (NBND) after 'ntyp' line
awk -v var="$NBND" '/ntyp/{print;print var;next}1' "$prefix.nscf.in" > tmp && mv tmp "$prefix.nscf.in"

# Truncate file after K_POINTS and append new k-grid
awk '/K_POINTS/ { print; exit } { print }' "$prefix.nscf.in" > tmp && mv tmp "$prefix.nscf.in"
echo "$kgrid $shifted_grid" >> "$prefix.nscf.in"

# Cleanup
rm -f tmp

# Run NSCF calculation
$MPIRUN_nscf "$prefix.nscf.in" > nscf.out

# Copy all outputs to bands directory
scp -r * "$PATH4"

				    	echo ">>> Ending QE NSCF [With Spin/Non-Spin] Computation."
				    	echo "============================================================="
       {
                            echo "    → NSCF done."
                            echo "    → Grid Used               : = $kgrid"
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of NSCF report."
                            echo ""
                           } >> "$DIR/work_report.txt"


exit 0
