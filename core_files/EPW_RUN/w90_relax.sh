#!/bin/bash




cd $DIR

					echo "====================================================="
				    	echo ">>> Running QE Geometry Relaxation For W90 & EPW:."
			

mkdir -p epw_$prefix

cd epw_$prefix
mkdir -p w90_calc

cd w90_calc
mkdir -p relax_w90

PATH1=$DIR/epw_$prefix/w90_calc/relax_w90


# Path for the log file
log_file="$DIR/check_folders.log"

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] QE geometry relax report:"
                           } >> "$DIR/EPW_work_report.txt"

# Create epw directory
#---------------------- DFT Geometric Relaxation Part ------------------------#




# Copy input
scp -r "$DIR/$prefix.relax" "$PATH1"

# Go to relax folder and run relaxation
cd "$PATH1" || exit 1
$MPIRUN_epw_PW "$prefix.relax" > relax.out

# Extract final atomic coordinates
awk '/Begin final coordinates/,/End final coordinates/' relax.out \
    | awk 'NR > 3 { print }' \
    | awk 'NR > 2 { print last } { last = $0 }' > 1

# Append last 3 lines of original input
tail -2 "$prefix.relax" >> 1

# Extract everything from &CONTROL to just before ATOMIC_POSITIONS
awk '/&CONTROL/,/ATOMIC_POSITIONS/' "$prefix.relax" \
    | awk 'NR > 1 { print last } { last = $0 }' >> 2

# Append new coordinates
cat 1 >> 2

# Replace vc-relax with scf (portable)
sed 's/vc-relax/scf/g' 2 > tmp && mv tmp 2

# Replace ibrav with 0 (portable)
sed -E "s/ibrav[[:space:]]*=[[:space:]]*$ibrav/ibrav = 0/g" 2 > tmp && mv tmp 2

# Inject force_symmorphic, lspinorb, and noncolin options
awk '/ntyp/{print;print "'"                                  force_symmorphic = $force_symmorphic"'";next}1' 2 > 3

awk '/ntyp/{print;print "'"                                          lspinorb = $lspinorb"'";next}1' 3 > 4

awk '/ntyp/{print;print "'"                                          noncolin = $noncolin"'";next}1' 4 > 5

# Clean and finalize
rm -f 1 2 3 4
mv 5 2

# Rename final input file in scf_soc directory
cd "$PATH2" || exit 1
mv 2 "$prefix.scf.in"

				    	echo ">>> Ending QE Geometry Relaxation For W90 & EPW."
					echo "====================================================="



                            {
                            echo "    → Relaxation done."
                            echo "    → Grid Used               : = $kgrid"
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of QE geometry relax report."
                            echo ""
                           } >> "$DIR/EPW_work_report.txt"

exit 0
