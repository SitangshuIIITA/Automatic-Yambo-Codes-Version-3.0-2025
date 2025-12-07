#!/bin/bash





cd $DIR
					echo "====================================================="
				    	echo ">>> Running QE SCF For W90 & EPW:."
			
                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] QE EPW report:"
                                } >> "$DIR/EPW_work_report.txt"

PATH1="$DIR/epw_$prefix/w90_calc"
PATH2="$DIR/epw_$prefix/w90_calc/scf_w90"
PATH3="$DIR/epw_$prefix/w90_calc/nscf_w90"

cd "$PATH1" || exit 1

mkdir -p scf_w90

cd "$PATH2" || exit 1

scp -r $PATH1/relax_w90/* ./. 

#--------------------------- DFT SCF-SOC Calculations ---------------------------#

cd "$PATH2" || exit 1

# Remove spin-related flags
sed '/lspinorb/d' "$prefix.scf.in" > tmp && mv tmp "$prefix.scf.in"
sed '/noncolin/d' "$prefix.scf.in" > tmp && mv tmp "$prefix.scf.in"

# Truncate at K_POINTS line
sed '/K_POINTS/q' "$prefix.scf.in" > tmp && mv tmp "$prefix.scf.in"

# Append new grid
echo "$kgrid_nonsoc $shifted_grid" >> "$prefix.scf.in"

# Run non-SOC SCF
$MPIRUN_epw_PW "$prefix.scf.in" > scf.out

				
				    	echo ">>> Ending QE SCF For W90 & EPW:."
					echo "====================================================="

                            {
                            echo "    → SCF done."
                            echo "    → Grid Used               : = $kgrid"
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of SCF report."
                            echo ""
                           } >> "$DIR/EPW_work_report.txt"


exit 0
