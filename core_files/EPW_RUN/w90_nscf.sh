#!/bin/bash





cd $DIR

					echo "====================================================="
				    	echo ">>> Running QE NSCF At Uniform Grids For W90 & EPW:."
			

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] QE EPW report:"
                                } >> "$DIR/EPW_work_report.txt"


mkdir -p epw_$prefix/w90_calc/nscf_w90

mkdir -p $DIR/epw_$prefix/W90_plots


PATH1="$DIR/epw_$prefix/w90_calc"
PATH2="$DIR/epw_$prefix/w90_calc/scf_w90"
PATH3="$DIR/epw_$prefix/w90_calc/nscf_w90"

cd "$PATH3" || exit 1

scp -r $PATH2/* ./. 

#--------------------------- DFT SCF-SOC Calculations ---------------------------#


perl $PERL_PATH $kgrid_nonsoc > QE_k_list.in

perl $PERL_PATH $kgrid_nonsoc wan > epw_k_list.in


# Rename SCF input to NSCF input
scp -r "$prefix.scf.in" "$prefix.epw_nscf.in"

# Replace calculation = 'scf' with 'nscf'
awk '{ if ($0 ~ /calculation *= *'\''scf'\''/) gsub(/scf/, "nscf"); print }' "$prefix.epw_nscf.in" > tmp && mv tmp "$prefix.epw_nscf.in"

# Delete line containing "force_symmorphic = .TRUE."
grep -v "force_symmorphic = .TRUE." "$prefix.epw_nscf.in" > tmp && mv tmp "$prefix.epw_nscf.in"


# Insert number of bands (NBND) after 'ntyp' line
awk -v var="$NBND" '/ntyp/{print;print var;next}1' "$prefix.epw_nscf.in" > tmp && mv tmp "$prefix.epw_nscf.in"

# Truncate file after K_POINTS and append new k-grid
awk '/K_POINTS {automatic}/ {exit} {print}' "$prefix.epw_nscf.in" > tmp && mv tmp "$prefix.epw_nscf.in"

echo "$(cat QE_k_list.in)" >> "$prefix.epw_nscf.in"

# Run NSCF calculation
$MPIRUN_epw_nscf "$prefix.epw_nscf.in" > epw_nscf.out

scp -r ../../../../core_files/EPW_RUN/analyze_windows.py ./.

python3 analyze_windows.py --nscf epw_nscf.out --num_wann $num_wann > wannier_dis_details.txt
scp -r wannier_dis_details.txt ../../W90_plots/.

# Cleanup
rm -f tmp

					echo ">>> Ending QE NSCF At Uniform Grids For W90 & EPW:."
					echo "====================================================="
				    	

                            {
                            echo "    → NSCF for Wannierization done."
                            echo "    → Grid Used               : = $kgrid"
                            echo "    → See W90_plots folder in $DIR/epw_{"$prefix"}."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of Wannier report."
                            echo ""
                           } >> "$DIR/EPW_work_report.txt"


exit 0
