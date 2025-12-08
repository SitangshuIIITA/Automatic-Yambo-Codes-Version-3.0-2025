#!/bin/bash





cd $DIR
                         		echo "============================================================="
				    	echo ">>> Running QE Thermo PW Computation."
				    	

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] QE Thermo-pw report:"
                           } >> "$DIR/work_report.txt"



mkdir -p Thermo_PW

PATH1="$DIR/Thermo_PW"
PATH2="$DIR/scf_soc"
PATH6="$DIR/scf_nonsoc"

#--------------------------- Thermo_pw input file generation -------------------------------------#

# Copy and rename relax file
cp "$prefix.relax" "$PATH1"
cd "$PATH1"

# Replace 'vc-relax' with 'scf' (cross-platform compatible way)
sed 's/vc-relax/scf/g' "$prefix.relax" > "$prefix.mur_lc_t.in"

# Delete line containing "assume_isolated = '2D'"
grep -v "assume_isolated *= *'2D'" "$prefix.mur_lc_t.in" > tmp && mv tmp "$prefix.mur_lc_t.in"

# Generate ph_control input file
cat > ph_control << EOF
ph_title
 &inputph
  tr2_ph=$TRN_PH,
  prefix='$prefix',
  fildyn='$prefix.dyn.xml',
  ldisp=.TRUE.
  nq1=$NQ1, nq2=$NQ2, nq3=$NQ3,
/
EOF

# Generate thermo_control input file
cat > thermo_control << EOF
&INPUT_THERMO
  what='mur_lc_t',
  deltat=3.,
  step_ngeo(1) = 0.05,
  enhance_plot=.TRUE.
/
EOF

# Run Thermo_pw calculation
$MPIRUN_THERMO_PW < "$prefix.mur_lc_t.in" > "$prefix.mur_lc_t.out"

                         		
				    	echo ">>> Ending QE Thermo PW Computation."
				    	echo "============================================================="

                            {
                            echo "    → Thermo-pw done."
                            echo "    → Grid Used               : = $kgrid"
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of Thermo-pw report."
                            echo ""
                           } >> "$DIR/work_report.txt"

exit
