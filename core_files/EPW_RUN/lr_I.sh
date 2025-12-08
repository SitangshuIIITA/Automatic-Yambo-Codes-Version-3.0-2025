#!/bin/bash





source ../../core_files/EPW_VAR/EPW_variables.sh
source ../../core_files/QE_VAR/QE_variables.sh

cd $DIR

					echo "=========================================================================="
				    	echo ">>> Running Dynamical Matrices For Quadrupole Correction File Generation: quadrupole.fmt"
			

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] Quadrupole report:"
                                } >> "$DIR/EPW_work_report.txt"


mkdir -p epw_$prefix/various_lr


PATH1="$DIR/epw_$prefix/w90_calc"
PATH2="$PATH1/relax_w90"
PATH3="$DIR/epw_$prefix/Quadrupole"
PATH4="$DIR/epw_$prefix/various_lr"


scp -r ../core_files/EPW_RUN/prepare_q.py $PATH4/.


cd "$PATH4" || exit 1

python3 prepare_q.py --prepare-q


###################	Relax computation 	############################################

scp -r $PATH2/$prefix.relax ./. 
sed -i "/occupations *= *'fixed'/d" "$prefix.relax"

awk -v var="					   outdir  = '"work"' " '/calculation/{print;print var;next}1' "$prefix.relax" > tmp && mv tmp "$prefix.relax"



		$MPIRUN_epw_PW "$prefix.relax" > relax.out





###################	scf computation 	############################################

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

# Inject lspinorb, and noncolin options

awk '/ntyp/{print;print "'"                                          lspinorb = $lspinorb"'";next}1' 2 > 3

awk '/ntyp/{print;print "'"                                          noncolin = $noncolin"'";next}1' 3 > 4


mv 4 2
mv 2 "$prefix.scf.in"

# Clean and finalize
rm -f 1 2 3 4


		$MPIRUN_epw_PW "$prefix.scf.in" > scf.out






###################	phonon calculation  ##################################
# Create phonon input file
cat > "ph.in" << EOF
&inputph
    verbosity = 'high',
       outdir = 'work'
       tr2_ph = $TRN_PH,
    alpha_mix = $alpha_mx,
       prefix = '$prefix',
       fildyn = 'dyn',
     fildvscf = 'dvscf',
        ldisp = .true.,
        epsil = .true.
    nq1 = $NQ1, nq2 = $NQ2, nq3 = $NQ3,
/
EOF

		$EPW_PH_DISP_RUN -inp "ph.in" > ph.out


cat > "q2r_3d.in" << EOF
&INPUT
    fildyn = 'dyn'
    flfrc = '3d.ifc'

    zasr = 'crystal'

    loto_2d = .false.
/
EOF

cat > "q2r_gaussian.in" << EOF
&INPUT
    fildyn = 'dyn'
    flfrc = 'gaussian.ifc'

    zasr = 'crystal'

    loto_2d = .true.
/
EOF

cat > "matdyn_3d.in" << EOF
&INPUT
flfrc = '3d.ifc'
flfrq = '3d.freq'
fldos = ' '
fleig = ' '
flvec = ' '
asr = 'crystal'
loto_2d = False
q_in_cryst_coord = True
/
EOF
echo "$(cat q.dat)" >> matdyn_3d.in


cat > "matdyn_gaussian.in" << EOF
&INPUT
flfrc = 'gaussian.ifc'
flfrq = 'gaussian.freq'
fldos = ' '
fleig = ' '
flvec = ' '
asr = 'crystal'
loto_2d = True
q_in_cryst_coord = True
/
EOF
echo "$(cat q.dat)" >> matdyn_gaussian.in


for lr in '3d' 'gaussian'
do

    $EPW_Q2R_DISP_RUN -inp q2r_$lr.in > q2r_$lr.out
    $EPW_MATDYN_DISP_RUN -inp matdyn_$lr.in > matdyn_$lr.out

done

		ph2epw



perl $PERL_PATH $kgrid_nonsoc > QE_k_list.in

# Rename SCF input to NSCF input
scp -r "$prefix.scf.in" "$prefix.nscf.in"

# Replace calculation = 'scf' with 'nscf'
awk '{ if ($0 ~ /calculation *= *'\''scf'\''/) gsub(/scf/, "nscf"); print }' "$prefix.nscf.in" > tmp && mv tmp "$prefix.nscf.in"

# Delete line containing "force_symmorphic = .TRUE."
grep -v "force_symmorphic = .TRUE." "$prefix.nscf.in" > tmp && mv tmp "$prefix.nscf.in"


# Insert number of bands (NBND) after 'ntyp' line
awk -v var="$NBND" '/ntyp/{print;print var;next}1' "$prefix.nscf.in" > tmp && mv tmp "$prefix.nscf.in"

# Truncate file after K_POINTS and append new k-grid
awk '/K_POINTS {automatic}/ {exit} {print}' "$prefix.nscf.in" > tmp && mv tmp "$prefix.nscf.in"

echo "$(cat QE_k_list.in)" >> "$prefix.nscf.in"


		$MPIRUN_epw_nscf "$prefix.nscf.in" > nscf.out

					
				    	echo ">>> Ending Dynamical Matrices Calculations For quadrupole.fmt File Generation."
					echo "=========================================================================="


                       {    echo "    → Quadrupole generated."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of Quadrupole.fmt report."
			    echo ""
                           } >> "$DIR/EPW_work_report.txt"
	

exit 0
