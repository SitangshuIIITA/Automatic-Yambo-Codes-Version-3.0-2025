#!/bin/bash





cd $DIR
					echo "====================================================="
				    	echo ">>> Running QE Phonon Dynamical Matrices For EPW: Don't Enter Anything. I am Automatic."
							    	

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] QE phonon-dispersion report:"
                           } >> "$DIR/EPW_work_report.txt"


mkdir -p $DIR/epw_$prefix/w90_calc/w90_phonon_disp

PATH2="$DIR/epw_$prefix/w90_calc/scf_w90"
PATH12="$DIR/epw_$prefix/w90_calc/w90_phonon_disp"

cd "$PATH12"

# Create phonon input file
cat > "$prefix.phonon.in" << EOF
&inputph
    verbosity = 'high',
    outdir = './'
    tr2_ph = $TRN_PH,
    alpha_mix = $alpha_mx,
    prefix = '$prefix',
    fildyn = '$prefix.dyn.xml',
    fildvscf = 'dvscf',
    ldisp = .true.,
    epsil = .true.,
    nq1 = $NQ1, nq2 = $NQ2, nq3 = $NQ3,
/
EOF

cp -r "$PATH2/$prefix.save" .


$EPW_PH_DISP_RUN -inp "$prefix.phonon.in" > output_phonon

cp -r "$DIR/bands/bands.out" .


# Extract q-points from bands.out
awk '/cryst. coord./,/Dense/' bands.out | awk 'NR > 2 { print }' | awk 'NR>2 {print last} {last=$0}' > tmp1
awk '/cryst. coord./,/Dense/' tmp1 | awk 'NR > 1 { print }' | awk 'NR>1 {print last} {last=$0}' > tmp2
awk '//{printf "%10s %15s %15s %15s\n", $5, $6, $7, $10}' tmp2 > tmp3

# Remove trailing comma if any (macOS/Linux compatible)
clean_sed() {
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' -e "$1" "$2"
    else
        sed -i -e "$1" "$2"
    fi
}
clean_sed 's/),//' tmp3
cat tmp3 | wc -l > tmp4.txt
echo "$(wc -l < tmp3) crystal" > tmp4.txt
cat tmp3 >> tmp4.txt
mv tmp4.txt filkf.txt

python3 $EPW_PATH_py/pp.py << EOF
$prefix
EOF

				    	echo ">>> Ending QE Phonon Dynamical Matrices For EPW."
				    	echo "====================================================="

					
                            {
                            echo "    → Phonon dispersion done."
                            echo "    → q-Grid Used             : = $NQ1 x $NQ2 x $NQ3"
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of Phonon dispersion report."
                            echo ""
                           } >> "$DIR/EPW_work_report.txt"


exit 0
