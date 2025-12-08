#!/bin/bash





cd $DIR

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] QE phonon-dispersion report:"
                           } >> "$DIR/work_report.txt"


mkdir -p phonon_disp

PATH2="$DIR/scf_nonsoc"
PATH12="$DIR/phonon_disp"

cd "$PATH12"

# Create phonon input file
cat > "$prefix.phonon.in" << EOF
&inputph
    verbosity = 'high',
    outdir = './'
    tr2_ph = $TRN_PH,
    alpha_mix = $alpha_mx,
    prefix = '$prefix',
    fildyn = '$prefix.dyn',
    fildvscf = 'dvscf',
    ldisp = .true.,
    nq1 = $NQ1, nq2 = $NQ2, nq3 = $NQ3,
/
EOF

cp -r "$PATH2/$prefix.save" .

    					echo "============================================================="
				    	echo ">>> Running Phonon Dynamical Matrices & Dispersion & DOS"
				    	

$PH_DISP_RUN -inp "$prefix.phonon.in" > output_phonon

# Create q2r input
cat > "$prefix.q2r.in" << EOF
&input
    zasr = 'crystal',
    flfrc = '$prefix.fc',
    fildyn = '$prefix.dyn',
    loto_2d = .true.,
/
EOF

$Q2R_DISP_RUN -inp "$prefix.q2r.in" > q2r.out

cp -r "$DIR/bands/bands.out" .





# Copy the band route file into the bands directory
cp "$DIR/$prefix.band_route" "$PATH12/"

# Prepare the K_POINTS file
mv "$prefix.band_route" "$prefix.k_points"

# Clean and modify the k_points file (use portable method)
sed 's/|//g' "$prefix.k_points" > "$prefix.k_points.clean"
sed 's/$/ 50/' "$prefix.k_points.clean" > "$prefix.k_points.final"

# Add header with number of points and crystal type
n_lines=$(wc -l < "$prefix.k_points.final")
{
    #echo "/"
    echo "$n_lines"
    cat "$prefix.k_points.final"
} > kpoints_block.in



# Extract q-points from bands.out
awk '/cryst. coord./,/Dense/' bands.out | awk 'NR > 2 { print }' | awk 'NR>2 {print last} {last=$0}' > tmp1
awk '/cryst. coord./,/Dense/' tmp1 | awk 'NR > 1 { print }' | awk 'NR>1 {print last} {last=$0}' > tmp2
awk '//{printf "%10s %15s %15s\n", $5, $6, $7}' tmp2 > tmp3

# Remove trailing comma if any (macOS/Linux compatible)
clean_sed() {
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' -e "$1" "$2"
    else
        sed -i -e "$1" "$2"
    fi
}
clean_sed 's/),//' tmp3
cat tmp3 | wc -l > tmp4
cat tmp3 >> tmp4

# Create matdyn input
cat > "$prefix.matdyn.in" << EOF
&input
    asr = 'crystal',
    flfrc = '$prefix.fc',
    flfrq = '$prefix.freq',
    q_in_cryst_coord = .true.,
    q_in_band_form = .true.,
    loto_2d = .true.,
/
EOF

cat kpoints_block.in >> "$prefix.matdyn.in"
$MATDYN_DISP_RUN -inp "$prefix.matdyn.in" > output_matdyn

# Create DOS input
cat > "$prefix.dos.matdyn.in" << EOF
&input
    asr = 'crystal',
    dos = .true.,
    loto_2d = .true.,
    fldos = '$prefix.dos',
    flfrc = '$prefix.fc',
    flfrq = '$prefix.dos.freq',
    flvec = '$prefix.dos.modes',
    q_in_cryst_coord = .true.,
    deltaE = 2.d0,
    nk1 = 90, nk2 = 90, nk3 = 1
/
EOF

$MATDYN_DISP_RUN -inp "$prefix.dos.matdyn.in" > output_dos_matdyn

# Clean up
rm -f tmp1 tmp2 tmp3 tmp4 bands.out
					
				    	echo ">>> Phonon Dynamical Matrices & Dispersion & DOS Completed."
				    	echo "============================================================="



                            {
                            echo "    → Phonon dispersion done."
                            echo "    → q-Grid Used             : = $NQ1 x $NQ2 x $NQ3"
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of Phonon dispersion report."
                            echo ""
                           } >> "$DIR/work_report.txt"


exit 0
