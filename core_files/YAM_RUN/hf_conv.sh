#!/bin/bash



cd $DIR

mkdir -p HF_folder

PATH3="$DIR/nscf"
PATH5="$DIR/band_structure"
PATH1="$DIR/HF_folder"
PATH2="$DIR/nscf/$prefix.save"

# Note: double-grid.ndb is inside the SAVE folder. But it is ineffective in HF computation. However, do not delete this.

# ------------------------------ HF Calculation ----------------------------- #

cd "$DIR/graph_data"
mkdir -p HF_band_structure
cd HF_band_structure

cat > hf_conv.txt << EOF
# Ry    direct-gap    indirect-gap
EOF

cd "$PATH1"

for line in $exx_value; do
    mkdir -p "HF-${line}${ryd}"
    cd "HF-${line}${ryd}"
    
    scp -r "$PATH2/SAVE" .
        #cd SAVE
        #    rm -f ndb.Double_Grid
        #cd ..

    rm -rf "${line}${ryd}.${hf_filename}"
    
    $YAMBO 2> out
    $YAMBO -r -rw -x -V par -F "${line}${ryd}.${hf_filename}" 2> out

    # ----------------------------- COULOMB cut-off & RIM ----------------------- #

    sed -i.bak 's/rimw_type= "default"/rimw_type= "'"$rimw_type"'"/g' "${line}${ryd}.${hf_filename}"
    sed -i.bak 's/RandGvecW= 1               RL/RandGvecW= '"$randgvecw $rl"'/g' "${line}${ryd}.${hf_filename}"
    sed -i.bak 's/RandQpts=0/RandQpts='"$randqpts"'/g' "${line}${ryd}.${hf_filename}"
    sed -i.bak 's/RandGvec= 1                RL/RandGvec= '"$randgvec $rl"'/g' "${line}${ryd}.${hf_filename}"
    sed -i.bak 's/CUTGeo= "none"/CUTGeo= "'"$slab"'"/g' "${line}${ryd}.${hf_filename}"

    # ------------------------- CPU parallel structure --------------------------- #

    sed -i.bak 's/SE_CPU= ""/SE_CPU="'"$se_cpu"'"/g' "${line}${ryd}.${hf_filename}"
    sed -i.bak 's/SE_ROLEs= ""/SE_ROLEs= "q qp b"/g' "${line}${ryd}.${hf_filename}"

    # ------------------------------- HF cut-offs -------------------------------- #

# Delete line 35 twice (if needed)
sed  '36d' "${line}${ryd}.${hf_filename}" > tmp1 && mv tmp1 "${line}${ryd}.${hf_filename}"
sed  '36d' "${line}${ryd}.${hf_filename}" > tmp1 && mv tmp1 "${line}${ryd}.${hf_filename}"

# Insert the two lines above line 35
awk -v exx="EXXRLvcs= ${line} $ryd" -v vxc="VXCRLvcs= ${line} $ryd" 'NR==34{print exx; print vxc} {print}' "${line}${ryd}.${hf_filename}" > tmp2 && mv tmp2 "${line}${ryd}.${hf_filename}"


    # ----------- Reading the k-point limits, max valence band, min conduction band --------------- #

# Extract k-points and prepare
awk -F':' '/IBZ K-points/ {gsub(/^ +/, "", $2); print "1|"$2}' r_setup > kpts_temp

# Extract filled bands
awk '/Filled Bands/ {print $5}' r_setup > filled_temp
awk -v v="$nv" '{print $1 - v, $2}' filled_temp > filled_parsed

# Extract empty bands
awk '/Empty Bands/ {print $5}' r_setup > empty_temp
awk -v v="$nc" '{print $1 + v, $2}' empty_temp > empty_parsed

# Now paste together without extra sed tricks
paste -d'|' kpts_temp filled_parsed empty_parsed | awk '{print $0"|"}' > insert_block

# Make sure no blank lines are present
sed '/^[[:space:]]*$/d' insert_block > tmp && mv tmp insert_block


awk -v n=39 -v f="insert_block" 'NR==n {
    while ((getline x < f) > 0) print x
    next
} {print}' "${line}${ryd}.${hf_filename}" > temp_file && mv temp_file "${line}${ryd}.${hf_filename}"

# Clean temporary files
rm -f kpts_temp filled_temp filled_parsed empty_temp empty_parsed insert_block *.bak


    # ------------------------- Running Yambo --------------------------- #

$MPIRUN_YAMBO yambo -F ${line}$ryd.$hf_filename 2> out -J output -C Report

cd Report

# Move the Yambo output file with proper naming
mv r-output_HF_and_locXC_rim_cut_rim_w "${prefix}.${line}${ryd}.r-output_HF_and_locXC_rim_cut_rim_w"

# Extract the section between 'Hartree-Fock occupations report' and 'Timing Overview'
awk '/Hartree-Fock occupations report/,/Timing Overview/' "${prefix}.${line}${ryd}.r-output_HF_and_locXC_rim_cut_rim_w" > report_section

# Extract Direct Gap and Indirect Gap
awk '/Direct Gap/ {print $5}' report_section | sed '$d' > direct_gap
awk '/Indirect Gap/ {print $5}' report_section | sed '$d' > indirect_gap

# Save the line value into rydberg file
echo "${line}" > rydberg

# Prepare the final output
paste rydberg direct_gap indirect_gap >> "${DIR}/graph_data/HF_band_structure/hf_conv.txt"

# Cleanup temporary files
rm -f report_section direct_gap indirect_gap rydberg

# Copy the renamed Yambo output to the HF_band_structure folder
scp "${prefix}.${line}${ryd}.r-output_HF_and_locXC_rim_cut_rim_w" "${DIR}/graph_data/HF_band_structure/."


    cd ..
    cd ..
done

echo "$exx_value" > exch0
awk '{ print $(NF-0) }' exch0 > con_exx
exx_lim="$(cat con_exx)"

{
echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Converged HF report:"
echo "    → Converged EXXRLvcs                : EXXRLvcs = $exx_lim Ry"
echo "    → Converged EXXRLvcs                : VXCRLvcs = $exx_lim Ry"
echo "    → See detailed gap report in $DIR/graph_data/HF_band_structure."
echo "    → Status                            : HF calculations completed."
echo "End of converged HF report."
echo ""
} >> "$DIR/work_report.txt"


exit 0
