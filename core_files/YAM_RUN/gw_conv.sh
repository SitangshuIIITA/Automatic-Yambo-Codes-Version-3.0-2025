#!/bin/bash



cd $DIR

mkdir -p GW_folder

PATH3="$DIR/nscf"
PATH5="$DIR/band_structure"
PATH6="$DIR/GW_folder"
PATH7="$DIR/HF_folder"

#-------------------------------GW Calculation-----------------------------------#
cd "$DIR/graph_data" || exit 1
mkdir -p GW_band_structure
cd GW_band_structure || exit 1

cat > gw_conv.txt << EOF
# Ry    direct-gap      indirect-gap
EOF

cd "$PATH6" || exit 1

echo "$exx_value" > exch1
awk '{ print $(NF-0) }' exch1 > exch2
exx="$(cat exch2)"

echo "$ngsblpp_value" > exch3
ngs="$(cat exch3)"

echo "$bndrnge" > exch4
bndr="$(cat exch4)"

set -- $bndrnge  # Use the set command properly

for i in $ngsblpp_value; do
    mkdir -p "GW-$i$ryd"
    cd "GW-$i$ryd" || exit 1

    scp -r "$PATH7/HF-$exx$ryd/SAVE" ./.
    scp -r "$PATH7/HF-$exx$ryd/output" ./.

    $YAMBO 2> out
    $YAMBO -r -rw -gw0 p -g n -V all -F "$gw_filename.$i$ryd" 2> out

    #------------------------------ COULOMB cut-off & RIM ------------------------------#
    sed -i.bak 's/rimw_type= "default"/rimw_type= "'"$rimw_type"'"/g'               "$gw_filename.$i$ryd"
    sed -i.bak 's/RandGvecW= 1               RL/RandGvecW= '"$randgvecw $rl"'/g'    "$gw_filename.$i$ryd"
    sed -i.bak 's/RandQpts=0/RandQpts='"$randqpts"'/g'                              "$gw_filename.$i$ryd"
    sed -i.bak 's/RandGvec= 1                RL/RandGvec= '"$randgvec $rl"'/g'      "$gw_filename.$i$ryd"
    sed -i.bak 's/CUTGeo= "none"/CUTGeo= "'"$slab"'"/g'                             "$gw_filename.$i$ryd"
    
    #------------------------------ Self Consistent GW ----------------------------------#
    # Uncomment to enable iterative GW:
    # sed -i '126s/GWIter=0/GWIter=1/g' "${gw_filename}.${ngs}${ryd}"

    #------------------------------ Cutoffs ---------------------------------------------#

    # Insert EXXRLvcs and VXCRLvcs at line 60 (after deleting lines 60 and 61)
    awk -v exx="EXXRLvcs= $exx $ryd" -v vxc="VXCRLvcs= $exx $ryd" \
    'NR==60{print exx; print vxc; next} NR==61{next} {print}' "$gw_filename.$i$ryd" > tmp && mv tmp "$gw_filename.$i$ryd"

    # Replace line 90 with NGsBlkXp
    awk -v val="NGsBlkXp= $i $ryd" \
    'NR==90{print val; next} {print}' "$gw_filename.$i$ryd" > tmp && mv tmp "$gw_filename.$i$ryd"



    # Delete line 88 (works on both systems)
    sed '88d' "$gw_filename.$i$ryd" > tmp1 && mv tmp1 "$gw_filename.$i$ryd"

    # Insert "1| $1|" after line 88 (compatible on both macOS and Linux)
    # Use a temporary file to handle the insert
    awk "NR==88{print \"1| $1|\"}1" "$gw_filename.$i$ryd" > tmp && mv tmp "$gw_filename.$i$ryd"


    # Delete line 88 (works on both systems)
    sed '123d' "$gw_filename.$i$ryd" > tmp1 && mv tmp1 "$gw_filename.$i$ryd"

    # Insert "1| $1|" after line 123 (compatible on both macOS and Linux)
    # Use a temporary file to handle the insert
    awk "NR==123{print \"1| $1|\"}1" "$gw_filename.$i$ryd" > tmp && mv tmp          "$gw_filename.$i$ryd"


    #------------------------------ CPU parallel structure ------------------------------#
    sed -i.bak "s/SE_CPU= \"\"/SE_CPU=\"$se_cpu\"/g"                                "$gw_filename.$i$ryd"
    sed -i.bak "s/SE_ROLEs= \"\"/SE_ROLEs= \"q qp b\"/g"                            "$gw_filename.$i$ryd"

    sed -i.bak "s/X_and_IO_CPU= \"\"/X_and_IO_CPU=\"$x_and_io\"/g"                  "$gw_filename.$i$ryd"
    sed -i.bak "s/X_and_IO_ROLEs= \"\"/X_and_IO_ROLEs= \"q g k c v\"/g"             "$gw_filename.$i$ryd"

    sed -i.bak "s/DIP_CPU= \"\"/DIP_CPU=\"$dip_cpu\"/g"                             "$gw_filename.$i$ryd"
    sed -i.bak "s/DIP_ROLEs= \"\"/DIP_ROLEs= \"k c v\"/g"                           "$gw_filename.$i$ryd"

    #------------------------------ GW cut-offs ----------------------------------------#
    sed -i.bak '99s/XTermKind= "none"/XTermKind= "BG"/g'                            "$gw_filename.$i$ryd"
    sed -i.bak '127s/GTermKind= "none"/GTermKind= "BG"/g'                           "$gw_filename.$i$ryd"
    sed -i.bak '132s/#ExtendOut/ExtendOut/g'                                        "$gw_filename.$i$ryd"

    #------------------------------ Band range parsing from r_setup ---------------------#
    # Extract k-points and prepare (assuming colon-delimited format)
    awk -F':' '/IBZ K-points/ {gsub(/^ +/, "", $2); print "1|"$2}' r_setup > kpts_temp

    # Extract filled bands and compute modified VB
    awk '/Filled Bands/ {print $5}' r_setup > filled_temp
    awk -v v="$nv" '{print $1 - v}' < filled_temp > filled_parsed

    # Extract empty bands and compute modified CB
    awk '/Empty Bands/ {print $5}' r_setup > empty_temp
    awk -v v="$nc" '{print $1 + v}' < empty_temp > empty_parsed

    # Combine all into one block, format as "1|kpt|vb|cb|"
    paste -d'|' kpts_temp filled_parsed empty_parsed | awk '{print $0"|"}' > insert_block

    # Remove any blank lines (safe cross-platform way)
    grep -v '^$' insert_block > insert_block.clean && mv insert_block.clean insert_block

    # Insert into line 136 (replacing that line)
    awk -v n=136 -v f="insert_block" 'NR==n {
        while ((getline x < f) > 0) print x
        next
    } {print}' "$gw_filename.$i$ryd" > temp_file && mv temp_file "$gw_filename.$i$ryd"

    cat "insert_block" > ../correction.txt

    rm -f kpts_temp filled_temp empty_temp filled_parsed empty_parsed insert_block temp_file
    
    sed -i.bak '132s/#ExtendOut/ExtendOut/g' "$gw_filename.$i$ryd"

    $MPIRUN_YAMBO $YAMBO $YAMBO -F "$gw_filename.$i$ryd" -J output -C Report 2> out


        cd Report

            # Remove comment lines and extract columns 3 and 4
            awk '!/^#/' o-output.qp | awk '{print $3 "\t" $4}' | sort -n > tmp_sorted

            rm -f tmp_sorted

            # Create 'rydberg' file
            echo "${line}" > rydberg

            # Extract variables from r-output_HF_and_locXC_gw0_dyson_rim_cut_rim_w_em1d_ppa_el_el_corr
                awk '/Hartree-Fock occupations report/,/Dyson equation: Newton solver/' \
    r-output_HF_and_locXC_gw0_dyson_rim_cut_rim_w_em1d_ppa_el_el_corr > tmp_extract

                var1=$(awk '/Direct Gap localized at k/ {print $8}' tmp_extract)
                var2=$(awk '/Indirect Gap between kpts/ {print $7}' tmp_extract)
                var3=$(awk '/Indirect Gap between kpts/ {print $8}' tmp_extract)
                var4=$(awk '/Filled Bands/ {print $5}' tmp_extract)
                var5=$(awk '/Empty Bands/ {print $5}' tmp_extract)

            # Extract GW corrected eigenvalues
                awk -v v1="$var1" -v v4="$var4" '$1 == v1 && $2 == v4 {print $4}' o-output.qp > tmp7
                awk -v v1="$var1" -v v5="$var5" '$1 == v1 && $2 == v5 {print $4}' o-output.qp > tmp8
                awk -v v1="$var1" -v v4="$var4" '$1 == v1 && $2 == v4 {print $4}' o-output.qp > tmp9
                awk -v v1="$var3" -v v5="$var5" '$1 == v1 && $2 == v5 {print $4}' o-output.qp > tmp10

            # Combine results
                paste tmp7 tmp8 tmp9 tmp10 > tmp11
                awk '{print $2 - $1}' tmp11 > tmp12
                awk '{print $4 - $3}' tmp11 > tmp13
                paste tmp12 tmp13 > tmp14
                paste rydberg tmp14 > tmp15

            # Append to the final file
                cat tmp15 >> "$DIR/graph_data/GW_band_structure/gw_conv.txt"

            # Cleanup
                #rm -f tmp{7..15} tmp_extract rydberg

    cd ..
    cd ..
    shift
done

echo "$ngsblpp_value" > exch3
awk '{ print $(NF-0) }' exch3 > con_ryd
ngs="$(cat con_ryd)"

echo "$bndrnge" > exch4
awk '{ print $(NF-0) }' exch4 > con_bnd
bnd_lim="$(cat con_bnd)"

# Step 1: Locate the second most recent GW-* folder as input to G0W0
target_folder=$(ls -td ./GW-* | sed -n 1p)
rename_folder=$(basename "$target_folder")

{
echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Converged G0W0 report:"
echo "    → Converged G0W0 folder is from     		: $rename_folder"
echo "    → Converged response block size               : NGsBlkXp = $ngs Ry"
echo "    → Converged polarization/Screening Bands      : Scheme = 1 | $bnd_lim"
echo "    → G0W0 corrected K-points|Bands     		: Scheme = $(cat correction.txt)"
echo "    → See detailed gap report in $DIR/graph_data/GW_band_structure."
echo "    → Status                            		: G0W0 convergence calculations completed."
echo "End of converged G0W0 report."
echo ""
} >> "$DIR/work_report.txt"

exit
