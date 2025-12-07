#!/bin/bash





cd "$DIR"

                         		echo "============================================================="
				    	echo ">>> Running QE Electronic Band Structures Computation."
				    	



{
    echo ""
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] QE bands report:"
} >> "$DIR/work_report.txt"

PATH4="$DIR/bands"
PATH3="$DIR/nscf"

# ------------------------------------ QE Band Structure Calculations -------------------------

echo "[$(date)] Starting QE Band Structure Calculation for $prefix" >> "$DIR/qe_band_calc.log"

# Ensure the bands folder is empty
rm -rf "$PATH4"/*
 
# Copy everything from nscf to bands folder
cp -r "$PATH3"/* "$PATH4"/

cd "$PATH4"

# Move the nscf input file to the bands input file
mv "$prefix.nscf.in" "$prefix.bands.in"

# Replace 'nscf' with 'bands' in the bands input file safely
awk '{ if ($0 ~ /calculation *= *'\''nscf'\''/) gsub(/nscf/, "bands"); print }' "$prefix.bands.in" > tmp && mv tmp "$prefix.bands.in"

# Backup original file
cp "$prefix.bands.in" "$prefix.bands.in.bak"

# Delete the 'nbnd' line if it exists
sed '' "$prefix.bands.in.bak" > "$prefix.bands.in.tmp"

# Remove lines after 'K_POINTS' (inclusive)
sed '/K_POINTS/,$d' "$prefix.bands.in.tmp" > "$prefix.bands.in"

# Cleanup temporary file
rm -f "$prefix.bands.in.tmp"

# Copy the band route file into the bands directory
cp "$DIR/$prefix.band_route" "$PATH4/"

# Prepare the K_POINTS file
mv "$prefix.band_route" "$prefix.k_points"

# Clean and modify the k_points file (use portable method)
sed 's/|//g' "$prefix.k_points" > "$prefix.k_points.clean"
sed 's/$/ 50/' "$prefix.k_points.clean" > "$prefix.k_points.final"

# Add header with number of points and crystal type
n_lines=$(wc -l < "$prefix.k_points.final")
{
    echo "K_POINTS {crystal_b}"
    echo "$n_lines"
    cat "$prefix.k_points.final"
} > kpoints_block.in

# Append the modified K_POINTS block to the bands input
cat kpoints_block.in >> "$prefix.bands.in"

# Clean up temporary files
rm -f "$prefix.k_points" "$prefix.k_points.clean" "$prefix.k_points.final" kpoints_block.in

# Run the band structure calculation
$MPIRUN_Bnds "$prefix.bands.in" > bands.out

                         		
				    	echo ">>> Ending QE Electronic Band Structures Computation."
				    	echo "============================================================="


{
    echo "    → Bands done."
    echo "    → Grid Used               : = $kgrid"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of bands report."
    echo ""
} >> "$DIR/work_report.txt"

exit

