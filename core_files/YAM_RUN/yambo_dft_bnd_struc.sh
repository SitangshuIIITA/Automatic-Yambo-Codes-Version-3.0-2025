#!/bin/bash



cd $DIR

# Create necessary folders
mkdir -p band_structure graph_data

# Define paths
PATH3="$DIR/nscf"
PATH4="$DIR/bands"
PATH5="$DIR/band_structure"

# ------------------ YAMBO Interpolator: DFT Band Structure Calculations ------------------ #
cd "$DIR/graph_data"
mkdir -p DFT_band_structure

cd "$PATH5"
mkdir -p work

cd "$PATH3/$prefix.save"
$P2Y 2> out
scp -r SAVE "$PATH5/work"

cd "$PATH5/work"
$YAMBO 2> out

$YPP -s b -V qp -F "$prefix.DFT_bands.in" 2> out

# ------------------ Reading the k-point limits, max valence band, min conduction band from r_setup file ------------------ #

# Step 1: Extract Filled and Empty Bands
awk '/Filled Bands/ {print $5}' r_setup > 5
awk -v v=$nv '{print $1 - v, $2}' 5 > 6

awk '/Empty Bands/ {print $5}' r_setup > 7
awk -v v=$nc '{print $1 + v, $2}' 7 > 8

# Step 2: Create file 9 with '|'
paste -d'|' 6 8 | sed 's/$/|/' > 9

# Step 3: Insert contents of 9 at line 27, replacing original line 27
awk -v n=27 -v file=9 '
NR==n {
    while ((getline line < file) > 0) print line
    next  # skip original line 27
}
{print}
' "$prefix.DFT_bands.in" > tmp && mv tmp "$prefix.DFT_bands.in"

# Step 4: Modify BZ route for DFT band diagram
sed -e 's/NN/BOLTZ/g' -e 's/BANDS_steps= 10/BANDS_steps= 100/g' "$prefix.DFT_bands.in" > tmp && mv tmp "$prefix.DFT_bands.in"
sed "s|BANDS_path= \"\"|BANDS_path=\"$bands_path\"|g" "$prefix.DFT_bands.in" > tmp && mv tmp "$prefix.DFT_bands.in"

# Step 5: Clean-up
rm -f 1 2 3 4 5 6 7 8 9 new_file.txt 1new_file.txt

# Step 1: Copy band_route file (scp not needed if local)
cp "$DIR/$prefix.band_route" .

# Step 2: Insert the contents of prefix.band_route at line 52, replacing it
awk -v n=52 -v file="$prefix.band_route" '
NR==n {
    while ((getline line < file) > 0) print line
    next  # skip original line 52
}
{print}
' "$prefix.DFT_bands.in" > tmp && mv tmp "$prefix.DFT_bands.in"

# Step 3: Append '%' at the end of DFT_bands.in
echo '%' >> "$prefix.DFT_bands.in"

$YPP -F "$prefix.DFT_bands.in" 2> out
mv o.bands_interpolated "$prefix.DFT.o.bands_interpolated"
scp -r "$prefix.DFT.o.bands_interpolated" "$DIR/graph_data/DFT_band_structure/."

rm -f l_electrons_bnds l_setup

# ------------------ Convert rows with line breaks into multiple columns (if necessary) ------------------ #
# This block can be uncommented if it's strongly required, but it is time-consuming.
#
# cd "$PATH4"
# scp -r Bandx.dat.gnu "$PATH5/work/FixSymm"
#
# cd "$PATH5/work/FixSymm"
#
# mv Bandx.dat.gnu "$prefix.Bandx.dat.gnu"
#
# sed 's/^ *//' "$prefix.Bandx.dat.gnu" > tmp && mv tmp "$prefix.Bandx.dat.gnu"
#
# sed -i '$ d' "$prefix.Bandx.dat.gnu"
#
# sed -E '/\S/!d;H;g;s/((\n\S+)(\s+)[^\n]*)(.*)\2\3(.*)$/\1\3\5\4/;h;$!d;x;s/.//' "$prefix.Bandx.dat.gnu" > 1
#
# mv 1 "$prefix.Bandx.dat.gnu"
#
# scp -r "$prefix.Bandx.dat.gnu" "$DIR/graph_data/DFT_band_structure/."
#
# rm -rf "$prefix.Bandx.dat.gnu"

exit


#$YPP -y -F rmv_symm.in 2> out

# Use a temporary file safely
#tmpfile=$(mktemp /tmp/rmv_symm.XXXXXX)

#sed 's/#RmAllSymm/RmAllSymm/g' "$PATH5/work/rmv_symm.in" > "$tmpfile" && mv "$tmpfile" "$PATH5/work/rmv_symm.in"

#$YPP -F rmv_symm.in 2> out

#cd "$PATH5/work/FixSymm"
#$YAMBO 2> out
#$YPP -s b -V qp -F "$prefix.DFT_bands.in" 2> out
