#!/bin/bash



cd $DIR

PATH5="$DIR/band_structure"
PATH15="$DIR/BSE_folder/BSE_slepC"

                           {
                            echo ""
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] BSE finite-Q exciton dispersion report:"
                           } >> "$DIR/work_report.txt"


# Output file for band path
BAND_ROUTE_FILE="$prefix.band_route"
> "$BAND_ROUTE_FILE"  # Clear existing file or create a new one

# Function to get k-point coordinates
get_kpoint_coord() {
  case "$1" in
    G) echo "0.00000 | 0.00000 | 0.00000 |" ;;
    M) echo "0.50000 | 0.00000 | 0.00000 |" ;;
    K) echo "0.33333 | 0.33333 | 0.00000 |" ;;
    *)
      echo "Error: Unknown high-symmetry point '$1'" >&2 >> "$DIR/work_report.txt"
      exit 1
      ;;
  esac
}

echo "            CALCULATING: QPATH is       : '$QPATH'" >> "$DIR/work_report.txt"
for point in $QPATH; do
  echo "            DEBUG: Processing point     : $point" >> "$DIR/work_report.txt"
  coord=$(get_kpoint_coord "$point")
  echo "            $coord" >> "$DIR/work_report.txt"
  echo "$coord" >> "$BAND_ROUTE_FILE"
done

# ------------------------- Excitonic Dispersion -------------------------
cd "$PATH15" || { echo "Directory $PATH15 not found"; exit 1; }

rm -rf $prefix.ypp_exciton_disp.in l-output* r-output_exciton* r-output_ypp*

# Run YPP to generate input file
$YPP -e i -F "$prefix.ypp_exciton_disp.in" 2> out

# Modify the input file with sed
sed -i.bak '18s/INTERP_mode= "NN"/INTERP_mode= "BOLTZ"/' "$prefix.ypp_exciton_disp.in"
sed -i.bak '19s/INTERP_Shell_Fac= 20.00000 /INTERP_Shell_Fac= 50.00000 /' "$prefix.ypp_exciton_disp.in"
sed -i.bak '21s/BANDS_steps= 10/BANDS_steps= 100/' "$prefix.ypp_exciton_disp.in"
sed -i.bak 's/States= "0 - 0"/States="'"$exciton_states"'"/' "$prefix.ypp_exciton_disp.in"

# Get number of k-points from setup file
awk '/Grid dimensions/ {print $4*3}' r_setup > kcount.tmp
kcountvar=$(cat kcount.tmp)
rm -f kcount.tmp

# Copy the generated band route into the working directory
cp "$DIR/$BAND_ROUTE_FILE" .

# Insert band route into YPP input file
sed -i.bak '35s/.*//' "$prefix.ypp_exciton_disp.in"
cat "$BAND_ROUTE_FILE" >> "$prefix.ypp_exciton_disp.in"
sed -i.bak '35d' "$prefix.ypp_exciton_disp.in"
echo "%" >> "$prefix.ypp_exciton_disp.in"

# Run YPP with modified input
$YPP -F "$prefix.ypp_exciton_disp.in" -J output 2> out

# Rename and copy the output
mv o-output.excitons_interpolated "$prefix.o-output.excitons_interpolated"
cp "$prefix.o-output.excitons_interpolated" "$DIR/graph_data/BSE_slepC_data/."

                                            {
                                                        echo "    → See detailed report in $DIR/graph_data/BSE_slepC_data."
                                                        echo "End of BSE finite-Q exciton dispersion."
                                                        echo ""
                                            } >> "$DIR/work_report.txt"
   

exit 0
