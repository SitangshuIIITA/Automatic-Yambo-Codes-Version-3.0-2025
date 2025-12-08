#!/bin/bash



cd $DIR



PATH6="$DIR/GW_folder"

cd "$PATH6" || exit 1

mkdir -p GW_selfC
scp -r exch3 GW_selfC
cd GW_selfC || exit 1
awk '{ print $(NF-0) }' exch3 > con_ryd
con="$(cat con_ryd)"

echo "$bndrnge" > exch4
awk '{ print $(NF-0) }' exch4 > con_bnd
bnd_lim="$(cat con_bnd)"

finmane="$gw_filename.$con$ryd"

# Step 1: Locate the second most recent GW-* folder as input to G0W0
target_folder=$(ls -td ../GW-* | sed -n 1p)

if [ ! -d "$target_folder" ]; then
    echo "Error: Could not find a valid target folder (../GW-*) to initialize G0W0."
    exit 1
fi

echo "Initializing G0W0 from $target_folder..."
cp -r "$target_folder" G0W0
cd G0W0

        cd Report

            # Remove comment lines and extract columns 3 and 4
            awk '!/^#/' o-output.qp | awk '{print $3 "\t" $4}' | sort -n > tmp_sorted

            # Split into positive and negative energy files
            awk '$2 !~ /^-/' tmp_sorted > G0W0.positiveE.txt
            awk '$2  ~ /^-/' tmp_sorted > G0W0.negativeE.txt
            
            cp -r G0W0.positiveE.txt G0W0.negativeE.txt $DIR/graph_data/GW_band_structure/.

	    
        cd ..
        
            mv output G0W0

            cp -r $DIR/band_structure/work/$prefix.DFT_bands.in ./.
  
            mv $prefix.DFT_bands.in $prefix.G0W0.bands.in
    
            sed -i.bak "s|.*GfnQPdb.*|GfnQPdb= \"E < ./output/ndb.QP\"|" $prefix.G0W0.bands.in
        
            $YPP -F "$prefix.G0W0.bands.in" 2> out
    
            mv o.bands_interpolated "$prefix.G0W0.o.bands_interpolated"

            cp -r "$prefix.G0W0.o.bands_interpolated" "$DIR/graph_data/GW_band_structure/."
            rm -f *.bak out l_setup

cd ..

# Step 3: Self-consistent GW iterations (evGW)
for ((i=1; i<=$gw_self_inter; i++)); do
    prev_i=$((i - 1))
    prev_folder="G${prev_i}W${prev_i}"
    curr_folder="G${i}W${i}"

    echo "Starting $curr_folder from $prev_folder..."

    cp -r "$prev_folder" "$curr_folder"

    cd "$curr_folder" || exit 1

    # Prepare the next GW input
    cp $finmane yambo_g${i}w${i}_input.in

    # Update input file to read previous QP corrections
    sed -i.bak "s|.*GfnQPdb.*|GfnQPdb= \"E < ../$prev_folder/$prev_folder/ndb.QP\"|" yambo_g${i}w${i}_input.in
    sed -i.bak "s|.*XfnQPdb.*|XfnQPdb= \"E < ../$prev_folder/$prev_folder/ndb.QP\"|" yambo_g${i}w${i}_input.in

    echo "Running Yambo for $curr_folder..."
    $MPIRUN_YAMBO $YAMBO -F yambo_g${i}w${i}_input.in -J $curr_folder 2> out
    
    cp -r o-$curr_folder.qp $DIR/graph_data/GW_band_structure
    cp -r o-$curr_folder.qp $curr_folder
    
            #cd Report

            # Remove comment lines and extract columns 3 and 4
            awk '!/^#/' o-$curr_folder.qp | awk '{print $3 "\t" $4}' | sort -n > tmp_sorted

            # Split into positive and negative energy files
            awk '$2 !~ /^-/' tmp_sorted > $curr_folder.positiveE.txt
            awk '$2  ~ /^-/' tmp_sorted > $curr_folder.negativeE.txt
            
            cp -r $curr_folder.positiveE.txt $curr_folder.negativeE.txt $DIR/graph_data/GW_band_structure/.
            rm -f tmp_sorted
            
            #cd ..
    
    cp -r $DIR/band_structure/work/$prefix.DFT_bands.in ./.
    
    mv $prefix.DFT_bands.in $prefix.$curr_folder.bands.in
    
    sed -i.bak "s|.*GfnQPdb.*|GfnQPdb= \"E < ./$curr_folder/ndb.QP\"|" $prefix.$curr_folder.bands.in
        
    $YPP -F "$prefix.$curr_folder.bands.in" 2> out
    
    mv o.bands_interpolated "$prefix.$curr_folder.o.bands_interpolated"
    cp -r "$prefix.$curr_folder.o.bands_interpolated" "$DIR/graph_data/GW_band_structure/."
    
    cd ..

done

{
echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Self-consistent GW report: $curr_folder"
echo "    → Initialization source   : Converged G0W0 results."
echo "    → Purpose                 : Finalized parameters for ev-GW loop."
echo "    → ev-GW Data location     : $DIR/graph_data/GW_band_structure."
echo "    → Status                  : ev-GW loop for $curr_folder completed."
echo "End of self-consistent GW report."
echo ""
} >> "$DIR/work_report.txt"
