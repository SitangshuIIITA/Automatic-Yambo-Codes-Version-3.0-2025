#!/bin/bash



cd $DIR

PATH3="$DIR/nscf"
PATH5="$DIR/band_structure"
PATH6="$DIR/GW_folder"
PATH7="$DIR/BSE_GW"
PATH8="$DIR/nscf-dg"
PATH14="$DIR/BSE_folder/BSE_Temperature"
PATH15="$DIR/graph_data/BSE_Temperature_Data"

cd "$PATH14"

                          {
                            echo ""
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Exciton intensities with temperature calculation report:"
                           } >> "$DIR/work_report.txt"


#--------- Exciton Intensity sorting: all positive energies and intensities between & including 0.3 and 1 --------#
for line in $temp_value; do
    cd "BSE-${line}$Kelvin" || continue

    awk '(NR>1) && ($1 > 0) && ($2 > 0.3 )' "$prefix.${line}$Kelvin.o-output.exc_qpt1_E_sorted" > tmp_energy_sorted.txt

    #--------- List of exciton states which are degenerate by 0.01 eV --------#
    awk 'NR>1 && ($1 - prev[1]) <= 0.01 { printf "%d - %d\n", prev[3], $3 } { split($0, prev) }' tmp_energy_sorted.txt > tmp_degenerate.txt

    awk '
        {
            first[NR]=$1;
            last[NR]=$NF;
            row[NR]=$0
        }
        END {
            for (i = 1; i <= NR; i++)
                if (last[i] != first[i+1])
                    print row[i]
        }
    ' tmp_degenerate.txt > exc_states.txt

    mv tmp_energy_sorted.txt "$prefix.${line}$Kelvin.energy_state_linewidth.txt"

    sed -i.bak '1 i\
#    E [ev]             Strength           Index              W [meV]' "$prefix.${line}$Kelvin.energy_state_linewidth.txt"

    sed -i.bak '1 i\
# Exciton Intensity sorting: all positive energies and intensities between & including 0.3 and 1' "$prefix.${line}$Kelvin.energy_state_linewidth.txt"

    sed -i.bak '1 i\
# List of exciton states which are degenerate by 0.01 eV' exc_states.txt

    cp exc_states.txt exc_states_list.txt

    sed 1d exc_states.txt

    rm -f tmp_degenerate.txt *.bak

    while IFS= read -r states; do
        "$YPP_PH" -e -e -F "exciton_eliash.$line$Kelvin.${states}" 2> out

        sed -i '' 's/States= "0 - 0"/States= "'"$states"'"/g' "exciton_eliash.$line$Kelvin.${states}"
        sed -i '' 's/PhStps= 200/PhStps= 1000/g' "exciton_eliash.$line$Kelvin.${states}"
        sed -i '' 's/PhBroad= 0.010000/PhBroad= 0.0050000/g' "exciton_eliash.$line$Kelvin.${states}"

        "$YPP_PH" -F "exciton_eliash.$line$Kelvin.${states}" -J output 2> out

        cd "$PATH14/BSE-${line}$Kelvin" || continue

        cp "$prefix.${line}$Kelvin.energy_state_linewidth.txt" exc_states_list.txt o-output.g_sq_F* "$PATH15/"

    done < exc_states.txt

    cd ..
done

                                         {
                                                        echo "    → Exciton Intensity sorting: all positive energies and intensities between & including 0.3 and 1."
                                                        echo "    → Exciton states which are degenerate less than 0.01 eV are also listed."
                                                        echo "    → See detailed report in $DIR/graph_data/BSE_Temperature_Data."
                                                        echo "End of exciton intensities with temperature calculation."
                                                        echo ""
                                            } >> "$DIR/work_report.txt"






exit
