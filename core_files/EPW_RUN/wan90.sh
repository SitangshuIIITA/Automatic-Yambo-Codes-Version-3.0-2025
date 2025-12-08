#!/bin/bash





PATH1="$DIR/relax"
PATH2="$DIR/epw_$prefix/w90_calc/nscf_w90"


cd $DIR

						echo "====================================================="
						echo ">>> Running W90 Calculations:."

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] QE EPW report:"
                                } >> "$DIR/EPW_work_report.txt"



# Create epw directory
cd $PATH2
scp -r ../../../bands/Bandx.dat.gnu ./.

rm -f UNK* epw_$prefix* pw2*


#  Writing a Wannier90 input script:

cat > $prefix.win << EOF
		  use_ws_distance = .true.
			num_bands = $num_bands
			 num_wann = $num_wann
			 num_iter = $num_iter
			   iprint = $iprint 
		  num_dump_cycles = $num_dump_cycles
                 num_print_cycles = $num_print_cycles
                     dis_num_iter = $dis_num_iter

!! To plot the WFs
			! restart = plot
		     wannier_plot = true
	   wannier_plot_supercell = 3
	      ! wannier_plot_list = 1,5
	      		  
	      
	      
!! To plot the WF interpolated bandstructure

		       bands_plot = true
			  mp_grid = $kgrid
!! Disentanglement windows: energies in eV, relative to QE Fermi level			  
		      dis_win_min = $dis_win_min
		      dis_win_max = $dis_win_max
		     dis_froz_min = $dis_froz_min
		     dis_froz_max = $dis_froz_max		  
EOF
scp -r $DIR/$prefix.band_route ./.
# Insert band path
echo "begin kpoint_path" >> $prefix.win
labels=($bands_path)

coords=()
while IFS="|" read -r x y z _; do
    x=$(echo $x | xargs)  # trim spaces
    y=$(echo $y | xargs)
    z=$(echo $z | xargs)
    coords+=("$x $y $z")
done < "$prefix.band_route"

# Loop over segments (pairs of points)
for ((i=0; i<${#labels[@]}-1; i++)); do
    echo "${labels[$i]} ${coords[$i]} ${labels[$((i+1))]} ${coords[$((i+1))]}" >> $prefix.win
done

echo "end kpoint_path" >> $prefix.win

echo "begin projections" >> $prefix.win


source ../../../../core_files/EPW_VAR/EPW_variables.sh
for orb in "${orbitals[@]}"; do
    echo "  $orb" >> $prefix.win
done    
echo "end projections" >> $prefix.win

echo "begin atoms_frac" >> $prefix.win
source ../../../../core_files/QE_VAR/QE_variables.sh
for pos in "${positions[@]}"; do
    echo "  $pos" >> $prefix.win
done
echo "end atoms_frac" >> $prefix.win

echo "begin unit_cell_cart" >> $prefix.win
echo bohr >> $prefix.win
# Extract final atomic coordinates
awk '/Begin final coordinates/,/End final coordinates/' relax.out \
    | awk 'NR > 3 { print }' \
    | awk 'NR > 2 { print last } { last = $0 }' > 1

awk '/CELL_PARAMETERS/,/ATOMIC_POSITIONS/' 1  > 2
sed '1d;$d' "2" > tmp && mv tmp "2"
sed '$d' "2" > tmp && mv tmp "2"

# Convert alat to Bohr
alat=$(awk -F"=" '/CELL_PARAMETERS/ {gsub(/\)/,"",$2); print $2; exit}' 1)
alat_bohr=$(echo "$alat" | bc -l)

# Multiply all numbers in file 2 by alat_bohr
awk -v a="$alat_bohr" '{printf "%12.6f %12.6f %12.6f\n", $1*a, $2*a, $3*a}' 2 > 2_bohr
echo "$(cat 2_bohr)" >> "$prefix.win"
echo "end unit_cell_cart" >> $prefix.win
echo "begin kpoints" >> $prefix.win
echo "$(cat epw_k_list.in)" >> "$prefix.win"
echo "end kpoints" >> $prefix.win
rm -rf 1 2 2_bohr


		$MPIRUN_wan -pp $prefix

#  Writing a Wannier90 input script:

cat > pw2wannier90.in << EOF
&inputpp
   outdir = '.'
   prefix = '$prefix'
   seedname = '$prefix'
   write_amn = .true.
   write_mmn = .true.
   write_unk = .true.
/
EOF

		$MPIRUN_pw2wan < pw2wannier90.in > pw2wannier90.out

		$MPIRUN_wan $prefix

scp -r *.xsf "$prefix"_band.dat $prefix.wout ../../W90_plots/.


					echo ">>> Ending W90 Calculations. Check W90_plots Folder."
					echo "====================================================="
					

                         {  echo "    → Wannierization done."
                            echo "    → Grid Used               : = $kgrid"
                            echo "    → See W90_plots folder in $DIR/epw_{"$prefix"}."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of Wannier report."
			    echo ""
                           } >> "../../../EPW_work_report.txt"

exit 0
