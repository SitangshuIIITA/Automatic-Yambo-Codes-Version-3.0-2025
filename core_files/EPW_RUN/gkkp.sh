#!/bin/bash





source ../../core_files/EPW_VAR/EPW_variables.sh
source ../../core_files/QE_VAR/QE_variables.sh

cd $DIR

					echo "============================================================"
				    	echo ">>> Running gkkp Matrix Elements Of Electron + Holes For: no_lr, 3d, gaussian, dipole_sp & quadrupole Conditions:"
				    	echo ">>> Note: The Matrix Elements Will Be Printed At allk and q of KGM route."
				    	

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] Quadrupole report:"
                                } >> "$DIR/EPW_work_report.txt"


mkdir -p epw_$prefix/gkkp-k-q-BZ


PATH1="$DIR/epw_$prefix/w90_calc/nscf_w90"
PATH2="$PATH1/relax_w90"
PATH3="$DIR/epw_$prefix/Quadrupole"
PATH4="$DIR/epw_$prefix/gkkp-k-q-BZ"
PATH5="$DIR/epw_$prefix/various_lr"


cd "$PATH4" || exit 1


scp -r $PATH5/* ./. 

scp -r $PATH1/wannier_dis_details.txt ./.


##################### get the CBM k #####################################
nscf="nscf.out"
outfile="k.txt"
cartfile="cart_coords.txt"
crystfile="cryst_coords.txt"

# --- 1. Get VBM energy ---
e_vbm=$(awk '/highest occupied, lowest unoccupied level/ {print $7}' "$nscf")

# --- 2. Extract Cartesian coordinates (2pi/alat) ---
awk '/^[[:space:]]*cart\. coord/ { in_cart=1; next }
     in_cart && /^[[:space:]]*k\(/ { print }' "$nscf" > "$cartfile"

# --- 3. Extract crystal coordinates ---
awk '/^[[:space:]]*cryst\. coord/ { in_cryst=1; next }
     in_cryst && /^[[:space:]]*k\(/ { print }' "$nscf" > "$crystfile"

# --- 4. Get fractional coordinates of CBM k-point ---
read kx ky kz <<< $(awk -v e="$e_vbm" '
/k =/ { last_kx=$3; last_ky=$4; last_kz=$5 }
/'"$e_vbm"'/ { print last_kx, last_ky, last_kz; exit }
' "$nscf")

# --- 5. Find closest k-point in NSCF and get its index ---
closest_line=$(awk -v kx="$kx" -v ky="$ky" -v kz="$kz" '
BEGIN { min_dist=1e10; closest_line="" }
/^[[:space:]]*k\(/ {
    gsub(/[(),]/,"",$5); gsub(/[(),]/,"",$6); gsub(/[(),]/,"",$7)
    cx=$5+0; cy=$6+0; cz=$7+0
    dx=cx - kx; dy=cy - ky; dz=cz - kz
    dist=dx*dx + dy*dy + dz*dz
    if (dist < min_dist) {
        min_dist=dist
        closest_line=$0
    }
}
END { print closest_line }
' "$cartfile")

# --- 6. Extract k-index from closest_line ---
kidx=$(echo "$closest_line" | sed -n 's/.*k(\s*\([0-9]\+\)\s*).*/\1/p')

# --- 7. Get matching line from crystal coords file ---
cryst_line=$(awk -v idx="$kidx" '
/^[[:space:]]*k\(/ {
    if (match($0,/k\([[:space:]]*([0-9]+)[[:space:]]*\)/,arr)) {
        if (arr[1]==idx) print
    }
}
' "$crystfile")

# --- 8. Extract only coords + weight ---
coords=$(echo "$cryst_line" | sed -E 's/.*=\s*\(([^)]*)\), wk =\s*([0-9.Ee+-]+)/\1 \2/')

# --- 9. Write clean output ---
{
  echo "1 crystal"
  echo "$coords"
} > "$outfile"

echo "CBM k-point written to $outfile"
read NK1 NK2 NK3 <<< "$kgrid"

cat > epw_no_lr.in << EOF
&inputepw
		  	   prefix = '$prefix',
		  	   outdir = 'work'
		  	dvscf_dir = 'save'
		  	  etf_mem = 1
		       
		  	     elph = .true.
		  	 epbwrite = .true.
		  	  epbread = .false.
		  	 epwwrite = .true.
		  	  epwread = .false.
		  	   lpolar = .false.
		  	  
		  	   use_ws = .true.
		  	     lifc = .false.
			
		  	degaussw  = $degaussw
		  	

		  	   prtgkk = .true.

		  	band_plot = .true.
			system_2d = 'no'
			  

			  restart = .true.
		       selecqread = .false.
		     restart_step = 50
		      efermi_read = .true.
		     !fermi_energy = $fermi_energy
		  	
			    filqf = 'q.dat'

			    filkf = 'q.dat'
			
		  	     nk1  = $NK1
		  	     nk2  = $NK2
		  	     nk3  = $NK3
		  	
		  	     nq1  = $NQ1
		  	     nq2  = $NQ2
		  	     nq3  = $NQ3
  	   
		  ! Converged Wannier inputs:
		       wannierize = .true.
		          nbndsub = $num_wann
		         num_iter = $num_iter
		     dis_froz_min = $dis_froz_min
		     dis_froz_max = $dis_froz_max
			   iprint = $iprint 
		    bands_skipped = 'exclude_bands = $exclude_bands'


EOF

# Write orbital projections
for i in "${!orbitals[@]}"; do
    idx=$((i+1))
    echo "proj($idx) = '${orbitals[$i]}'" >> epw_no_lr.in
done


# --- wdata block ---

# Loop over path segments
for ((i=0; i<${#labels[@]}-1; i++)); do
    echo "wdata($((i+3))) = '${labels[$i]} ${coords[$i]} ${labels[$((i+1))]} ${coords[$((i+1))]}'" >> epw_no_lr.in
done

last_idx=$(( ${#labels[@]} + 1 ))
echo "wdata($((last_idx))) = 'dis_win_min = $dis_win_min'" >> epw_no_lr.in
echo "wdata($((last_idx+1))) = 'dis_win_max = $dis_win_max'" >> epw_no_lr.in
echo "wdata($((last_idx+2))) = 'dis_num_iter = $dis_num_iter'" >> epw_no_lr.in
echo "/" >> epw_no_lr.in



scp -r epw_no_lr.in epw_3d.in
sed 's/wannierize = .true./wannierize = .false./g' 		epw_3d.in > tmp && mv tmp epw_3d.in
sed 's/epbread = .false./epbread = .true./g' 			epw_3d.in > tmp && mv tmp epw_3d.in
sed 's/epbwrite = .true./epbwrite = .false./g' 			epw_3d.in > tmp && mv tmp epw_3d.in
sed 's/lpolar = .false./lpolar = .true./g' 			epw_3d.in > tmp && mv tmp epw_3d.in


scp -r epw_3d.in epw_gaussian.in
sed 's/wannierize = .true./wannierize = .false./g' 		epw_gaussian.in > tmp && mv tmp epw_gaussian.in
sed 's/lpolar = .false./lpolar = .true./g' 			epw_gaussian.in > tmp && mv tmp epw_gaussian.in
sed 's/system_2d = '"'no'"'/system_2d = '"'gaussian'"'/g' 	epw_gaussian.in > tmp && mv tmp epw_gaussian.in


scp -r epw_gaussian.in epw_dipole_sp.in
sed 's/system_2d = '"'gaussian'"'/system_2d = '"'dipole_sp'"'/g' epw_dipole_sp.in > tmp && mv tmp epw_dipole_sp.in


scp -r epw_dipole_sp.in epw_quadrupole.in
sed 's/system_2d = '"'dipole_sp'"'/system_2d = '"'quadrupole'"'/g' 				epw_quadrupole.in > tmp && mv tmp epw_quadrupole.in
awk -v var="		             vme  = '"wannier"' " '/lpolar/{print;print var;next}1' 	epw_quadrupole.in > tmp && mv tmp epw_quadrupole.in

scp -r $PATH3/quadrupole.fmt ./.
mv quadrupole.fmt _quadrupole.fmt

for lr in 'no_lr' '3d' 'gaussian' 'dipole_sp' 'quadrupole'
do
    test $lr = 'quadrupole' && mv _quadrupole.fmt quadrupole.fmt

    		$MPIRUN_EPW -input epw_$lr.in > epw_$lr.out

    test $lr = 'quadrupole' && mv quadrupole.fmt _quadrupole.fmt

    mv work/$prefix.epmatwp $lr.epmatwp
done


				    	echo ">>> Ending gkkp Matrix Elements for: Electron + Holes."
				    	echo "============================================================"

cd $DIR
		
                         {  echo "    → Electronic gkkp matrix elements are evaluated."
                            echo "    → Plottable folder is inside files $DIR/epw_{"$prefix"}/EPW_plots"
                            echo "    → Before final production, cross check the convergences."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of gkkp report."
			    echo ""
                           } >> "$DIR/EPW_work_report.txt"
		
		



exit 0
