#!/bin/bash





source ../../core_files/EPW_VAR/EPW_variables.sh
source ../../core_files/QE_VAR/QE_variables.sh

cd $DIR


					echo "=========================================================================="
				    	echo ">>> Running EPW Wannierization for CB and VB:"
			

PATH1="$DIR/relax"
PATH2="$DIR/epw_$prefix/w90_calc"
PATH3="$DIR/epw_$prefix/EPW_prelim"
PATH4="$DIR/epw_$prefix/EPW_calc/scf_nonsoc"
PATH5="$DIR/epw_$prefix/EPW_calc/nscf"

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] EPW-priliminary report:"
                                } >> "$DIR/EPW_work_report.txt"


mkdir -p $DIR/epw_$prefix/EPW_prelim
mkdir -p $DIR/epw_$prefix/EPW_prelim/work

mkdir -p $DIR/epw_$prefix/EPW_plots

# Create epw directory
cd $PATH3


scp -r $PATH2/nscf_w90/$prefix.save ./work/.
scp -r $PATH2/nscf_w90/$prefix.xml ./work/.
scp -r $PATH2/nscf_w90/relax.out ./.

scp -r $DIR/$prefix.band_route ./.

rm -rf *amn *chk *eig *mmn *nnkp UNK* pw2* *.xsf *.win $prefix.wout *.dat *band.gnu *band.kpt 

scp -r $PATH2/w90_phonon_disp/save ./.
scp -r $PATH2/w90_phonon_disp/_ph0 ./work/.

scp -r $PATH2/w90_phonon_disp/filkf.txt ./.
scp -r $PATH2/nscf_w90/wannier_dis_details.txt ./.

read NK1 NK2 NK3 <<< "$kgrid"

cat > $prefix.in << EOF
&inputepw
		  	   prefix = '$prefix',
		  	   outdir = 'work'
		  	dvscf_dir = 'save'
		  	  etf_mem = 0
			
		  	     elph = .true.
		  	 epbwrite = .true.
		  	  epbread = .false.
		  	 epwwrite = .true.
		  	  epwread = .false.
		  	
		  	   lpolar = .true.
		  	      vme = 'dipole'
		  	  fsthick = 100
		  	!system_2d = 'quadrupole'
		  	
		  	degaussw  = 0.001
		  	band_plot = .true.
		      efermi_read = .true.
		          scissor = 0.0! $scissor		  	
		  	
		  	
		  	    filkf = './filkf.txt'
		  	    filqf = './filkf.txt'
		  	
		  	     nk1  = $NK1
		  	     nk2  = $NK2
		  	     nk3  = $NK3
		  	
		  	     nq1  = $NQ1
		  	     nq2  = $NQ2
		  	     nq3  = $NQ3
  	   
		  ! Converged Wannier inputs:
		       wannierize = .true.
		     wannier_plot = .true.
		          nbndsub = $num_wann
		         num_iter = $num_iter
		     dis_froz_min = $dis_froz_min
		     dis_froz_max = $dis_froz_max
			   iprint = 2 
		    bands_skipped = 'exclude_bands = $exclude_bands'


EOF

# Write orbital projections
for i in "${!orbitals[@]}"; do
    idx=$((i+1))
    echo "proj($idx) = '${orbitals[$i]}'" >> $prefix.in
done


# --- wdata block ---
echo "wdata(1) = 'bands_plot = .true.'" >> ${prefix}.in
echo "wdata(2) = 'begin kpoint_path'" >> ${prefix}.in

# --- generate k-point path automatically ---
labels=($bands_path)

coords=()
while IFS="|" read -r x y z _; do
    x=$(echo $x | xargs)  # trim spaces
    y=$(echo $y | xargs)
    z=$(echo $z | xargs)
    coords+=("$x $y $z")
done < "$prefix.band_route"

# Loop over path segments
for ((i=0; i<${#labels[@]}-1; i++)); do
    echo "wdata($((i+3))) = '${labels[$i]} ${coords[$i]} ${labels[$((i+1))]} ${coords[$((i+1))]}'" >> ${prefix}.in
done

last_idx=$(( ${#labels[@]} + 2 ))
echo "wdata(${last_idx}) = 'end kpoint_path'" >> ${prefix}.in
echo "wdata($((last_idx+1))) = 'dis_win_min = $dis_win_min'" >> ${prefix}.in
echo "wdata($((last_idx+2))) = 'dis_win_max = $dis_win_max'" >> ${prefix}.in
echo "wdata($((last_idx+3))) = 'num_dump_cycles = $num_dump_cycles'" >> ${prefix}.in
echo "wdata($((last_idx+4))) = 'num_print_cycles = $num_print_cycles'" >> ${prefix}.in
echo "wdata($((last_idx+5))) = 'guiding_centres = .true.'" >> ${prefix}.in
echo "wdata($((last_idx+6))) = 'dis_num_iter = $dis_num_iter'" >> ${prefix}.in
echo "wdata($((last_idx+7))) = 'use_ws_distance = T'" >> ${prefix}.in
echo "wdata($((last_idx+8))) = 'conv_tol = $conv_tol'" >> ${prefix}.in
echo "wdata($((last_idx+9))) = 'bands_plot_format = gnuplot'" >> ${prefix}.in

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

echo "wdata($((last_idx+10))) = 'begin unit_cell_cart'" >> ${prefix}.in
echo "wdata($((last_idx+11))) = 'bohr'" >> ${prefix}.in

# Print 3 lattice vectors from 2_bohr
lineno=0
while read -r x y z; do
    lineno=$((lineno+1))
    echo "wdata($((last_idx+11+lineno))) = '$x     $y     $z'" >> ${prefix}.in
done < 2_bohr

echo "wdata($((last_idx+15))) = 'end unit_cell_cart'" >> ${prefix}.in
echo "/" >> ${prefix}.in

		$MPIRUN_EPW -input ${prefix}.in > ${prefix}.out 


scp -r $DIR/bands/Bandx.dat.gnu ../EPW_plots/.
scp -r $prefix.out "$prefix"_band.dat phband.freq decay* ../EPW_plots/.
scp -r $DIR/phonon_disp/$prefix.freq.gp ../EPW_plots/.


						    	echo ">>> Ending EPW Wannierization For CB and VB. Before Moving Ahead, Check Relevent Electron and Phonon band Structure Files: *band.dat and phband.freq. Check EPW_plots Folder For Plottable Data."
							echo "=========================================================================="			


                         {  echo "    → EPW Bandstructures (ELectron & Phonon) are done."
                            echo "    → Essestial files generated: *eig, decay.* *.mmn, *.nnkp, etc."
                            echo "    → Before moving ahead, check the relevent electron and phonon band structure files: *band.dat and phband.freq."
                            echo "    → See EPW_plots folder in $DIR/epw_{"$prefix"}."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of EPW_priliminary report."
			    echo ""
                           } >> "$DIR/EPW_work_report.txt"

exit 0
