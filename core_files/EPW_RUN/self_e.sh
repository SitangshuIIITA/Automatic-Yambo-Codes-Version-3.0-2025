#!/bin/bash





source ../../core_files/EPW_VAR/EPW_variables.sh
source ../../core_files/QE_VAR/QE_variables.sh

cd $DIR 

Temp="77 100 300 500"

					echo "====================================================="
				    	echo ">>> Running Electronic CBM Self-Energy at $Temp K."
			

PATH2="$DIR/epw_$prefix"
PATH3="$DIR/epw_$prefix/EPW_elec_self"
PATH4="$DIR/epw_$prefix/gkkp-k-q-BZ"

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] Electronic Self Energy at 300 K report:"
                                } >> "$DIR/EPW_work_report.txt"

                               

# Create epw directory

cd $PATH2
mkdir -p EPW_elec_self

cd $PATH3

scp -r $PATH4/* ./.

scp -r ../../../core_files/EPW_RUN/prepare_q.py ./.
python3 prepare_q.py --prepare-q

#cp -r $PATH2/Quadrupole/_quadrupole.fmt ./.

#mv quadrupole.fmt _quadrupole.fmt

read NK1 NK2 NK3 <<< "$kgrid"
rm -rf restart_0.fmt

# Extract Fermi energy from output and apply shift
fermi_energy=$(grep "Detected VBM, CBM from nscf" "wannier_dis_details.txt" | awk -v shift="$fermi_shift_cond" '{printf "%.6f", $7+shift}')


cat > eself_no_lr.in << EOF
&inputepw
		  	   prefix = '$prefix',
		  	   outdir = 'work'
		  	dvscf_dir = 'save'
		  	  etf_mem = 1
		       iverbosity = 3
		       
		  	     elph = .true.
		  	 epbwrite = .false.
		  	  epbread = .false.
		  	 epwwrite = .true.
		  	  epwread = .false.
		  	   lpolar = .false.
		  	  
		  	  ! use_ws = .true.
		  	   !  lifc = .false.
			

		  	  fsthick = $fsthick
		  	degaussw  = $degaussw
		  	
		  	!band_plot = .true.
			system_2d = 'no'
	
		      epmatkqread = .false.
	
			   nstemp = 4
			    temps = $Temp
			
			  restart = .true.
		       selecqread = .false.
		     restart_step = 1000
		      efermi_read = .true.
		     fermi_energy = $fermi_energy
 	                 !scissor = 0 !$scissor

		  	   
 		       elecselfen = .true.

		          
	        	    filkf = 'q.dat'
	  	   
		  	    nqf1  = $nqf1
		  	    nqf2  = $nqf2
		  	    nqf3  = $nqf3


		  	     nk1  = $NK1
		  	     nk2  = $NK2
		  	     nk3  = $NK3
		  	  
		  	     nq1  = $NQ1
		  	     nq2  = $NQ2
		  	     nq3  = $NQ3


  	   
		  ! Converged Wannier inputs:
		       wannierize = .false.
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
    echo "proj($idx) = '${orbitals[$i]}'" >> eself_no_lr.in
done


# --- wdata block ---

# Loop over path segments
for ((i=0; i<${#labels[@]}-1; i++)); do
    echo "wdata($((i+3))) = '${labels[$i]} ${coords[$i]} ${labels[$((i+1))]} ${coords[$((i+1))]}'" >> eself_no_lr.in
done

last_idx=$(( ${#labels[@]} + 1 ))
echo "wdata($((last_idx))) = 'dis_win_min = $dis_win_min'" >> eself_no_lr.in
echo "wdata($((last_idx+1))) = 'dis_win_max = $dis_win_max'" >> eself_no_lr.in
echo "wdata($((last_idx+2))) = 'dis_num_iter = $dis_num_iter'" >> eself_no_lr.in
echo "/" >> eself_no_lr.in



scp -r eself_no_lr.in eself_3d.in
sed 's/lpolar = .false./lpolar = .true./g' 			eself_3d.in > tmp && mv tmp eself_3d.in

scp -r eself_3d.in eself_gaussian.in
sed 's/system_2d = '"'no'"'/system_2d = '"'gaussian'"'/g' 	eself_gaussian.in > tmp && mv tmp eself_gaussian.in


scp -r eself_gaussian.in eself_dipole_sp.in
sed 's/system_2d = '"'gaussian'"'/system_2d = '"'dipole_sp'"'/g' eself_dipole_sp.in > tmp && mv tmp eself_dipole_sp.in


scp -r eself_dipole_sp.in eself_quadrupole.in
sed 's/system_2d = '"'dipole_sp'"'/system_2d = '"'quadrupole'"'/g' 				eself_quadrupole.in > tmp && mv tmp eself_quadrupole.in
awk -v var="		             vme  = '"dipole"' " '/lpolar/{print;print var;next}1' 	eself_quadrupole.in > tmp && mv tmp eself_quadrupole.in


	for lr in 'no_lr' '3d' 'gaussian' 'dipole_sp' 'quadrupole'
		do
		
    			test $lr = 'quadrupole' && mv _quadrupole.fmt quadrupole.fmt 

    				$MPIRUN_EPW -input eself_$lr.in > eself_$lr.out
    
    					for line in $Temp; do
    						mv linewidth.elself.${line}.000K $lr.linewidth.elself.${line}.000K
					done

    			test $lr = 'quadrupole' && mv quadrupole.fmt _quadrupole.fmt
			rm -rf restart_0.fmt
    			mv  work/$prefix.epmatwp $lr.epmatwp
    
	done    
    	


					mkdir -p "$PATH2/EPW_plots/Elself_LW"
						cd "$PATH2/EPW_plots/Elself_LW" || exit 1

							for line in $Temp; do
    								mkdir -p "Elself_LW-${line}"
    									cd "Elself_LW-${line}" || exit 1

    										for lr in 'no_lr' '3d' 'gaussian' 'dipole_sp' 'quadrupole'; do
        									cp -r "$PATH3/${lr}.linewidth.elself.${line}.000K" .
        									cp -r "$PATH3/eself_$lr.out" .
        									grep 'ik =' eself_$lr.out > klist.dat
        									cp -r klist.dat ../../.
    										done

    									cd ..
							done



					
cd $DIR
				    	echo ">>> Ending Electronic CBM Self-Energy at $Temp K. Check EPW_plots Folder."
					echo "====================================================="

		
                         {  echo "    → Electronic Self Energy at 300 K for each phonon mode along q-BZ is done."
                            echo "    → Essestial files generated: *eig, decay.* *.mmn, *.nnkp, etc."
                            echo "    → Before moving ahead, check the relevent electron and phonon band structure files: *band.dat and phband.freq."
                            echo "    → See EPW_plots folder in $DIR/epw_{"$prefix"}."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of EPW_priliminary report."
			    echo ""
                           } >> "$DIR/EPW_work_report.txt"
		
		

exit 0


















