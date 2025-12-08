#!/bin/bash





source ../../core_files/EPW_VAR/EPW_variables.sh
source ../../core_files/QE_VAR/QE_variables.sh

cd $DIR 


					echo "=========================================================================="
				    	echo ">>> Running mobility for: electrons at $ncarrier cm^-2 and $temps K:"
				    	


PATH1="$DIR/epw_$prefix"

PATH2="$PATH1/gkkp-k-q-BZ"
PATH3="$PATH1/E-Mobility_vs_T"


                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] QE Electronic Mobility report:"
                                } >> "$DIR/EPW_work_report.txt"

# Create epw directory

mkdir -p $PATH1/EPW_plots/e-Mobility_vs_T



mkdir -p $PATH1/E-Mobility_vs_T


cd $PATH3

cp -r $PATH2/* ./.

read NK1 NK2 NK3 <<< "$kgrid"
rm -rf restart_0.fmt

# Extract Fermi energy from output and apply shift
fermi_energy=$(grep "Detected VBM, CBM from nscf" "wannier_dis_details.txt" | awk -v shift="$fermi_shift_cond" '{printf "%.6f", $7+shift}')

cat > emob_dipole_sp.in << EOF
&inputepw
		  	   prefix = '$prefix',
		  	   outdir = 'work'
		  	dvscf_dir = 'save'
		  	  etf_mem = 3
		       iverbosity = 3
		       
		  	     elph = .true.
		  	 epbwrite = .false.
		  	  epbread = .false.
		  	 epwwrite = .true.
		  	  epwread = .false.
		  	   lpolar = .true.
		  	mp_mesh_k = .true.
			

		  	  fsthick = $fsthick
		  	degaussw  = 0.0
		       scattering = .true.
		 scattering_serta = .true.
			  int_mob = .false.
			  carrier = .true.
			 ncarrier = $ncarrier
		    iterative_bte = .true.
		      epmatkqread = .false.
		      mob_maxiter = 500
		     broyden_beta = 1.0
		          bfieldx = 0.0d0
			  bfieldy = 0.0d0
			  bfieldz = 1.0d-10
		  	
		  	   nstemp = $nstemp
			    temps = $temps
			system_2d = 'dipole_sp'

			  restart = .true.
		       selecqread = .false.
		     restart_step = 1000
		      efermi_read = .true.
		     fermi_energy = $fermi_energy
 	                 !scissor = 0 !$scissor

		  	   
 		            nkf1  = $nkf1
		  	    nkf2  = $nkf2
		  	    nkf3  = $nkf3

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
    echo "proj($idx) = '${orbitals[$i]}'" >> emob_dipole_sp.in
done


# --- wdata block ---

# Loop over path segments
for ((i=0; i<${#labels[@]}-1; i++)); do
    echo "wdata($((i+3))) = '${labels[$i]} ${coords[$i]} ${labels[$((i+1))]} ${coords[$((i+1))]}'" >> emob_dipole_sp.in
done

last_idx=$(( ${#labels[@]} + 1 ))
echo "wdata($((last_idx))) = 'dis_win_min = $dis_win_min'" >> emob_dipole_sp.in
echo "wdata($((last_idx+1))) = 'dis_win_max = $dis_win_max'" >> emob_dipole_sp.in
echo "wdata($((last_idx+2))) = 'dis_num_iter = $dis_num_iter'" >> emob_dipole_sp.in
echo "/" >> emob_dipole_sp.in



awk -v var="		             vme  = '"dipole"' " '/lpolar/{print;print var;next}1' 	emob_dipole_sp.in > tmp && mv tmp emob_dipole_sp.in


cp -r emob_dipole_sp.in emob_quadrupole.in
sed 's/system_2d = '"'dipole_sp'"'/system_2d = '"'quadrupole'"'/g' 				emob_quadrupole.in > tmp && mv tmp emob_quadrupole.in


cp -r emob_quadrupole.in spec_decom_quadrupole.in

sed 's/nstemp = '"$nstemp"'/nstemp = 1/g'			 			spec_decom_quadrupole.in > tmp && mv tmp spec_decom_quadrupole.in
sed 's/temps = '"$temps"'/temps = 300/g'			 			spec_decom_quadrupole.in > tmp && mv tmp spec_decom_quadrupole.in
awk -v var="		     mob_maxfreq  = 2500 " '/lpolar/{print;print var;next}1' 	spec_decom_quadrupole.in > tmp && mv tmp spec_decom_quadrupole.in
awk -v var="		       mob_nfreq  = 300 " '/lpolar/{print;print var;next}1' 	spec_decom_quadrupole.in > tmp && mv tmp spec_decom_quadrupole.in
    	



			for lr in 'dipole_sp' 'quadrupole'
				do
				    	
				    echo ">>> Running mobility for: electron-$lr."
				    


					    # For quadrupole, temporarily rename the format file
    						test $lr = 'quadrupole' && mv _quadrupole.fmt quadrupole.fmt

							
					    # Run EPW
						    $MPIRUN_EPW -input emob_$lr.in > emob_$lr.out
						
								mkdir -p "e-Mobility_$lr"
    					
    					      				cp -r emob_$lr.out ./"e-Mobility_$lr"
    					      				cp -r inv_*.fmt ./"e-Mobility_$lr"
    					      				cp -r m_effective.fmt ./"e-Mobility_$lr"
    					      				cp -r mobility_*.fmt ./"e-Mobility_$lr"
    					      				cp -r IBTE*_*.fmt ./"e-Mobility_$lr"
    					      
    					      				cp -r "e-Mobility_$lr" $PATH1/EPW_plots/e-Mobility_vs_T/.

						    # Restore quadrupole format file
						    test $lr = 'quadrupole' && mv quadrupole.fmt _quadrupole.fmt
						  
						    rm -rf restart_0.fmt
    						    mv  work/$prefix.epmatwp $lr.epmatwp

   							 echo ">>> Finished processing $lr"
							 echo ""
				done

 					
    	
    	mv _quadrupole.fmt quadrupole.fmt
    						
    						
    					
				    			echo ">>> Running spectral decomposition for: electron-quadrupole at 300 K."

								rm -rf restart_0.fmt
    						    		
								$MPIRUN_EPW -input spec_decom_quadrupole.in > spec_decom_quadrupole.out

							echo ">>> Finished processing spectral decomposition at 300 K."


	
	mv quadrupole.fmt _quadrupole.fmt
	rm -rf restart_0.fmt restart_0.fmt.bak 					
 					
 					
 					
 					
 					
 					
cd $DIR

					echo ">>> Ending mobility for: electrons."
				    	echo "=========================================================================="
				
				
				
		
                         {  echo "    → Electronic Mobility at "$temps" for dipole and dipole+quadrupole corrections."
                            echo "    → Plottable folder is inside files $DIR/epw_{"$prefix"}/EPW_plots"
                            echo "    → Before final production, cross check the convergences."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of Electronic mobility report."
			    echo ""
                           } >> "$DIR/EPW_work_report.txt"
		
		


exit 0


















