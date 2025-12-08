#!/bin/bash





source ../../core_files/EPW_VAR/EPW_variables.sh
source ../../core_files/QE_VAR/QE_variables.sh

cd $DIR 


				    echo "============================================================="
				    echo ">>> Running Ionized and Grain Boundary Impurities Mobility For: Electron (quadrupole) at $ii_imp cm^-2 at $t K:"
				    



PATH2="$DIR/epw_$prefix"
PATH3="$DIR/epw_$prefix/E-mu_vs_ii"
PATH4="$DIR/epw_$prefix/gkkp-k-q-BZ"

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] QE Electron Mobility (Ionized +Grain Boundary) report:"
                                } >> "$DIR/EPW_work_report.txt"

# Create epw directory

cd $PATH2
mkdir -p E-mu_vs_ii

#mkdir -p $PATH2/EPW_plots/e-mu_vs_ii

cd $PATH3

scp -r $PATH4/* ./.
read NK1 NK2 NK3 <<< "$kgrid"
rm -rf restart_0.fmt

# Extract Fermi energy from output and apply shift
fermi_energy=$(grep "Detected VBM, CBM from nscf" "wannier_dis_details.txt" | awk -v shift="$fermi_shift_cond" '{printf "%.6f", $7+shift}')

mv _quadrupole.fmt quadrupole.fmt

for ii in $ii_imp; do

carrier_d=$(awk -v n="$ii" 'BEGIN { printf("%.1E", -n) }' | sed 's/+0*/ /; s/ //g')

cat > ${ii}_quad.in << EOF
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

		  	     vme  = 'dipole'
		  	
			
		  	  fsthick = $fsthick
		  	degaussw  = 0.0
		       scattering = .true.
		 scattering_serta = .true.
			  int_mob = .false.
			  carrier = .true.
			 ncarrier = $ii
		    iterative_bte = .true.
		      epmatkqread = .false.
		      mob_maxiter = 500
		     broyden_beta = 1.0
		          bfieldx = 0.0d0
			  bfieldy = 0.0d0
			  bfieldz = 1.0d-10
		  	
		  	   nstemp = 1
			    temps = $t_i
			system_2d = 'quadrupole'

 			     ii_g = .true.
                        ii_charge = 1.0d0
                             ii_n = $ii
                    ii_scattering = .true.
                          ii_only = .false.
                          
                    gb_scattering = .true.
                          gb_size = $gb_size	!grain size in nm units
                          gb_only = .false.

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
			   iprint = 2 
		    bands_skipped = 'exclude_bands = $exclude_bands'


EOF

# Write orbital projections
for i in "${!orbitals[@]}"; do
    idx=$((i+1))
    echo "proj($idx) = '${orbitals[$i]}'" >> ${ii}_quad.in
done


# --- wdata block ---

# Loop over path segments
for ((i=0; i<${#labels[@]}-1; i++)); do
    echo "wdata($((i+3))) = '${labels[$i]} ${coords[$i]} ${labels[$((i+1))]} ${coords[$((i+1))]}'" >> ${ii}_quad.in
done

last_idx=$(( ${#labels[@]} + 1 ))
echo "wdata($((last_idx))) = 'dis_win_min = $dis_win_min'" >> ${ii}_quad.in
echo "wdata($((last_idx+1))) = 'dis_win_max = $dis_win_max'" >> ${ii}_quad.in
echo "wdata($((last_idx+2))) = 'dis_num_iter = $dis_num_iter'" >> ${ii}_quad.in
echo "/" >> ${ii}_quad.in


				    

		
					mkdir -p II_E-mobility/"e-Mobility_${ii}"
					
					
					    # Run EPW
						    $MPIRUN_EPW -input ${ii}_quad.in > ${ii}_quad.out
						  
						    	rm -rf restart_0.fmt restart_0.fmt.bak

						    		
    									cp -r ${ii}_quad.out ./II_E-mobility/"e-Mobility_${ii}"
    					      				cp -r inv_*.fmt ./II_E-mobility/"e-Mobility_${ii}"
    					      				cp -r m_effective.fmt ./II_E-mobility/"e-Mobility_${ii}"
    					      				cp -r mobility_*.fmt ./II_E-mobility/"e-Mobility_${ii}"
    					      				cp -r IBTE*_*.fmt ./II_E-mobility/"e-Mobility_${ii}"
    					      				
    					      				cp -r "II_E-mobility" $PATH2/EPW_plots
    					      		
    					      		mv  work/$prefix.epmatwp $ii.epmatwp

    					      				
   							 echo ">>> Finished processing $ii cm^-2."
							 echo ""
done

				    
				    echo ">>> Ending Ionized and Grain Boundary Impurities For Electron Mobility."
				    echo "============================================================="

					
cd $DIR
		
                         {  echo "    → Electron Mobility at ionized Impurities+Grain Boundaries for dipole and dipole+quadrupole corrections."
                            echo "    → Plottable folder is inside files $DIR/epw_{"$prefix"}/EPW_plots"
                            echo "    → Before final production, cross check the convergences."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of Electronic mobility report."
			    echo ""
                           } >> "$DIR/EPW_work_report.txt"
		
		

exit 0














