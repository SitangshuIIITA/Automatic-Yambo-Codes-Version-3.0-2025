#!/bin/bash


source ../../core_files/EPW_VAR/EPW_variables.sh
source ../../core_files/QE_VAR/QE_variables.sh

cd $DIR 

PATH2="$DIR/epw_$prefix"
PATH3="$DIR/epw_$prefix/E-scattering"
PATH4="$DIR/epw_$prefix/gkkp-E-Matrix"

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] QE Electronic Mobility report:"
                                } >> "$DIR/EPW_work_report.txt"

# Create epw directory

cd $PATH2
mkdir E-scattering

mkdir $PATH2/EPW_plots/e-scattering

cd $PATH3

scp -r $PATH4/* ./.
#rm -rf *.epb*  *.epmatwp *.bvec *cube  decay* *amn *chk *eig *mmn *nnkp UNK* pw2* *.xsf *.win $prefix.wout *.dat *band.gnu *band.kpt *.k* *.png *.freq* epw_* dyn* q2r* EPW* matdyn* *.py *.ifc *.in *.ukk 

read NK1 NK2 NK3 <<< "$kgrid"
rm -rf restart_0.fmt
scp -r dipole_sp.epmatwp ./work/.
mv ./work/dipole_sp.epmatwp ./work/$prefix.epmatwp


# Extract Fermi energy from output and apply shift
fermi_energy=$(grep "Detected VBM, CBM from nscf" "wannier_dis_details.txt" | awk -v shift="$fermi_shift_cond" '{printf "%.6f", $7+shift}')

cat > esca_dipole_sp.in << EOF
&inputepw
		  	   prefix = '$prefix',
		  	   outdir = 'work'
		  	dvscf_dir = 'save'
		  	  etf_mem = 1
		       iverbosity = 3
		       
		  	     elph = .true.
		  	 epbwrite = .false.
		  	  epbread = .false.
		  	 epwwrite = .false.
		  	  epwread = .true.
		  	   lpolar = .true.
		  	  
		  	   use_ws = .true.
		  	     lifc = .false.
		  	mp_mesh_k = .true.
			

		  	  fsthick = $fsthick
		  	degaussw  = 0.0
		  	!elecselfen = .true.
		       scattering = .true.
		 scattering_serta = .true.
			  int_mob = .false.
			  carrier = .true.
			 ncarrier = $ncarrier
		    !iterative_bte = .true.
		      epmatkqread = .false.
		     ! mob_maxiter = 500
		     !broyden_beta = 1.0
		      !    bfieldx = 0.0d0
			!  bfieldy = 0.0d0
			 ! bfieldz = 1.0d-10
		  	
		  	   nstemp = $nstemp
			    temps = $temps
			system_2d = 'dipole_sp'

			  restart = .true.
		       selecqread = .false.
		     restart_step = 50
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
    echo "proj($idx) = '${orbitals[$i]}'" >> esca_dipole_sp.in
done


# --- wdata block ---

# Loop over path segments
for ((i=0; i<${#labels[@]}-1; i++)); do
    echo "wdata($((i+3))) = '${labels[$i]} ${coords[$i]} ${labels[$((i+1))]} ${coords[$((i+1))]}'" >> esca_dipole_sp.in
done

last_idx=$(( ${#labels[@]} + 1 ))
echo "wdata($((last_idx))) = 'dis_win_min = $dis_win_min'" >> esca_dipole_sp.in
echo "wdata($((last_idx+1))) = 'dis_win_max = $dis_win_max'" >> esca_dipole_sp.in
echo "wdata($((last_idx+2))) = 'dis_num_iter = $dis_num_iter'" >> esca_dipole_sp.in
echo "/" >> esca_dipole_sp.in



awk -v var="		             vme  = '"wannier"' " '/lpolar/{print;print var;next}1' 	esca_dipole_sp.in > tmp && mv tmp esca_dipole_sp.in


scp -r esca_dipole_sp.in esca_quadrupole.in
sed 's/system_2d = '"'dipole_sp'"'/system_2d = '"'quadrupole'"'/g' 				esca_quadrupole.in > tmp && mv tmp esca_quadrupole.in





			for lr in 'dipole_sp' 'quadrupole'
				do
				    echo "==============================="
				    echo ">>> Running EPW for: electron-$lr"
				    echo "==============================="


					    # For quadrupole, temporarily rename the format file
    						test $lr = 'quadrupole' && mv _quadrupole.fmt quadrupole.fmt

					    # Run EPW
						    $MPIRUN_EPW -input esca_$lr.in > esca_$lr.out

								mkdir -p "e-scattering_$lr"
    					
    					      				scp -r esca_$lr.out ./"e-scattering_$lr"
    					      				scp -r inv_*.fmt ./"e-scattering_$lr"
    					      				scp -r m_effective.fmt ./"e-scattering_$lr"
    					      				scp -r mobility_*.fmt ./"e-scattering_$lr"
    					      				scp -r IBTE*_*.fmt ./"e-scattering_$lr"
    					      
    					      				scp -r "e-scattering_$lr" $PATH2/EPW_plots/e-scattering/.


						    rm -f restart_0.fmt	

						    # Restore quadrupole format file
						    test $lr = 'quadrupole' && mv quadrupole.fmt _quadrupole.fmt

							    # Move the epmatwp file if it exists
								    if [ -f work/$prefix.epmatwp ]; then
								        #mv work/$prefix.epmatwp $lr.epmatwp
								        echo "Saved $lr.epmatwp successfully."
								    else
								        echo "Error: work/$prefix.epmatwp not found! Skipping $lr..."
								    fi

   							 echo ">>> Finished processing $lr"
							    echo ""
				done



					
cd $DIR
		
                         {  echo "    → Electronic Mobility at "$temps" for dipole and dipole+quadrupole corrections."
                            echo "    → Plottable folder is inside files $DIR/epw_{"$prefix"}/EPW_plots"
                            echo "    → Before final production, cross check the convergences."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of Electronic mobility report."
			    echo ""
                           } >> "./EPW_work_report.txt"
		
		

exit 0


















