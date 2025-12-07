#!/bin/bash





PATH1="$DIR/epw_$prefix/MC_Electron"


cd $PATH1

						echo "====================================================="
						echo ">>> Parsing EPW Quadrupole Velocity, Energies and Relaxation Time for Low --> High Field Ensemble Monte-Carlo Simulations."

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] Monte-Carlo report:"
                                } >> "$DIR/EPW_work_report.txt"


cp -r ../../../core_files/EPW_RUN/MC_Electron/drift_parser.py ./.
cp -r ../../../core_files/EPW_RUN/MC_Electron/analyze_abs_em.py ./.
cp -r ../../../core_files/EPW_RUN/MC_Electron/analyze_BZ_Scattering.py ./.
cp -r ../../../core_files/EPW_RUN/MC_Electron/compute_scattering_stats.py ./.
cp -r ../../../core_files/EPW_RUN/MC_Electron/update_velocities_from_gradient.py ./.
cp -r ../../../core_files/EPW_RUN/MC_Electron/mc-1_new.py ./.
cp -r ../../../core_files/EPW_RUN/MC_Electron/cell_area.py ./.
cp -r ../../../core_files/EPW_RUN/MC_Electron/analyze_BZ_Scattering.py ./.
cp -r ../../../core_files/EPW_RUN/MC_Electron/map_final_pop_to_epw_grid.py ./.
cp -r ../../../core_files/EPW_RUN/MC_Electron/mc-2.py ./.
cp -r ../../../core_files/EPW_RUN/MC_Electron/plot_pop.py ./.
cp -r ../../../core_files/EPW_RUN/MC_Electron/plot_drift.py ./.
cp -r ../gkkp-k-q-BZ/epw_quadrupole.out ./.
cp -r ../gkkp-k-q-BZ/q.dat ./.
cp -r ../w90_calc/scf_w90/scf.out ./.


					#python3 drift_parser.py \
  					#	--quad epw_quadrupole.out \
    					#	--qdat q.dat \
    					#	--out epw_quadrupole_parsed.h5


					echo ">>> Ending Parsing of Data Calculations. Check MC_Electron Folder for HDF5 parsed file."
					echo "====================================================="
					echo ""
					echo ""
					echo ""
					echo ">>> Running Monte-Carlo Calculations for field dependent mobility [Ensemble MonteCarlo]:."
					echo "====================================================="
					
				
				
				#python3 compute_scattering_stats.py 
			
				#python3 update_velocities_from_gradient.py

	
				#mpirun -np 4 python3 mc-2.py
				
				#python3 plot_drift.py
				
				#python3 map_final_pop_to_epw_grid.py
				
				#python3 plot_pop.py
				
				
				

				#rm -rf drift_parser.py single_file_mc_fast.py


				


				
						
				echo "====================================================="
                         {  echo "    → Drift velocity vs Electrinc Field (along x direction) is done."
                            echo "    → See MC_Electrons folder in $DIR/epw_{"$prefix"} for results."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of Drift-velcoty MC report."
			    echo ""
                           } >> "../../../EPW_work_report.txt"

exit 0
