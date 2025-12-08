#!/bin/bash





PATH1="$DIR/epw_$prefix/"


cd $DIR

						echo "====================================================="
						echo ">>> Parsing EPW Quadrupole Velocity, Energies and Relaxation Time for Zero-Field Monte-Carlo Simulations."

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] Monte-Carlo report:"
                                } >> "$DIR/EPW_work_report.txt"



# Create epw directory
cd $PATH1
mkdir -p Einstien_e-MC
cd Einstien_e-MC


cp -r ../../../core_files/MC_RUN/epw_e_mc_parser.py ./.
cp -r ../../../core_files/MC_RUN/mc_e_driver.py ./.
cp -r ../../../core_files/MC_RUN/plot_e_mu.py ./.


cp -r ../EPW_plots/e-Mobility_vs_T/e-Mobility_quadrupole/IBTEvel_* ./.
cp -r ../EPW_plots/e-Mobility_vs_T/e-Mobility_quadrupole/m_effective.fmt ./.
cp -r ../EPW_plots/e-Mobility_vs_T/e-Mobility_quadrupole/inv_taucb_0.fmt ./.
cp -r ../w90_calc/scf_w90/scf.out ./.
cp -r ../EPW_plots/mobility_summary_electron.txt ./.


mv IBTEvel_sup_0.fmt electron_IBTEvel_sup_0.fmt
mv m_effective.fmt electron_m_effective.fmt
mv inv_taucb_0.fmt electron_inv_taucb_0.fmt

				
				python3 epw_e_mc_parser.py
				
				
					echo ">>> Ending Parsing of Data Calculations. Check MC_Electron Folder for HDF5 parsed file."
					echo "====================================================="
					echo ""
					echo ">>> Running Monte-Carlo Calculations for zero field electronic mobility [Einstien Mobility]:."
					
					
				python3 mc_e_driver.py	

				python plot_e_mu.py
						
				echo "====================================================="


				rm -rf IBTE* *.fmt *.py	 





                         {  echo "    → See MC_Electrons/ log and MC_OUTPUT folders in $DIR/epw_{"$prefix"} for results."
                            echo "    → ZeroField Calculation on mobilities are done."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of Einstein mobility report."
			    echo ""
                           } >> "../../../EPW_work_report.txt"

exit 0
