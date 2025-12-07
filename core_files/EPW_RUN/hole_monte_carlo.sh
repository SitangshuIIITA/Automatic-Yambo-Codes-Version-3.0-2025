#!/bin/bash





PATH1="$DIR/epw_$prefix/"


cd $DIR

						echo "====================================================="
						echo ">>> Parsing EPW Quadrupole Velocity, Energies and Relaxation Time for Monte-Carlo Simulations."

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] Monte-Carlo report:"
                                } >> "$DIR/EPW_work_report.txt"



# Create epw directory
cd $PATH1

cp -r ../../core_files/EPW_RUN/MC_Holes ./.
cd MC_Holes



cp -r ../EPW_plots/h-Mobility_vs_T/h-Mobility_quadrupole/IBTEvel_* ./.
cp -r ../EPW_plots/h-Mobility_vs_T/h-Mobility_quadrupole/m_effective.fmt ./.
cp -r ../EPW_plots/h-Mobility_vs_T/h-Mobility_quadrupole/inv_tau_0.fmt ./.
cp -r ../w90_calc/scf_w90/scf.out ./.
cp -r ../EPW_plots/mobility_summary_hole.txt ./.


mv IBTEvel_sup_0.fmt hole_IBTEvel_sup_0.fmt
mv m_effective.fmt hole_m_effective.fmt
mv inv_tau_0.fmt hole_inv_tau_0.fmt

				
				python3 epw_mc_parser.py
				
				
					echo ">>> Ending Parsing of Data Calculations. Check MC_NEW Folder for HDF5 parsed file."
					echo "====================================================="
					echo ""
					echo ">>> Running Monte-Carlo Calculations for zero field mobility [Einstien Mobility]:."
					
					
				python3 mc_driver.py	

				python plot_mu.py
						
				echo "====================================================="












                         {  echo "    → Wannierization done."
                            echo "    → Grid Used               : = $kgrid"
                            echo "    → See W90_plots folder in $DIR/epw_{"$prefix"}."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of Wannier report."
			    echo ""
                           } >> "../../../EPW_work_report.txt"

exit 0
