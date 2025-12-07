#!/bin/bash





source ../../core_files/EPW_VAR/EPW_variables.sh
source ../../core_files/QE_VAR/QE_variables.sh

cd $DIR 

					echo "===================================================="
				    	echo ">>> Running Mobility (E+H) Extraction Process at "$temps" and Impurities for dipole and dipole+quadrupole corrections:"
				    	

PATH2="$DIR/epw_$prefix/EPW_plots"
PATH3="$DIR/epw_$prefix/E-Mobility_vs_T"
PATH4="$DIR/epw_$prefix/gkkp-E-Matrix"

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] Mobility extractions:"
                                } >> "$DIR/EPW_work_report.txt"

# Create epw directory

cd $PATH2

	scp -r ../../../core_files/EPW_RUN/mobility_extract.py .
	
	scp -r ../../../core_files/EPW_RUN/mobility_vs_T_Plot.py .
	
	scp -r ../../../core_files/EPW_RUN/gaussian-h.py .
	
	#scp -r ../../../core_files/EPW_RUN/scattering_rate.py .
	
	scp -r ../../../core_files/EPW_RUN/mobility_extract_ii.py .
	
	
	python3 mobility_extract.py
rm -rf mobility_extract.py	

	python3 gaussian-h.py
rm -rf gaussian-h.py

	python3 mobility_vs_T_Plot.py
rm -rf mobility_vs_T_Plot.py


	python3 mobility_extract_ii.py
rm -rf mobility_vs_T_Plot.py
	

	
	
				    	echo ">>> Ending Mobility (E+H) Extraction Process. Check Folder: EPW_plots."
				    	echo "===================================================="
	

         {  echo "    → Electronic & Hole Mobility at "$temps" and Impurities for dipole and dipole+quadrupole corrections."
                            echo "    → Plottable folder is inside files $DIR/epw_"{$prefix}"/EPW_plots"
                            echo "    → Before final production, cross check the convergences."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of Mobility extraction report."
			    echo ""
                           } >> "$DIR/EPW_work_report.txt"
		
		

exit 0


















