#!/bin/bash





source ../../core_files/EPW_VAR/EPW_variables.sh
source ../../core_files/QE_VAR/QE_variables.sh

cd $DIR

					echo "====================================================="
				    	echo ">>> Running & Plotting Electronic Self-Energy."
			

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] Plotting of CB Self-Energy report:"
                                } >> "$DIR/EPW_work_report.txt"


PATH4="$DIR/epw_$prefix/EPW_plots"

cd "$PATH4" || exit 1

	scp -r "../../../core_files/EPW_RUN/plot_self_e_KGM.py" .
	INPUT="klist.dat"
	OUTPUT="klist_clean.dat"

		if [[ ! -f "$INPUT" ]]; then
    			echo "Error: $INPUT not found!"
    			exit 1
		fi

		# Use awk to keep only first occurrence of ik
			awk '
				/coord.:/ {
  					  ik = $3
  						  if (!seen[ik]++) {
 						       print $0
 						   }
						}
						' "$INPUT" > "$OUTPUT"

		#echo "Cleaned klist written to $OUTPUT"
		mv "$OUTPUT" "$INPUT"

		scp -r "$DIR/epw_$prefix/EPW_prelim/wannier_dis_details.txt" .

python3 plot_self_e_KGM.py

#scp -r *.png $DIR/epw_${prefix}/EPW_plots/.

scp -r "../../../core_files/EPW_RUN/plot_e_quadrupole_mode_bar.py" .
python3 plot_e_quadrupole_mode_bar.py


scp -r "../../../core_files/EPW_RUN/plot_e_scattering.py" .
python3 plot_e_scattering.py




				    	echo ">>> Ending & Plotting Electronic Self-Energy. Check EPW_plots/self_energy Folders."
					echo "====================================================="



   			 {  echo "    → Electronic Self energies are plotted."
                            echo "    → Plottable folder is inside files $DIR/epw_"{$prefix}"/EPW_plots"
                            echo "    → Before final production, cross check the convergences."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of Self energy report."
			    echo ""
                           } >> "$DIR/EPW_work_report.txt"		




























exit 0
