#!/bin/bash





source ../../core_files/EPW_VAR/EPW_variables.sh
source ../../core_files/QE_VAR/QE_variables.sh

cd $DIR

					echo "=========================================================================="
				    	echo ">>> Running & Plotting no_lr and quadrupole electronic gkkp [intra and interband] along ph-BZ route KGM for k = CBM."
			

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] Quadrupole report:"
                                } >> "$DIR/EPW_work_report.txt"


PATH4="$DIR/epw_$prefix/gkkp-E-Matrix"
scp -r ../core_files/EPW_RUN/comp.py $PATH4/.

cd "$PATH4" || exit 1

python3 comp.py
scp -r *.png $DIR/epw_${prefix}/EPW_plots/.


					echo ">>> Ending  Plots. Check EPW_plots Folder."
					echo "=========================================================================="
			

         		 {  echo "    → This computes and plots the no_lr and quadrupole electronic gkkp [intra and interband] along ph-BZ route KGM for k = CBM. See EPW_variable list to use desired CBM and VBM ith and jth band numbers."
                            echo "    → Plottable folder is inside files $DIR/epw_"{$prefix}"/EPW_plots"
                            echo "    → Before final production, cross check the convergences."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of Matrix elements are plotted."
			    echo ""
                           } >> "$DIR/EPW_work_report.txt"




























exit 0
