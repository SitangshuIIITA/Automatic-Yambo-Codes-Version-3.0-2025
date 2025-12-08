#!/bin/bash



cd $DIR

                      {
                            echo ""
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Excitonic oscillator strength at \Gamma [Q=0] report:"
                           } >> "$DIR/work_report.txt"
     

        PATH1=$DIR/GW_folder
       	PATH2=$DIR/BSE_folder
        PATH5=$DIR/graph_data/BSE_data_files
        
        
        scp -r $PATH1/exch3 $PATH2/.
        cd $PATH2 || exit 1
        awk '{ print $(NF-0) }' exch3 > con_ngs
        con="$(cat con_ngs)"
        
        cd "BSE-$con$ryd"
        
        
#			1. Run first the BSE-GW script: bse-gw.sh


 				#-------------Excitonic Oscillator Strength & Amplitude Calculation--------------------------#

               					$YPP -J output -e a 1 2> out
               				
                                mv o-output.exc_qpt1_amplitude* o-output.exc_qpt1_weights* $PATH5/.

                                                {
                                                        echo "    → See detailed report in $DIR/graph_data/BSE_data_files."
                                                        echo "End of BSE [Q=0] oscillator report."
                                                        echo ""
                                                } >> "$DIR/work_report.txt"


exit;
