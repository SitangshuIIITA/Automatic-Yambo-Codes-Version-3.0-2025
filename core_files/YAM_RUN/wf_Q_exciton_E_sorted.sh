#!/bin/bash



cd $DIR

        PATH1=$DIR/GW_folder
       	PATH2=$DIR/BSE_folder/BSE_slepC
        PATH5=$DIR/graph_data/BSE_slepC_data
        
        
        cd $PATH2 || exit 1
              
                           {
                            echo ""
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] BSE finite-Q post-processing exciton E-sorted weight report:"
                           } >> "$DIR/work_report.txt"



#-----------------------------------FInite-Q Excitonic distribution plot in 3D xcrysden------------------------------#

                        for line in $QPT; do
                            
                                $YPP -J output -e s -b "${line}"  2> out
                                
                                 mv *amplitude* *weights* $PATH5/.
 
                                            {
                                                        echo "    → Exciton finite Q used                                : Q = ${line}"
                                                        echo ""
                                            } >> "$DIR/work_report.txt"
   
        


                        done
                                            {
                                                        echo "    → See detailed report in $DIR/graph_data/BSE_slepC_data."
                                                        echo "End of exciton E-sorted weight at finite-Q calculation."
                                                        echo ""
                                            } >> "$DIR/work_report.txt"
   
        

exit;
