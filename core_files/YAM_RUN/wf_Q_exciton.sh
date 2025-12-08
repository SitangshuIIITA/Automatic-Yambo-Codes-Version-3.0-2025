#!/bin/bash



cd $DIR

        PATH1=$DIR/GW_folder
       	PATH2=$DIR/BSE_folder/BSE_slepC
        PATH5=$DIR/graph_data/BSE_slepC_data
        
        
        cd $PATH2 || exit 1
              
                           {
                            echo ""
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Excitonic pdf in xsf format at Q!=0 report:"
                           } >> "$DIR/work_report.txt"



#-----------------------------------FInite-Q Excitonic distribution plot in 3D xcrysden------------------------------#

                        for line in $QPT; do
                            
                                $YPP -J output -e w -b "${line}" -F $prefix.finite_Q_${line}.ypp_WF.in 2> out
                                
                                # Step 1: Extract Filled and Empty Bands
                                        awk '/Grid dimensions/ {print $4}' r_setup > 7
                                        awk '{print $1 *2}' 7 > grid
                                        exc_grid="$(cat grid)"
                                        
                                sed -i.bak 's/Format= "g"/Format= "x"/g'                                                                                $prefix.finite_Q_${line}.ypp_WF.in
                                sed -i.bak 's/Direction= "1"/Direction= "123"/g'                                                                        $prefix.finite_Q_${line}.ypp_WF.in
                                sed -i.bak '27s/1 | 1 | 1 |/'$exc_grid' | '$exc_grid' | 1 |/g'                                                          $prefix.finite_Q_${line}.ypp_WF.in
                                sed -i.bak '30s/0.000000 | 0.000000 | 0.000000 |/'$at_coord_x' | '$at_coord_y' | '$at_coord_z' |/g'                     $prefix.finite_Q_${line}.ypp_WF.in
                                sed -i.bak 's/States= "0 - 0"/States= "'"$exciton_states"'"/g'                                                          $prefix.finite_Q_${line}.ypp_WF.in
                                sed -i.bak 's/Degen_Step= 0.010000/Degen_Step= 0.00000/g'                                                               $prefix.finite_Q_${line}.ypp_WF.in

                                
                                $MPIRUN_YAMBO $YPP -F $prefix.finite_Q_$QPT.ypp_WF.in -J output 2> out
                                
                                rm -rf $prefix.finite_Q_${line}.ypp_WF.in.bak
                                
                                mv *.xsf $PATH5/.
 
                                            {
                                                        echo "    → Exciton finite Q used                                : Q = ${line}"
                                                        echo ""
                                            } >> "$DIR/work_report.txt"
   
        


                        done
                                            {
                                                        echo "    → See detailed report in $DIR/graph_data/BSE_slepC_data."
                                                        echo "End of excitonic pdf in xsf format at Q!=0 calculation."
                                                        echo ""
                                            } >> "$DIR/work_report.txt"
   
        

exit;
