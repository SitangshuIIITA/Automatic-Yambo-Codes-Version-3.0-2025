#!/bin/bash



cd $DIR

        PATH1=$DIR/GW_folder
       	PATH2=$DIR/BSE_folder
        PATH5=$DIR/graph_data/BSE_data_files
        
                      {
                            echo ""
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Excitonic pdf in xsf format at \Gamma [Q=0] report:"
                           } >> "$DIR/work_report.txt"
     

        cd $PATH2 || exit 1
        con="$(cat con_ngs)"
        
        cd "BSE-$con$ryd"
        
        

				#-----------------------------------Excitonic distribution plot in 3D xcrysden------------------------------#

                                $YPP -F $prefix.ypp_WF.in -J output  -e w  1 2> out
                                
                                # Step 1: Extract Filled and Empty Bands
                                        awk '/Grid dimensions/ {print $4}' r_setup > 7
                                        awk '{print $1 *2}' 7 > grid
                                        exc_grid="$(cat grid)"
                                        
                                sed -i.bak 's/Format= "g"/Format= "x"/g'                                                                                $prefix.ypp_WF.in
                                sed -i.bak 's/Direction= "1"/Direction= "123"/g'                                                                        $prefix.ypp_WF.in
                                sed -i.bak '26s/1 | 1 | 1 |/'$exc_grid' | '$exc_grid' | 1 |/g'                                                          $prefix.ypp_WF.in
                                sed -i.bak '29s/0.000000 | 0.000000 | 0.000000 |/'$at_coord_x' | '$at_coord_y' | '$at_coord_z' |/g'                     $prefix.ypp_WF.in
                                sed -i.bak 's/States= "0 - 0"/States= "'"$exciton_states"'"/g'                                                          $prefix.ypp_WF.in
                                sed -i.bak 's/Degen_Step= 0.010000/Degen_Step= 0.00000/g'                                                               $prefix.ypp_WF.in

               					$MPIRUN_YAMBO $YPP -F $prefix.ypp_WF.in -J output 2> out
                                mv *.xsf $PATH5/.
 
                                                 {
                                                        echo "    → See detailed report in $DIR/graph_data/BSE_data_files."
                                                        echo "End of BSE [Q=0] oscillator report."
                                                        echo ""
                                                } >> "$DIR/work_report.txt"

 
exit;
