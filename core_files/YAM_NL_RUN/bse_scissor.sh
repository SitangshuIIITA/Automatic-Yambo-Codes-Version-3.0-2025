#!/bin/bash



cd $DIR

        mkdir -p Real-Time
        
        PATH1=$DIR/GW_folder
        PATH2=$DIR/BSE_folder
        PATH3=$DIR/Real-Time
        PATH4=$DIR/GW_folder/GW_selfC
        
        mkdir -p $DIR/graph_data/Real_Time_data
        mkdir -p $DIR/graph_data/Real_Time_data/BSE_scissor_data

        cd $PATH3
        		rm -rf Real-Time
                          {
                            echo ""
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] BSE-scissor with diagonalization calculation report:"
                           } >> "$DIR/work_report.txt"


        mkdir -p BSE_scissor_slab_z BSE_scissor_box_z
        
        cd BSE_scissor_slab_z
        
        scp -r $PATH1/exch1 ./.
        awk '{ print $(NF-0) }' exch1 > bsengexx
        bseexx="$(cat bsengexx)"

        scp -r $PATH1/exch3 ./.
        awk '{ print $(NF-0) }' exch3 > con_ngs
        con="$(cat con_ngs)"
        
        scp -r $PATH1/exch4 ./.
        echo "$bndrnge" > exch4
        awk '{ print $(NF-0) }' exch4 > dip_bnd
        dip_b="$(cat dip_bnd)"


        
        scp -r $PATH2/"BSE-$con$ryd"/SAVE ./.
              
        scp -r $PATH2/"BSE-$con$ryd"/output ./.
        
              					cd output
                                        	
                                        		rm -rf ndb.BS* ndb.cut* ndb.dip* ndb.em1s* ndb.QP_i* ndb.RIM*
                                        	cd ..
                
        scp -r $PATH2/"BSE-$con$ryd"/Report/r-output_rim_cut_rim_w_optics_dipoles_bss_bse_em1s ./.
                       
                    
        scp -r $PATH2/"BSE-$con$ryd"/r_setup ./.
        
                                            {
                                                        echo "    → Response block size               : NGsBlkXs = $con Ry"
                                                        echo "    → Response block size               : BSENGexx = $bseexx Ry"
                                                        echo "    → Response block size               : BSENGBlk = $con Ry"
                                                        echo "    → Polarization/Screening Bands      : Scheme = 1 | $dip_b"
                                            } >> "$DIR/work_report.txt"
        
#-----------------------------------------------BSE + scissor : diagonalization solver--------------------------------------------#

            rm -rf "$sci_filename"
        
            $YAMBO -r -rw -X s -optics b -kernel sex -Ksolver d -V all -F $sci_filename 2> out


    #-----------------------------------------------COULOMB cut-off & RIM-------------------------------------------------#
                        
                                            sed -i.bak 's/rimw_type= "default"/rimw_type= "'"$rimw_type"'"/g'                               "$sci_filename"
                                            sed -i.bak 's/RandGvecW= 1               RL/RandGvecW= '"$randgvecw $rl"'/g'                    "$sci_filename"
                                            sed -i.bak 's/RandQpts=0/RandQpts='$randqpts'/g'                                                "$sci_filename"
                                            sed -i.bak 's/RandGvec= 1                RL/RandGvec= '"$randgvec $rl"'/g'                      "$sci_filename"
                                            sed -i.bak 's/CUTGeo= "none"/CUTGeo= "'"$slab"'"/g'                                             "$sci_filename"
                                            sed -i.bak 's/DBsIOoff= "none"/DBsIOoff= "BS"/g'                                                "$sci_filename"
                                            sed -i.bak 's/DBsFRAGpm= "none"/DBsFRAGpm= "+BS"/g'                                             "$sci_filename"
    

         #--------------------------------------------CPU parallel structure----------------------------------------------#

                                            sed -i.bak "s/BS_CPU= \"\"/BS_CPU=\"$bs_cpu\"/g"                                                "$sci_filename"
                                            sed -i.bak "s/BS_ROLEs= \"\"/BS_ROLEs= \"k eh t\"/g"                                            "$sci_filename"

                                            sed -i.bak "s/X_and_IO_CPU= \"\"/X_and_IO_CPU=\"$x_and_io\"/g"                                  "$sci_filename"
                                            sed -i.bak "s/X_and_IO_ROLEs= \"\"/X_and_IO_ROLEs= \"q g k c v\"/g"                             "$sci_filename"

                                            sed -i.bak "s/DIP_CPU= \"\"/DIP_CPU=\"$dip_cpu\"/g"                                             "$sci_filename"
                                            sed -i.bak "s/DIP_ROLEs= \"\"/DIP_ROLEs= \"k c v\"/g"                                           "$sci_filename"


#---------------------------------------------FFT-GVecs and other Response size block------------------------------------#

                                            sed '174d' "$sci_filename" > tmp1 && mv tmp1 "$sci_filename"
                                            awk "NR==174{print \"NGsBlkXs=$con $ryd\"}1" "$sci_filename" > tmp && mv tmp                   "$sci_filename"
                                                
                                            sed '88d' "$sci_filename" > tmp1 && mv tmp1 "$sci_filename"
                                            awk "NR==88{print \"BSENGBlk=$con $ryd\"}1" "$sci_filename" > tmp && mv tmp                    "$sci_filename"
 

                                            sed '82d' "$sci_filename" > tmp1 && mv tmp1 "$sci_filename"
                                            awk "NR==82{print \"BSENGexx=$bseexx $ryd\"}1" "$sci_filename" > tmp && mv tmp                 "$sci_filename"
 
 
#---------------------------------------------Includes el-hole coupling--------------------------------------------------#

                                            sed -i.bak 's/BSEmod= "resonant"/BSEmod= "resonant"/g'                                          "$sci_filename"
                                            sed -i.bak 's/Lkind= "BAR"/Lkind= "'"$lkind"'"/g'                                               "$sci_filename"
                                            #sed -i.bak 's/#WehCpl/WehCpl/g'                                                                "$sci_filename"

#---------------------------------------------Includes G0W0 or self_C_GW + gap correction----------------------------------------#

                                                             if [ "$BSE_gap" = "g0w0" ]; then
                                                             
                                    				    #scp -r $PATH1/con_ryd ./.
                                    				    #con_r="$(cat con_ryd)"
                                    	
                                    					#    scp -r "$PATH1/GW-$con_r$ryd/SAVE" ./.
                                    					#    scp -r "$PATH1/GW-$con_r$ryd/output" ./.
                                    					        {
                                    					                    echo "    → QP gap is fixed from              : converged G0W0 = GW-$con$ryd."
                                                        
                                    					        } >> "$DIR/work_report.txt"

                                                             
                                                                    scp -r "$PATH1/GW-$con$ryd/Report/o-output.qp" ./.
                                                                    
                                                                   # Remove comment lines and extract columns 3 and 4
                                                                    awk '!/^#/' o-output.qp | awk '{print $3 "\t" $4}' | sort -n > tmp_sorted

                                                                   # Split into positive and negative energy files
                                                                    awk '$2 !~ /^-/' tmp_sorted > GW_positiveE.txt
                                                                    awk '$2  ~ /^-/' tmp_sorted > GW_negativeE.txt
      
                                                                    scp -r $DIR/lin_reg.py ./.
                                                                    python3 lin_reg.py
                                                                    
                                                                    {
                                                                        echo "    → GW gap is fixed from              : G0W0 result."
                                                                    } >> "$DIR/work_report.txt"
                                                                    


                                                            elif [ "$BSE_gap" = "gw_selfC" ]; then
                                                                    scp -r "$PATH4/G${gw_self_inter}W${gw_self_inter}/o-G${gw_self_inter}W${gw_self_inter}.qp" ./.
                                                                     # Remove comment lines and extract columns 3 and 4
                                                                    awk '!/^#/' o-G${gw_self_inter}W${gw_self_inter}.qp | awk '{print $3 "\t" $4}' | sort -n > tmp_sorted

                                                                   # Split into positive and negative energy files
                                                                    awk '$2 !~ /^-/' tmp_sorted > GW_positiveE.txt
                                                                    awk '$2  ~ /^-/' tmp_sorted > GW_negativeE.txt
      
                                                                    scp -r $DIR/lin_reg.py ./.
                                                                    python3 lin_reg.py
                                                                    {
                                                                        echo "    → GW gap is fixed from              : sc-GW result."
                                                                    } >> "$DIR/work_report.txt"
                                                                  
                                                            fi
                                                                                                                                    
                                                                    awk '/coarse-grid/ {print $6}' r-output_rim_cut_rim_w_optics_dipoles_bss_bse_em1s > 1
                                                                    sed '1d' 1 > 2
                                                                    echo "$(cat 2) |$(cat slope_c.txt) |$(cat slope_v.txt) |" > scissor_block
                                                                    grep -v '^$' scissor_block > scissor_block.clean && mv scissor_block.clean scissor_block

                                                                    awk -v n=96 -v f="scissor_block" 'NR==n {
                                                                        while ((getline x < f) > 0) print x
                                                                            next
                                                                            } {print}' "$sci_filename" > temp_file && mv temp_file "$sci_filename"

                                                                    rm -rf 1 2

                                                    # Normal BSE band selection
                                                            filled_val=$(awk '/Filled Bands/ {print $5}' r_setup)
                                                            empty_val=$(awk '/Empty Bands/ {print $5}' r_setup)

                                                            vb=$((filled_val - $nv))
                                                            cb=$((empty_val + $nc))

                                                            echo "${vb} |${cb} |" > insert_block
                                                            grep -v '^$' insert_block > insert_block.clean && mv insert_block.clean insert_block

                                                            awk -v n=121 -v f="insert_block" 'NR==n {
                                                                while ((getline x < f) > 0) print x
                                                                    next
                                                                } {print}' "$sci_filename" > temp_file && mv temp_file "$sci_filename"
                                                                
                                                                
     
#----------------------------Includes photon energy range, W screening bands, writes exciton wavefunctions---------------#


                                            sed -i.bak '136s/0.100000 | 0.100000 |/'"$Bdmrange"' | '"$Bdmrange"' |/g'                       "$sci_filename"
                                            sed -i.bak '133s/  0.00000 | 10.00000 |/  0.00000 | '"$phot_max_en_lim"' |/g'                   "$sci_filename"

                                            sed -i.bak 's/BEnSteps= 100/BEnSteps= 1000/g'                                                   "$sci_filename"
                                            sed -i.bak 's/#WRbsWF/WRbsWF/g'                                                                 "$sci_filename"
 
                                             # Delete line 172 (works on both systems): Polarization bands
                                            sed '172d' "$sci_filename" > tmp1 && mv tmp1                                                    "$sci_filename"
                                            awk "NR==172{print \"1| $dip_b|\"}1" "$sci_filename" > tmp && mv tmp                            "$sci_filename"
 
                                              # Delete line 72 (works on both systems): Dip bands
                                            sed '72d' "$sci_filename" > tmp1 && mv tmp1                                                     "$sci_filename"
                                            awk "NR==72{print \"1| $dip_b|\"}1" "$sci_filename" > tmp && mv tmp                             "$sci_filename"
 
 
 #---------------------------- Define electric field direction vector based on mode_electric_field-----------------------#

                                                        case "$mode_electric_field" in
                                                            x)
                                                                vec="1.000000 | 0.000000 | 0.000000 |"
                                                                ;;
                                                            y)
                                                                vec="0.000000 | 1.000000 | 0.000000 |"
                                                                ;;
                                                            z)
                                                                vec="0.000000 | 0.000000 | 1.000000 |"
                                                                ;;
                                                            xy|yx)
                                                                vec="1.000000 | 1.000000 | 0.000000 |"
                                                               ;;
                                                            yz|zy)
                                                                vec="0.000000 | 1.000000 | 1.000000 |"
                                                               ;;
                                                            zx|xz)
                                                                vec="1.000000 | 0.000000 | 1.000000 |"
                                                                ;;
                                                            *)
                                                                echo "Error: Unsupported mode_electric_field='$mode_electric_field'"
                                                                exit 1
                                                                ;;
                                                        esac

                                                    # Use sed for in-place editing without the backup extension

                                                        sed -i.bak '141s/1.000000 | 0.000000 | 0.000000 |/'"$vec"'/g'                         "$sci_filename"
                                                        sed -i.bak '187s/1.000000 | 0.000000 | 0.000000 |/'"$vec"'/g'                         "$sci_filename"



#-----------------------------------------------------Running Yambo:----------------------------------------------------#

                                                    $MPIRUN_YAMBO $YAMBO -F "$sci_filename" -J output -C Report 2> out
                                                    
                                                    cd Report
                                                    
                                                            mv o-output.alpha_q1_diago_bse $prefix.scissor.slab_z.o-output.alpha_q1_diago_bse
                                                            

                                                            scp -r $prefix.scissor.slab_z.o-output.alpha_q1_diago_bse $DIR/graph_data/Real_Time_data/BSE_scissor_data

                                                    cd ..
                                                    
                                             
                                         {
                                                        echo "    → Scissor-parameters                : Scheme = $(cat scissor_block)"
                                                        echo "    → BSE Transition Bands              : Scheme = $(cat insert_block)"
                                                        echo "    → See detailed report in $DIR/graph_data/BSE_scissor_data."
                                            } >> "$DIR/work_report.txt"
                                             
                                             
                            cd ../BSE_scissor_box_z
                            
                            	scp -r ../BSE_scissor_slab_z/SAVE ./.
                            	cp -r ../BSE_scissor_slab_z/$sci_filename ./.
                            	scp -r ../BSE_scissor_slab_z/r_setup ./.
                            	
                            				sed -i.bak 's/rim_w/#rim_w/g'                               				"$sci_filename"
                            				sed -i.bak 's/RandGvecW/#RandGvecW/g'                       				"$sci_filename"
                            				sed -i.bak 's/rimw_type/#rimw_type/g'                       				"$sci_filename"
                            				sed -i.bak 's/rimw_type/#rimw_type/g'                       				"$sci_filename"
                            				sed -i.bak 's/CUTGeo= "slab z"/CUTGeo= "box z"/g'             				"$sci_filename"
                    					
                    				
                                        		# Extract the value from r_setup using awk and calculate (value - 1)
                                            		    awk '/Alat factors/ {printf "%.6f|\n", $6 - 1}' r_setup > tmp1

                                       			# Add leading zeros to match "0.000000|0.000000|value|"
                                            		    awk '{print "0.000000|0.000000|" $0}' tmp1 > tmp2

		                                        # Remove line 59 from the original file
                			                    awk 'NR != 59' "$sci_filename" > tmp3

                		                        # Insert the new line at line 59
                		                            awk 'NR==58 {print; system("cat tmp2")} NR!=58' tmp3 > tmp4

                                		        # Replace the original file
                                		            mv tmp4 "$sci_filename"

                                		        # Clean up temporary files
                                		            rm -f tmp1 tmp2 tmp3
                    
                    
#-----------------------------------------------------Running Yambo:----------------------------------------------------#

                                                    $MPIRUN_YAMBO $YAMBO -F "$sci_filename" -J output -C Report 2> out
                                                    
                                                    cd Report
                                                    
                                                            mv o-output.alpha_q1_diago_bse $prefix.box_z.scissor.o-output.alpha_q1_diago_bse
                                                            

                                                            scp -r $prefix.box_z.scissor.o-output.alpha_q1_diago_bse $DIR/graph_data/Real_Time_data/BSE_scissor_data

                                                    cd ..
                                                    
                                        cd ..    
                                             
                                         {
                                                        
                                                        echo "End of BSE-scissor with slab_z cut-off report."
                                                        echo "End of BSE-scissor with box_z cut-off report."
                                                        echo ""
                                                        echo ""
                                            } >> "$DIR/work_report.txt"
                     
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                            
        
exit;
