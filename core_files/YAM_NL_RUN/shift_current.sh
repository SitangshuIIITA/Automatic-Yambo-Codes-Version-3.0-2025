#!/bin/bash


        cd $DIR
        
        PATH0=$DIR/Real-Time/BSE_scissor_slab_z
        PATH1=$DIR/GW_folder
        PATH2=$DIR/Real-Time/HXC_coll
        PATH3=$DIR/Real-Time

                          {
                            echo ""
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] SHift-current report:"
                           } >> "$DIR/work_report.txt"

        cd $PATH3
        
            mkdir -p Shift
        
        cd Shift

        scp -r $PATH0/insert_block ./.
        scp -r $PATH0/scissor_block ./.
        scp -r $PATH0/bsengexx ./.
                bseexx="$(cat bsengexx)"

        scp -r $PATH0/exch3 ./.
        awk '{ print $(NF-0) }' exch3 > con_ngs
        con="$(cat con_ngs)"

        
                            scp -r $PATH2/SAVE ./.
                            scp -r $PATH2/r_setup ./.
a=5
b=2

min_shg_NLEnRange=$(printf "%.4f" "$(echo "scale=4; $min_NLEnRange * $a" | bc)")
max_shg_NLEnRange=$(printf "%.4f" "$(echo "scale=4; $phot_max_en_lim / $b" | bc)")



                                         $YAMBO_NL -r -u n -V all -F $prefix.shift.in 2> out

    #-----------------------------------------------COULOMB cut-off & RIM-------------------------------------------------#


                                        # Extract the value from r_setup using awk and calculate (value - 1)
                                            awk '/Alat factors/ {printf "%.6f|\n", $6 - 1}' r_setup > tmp1

                                        # Add leading zeros to match "0.000000|0.000000|value|"
                                            awk '{print "0.000000|0.000000|" $0}' tmp1 > tmp2

                                        # Remove line 45 from the original file
                                            awk 'NR != 45' "$prefix.shift.in" > tmp3

                                        # Insert the new line at line 45
                                            awk 'NR==44 {print; system("cat tmp2")} NR!=44' tmp3 > tmp4

                                        # Replace the original file
                                            mv tmp4 "$prefix.shift.in"

                                        # Clean up temporary files
                                            rm -f tmp1 tmp2 tmp3

                                            sed -i.bak 's/RandQpts=0/RandQpts='$randqpts'/g'                                            	"$prefix.shift.in"
                                            sed -i.bak 's/RandGvec= 1                RL/RandGvec= '"$randgvec $rl"'/g'                      	"$prefix.shift.in"
                                            sed -i.bak 's/CUTGeo= "none"/CUTGeo= "'"$slab_nl"'"/g'                                             	"$prefix.shift.in"
                                            sed -i.bak 's/DBsIOoff= "none"/DBsIOoff= "BS"/g'                                                	"$prefix.shift.in"
                                            sed -i.bak 's/DBsFRAGpm= "none"/DBsFRAGpm= "+BS"/g'                                             	"$prefix.shift.in"
    

         #--------------------------------------------CPU parallel structure----------------------------------------------#

                                            sed -i.bak "s/NL_CPU= \"\"/NL_CPU=\"$nl_cpu\"/g"                         	                    	"$prefix.shift.in"
                                            sed -i.bak "s/NL_ROLEs= \"\"/NL_ROLEs= \"w k\"/g"                                               	"$prefix.shift.in"

                                            sed -i.bak "s/OSCLL_CPU= \"\"/OSCLL_CPU=\"$oscll_cpus\"/g"                                   	"$prefix.shift.in"
                                            sed -i.bak "s/OSCLL_ROLEs= \"\"/OSCLL_ROLEs= \"k b\"/g"                                         	"$prefix.shift.in"

                                            sed -i.bak "s/DIP_CPU= \"\"/DIP_CPU=\"$dip_cpu\"/g"                                             	"$prefix.shift.in"
                                            sed -i.bak "s/DIP_ROLEs= \"\"/DIP_ROLEs= \"k c v\"/g"                                           	"$prefix.shift.in"




#---------------------------------------------Includes G0W0 or self_C_GW + gap correction----------------------------------------#

                                                             if [ "$BSE_gap" = "g0w0" ]; then
                                                                    
                                                                    {
                                                                        echo "    → GW gap is fixed from              : G0W0 result [scissor corrected]."
                                                                    } >> "$DIR/work_report.txt"
                                                                    


                                                            elif [ "$BSE_gap" = "gw_selfC" ]; then
                                                                    {
                                                                        echo "    → GW gap is fixed from              : sc-GW result [scissor corrected]."
                                                                    } >> "$DIR/work_report.txt"
                                                                  
                                                            fi
                                                                                                                                    
                                                                    awk -v n=86 -v f="scissor_block" 'NR==n {
                                                                        while ((getline x < f) > 0) print x
                                                                            next
                                                                            } {print}' "$prefix.shift.in" > temp_file && mv temp_file "$prefix.shift.in"

                                                                    rm -rf 1 2

                                                            awk -v n=52 -v f="insert_block" 'NR==n {
                                                                while ((getline x < f) > 0) print x
                                                                    next
                                                                } {print}' "$prefix.shift.in" > temp_file && mv temp_file "$prefix.shift.in"


#---------------------------------------------Nonlinear simulation parameters-------------------------------------------#

                                            sed -i.bak 's/NLtime=-1.000000/NLtime='"$NLtime"'/g'                                           	"$prefix.shift.in"

                                            sed -i.bak 's/NLintegrator= "INVINT"/NLintegrator= "'"$nl_integrator"'"/g'                     	"$prefix.shift.in"

                                            sed -i.bak 's/NLCorrelation= "IPA"/NLCorrelation= "'"$nl_Correlation"'"/g'                     	"$prefix.shift.in"

                                            sed -i.bak 's/NLDamping= 0.00000/NLDamping= '"$nl_Damping"'/g'                            		"$prefix.shift.in"
                                            
                                            sed -i.bak 's/NLEnSteps=0/NLEnSteps= '"$nl_ensteps"'/g'                                 		"$prefix.shift.in"
                                            
                                            sed -i.bak 's/#EvalCurrent/EvalCurrent/g'		                                 		"$prefix.shift.in"
                                            
                                            sed -i.bak "62s/-1.000000 |-1.000000 |/${min_shg_NLEnRange}|${max_shg_NLEnRange}|/g"               	"$prefix.shift.in"

                                            
#---------------------------------------------Exchange and correlation cut-offs----------------------------------------#

                                            sed '78d' "$prefix.shift.in" > tmp1 && mv tmp1 "$prefix.shift.in"
                                            awk "NR==78{print \"HARRLvcs= $bseexx $ryd\"}1" "$prefix.shift.in" > tmp && mv tmp            	"$prefix.shift.in"
                                                
                                            sed '79d' "$prefix.shift.in" > tmp1 && mv tmp1 "$prefix.shift.in"
                                            awk "NR==79{print \"EXXRLvcs= $bseexx $ryd\"}1" "$prefix.shift.in" > tmp && mv tmp            	"$prefix.shift.in"
 

                                            sed '80d' "$prefix.shift.in" > tmp1 && mv tmp1 "$prefix.shift.in"
                                            awk "NR==80{print \"CORRLvcs= $con $ryd\"}1" "$prefix.shift.in" > tmp && mv tmp            		"$prefix.shift.in"
 
 
#---------------------------------------------Choice of field----------------------------------------#

                                            sed -i.bak 's/"SOFTSIN"/"'"SOFTSIN"'"/g'                                                         	"$prefix.shift.in"

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

                                                        sed -i.bak '110s/0.000000 | 0.000000 | 0.000000 |/'"$vec"'/g'                     	"$prefix.shift.in"
                                                   
                                                   
                                        $MPIRUN_YAMBO $YAMBO_NL -F $prefix.shift.in 2> out



    
                                            {
                                                        echo "End of Shift-current calculations."
                                                        echo ""
                                            } >> "$DIR/work_report.txt"
                                            
                                            
                                            


exit;
