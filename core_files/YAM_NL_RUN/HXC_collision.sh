#!/bin/bash


        cd $DIR
        
        PATH0=$DIR/Real-Time/BSE_scissor_slab_z
        PATH1=$DIR/GW_folder
        PATH2=$DIR/Real-Time/NL_setup/FixSymm
        PATH3=$DIR/Real-Time

                          {
                            echo ""
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] HX collision report:"
                           } >> "$DIR/work_report.txt"

        cd $PATH3
        
            mkdir -p HXC_coll
        
        cd HXC_coll

        scp -r $PATH0/insert_block ./.
        scp -r $PATH0/scissor_block ./.
        
        scp -r $PATH0/exch1 ./.
        awk '{ print $(NF-0) }' exch1 > bsengexx
        bseexx="$(cat bsengexx)"

        scp -r $PATH0/exch3 ./.
        awk '{ print $(NF-0) }' exch3 > con_ngs
        con="$(cat con_ngs)"
        
        scp -r $PATH0/exch4 ./.
        echo "$bndrnge" > exch4
        awk '{ print $(NF-0) }' exch4 > dip_bnd
        dip_b="$(cat dip_bnd)"



        
                            scp -r $PATH2/SAVE ./.
                            scp -r $PATH2/r_setup ./.

                                                        $YAMBO_NL -r -e -v h+sex -V par -F $prefix.hxc.in 2> out

    #-----------------------------------------------COULOMB cut-off & RIM-------------------------------------------------#


                                        # Extract the value from r_setup using awk and calculate (value - 1)
                                            awk '/Alat factors/ {printf "%.6f|\n", $6 - 1}' r_setup > tmp1

                                        # Add leading zeros to match "0.000000|0.000000|value|"
                                            awk '{print "0.000000|0.000000|" $0}' tmp1 > tmp2

                             		# Replace line 34 of $prefix.hxc.in with contents of tmp2
					    awk 'NR==34 { while ((getline line < "tmp2") > 0) print line; next } { print }' "$prefix.hxc.in" > tmp4

					# Overwrite the original file
					    mv tmp4 "$prefix.hxc.in"

                                        # Clean up temporary files
                                            rm -f tmp1 tmp2 tmp3

                                            sed -i.bak 's/RandQpts=0/RandQpts='$randqpts'/g'                                            "$prefix.hxc.in"
                                            sed -i.bak 's/RandGvec= 1                RL/RandGvec= '"$randgvec $rl"'/g'                      "$prefix.hxc.in"
                                            sed -i.bak 's/CUTGeo= "none"/CUTGeo= "'"$slab_nl"'"/g'                                          "$prefix.hxc.in"
                                            sed -i.bak 's/DBsIOoff= "none"/DBsIOoff= "BS"/g'                                                "$prefix.hxc.in"
                                            sed -i.bak 's/DBsFRAGpm= "none"/DBsFRAGpm= "+BS"/g'                                             "$prefix.hxc.in"
    

         #--------------------------------------------CPU parallel structure----------------------------------------------#

                                            sed -i.bak "s/RT_CPU= \"\"/RT_CPU=\"$rt_cpu\"/g"                                                "$prefix.hxc.in"
                                            sed -i.bak "s/RT_ROLEs= \"\"/RT_ROLEs= \"k b q qp\"/g"                                          "$prefix.hxc.in"

                                            sed -i.bak "s/X_and_IO_CPU= \"\"/X_and_IO_CPU=\"$x_and_io\"/g"                                  "$prefix.hxc.in"
                                            sed -i.bak "s/X_and_IO_ROLEs= \"\"/X_and_IO_ROLEs= \"q g k c v\"/g"                             "$prefix.hxc.in"

                                            sed -i.bak "s/DIP_CPU= \"\"/DIP_CPU=\"$dip_cpu\"/g"                                             "$prefix.hxc.in"
                                            sed -i.bak "s/DIP_ROLEs= \"\"/DIP_ROLEs= \"k c v\"/g"                                           "$prefix.hxc.in"



#---------------------------------------------Nonlinear simulation parameters-------------------------------------------#

                                            # Delete line 41 (works on both systems): Polarization bands
                                            sed '41d' "$prefix.hxc.in" > tmp1 && mv tmp1                                                    "$prefix.hxc.in"
                                            awk "NR==41{print \"1| $dip_b|\"}1" "$prefix.hxc.in" > tmp && mv tmp                            "$prefix.hxc.in"
 
                                              # Delete line 48 (works on both systems): Collision bands
                                                           awk -v n=48 -v f="insert_block" 'NR==n {
                                                                while ((getline x < f) > 0) print x
                                                                    next
                                                                } {print}' "$prefix.hxc.in" > temp_file && mv temp_file "$prefix.hxc.in"
   
 

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

                                                        sed -i.bak '45s/1.000000 | 0.000000 | 0.000000 |/'"$vec"'/g'                     	"$prefix.hxc.in"
                                                   
#---------------------------------------------Exchange and correlation cut-offs----------------------------------------#

                                            sed '43d' "$prefix.hxc.in" > tmp1 && mv tmp1 "$prefix.hxc.in"
                                            awk "NR==43{print \"NGsBlkXs= $con    $ryd\"}1" "$prefix.hxc.in" > tmp && mv tmp        		"$prefix.hxc.in"	
 
                                            sed '51d' "$prefix.hxc.in" > tmp1 && mv tmp1 "$prefix.hxc.in"
                                            awk "NR==51{print \"HARRLvcs= $bseexx $ryd\"}1" "$prefix.hxc.in" > tmp && mv tmp            	"$prefix.hxc.in"
                                                
                                            sed '52d' "$prefix.hxc.in" > tmp1 && mv tmp1 "$prefix.hxc.in"
                                            awk "NR==52{print \"EXXRLvcs= $bseexx $ryd\"}1" "$prefix.hxc.in" > tmp && mv tmp           		"$prefix.hxc.in"
 

                                   	    awk -v val1="$con" -v val2="$ryd" 'NR==53 { print "CORRLvcs= " val1, val2; next } { print }' "$prefix.hxc.in" > tmp && mv tmp "$prefix.hxc.in"

                                                   
                                        $MPIRUN_YAMBO $YAMBO_NL -F $prefix.hxc.in 2> out



    
                                            {
                                                        echo "End of nonlinear/HX collision for Real-time BSE calculations."
                                                        echo ""
                                            } >> "$DIR/work_report.txt"


exit;
