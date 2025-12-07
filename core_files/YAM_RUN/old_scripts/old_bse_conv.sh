#!/bin/bash

    cd $DIR
    
        
	PATH6=$DIR/GW_folder/GW_selfC
	mkdir BSE_folder	


	cd BSE_folder

				# Extract and process variables
				bseexx=($bsengexx)
				ngsblk=($ngsblks)
				bseng=($bsengblk)
				NC="$nc" # Number of lowest conduction bands

				# Compute cbn and vbn
				cbn=($(seq 1 $NC))
				vbn=($(seq 1 $NC))

					# Write cbn to file
					echo "${cbn[*]}" > cbn

					# Clean up any stray "cbn" files
					rm -f cbn

					# Loop through each element of the arrays
					for i in "${!ngsblk[@]}"; do
					    char_ngsblk=${ngsblk[$i]}
					    char_bseexx=${bseexx[$i]}
					    char_bseng=${bseng[$i]}
					    dir_name="BSE-${char_ngsblk}$ryd"

				 	# Use "16 1 1" CPUs for the 1 1 transition bandss, otherwise use bs_cpu values
                                                if [ $i -lt 1 ]; then
                                                        char_bs_cpu="4 1 1"
                                                else
                                                        char_bs_cpu="$bs_cpu"
                                                fi	

    
				# Create directory if not exists
				mkdir -p "$dir_name"
    	
				cd $dir_name


	 #-----------------------------------------------BSE+GW : diago solver--------------------------------------------#

					# Selecting a converged GW folder
					
					echo "$bndrnge" > exch1
					awk '{ print $(NF-1) }' exch1 > exch2
					bndr="$(cat exch2)"
					set -- $bndrnge
					
                	
						scp -r $PATH6/G"$gw_self_inter"W"$gw_self_inter"/SAVE ./.
                        
                        scp -r $PATH6/G"$gw_self_inter"W"$gw_self_inter"/G"$gw_self_inter"W"$gw_self_inter" ./.	

						mv G"$gw_self_inter"W"$gw_self_inter" output



							$YAMBO 2> out
		                                $YAMBO -r -rw -X s -optics b -kernel sex -Ksolver d -V all -F $bse_gw_filename  2> out


	#-----------------------------------------------COULOMB cut-off & RIM-------------------------------------------------#
						
										   	sed -i.bak 's/rimw_type= "default"/rimw_type= "'"$rimw_type"'"/g'               				"$bse_gw_filename"
    										sed -i.bak 's/RandGvecW= 1               RL/RandGvecW= '"$randgvecw $rl"'/g'    				"$bse_gw_filename"
    										sed -i.bak 's/RandQpts=0/RandQpts="'"$randqpts"'"/g'                            				"$bse_gw_filename"
    										sed -i.bak 's/RandGvec= 1                RL/RandGvec= '"$randgvec $rl"'/g'      				"$bse_gw_filename"
    										sed -i.bak 's/CUTGeo= "none"/CUTGeo= "'"$slab"'"/g'                             				"$bse_gw_filename"
    
	                                       	sed -i.bak 's/DBsIOoff= "none"/DBsIOoff= "BS"/g'                           						"$bse_gw_filename"
                                        	sed -i.bak 's/DBsFRAGpm= "none"/DBsFRAGpm= "+BS"/g'                        						"$bse_gw_filename"
	

         #--------------------------------------------CPU parallel structure----------------------------------------------#

    										sed -i.bak "s/BS_CPU= \"\"/BS_CPU=\"$char_bs_cpu\"/g"                                			"$bse_gw_filename"
    										sed -i.bak "s/BS_ROLEs= \"\"/BS_ROLEs= \"k eh t\"/g"                            				"$bse_gw_filename"

   				 							sed -i.bak "s/X_and_IO_CPU= \"\"/X_and_IO_CPU=\"$x_and_io\"/g"                  				"$bse_gw_filename"
    										sed -i.bak "s/X_and_IO_ROLEs= \"\"/X_and_IO_ROLEs= \"q g k c v\"/g"             				"$bse_gw_filename"

    										sed -i.bak "s/DIP_CPU= \"\"/DIP_CPU=\"$dip_cpu\"/g"                             				"$bse_gw_filename"
    										sed -i.bak "s/DIP_ROLEs= \"\"/DIP_ROLEs= \"k c v\"/g"                           				"$bse_gw_filename"


#---------------------------------------------FFT-GVecs and other Response size block------------------------------------#

											sed '174d' "$bse_gw_filename" > tmp1 && mv tmp1 "$bse_gw_filename"
											awk "NR==174{print \"NGsBlkXs= $char_ngsblk $ryd\"}1" "$bse_gw_filename" > tmp && mv tmp		"$bse_gw_filename"	
                                                
											sed '88d' "$bse_gw_filename" > tmp1 && mv tmp1 "$bse_gw_filename"
											awk "NR==88{print \"BSENGBlk= $char_bseng $ryd\"}1" "$bse_gw_filename" > tmp && mv tmp			"$bse_gw_filename"	
 

											sed '82d' "$bse_gw_filename" > tmp1 && mv tmp1 "$bse_gw_filename"
											awk "NR==82{print \"BSENGexx= $char_bseexx $ryd\"}1" "$bse_gw_filename" > tmp && mv tmp			"$bse_gw_filename"	
 
 
#---------------------------------------------Includes el-hole coupling--------------------------------------------------#

                                            sed -i.bak 's/BSEmod= "resonant"/BSEmod= "coupling"/g'                      					"$bse_gw_filename"
                                            sed -i.bak 's/Lkind= "BAR"/Lkind= "'"$lkind"'"/g'				                      			"$bse_gw_filename"
                                            sed -i.bak 's/#WehCpl/WehCpl/g'                                             					"$bse_gw_filename"

#---------------------------------------------Includes self_C_GW + gap correction----------------------------------------#

                                            sed -i.bak 's/KfnQPdb= "none"/KfnQPdb= "E < .\/output\/ndb.QP"/g'       						"$bse_gw_filename"

       
#----------------------------Includes photon energy range, W screening bands, writes exciton wavefunctions---------------#


                                            sed -i.bak '136s/0.100000 | 0.100000 |/'"$Bdmrange"' | '"$Bdmrange"' |/g'             			"$bse_gw_filename"
											sed -i.bak '133s/  0.00000 | 10.00000 |/  0.00000 | '"$phot_max_en_lim"' |/g'	        		"$bse_gw_filename"

                                            sed -i.bak 's/BEnSteps= 100/BEnSteps= 1000/g'                                					"$bse_gw_filename"
                                            sed -i.bak 's/#WRbsWF/WRbsWF/g'                                           						"$bse_gw_filename"
 
     										# Delete line 172 (works on both systems): Polarization bands
											sed '172d' "$bse_gw_filename" > tmp1 && mv tmp1 												"$bse_gw_filename"
											awk "NR==172{print \"1| $bndr|\"}1" "$bse_gw_filename" > tmp && mv tmp 							"$bse_gw_filename"
 
									 	     # Delete line 72 (works on both systems): Dip bands
											sed '72d' "$bse_gw_filename" > tmp1 && mv tmp1 													"$bse_gw_filename"
											awk "NR==72{print \"1| $bndr|\"}1" "$bse_gw_filename" > tmp && mv tmp 							"$bse_gw_filename"
 
 
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

														sed -i.bak '141s/1.000000 | 0.000000 | 0.000000 |/'"$vec"'/g'             			"$bse_gw_filename"
														sed -i.bak '184s/1.000000 | 0.000000 | 0.000000 |/'"$vec"'/g'             			"$bse_gw_filename"


 #------------------Reading the k point limits, max valence band, min conduction band from r_setup file-------------#

                                                    # Extract original filled and empty band indices
                                                        filled_val=$(awk '/Filled Bands/ {print $5}' r_setup)
                                                        empty_val=$(awk '/Empty Bands/ {print $5}' r_setup)

                                                    # Calculate valence and conduction band indices for this iteration
                                                        vb=$((filled_val - i))
                                                        cb=$((empty_val + i))

                                                    # Write to insert_block
                                                        echo "${vb} |${cb} |" > insert_block

                                                        echo "BSE transition bands in $i convergence:"
                                                        cat insert_block
    
                                                    # Remove any blank lines (safe cross-platform way)
                                                        grep -v '^$' insert_block > insert_block.clean && mv insert_block.clean insert_block

                                                    # Insert into line 136 (replacing that line)
                                                            awk -v n=121 -v f="insert_block" 'NR==n {
                                                                    while ((getline x < f) > 0) print x
                                                                        next
                                                                        } {print}' "$bse_gw_filename" > temp_file && mv temp_file "$bse_gw_filename"


 #-----------------------------------------------------Running Yambo:----------------------------------------------------#

                                                    $MPIRUN_YAMBO $YAMBO -F "$bse_gw_filename" -J output -C Report 2> out


					cd ..

					done

exit 0

