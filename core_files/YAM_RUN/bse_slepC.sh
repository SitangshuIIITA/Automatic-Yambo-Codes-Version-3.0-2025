#!/bin/bash



cd $DIR

        PATH1=$DIR/GW_folder
        PATH2=$DIR/BSE_folder
        PATH3=$DIR/BSE_folder/BSE_slepC
        PATH4=$DIR/GW_folder/GW_selfC
        
        mkdir -p $PATH2/BSE_slepC
        mkdir -p $DIR/graph_data/BSE_slepC_data

        cd $PATH3
        
                
                          {
                            echo ""
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] BSE-SlepC Q!=0 calculation report:"
                           } >> "$DIR/work_report.txt"

        
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


#-----------------------------------------------Optical Gap is fixed from converged BSE--------------------------------------------#

                                    if [ "$BSE_gap" = "g0w0" ]; then
                                    
                                        {
                                              echo "    → This is BSE-SlepC. QP gap is fixed from              : converged G0W0 result."
                                        } >> "$DIR/work_report.txt"
                                                                    
                                        
                                    elif [ "$BSE_gap" = "gw_selfC" ]; then
                                        {
                                              echo "    → This is BSE-SlepC. QP gap is fixed from              : converged self-consistent G"$gw_self_inter"W"$gw_self_inter"."
                                            } >> "$DIR/work_report.txt"
                                       
                                    fi

        
        scp -r $PATH2/"BSE-$con$ryd"/SAVE ./.
                
        scp -r $PATH2/"BSE-$con$ryd"/output ./.
        
                cd output
                    rm -rf ndb.BS*
                cd ..
                
        scp -r $PATH2/"BSE-$con$ryd"/Report/r-output_rim_cut_rim_w_optics_dipoles_bss_bse_em1s ./.
                       
                    
        scp -r $PATH2/"BSE-$con$ryd"/r_setup ./.
        
        
#-----------------------------------------------BSE + SlepC-Diagonalization solver--------------------------------------------#
        
            rm -rf "$bse_slepC_filename"
            $YAMBO -r -rw -X s -optics b -kernel sex -Ksolver s -V all -F $bse_slepC_filename 2> out


#-----------------------------------------------COULOMB cut-off & RIM-------------------------------------------------#
                        
                                            sed -i.bak 's/rimw_type= "default"/rimw_type= "'"$rimw_type"'"/g'                               "$bse_slepC_filename"
                                            sed -i.bak 's/RandGvecW= 1               RL/RandGvecW= '"$randgvecw $rl"'/g'                    "$bse_slepC_filename"
                                            sed -i.bak 's/RandQpts=0/RandQpts='$randqpts'/g'                                            "$bse_slepC_filename"
                                            sed -i.bak 's/RandGvec= 1                RL/RandGvec= '"$randgvec $rl"'/g'                      "$bse_slepC_filename"
                                            sed -i.bak 's/CUTGeo= "none"/CUTGeo= "'"$slab"'"/g'                                             "$bse_slepC_filename"
                                            sed -i.bak 's/DBsIOoff= "none"/DBsIOoff= "BS"/g'                                                "$bse_slepC_filename"
                                            sed -i.bak 's/DBsFRAGpm= "none"/DBsFRAGpm= "+BS"/g'                                             "$bse_slepC_filename"
    

#--------------------------------------------CPU parallel structure---------------------------------------------------#

                                            sed -i.bak "s/BS_CPU= \"\"/BS_CPU=\"$bs_cpu\"/g"                                                "$bse_slepC_filename"
                                            sed -i.bak "s/BS_ROLEs= \"\"/BS_ROLEs= \"k eh t\"/g"                                            "$bse_slepC_filename"

                                            sed -i.bak "s/X_and_IO_CPU= \"\"/X_and_IO_CPU=\"$x_and_io\"/g"                                  "$bse_slepC_filename"
                                            sed -i.bak "s/X_and_IO_ROLEs= \"\"/X_and_IO_ROLEs= \"q g k c v\"/g"                             "$bse_slepC_filename"

                                            sed -i.bak "s/DIP_CPU= \"\"/DIP_CPU=\"$dip_cpu\"/g"                                             "$bse_slepC_filename"
                                            sed -i.bak "s/DIP_ROLEs= \"\"/DIP_ROLEs= \"k c v\"/g"                                           "$bse_slepC_filename"

#---------------------------------------------FFT-GVecs and other Response size block------------------------------------#

                                            sed '184d' "$bse_slepC_filename" > tmp1 && mv tmp1 "$bse_slepC_filename"
                                            awk "NR==184{print \"NGsBlkXs= $con $ryd\"}1" "$bse_slepC_filename" > tmp && mv tmp             "$bse_slepC_filename"
                                                
                                            sed '88d' "$bse_slepC_filename" > tmp1 && mv tmp1 "$bse_slepC_filename"
                                            awk "NR==88{print \"BSENGBlk= $con $ryd\"}1" "$bse_slepC_filename" > tmp && mv tmp              "$bse_slepC_filename"
 

                                            sed '82d' "$bse_slepC_filename" > tmp1 && mv tmp1 "$bse_slepC_filename"
                                            awk "NR==82{print \"BSENGexx= $bseexx $ryd\"}1" "$bse_slepC_filename" > tmp && mv tmp           "$bse_slepC_filename"
 

#---------------------------------------------Includes el-hole coupling--------------------------------------------------#

                                            sed -i.bak 's/BSEmod= "resonant"/BSEmod= "coupling"/g'                                          "$bse_slepC_filename"
                                            sed -i.bak 's/Lkind= "BAR"/Lkind= "'"$lkind"'"/g'                                               "$bse_slepC_filename"
                                            sed -i.bak 's/#WehCpl/WehCpl/g'                                                                 "$bse_slepC_filename"


#---------------------------------------------Includes location of QP gap correction----------------------------------------#

                                            sed -i.bak 's/KfnQPdb= "none"/KfnQPdb= "E < .\/output\/ndb.QP"/g'                               "$bse_slepC_filename"


#---------------------------------------------Includes optical transition bands ----------------------------------------#
                                                                                                                                    
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
                                                                } {print}' "$bse_slepC_filename" > temp_file && mv temp_file "$bse_slepC_filename"
                                                                
#----------------------------Includes photon energy range, W screening bands, writes exciton wavefunctions---------------#


                                            sed -i.bak '136s/0.100000 | 0.100000 |/'"$Bdmrange"' | '"$Bdmrange"' |/g'                       "$bse_slepC_filename"
                                            sed -i.bak '133s/  0.00000 | 10.00000 |/  0.00000 | '"$phot_max_en_lim"' |/g'                   "$bse_slepC_filename"

                                            sed -i.bak 's/BEnSteps= 100/BEnSteps= 1000/g'                                                   "$bse_slepC_filename"
                                            sed -i.bak 's/#WRbsWF/WRbsWF/g'                                                                 "$bse_slepC_filename"
 
                                             # Delete line 175 (works on both systems): Polarization bands
                                            sed '182d' "$bse_slepC_filename"> tmp1 && mv tmp1                                               "$bse_slepC_filename"
                                            awk "NR==182{print \"1| $dip_b|\"}1" "$bse_slepC_filename" > tmp && mv tmp                      "$bse_slepC_filename"
 
                                              # Delete line 72 (works on both systems): Dip bands
                                            sed '72d' "$bse_slepC_filename" > tmp1 && mv tmp1                                               "$bse_slepC_filename"
                                            awk "NR==72{print \"1| $dip_b|\"}1" "$bse_slepC_filename" > tmp && mv tmp                       "$bse_slepC_filename"
 
 
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

                                                        sed -i.bak '141s/1.000000 | 0.000000 | 0.000000 |/'"$vec"'/g'                        "$bse_slepC_filename"
                                                        sed -i.bak '194s/1.000000 | 0.000000 | 0.000000 |/'"$vec"'/g'                        "$bse_slepC_filename"

							sed -i.bak 's/BSSSlepcApproach= "none"/BSSSlepcApproach= "Krylov-Schur"/g'           "$bse_slepC_filename"
							



#----------------------------BSE Finite Q points--------------------------------------------------------------------------#

                                                    # Normal BSE band selection
                                                            Last_Q=$(awk '/IBZ Q-points/ {print $4}' r_setup)
                                                            echo "1 |${Last_Q} |" > insert_Q
                                                            grep -v '^$' insert_Q> insert_Q.clean && mv insert_Q.clean insert_Q

                                                            awk -v n=118 -v f="insert_Q" 'NR==n {
                                                                while ((getline x < f) > 0) print x
                                                                    next
                                                                } {print}' "$bse_slepC_filename" > temp_file && mv temp_file "$bse_slepC_filename"



#-------------------------------SLEPC solver: first n excitons------------------------------------------------------------#

                                                        sed -i.bak 's/BSSNEig=0/BSSNEig='$bssneig'/g'                                       "$bse_slepC_filename"
                                                        sed -i.bak 's/BSSEnTarget= 0.000000/BSSEnTarget= '$target'/g'                       "$bse_slepC_filename"



#-----------------------------------------------------Running Yambo:------------------------------------------------------#

                                                        $MPIRUN_YAMBO $YAMBO -F "$bse_slepC_filename" -J output -C Report 2> out
                                                    
                                                    	cd Report
                                                    	
                                                    	scp -r o-output.alpha_q* $DIR/graph_data/BSE_slepC_data/.
                                                    	
                                                    	
                                                    	cd ..
                                                    
                                            {
                                                        echo "    → Response block size                                  : NGsBlkXs = $con Ry"
                                                        echo "    → Response block size                                  : BSENGexx = $bseexx Ry"
                                                        echo "    → Response block size                                  : BSENGBlk = $con Ry"
                                                        echo "    → Polarization/Screening Bands                         : Scheme   = 1 | $dip_b"
                                                        echo "    → Exciton Q points                                     : Scheme   = 1 | $insert_Q"
                                                        echo "    → BSE Transition Bands                                 : Scheme = $(cat insert_block)"
                                                        echo "    → See detailed report in $DIR/graph_data/BSE_slepC_data."
                                                        echo ""
                                            } >> "$DIR/work_report.txt"
   
                                                 
                                                    
                            cd ..
        
exit;
