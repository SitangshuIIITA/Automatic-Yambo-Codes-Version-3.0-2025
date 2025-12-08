#!/bin/bash


        cd $DIR
        
        PATH1=$DIR/nscf/$prefix.save
        
        
        mkdir -p Nonlinear
        mkdir -p $DIR/graph_data/Nonlinear_data

                          {
                            echo ""
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Nonlinear initialization report:"
                           } >> "$DIR/work_report.txt"

        
                            scp -r $PATH1/SAVE ./.
                            
                                                                            $YPP_NL -y 2> out
                                                                            
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

                                                        sed -i.bak '141s/1.000000 | 0.000000 | 0.000000 |/'"$vec"'/g'                                                           "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        sed -i.bak '184s/1.000000 | 0.000000 | 0.000000 |/'"$vec"'/g'                                                           "$prefix.BSE-${line}$Kelvin.diago.in"

        
        
        rm -rf BSE-*

        
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




            for line in $temp_value
                        do
                            mkdir "BSE-${line}$Kelvin"


                cd BSE-${line}$Kelvin

                    scp -r $PATH2/"BSE-$con$ryd"/SAVE ./.

                    scp -r $PATH2/"BSE-$con$ryd"/output ./.
                    
                    scp -r $PATH2/"BSE-$con$ryd"/r_setup ./.
        
                cd output
                    rm -rf ndb.BS* ndb.QP
                cd ..

                cd SAVE
                    rm -rf ndb.qp
                cd ..
                
                                           
#-----------------------------------------------Choise of QP gap: G0W0 or self-consistent GW and GW + temperature ELPH database merging ----------------------#

                                                $YPP -qpdb -m -F merge.in 2> out
                                                
                                                    # Delete line 19
                                                    sed '19d' merge.in > tmp && mv tmp merge.in


                                                        if [ "$BSE_gap" = "g0w0" ]; then
                                                                # Insert at line 19
                                                                    sed '19i\
"E" |"+" |"1" |"../../../GW_folder/GW_selfC/G0W0/output/ndb.QP" |
                                                                    ' merge.in > tmp && mv tmp merge.in
                                                               {
                                                                    echo "    → QP gap is fixed from              : converged G0W0."
                                                        
                                                                } >> "$DIR/work_report.txt"

                                        
                                                        elif [ "$BSE_gap" = "gw_selfC" ]; then
                                                               # Insert at line 19 (portable on macOS and Linux)
sed "19i\\
\"E\" |\"+\" |\"1\" |\"../../../GW_folder/GW_selfC/G${gw_self_inter}W${gw_self_inter}/G${gw_self_inter}W${gw_self_inter}/ndb.QP\" |
" merge.in > tmp && mv tmp merge.in

                                                                {
                                                                    echo "    → QP gap is fixed from              : converged self-consistent G"$gw_self_inter"W"$gw_self_inter"."
                                                        
                                                                } >> "$DIR/work_report.txt"

                                                        fi

                                                # Insert at line 20
                                                sed "20i\\
\"E W\" |\"+\" |\"1\" |\"../../../ELPH/QP-${line}${Kelvin}/qp-${line}${Kelvin}/ndb.QP\" |
                                                " merge.in > tmp && mv tmp merge.in
                                                
                                                $YPP -F merge.in 2> out


                                               $YAMBO_PH 2> out
                                               $YAMBO_PH  -r -rw -X s -optics b -kernel sex -Ksolver d -V all -F $prefix.BSE-${line}$Kelvin.diago.in 2> out

#-----------------------------------------------COULOMB cut-off & RIM-------------------------------------------------#
                        
                                                        sed -i.bak 's/rimw_type= "default"/rimw_type= "'"$rimw_type"'"/g'                                                       "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        sed -i.bak 's/RandGvecW= 1               RL/RandGvecW= '"$randgvecw $rl"'/g'                                            "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        sed -i.bak 's/RandQpts=0/RandQpts="'"$randqpts"'"/g'                                                                    "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        sed -i.bak 's/RandGvec= 1                RL/RandGvec= '"$randgvec $rl"'/g'                                              "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        sed -i.bak 's/CUTGeo= "none"/CUTGeo= "'"$slab"'"/g'                                                                     "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        sed -i.bak 's/DBsIOoff= "none"/DBsIOoff= "BS"/g'                                                                        "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        sed -i.bak 's/DBsFRAGpm= "none"/DBsFRAGpm= "+BS"/g'                                                                     "$prefix.BSE-${line}$Kelvin.diago.in"
    

#--------------------------------------------CPU parallel structure----------------------------------------------#

                                                        sed -i.bak "s/BS_CPU= \"\"/BS_CPU=\"$bs_cpu\"/g"                                                                        "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        sed -i.bak "s/BS_ROLEs= \"\"/BS_ROLEs= \"k eh t\"/g"                                                                    "$prefix.BSE-${line}$Kelvin.diago.in"

                                                        sed -i.bak "s/X_and_IO_CPU= \"\"/X_and_IO_CPU=\"$x_and_io\"/g"                                                          "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        sed -i.bak "s/X_and_IO_ROLEs= \"\"/X_and_IO_ROLEs= \"q g k c v\"/g"                                                     "$prefix.BSE-${line}$Kelvin.diago.in"

                                                        sed -i.bak "s/DIP_CPU= \"\"/DIP_CPU=\"$dip_cpu\"/g"                                                                     "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        sed -i.bak "s/DIP_ROLEs= \"\"/DIP_ROLEs= \"k c v\"/g"                                                                   "$prefix.BSE-${line}$Kelvin.diago.in"

#---------------------------------------------FFT-GVecs and other Response size block------------------------------------#

                                                        sed '174d' "$prefix.BSE-${line}$Kelvin.diago.in" > tmp1 && mv tmp1 "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        awk "NR==174{print \"NGsBlkXs= $con $ryd\"}1" "$prefix.BSE-${line}$Kelvin.diago.in" > tmp && mv tmp                     "$prefix.BSE-${line}$Kelvin.diago.in"
                                                
                                                        sed '88d' "$prefix.BSE-${line}$Kelvin.diago.in" > tmp1 && mv tmp1 "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        awk "NR==88{print \"BSENGBlk= $con $ryd\"}1" "$prefix.BSE-${line}$Kelvin.diago.in" > tmp && mv tmp                      "$prefix.BSE-${line}$Kelvin.diago.in"
 

                                                        sed '82d' "$prefix.BSE-${line}$Kelvin.diago.in" > tmp1 && mv tmp1 "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        awk "NR==82{print \"BSENGexx= $bseexx $ryd\"}1" "$prefix.BSE-${line}$Kelvin.diago.in" > tmp && mv tmp                   "$prefix.BSE-${line}$Kelvin.diago.in"
 
 
#---------------------------------------------Includes el-hole coupling--------------------------------------------------#

                                                        sed -i.bak 's/BSEmod= "resonant"/BSEmod= "coupling"/g'                                                                  "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        sed -i.bak 's/Lkind= "BAR"/Lkind= "'"$lkind"'"/g'                                                                       "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        sed -i.bak 's/#WehCpl/WehCpl/g'                                                                                         "$prefix.BSE-${line}$Kelvin.diago.in"

                                           
#---------------------------------------------Includes self_C_GW + gap correction----------------------------------------#

                                                        sed -i.bak 's/KfnQPdb= "none"/KfnQPdb= "E W < SAVE\/ndb.QP_merged_1_gw_ppa_el_ph"/g'                                        "$prefix.BSE-${line}$Kelvin.diago.in"

#----------------------------Includes photon energy range, W screening bands, writes exciton wavefunctions---------------#


                                                        sed -i.bak '136s/0.100000 | 0.100000 |/0.00100000 | 0.00100000 |/g'                                                     "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        sed -i.bak '133s/  0.00000 | 10.00000 |/  0.00000 | '"$phot_max_en_lim"' |/g'                                           "$prefix.BSE-${line}$Kelvin.diago.in"

                                                        sed -i.bak 's/BEnSteps= 100/BEnSteps= 1000/g'                                                                           "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        sed -i.bak 's/#WRbsWF/WRbsWF/g'                                                                                         "$prefix.BSE-${line}$Kelvin.diago.in"
 
                                             # Delete line 172 (works on both systems): Polarization bands
                                                        sed '172d' "$prefix.BSE-${line}$Kelvin.diago.in" > tmp1 && mv tmp1                                                      "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        awk "NR==172{print \"1| $dip_b|\"}1" "$prefix.BSE-${line}$Kelvin.diago.in" > tmp && mv tmp                               "$prefix.BSE-${line}$Kelvin.diago.in"
 
                                              # Delete line 72 (works on both systems): Dip bands
                                                        sed '72d' "$prefix.BSE-${line}$Kelvin.diago.in" > tmp1 && mv tmp1                                                       "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        awk "NR==72{print \"1| $dip_b|\"}1" "$prefix.BSE-${line}$Kelvin.diago.in" > tmp && mv tmp                                "$prefix.BSE-${line}$Kelvin.diago.in"
 
 
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

                                                        sed -i.bak '141s/1.000000 | 0.000000 | 0.000000 |/'"$vec"'/g'                                                           "$prefix.BSE-${line}$Kelvin.diago.in"
                                                        sed -i.bak '184s/1.000000 | 0.000000 | 0.000000 |/'"$vec"'/g'                                                           "$prefix.BSE-${line}$Kelvin.diago.in"

 #------------------Reading the k point limits, max valence band, min conduction band from r_setup file-------------#


                                                    # Normal BSE band selection
                                                            filled_val=$(awk '/Filled Bands/ {print $5}' r_setup)
                                                            empty_val=$(awk '/Empty Bands/ {print $5}' r_setup)

                                                            vb=$((filled_val - $nv))
                                                            cb=$((empty_val + $nc))

                                                            echo "${vb} |${cb} |" > insert_block
                                                            grep -v '^$' insert_block > insert_block.clean && mv insert_block.clean insert_block
                                                            scp -r insert_block ../.

                                                            awk -v n=121 -v f="insert_block" 'NR==n {
                                                                while ((getline x < f) > 0) print x
                                                                    next
                                                                } {print}' "$prefix.BSE-${line}$Kelvin.diago.in" > temp_file && mv temp_file "$prefix.BSE-${line}$Kelvin.diago.in"

#-----------------------------------------------------Running Yambo:----------------------------------------------------#

                                                    $MPIRUN_YAMBO $YAMBO_PH -F "$prefix.BSE-${line}$Kelvin.diago.in" -J output -C Report 2> out
   
                                                        cd Report
                                                    
                                                                mv o-output.alpha_q1_diago_bse $prefix.${line}$Kelvin.BSE_GW.o-output.alpha_q1_diago_bse
                                                    
                                                                scp -r $prefix.${line}$Kelvin.BSE_GW.o-output.alpha_q1_diago_bse $DIR/graph_data/BSE_Temperature_Data/.
                                                    
                                                        cd ..
                                                    
                                                    $YPP -J output -e s 1 2> out
                                                    
                                                                mv o-output.exc_qpt1_E_sorted $prefix.${line}$Kelvin.o-output.exc_qpt1_E_sorted
                                                                scp -r $prefix.${line}$Kelvin.o-output.exc_qpt1_E_sorted $DIR/graph_data/BSE_Temperature_Data/.
                                                    

                                           
                                           
                                            echo "$PWD"
                                             
                                         {
                                                        echo "    → Temperature used                              : T = ${line} K"
                                                        echo "    → See detailed report in $DIR/graph_data/BSE_Temperature_Data."
                                                        echo ""
                                            } >> "$DIR/work_report.txt"
                                             
                                             
                            cd ..
                            
                     done
                                         {
                                                        echo "    → Response block size used                : NGsBlkXs = $con Ry"
                                                        echo "    → Response block size used                : BSENGexx = $bseexx Ry"
                                                        echo "    → Response block size used                : BSENGBlk = $con Ry"
                                                        echo "    → Polarization/Screening Bands used       : Scheme = 1 | $dip_b"
                                                        echo "    → BSE Transition Bands used               : Scheme = $(cat insert_block)"
                                                        echo ""
                                                        echo "End of BSE-Tempreature calculations."
                                                        echo ""
                                            } >> "$DIR/work_report.txt"


exit;
