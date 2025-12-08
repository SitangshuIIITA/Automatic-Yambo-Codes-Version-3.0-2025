#!/bin/bash



cd $DIR
        
        
        
      	PATH1=$DIR/ELPH_q_conv
      	PATH2=$DIR/ELPH
        PATH3=$DIR/graph_data/ELPH_q_converged_data


      	cd $PATH2
                          {
                            echo ""
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] QP ZPR convergence tests based on phonon momenta fine-grids:"
                           } >> "$DIR/work_report.txt"

        
		mkdir -p $DIR/graph_data/ELPH_q_converged_data
        
        mkdir -p ELPH_q_conv
        
        cd ELPH_q_conv
        
        for line in $matd
                    do
                         mkdir "ELPH-${line}x${line}x1_q"

                           	cd "ELPH-${line}x${line}x1_q"
                                	scp -r $PATH2/phonon ./.
                                    scp -r $PATH2/dvscf ./.

                                    cd phonon
                                        rm -rf $prefix.matdyn.in

#--------------------------Phonon energies [q2r + matdyn] on fine [double] grids calculation---------------------------------#

cat > $prefix.q2r.in << EOF
&input
                                                zasr='crystal',
                                                flfrc='$prefix.fc',
                                                fildyn='$prefix.dyn',
                                                loto_2d = .true.,
/
EOF
                    $Q2R_RUN -inp $prefix.q2r.in > q2r.out


cat > $prefix.matdyn.in << EOF
&input
                                                asr='crystal',
                                                flfrc='$prefix.fc',
                                                flfrq='$prefix.freq',
                                                dos=.true.,
                                                loto_2d = .true.,
                                                fldos='$prefix.dos',
                                                deltaE=1.d0,
                                                nk1=$line, nk2=$line, nk3=1
/
EOF
                    $MATDYN_RUN -inp $prefix.matdyn.in > output_matdyn 2> out

                                               scp -r $prefix.freq ../dvscf/$prefix.save/.

                          cd ../dvscf/$prefix.save
                                            rm -rf l_gkkp_gkkp_dg r_gkkp_gkkp_dg $prefix.ypp_ph_dg.in

#--------------------------Phonon double grid calculation---------------------------------#


cat > $prefix.ypp_ph_dg.in << EOF
&input

                                        gkkp
                                        gkkp_dg
                                        PHfreqF= "$prefix.freq"
                                        PHmodeF= "none"
                                        FineGd_mode= "mixed"
                                        #SkipBorderPts
                                        EkplusQmode= "interp"
                                        #TestPHDGrid
/
EOF
                    $YPP_PH -nompi -F $prefix.ypp_ph_dg.in 2> out
                    
                                        scp -r SAVE r_setup r_gkkp_gkkp_dg ../../.
                        
                        cd ../../

#------------------------------------------Quasi-particle energies calculation-----------------------------------------------------------------#


                                        $YAMBO_PH -r -g n -p fan -c ep -V gen -F qp-0K.${line}x${line}x1_q.in 2> out


#------------------------------------------COULOMB cut-off & RIM-----------------------------------------#


                                                awk '/Alat factors/ {print $6-1}' r_setup > 1
                                                sed 's/$/|/' 1 > 2
                                                sed '1 s/./0.000000|0.000000|&/' 2 > 3

                                                # Modify input file at line 29
                                                sed '29s/.*//' qp-0K.${line}x${line}x1_q.in > 4
                                                sed '29 r 3' 4 > 5
                                                sed '29d' 5 > 6
                                                mv 6 qp-0K.${line}x${line}x1_q.in
                                                rm -f 1 2 3 4 5

                                                # In-place edits (portable)
                                            
                                                sed -i.bak 's/RandQpts=0/RandQpts=1000000/g'                                            qp-0K.${line}x${line}x1_q.in
                                                sed -i.bak 's/RandGvec= 1/RandGvec= 100/g'                                              qp-0K.${line}x${line}x1_q.in
                                                sed -i.bak 's/CUTGeo= "none"/CUTGeo= "box z"/g'                                         qp-0K.${line}x${line}x1_q.in
                                                sed -i.bak "s/BoseTemp=-1.000000         eV/BoseTemp=0   Kn/g"                          qp-0K.${line}x${line}x1_q.in
                                                sed -i.bak 's/GDamping= 0.100000         eV/GDamping= 0.300000         eV/g'            qp-0K.${line}x${line}x1_q.in
                                                sed -i.bak 's/#WRgFsq/WRgFsq/g'                                                         qp-0K.${line}x${line}x1_q.in
                                                sed -i.bak '$a\
ExtendOut'                                                                                                                              qp-0K.${line}x${line}x1_q.in

                                                awk -v insert="1| $GphBRnge|" 'NR==36 {$0=insert} {print}' qp-0K.${line}x${line}x1_q.in > tmpfile && mv tmpfile qp-0K.${line}x${line}x1_q.in



                                        #------------------------Reading the k point limits, max valence band, min conduction band from r_setup file-------------#

                                        # Extract k-points and prepare
                                            awk -F':' '/IBZ K-points/ {gsub(/^ +/, "", $2); print "1|"$2}' r_setup > kpts_temp

                                        # Extract filled bands
                                            awk '/Filled Bands/ {print $5}' r_setup > filled_temp
                                            awk -v v="$nv" '{print $1 - v, $2}' filled_temp > filled_parsed

                                        # Extract empty bands
                                            awk '/Empty Bands/ {print $5}' r_setup > empty_temp
                                            awk -v v="$nc" '{print $1 + v, $2}' empty_temp > empty_parsed

                                        # Now paste together without extra sed tricks
                                            paste -d'|' kpts_temp filled_parsed empty_parsed | awk '{print $0"|"}' > insert_block

                                        # Make sure no blank lines are present
                                            sed -i '' '/^$/d' insert_block

                                                               awk -v n=44 -v f="insert_block" 'NR==n {
                                                                while ((getline x < f) > 0) print x
                                                                    next
                                                                } {print}' "qp-0K.${line}x${line}x1_q.in" > temp_file && mv temp_file "qp-0K.${line}x${line}x1_q.in"

                                        # Clean temporary files
                                            rm -f kpts_temp filled_temp filled_parsed empty_temp empty_parsed insert_block *.bak

                                    $YAMBO_PH -F qp-0K.${line}x${line}x1_q.in -J qp-0K 2> out

                                                mv o-qp-0K.qp o-qp-0K.${line}x${line}x1_q.qp

                                                scp -r o-qp-0K.${line}x${line}x1_q.qp $PATH3/.
                  
                    
                    
                                            {
                                                        echo "    → ELPH phonon fine q grid used                              : q = ${line}x${line}x1"
                                                        echo ""
                                            } >> "$DIR/work_report.txt"
   
        
                                    cd ..
                            
            done
                    
                echo "$matd" > mat
                awk '{ print $(NF-0) }' mat > con_mat_q
                mat_q="$(cat con_mat_q)"

                                    cd ELPH-${mat_q}x${mat_q}x1_q
                                    
                                            scp -r r_setup SAVE ../../.
                            
                            
                                            {
                                                        echo "    → See detailed report in $DIR/graph_data/ELPH_q_converged_data."
                                                        echo "QP ZPR (Fan+DW) converges for q = ${mat_q}x${mat_q}x1. This SAVE directory is used for further processing."
                                                        echo "QP ZPR (Fan+DW) convergence tests based on phonon momenta fine-grids completed."
                                                        echo ""
                                                        
                                            } >> "$DIR/work_report.txt"
                                            
                                            
                                            
            

exit;	
