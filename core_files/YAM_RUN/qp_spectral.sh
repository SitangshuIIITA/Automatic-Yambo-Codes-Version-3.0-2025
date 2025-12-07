#!/bin/bash



cd $DIR


    PATH1=$DIR/ELPH
    mkdir $DIR/graph_data/Spectral_function

                          {
                            echo ""
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Quasi-particle Spectral Function Calculations:"
                           } >> "$DIR/work_report.txt"

  
	
	cd $PATH1


		for line in $temp_value
			do
   				cd QP-${line}$Kelvin
				
#------------------------------------------Spectral Function for Indirect and Direct Band-Gap-----------------------------------------------------------------#


                                                $YAMBO_PH -r -g g -p fan -c ep -V gen -F sf-${line}$Kelvin.ind_gap.in 2> out


#------------------------------------------COULOMB cut-off & RIM-------------------------------------------------------------------------------#

                                                awk '/Alat factors/ {print $6-1}' r_setup > 1
                                                sed 's/$/|/' 1 > 2
                                                sed '1 s/./0.000000|0.000000|&/' 2 > 3

                                                # Modify input file at line 29
                                                sed '29s/.*//' sf-${line}$Kelvin.ind_gap.in > 4
                                                sed '29 r 3' 4 > 5
                                                sed '29d' 5 > 6
                                                mv 6 sf-${line}$Kelvin.ind_gap.in
                                                rm -f 1 2 3 4 5

                                                # In-place edits (portable)
                                            
                                                sed -i.bak 's/RandQpts=0/RandQpts=1000000/g'                                            sf-${line}$Kelvin.ind_gap.in
                                                sed -i.bak 's/RandGvec= 1/RandGvec= 100/g'                                              sf-${line}$Kelvin.ind_gap.in
                                                sed -i.bak 's/CUTGeo= "none"/CUTGeo= "box z"/g'                                         sf-${line}$Kelvin.ind_gap.in
                                                sed -i.bak "s/BoseTemp=-1.000000         eV/BoseTemp=${line}   Kn/g"                    sf-${line}$Kelvin.ind_gap.in
                                                sed -i.bak 's/GEnSteps= 100/GEnSteps= 1000/g'                                           sf-${line}$Kelvin.ind_gap.in

                                                awk 'NR==44 { print "-2| 2|  eV"; next } { print }' sf-${line}$Kelvin.ind_gap.in > tmp && mv tmp sf-${line}$Kelvin.ind_gap.in
                                                awk "NR==47 { print \"$GDmRnge| $GDmRnge|  eV\"; next } { print }" sf-${line}$Kelvin.ind_gap.in > tmp && mv tmp sf-${line}$Kelvin.ind_gap.in

                                                sed -i.bak 's/#WRgFsq/WRgFsq/g'                                                         sf-${line}$Kelvin.ind_gap.in
                                                sed -i.bak '$a\
ExtendOut'                                                                                                                              sf-${line}$Kelvin.ind_gap.in

                                                awk -v insert="1| $GphBRnge|" 'NR==36 {$0=insert} {print}' sf-${line}$Kelvin.ind_gap.in > tmpfile && mv tmpfile sf-${line}$Kelvin.ind_gap.in


#----------------------------------Reading the k point limits, max valence band, min conduction band from r_setup file-------------------------#

                                                # Step 1: Extract values from `r_setup`
                                                val1=$(awk '/Indirect Gap between kpts/ {print $7}' r_setup)
                                                val2=$(awk '/Indirect Gap between kpts/ {print $8}' r_setup)
                                                val3=$(awk '/Filled Bands/ {print $5}' r_setup)
                                                val4=$(awk '/Empty Bands/ {print $5}' r_setup)

                                                # Step 2: Format the line for indirect k-points
                                                new_line1="${val2}|${val2}|${val4}|${val4}|"

                                                # Step 3: Replace line 50 in target file
                                                # - Use awk for line replacement (portable across platforms)
                                                awk -v newline1="$new_line1" 'NR==50 {$0=newline1} {print}' sf-${line}$Kelvin.ind_gap.in > tmpfile && mv tmpfile sf-${line}$Kelvin.ind_gap.in


                                            # Clean up
                                                rm -f tmp1 tmp2 tmp3 tmp4 tmp5 tmp_new tmp_new2 *.bak

                                                $YAMBO_PH -F sf-${line}$Kelvin.ind_gap.in -J sf-${line}$Kelvin 2> out
                    
                                                rm -rf sf-${line}$Kelvin
                                                
                                                new_line2="${val1}|${val1}|${val3}|${val4}|"
                                                
                                                # Step 3: Replace line 50 in target file
                                                # - Use awk for line replacement (portable across platforms)
                                                awk -v newline2="$new_line2" 'NR==50 {$0=newline2} {print}' sf-${line}$Kelvin.ind_gap.in > tmpfile && mv tmpfile sf-${line}$Kelvin.ind_gap.in


                                            # Clean up
                                                rm -f tmp1 tmp2 tmp3 tmp4 tmp5 tmp_new tmp_new2 *.bak

                                                $YAMBO_PH -F sf-${line}$Kelvin.ind_gap.in -J sf-${line}$Kelvin 2> out
                                                
                                                scp -r o-sf* $DIR/graph_data/Spectral_function/.
                    
                                            {
                                                        echo "    → Temperature used                              : T = ${line} K"
                                                        echo ""
                                            } >> "$DIR/work_report.txt"




  				cd ..
			done	 

                                            {
                                                        echo "    → See output files in $DIR/graph_data/Spectral_function/."
                                                        echo ""
                                                        echo "End of Spectral function calculations."
                                            } >> "$DIR/work_report.txt"


exit;
