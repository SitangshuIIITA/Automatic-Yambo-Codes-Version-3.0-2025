#!/bin/bash



cd $DIR

    PATH1=$DIR/ELPH
		
	
	cd $PATH1

                          {
                            echo ""
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Quasi-particle Electronic Eliashberg Function Calculations:"
                           } >> "$DIR/work_report.txt"


echo "$temp_value" > Tempera
awk '{ print $(NF-0) }' Tempera > Temp_fixed
TEMP_used="$(cat Temp_fixed)"


		for line in $TEMP_used
			do
   				cd QP-${line}$Kelvin
				
#------------------------------------------ELectronic Eliashberg Function for Indirect and Direct Band-Gap-----------------------------------------------------------------#


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
  
#--------------------------------------------------Eliashberg Function: Indirect and direct gap----------------------------------------------#


                                                $YPP_PH -s e -F eliash.${line}$Kelvin.ind.in 2> out

                                                # Replace patterns in input file
                                                sed 's/PhBroad= 0.010000/PhBroad= 0.0010000/g; s/PhStps= 200/PhStps= 1000/g' eliash.${line}$Kelvin.ind.in > tmp.in && mv tmp.in eliash.${line}$Kelvin.ind.in

                                              
                                                # Replace line 21 in input file : indirect line points
                                                awk -v newline1="$new_line1" 'NR==21 {$0=newline1} {print}' eliash.${line}$Kelvin.ind.in > new_eliash.in && mv new_eliash.in eliash.${line}$Kelvin.ind.in
                                                
                                                #Run YPP_PH
                                                $YPP_PH -F eliash.${line}$Kelvin.ind.in -J qp-${line}$Kelvin 2> out
                                                
                                                
                                                 # Replace line 21 in input file : indirect line points
                                                awk -v newline2="$new_line2" 'NR==21 {$0=newline2} {print}' eliash.${line}$Kelvin.ind.in > new_eliash.in && mv new_eliash.in eliash.${line}$Kelvin.ind.in

                                                #Run YPP_PH
                                                $YPP_PH -F eliash.${line}$Kelvin.ind.in -J qp-${line}$Kelvin 2> out
                                            
                                                
                                                # Compute difference and write to file
                                                file1="o-qp-${line}${Kelvin}.g_sq_F_b_${val4}_k_${val2}"        # indirect conduction.
                                                file2="o-qp-${line}${Kelvin}.g_sq_F_b_${val3}_k_${val1}"        # top valence.
                                                file3="o-qp-${line}${Kelvin}.g_sq_F_b_${val4}_k_${val1}"        # direct conduction.

                                                outfile1="${prefix}.${line}${Kelvin}.eliash_diff_direct.txt"
                                                outfile2="${prefix}.${line}${Kelvin}.eliash_diff_indirect.txt"

                                                paste "$file1" "$file2" | awk '{print $1, $2 - $10}' > "$outfile2"
                                                paste "$file3" "$file2" | awk '{print $1, $2 - $10}' > "$outfile1"

                                                # Remove empty and comment lines
                                                sed '/^\s*$/d; /^\s*#/d' $prefix.${line}$Kelvin.eliash_diff_indirect.txt > tmp && mv tmp $prefix.${line}$Kelvin.eliash_diff_indirect.txt

                                                # Cleanup
                                                rm -f f1 f2 f3 f4 insert_line

                                            

                  cd ..
            done

                                            {
                                                        echo "    → See output files in $DIR/graph_data/Spectral_function/."
                                                        echo "    → Note: Eliashberg difference file is created at both direct and indirect gaps. Difference taken is g^2F(conduction)-g^2F(valence)."
                                                        echo "End of Electronic Eliashberg function calculations."
                                            } >> "$DIR/work_report.txt"


exit;
