#!/bin/bash



cd $DIR

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] QP [Fan + Debye Waller] ZPR energies at various temperature report:"
                                } >> "$DIR/work_report.txt"



    PATH1=$DIR/ELPH
		
	mkdir $DIR/graph_data/QP_Temperature_data

	cd $PATH1


		for line in $temp_value
			do
   				mkdir "QP-${line}$Kelvin"
				scp -r SAVE ./QP-${line}$Kelvin
				cp -r r_setup ./QP-${line}$Kelvin
   				cd QP-${line}$Kelvin
				
#------------------------------------------Quasi-particle energies calculation-----------------------------------------------------------------#


					$YAMBO_PH -r -g n -p fan -c ep -V gen -F qp-${line}$Kelvin.in 2> out


#------------------------------------------COULOMB cut-off & RIM-------------------------------------------------------------------------------#

                                                awk '/Alat factors/ {print $6-1}' r_setup > 1
                                                sed 's/$/|/' 1 > 2
                                                sed '1 s/./0.000000|0.000000|&/' 2 > 3

                                                # Modify input file at line 29
                                                sed '29s/.*//' qp-${line}$Kelvin.in > 4
                                                sed '29 r 3' 4 > 5
                                                sed '29d' 5 > 6
                                                mv 6 qp-${line}$Kelvin.in
                                                rm -f 1 2 3 4 5

                                                # In-place edits (portable)
                                            
                                                sed -i.bak 's/RandQpts=0/RandQpts=1000000/g'                                            qp-${line}$Kelvin.in
                                                sed -i.bak 's/RandGvec= 1/RandGvec= 100/g'                                              qp-${line}$Kelvin.in
                                                sed -i.bak 's/CUTGeo= "none"/CUTGeo= "box z"/g'                                         qp-${line}$Kelvin.in
                                                sed -i.bak "s/BoseTemp=-1.000000         eV/BoseTemp=${line}   Kn/g"                    qp-${line}$Kelvin.in
                                                sed -i.bak 's/GDamping= 0.100000         eV/GDamping= 0.300000         eV/g'            qp-${line}$Kelvin.in
                                                sed -i.bak 's/#WRgFsq/WRgFsq/g'                                                         qp-${line}$Kelvin.in
                                                sed -i.bak '$a\
ExtendOut'                                                                                                                              qp-${line}$Kelvin.in

                                                awk -v insert="1| $GphBRnge|" 'NR==36 {$0=insert} {print}' qp-${line}$Kelvin.in > tmpfile && mv tmpfile qp-${line}$Kelvin.in


#----------------------------------Reading the k point limits, max valence band, min conduction band from r_setup file-------------------------#

 
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
                                                                } {print}' "qp-${line}$Kelvin.in" > temp_file && mv temp_file "qp-${line}$Kelvin.in"

                                        # Clean temporary files
                                            rm -f kpts_temp filled_temp filled_parsed empty_temp empty_parsed insert_block *.bak

					$YAMBO_PH -F qp-${line}$Kelvin.in -J qp-${line}$Kelvin 2> out

                                                mv o-qp-${line}$Kelvin.qp $prefix.o-qp-${line}$Kelvin.qp

                                                scp -r $prefix.o-qp-${line}$Kelvin.qp $DIR/graph_data/QP_Temperature_data/.
		
						
  				cd ..
			done	 

            
                            {
                            echo "    → ZPR energies written."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of ZPR report."
                            echo "    → See detailed gap report in $DIR/graph_data/QP_Temperature_data for details."
                            echo ""
                           } >> "$DIR/work_report.txt"


exit;
