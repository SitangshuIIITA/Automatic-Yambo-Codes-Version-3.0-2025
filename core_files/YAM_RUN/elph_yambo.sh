#!/bin/bash



cd $DIR


                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] Computation of initial and scattered electronic states |nk> and |n'k+q> report:"
                                } >> "$DIR/work_report.txt"



    PATH2=$DIR/scf_nonsoc
    PATH3=$DIR/nscf
    PATH5=$DIR/band_structure
    PATH12=$DIR/phonon
    PATH13=$DIR/dvscf
	
		
                          # The self-consistent file and this file defines the number of k-points for the density and the plane-wave cutoff
                          
                               SCF_FILE="$prefix.scf.in"
                               PREFIX="$prefix"
                               DYN_FILE="$prefix.dyn"

                               # The files for nscf, phonons, and dvscf without the q-points
                               # in the nscf file we define the number of bands
                               NSCF_FILE="$prefix.nscf.in"
                               PHONON_FILE="$prefix.phonon.in"
                               DVSCF_FILE="$prefix.dvscf.in"
                               DYN0_FILE="$prefix.dyn0"

        STARTPWD=$DIR

            Q_GRID_NAME="ELPH" # Name of the folder where calculations are performed

			cd $Q_GRID_NAME

#			# Run the dvscf
			mkdir dvscf
			cd dvscf
				cat ../../dvscf/$DVSCF_FILE > $DVSCF_FILE
				cat $STARTPWD/qpoints  >> $DVSCF_FILE
				scp -r ../nscf/$PREFIX.save ./
				scp -r ../phonon/_ph0 ./
				scp ../phonon/$DYN_FILE* ./
			echo "Running DVSCF ... "
			$PH_RUN_dvscf -inp $DVSCF_FILE > output_dvscf 2> out




#			#Read gkkp matrix elements

			cd $PREFIX.save
				echo "gkkp"> ypp.in_dbs
				echo "gkkp_db" >> ypp.in_dbs
				echo 'DBsPATH= "../elph_dir/"' >> ypp.in_dbs
				echo "Running YPP_PH gkkp ... "
				$YPP_PH -nompi -F ypp.in_dbs 2> output_dbs 2> out


            
                            {
                            echo "    → Scattering states written."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of scattering states report."
                            echo ""
                           } >> "$DIR/work_report.txt"


			cd $STARTPWD
 
exit;
