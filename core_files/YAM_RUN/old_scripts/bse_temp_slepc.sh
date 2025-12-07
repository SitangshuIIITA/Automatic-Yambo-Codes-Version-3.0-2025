#!/bin/bash
#SBATCH -N 16
#SBATCH --ntasks-per-node=40
#SBATCH --job-name=sitangshu
#SBATCH --error=err.out
#SBATCH --output=out.out
#SBATCH --time=96:00:00
#SBATCH --partition=standard

        cd $SLURM_SUBMIT_DIR
        export I_MPI_FALLBACK=disable
        export I_MPI_FABRICS=shm:dapl
	module load gnu8/8.3.0 openmpi3/3.1.4

	mkdir BSE_Tempreature_slpec

	PATH14=$DIR/BSE_Tempreature_slpec

	cd $DIR/graph_data/

		mkdir BSE_Tempreature_data_slpec
	
	cd $PATH14

		for line in $temp_value
                        do
                                mkdir "BSE-${line}$Kelvin"

 
#       # 1. BSE + GW  spectra using slepC solver is computed at different temperatures. 

	#-----------------------------------------------BSE + GW + phonon Calculation: SlepC solver------------------------------#


						scp -r $DIR/database/SAVE ./BSE-${line}$Kelvin/.

						#scp -r $PATH7/Report/r-output_rim_cut_optics_dipoles_bss_bse_em1s ./BSE-${line}$Kelvin/.


                                       cd BSE-${line}$Kelvin
                                               cd SAVE
                                               rm -rf ndb.QP
                                               cd ..
		
		#----------------------------GW + temperature ELPH data merging-----------------------------------------------------#

                                                ypp -qpdb -m -F merge.in
                                                sed '19d' merge.in > 1
                                                mv 1 merge.in
                                                rm -rf 1
                                                sed -i '19i "E" |"+" |"1" |"../../GW_folder/output/ndb.QP" |' merge.in

                                                sed -i '20i "E W" |"+" |"1" |"'"../../ELPH/QP-${line}$Kelvin/qp-${line}$Kelvin/ndb.QP"'" |' merge.in

                                                ypp -F merge.in



                                               yambo_ph
                                               yambo_ph -r -rw -X s -optics b -kernel sex -Ksolver s -V all -F $prefix.BSE-${line}$Kelvin.slpc.in
		
				

  		#-----------------------------------------------COULOMB cut-off & RIM---------------------------------------------#

                                        sed -i 's/rimw_type= "default"/rimw_type= "'"$rimw_type"'"/g'                   $prefix.BSE-${line}$Kelvin.slpc.in
                                        sed -i 's/RandGvecW= 1               RL/RandGvecW= '"$randgvecw $rl"'/g'        $prefix.BSE-${line}$Kelvin.slpc.in
                                        sed -i 's/RandQpts=0/RandQpts='"$randqpts"'/g'                                  $prefix.BSE-${line}$Kelvin.slpc.in
                                        sed -i 's/RandGvec= 1                RL/RandGvec= '"$randgvec $rl"'/g'          $prefix.BSE-${line}$Kelvin.slpc.in
                                        sed -i 's/CUTGeo= "none"/CUTGeo= "'"$slab"'"/g'                                 $prefix.BSE-${line}$Kelvin.slpc.in

                                        sed -i 's/DBsIOoff= "none"/DBsIOoff= "BS"/g'                                    $prefix.BSE-${line}$Kelvin.slpc.in
                                        sed -i 's/DBsFRAGpm= "none"/DBsFRAGpm= "+BS"/g'                                 $prefix.BSE-${line}$Kelvin.slpc.in



                #----------------------------------------------CPU parallel structure---------------------------------------------#


                                                sed -i 's/X_and_IO_CPU= ""/X_and_IO_CPU="'"$x_and_io"'"/g'              $prefix.BSE-${line}$Kelvin.slpc.in
                                                sed -i 's/X_and_IO_ROLEs= ""/X_and_IO_ROLEs= "q g k c v"/g'             $prefix.BSE-${line}$Kelvin.slpc.in

                                                sed -i 's/DIP_CPU= ""/DIP_CPU="'"$dip_cpu"'"/g'                         $prefix.BSE-${line}$Kelvin.slpc.in
                                                sed -i 's/DIP_ROLEs= ""/DIP_ROLEs= "k c v"/g'                           $prefix.BSE-${line}$Kelvin.slpc.in

                                                sed -i 's/BS_CPU= ""/BS_CPU= "'"$bs_cpu"'"/g'                           $prefix.BSE-${line}$Kelvin.slpc.in
                                                sed -i 's/BS_ROLEs= ""/BS_ROLEs= "k eh t"/g'                            $prefix.BSE-${line}$Kelvin.slpc.in

               #----------------------------------------------FFT-GVecs and other Response size block-----------------------------#

                                                sed -i.bak -e '181d'                                                    $prefix.BSE-${line}$Kelvin.slpc.in
                                                sed -i.bak "181i NGsBlkXs= $ngsblks"                                    $prefix.BSE-${line}$Kelvin.slpc.in

                                                sed -i.bak -e '89d'                                                     $prefix.BSE-${line}$Kelvin.slpc.in
                                                sed -i.bak "89i BSENGBlk= $bsengblk"                                    $prefix.BSE-${line}$Kelvin.slpc.in

                                                sed -i.bak -e '83d'                                                     $prefix.BSE-${line}$Kelvin.slpc.in
                                                sed -i.bak "83i BSENGexx= $bsengexx"                                    $prefix.BSE-${line}$Kelvin.slpc.in



                #----------------------------In-plane electric field directions-------------------#             

                                                sed -i '142s/ 1.000000 | 0.000000 | 0.000000 |/ 1.000000 | 1.000000 | 0.000000 |/g'     $prefix.BSE-${line}$Kelvin.slpc.in
                                                sed -i '191s/ 1.000000 | 0.000000 | 0.000000 |/ 1.000000 | 1.000000 | 0.000000 |/g'     $prefix.BSE-${line}$Kelvin.slpc.in

                #----------------------------Includes el-hole coupling----------------------------#

                                                sed -i 's/BSEmod= "resonant"/BSEmod= "coupling"/g'                      $prefix.BSE-${line}$Kelvin.slpc.in
                                                sed -i 's/#WehCpl/WehCpl/g'                                             $prefix.BSE-${line}$Kelvin.slpc.in

		#----------------------------Include Bose Temperature-----------------------------#


						sed -i 's/BoseTemp=-1.000000         eV/BoseTemp='${line}'     K/g' $prefix.BSE-${line}$Kelvin.slpc.in							
						sed -i.bak -e '137d'                                                $prefix.BSE-${line}$Kelvin.slpc.in
                                                sed -i.bak "137i 0.0050000 | 0.0050000 |         eV"                $prefix.BSE-${line}$Kelvin.slpc.in

		#---------------------Includes GW + temperature dependent band gap energies--------#

                                                sed -i 's/KfnQPdb= "none"/KfnQPdb= "E W < SAVE\/ndb.QP_merged_1_gw_ppa_el_ph"/g'       $prefix.BSE-${line}$Kelvin.slpc.in

               	#------------------------Reading the k point limits, max valence band, min conduction band from r_setup file-------------#

                                                awk '/Filled Bands/ {print $5}' r_setup > 5
                                                awk -v v=$nv '{print $1 - v, $2}' 5 > 6


                                                awk '/Empty Bands/ {print $5}' r_setup > 7
                                                awk -v v=$nc '{print $1 + v, $2}' 7 > 8
                                                paste --delimiter='|' 6 8 > 9
                                                sed -i 's/$/|/' 9
                                                sed -i.bak '12s/.*//' $prefix.BSE-${line}$Kelvin.slpc.in
                                                sed -i.bak '122 r 9' $prefix.BSE-${line}$Kelvin.slpc.in
                                                sed -i.bak '122d' $prefix.BSE-${line}$Kelvin.slpc.in
                                                rm -rf 1 2 3 4 5 6 7 8 9  new_file.txt 1new_file.txt *.bak

                #----------------------Includes photon energy range, W screening bands, writes exciton wavefunctions--------------------#

                                                sed -i '134s/0.00000 | 10.00000 |/0.00000 | 7.00000 |/g'                	$prefix.BSE-${line}$Kelvin.slpc.in
                                                sed -i 's/BEnSteps= 100/BEnSteps= 1000/g'                                	$prefix.BSE-${line}$Kelvin.slpc.in
                                                sed -i 's/#WRbsWF/WRbsWF/g'                                             	$prefix.BSE-${line}$Kelvin.slpc.in


		#-------------------------SLEPC solver: first 5 excitons----------------------------------------------------------------#

                                                sed -i 's/BSSNEig=0/BSSNEig='$bssneig'/g'                               	$prefix.BSE-${line}$Kelvin.slpc.in
                                                sed -i 's/BSSEnTarget= 0.000000/BSSEnTarget= '$target'/g'                    	$prefix.BSE-${line}$Kelvin.slpc.in


                                                $MPIRUN_YAMBO yambo_ph -F $prefix.BSE-${line}$Kelvin.slpc.in -J output -C Report


                              	cd Report

                                               mv o-output.alpha_q1_slepc_bse $prefix.${line}$Kelvin.o-output.alpha_q1_slepc_bse

                                               scp -r $prefix.${line}$Kelvin.o-output.alpha_q1_slepc_bse $DIR/graph_data/BSE_Tempreature_data_slepc



       				cd ..

						ypp_ph -J output -e s 1

               					mv o-output.exc_qpt1_I_sorted $prefix.${line}$Kelvin.BSE_GW.o-output.exc_qpt1_I_sorted
               					mv o-output.exc_qpt1_E_sorted $prefix.${line}$Kelvin.BSE_GW.o-output.exc_qpt1_E_sorted
						scp -r $prefix.${line}$Kelvin.BSE_GW.o-output.exc_qpt1_E_sorted $DIR/graph_data/BSE_Tempreature_data_slepc
						scp -r $prefix.${line}$Kelvin.BSE_GW.o-output.exc_qpt1_I_sorted $DIR/graph_data/BSE_Tempreature_data_slepc

				cd ..
                        done


exit;
