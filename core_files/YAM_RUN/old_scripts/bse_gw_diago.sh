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

	mkdir BSE_GW_diago

	PATH3=$DIR/nscf
        PATH6=$DIR/GW_folder
	PATH14=$DIR/BSE_GW_diago

	cd $DIR/graph_data/

		mkdir BSE_GW_diago_data
	 


	cd $PATH6
                                                scp -r SAVE $PATH14

	
        cd $PATH14

#       # 1. BSE spectra using diagonalization solver is computed at different temperatures. 

	#-----------------------------------------------BSE + GW + phonon Calculation: diago solver------------------------------#
	        				yambo
                                                yambo -r -X s -optics b -kernel sex -Ksolver d -V all -F $bse_gw_filename

				
        #------------------------------COULOMB cut-off & RIM------------------------------#


                                                awk '/Alat factors/ {print $6-1}' r_setup > 1
                                                sed 's/$/|/' 1 > 2
                                                sed '1 s/./0.000000|0.000000|&/' 2 > 3
                                                sed '55s/.*//' $bse_gw_filename > 4
                                                sed '55 r 3' 4 > 5
                                                sed '55d' 5 > 6
                                                mv 6 $bse_gw_filename
                                                rm -rf 1 2 3 4 5

                                                sed -i 's/RandQpts=0/RandQpts=1000000/g'                               $bse_gw_filename
                                                sed -i 's/RandGvec= 1/RandGvec= 100/g'                                 $bse_gw_filename
                                                sed -i 's/CUTGeo= "none"/CUTGeo= "box z"/g'                            $bse_gw_filename

                                                sed -i 's/DBsIOoff= "none"/DBsIOoff= "BS"/g'                           $bse_gw_filename
                                                sed -i 's/DBsFRAGpm= "none"/DBsFRAGpm= "+BS"/g'                        $bse_gw_filename



                #------------------------------CPU parallel structure-----------------------------#


                                                sed -i 's/X_and_IO_CPU= ""/X_and_IO_CPU="'"$x_and_io"'"/g'              $bse_gw_filename
                                                sed -i 's/X_and_IO_ROLEs= ""/X_and_IO_ROLEs= "q g k c v"/g'             $bse_gw_filename

                                                sed -i 's/DIP_CPU= ""/DIP_CPU="'"$dip_cpu"'"/g'                        	$bse_gw_filename
                                                sed -i 's/DIP_ROLEs= ""/DIP_ROLEs= "k c v"/g'                           $bse_gw_filename

                                                sed -i 's/BS_CPU= ""/BS_CPU= "'"$bs_cpu"'"/g'                           $bse_gw_filename
                                                sed -i 's/BS_ROLEs= ""/BS_ROLEs= "k eh t"/g'                            $bse_gw_filename

               #---------------------------FFT-GVecs and other Response size block---------------#

                                                sed -i.bak -e '30d'                                                     $bse_gw_filename
                                                sed -i.bak "30i FFTGvecs= $fftgvecs"                                    $bse_gw_filename
                                                
                                                sed -i.bak -e '164d'                                                    $bse_gw_filename
                                                sed -i.bak "164i NGsBlkXs= $ngsblks"                                    $bse_gw_filename

                                                sed -i.bak -e '81d'                                                     $bse_gw_filename
                                                sed -i.bak "81i BSENGBlk= $bsengblk"                                    $bse_gw_filename

                                                sed -i.bak -e '79d'                                                     $bse_gw_filename
                                                sed -i.bak "79i BSENGexx= $bsengexx"                                    $bse_gw_filename



                #----------------------------In-plane electric field directions-------------------#             

                                                sed -i '136s/ 1.000000 | 0.000000 | 0.000000 |/ 1.000000 | 1.000000 | 0.000000 |/g'    	$bse_gw_filename
                                                sed -i '174s/ 1.000000 | 0.000000 | 0.000000 |/ 1.000000 | 1.000000 | 0.000000 |/g'     $bse_gw_filename

                #----------------------------Includes el-hole coupling----------------------------#

                                                sed -i 's/BSEmod= "resonant"/BSEmod= "resonant"/g'                      $bse_gw_filename
                                                sed -i 's/#WehCpl/WehCpl/g'                                             $bse_gw_filename
						sed -i '131s/0.100000 | 0.100000 |/0.070000 | 0.070000 |/g'             $bse_gw_filename

		#---------------------Includes GW + temperature dependent band gap energies--------#

                                                sed -i 's/KfnQPdb= "none"/KfnQPdb= "E < ..\/GW_folder\/output\/ndb.QP"/g'       $bse_gw_filename

               	#------------------------Reading the k point limits, max valence band, min conduction band from r_setup file-------------#

                                                awk '/Filled Bands/ {print $5}' r_setup > 5
                                                awk -v v=$nv '{print $1 - v, $2}' 5 > 6


                                                awk '/Empty Bands/ {print $5}' r_setup > 7
                                                awk -v v=$nc '{print $1 + v, $2}' 7 > 8
                                                paste --delimiter='|' 6 8 > 9
                                                sed -i 's/$/|/' 9
                                                sed -i.bak '116s/.*//' $bse_gw_filename
                                                sed -i.bak '116 r 9' $bse_gw_filename
                                                sed -i.bak '116d' $bse_gw_filename
                                                rm -rf 1 2 3 4 5 6 7 8 9  new_file.txt 1new_file.txt *.bak

                #----------------------Includes photon energy range, W screening bands, writes exciton wavefunctions--------------------#

                                                sed -i '128s/0.00000 | 10.00000 |/0.00000 | 7.00000 |/g'                $bse_gw_filename
                                                sed -i 's/BEnSteps= 100/BEnSteps= 500/g'                                $bse_gw_filename
                                                sed -i 's/#WRbsWF/WRbsWF/g'                                             $bse_gw_filename



                                                $MPIRUN_YAMBO yambo -F $bse_gw_filename -J output -C Report


                              	cd Report

                                               	mv o-output.alpha_q1_diago_bse $prefix.BSE_GW.o-output.alpha_q1_diago_bse
                                               	mv o-output.eps_q1_diago_bse $prefix.BSE_GW.o-output.eps_q1_diago_bse
                                              	mv o-output.jdos_q1_diago_bse $prefix.BSE_GW.o-output.jdos_q1_diago_bse

						scp -r $prefix.BSE_GW.o-output.alpha_q1_diago_bse $prefix.BSE_GW.o-output.jdos_q1_diago_bse $DIR/graph_data/BSE_GW_diago_data


       				cd ..

						ypp -J output -e s 1

               					mv o-output.exc_qpt1_I_sorted $prefix.BSE_GW.diago.o-output.exc_qpt1_I_sorted
               					mv o-output.exc_qpt1_E_sorted $prefix.BSE_GW.diago.o-output.exc_qpt1_E_sorted
						scp -r $prefix.BSE_GW.diago.o-output.exc_qpt1_E_sorted $DIR/graph_data/BSE_GW_diago_data/.
						scp -r $prefix.BSE_GW.diago.o-output.exc_qpt1_I_sorted $DIR/graph_data/BSE_GW_diago_data/.


exit;
