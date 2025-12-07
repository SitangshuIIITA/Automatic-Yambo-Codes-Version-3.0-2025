#!/bin/bash

#--------1. Working & Pseudo QE bin Directory Location-------------------------------#
currentlocation=$(pwd)
DIR=${currentlocation%/*\/*}/work_dir
PSEUDO_DIR=$DIR/pseudo
#---------------------------------------------------------------------------------#

#---------2. HPC-CPU parallel distribution and prefix------------------------------#

QE_PATH=/home/sitangshu/Software/q-e-EPW-6.0/bin
#QE_PATH=/home/sitangshu/Software/q-e-qe-6.7MaX-Release/bin

PARALLEL="mpirun"

  MPIRUN_PW="$PARALLEL -np 4  $QE_PATH/pw.x -ndiag 1 -inp"
MPIRUN_nscf="$PARALLEL -np 4 $QE_PATH/pw.x -inp"
MPIRUN_Bnds="$PARALLEL -np 4 $QE_PATH/pw.x  -inp"
MPIRUN_Bndx="$PARALLEL -np 4 $QE_PATH/bands.x -inp"
MPIRUN_proj="$PARALLEL -np 4 $QE_PATH/projwfc.x -inp"

#For Thermo_pw calculations:
MPIRUN_THERMO_PW="$PARALLEL -np 4  $QE_PATH/thermo_pw.x -npool 2"

#For EPW calculations:
MPIRUN_EPW="$PARALLEL -np 4  $EPW_PATH/epw.x -npool 4"



#For phonon dispersion and dos

PH_DISP_RUN="$PARALLEL -np 4 $QE_PATH/ph.x -npool 2"  #k-point parallelization is supported.
MATDYN_DISP_RUN="$PARALLEL -np 8 $QE_PATH/matdyn.x"
Q2R_DISP_RUN="$PARALLEL -np 8 $QE_PATH/q2r.x"


#For electron-phonon:
PW_RUN="$PARALLEL -np 4 $QE_PATH/pw.x"
PH_RUN_phon="$PARALLEL -np 4 $QE_PATH/ph.x -npool 2"	#k-point (npool) parallelization is supported.
PH_RUN_dvscf="$PARALLEL -np 4 $QE_PATH/ph.x"		#k-point (npool) parallelization is not supported.
MATDYN_RUN="$PARALLEL -np 4 $QE_PATH/matdyn.x"
Q2R_RUN="$PARALLEL -np 4 $QE_PATH/q2r.x"
 

#----------3. DFT input prefix----------------------------------------------------#
prefix=bn                           # Use your prefix
forc_conv_thr=1e-05
etot_conv_thr=1e-05
ibrav=4
ecutwfc=60 	                        #Kinetic cut-off for upf pseudo.
celldm_1=4.716
celldm_3=7
nat=2
ntyp=2
occupation=fixed
force_symmorphic=.FALSE.
lspinorb=.FALSE.
noncolin=.FALSE.

# Truncation of the Coulomb interaction in the z direction for structures periodic in the x-y plane.
#Total energy, forces and stresses are computed in a two-dimensional framework:
assume_isolated="2D"                # "2D" or "none".
cell_dofree="2Dxy"                  #  only x and y components are allowed to change

kgrid_nonsoc="12 12 1"                #   kpoint grid in absence of SOC for nonsoc-scf, phonon dispersion & elph.
kgrid="12 12 1"	                    #   kpoint grid- for nscf, GW, BSE, etc.
shifted_grid="0 0 0"

#----------4. Number of atoms & positions----------------------------------------#
# Don't put extra commas or inverted commas, etc.

atoms="B N"                         #   Write the atomic species to be used in electronic stability calculation.
mass_atom="10.811000 14.0067"       #   Atomic mass of atom. Put them in the order of atoms.
pseudo_prefix="B N"   #   pseudopotential prefix name. Put them in the order of atoms. The extension is already in UPF.

# Atomic positions in crystal coordinates
positions=(
  "B   0.66666667  0.33333333 0.00000000"
  "N  -0.66666667 -0.33333333 0.00000000"
)

#----------5. Phonon dispersion q-grid--------------------------------------------#

NQ1="12"                         #Give these nq1, nq2, nq3 same as in kgrid_nonsoc for phonon dispersion.
NQ2="12"
NQ3="1"
TRN_PH="1.0d-16"
alpha_mx="0.7"


#----------6. Double-grid k points: to be used in BSE scissor---------------------#

dgkpoint="8 8 1 1 1 0" 	        #Put as double the k-grid on shifted 1 1 0


#----------7. nscf no. of bands (occupied+unoccupied)-----------------------------#
NBND="                                              nbnd = 20"
NBND_ELPH=20			        #This nbnd is used in elph calculations.
#---------------------------------------------#



#----------8. Choice of BZ route for band diagram---------------------------------#
bands_path="G M K G"
#-----------------End of DFT Script-----------------------------------------------#

