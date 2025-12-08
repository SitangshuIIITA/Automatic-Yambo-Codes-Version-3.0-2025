#!/bin/bash

#--------1. Working & Pseudo QE bin Directory Location-------------------------------#
currentlocation=$(pwd)
DIR=${currentlocation%/*\/*}/work_dir
PSEUDO_DIR=$DIR/pseudo
#---------------------------------------------------------------------------------#

#---------2. HPC-CPU parallel distribution and prefix------------------------------#

EPW_PATH=/home/sitangshu/Software/q-e-EPW-6.0/bin
EPW_PATH_py=/home/sitangshu/Software/q-e-EPW-6.0/EPW/bin
PERL_PATH=/home/sitangshu/Software/q-e-EPW-6.0/external/wannier90/utility/kmesh.pl


QUAD_EPW_PATH=/home/sitangshu/Software/q-e/bin
QUAD_EPW_PATH_py=/home/sitangshu/Software/q-e/EPW/bin

PARALLEL="mpirun"

  MPIRUN_epw_PW="$PARALLEL -np 4  $EPW_PATH/pw.x -ndiag 1 -inp"
MPIRUN_epw_nscf="$PARALLEL -np 4 $EPW_PATH/pw.x -inp"
MPIRUN_epw_Bnds="$PARALLEL -np 4 $EPW_PATH/pw.x  -inp"
MPIRUN_epw_Bndx="$PARALLEL -np 4 $EPW_PATH/bands.x -inp"
MPIRUN_epw_proj="$PARALLEL -np 4 $EPW_PATH/projwfc.x -inp"

#For phonon dispersion and dos

EPW_PH_DISP_RUN="$PARALLEL -np 4 $EPW_PATH/ph.x -nk 4"  #k-point parallelization is supported.
EPW_MATDYN_DISP_RUN="$PARALLEL -np 4 $EPW_PATH/matdyn.x"
EPW_Q2R_DISP_RUN="$PARALLEL -np 4 $EPW_PATH/q2r.x"


#For EPW calculations:
MPIRUN_wan="$EPW_PATH/wannier90.x"
MPIRUN_pw2wan="$EPW_PATH/pw2wannier90.x"
MPIRUN_EPW="$PARALLEL -np 4 $EPW_PATH/epw.x -npool 4"





#--------Quadrupole-EPW---------------------------------------------------------#

  QUAD_MPIRUN_epw_PW="$PARALLEL -np 4 $QUAD_EPW_PATH/pw.x -ndiag 1 -inp"
QUAD_MPIRUN_epw_nscf="$PARALLEL -np 4 $QUAD_EPW_PATH/pw.x -inp"
QUAD_MPIRUN_epw_Bnds="$PARALLEL -np 4 $QUAD_EPW_PATH/pw.x  -inp"
QUAD_MPIRUN_epw_Bndx="$PARALLEL -np 4 $QUAD_EPW_PATH/bands.x -inp"
QUAD_MPIRUN_epw_proj="$PARALLEL -np 4 $QUAD_EPW_PATH/projwfc.x -inp"

#For phonon dispersion and dos

QUAD_EPW_PH_DISP_RUN="$PARALLEL -np 4 $QUAD_EPW_PATH/ph.x -nk 4"  #k-point parallelization is supported.
QUAD_EPW_MATDYN_DISP_RUN="$PARALLEL -np 4 $QUAD_EPW_PATH/matdyn.x"
QUAD_EPW_Q2R_DISP_RUN="$PARALLEL -np 4 $QUAD_EPW_PATH/q2r.x"


#For EPW calculations:
QUAD_MPIRUN_wan="$QUAD_EPW_PATH/wannier90.x"
QUAD_MPIRUN_pw2wan="$QUAD_EPW_PATH/pw2wannier90.x"
QUAD_MPIRUN_EPW="$PARALLEL -np 4 $QUAD_EPW_PATH/epw.x -npool 4"




#----------1. W90 input prefix----------------------------------------------------#

# Orbital selections:
orbitals=(
  "B: s; pz"
  "N: px; py; pz"
)

num_bands=20		# Total number of bands in nscf calculation: MUST.
num_wann=5		# Total number of orbitals (WF) that make the VB and CB bands.
dis_num_iter=30000
num_iter=30000		# Iteration criteria for DLTA convergence. Make it large so that DLTA~10^-15.	
iprint=2 		# Verbosity level of Wannier90 code.
num_dump_cycles=5	# Steps of iterations
num_print_cycles=5	# Print steps of iterations
conv_tol=1E-14		# W90 convergence criteria for iteration.



#-----Converge these values:
dis_win_min=-15.5	# Min outer window for W90.	
dis_win_max=10.6		# Max outer window for W90.
dis_froz_min=-15.5	# Min inner window for W90.
dis_froz_max=-1.0	# Max inner window for W90.

#----------2. EPW input prefix----------------------------------------------------#

epbwrite=.true.		# If epbwrite = .true., the electron-phonon matrix elements in the coarse Bloch representation and relevant data (dyn matrices) are written to disk. If epbread = .true. the above quantities are read from the ‘prefix.epb’ files. Pool dependent files.
epbread=.false.

epwwrite=.true.		# If epwread = .true., the electron-phonon matrix elements in the coarse Wannier representation are read from the ‘epwdata.fmt’ and ‘XX.epmatwpX’ files. Each pool reads the same file. It is used for a restart calculation and requires kmaps = .true. A prior calculation with epwwrite = .true is also required.
epwread=.false.

etf_mem=1		# If etf_mem = 1, then all the fine Bloch-space el-ph matrix elements are stored in memory (faster). When etf_mem = 1, more IO (slower) but less memory is required. When etf_mem = 2, an additional loop is done on mode for the fine grid interpolation part. This reduces the memory further by a factor “nmodes”. The etf_mem = 3 is like etf_mem = 1 but the fine k-point grid is generated with points within the fsthick window using k-point symmetry (mp_mesh_k = .true. is needed) and the fine q-grid is generated on the fly. The option etf_mem = 3 is used for transport calculations with ultra dense fine momentum grids.

wannierize=.true.	# Calculate the Wannier functions using W90 library calls and write rotation matrix to file ‘filukk’. If .false., filukk is read from disk.

lpolar=.true.		# If .true. enable the correct Wannier interpolation in the case of polar material.

vme=dipole		# if ‘dipole’ then computes the velocity as dipole+commutator = <psi_{mk} |p+i[V_{NL},r]| psi_{nk}>`. If ‘wannier’ then computes the velocity as dH_nmk/dk - i(e_nk-e_mk)A_nmk where A is the Berry connection. Note: Before v5.4, vme = .FALSE. was the velocity in the local approximation as <psi_mk|p|psi_nk>. Before v5.4, vme = .TRUE. was the same as ‘wannier’.
 
exclude_bands=1,7:20	# Number of bands you dont want to compute scatttering for, example: 1-3,6-20

 
degaussw=0.005 	# Smearing for sum over q in the e-ph coupling in [meV]

#-----Fine k and q grids for interpolations: Converge these by increasing them.
nkf1=80
nkf2=80
nkf3=1
nqf1=80
nqf2=80
nqf3=1

system_2d=dipole_sp
scissor=1.88				# This is GW gap minus DFT gap. Run GW for this.


#----------Mobility for the electron-phonon part only -----------------------#


nstemp=10 						# Number of temperatures
temps="77 100 150 200 250 300 350 400 450 500"		# Temperatures for mobility calculation [Careful: no space between and no 0 temperature].
ncarrier=1E9						# Carrier concentration in (cm^-2) in 2D. If positive [ncarrier=1E10] then electron mobility, if negative [ncarrier=-1E10] then hole mobility.

#----------------------------------------------------------------------------#



#----Mobility for the ionized impurity + grain Boundary Scattering part -----#



ii_imp="1.0d9 1.0d10 1.0d11 1.0d12 1.0d13 1.0d14" # Inoized impurities concentration for mobilities in cm^-2. 
gb_size=100				# Size of the grains in grain-boundary dominated mobilities.
t_i=77					# Temperature (in K) at which mobilities Vs ionized impurities + grain boundaries will be carried out.

#----------------------------------------------------------------------------#





fsthick=1.5				# Width of the Fermi surface window to take into account states in the self-energy delta functions in [eV]. Narrowing this value reduces the number of bands included in the selfenergy calculations.

fermi_shift_cond=0.15			# If "0 eV", the it is Ec=Ef. If positive, then Ef enters in CB. 	
fermi_shift_vale=0.15			# If "0 eV", the it is Ev=Ef. If positive, then Ef enters in VB.

ibndprt=2		# DFPT band number at k=0 [which electronic band you want to print]: Used in Quadrupole part of the reference DFPT computation. NOte Q do not depend on n,m, k. ibndprt could be one of the deep valence band [wannier90 is more localized here].


n_band=4   # This is initial band number. Requierd in |g_nm,k=\Gamma| plot for no_lr and dipole+Quadrupole comparison plot. Chooose the band indices properly. These should be within "exclude_bands". 5 here is CBM. If n_band=m_band, its intraband.
m_band=4   # This is final band number. Requierd in |g_nm,k=\Gamma| plot for no_lr and dipole+Quadrupole comparison plot. Chooose the band indices properly. These should be within "exclude_bands". 5 here is CBM. If n_band=m_band, its intraband.




























