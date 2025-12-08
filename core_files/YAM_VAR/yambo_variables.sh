#!/bin/bash

#------------------------------------- This works for version Yambo 5.3.0 ----------------------------------------#

#--------------- Variables of YAMBO Calculations from this section onwards--------#
#---------------- Yambo & CPU Parallel structure ---------- #

PARALLEL_YAMBO="mpirun"
MPIRUN_YAMBO="$PARALLEL_YAMBO -np 8"
YAMBO_PH="/home/sitangshu/Software/yambo-5.3.0/bin/yambo_ph"
YPP_PH="/home/sitangshu/Software/yambo-5.3.0/bin/ypp_ph"
P2Y="p2y"
YAMBO="yambo"
YPP="ypp"
YPP_NL="ypp_nl"
YAMBO_NL="yambo_nl"

se_cpu="1 4 2"
x_and_io="1 1 1 4 2"
dip_cpu="1 4 2"
bs_cpu="1 4 2"
#---------------------------------------------------#

#---------------- Hartree-Fock Convergence ---------#
exx_value="10 20 30"	        #   This value is both exxrlvecs and vxcrlvecs and is converged at the Hartree-Fock Level.
                                #   Convergence is obtained by varying both EXXRLvcs = VXCRLvcs = "10 30 50" at the samte time.
ryd="Ry"
rl="RL"
rimw_type="semiconductor"	    #   This is RIM_W type: default/semiconductor/metal/graphene
randqpts="3000000"
randgvec="100"		            #   This is RIM value in RL.
randgvecw="5"		            #   This is RIM_W value in RL.
slab="slab z"		            #   box/cylinder/sphere/ws/slab X/Y/Z/XY..

#---------------------------------------------------#

# Converge the followings for GW:
#---------------------------------------------------#

# Note: both ngsblpp_value and bndrnge must be of same array size. Note also: The total number of arrays must be same and should exist with the highest number of valence bands. This means if array size is four, then at least four top valence bands should exist in both GW and BSE.


ngsblpp_value="7 9 10"           #   This is response block size NGsBlkXp in Ry. Converge this value here by increasing such as 1 5 9 12, etc.
bndrnge="10 12 15"               #   Converge this screening bands. This is [Xp] Polarization function bands. See line 29 in QE_variables. It should be less than equal to NBND.
                                #   Convergence is obtained by varying both BndsRnXp = GbndRnge = "5 10 14 17"  at the samte time (four entities in this example). You can increase more.
                        
gw_self_inter="3"               #   GW self consistent iteration number: 3 means G0W0 to G1W1 to G2W2 to G3W3.

#-------------------------------------------------- #

# Converge the followings for BSE:
#---------------------------------------------------#
BSE_gap="gw_selfC"                  #   This is type and location of converged QP gap. Use it as "g0w0" if G0W0, or "gw_selfC" if self-consistent GW. Run the previouis gw_converge and gw_selfC first.
ngsblks="$ngsblpp_value"            #   Converge this value here. This is response block size. Keep this same as ngsblkp.
bsengblk="$ngsblpp_value"           #   Converge this value here. This is response block size. Keep this same as ngsblks.
bsengexx="$exx_value"               #   Converge this value here. This is response block size. Keep this same as ngsblks.

lkind="BAR"                        #   full or BAR in BSE kernel.
Bdmrange="0.1"                      #   This  is BDmRange when fixed atom condition is used. BSS damping range: Increase this or decrease to get experimental broadening of the spectrum around each excitons.
phot_max_en_lim="9.0000"            #   This is BEnRange. Max. Photon energy (in eV) in the absorption spectrum.
mode_electric_field="x"            #   This is LongDrXs and BLongDir. Put it either x, or y or z or xy or yz or zx. For 2D: it is xy.

#-----------GW/BSE Choice of transition bands-------#
nc="2"                              #   This is the total no. of conduction bands ABOVE THE CBM. If nc = 0, then the code takes only the lowest CB.
nv="2"                              #   This is the total no. of valence bands BELOW THE VBM. If nv = 0, then the code takes only the highest VB.
#-------------------------------------------------- #


#-----------Exciton PDF 3D Xcrysden plotting at Q=0 [optical dipole limit]---------------#
exciton_states="1 - 3"             #   These are the 10 lowest exciton states whose dispersion diagram will be plotted along the BZ. MIND THE GAPS BETWEEN NUMBERS.
at_coord_x="4.67428"                #   Hole position (x) in unit cell (positive) in crystal coordinate. Look into the *.E-sorted file.
at_coord_y="4.67428"                #   Hole position (y) in unit cell (positive) in crystal coordinate. Look into the *.E-sorted file.
at_coord_z="4.67428"                #   Hole position (z) in unit cell (positive) in crystal coordinate. Look into the *.E-sorted file.
#-------------------------------------------------- #

#-----------Exciton Finite Q dispersion states---------------#
bssneig="3"                        #  This is BSSNEig [Number of eigenvalues to compute] and is the total number of lowest excitons which will be evaluated in BSE/exciton dispersion using SlepC. See exciton_states above.
target="5.00"		                #  This is BSSEntarget. This value is in eV and the SlepC solver look around this number to find the excitons.
#-------------------------------------------------- #

#-----------Exciton FInite-Q dispersion states---------------#
QPT="1 4 6"                         #   This is BSE Finite Q (momenta) point: BSQindex [Q-Index of the BS state(s)].: Required to sketch exciton 3D PDF xsf plots at these points.
#-------------------------------------------------- #

#------------- ELPH double grid convergence: Convergence of el-ph QP energies on phonon q double grid:---------------#
matd="18 24 32"              #   These are the matdyn phonon grids: It is read like 18x18x1, 24x24x1, 32x32x1, 60x60x1, 120x120x1, etc.
GphBRnge="$NBND_ELPH"               #   This should be less or equal to nbnd in nscf and is the [ELPH] G[W] bands range.
#----------------------------------------------------#


#-----Temperature file names------------------------#
temp_value="0 300"
Kelvin="K"
GDmRnge="0.2"                       # This is G_gw damping for spectral function in eV.

#----------------------------------------------------#

#-----Exciton Phonon Luminiousity-------------------#
dos_broad_exc_phon=0.01                                # The broadening is in eV
dos_energy_range="5.000000 |5.200000 |         eV"     #Energy steps
#----------------------------------------------------#


#-----Exciton dispersion BZ path---------------------#
QPATH="G M K"                     # Excitonic BZ: Allowed values: G, M, K — separated by spaces. Path to exciton dispersion plot. You can chose any order: G K G or M G K, etc.

#----------------------------------------------------#


#------------ File names --------------------#
hf_filename=$prefix.yambo_hf_input.in
hf_ypp_bands=$prefix.hf_bands.in
gw_filename=$prefix.yambo_GW_input.in
gw_ypp_bands=$prefix.GW_bands.in
bse_filename=$prefix.yambo_bse_input.in
bse_slepC_filename=$prefix.slepC.yambo_bse_input.in
bse_gw_filename=$prefix.yambo_bse_diago_input.in
sci_filename=$prefix.yambo_bse_scissor_input.in
ex_prefix=$prefix.bse_scissor
exc_disp_filename=$prefix.yambo_exc_disp_input.in
ypp_exciton_disp=$prefix.ypp_exciton_disp.in
#---------------------------------------------#


#-----Temperature file names------------------#
temp_prefix="$temp_value Kn"
bse_temp_filename=$prefix.yambo_bse.$temp_value$Kelvin.input.in
exprefix0K=$prefix.$temp_value$Kelvin
exc_disp_temp=$prefix.yambo_bse_exc_disp.$temp_value$Kelvin.input.in
#---------------------------------------------#

