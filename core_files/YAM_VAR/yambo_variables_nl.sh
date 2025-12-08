#!/bin/bash

#------------------------------------- This works for version Yambo 5.3.0 ----------------------------------------#

#--------------- Variables of YAMBO Calculations from this section onwards--------#
#---------------- Yambo & CPU Parallel structure ---------- #

PARALLEL_YAMBO="mpirun"
MPIRUN_YAMBO="$PARALLEL_YAMBO -np 8"
YAMBO_PH="/Users/sitangshubhattacharya/Softwares/yambo-5.3.0/bin/yambo_ph"
YPP_PH="/Users/sitangshubhattacharya/Softwares/yambo-5.3.0/bin/ypp_ph"
P2Y="p2y"
YAMBO="yambo"
YPP="ypp"
YPP_NL="ypp_nl"
YAMBO_NL="yambo_nl"

nl_cpu_rt_bse="1 8"			# linear response BSE: For max performance: Put max cpus in k and 1 for w.
nl_cpu="8 1"				# SHG, THG: For max performance: Put max cpus in w and 1 for k.
dip_cpu="1 4 2"
oscll_cpus="4 2"
rt_cpu="1 2 2 2"
x_and_io="1 1 1 4 2"
#---------------------------------------------------#


#---------------- Nonliner RIM cut-off [box z only] ---------#
#Note: Yambo 5.3.0: slab z dont work on nl optics yet.

slab_nl="box z"                    #   box/cylinder/sphere/ws/slab X/Y/Z/XY..


# Common parameters for nonlinear calculations: real-time BSE and SHG:
#---------------------------------------------------#
correlation="SEX"		    #	This is HXC correlation type: Make it "h" or "SEX" or "h+SEX", etc. See yambo_nl -h -potential.
NLtime="55"                         #   This is total simulation time in "fs" used in Cranck-Nicholson iterative scheme. Keep this in 55-60 for small systems.
nl_integrator="CRANKNIC"            #   This is time step solver for TD-DFT equation using crank-nicholson algorithm. Time step used=0.01 fs.
nl_Correlation="SEX"                #   This is type of electron-hole correlation. SEX is equivalent to BSE.
min_NLEnRange="0.1"                 #   Min Laser frequency [x-axis energy range max value in eV] for SHG calculations or choose other in numbers. Note: This should be in half (w/2) eV.
max_NLEnRange="$phot_max_en_lim"  #   	Max Laser frequency [x-axis energy range max value in eV] for SHG calculations or choose other in numbers. Note: This should be in half (w/2) eV.

# Parameters for real-time BSE:
#---------------------------------------------------#
nl_Damping="0.100000"               #   This is the damping or the dephasing in "eV". It is required in both real-time BSE and SHG.

# Parameters for nonlinear SHG:
#---------------------------------------------------#
nl_ensteps="80"                    #   For SHG only. Value of energy steps [laser frequency w] for SHG calculations. A suitable value: 120.

#---------------------------------------------#

