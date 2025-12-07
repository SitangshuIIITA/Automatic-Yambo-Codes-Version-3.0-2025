#!/bin/bash






# Note: What to run: Make sure you installed them. 

#	1. You need two QEs: 6.8 [q-e] and 7.5 [qe-epw-6.0]. Both patched.
#	2. To install patched q-e [QE-6.8] use: https://github.com/janberges/qe2screen
#	3. To install patched q-e-epw-6.0 [QE-7.5], use: https://github.com/janberges/elphmod

# AFter installation, go to core_files/EPW_VAR and carefully make changes in EPW_variables.sh.

# 1. Run this script run_EPW.sh by uncommenting (one by one or all of them) the lines.
#Suggestion: Before initiating Wannierization and EPW, run a converged QE relax, scf, nscf, bands and orbital projections first.

cd ../core_files/QE_VAR
set -a
source ./QE_variables.sh
set +a

cd ../../work_dir

cd ../core_files/EPW_VAR
set -a
source ./EPW_variables.sh
set +a

cd ../../work_dir


cd ../core_files/EPW_RUN 

# Run sequentially or at one shot, the following lines :

#---------------A0. W90 based tasks [Quantum Espresso] Converge these first after running run_QE_Yambo.sh :---------------------#
# Serial used: patched QE-7.5

#		bash w90_relax.sh		 	# calculates relax.

#		bash w90_scf.sh			 	# calculates scf.
		
#		bash w90_phon_disp.sh			# calculates phonon dispersion and dvscfs.
		
#		bash w90_nscf.sh			# calculates nscf by given nband by breaking crystal symmetry.
		
#		bash wan90.sh				# generates and compute the wan90 input and wannier files for convergence study and subsequent epw calculation.


#---------------A1. EPW based tasks--------------------------------------------------------------------------------------------#
# To Run a succesful EPW:
#  1. Converge first W90 in previous steps using suitable disentanglement windows. Edit core_files/EPW_VAR/EPW_variables.sh to set them.
#  2. Run first the above scripts. This will set dvscf + phonon dispersion.


# Serial used: patched QE-7.5

#		bash epw_prelim.sh			# This generates and compute the EPW input and wannier electron and phonon dispersions for subsequent mobility, etc. calculation.

# Serial used: patched QE-6.8

#		bash quadrupole.sh			# This generates quadrupole.fmt by fiiting bare DFPT results.

# Serial used: patched QE-7.5

#		bash lr_I.sh				# This computes and prepares the gkkp folder for 'no_lr' '3d' 'gaussian' 'dipole_sp' 'quadrupole' conditions along KGM route. 

#		bash gkkp_matrix.sh			# This generates and compute the electron/hole inter/intraband matrix elements at CBM/VBM along q-BZ from the previous run.

#		bash plot_matrix.sh			# This computes and plots the no_lr and quadrupole electronic gkkp [intra and interband] along ph-BZ route "KGM" for k = CBM. See EPW_variable list to use desired CBM and VBM ith and jth band numbers.


#		bash gkkp.sh				# This generates and compute the electron/hole inter/intraband matrix elements along all k and q-BZ from the previous runs.


#		bash self_e.sh				# This compute the electronic self energy [Im \Sigma] for the conduction bands KGM route at various temperatures: 77, 100, 300, 500 K.


#		bash plot_self_e_KGM.sh			# Check the self energy plots in EPW_plots folder.


#		bash self_h.sh				# This compute the hole self energy [Im \Sigma] for the valence bands KGM route at various temperatures: 77, 100, 300, 500 K.


#		bash plot_self_h_KGM.sh			# Check the self energy plots in EPW_plots folder.


#		bash e_mobility.sh			# Finally! this computes the electronic mobility [both dipole_sp and quadrupole] at various temperatures specified by ntemps and a fixed carrier density [ncarrier]. Check the folder EPW_plots to plot the data.

#		bash h_mobility.sh			# Finally! this computes the hole mobility [both dipole_sp and quadrupole] at various temperatures specified by ntemps and a fixed carrier density [ncarrier]. Check the folder EPW_plots to plot the data. NOTE: You do not have to make negative of ncarrier. Its taken care of.

#		bash e_mobility_vs_iicarrier.sh		# This computes the electronic mobility [quadrupole] vs ionized impurity densities at 300 K. Check the folder EPW_plots to plot the data. The total mobility has both mobility_ph + mobility_ionized_imp values.

#		bash h_mobility_vs_iicarrier.sh		# This computes the hole mobility [quadrupole] vs ionized impurity densities at 300 K. Check the folder EPW_plots to plot the data. The total mobility has both mobility_ph + mobility_ionized_imp values.

#		bash extract_mobility.sh		# This extracts the data files to plot. Check EPW_plots folder. Change the variable in gaussian-h.py inside core_files/EPW_RUN for spectral decomposition plots. 

#		bash electron_monte_carlo.sh		# This extracts Einstien electronic mobility and diffucion coefficient using zero field algorithm of MonteCarlo.
		
#		bash hole_monte_carlo.sh		# This extracts Einstien hole mobility and diffucion coefficient using zero field algorithm of MonteCarlo.


		bash electron_drift.sh			# This extracts electronic drift velocity vs electric field [along x direction] algorithm of MonteCarlo.



######################################################  GAME OVER ###############################################################

	
rm -rf count

exit;
