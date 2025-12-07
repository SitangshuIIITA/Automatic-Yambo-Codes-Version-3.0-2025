#!/bin/bash






# Note: How to run:

# 1. Run this script run.sh by uncommenting input_relax.sh line below to generate the input script. Once the input is generated, comment it and go to the QE and Yambo shell script lines below. Uncomment (one by one or all of them) them and run this script again:
# 2. In the Yambo BSE command runs, SlepC module for yambo is required. Install them if you dont have.

cd ../core_files/QE_VAR
set -a
source ./QE_variables.sh
set +a

cd ../../work_dir

#---------------A0 DFT Input File Generation [Quantum Espresso] :------------------------------------------#
	       	bash input_relax.sh                                       # Relax input file generation

cd ../core_files/QE_VAR

cd ../YAM_VAR
set -a
source ./yambo_variables.sh
set +a

set -a
source ./yambo_variables_nl.sh
set +a


cd ../QE_RUN 

# Run sequentially or at one shot, the following lines :

#---------------A1. DFT based tasks [Quantum Espresso] :---------------------------------------------------#

#		bash QE_relax.sh				# calculates energy/geometry minimization/relax.

#		bash QE_electron_stability.sh			# calculates formation energy of the system.

#		bash QE_scf.sh					# calculates scf (both with spin and no-spin).

#		bash QE_nscf.sh					# calculates non self consistent part.

#		bash QE_bands.sh				# calculates bands using QE band interpolator.
  
 # 		bash QE_bands_projection.sh           		#calculates bands ldos/pdos/. etc with spin.

#		bash QE_double_nscf.sh				# calculates nscf on denser grid.

#		bash QE_phon_disp.sh				# calculates phonon dispersion on regilar grid. Won't run unless bands are done properly.

#		bash QE_thermo_pw.sh                  		# calculates vibrational zero-point energies under QHA using thermo_pw.x (lattice thermal expansion).





	
cd ../YAM_RUN

# The following tasks depend on the previous QE folders. So, run all of them first.

# Linear properties: Run sequentially or at one shot, the following lines:

#	       	bash yambo_dft_bnd_struc.sh			# calculates DFT band structure using yambo interpolator.


#---------------A2. CPU intense tasks [GW+BSE]:----------------------------------------------------------#

#	       	bash hf_conv.sh					# Hartree-Fock energy convergence computation. The convergences are for exx and vxx.
        
#	       	bash gw_conv.sh                      		# G0W0 energy convergence computation. The convergences are for nglx, and G_band and .
        
#	       	bash gw_selfC.sh                     		# Self_consistent GW computation: G0W0, G1W1, G2W2,.. Conduction and valence slopes using regression method is also computed.
        
#	       	bash bse_conv.sh                     		# GW-selfC (or G0W0) + BSE DIAGONALIZATION SOLVER with automatic convergence on transition bands, bsengexx, ngsblks, bsengblk, and transition bands at same time.

#		bash bse_bands_conv.sh
        
#        	bash bse_slepC.sh                    		# GW-selfC (or G0W0) + BSE finite Q calculation using DIAGONALIZATION SOLVER and SLEPC interation on a converged BSE parameters.


#---------------A3. CPU intense tasks [Temperature-Dependent el-ph + BSE]:-------------------------------#

#        	bash elph_dvscf.sh                   		# Computation of perturbed potential [dV(q)/d\xi] and phonon dynamical matrix elements for el-ph interaction.
  
#        	bash elph_yambo.sh                   		# Computation of initial and scattered electronic states |nk> and |n'k+q>.
        
#        	bash elph_converg_q.sh               		# convergence of el-ph QP energies [Fan+Debye Waller] at T = 0 K on various dense phonon double grids.

#        	bash qp_energies.sh                  		# QP [Fan + Debye Waller] ZPR energies at various temperatures: This gives band-gap Vs Temperature data.

#        	bash qp_spectral.sh                  		# QP spectral functions at various temperatures.

#        	bash qp_eliashberg.sh                		# QP electronic Eliashberg function g^F(\omega) calcultaion.

#        	bash bse_temperature.sh              		# BSE temperature calculation (in loop) using full daigonalization solver with merged GW + elph database.



#--------------- A4. YPP Post-processing tasks:---------------------------------------------------------#

#        	bash wf_E_sorted.sh                   		# Output data-file of exciton E-sorted weights, etc. at Q =0.

#        	bash wf_exciton.sh                    		# Xcrysden post-processing of 3D-Plotting data-file of exciton wavefunction/amplitude, etc. from BSE diagonalization at Q =0.
       
#        	bash wf_Q_exciton.sh                  		# Finite-Q 3D PDF exciton wavefunction using Xcrysden post-processing from BSE-SlepC diagonalization.

#        	bash wf_Q_exciton_E_sorted.sh         		# Finite-Q output data-file of exciton E-sorted weights, etc.

#        	bash exciton_disp.sh				# Exciton dispersion calculation. Run BSE at all Q points first [bse_slepC.sh].

#        	bash exciton_eliash.sh                		# Exciton-eliashberg function calculation at various temperatures. Run bse_tempeature.sh first to obtain the result.


cd ../YAM_NL_RUN

#----------------A5. NonLinear Optical Properties: CPU intense tasks:-----------------------------------#
# Run first converged GW and BSE for the following tasks: WORKS FOR BOX Z TYPE TRUNCATION.

#      		bash bse_scissor.sh                   		# BSE calculation using INVERSION SOLVER with the automatic G0W0 slope regression evaluation and a converged BSE parameters.

#       	bash nl_setup.sh                      		# Setting up environment for real-Time BSE/non-linear calculations.

#       	bash HXC_collision.sh		       		# Setting of HX [Hartree+SEX] collision matrices.
  	  	
#     		bash real_time_bse.sh                 		# Real-Time BSE calculations.

#     		bash real_time_bse_ypp.sh			# post processing real-time BSE for "delta-like" input impulse.  
		
#     		bash shg.sh					# Second Harmonic Generation for "Quasi-sine-like" input.	
	
#     		bash shg_ypp.sh					# post-processing of Second Harmonic Generation for "Quasi-sine-like" input.

#     		bash shift_current.sh				# Shift-Current Generation for "Soft-Sine" type input.	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	


cd ..
#		rm -rf count


exit;
