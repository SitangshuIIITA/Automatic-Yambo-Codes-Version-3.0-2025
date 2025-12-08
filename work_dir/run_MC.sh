#!/bin/bash






# Note: What to run: Make sure you installed them. 

#	1. You need two EPW run calculation files: epw.out that stores gkkp elements.
#	2. python3 with mpi & gif generators library. 

#	Suggestion: Before running monte carlo engines, run a converged QE relax, scf, nscf, bands and orbital projections and EPW matrix elements and mobility out files first.


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


cd ../core_files/MC_RUN 

# Run sequentially or at one shot, the following lines :


		#bash electron_monte_carlo.sh		# This extracts Einstien electronic mobility and diffucion coefficient using zero field algorithm of MonteCarlo.
		
		#bash hole_monte_carlo.sh		# This extracts Einstien hole mobility and diffucion coefficient using zero field algorithm of MonteCarlo.


		bash electron_drift.sh			# This extracts electronic drift velocity vs electric field [along x direction] algorithm of MonteCarlo.



######################################################  GAME OVER ###############################################################

	
rm -rf count

exit;
