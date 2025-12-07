#!/bin/bash





# This will clean all previous folders. BEWARE!

cd ../core_files/QE_VAR
set -a
source ./QE_variables.sh
set +a


cd ../YAM_VAR
set -a
source ./yambo_variables.sh
set +a


cd $DIR

    rm -rf atomic* band* $prefix.band_route databse* $prefix.relax electronic* nscf* database* phonon_disp qe_band* relax scf* variables* work_report* Thermo_PW orbital_* qe_calc* variables.txt work_report.txt atomic species.txt bn.band_route




exit;
