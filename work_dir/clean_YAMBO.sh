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

    rm -rf BSE_folder dvscf ELPH graph_data GW_folder HF_folder phonon qpoints work_report*




exit;
