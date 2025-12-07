#!/bin/bash





cat > variables.txt << EOF
# Atoms array
atoms=($atoms)

# Mass array
mass_atom=($mass_atom)

# Pseudo-prefix array
pseudo_prefix=($pseudo_prefix)
EOF

source ./variables.txt

# Create the atomic_species.txt file
cat > atomic_species.txt << EOF
ATOMIC_SPECIES
EOF

# Loop through the atoms array and generate the required format
for i in "${!atoms[@]}"; do
  
  echo "  ${atoms[$i]} ${mass_atom[$i]} ${pseudo_prefix[$i]}.upf" >> atomic_species.txt
    done


#------------------------------1. DFT QE input-----------------------------------#
#--------------------------------------------------------------------------------#
cat > $prefix.relax << EOF

&CONTROL
                                        wf_collect = .true.
                                       calculation = 'vc-relax'
                                         verbosity = 'high'
                                        pseudo_dir = '$PSEUDO_DIR'
                                     forc_conv_thr = $forc_conv_thr
                                            prefix = '$prefix'
                                     etot_conv_thr = $etot_conv_thr
/&end

&SYSTEM
                                           ecutwfc = $ecutwfc
                                       occupations = '$occupation'
                                         celldm(1) = $celldm_1
                                             ibrav = $ibrav
                                         celldm(3) = $celldm_3
                                               nat = $nat
                                              ntyp = $ntyp
                                   assume_isolated = '$assume_isolated'

/&end
/

&ELECTRONS
                                   diagonalization = 'david'
                                          conv_thr = 1e-10
                                  electron_maxstep = 400
/&end
/ 
 
&ions
                                      ion_dynamics = 'bfgs'
/&end
&cell
                                       cell_dofree = '$cell_dofree'
                                     cell_dynamics = 'bfgs'
/&end

EOF

cat  atomic_species.txt >> $prefix.relax
source ../core_files/QE_VAR/QE_variables.sh
echo "ATOMIC_POSITIONS {crystal}" >> $prefix.relax
for pos in "${positions[@]}"; do
    echo "  $pos" >> $prefix.relax
done
echo "K_POINTS {automatic}" >> $prefix.relax
echo "$kgrid_nonsoc $shifted_grid" >> $prefix.relax
#---------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------#


#------2. Choice of BZ route for band diagram: Change as per requirement----------------------#

#------------------------Do not change the format below---------------------------------------#
#----Change the BZ route name order (example: G M K G) in the QE_variables.sh file------------#
cat > $prefix.band_route << EOF
0.00000 | 0.00000 | 0.00000 |
0.50000 | 0.00000 | 0.00000 |
0.33333 | 0.33333 | 0.00000 |
0.00000 | 0.00000 | 0.00000 |
EOF
#---------------------------------------------------------------------------------------------#


