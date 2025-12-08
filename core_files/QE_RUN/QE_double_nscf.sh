#!/bin/bash





cd $DIR
                         		echo "============================================================="
				    	echo ">>> Running QE NSCF Computation at Double k-grid."
				    	


                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] QE NSCF-double-grid report:"
                           } >> "$DIR/work_report.txt"

# Create necessary folders
mkdir -p database database_double nscf-dg

# Define paths
PATH2="$DIR/scf_soc"
PATH3="$DIR/nscf"
PATH8="$DIR/nscf-dg"

# ------------------ NSCF Double-Grid Calculation ------------------ #
cd "$PATH2"
cp -r * "$PATH8"

cd "$PATH8"
rm -f scf.out

# Rename input file
mv "$prefix.scf.in" "$prefix.nscf_dg.in"

# Replace 'scf' with 'nscf'
sed -i.bak "s/calculation *= *'scf'/calculation = 'nscf'/g" "$prefix.nscf_dg.in"

# Truncate everything after (and including) K_POINTS
awk '/^[ ]*K_POINTS/ { print; exit } { print }' "$prefix.nscf_dg.in" > tmp1

# Append your dgkpoint
echo "$dgkpoint" >> tmp1

# Insert NBND after 'ntyp'
awk -v var="$NBND" '/ntyp/{print;print var;next}1' tmp1 > "$prefix.nscf_dg.in"

# Clean
rm -f tmp1 *.bak


# Run NSCF calculation
$MPIRUN_nscf "$prefix.nscf_dg.in" > "$prefix.nscf_dg.out"


# ------------------ Making Double Grid ------------------ #
cd "$PATH8/$prefix.save"
$P2Y 2> out
scp -r SAVE "$DIR/database_double/."
$YAMBO 2> out
rm -rf SAVE l_p* out r_setup*


cd "$PATH3/$prefix.save"
$P2Y 2> out
scp -r SAVE "$DIR/database/."
rm -rf SAVE l_p* out


# Navigate to the database directory
cd "$DIR/database"

# Run YAMBO with error output redirection (for logging)
$YAMBO 2> out

# Run YPP for Double Grid (with necessary modifications)
$YPP -m -nompi -F ypp_dg.in 2> out

# Modify ypp_dg.in file using sed (both macOS and Linux/Unix compatible)
sed 's/#SkipCheck/SkipCheck/g' ypp_dg.in > tmp && mv tmp ypp_dg.in
sed '22s/"none" |/"..\/database_double" |/g' ypp_dg.in > tmp && mv tmp ypp_dg.in

# Run YPP again with no MPI for non-parallel execution
$YPP -nompi -F ypp_dg.in 2> output_dbs

                         		
				    	echo ">>> Ending QE NSCF Computation at Double k-grid."
				    	echo "============================================================="


                            {
                            echo "    → NSCF on double-grid done."
                            echo "    → Grid Used               : = $dgkpoint"
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of NSCF report."
                            echo ""
                           } >> "$DIR/work_report.txt"


exit
