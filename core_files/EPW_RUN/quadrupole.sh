#!/bin/bash





source ../../core_files/EPW_VAR/EPW_variables.sh
source ../../core_files/QE_VAR/QE_variables.sh

cd $DIR

					echo "=========================================================================="
				    	echo ">>> Running & Optimizing DFPT q-points to Fit For File Generation: quadrupole.fmt"
			

                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] Quadrupole report:"
                                } >> "$DIR/EPW_work_report.txt"


mkdir -p epw_$prefix/Quadrupole


PATH1="$DIR/epw_$prefix/w90_calc"
PATH2="$PATH1/relax_w90"
PATH3="$DIR/epw_$prefix/Quadrupole"

cd "$PATH3" || exit 1


###################	Relax computation 	############################################

scp -r $PATH2/$prefix.relax ./. 
sed -i "/occupations *= *'fixed'/d" "$prefix.relax"

awk -v var="					   outdir  = '"work"' " '/calculation/{print;print var;next}1' "$prefix.relax" > tmp && mv tmp "$prefix.relax"

awk -v var="  				       occupations = '"smearing"' " '/ecutwfc/{print;print var;next}1' "$prefix.relax" > tmp && mv tmp "$prefix.relax"
awk -v var="  					  smearing = '"marzari-vanderbilt"'" '/ecutwfc/{print;print var;next}1' "$prefix.relax" > tmp && mv tmp "$prefix.relax"
awk -v var="                                           degauss = 0.02" '/ecutwfc/{print;print var;next}1' "$prefix.relax" > tmp && mv tmp "$prefix.relax"



		$QUAD_MPIRUN_epw_PW "$prefix.relax" > relax.out







###################	scf computation 	############################################

# Extract final atomic coordinates
awk '/Begin final coordinates/,/End final coordinates/' relax.out \
    | awk 'NR > 3 { print }' \
    | awk 'NR > 2 { print last } { last = $0 }' > 1

# Append last 3 lines of original input
tail -2 "$prefix.relax" >> 1

# Extract everything from &CONTROL to just before ATOMIC_POSITIONS
awk '/&CONTROL/,/ATOMIC_POSITIONS/' "$prefix.relax" \
    | awk 'NR > 1 { print last } { last = $0 }' >> 2

# Append new coordinates
cat 1 >> 2

# Replace vc-relax with scf (portable)
sed 's/vc-relax/scf/g' 2 > tmp && mv tmp 2

# Replace ibrav with 0 (portable)
sed -E "s/ibrav[[:space:]]*=[[:space:]]*$ibrav/ibrav = 0/g" 2 > tmp && mv tmp 2

# Inject lspinorb, and noncolin options

awk '/ntyp/{print;print "'"                                          lspinorb = $lspinorb"'";next}1' 2 > 3

awk '/ntyp/{print;print "'"                                          noncolin = $noncolin"'";next}1' 3 > 4


mv 4 2
mv 2 "$prefix.scf.in"

# Clean and finalize
rm -f 1 2 3 4


		$QUAD_MPIRUN_epw_PW "$prefix.scf.in" > scf.out





###################	bare phonon calculation  ##################################
# Create phonon input file
cat > "ph.in" << EOF
&inputph
    verbosity = 'high',
       outdir = 'work'
       tr2_ph = $TRN_PH,
    alpha_mix = $alpha_mx,
       prefix = '$prefix',
       fildyn = 'dyn',
     fildvscf = 'dvscf',
        ldisp = .true.,
         bare = .true.
    nq1 = $NQ1, nq2 = $NQ2, nq3 = $NQ3,
/
EOF



		$QUAD_EPW_PH_DISP_RUN -inp "ph.in" > ph.out

		ph2epw





###################	nscf computation at uniform grid  ##################################

perl $PERL_PATH $kgrid_nonsoc > QE_k_list.in

# Rename SCF input to NSCF input
scp -r "$prefix.scf.in" "$prefix.nscf.in"

# Replace calculation = 'scf' with 'nscf'
awk '{ if ($0 ~ /calculation *= *'\''scf'\''/) gsub(/scf/, "nscf"); print }' "$prefix.nscf.in" > tmp && mv tmp "$prefix.nscf.in"

# Delete line containing "force_symmorphic = .TRUE."
grep -v "force_symmorphic = .TRUE." "$prefix.nscf.in" > tmp && mv tmp "$prefix.nscf.in"


# Insert number of bands (NBND) after 'ntyp' line
awk -v var="$NBND" '/ntyp/{print;print var;next}1' "$prefix.nscf.in" > tmp && mv tmp "$prefix.nscf.in"

# Truncate file after K_POINTS and append new k-grid
awk '/K_POINTS {automatic}/ {exit} {print}' "$prefix.nscf.in" > tmp && mv tmp "$prefix.nscf.in"

echo "$(cat QE_k_list.in)" >> "$prefix.nscf.in"


		$QUAD_MPIRUN_epw_nscf "$prefix.nscf.in" > nscf.out


read NK1 NK2 NK3 <<< "$kgrid"

cat > epw.in << EOF
&inputepw
		  	   prefix = '$prefix',
		  	   outdir = 'work'
		  	dvscf_dir = 'save'
		  	     elph = .true.
		  	 epbwrite = .true.
		  	  epbread = .false.
		  	 epwwrite = .true.
		  	  epwread = .false.
		  	
		  	     nkf1  = $nkf1
		  	     nkf2  = $nkf2
		  	     nkf3  = $nkf3
		  	
		  	     nqf1  = $nqf1
		  	     nqf2  = $nqf2
		  	     nqf3  = $nqf3
		  	
		  	     nk1  = $NK1
		  	     nk2  = $NK2
		  	     nk3  = $NK3
		  	
		  	     nq1  = $NQ1
		  	     nq2  = $NQ2
		  	     nq3  = $NQ3

		       wannierize = $wannierize
		          nbndsub = $num_wann
		         num_iter = $num_iter
		     dis_froz_min = $dis_froz_min
		     dis_froz_max = $dis_froz_max
			   iprint = $iprint 
		    bands_skipped = 'exclude_bands = $exclude_bands'


EOF

# Write orbital projections
for i in "${!orbitals[@]}"; do
    idx=$((i+1))
    echo "proj($idx) = '${orbitals[$i]}'" >> epw.in
done


# --- wdata block ---

# Loop over path segments
for ((i=0; i<${#labels[@]}-1; i++)); do
    echo "wdata($((i+3))) = '${labels[$i]} ${coords[$i]} ${labels[$((i+1))]} ${coords[$((i+1))]}'" >> epw.in
done

last_idx=$(( ${#labels[@]} + 1 ))
echo "wdata($((last_idx))) = 'dis_win_min = $dis_win_min'" >> epw.in
echo "wdata($((last_idx+1))) = 'dis_win_max = $dis_win_max'" >> epw.in
echo "wdata($((last_idx+2))) = 'dis_num_iter = $dis_num_iter'" >> epw.in
echo "wdata($((last_idx+3))) = 'conv_tol = $conv_tol'" >> epw.in

echo "/" >> epw.in


		$QUAD_MPIRUN_EPW -input epw.in > epw.out 



		$QUAD_MPIRUN_epw_PW "$prefix.scf.in" > scf.out




###################	bare DFPT phonon calculation on uniform grid  ##################################
# Create phonon input file
cat > "phref.in" << EOF
&inputph
        verbosity = 'high',
           outdir = 'work'
           tr2_ph = $TRN_PH,
        alpha_mix = $alpha_mx,
           prefix = '$prefix',
           fildyn = 'dynref',
         fildvscf = 'dvscfref',
            ldisp = .true.,
             bare = .true.
            qplot = .true.
  electron_phonon = 'simple'
            ibndprt = $ibndprt
/
EOF

scp -r ../../../core_files/EPW_RUN/qpath_M_Gamma_K.py ./.

python3 qpath_M_Gamma_K.py
cat qpath_GM.txt >> phref.in  

		$QUAD_EPW_PH_DISP_RUN -inp "phref.in" > phref.out
		
		

scp -r ../../../core_files/EPW_RUN/Quad_optimize.py ./.

		python3 Quad_optimize.py > log.txt 2>&1

scp -r ../../../core_files/EPW_RUN/g_plot.py ./.
		
		python3 g_plot.py 
scp -r fitQa.png "$DIR/epw_$prefix/EPW_plots/"
scp -r fitQb.png "$DIR/epw_$prefix/EPW_plots/"
scp -r fitQ.png "$DIR/epw_$prefix/EPW_plots/"

				    	echo ">>> Ending Plots For Optimized DFPT q-points to Fit: quadrupole.fmt. Check EPW_plots Folder. "
					echo "=========================================================================="
	

          		 {  echo "    → Quadrupole.fmt and validation images are generated."
                            echo "    → Plottable folder is inside files $DIR/epw_{"$prefix"}/EPW_plots"
                            echo "    → Before final production, cross check the convergences."
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of fitQ report."
			    echo ""
                           } >> "$DIR/EPW_work_report.txt"
		

exit 0
