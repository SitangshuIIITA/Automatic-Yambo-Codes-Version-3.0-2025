#!/bin/bash





cd $DIR


                            {   echo ""
                                echo "[$(date '+%Y-%m-%d %H:%M:%S')] QE bands projection & DOS report:"
                           } >> "$DIR/work_report.txt"


                         		echo "============================================================="
				    	echo ">>> Running QE Spin Resolved DOS Computation."
				    	


PATH4=$DIR/bands

# ------------------- QE bands: bandx and projwfc calculations -------------------

cd $PATH4

# Logging
echo "[$(date)] Starting QE Band Structure Post-Processing for $prefix" >> bands.log

# ---------------------- Generate bandx input file -----------------------
cat > "$prefix.bandx.in" << EOF
&BANDS
    prefix = "$prefix"
    outdir = "."
    filband = "Bandx.dat"
/
EOF

# Run bandx calculation
echo "[$(date)] Running bandx calculation..." >> bands.log
$MPIRUN_Bndx "$prefix.bandx.in" > bandx.out
echo "[$(date)] Finished bandx calculation." >> bands.log

# ---------------------- Determine kresolveddos from lspinorb -----------------------
kresolveddos_value=".true."  # default

if grep -qi "lspinorb *= *.true." "$prefix.bands.in"; then
    kresolveddos_value=".true."
    echo "[$(date)] Detected lspinorb = .true. → Enabling kresolveddos." >> bands.log
else
    echo "[$(date)] Detected lspinorb = .false. or not set → Disabling kresolveddos." >> bands.log
fi

# ---------------------- Generate projwfc input file -----------------------
cat > "$prefix.projwfc.in" << EOF
&projwfc
    prefix = '$prefix'
    outdir = './'
    ngauss = 0,
    degauss = 0.016748,
    DeltaE = 0.030000,
    kresolveddos = $kresolveddos_value
/
EOF

# Run projwfc calculation
echo "[$(date)] Running projwfc calculation with kresolveddos = $kresolveddos_value..." >> bands.log
$MPIRUN_proj "$prefix.projwfc.in" > projwfc.out
echo "[$(date)] Finished projwfc calculation." >> bands.log

mkdir $DIR/orbital_plot_data
mkdir $DIR/orbital_plot_data/orbital_data

scp -r band* Bandx* $prefix.pdos* $prefix.pdos_tot nscf.out projwfc.out $DIR/orbital_plot_data/orbital_data
scp -r ../../core_files/QE_RUN/file_analyzer.py $DIR/orbital_plot_data
scp -r ../../core_files/QE_RUN/orbital_plot.py $DIR/orbital_plot_data

cd $DIR/orbital_plot_data
python3 orbital_plot.py

                         		
				    	echo ">>> Ending QE Spin Resolved DOS Computation."
				    	echo "============================================================="

                            {
                            echo "    → Projection and DOS done. Check orbital_plot_data folder to plot BandStruc and orbitals."
                            echo "    → Grid Used               : = $kgrid"
                            echo "[$(date '+%Y-%m-%d %H:%M:%S')] End of bands projection & DOS report."
                            echo ""
                           } >> "$DIR/work_report.txt"


exit
