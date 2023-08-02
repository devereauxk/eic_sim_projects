#!/bin/bash
IN_DIR=/eic/u/kdevereaux/work/eHIJING/eHIJING-pythia/eHIJING-examples/Events/eAu_10_100_K4_density
OUT_DIR=/eic/u/kdevereaux/work/e3c/analysis/eAu_10_100_K4_density_pow1_byhand

# Go into scratch directory (if needed, you can specify a different temporary directory)
# For this example, we'll use the script's current working directory.
# Uncomment the line below if you need to use a specific scratch directory.
# chmod g+rx ${_CONDOR_SCRATCH_DIR}
cd "$(dirname "${BASH_SOURCE[0]}")"

# Create the output directory if it doesn't exist
mkdir -p "$OUT_DIR"

for ((INPUT=0; INPUT<=1000; INPUT++))
do
    echo "processing file: $INPUT"
    
    DIR=$(printf "%04d" "$INPUT")
    mkdir -p "$DIR"

    # Symlink the input file and the processing script to the temporary directory
    ln -s "$IN_DIR/eAu_${INPUT}.dat" "$DIR/eAu_${INPUT}.dat"
    ln -s /eic/u/kdevereaux/work/e3c/e3c_hists.C "$DIR/e3c_hists.C"

    # Process the current input file
    #                                                                                     DOUBLE CHECK THESE
    #                                                                                     power, boost, Q2x
    root -l -b -q "$DIR/e3c_hists.C(\"$DIR/eAu_${INPUT}.dat\",\"$OUT_DIR/hists_eec_${INPUT}.root\", 1, 2147.95, 100, 1, 1.0, 1, 1)"

    # Optionally, remove the temporary directory if you don't need it for any other purposes.
    # rm -rf "$DIR"
done
