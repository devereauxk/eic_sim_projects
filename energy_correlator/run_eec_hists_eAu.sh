#!/bin/bash
IN_DIR=/eic/u/kdevereaux/work/eHIJING/eHIJING-pythia/eHIJING-examples/Events/eAu_1E8_K0_condor
OUT_DIR=/eic/u/kdevereaux/work/energy_correlator/eHIJING/eAu_1E8_K0_pow2

if [ -z "$1" ]
then
        echo "No job number set."
        echo "Please run as ./run_ep.sh jobnumber"
        echo "Exiting..."
        exit 1
fi

#Go into scratch directory
chmod g+rx ${_CONDOR_SCRATCH_DIR}
cd ${_CONDOR_SCRATCH_DIR}

#Make subdirectory and move there
INPUT=$(( 0 + $1 ))
echo $INPUT
DIR=`printf "%04d" $INPUT`
mkdir $DIR

ln -s $IN_DIR/eAu_${INPUT}.dat
ln -s /eic/u/kdevereaux/work/energy_correlator/eec_hists.C

mkdir $OUT_DIR

root -l -b -q "eec_hists.C(\"eAu_${INPUT}.dat\",\"$OUT_DIR/hists_eec_${INPUT}.root\", 1, 2147.95, 100, 1, 2)"

# 10 on 100: 2147.95 GeV
# 18 on 110: 4252.93 GeV
# 5 on 41: 440.378 GeV
