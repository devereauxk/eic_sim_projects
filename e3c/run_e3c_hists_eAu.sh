#!/bin/bash
IN_DIR=/eic/u/kdevereaux/work/eHIJING/eHIJING-pythia/eHIJING-examples/Events/eAu_1E8_K4
OUT_DIR=/eic/u/kdevereaux/work/e3c/analysis/eAu_1E8_K4_pow05

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
ln -s /eic/u/kdevereaux/work/e3c/e3c_hists.C

mkdir $OUT_DIR

#                                                                                                        DOUBLE CHECK THESE
#                                                                                                        power, boost, Q2x
root -l -b -q "e3c_hists.C(\"eAu_${INPUT}.dat\",\"$OUT_DIR/hists_eec_${INPUT}.root\", 1, 2130.17, 100, 1, 0.5, 1, 1)"