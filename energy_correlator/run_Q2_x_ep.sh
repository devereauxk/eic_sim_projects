#!/bin/bash
IN_DIR=/eic/u/kdevereaux/work/eHIJING/eHIJING-pythia/eHIJING-examples/Events/ep_10_100_K0_Q2x
OUT_DIR=/eic/u/kdevereaux/work/energy_correlator/eHIJING/ep_10_100_K0_Q2x

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

ln -s $IN_DIR/ep_${INPUT}.dat
ln -s /eic/u/kdevereaux/work/energy_correlator/Q2_x.C

mkdir $OUT_DIR

root -l -b -q "Q2_x.C(\"ep_${INPUT}.dat\",\"$OUT_DIR/hists_Q2_x_${INPUT}.root\")"

#root -l -q "Q2_x.C(\"$EVENTS/ep_10_100_Q2x/ep_1.dat\",\"eHIJING/ep_10_100_K0_Q2x/test_1.root\")"
