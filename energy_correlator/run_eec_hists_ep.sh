#!/bin/bash
IN_DIR=/eic/u/kdevereaux/work/eHIJING/eHIJING-pythia/eHIJING-examples/Events/ep_1E6_condor
OUT_DIR=/eic/u/kdevereaux/work/energy_correlator/eHIJING/ep_1E6_condor

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
ln -s /eic/u/kdevereaux/work/energy_correlator/eec_hists.C

root -l -b -q "eec_hists.C(\"ep_${INPUT}.dat\",\"$OUT_DIR/hists_eec_${INPUT}.root\", 1, 2131.56, 100, 0, 0)"
