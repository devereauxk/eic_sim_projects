#!/bin/bash
IN_DIR=/eic/u/kdevereaux/work/eHIJING/eHIJING-pythia/eHIJING-examples/Events/ep_10_100_K0_Q2x
OUT_DIR=/eic/u/kdevereaux/work/energy_correlator/eHIJING/ep_10_100_K0_noboost

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

mkdir $OUT_DIR

#                                                                                                        DOUBLE CHECK THESE
#                                                                                                        power, boost, Q2x
root -l -b -q "eec_hists.C(\"ep_${INPUT}.dat\",\"$OUT_DIR/hists_eec_${INPUT}.root\", 1, 2131.56, 100, 0, 0.5, 0, 1)"

#root -l -b -q "eec_hists.C(\"/eic/u/kdevereaux/work/eHIJING/eHIJING-pythia/eHIJING-examples/Events/ep_10_100_K0/ep_1.dat\",\"temp.root\", 1, 2131.56, 100, 0, 0.5, 1, 0)"
