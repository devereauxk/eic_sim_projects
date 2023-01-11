#!/bin/bash
IN_DIR=/eic/u/kdevereaux/work/reconstruction/ep_10_100/outfiles
OUT_DIR=/eic/u/kdevereaux/work/energy_correlator/eHIJING/ep_10_100_Q2x_pythia

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

ln -s $IN_DIR/ep_${INPUT}.root
ln -s /eic/u/kdevereaux/work/energy_correlator/Q2_x.C

mkdir $OUT_DIR

root -l -b -q "Q2_x.C(\"ep_${INPUT}.root\",\"$OUT_DIR/hists_Q2_x_${INPUT}.root\", 0)"

#root -l -q "Q2_x.C(\"/eic/u/kdevereaux/work/reconstruction/ep_10_100/outfiles/ep_1.root\",\"./test.root\", 0)"
