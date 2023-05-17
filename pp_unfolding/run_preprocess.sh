#!/bin/bash

IN_DIR=$2
OUT_DIR=$3
POW=$4

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

ln -s $IN_DIR/hepmcout_${INPUT}.root
ln -s /eic/u/kdevereaux/work/pp_unfolding/preprocess.C

mkdir $OUT_DIR

#                                                                                         DOUBLE CHECK THESE
root -l -b -q "preprocess.C(\"hepmcout_${INPUT}.root\",\"$OUT_DIR/eec_${INPUT}.root\", $POW)"
