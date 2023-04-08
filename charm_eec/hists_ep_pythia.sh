#!/bin/bash
#IN_DIR=/gpfs/mnt/gpfs02/eic/wfan/data/pythia8HepMC_e10p100_ft_MB
#IN_DIR=/gpfs/mnt/gpfs02/eic/wfan/data/pythia8HepMC_e10p100_MB
#IN_DIR=/eic/u/kdevereaux/work/pythia8/ep_18_275/outfiles

POW=$2
DO_BOOST=$3
DO_FORCE_INJET=$4
DO_FORCE_INPAIR=$5
FIXED_ID=$6
IN_DIR=$7
OUT_DIR=$8

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
ln -s /eic/u/kdevereaux/work/charm_eec/charm_eec_hists.C

mkdir $OUT_DIR

#                                                                                                        DOUBLE CHECK THESE
#                                                                                                                            pow,boost,Q2x,jet,pair
root -l -b -q "charm_eec_hists.C(\"hepmcout_${INPUT}.root\",\"$OUT_DIR/hists_eec_${INPUT}.root\", -1, 2131, 100, 0, $POW, $DO_BOOST, 1, $DO_FORCE_INJET, $DO_FORCE_INPAIR, $FIXED_ID)"

#root -l -b -q "charm_eec_hists.C(\"/eic/u/kdevereaux/work/pythia8/ep_18_275_nodecay/outfiles/hepmcout_1.root\",\"temp.root\", -1, 2131, 100, 0, 1, 0, 1, 1, 0, 421)"
