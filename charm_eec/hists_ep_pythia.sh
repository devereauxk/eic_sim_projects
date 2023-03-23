#!/bin/bash
#IN_DIR=/gpfs/mnt/gpfs02/eic/wfan/data/pythia8HepMC_e10p100_ft_MB
IN_DIR=/gpfs/mnt/gpfs02/eic/wfan/data/pythia8HepMC_e10p100_MB
OUT_DIR=/eic/u/kdevereaux/work/charm_eec/analysis/ep_10_100_pythia8_D0injet

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

ln -s $IN_DIR/ep_minbias_highQ2_${INPUT}.root
ln -s /eic/u/kdevereaux/work/charm_eec/charm_eec_hists.C

mkdir $OUT_DIR

POW=$2
DO_BOOST=$3
DO_Q2X=$4
DO_FORCE_INJET=$5
DO_FORCE_INPAIR=$6
FIXED_ID=$7

#                                                                                                        DOUBLE CHECK THESE
#                                                                                                                            pow,boost,Q2x,jet,pair
root -l -b -q "charm_eec_hists.C(\"ep_minbias_highQ2_${INPUT}.root\",\"$OUT_DIR/hists_eec_${INPUT}.root\", -1, 2131, 100, 0, $POW, $DO_BOOST, $DO_Q2X, $DO_FORCE_INJET, $DO_FORCE_INPAIR, $FIXED_ID)"
#root -l -b -q "eec_hists.C(\"hepmcout-${INPUT}.root\",\"$OUT_DIR/hists_eec_${INPUT}.root\", -1, 2131, 100, 0, 0.5, 1, 1)"

#root -l -b -q "charm_eec_hists.C(\"/gpfs/mnt/gpfs02/eic/wfan/data/pythia8HepMC_e10p100_MB/ep_minbias_highQ2_99.root\",\"temp.root\", -1, 2131, 100, 0, 0.5, 0, 1, 1, 1, 421)"
