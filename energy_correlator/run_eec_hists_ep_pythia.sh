#!/bin/bash
IN_DIR=/gpfs/mnt/gpfs02/eic/wfan/data/pythia8HepMC_e10p100_ft_MB
OUT_DIR=/eic/u/kdevereaux/work/energy_correlator/eHIJING/ep_10_100_pythia8_ft

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

ln -s $IN_DIR/hepmcout-${INPUT}.root
ln -s /eic/u/kdevereaux/work/energy_correlator/eec_hists.C

mkdir $OUT_DIR

root -l -b -q "eec_hists.C(\"hepmcout-${INPUT}.root\",\"$OUT_DIR/hists_eec_${INPUT}.root\", -1, 2131, 100, 0, 0.5, 1)"

#root -l -b -q "eec_hists.C(\"/gpfs/mnt/gpfs02/eic/wfan/data/pythia8HepMC_e10p100_ft_MB/hepmcout-1.root\",\"temp.root\", -1, 2131, 100, 0, 0.5, 1)"
