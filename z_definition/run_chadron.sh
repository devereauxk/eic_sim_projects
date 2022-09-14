#!/bin/sh

FOLDER=/eic/u/kdevereaux/work/reconstruction/ep_10_100
WORKING_DIR=/eic/u/kdevereaux/work/z_definition

chmod g+rx ${_CONDOR_SCRATCH_DIR}
cd ${_CONDOR_SCRATCH_DIR}

INPUT=$(( 0 + $1 ))
echo $INPUT

DIR=`printf "%05d" $INPUT`
mkdir -p $DIR
pushd $DIR
echo start running in directory $DIR

NEVT=$(( $2 ))
SMEAR_OPT=$3 # 0--no smearing, 1--DM smearing, 2--LBL smearing, 3--Hybrid smearing, 4--ATHENA smearing
BFIELD=$4 # 0--Barbar, 1--Beast
ID_OPT=$5 # 0--no hID (but with eID), 1--PID with no low momentum cutoff, 2--PID with low momentum cutoff & some mis-identified pi, K, 3--PID with low momentum cutoff & all identified pi, K

mkdir $WORKING_DIR/test/S${SMEAR_OPT}_B${BFIELD}_ID${ID_OPT}

ln -s $WORKING_DIR/fast_sim.h .
ln -s $WORKING_DIR/bins_fine.h .
ln -s $WORKING_DIR/D0_tree_patch_zmod.C .
#ln -s /gpfs/mnt/gpfs02/eic/kdevereaux/reconstruction/ATHENA_Resolutions_r.root .
#ln -s /gpfs/mnt/gpfs02/eic/kdevereaux/reconstruction/VertexRes_ATHENA.root .

ln -s $FOLDER/outfiles/ep_{$1}.root
hname=`echo hists-$1.root`

root -b -q 'D0_tree_patch_zmod.C("ep_'$1'.root","'$hname'",'$NEVT','$SMEAR_OPT','$BFIELD','$ID_OPT')'

mv $hname $FOLDER/S${SMEAR_OPT}_B${BFIELD}_ID${ID_OPT}
