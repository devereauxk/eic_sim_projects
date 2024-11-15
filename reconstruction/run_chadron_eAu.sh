#!/bin/sh

FOLDER=/gpfs/mnt/gpfs02/eic/kdevereaux/reconstruction/eAu_10_110

LIST=`ls -lhtr $FOLDER/outForPythiaMode/*.root | awk '{printf("%s\n",$9)}'` # since file full dir given as 9th output paran for ls -lhtr
NUM=0

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

mkdir $FOLDER/S${SMEAR_OPT}_B${BFIELD}_ID${ID_OPT}

ln -s /gpfs/mnt/gpfs02/eic/kdevereaux/reconstruction/fast_sim.h .
ln -s /gpfs/mnt/gpfs02/eic/kdevereaux/reconstruction/bins.h .
ln -s /gpfs/mnt/gpfs02/eic/kdevereaux/reconstruction/D0_tree_patch.C .
ln -s /gpfs/mnt/gpfs02/eic/kdevereaux/reconstruction/ATHENA_Resolutions_r.root .
ln -s /gpfs/mnt/gpfs02/eic/kdevereaux/reconstruction/VertexRes_ATHENA.root .

for file in $LIST
do
  if (( $NUM == $1 ))
  then
    fname=`echo $file | awk -F \/ '{printf("%s\n",$10)}'` # $10 since files are ten levels down from root directory
    ln -s $FOLDER/outForPythiaMode/$fname .
    fno=`echo $fname | awk -F \_ '{printf("%s\n",$2)}' | awk -F \. '{printf("%s\n",$1)}'`    # $2 since file of form ep_<some number>.root
    echo $file
    echo file name is $fname with file number $fno
    hname=`echo hists-$fno.root`
    echo histogram name is $hname

    if [ -a $FOLDER/S${SMEAR_OPT}_B${BFIELD}_ID${ID_OPT}/$hname ]
    then
      echo "File already exists"
    else
      root -b -q 'D0_tree_patch.C("'$file'","'$hname'",'$NEVT','$SMEAR_OPT','$BFIELD','$ID_OPT')'
      echo "job done...move output file: ${hname}"
      mv $hname $FOLDER/S${SMEAR_OPT}_B${BFIELD}_ID${ID_OPT}
    fi
  fi
  NUM=$(( $NUM + 1 ))
done

popd
