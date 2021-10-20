#!/bin/sh

FOLDER=/gpfs/mnt/gpfs02/eic/kdevereaux/reconstruction/eAu_18_110/S0_B0_ID0_DCA0
OUT_FOLDER=/gpfs/mnt/gpfs02/eic/kdevereaux/reconstruction/hist_output/eAu_5E7_0_0_0_0

LIST=`ls -lhtr $FOLDER/*.root | awk '{printf("%s\n",$9)}'` # since file full dir given as 9th output paran for ls -lhtr
NUM=0

chmod g+rx ${_CONDOR_SCRATCH_DIR}
cd ${_CONDOR_SCRATCH_DIR}

INPUT=$(( 0 + $1 ))
echo $INPUT

DIR=`printf "%05d" $INPUT`
mkdir -p $DIR
pushd $DIR
echo start running in directory $DIR

ls -s /gpfs/mnt/gpfs02/eic/kdevereaux/reconstruction/plot_histogram.C .

for file in $LIST
do
  if (( $NUM == $1 ))
  then

    root -b -q 'plot_histogram.C("'$file'","'$OUT_FOLDER'","e + Au @ 18, 110 GeV", "50000000")'

  fi
  NUM=$(( $NUM + 1 ))
done

popd
