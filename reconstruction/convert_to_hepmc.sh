#!/bin/sh

FOLDER=/gpfs/mnt/gpfs02/eic/kdevereaux/reconstruction/eAu_10_110

LIST=`ls -lhtr $FOLDER/outForPythiaMode/*.root | awk '{printf("%s\n",$9)}'`
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

mkdir $FOLDER/outHepMC

for file in $LIST
do
  if (( $NUM == $1 ))
  then
    echo $file
    fname=`echo $file | awk -F \/ '{printf("%s\n",$10)}'` # $10 since files are ten levels down from root directory
    ln -s $FOLDER/outForPythiaMode/$fname .
    fno=`echo $fname | awk -F \_ '{printf("%s\n",$2)}' | awk -F \. '{printf("%s\n",$1)}'`    # $2 since file of form ep_<some number>.root
    species=`echo $fname | awk -F \. 'printf("%s/n",$2)'`
    echo $species
    echo file name is $fname with file number $fno
    hname=`echo $species.hepmc`
    echo hepmc file name is $hname

    if [ -a $FOLDER/$fname ]
    then
      echo 'TreeToHepMC("'$FOLDER/$fname'.root")' | eic-smear
      mv $hname $FOLDER/outHepMC
    fi
  fi
  NUM=$(( $NUM + 1 ))
done

popd
