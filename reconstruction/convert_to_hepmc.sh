#!/bin/sh

FOLDER=/gpfs/mnt/gpfs02/eic/wfan/event_gen/pythiaeRHIC/10_100/
FOLDER_OUT=/gpfs/mnt/gpfs02/eic/kdevereaux/reconstruction/hepmc_out

LIST=`ls -lhtr $FOLDER/ep_filtered_highQ2_*.root | awk '{printf("%s\n",$9)}'`
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

#mkdir $FOLDER/outHepMC

for file in $LIST
do
  if (( $NUM == $1 ))
  then
    echo $file
    fname=`echo $file | awk -F \/ '{printf("%s\n",$10)}'` # $10 since files are ten levels down from root directory
    ln -s $FOLDER/$fname .
    fno=`echo $fname | awk -F \_ '{printf("%s\n",$2)}' | awk -F \. '{printf("%s\n",$1)}'`    # $2 since file of form ep_<some number>.root
    species=`echo $fname | awk -F \. '{printf("%s\n",$1)}'`
    echo $species
    echo file name is $fname with file number $fno
    hname=`echo $species.hepmc`
    echo hepmc file name is $hname

    if [ -a $FOLDER_OUT/$hname ]
    then
      echo "File already exists"
    else
      echo 'TreeToHepMC("'$fname'")' | eic-smear
      echo "job done...move output file: ${hname}"
      mv $hname $FOLDER_OUT/
    fi
  fi
  NUM=$(( $NUM + 1 ))
done

popd
