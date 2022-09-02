#!/bin/sh

#FOLDER=BeAGLE_e10Au100_MB # 77 files
FOLDER=pythiaeRHIC_e10Au100_MB_highQ2_small
#FOLDER=pythiaeRHIC_e10p100_MB # 101 high Q2 files
#FOLDER=BeAGLE_e10p100_MB # old BeAGLE files
#FOLDER=BeAGLE102_INCoff_NLO_e10Au100_MB

LIST=`ls -lhtr /gpfs/mnt/gpfs02/eic/wfan/data/$FOLDER/eAu_2.root | awk '{printf("%s\n",$9)}'`
#LIST=`ls -lhtr /gpfs/mnt/gpfs02/eic/wfan/data/$FOLDER/ep_minbias_highQ2*.root | awk '{printf("%s\n",$9)}'`
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

GENTYPE=0 # 0--Pythia6, 1--BeAGLE

ln -s /gpfs/mnt/gpfs02/eic/wfan/HF_vtx_study/EventCounts_EIC.C .

for file in $LIST
do
  if (( $NUM == $1 ))
  then
    fname=`echo $file | awk -F \/ '{printf("%s\n",$9)}'`
    ln -s /gpfs/mnt/gpfs02/eic/wfan/data/$FOLDER/$fname .
    fno=`echo $fname | awk -F \_ '{printf("%s\n",$2)}' | awk -F \. '{printf("%s\n",$1)}'`
    #fno=`echo $fname | awk -F \_ '{printf("%s\n",$4)}' | awk -F \. '{printf("%s\n",$1)}'`
    echo $file
    echo file name is $fname with file number $fno
    hname=`echo hists-$fno.root`
    echo histogram name is $hname

    if [ -a /gpfs/mnt/gpfs02/eic/wfan/HF_vtx_study/hist_output/$hname ]
    then
      echo "File already exists"
    else
      root -b -q 'EventCounts_EIC.C("'$file'","'$hname'",'$NEVT','$GENTYPE')'
      echo "job done...move output file: ${hname}"
      mv $hname /gpfs/mnt/gpfs02/eic/wfan/HF_vtx_study/hist_output/
    fi
  fi
  NUM=$(( $NUM + 1 ))
done

popd

