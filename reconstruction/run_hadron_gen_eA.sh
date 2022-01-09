#!/bin/sh

#FOLDER=BeAGLE_e10Au110_MB # new BeAGLE files with the vertex fix (496 root output)
#FOLDER=BeAGLE_e5Au41_MB # new BeAGLE files with the vertex fix (499 root output)
#FOLDER=BeAGLE_e10p100_MB # old BeAGLE files
FOLDER=BeAGLE_e18C110_MB
#FOLDER=pythiaeRHIC_e18p110_MB

#LIST=`ls -lhtr /gpfs/mnt/gpfs02/eic/wfan/data/$FOLDER/eAu_*.root | awk '{printf("%s\n",$9)}'`
LIST=`ls -lhtr /gpfs/mnt/gpfs02/eic/wfan/data/$FOLDER/merged.root | awk '{printf("%s\n",$9)}'`
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

mkdir /gpfs/mnt/gpfs02/eic/wfan/HF_vtx_study/hist_hadron_gen_eC

ln -s /gpfs/mnt/gpfs02/eic/wfan/HF_vtx_study/bins.h .
ln -s /gpfs/mnt/gpfs02/eic/wfan/HF_vtx_study/ana_hadron_gen.C .

for file in $LIST
do
  if (( $NUM == $1 ))
  then
    fname=`echo $file | awk -F \/ '{printf("%s\n",$9)}'`
    ln -s /gpfs/mnt/gpfs02/eic/wfan/data/$FOLDER/$fname .
    fno=`echo $fname | awk -F \_ '{printf("%s\n",$2)}' | awk -F \. '{printf("%s\n",$1)}'`
    echo $file
    echo file name is $fname with file number $fno
    hname=`echo hists-$fno.root`
    echo histogram name is $hname

    if [ -a /gpfs/mnt/gpfs02/eic/wfan/HF_vtx_study/hist_hadron_gen_eC/$hname ]
    then
      echo "File already exists"
    else
      root -b -q 'ana_hadron_gen.C("'$file'","'$hname'",'$NEVT')'
      echo "job done...move output file: ${hname}"
      mv $hname /gpfs/mnt/gpfs02/eic/wfan/HF_vtx_study/hist_hadron_gen_eC/
    fi
  fi
  NUM=$(( $NUM + 1 ))
done

popd
