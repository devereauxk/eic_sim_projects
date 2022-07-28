#!/usr/bin/bash
#-----#-----#----#----#----#----#----#----#----#

WORKING_DIR=/eic/u/kdevereaux/work/BeAGLE-debug/eAu_10_100_qhat0_nlo

if [ -z "$1" ]
then
        echo "No job number set."
        echo "Please run as ./run_kdebug.sh jobnumber"
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
cd $DIR

#Environmental Variables
source /cvmfs/sphenix.opensciencegrid.org/gcc-8.3/opt/sphenix/core/gcc/8.3.0.1-0a5ad/x86_64-centos7/setup.sh

#Soft links to necessary files
rm ${WORKING_DIR}/logs/kdebug_bin_${INPUT}.txt
ln -s ${WORKING_DIR}/logs/eAu_${INPUT}.log


echo "start running in directory $PWD"

echo "Running Job Number $1"
while read line; do
  # if line of file starts with "@kdebug " then prints whole line to $fout
  first_word=`echo $line | awk '{print $1;}'`
  if [[ $first_word == "@kdebug" ]]
  then
    printf $line >> kdebug_bin.txt
  fi
done < eAu_${INPUT}.log
echo "Completed!!!"

echo ""

#Move output files and cleanup
echo "Cleaning Up..."
mv -v kdebug_bin.txt ${WORKING_DIR}/logs/kdebug_bin_${INPUT}.txt
echo "DONE!!!"
