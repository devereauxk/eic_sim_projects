#!/usr/bin/bash
#-----#-----#----#----#----#----#----#----#----#

WORKING_DIR=/eic/u/kdevereaux/work/pythia8/ep_18_275_nodecay

if [ -z "$1" ]
then
        echo "No job number set."
        echo "Please run as ./run_eAu.sh jobnumber"
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

ln -s ${WORKING_DIR}/DIS_18275
ln -s ${WORKING_DIR}/make_tree.C

#Run simulation
echo "start running in directory $PWD"
echo "Running Job Number $1"

./DIS_18275 > DIS_18275.log

echo "Completed Simulation!!!"
echo ""

echo "Making Output ROOT File..."
root -l -b -q 'make_tree.C("hepmcout.dat")'
echo "Done!!!"
echo ""

#Move output files and cleanup
echo "Cleaning Up..."
mv -v hempmcout.root ${WORKING_DIR}/outfiles/hempmcout_${INPUT}.root
mv -v DIS_18275.log ${WORKING_DIR}/logfiles/DIS_18275_${INPUT}.log
echo "DONE!!!"
