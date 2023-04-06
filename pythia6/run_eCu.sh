#!/usr/bin/bash
#-----#-----#----#----#----#----#----#----#----#

WORKING_DIR=/eic/u/kdevereaux/work/reconstruction/BeAGLE_v102/eCu_10_100_taufor05_qhat0_nlo

if [ -z "$1" ]
then
        echo "No job number set."
        echo "Please run as ./run_eCu.sh jobnumber"
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
export BEAGLESYS="${EICDIRECTORY}/PACKAGES/BeAGLE"
export LHAPDF5="${EICDIRECTORY}/gcc-8.3/lib"
export LD_LIBRARY_PATH="${LHAPDF5}:$LD_LIBRARY_PATH"
source /cvmfs/sphenix.opensciencegrid.org/gcc-8.3/opt/sphenix/core/gcc/8.3.0.1-0a5ad/x86_64-centos7/setup.sh

#Soft links to necessary files
ln -s ${WORKING_DIR}/inputFiles/eCu.inp
ln -s ${WORKING_DIR}/inputFiles/S1ALL003
ln -s ${WORKING_DIR}/nuclear.bin
ln -s ${WORKING_DIR}/make_tree.C

#Run simulation
echo "start running in directory $PWD"

echo "Running Job Number $1"
$BEAGLESYS/BeAGLE < eCu.inp > eCu.log
echo "Completed Simulation!!!"

echo ""

echo "Making Output ROOT File..."
root -l -b -q "make_tree.C(\"eCu.txt\")"
echo "Done!!!"
echo ""

#Move output files and cleanup
echo "Cleaning Up..."
mv -v eCu.txt ${WORKING_DIR}/outForPythiaMode/eCu_${INPUT}.txt
mv -v eCu.root ${WORKING_DIR}/outForPythiaMode/eCu_${INPUT}.root
mv -v eCu.log ${WORKING_DIR}/logs/eCu_${INPUT}.log
echo "DONE!!!"
