#!/usr/bin/bash
#-----#-----#----#----#----#----#----#----#----#

WORKING_DIR=/eic/u/kdevereaux/work/reconstruction/BeAGLE_v102/eAu_10_100_gen10_qhat0_nlo

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

#Environmental Variables
#export BEAGLESYS="${EICDIRECTORY}/PACKAGES/BeAGLE" # manual overide
export BEAGLESYS="${EICDIRECTORY}/PACKAGES/BeAGLE"
export LHAPDF5="${EICDIRECTORY}/gcc-8.3/lib"
export LD_LIBRARY_PATH="${LHAPDF5}:$LD_LIBRARY_PATH"
source /cvmfs/sphenix.opensciencegrid.org/gcc-8.3/opt/sphenix/core/gcc/8.3.0.1-0a5ad/x86_64-centos7/setup.sh

#Soft links to necessary files
#ln -s /eic/data/baraks/BeAGLE/inputFiles/eAu.inp
#ln -s /eic/data/baraks/BeAGLE/inputFiles/S3ALL003
#ln -s /eic/data/baraks/BeAGLE/nuclear.bin
#ln -s /eic/data/baraks/BeAGLE/make_tree.C
ln -s ${WORKING_DIR}/inputFiles/eAu.inp
ln -s ${WORKING_DIR}/inputFiles/S1ALL003
ln -s ${WORKING_DIR}/nuclear.bin
ln -s ${WORKING_DIR}/make_tree.C


#Run simulation
echo "start running in directory $PWD"

echo "Running Job Number $1"
$BEAGLESYS/BeAGLE < eAu.inp > eAu.log
#$BEAGLESYS/BeAGLE-v1.00.00 < eAu.inp > eAu.log
echo "Completed Simulation!!!"

echo ""

echo "Making Output ROOT File..."
root -l -b -q "make_tree.C(\"eAu.txt\")"
echo "Done!!!"
echo ""

#Move output files and cleanup
echo "Cleaning Up..."
#mv -v eA.txt /eic/data/baraks/BeAGLE/outForPythiaMode/5_41/eAu/eAu_${INPUT}.txt
#mv -v eA.root /eic/data/baraks/BeAGLE/outForPythiaMode/5_41/eAu/eAu_${INPUT}.root
#mv -v eAu.log /eic/data/baraks/BeAGLE/logs/5_41/eAu/eAu_${INPUT}.log
mv -v eAu.txt ${WORKING_DIR}/outForPythiaMode/eAu_${INPUT}.txt
mv -v eAu.root ${WORKING_DIR}/outForPythiaMode/eAu_${INPUT}.root
mv -v eAu.log ${WORKING_DIR}/logs/eAu_${INPUT}.log
echo "DONE!!!"
