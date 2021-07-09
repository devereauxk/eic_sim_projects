#!/usr/bin/bash
#-----#-----#----#----#----#----#----#----#----#

if [ -z "$1" ]
then
        echo "No job number set."
        echo "Please run as ./run_eD.sh jobnumber"
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
ln -s /eic/data/baraks/BeAGLE/inputFiles/eD.inp
ln -s /eic/data/baraks/BeAGLE/inputFiles/S1ALL003
ln -s /eic/data/baraks/BeAGLE/nuclear.bin
ln -s /eic/data/baraks/BeAGLE/make_tree.C

#Run simulation
echo "start running in directory $PWD"

echo "Running Job Number $1"
$BEAGLESYS/BeAGLE < eD.inp > eD.log
echo "Completed Simulation!!!"

echo ""

echo "Making Output ROOT File..."
root -l -b -q "make_tree.C(\"eD.txt\")"
echo "Done!!!"
echo ""

#Move output files and cleanup
echo "Cleaning Up..."
mv -v eD.txt /eic/data/baraks/BeAGLE/outForPythiaMode/5_41/eD/eD_${INPUT}.txt
mv -v eD.root /eic/data/baraks/BeAGLE/outForPythiaMode/5_41/eD/eD_${INPUT}.root
mv -v eD.log /eic/data/baraks/BeAGLE/logs/5_41/eD/eD_${INPUT}.log
echo "DONE!!!"
