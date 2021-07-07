#!/usr/bin/bash
#-----#-----#----#----#----#----#----#----#----#

if [ -z "$1" ]
then
        echo "No job number set."
        echo "Please run as ./run_eC.sh jobnumber"
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
ln -s /eic/data/baraks/BeAGLE/inputFiles/eC.inp
ln -s /eic/data/baraks/BeAGLE/inputFiles/S3ALL003
ln -s /eic/data/baraks/BeAGLE/nuclear.bin
ln -s /eic/data/baraks/BeAGLE/make_tree.C

#Run simulation
echo "start running in directory $PWD"

echo "Running Job Number $1"
$BEAGLESYS/BeAGLE < eC.inp > eC.log
echo "Completed Simulation!!!"

echo ""

echo "Making Output ROOT File..."
root -l -b -q "make_tree.C(\"eA.txt\")"
echo "Done!!!"
echo ""

#Move output files and cleanup
echo "Cleaning Up..."
mv -v eA.txt /eic/data/baraks/BeAGLE/outForPythiaMode/5_41/eC/eC_${INPUT}.txt
mv -v eA.root /eic/data/baraks/BeAGLE/outForPythiaMode/5_41/eC/eC_${INPUT}.root
mv -v eC.log /eic/data/baraks/BeAGLE/logs/5_41/eC/eC_${INPUT}.log
echo "DONE!!!"
