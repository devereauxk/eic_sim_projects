#!/usr/bin/bash
#-----#-----#----#----#----#----#----#----#----#

if [ -z "$1" ]
then
        echo "No job number set."
        echo "Please run as ./run_ep_BeAGLE.sh jobnumber"
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
ln -s /eic/u/kdevereaux/work/multiplicity/condor/ep_10_100_BeAGLE/inputFiles/ep.inp
ln -s /eic/u/kdevereaux/work/multiplicity/condor/ep_10_100_BeAGLE/inputFiles/S1ALL003
ln -s /eic/u/kdevereaux/work/multiplicity/condor/ep_10_100_BeAGLE/nuclear.bin
ln -s /eic/u/kdevereaux/work/multiplicity/condor/ep_10_100_BeAGLE/make_tree.C

#Run simulation
echo "start running in directory $PWD"

echo "Running Job Number $1"
$BEAGLESYS/BeAGLE < ep.inp > ep.log
echo "Completed Simulation!!!"

echo ""

echo "Making Output ROOT File..."
root -l -b -q "make_tree.C(\"ep.txt\")"
echo "Done!!!"
echo ""

#Move output files and cleanup
echo "Cleaning Up..."
mv -v ep.txt /eic/u/kdevereaux/work/multiplicity/condor/ep_10_100_BeAGLE/outForPythiaMode/ep_${INPUT}.txt
mv -v ep.root /eic/u/kdevereaux/work/multiplicity/condor/ep_10_100_BeAGLE/outForPythiaMode/ep_${INPUT}.root
mv -v ep.log /eic/u/kdevereaux/work/multiplicity/condor/ep_10_100_BeAGLE/logs/ep_${INPUT}.log
echo "DONE!!!"
