#!/bin/bash
OUT_DIR=/eic/u/kdevereaux/work/eHIJING/eHIJING-pythia/eHIJING-examples

if [ -z "$1" ]
then
        echo "No job number set."
        echo "Please run as ./run_ep.sh jobnumber"
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

export PYTHIA8DATA=/eic/u/kdevereaux/work/eHIJING/local/share/Pythia8/xmldoc
export PYTHIA8=/eic/u/kdevereaux/work/eHIJING/local

exe=/eic/u/kdevereaux/work/eHIJING/eHIJING-pythia/eHIJING-examples/build/ehijing-test
Configfile=/eic/u/kdevereaux/work/eHIJING/eHIJING-pythia/eHIJING-examples/s20.setting
Neve=400000

K=4 # default: 4.0 for EIC
M=1 # Generlizaed HT:1,  HIgher-twist:0, both in the soft gluon emission limit.

label="eHe4_10_100_pdf0"
folder=$OUT_DIR/Events/$label
TablePath=$OUT_DIR/Tables/$label
mkdir -p $folder
mkdir -p $TablePath
Z=2
A=4
$exe $Neve $Z $A $M $K $DIR $DIR $Configfile # > /dev/null 2>&1

mv -v $DIR/events.dat $folder/eHe4_${INPUT}.dat
mv -v $DIR/GHT.dat $TablePath/GHT_${INPUT}.dat
mv -v $DIR/Qs.dat $TablePath/Qs_${INPUT}.dat
