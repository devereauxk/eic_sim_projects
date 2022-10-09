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
cd $DIR

export PYTHIA8DATA=/eic/u/kdevereaux/work/eHIJING/local/share/Pythia8/xmldoc
export PYTHIA8=/eic/u/kdevereaux/work/eHIJING/local

exe=/eic/u/kdevereaux/work/eHIJING/eHIJING-pythia/eHIJING-examples/build/ehijing-test
Configfile=/eic/u/kdevereaux/work/eHIJING/eHIJING-pythia/eHIJING-examples/s20.setting
Neve=10000 #10000

K=2.0
M=1 # Generlizaed HT:1,  HIgher-twist:0, both in the soft gluon emission limit.

label="ep_1E6_condor"
folder=$OUT_DIR/Events/$label
TablePath=$OUT_DIR/Tables/$K
local_folder=./Events
local_TablePath=./Tables
mkdir -p $folder
mkdir -p $TablePath
Z=1
A=1
$exe $Neve $Z $A $M $K $TablePath $folder $Configfile # > /dev/null 2>&1
