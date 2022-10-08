#!/bin/bash
export PYTHIA8DATA=/eic/u/kdevereaux/work/eHIJING/local/share/Pythia8/xmldoc
export PYTHIA8=/eic/u/kdevereaux/work/eHIJING/local

exe=build/ehijing-test
Neve=10000 #10000
Configfile=s20.setting

K=2.0
M=1 # Generlizaed HT:1,  HIgher-twist:0, both in the soft gluon emission limit.

label="ep"
folder="Events"/$label
TablePath="Tables/$K"
mkdir -p $folder
mkdir -p $TablePath
Z=1
A=1
$exe $Neve $Z $A $M $K $TablePath $folder $Configfile # > /dev/null 2>&1
