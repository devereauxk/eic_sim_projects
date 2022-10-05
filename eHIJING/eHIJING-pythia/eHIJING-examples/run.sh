#!/bin/bash
export PYTHIA8DATA=/global/homes/w/wk42/miniconda3/envs/ehijing/share/Pythia8/xmldoc
export PYTHIA8=/global/homes/w/wk42/miniconda3/envs/ehijing

exe=ehijing-test
Neve=100000
Configfile=s20.setting

K=2.0
M=1 # Generlizaed HT:1,  HIgher-twist:0, both in the soft gluon emission limit.

label="ed"
folder="Events"/$label
TablePath="Tables/$K"
mkdir -p $folder
mkdir -p $TablePath
Z=1
A=2
$exe $Neve $Z $A $M $K $TablePath $folder $Configfile # > /dev/null 2>&1

label="eA"
folder="Events"/$label
TablePath="Tables/$K"
mkdir -p $folder
mkdir -p $TablePath
Z=82
A=208
$exe $Neve $Z $A $M $K $TablePath $folder $Configfile # > /dev/null 2>&1

