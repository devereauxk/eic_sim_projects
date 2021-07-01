#!/usr/bin/bash

echo "-----------------------------------"
echo "Running BeAGLE Simulation for eA Collider!!!"
echo "..."
echo ""

$BEAGLESYS/BeAGLE < inputFiles/eAu.inp > logs/eAu.log
#$BEAGLESYS/BeAGLE < inputFiles/eAu_noFermi.inp > logs/eAu.log

echo "Making Output ROOT File..."
root -l -b -q 'make_tree.C("eAu.txt")'
echo "-----------------------------------"


echo "Done!!!"


