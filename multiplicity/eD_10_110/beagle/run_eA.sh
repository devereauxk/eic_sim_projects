#!/usr/bin/bash

echo "-----------------------------------"
echo "Running BeAGLE Simulation for eA Collider!!!"
echo "..."
echo ""


$BEAGLESYS/BeAGLE < input/eD.inp > logs/eD.log

mv eD.txt outForPythiaMode/

echo "Making Output ROOT File..."
root -l -b -q 'make_tree.C("eD.txt")'
echo "-----------------------------------"


echo "Done!!!"
