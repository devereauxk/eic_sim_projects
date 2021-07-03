#!/usr/bin/bash

echo "-----------------------------------"
echo "Running analysis!!!"
echo "..."
echo ""


echo "Making Output ROOT File with histogram data for eD_10_110 events..."
root -l -b -q 'access_tree.C("./beagle/outForPythiaMode/eD.root", "output.root")'
echo "-----------------------------------"
echo ""

echo "Printing histograms for eAu_18_110 events..."
root -l -b -q 'plot_histogram.C("output.root")'
echo "-----------------------------------"
echo ""


echo "Done!!!"
