#!/usr/bin/bash

echo "-----------------------------------"
echo "Running analysis!!!"
echo "..."
echo ""


echo "Making Output ROOT File with histogram data for ep_10_100 events..."
root -l -b -q 'access_tree("../pythia/outfiles/ep_10_100_norad.root", "ep_10_100/output.root")'
echo "-----------------------------------"
echo ""

echo "Printing histograms for ep_10_100 events..."
root -l -b -q 'plot_histogram("ep_10_100/output.root", "ep_10_100/")'
echo "-----------------------------------"
echo ""

echo "Making Output ROOT File with histogram data for eAu_18_110 events..."
root -l -b -q 'access_tree("../beagle/outForPythiaMode/eAu.root", "eAu_18_110/output.root")'
echo "-----------------------------------"
echo ""

echo "Printing histograms for eAu_18_110 events..."
root -l -b -q 'plot_histogram("eAu_18_110/output.root", "eAu_18_110/")'
echo "-----------------------------------"
echo ""


echo "Done!!!"
