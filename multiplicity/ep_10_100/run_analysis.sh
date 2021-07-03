#!/usr/bin/bash

echo "-----------------------------------"
echo "Running analysis!!!"
echo "..."
echo ""


echo "Making Output ROOT File with histogram data for ep_10_100 events..."
root -l -b -q 'access_tree.C("./pythia/outfiles/ep_10_100_norad.root", "output.root")'
echo "-----------------------------------"
echo ""

echo "Printing histograms for ep_10_100 events..."
root -l -b -q 'plot_histogram.C("output.root")'
echo "-----------------------------------"
echo ""

echo "Done!!!"
