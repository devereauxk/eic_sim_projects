#!/usr/bin/bash

echo "-----------------------------------"
echo "Running analysis!!!"
echo "..."
echo ""


echo "Making Output ROOT File with histogram data for events..."
root -l -b -q 'access_tree_pt_eta_total_TH2D.C()'
echo "-----------------------------------"
echo ""

#echo "Printing histograms for events..."
#root -l -b -q 'plot_histogram.C("output.root")'
#echo "-----------------------------------"
#echo ""

echo "Done!!!"
