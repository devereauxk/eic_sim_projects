#!/usr/bin/bash

echo "-----------------------------------"
echo "Running analysis!!!"
echo "..."
echo ""


echo "Making Output ROOT File with histogram data for events..."
root -l -b -q 'access_tree_mult_eta_binned.C()'
echo "-----------------------------------"
echo ""

echo "Printing histograms for events..."
root -l -b -q 'plot_histograms_mult_eta_binned.C()'
echo "-----------------------------------"
echo ""

echo "Done!!!"
