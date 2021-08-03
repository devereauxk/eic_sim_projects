#!/usr/bin/bash

echo "-----------------------------------"
echo "Running analysis!!!"
echo "..."
echo ""

echo "Running plot_multiplicities_vs_atomic_number ... "
root -l -b -q 'plot_multiplicities_vs_atomic_number.C()'
echo "-----------------------------------"
echo ""


echo "Done!!!"
