root -l -q "plot_eec_hists.C(\"./eHIJING/eAu_1E8_K0_condor_v2/merged.root\", \"./eHIJING/eAu_1E8_K0_condor_v2/\", 1, \"./eHIJING/eAu_1E8_K0_condor_v2/merged.root\")"

root -l -q "plot_eec_hists.C(\"./eHIJING/eAu_1E8_condor_v2/merged.root\", \"./eHIJING/eAu_1E8_condor_v2/\", 1, \"./eHIJING/eAu_1E8_K0_condor_v2/merged.root\")"

root -l -q "plot_eec_hists.C(\"./eHIJING/eAu_1E8_K2/merged.root\", \"./eHIJING/eAu_1E8_K2/\", 1, \"./eHIJING/eAu_1E8_K0_condor_v2/merged.root\")"

root -l -q "plot_eec_hists.C(\"./eHIJING/eAu_1E8_K4/merged.root\", \"./eHIJING/eAu_1E8_K4/\", 1, \"./eHIJING/eAu_1E8_K0_condor_v2/merged.root\")"

root -l -q "plot_eec_hists.C(\"./eHIJING/eAu_1E8_K0_pow2/merged.root\", \"./eHIJING/eAu_1E8_K0_pow2/\")"

root -l -q "plot_eec_hists.C(\"./eHIJING/eAu_1E8_K4_pow2/merged.root\", \"./eHIJING/eAu_1E8_K4_pow2/\", 1, \"./eHIJING/eAu_1E8_K0_pow2/merged.root\")"

root -l -q "plot_eec_hists.C(\"./eHIJING/eC_1E8_K4/merged.root\", \"./eHIJING/eC_1E8_K4/\", 1, \"./eHIJING/eC_1E8_K0/merged.root\")"

root -l -q "plot_eec_hists.C(\"./eHIJING/eCu_1E8_K4/merged.root\", \"./eHIJING/eCu_1E8_K4/\", 1, \"./eHIJING/eCu_1E8_K0/merged.root\")"

# above deprecated since they use K=0 instead of ep as the baseline

# baseline comparisons
root -l -q "plot_eec_hists.C(\"./eHIJING/eD_10_100_K0/merged.root\", \"./eHIJING/eD_10_100_K0/\", 1, \"./eHIJING/ep_10_100_K0/merged.root\")"

root -l -q "plot_eec_hists.C(\"./eHIJING/eD_10_100_K4/merged.root\", \"./eHIJING/eD_10_100_K4/\", 1, \"./eHIJING/ep_10_100_K0/merged.root\")"

root -l -q "plot_eec_hists.C(\"./eHIJING/eAu_1E8_K0_condor_v2/merged.root\", \"./eHIJING/eAu_1E8_K0_condor_v2/\", 1, \"./eHIJING/ep_10_100_K0/merged.root\")"

# hists for eHIJIJNG and PYHTIA8 comparison
root -l -q "plot_eec_hists.C(\"./eHIJING/ep_10_100_K0/merged.root\", \"./eHIJING/ep_10_100_K0/\")"

root -l -q "plot_eec_hists.C(\"./eHIJING/ep_10_100_pythia8/merged.root\", \"./eHIJING/ep_10_100_pythia8/\")"

root -l -q "plot_eec_hists.C(\"./eHIJING/ep_10_100_pythia8_ft/merged.root\", \"./eHIJING/ep_10_100_pythia8/\")"
