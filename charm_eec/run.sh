
INJET_DIR=./analysis/ep_10_100_pythia8_D0injet

rm $INJET_DIR/merged.root
hadd -j $INJET_DIR/merged.root $INJET_DIR/*.root
root -l -q "plot_charm_eec_hists.C(\"./analysis/ep_10_100_pythia8_D0injet/merged.root\", \"./analysis/ep_10_100_pythia8_D0injet/\")"


#root -l -q "plot_charm_eec_hists.C(\"./analysis/ep_10_100_pythia8_D0inpair/merged.root\", \"./analysis/ep_10_100_pythia8_D0inpair/\")"

#root -l -q "plot_charm_eec_hists.C(\"../energy_correlator/eHIJING/ep_10_100_K0_pow05/merged.root\", \"./analysis/ep_10_100_pythia8_inclusive/\")"
