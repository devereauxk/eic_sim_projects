
INJET_DIR=./analysis/ep_10_100_pythia8_pow1_D0injet

rm $INJET_DIR/merged.root
hadd -j $INJET_DIR/merged.root $INJET_DIR/*.root
root -l -q "plot_charm_eec_hists.C(\"$INJET_DIR/merged.root\", \"$INJET_DIR/\")"

INC_DIR=./analysis/ep_10_100_pythia8_pow1_inclusive

rm $INC_DIR/merged.root
hadd -j $INC_DIR/merged.root $INC_DIR/*.root
root -l -q "plot_charm_eec_hists.C(\"$INC_DIR/merged.root\", \"$INC_DIR/\")"

#root -l -q "plot_charm_eec_hists.C(\"./analysis/ep_10_100_pythia8_D0inpair/merged.root\", \"./analysis/ep_10_100_pythia8_D0inpair/\")"
