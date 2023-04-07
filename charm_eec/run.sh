
DATA_DIR=../pythia8/ep_18_275/outfiles
ANALYSIS_DIR=./analysis/ep_18_275

mkdir $ANALYSIS_DIR
rm $DATA_DIR/merged.root
hadd -j $DATA_DIR/merged.root $DATA_DIR/*.root
root -l -q "plot_charm_eec_hists.C(\"$DATA_DIR/merged.root\", \"$ANALYSIS_DIR/\")"