EP_DIR=./BeAGLE_v101/ep_10_100_baseline_parp2
EA_DIR=./BeAGLE_v102/eC_10_100_taufor05_qhat0_nlo


echo "merging root eA root files..."
hadd -j $EA_DIR/outForPythiaMode/merged.root $EA_DIR/outForPythiaMode/*.root

#echo "generating hadron histograms for ep (writing to ana_merged.root file)..."
#root -l -q "ana_hadron_gen.C(\"$EP_DIR/outfiles/merged.root\", \"$EP_DIR/outfiles/ana_merged.root\", 0, 0, 0)"

echo "generating hadron histograms for eA (writing to ana_merged.root file)..."
root -l -q "ana_hadron_gen.C(\"$EA_DIR/outForPythiaMode/merged.root\", \"$EA_DIR/outForPythiaMode/ana_merged.root\", 0, 0, 1)"

echo "generating histogram plots for ep/eA (writing to figs/ folder)..."
mkdir figs
mkdir $EA_DIR/figs
root -l -q "plot_chadron_gen.C(\"$EP_DIR/outfiles/ana_merged.root\", 0, 1, \"$EA_DIR/outForPythiaMode/ana_merged.root\", 4, 1, \"$EA_DIR/figs/hists_gen.root\")"
mv figs/* $EA_DIR/figs
rm -r figs

#echo "converting plot files from pdf to png..."
#cp topng.sh $EA_DIR/figs
#bash $EA_DIR/figs/topng.sh
#rm $EA_DIR/figs/topng.sh
