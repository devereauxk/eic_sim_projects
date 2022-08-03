EP_DIR=./BeAGLE_v101/ep_10_100_baseline_parp2
EA_DIR=./BeAGLE_v102/eC_10_100_tauforOff_qhat0_nlo
#EA_DIR=./BeAGLE_v102/eXe_27.6_0_tauforOff_nlo


echo "merging root eA root files..."
hadd -f -j $EA_DIR/outForPythiaMode/merged.root $EA_DIR/outForPythiaMode/*.root

#echo "generating hadron histograms for ep (writing to ana_merged.root file)..."
#root -l -q "ana_hadron_gen.C(\"$EP_DIR/outfiles/merged.root\", \"$EP_DIR/outfiles/ana_merged.root\", 0, 0, 0)"

echo "generating hadron histograms for eA (writing to ana_merged.root file)..."
root -l -q "ana_hadron_gen.C(\"$EA_DIR/outForPythiaMode/merged.root\", \"$EA_DIR/outForPythiaMode/ana_merged.root\", 0, 0, 1)"

echo "generating histogram plots for eA/ep (writing to figs/ folder)..."
mkdir figs
mkdir $EA_DIR/figs_norm_0
root -l -q "plot_chadron_gen.C(\"$EP_DIR/outfiles/ana_merged.root\", 0, 1, \"$EA_DIR/outForPythiaMode/ana_merged.root\", 4, 1, \"$EA_DIR/figs_norm_0/hists_gen.root\", 0)"
mv figs/* $EA_DIR/figs_norm_0
rm -r figs

mkdir figs
mkdir $EA_DIR/figs_norm_1
root -l -q "plot_chadron_gen.C(\"$EP_DIR/outfiles/ana_merged.root\", 0, 1, \"$EA_DIR/outForPythiaMode/ana_merged.root\", 4, 1, \"$EA_DIR/figs_norm_1/hists_gen.root\", 1)"
mv figs/* $EA_DIR/figs_norm_1
rm -r figs


#const int sys_bins = 11;
#const char* sys_name[sys_bins] = {"e+p", "e+Au", "e+Cu", "e+Ca", "e+C", "e+p (BeAGLE)", "e+D", "e+Au (pythia)", "e+Xe", "e+C (pythia)", "e+Pb"};
#const char* sys_abbr[sys_bins] = {"ep", "eAu", "eCu", "eCa", "eC", "ep_BeAGLE", "eD", "eAu_pythia", "eXe", "eC_pythia", "ePb"};
#
#const int energy_bins = 6;
#const char* energy_name[energy_bins] = {"5x41 GeV", "10x100 GeV", "10x110 GeV", "18x110 GeV", "18x275 GeV", "27.6x0 GeV"};
#const char* energy_abbr[energy_bins] = {"5_41", "10_100", "10_110", "18_110", "18_275", "27.6_0"};
