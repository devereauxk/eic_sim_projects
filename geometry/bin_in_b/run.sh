#const int sys_bins = 5;
#const char* sys_name[sys_bins] = {"e+p", "e+Au [0,2.5]", "e+Au [2.5,5]", "e+Au [5,7.5]", "e+Au[7.5,10]"};
#const char* sys_abbr[sys_bins] = {"ep", "eAu_0_25", "eAu_25_5", "eAu_5_75", "eAu_75_10"};
#
#const int energy_bins = 6;
#const char* energy_name[energy_bins] = {"5x41 GeV", "10x100 GeV", "10x110 GeV", "18x110 GeV", "18x275 GeV", "27.6x0 GeV"};
#const char* energy_abbr[energy_bins] = {"5_41", "10_100", "10_110", "18_110", "18_275", "27.6_0"};

# RUN THIS LINE BEFORE TO MAKE APPROPRIATE ana_merged_in_thickness FILES
# root -l -q "ana_hadron_gen.C(\"../../reconstruction/BeAGLE_v102/eAu_10_100_taufor05_qhat0_nlo/outForPythiaMode/merged.root\", \"../../reconstruction/BeAGLE_v102/eAu_10_100_taufor05_qhat0_nlo/outForPythiaMode/ana_merged_in_b_0_25.root\", 0, 0, 1, 0, 2.5)"

EP_DIR=../../reconstruction/BeAGLE_v101/ep_10_100_baseline_parp2
EA_DIR=../../reconstruction/BeAGLE_v102/eAu_10_100_qhat0_nlo

mkdir $EA_DIR/figs_norm_1/figs_in_b

echo "generating histogram plots for eA/ep in appropriate b bin (writing to figs_in_b/<subfolder> folder)..."
mkdir figs
mkdir $EA_DIR/figs_norm_1/figs_in_b/0_25
root -l -q "plot_chadron_gen.C(\"$EP_DIR/outfiles/ana_merged.root\", 0, 1, \"$EA_DIR/outForPythiaMode/ana_merged_in_b_0_25.root\", 1, 1, \"$EA_DIR/figs_norm_1/figs_in_b/0_25/hists_gen_eAu_0_25.root\", 1)"
mv figs/* $EA_DIR/figs_norm_1/figs_in_b/0_25
rm -r figs
cp $EA_DIR/figs_norm_1/figs_in_b/0_25/hists_gen_eAu_0_25.root ./

echo "generating histogram plots for eA/ep in appropriate b bin (writing to figs_in_b/<subfolder> folder)..."
mkdir figs
mkdir $EA_DIR/figs_norm_1/figs_in_b/25_5
root -l -q "plot_chadron_gen.C(\"$EP_DIR/outfiles/ana_merged.root\", 0, 1, \"$EA_DIR/outForPythiaMode/ana_merged_in_b_25_5.root\", 2, 1, \"$EA_DIR/figs_norm_1/figs_in_b/25_5/hists_gen_eAu_25_5.root\", 1)"
mv figs/* $EA_DIR/figs_norm_1/figs_in_b/25_5
rm -r figs
cp $EA_DIR/figs_norm_1/figs_in_b/25_5/hists_gen_eAu_25_5.root ./

echo "generating histogram plots for eA/ep in appropriate b bin (writing to figs_in_b/<subfolder> folder)..."
mkdir figs
mkdir $EA_DIR/figs_norm_1/figs_in_b/5_75
root -l -q "plot_chadron_gen.C(\"$EP_DIR/outfiles/ana_merged.root\", 0, 1, \"$EA_DIR/outForPythiaMode/ana_merged_in_b_5_75.root\", 3, 1, \"$EA_DIR/figs_norm_1/figs_in_b/5_75/hists_gen_eAu_5_75.root\", 1)"
mv figs/* $EA_DIR/figs_norm_1/figs_in_b/5_75
rm -r figs
cp $EA_DIR/figs_norm_1/figs_in_b/5_75/hists_gen_eAu_5_75.root ./

echo "generating histogram plots for eA/ep in appropriate b bin (writing to figs_in_b/<subfolder> folder)..."
mkdir figs
mkdir $EA_DIR/figs_norm_1/figs_in_b/75_10
root -l -q "plot_chadron_gen.C(\"$EP_DIR/outfiles/ana_merged.root\", 0, 1, \"$EA_DIR/outForPythiaMode/ana_merged_in_b_75_10.root\", 4, 1, \"$EA_DIR/figs_norm_1/figs_in_b/75_10/hists_gen_eAu_75_10.root\", 1)"
mv figs/* $EA_DIR/figs_norm_1/figs_in_b/75_10
rm -r figs
cp $EA_DIR/figs_norm_1/figs_in_b/75_10/hists_gen_eAu_75_10.root ./

mkdir figs
mkdir $EA_DIR/figs_norm_1/figs_in_b/compare_figs
root -l -q "compare_chadron_gen_diff_sys.C(1)"
mv figs/* $EA_DIR/figs_norm_1/figs_in_b/compare_figs
rm -r figs
