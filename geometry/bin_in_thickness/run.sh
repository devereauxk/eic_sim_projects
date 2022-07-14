o#const int sys_bins = 5;
#const char* sys_name[sys_bins] = {"e+p", "e+Au [0,3]", "e+Au [3,6]", "e+Au [6,9]", "e+Au[9,13]"};
#const char* sys_abbr[sys_bins] = {"ep", "eAu_0_3", "eAu_3_6", "eAu_6_9", "eAu_9_13"};
#
#const int energy_bins = 6;
#const char* energy_name[energy_bins] = {"5x41 GeV", "10x100 GeV", "10x110 GeV", "18x110 GeV", "18x275 GeV", "27.6x0 GeV"};
#const char* energy_abbr[energy_bins] = {"5_41", "10_100", "10_110", "18_110", "18_275", "27.6_0"};

# RUN THIS LINE BEFORE TO MAKE APPROPRIATE ana_merged_in_thickness FILES
# root -l -q "ana_hadron_gen.C(\"../../reconstruction/BeAGLE_v102/eAu_10_100_taufor05_qhat0_nlo/outForPythiaMode/merged.root\", \"../../reconstruction/BeAGLE_v102/eAu_10_100_taufor05_qhat0_nlo/outForPythiaMode/ana_merged_in_thickness_0_3.root\", 0, 0, 1, 0, 3)"

EP_DIR=../../reconstruction/BeAGLE_v101/ep_10_100_baseline_parp2
EA_DIR=../../reconstruction/BeAGLE_v102/eAu_10_100_qhat0_nlo

SUB_DIR=$EA_DIR/figs_norm_1/figs_in_thickness_no_eta

mkdir $SUB_DIR

echo "generating histogram plots for eA/ep in appropriate thickness bin (writing to figs_in_thickness/<subfolder> folder)..."
mkdir figs
mkdir $SUB_DIR/0_3
root -l -q "plot_chadron_gen.C(\"$EP_DIR/outfiles/ana_merged.root\", 0, 1, \"$EA_DIR/outForPythiaMode/ana_merged_in_thickness_0_3.root\", 1, 1, \"$SUB_DIR/0_3/hists_gen_eAu_0_3.root\", 1)"
mv figs/* $SUB_DIR/0_3
rm -r figs
cp $SUB_DIR/0_3/hists_gen_eAu_0_3.root ./

echo "generating histogram plots for eA/ep in appropriate thickness bin (writing to figs_in_thickness/<subfolder> folder)..."
mkdir figs
mkdir $SUB_DIR/3_6
root -l -q "plot_chadron_gen.C(\"$EP_DIR/outfiles/ana_merged.root\", 0, 1, \"$EA_DIR/outForPythiaMode/ana_merged_in_thickness_3_6.root\", 2, 1, \"$SUB_DIR/3_6/hists_gen_eAu_3_6.root\", 1)"
mv figs/* $SUB_DIR/3_6
rm -r figs
cp $SUB_DIR/3_6/hists_gen_eAu_3_6.root ./

echo "generating histogram plots for eA/ep in appropriate thickness bin (writing to figs_in_thickness/<subfolder> folder)..."
mkdir figs
mkdir $SUB_DIR/6_9
root -l -q "plot_chadron_gen.C(\"$EP_DIR/outfiles/ana_merged.root\", 0, 1, \"$EA_DIR/outForPythiaMode/ana_merged_in_thickness_6_9.root\", 3, 1, \"$SUB_DIR/6_9/hists_gen_eAu_6_9.root\", 1)"
mv figs/* $SUB_DIR/6_9
rm -r figs
cp $SUB_DIR/6_9/hists_gen_eAu_6_9.root ./

echo "generating histogram plots for eA/ep in appropriate thickness bin (writing to figs_in_thickness/<subfolder> folder)..."
mkdir figs
mkdir $SUB_DIR/9_13
root -l -q "plot_chadron_gen.C(\"$EP_DIR/outfiles/ana_merged.root\", 0, 1, \"$EA_DIR/outForPythiaMode/ana_merged_in_thickness_9_13.root\", 4, 1, \"$SUB_DIR/9_13/hists_gen_eAu_9_13.root\", 1)"
mv figs/* $SUB_DIR/9_13
rm -r figs
cp $SUB_DIR/9_13/hists_gen_eAu_9_13.root ./

mkdir figs
mkdir $SUB_DIR/compare_figs
root -l -q "compare_chadron_gen_diff_sys.C(1)"
mv figs/* $SUB_DIR/compare_figs
rm -r figs
rm *.root
