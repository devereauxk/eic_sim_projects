#const char* sys_name[sys_bins] = {"e+p", "e+Au", "e+Cu", "e+Ca", "e+C", "e+p (BeAGLE)", "e+D", "e+Au (pythia)", "e+Xe", "e+C (pythia)", "e+Pb"};
#const char* sys_abbr[sys_bins] = {"ep", "eAu", "eCu", "eCa", "eC", "ep_BeAGLE", "eD", "eAu_pythia", "eXe", "eC_pythia", "ePb"};

OUT_DIR=./overlay_taufor05
NORM=0

FIG_DIR=$OUT_DIR/figs_norm_${NORM}
mkdir $FIG_DIR

EP_DIR=./BeAGLE_v101/ep_10_100_baseline_parp2
EA_DIR=./BeAGLE_v102/eC_10_100_taufor05_qhat0_nlo

mkdir figs
mkdir $EA_DIR/figs_norm_${NORM}
root -l -q "plot_chadron_gen.C(\"$EP_DIR/outfiles/ana_merged.root\", 0, 1, \"$EA_DIR/outForPythiaMode/ana_merged.root\", 4, 1, \"$EA_DIR/figs_norm_${NORM}/hists_gen.root\", 1)"
mv figs/* $EA_DIR/figs_norm_${NORM}
rm -r figs
cp $EA_DIR/figs_norm_${NORM}/hists_gen.root $OUT_DIR/hists_gen_eC.root


EP_DIR=./BeAGLE_v101/ep_10_100_baseline_parp2
EA_DIR=./BeAGLE_v102/eCu_10_100_taufor05_qhat0_nlo

mkdir figs
mkdir $EA_DIR/figs_norm_${NORM}
root -l -q "plot_chadron_gen.C(\"$EP_DIR/outfiles/ana_merged.root\", 0, 1, \"$EA_DIR/outForPythiaMode/ana_merged.root\", 2, 1, \"$EA_DIR/figs_norm_${NORM}/hists_gen.root\", 1)"
mv figs/* $EA_DIR/figs_norm_${NORM}
rm -r figs
cp $EA_DIR/figs_norm_${NORM}/hists_gen.root $OUT_DIR/hists_gen_eCu.root


EP_DIR=./BeAGLE_v101/ep_10_100_baseline_parp2
EA_DIR=./BeAGLE_v102/eAu_10_100_taufor05_qhat0_nlo

mkdir figs
mkdir $EA_DIR/figs_norm_${NORM}
root -l -q "plot_chadron_gen.C(\"$EP_DIR/outfiles/ana_merged.root\", 0, 1, \"$EA_DIR/outForPythiaMode/ana_merged.root\", 1, 1, \"$EA_DIR/figs_norm_${NORM}/hists_gen.root\", 1)"
mv figs/* $EA_DIR/figs_norm_${NORM}
rm -r figs
cp $EA_DIR/figs_norm_${NORM}/hists_gen.root $OUT_DIR/hists_gen_eAu.root


EP_DIR=./BeAGLE_v101/ep_10_100_baseline_parp2
EA_DIR=./BeAGLE_v102/ePb_10_100_taufor05_qhat0_nlo

mkdir figs
mkdir $EA_DIR/figs_norm_${NORM}
root -l -q "plot_chadron_gen.C(\"$EP_DIR/outfiles/ana_merged.root\", 0, 1, \"$EA_DIR/outForPythiaMode/ana_merged.root\", 10, 1, \"$EA_DIR/figs_norm_${NORM}/hists_gen.root\", 1)"
mv figs/* $EA_DIR/figs_norm_${NORM}
rm -r figs
cp $EA_DIR/figs_norm_${NORM}/hists_gen.root $OUT_DIR/hists_gen_ePb.root

mkdir figs
root -l -q "compare_chadron_gen_diff_sys.C(1)"
mv figs/* $FIG_DIR
rm -r figs
