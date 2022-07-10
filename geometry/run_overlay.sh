#const char* sys_name[sys_bins] = {"e+p", "e+Au", "e+Cu", "e+Ca", "e+C", "e+p (BeAGLE)", "e+D", "e+Au (pythia)", "e+Xe", "e+C (pythia)", "e+Pb"};
#const char* sys_abbr[sys_bins] = {"ep", "eAu", "eCu", "eCa", "eC", "ep_BeAGLE", "eD", "eAu_pythia", "eXe", "eC_pythia", "ePb"};

END_DIR=overlay_taufor10

mkdir $END_DIR

DIR=../reconstruction/BeAGLE_v102/eC_10_100_qhat0_nlo
cp $DIR/inc_figs/inc_hists_gen_eC.root $END_DIR

DIR=../reconstruction/BeAGLE_v102/eCu_10_100_qhat0_nlo
cp $DIR/inc_figs/inc_hists_gen_eCu.root $END_DIR

DIR=../reconstruction/BeAGLE_v102/eAu_10_100_qhat0_nlo
cp $DIR/inc_figs/inc_hists_gen_eAu.root $END_DIR

DIR=../reconstruction/BeAGLE_v102/ePb_10_100_qhat0_nlo
cp $DIR/inc_figs/inc_hists_gen_ePb.root $END_DIR
