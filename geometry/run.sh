#const int sys_bins = 11;
#const char* sys_name[sys_bins] = {"e+p", "e+Au", "e+Cu", "e+Ca", "e+C", "e+p (BeAGLE)", "e+D", "e+Au (pythia)", "e+Xe", "e+C (pythia)", "e+Pb"};
#const char* sys_abbr[sys_bins] = {"ep", "eAu", "eCu", "eCa", "eC", "ep_BeAGLE", "eD", "eAu_pythia", "eXe", "eC_pythia", "ePb"};
#
#const int energy_bins = 6;
#const char* energy_name[energy_bins] = {"5x41 GeV", "10x100 GeV", "10x110 GeV", "18x110 GeV", "18x275 GeV", "27.6x0 GeV"};
#const char* energy_abbr[energy_bins] = {"5_41", "10_100", "10_110", "18_110", "18_275", "27.6_0"};

DIR=../reconstruction/BeAGLE_v102/eC_10_100_qhat0_nlo

#echo "merging root files..."
#hadd -f -j $DIR/outForPythiaMode/merged.root $DIR/outForPythiaMode/*.root

#echo "generating Ninc/Nincch histogram (writing to inc_merged.root file)..."
#root -l -q "inc_hadron_gen.C(\"$DIR/outForPythiaMode/merged.root\", \"$DIR/outForPythiaMode/inc_merged.root\", 0, 0, 1)"

echo "generating histogram plots for species (writing to inc_figs/ folder)..."
mkdir $DIR/inc_figs
root -l -q "plot_inc_hadron.C(\"$DIR/outForPythiaMode/inc_merged.root\", 4, 1, \"$DIR/inc_figs/\")"

DIR=../reconstruction/BeAGLE_v102/eCu_10_100_qhat0_nlo

#echo "merging root files..."
#hadd -f -j $DIR/outForPythiaMode/merged.root $DIR/outForPythiaMode/*.root

#echo "generating Ninc/Nincch histogram (writing to inc_merged.root file)..."
#root -l -q "inc_hadron_gen.C(\"$DIR/outForPythiaMode/merged.root\", \"$DIR/outForPythiaMode/inc_merged.root\", 0, 0, 1)"

echo "generating histogram plots for species (writing to inc_figs/ folder)..."
mkdir $DIR/inc_figs
root -l -q "plot_inc_hadron.C(\"$DIR/outForPythiaMode/inc_merged.root\", 2, 1, \"$DIR/inc_figs/\")"

DIR=../reconstruction/BeAGLE_v102/eAu_10_100_qhat0_nlo

#echo "merging root files..."
#hadd -f -j $DIR/outForPythiaMode/merged.root $DIR/outForPythiaMode/*.root

#echo "generating Ninc/Nincch histogram (writing to inc_merged.root file)..."
#root -l -q "inc_hadron_gen.C(\"$DIR/outForPythiaMode/merged.root\", \"$DIR/outForPythiaMode/inc_merged.root\", 0, 0, 1)"

echo "generating histogram plots for species (writing to inc_figs/ folder)..."
mkdir $DIR/inc_figs
root -l -q "plot_inc_hadron.C(\"$DIR/outForPythiaMode/inc_merged.root\", 1, 1, \"$DIR/inc_figs/\")"

DIR=../reconstruction/BeAGLE_v102/ePb_10_100_qhat0_nlo

#echo "merging root files..."
#hadd -f -j $DIR/outForPythiaMode/merged.root $DIR/outForPythiaMode/*.root

#echo "generating Ninc/Nincch histogram (writing to inc_merged.root file)..."
#root -l -q "inc_hadron_gen.C(\"$DIR/outForPythiaMode/merged.root\", \"$DIR/outForPythiaMode/inc_merged.root\", 0, 0, 1)"

echo "generating histogram plots for species (writing to inc_figs/ folder)..."
mkdir $DIR/inc_figs
root -l -q "plot_inc_hadron.C(\"$DIR/outForPythiaMode/inc_merged.root\", 10, 1, \"$DIR/inc_figs/\")"
