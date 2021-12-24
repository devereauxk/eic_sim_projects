#!/bin/sh

echo "This script construct ratio of hadron production (truth level) for two systems."
echo "Numerator system is usually e+A and denominator sysmtem is usually e+p."
# Exe line to obtain e+Au @ 18x110 / e+p @ 10x100 is like below:
# bash check_chadron_vs_z_gen.sh 1 3 0 1

sys_name=("e+p" "e+Au" "e+Cu" "e+Ca" "e+C" "e+p (BeAGLE)" "e+d")
sys_abbr=("ep" "eAu" "eCu" "eCa" "eC" "ep_BeAGLE" "ed")

energy_name=("5x41" "10x100" "10x110" "18x110" "18x275")
energy_abbr=("5_41" "10_100" "10_110" "18_110" "18_275")

num_sys=$(( $1 ))
echo "Numerator collision system is ${sys_name[$num_sys]}"

num_energy=$(( $2 ))
echo "Numerator collision beam energy is ${energy_name[$num_energy]} GeV"

denom_sys=$(( $3 ))
echo "Denominator collision system is ${sys_name[$denom_sys]}"

denom_energy=$(( $4 ))
echo "Denominator collision beam energy is ${energy_name[$denom_energy]} GeV"

mkdir kin_gen
pushd kin_gen

mkdir ${sys_abbr[$num_sys]}_${energy_abbr[$num_energy]}_over_${sys_abbr[$denom_sys]}_${energy_abbr[$denom_energy]}
pushd ${sys_abbr[$num_sys]}_${energy_abbr[$num_energy]}_over_${sys_abbr[$denom_sys]}_${energy_abbr[$denom_energy]}

rm -rf figs
mkdir figs

# ln -sf ~/EIC/fast_sim/chadron/ep_data/10_100/hists_hadron_gen_e10p100_ycut_v2.root hists_gen_ep.root

ln -sf ~/EIC/fast_sim/chadron/${sys_abbr[$denom_sys]}_data/${energy_abbr[$denom_energy]}/hists_ycut_v3.root hists_gen_ep.root
ln -sf ~/EIC/fast_sim/chadron/${sys_abbr[$num_sys]}_data/${energy_abbr[$num_energy]}/hists_ycut_v3.root hists_gen_eA.root

# ln -sf ~/EIC/data/hists_gen_ep.root hists_gen_ep.root
# ln -sf ~/EIC/data/hists_gen_eAu.root hists_gen_eAu.root

# ln -sf ~/EIC/fast_sim/chadron/ep_data/5_41/hists_hadron_gen_e5p41.root hists_gen_ep.root
# ln -sf ~/EIC/fast_sim/chadron/eAu_data/5_41/hists_hadron_gen_e5Au41.root hists_gen_eAu.root


# ln -sf ~/EIC/fast_sim/chadron/ep_data/10_100/hists_hadron_BeAGLE_e10p100_ycut.root hists_gen_eAu.root
# ln -sf ~/EIC/fast_sim/chadron/eAu_data/10_110/hists_hadron_gen_e10Au110_ycut_v2.root hists_gen_eA.root

ln -sf ~/EIC/fast_sim/chadron/plot_chadron_gen.C .
ln -sf ~/EIC/fast_sim/chadron/bins.h .
ln -sf ~/Documents/stuff/Pad3x3.h .
ln -sf ~/Documents/stuff/Pad4x5.h .

root -b -q 'plot_chadron_gen.C("hists_gen_ep.root", '${denom_sys}', '${denom_energy}', "hists_gen_eA.root", '${num_sys}', '${num_energy}', "hists_gen.root")'

popd

popd

