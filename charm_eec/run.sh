
DIR=./analysis/ep_18_275_nodecay_pow1_D0inpair

rm $DIR/merged.root
hadd $DIR/merged.root $DIR/*.root
root -l -q "plot_charm_eec_hists.C(\"$DIR/merged.root\", \"$DIR/\")"

#DIR=./analysis/ep_18_275_nodecay_pow1_D0injet

#rm $DIR/merged.root
#hadd $DIR/merged.root $DIR/*.root
#root -l -q "plot_charm_eec_hists.C(\"$DIR/merged.root\", \"$DIR/\")"

#root -l -q "plot_charm_eec_hists.C(\"./analysis/ep_18_275_pow1_inclusive/merged.root\", \"./analysis/ep_18_275_pow1_inclusive/\")"